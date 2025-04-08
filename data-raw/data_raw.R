library(tidyverse)
library(magrittr)
library(rmapshaper)
library(sf)

# # Download HUC4 dataset from https://nrcs.app.box.com/v/gateway/file/233091823212
# unzip("./data-raw/wbdhu4_a_us_september2017.zip",
#   exdir = "./data-raw/wbdhu4_a_us_september2017"
# )
#
# nhd_huc4 <- sf::st_read("./data-raw/wbdhu4_a_us_september2017/wbdhu4_a_us_september2017.gdb",
#   layer = "WBDHU4"
# ) %>%
#   dplyr::filter(
#     !grepl("AK", STATES),
#     !grepl("HI", STATES)
#   ) %>%
#   dplyr::select(dplyr::starts_with("HUC")) %>%
#   dplyr::rename(geometry = Shape) %>%
#   tibble::as_tibble() %>%
#   sf::st_as_sf() %>%
#   sf::st_cast("MULTIPOLYGON")
#
# sf::st_write(nhd_huc4,
#   "./data-raw/nhd_huc4.gpkg",
#   delete_dsn = TRUE
# )
#
# zip("./data-raw/nhd_huc4.gpkg.zip", "./data-raw/nhd_huc4.gpkg")
#
# # devtools::use_data(nhd_huc4,
# #                    overwrite = TRUE)
#
# unlink("./data-raw/wbdhu4_a_us_september2017.zip")
#
# unlink("./data-raw/wbdhu4_a_us_september2017",
#   recursive = TRUE
# )
#
# unlink("./data-raw/nhd_huc4.gpkg",
#   recursive = TRUE
# )


## The NASS CDL RAT
rat <- tempfile(fileext = ".zip")
download.file("https://www.nass.usda.gov/Research_and_Science/Cropland/docs/generic_cdl_attributes.tif.vat.dbf.zip",
  destfile = rat
)
unzip(rat, exdir = tempdir())

nass <-
  foreign::read.dbf(paste0(tempdir(), "/ESRI_attribute_files/ArcGIS10.7.0_2019_30m_cdls.img.vat.dbf")) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(rgb = rgb(red = RED, green = GREEN, blue = BLUE, alpha = OPACITY, maxColorValue = 255)) %>%
  dplyr::transmute(
    ID = VALUE,
    Class = as.character(CLASS_NAME),
    Color = rgb
  ) %>%
  dplyr::mutate(
    Class = ifelse(ID == 131, "Barren (2)", Class),
    Class = ifelse(ID == 152, "Shrubland (2)", Class)
  ) %>%
  dplyr::filter(!is.na(Class))

tablesHeaders <-
  readr::read_rds("data-raw/tablesHeaders.rds")

nlcd <-
  readr::read_csv("data-raw/nlcd_rat.csv")

(httr::GET(
  "https://www.mrlc.gov/geoserver/ows",
  query = list(
    service = "WCS",
    version = "2.0.1",
    request = "GetCapabilities"
  )
) |>
  httr::content() |>
  xml2::as_list())$Capabilities$Contents |>
  purrr::map(\(x)
  list(
    Title = x$Title[[1]],
    id = x$CoverageId[[1]]
  )) |>
  magrittr::set_names(NULL) |>
  purrr::transpose() |>
  purrr::map(unlist) |>
  tibble::as_tibble() |>
  dplyr::filter(
    stringr::str_starts(id, "mrlc_download"),
    stringr::str_detect(id, "NLCD|nlcd"),
    stringr::str_detect(id, "Science_Product", negate = TRUE),
    stringr::str_detect(id, "Disturbance_Date", negate = TRUE),
    stringr::str_detect(id, "Fractional_Component", negate = TRUE),
    stringr::str_detect(id, "descriptor", negate = TRUE),
    stringr::str_detect(id, "Descriptor", negate = TRUE),
    stringr::str_detect(id, "Pixels", negate = TRUE),
    stringr::str_detect(id, "Annual_NLCD", negate = TRUE),
    stringr::str_detect(id, "Count", negate = TRUE),
    stringr::str_detect(id, "Index", negate = TRUE)
  ) |>
  dplyr::mutate(
    Title = stringr::str_replace(
      Title,
      "nlcd_tcc_conus",
      "NLCD"
    ),
    Title = stringr::str_replace(
      Title,
      "v2021-4",
      "Tree_Canopy_L48"
    ),
    Title = stringr::str_replace(
      Title,
      "Land_Cover",
      "landcover"
    ),
    Title = stringr::str_replace(
      Title,
      "Impervious",
      "impervious"
    ),
    Title = stringr::str_replace(
      Title,
      "Tree_Canopy",
      "canopy"
    )
  ) |>
  tidyr::separate_wider_delim(
    cols = Title,
    delim = "_",
    names = c(
      "NLCD",
      "year",
      "dataset",
      "landmass"
    )
  ) |>
  dplyr::select(!NLCD) |>
  dplyr::mutate(
    year = as.integer(year),
    dataset = factor(dataset,
      levels = c(
        "landcover",
        "impervious",
        "canopy"
      ),
      ordered = TRUE
    ),
    landmass = factor(landmass,
      levels = c(
        "L48",
        "AK",
        "HI",
        "PR"
      ),
      ordered = TRUE
    )
  ) |>
  dplyr::arrange(landmass, dataset, year) |>
  dplyr::mutate(
    available =
      purrr::map_lgl(
        id,
        \(x){
          out <-
            tryCatch(
              paste0(
                "WCS:",
                httr::modify_url(
                  "https://dmsdata.cr.usgs.gov/geoserver/ows",
                  query =
                    list(
                      version = "2.0.1",
                      coverageid = x
                    )
                )
              ) |>
                terra::rast(),
              error = function(e) FALSE
            )
          if (identical(out, FALSE)) {
            FALSE
          } else {
            TRUE
          }
        }
      )
  ) |>
  print(n = 1000)


bigger_dim <- "https://www.mrlc.gov/geoserver/rcmap_anhb/wcs?request=GetCoverage&service=WCS&version=2.0.1&coverageid=rcmap_anhb__rcmap_annual_herbaceous_2020&subset=X(-1054529.13659319,-1022600.05347719)&subset=Y(2020307.11783295,2056926.21900819)"

httr::GET(
  "https://www.mrlc.gov/geoserver/mrlc_download/NLCD_2001_Land_Cover_AK/wcs",
  query =
    list(
      service = "WCS",
      request = "DescribeCoverage",
      version = "2.0.1",
      coverageid = "NLCD_2001_Land_Cover_AK"
    )
) %>%
  httr::content() %>%
  xml2::as_list()

FedData::meve %>%
  sf::st_transform(5070) %>%
  sf::st_bbox()

test <-
  httr::modify_url(
    "https://www.mrlc.gov/geoserver/mrlc_download/NLCD_2001_Land_Cover_AK/wcs",
    query =
      list(
        service = "WCS",
        request = "GetCoverage",
        version = "2.0.1",
        coverageid = "NLCD_2001_Land_Cover_AK"
      )
  ) |>
  terra::rast()

paste0("WCS:", )

httr::modify_url(
  "https://www.mrlc.gov/geoserver/mrlc_download/NLCD_2001_Land_Cover_AK/wcs",
  query =
    list(
      version = "1.1.1",
      identifiers = "mrlc_download:NLCD_2001_Land_Cover_AK"
    )
) |>
  terra::rast(drivers = "WCS")


sf::st_crs(test)

##### National Park Spatial Polygon
meve <-
  FedData::get_padus(
    template = "Mesa Verde National Park",
    label = "meve",
    layer = "Proclamation_and_Other_Planning_Boundaries",
    force.redo = TRUE
  ) %$%
  Proclamation_and_Other_Planning_Boundaries %>%
  sf::st_geometry() %>%
  sf::st_transform(4326)



# A test dataset of the raw NLCD
# httr::GET(
#   url = "https://s3-us-west-2.amazonaws.com/mrlc/NLCD_2016_Land_Cover_L48_20190424.zip",
#   httr::write_disk(
#     path = paste0(tempdir(), "/NLCD_2016.zip"),
#     overwrite = TRUE
#   ),
#   httr::progress()
# )
#
# unzip(paste0(tempdir(), "/NLCD_2016.zip"),
#   exdir = tempdir()
# )
#
# nlcd <-
#   paste0(tempdir, "/NLCD_2016_Land_Cover_L48_20190424.img") %>%
#   terra::rast() %>%
#   terra::crop(
#     .,
#     sf::st_transform(
#       meve,
#       terra::crs(.)
#     )
#   )

## A 1x1 degree grid
# grid <- sp::GridTopology(
#   cellcentre.offset = c(-179.5, -89.5),
#   cellsize = c(1, 1),
#   cells.dim = c(360, 180)
# ) %>%
#   # sp::SpatialGrid(proj4string=CRS("+proj=longlat +datum=WGS84")) %>%
#   methods::as("SpatialPolygons") %>%
#   sf::st_as_sfc() %>%
#   sf::st_set_crs(4326)

usethis::use_data(tablesHeaders, nass, nlcd, grid, overwrite = TRUE, internal = TRUE)
usethis::use_data(meve, overwrite = TRUE)
