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
