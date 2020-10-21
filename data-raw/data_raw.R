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
# #                    overwrite = T)
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
  dplyr::select(
    ID = VALUE,
    `Land Cover` = CLASS_NAME,
    Color = rgb
  )

tablesHeaders <-
  readr::read_rds("data-raw/tablesHeaders.rds")



##### National Park Spatial Polygon

# mvnp <- sf::read_sf("https://gist.githubusercontent.com/bocinsky/2081c48598de66c3b3dd377a5bb4519d/raw/707483e40e4a4f72fd3797096af1c5d2a4dfd0bb/meve.geojson")

glac <-
  "~/Downloads/NPS_-_Land_Resources_Division_Boundary_and_Tract_Data_Service-shp/" %>%
  sf::read_sf() %>%
  dplyr::filter(UNIT_CODE == "GLAC") %>%
  dplyr::select(GNIS_ID, UNIT_CODE, UNIT_NAME, UNIT_TYPE, REGION) %>%
  sf::st_transform(4326)



usethis::use_data(tablesHeaders, nass, overwrite = TRUE, internal = TRUE)
# usethis::use_data(mvnp, overwrite = TRUE)
usethis::use_data(glac, overwrite = TRUE)
