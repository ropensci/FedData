library(tidyverse)
library(magrittr)
library(rmapshaper)
library(sf)

# Download HUC4 dataset from https://nrcs.app.box.com/v/gateway/file/233091823212
unzip("./data-raw/wbdhu4_a_us_september2017.zip",
      exdir = "./data-raw/wbdhu4_a_us_september2017")

nhd_huc4 <- sf::st_read("./data-raw/wbdhu4_a_us_september2017/wbdhu4_a_us_september2017.gdb",
                        layer = "WBDHU4") %>%
  dplyr::filter(!grepl("AK",STATES),
                !grepl("HI",STATES)) %>%
  dplyr::select(dplyr::starts_with("HUC")) %>%
  dplyr::rename(geometry = Shape) %>%
  tibble::as_tibble() %>%
  sf::st_as_sf() %>%
  sf::st_cast('MULTIPOLYGON')

sf::st_write(nhd_huc4,
             "./data-raw/nhd_huc4.gpkg",
             delete_dsn = TRUE)

zip("./data-raw/nhd_huc4.gpkg.zip","./data-raw/nhd_huc4.gpkg")

# devtools::use_data(nhd_huc4,
#                    overwrite = T)

unlink("./data-raw/wbdhu4_a_us_september2017.zip")

unlink("./data-raw/wbdhu4_a_us_september2017",
       recursive = TRUE)

unlink("./data-raw/nhd_huc4.gpkg",
       recursive = TRUE)
