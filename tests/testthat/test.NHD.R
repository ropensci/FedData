library(FedData)
library(httr)
context("National Hydrography Dataset tests")

test_that("The NHD Subregion boundaries dataset is availableat the correct URL", {
  url <- 'ftp://ftp.igsb.uiowa.edu/gis_library/USA/huc_04.zip'
  expect_false(suppressWarnings(httr::http_error(url)))
  
  url <- 'ftp://ftp.igsb.uiowa.edu/gis_library/USA/blah.zip'
  expect_error(suppressWarnings(httr::http_error(url)))
})

test_that("The NHD staged subregions are available at the correct URL", {
  area <- 1404
  url <- paste('ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/Hydrography/NHD/HU4/HighResolution/GDB/NHD_H_',area,'_GDB.zip',sep='')
  expect_false(suppressWarnings(httr::http_error(url)))
  
  area <- "blah"
  url <- paste('ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/Hydrography/NHD/HU4/HighResolution/GDB/NHD_H_',area,'_GDB.zip',sep='')
  expect_error(suppressWarnings(httr::http_error(url)))
})
