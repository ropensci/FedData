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
  url <- paste0("https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU4/HighResolution/GDB/NHD_H_",area,"_HU4_GDB.zip")
  expect_false(suppressWarnings(httr::http_error(url)))
  
  area <- "blah"
  url <- paste0("https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU4/HighResolution/GDB/NHD_H_",area,"_HU4_GDB.zip")
  expect_true(suppressWarnings(httr::http_error(url)))
})
