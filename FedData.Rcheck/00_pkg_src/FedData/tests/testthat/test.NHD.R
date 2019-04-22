library(FedData)
library(httr)
context("National Hydrography Dataset tests")

test_that("The NHD Subregion boundaries dataset is availableat the correct URL", {
  skip_on_cran()
  
  url <- "https://github.com/ropensci/FedData/raw/master/data-raw/nhd_huc4.gpkg.zip"
  expect_false(suppressWarnings(httr::http_error(url)))
  
  url <- "https://github.com/ropensci/FedData/raw/master/data-raw/blah.zip"
  expect_true(suppressWarnings(httr::http_error(url)))
})

test_that("The NHD staged subregions are available at the correct URL", {
  skip_on_cran()
  
  area <- 1404
  url <- paste0("https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU4/HighResolution/GDB/NHD_H_",area,"_HU4_GDB.zip")
  expect_false(suppressWarnings(httr::http_error(url)))
  
  area <- "blah"
  url <- paste0("https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU4/HighResolution/GDB/NHD_H_",area,"_HU4_GDB.zip")
  expect_true(suppressWarnings(httr::http_error(url)))
})
