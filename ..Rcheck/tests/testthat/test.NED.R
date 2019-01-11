library(FedData)
library(httr)
context("National Elevation Dataset tests")

test_that("The NED tiles are available at the correct URL", {
  res <- "1"
  tileNorthing <- "35"
  tileWesting <- "100"
  url <- paste0("https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/", res, "/ArcGrid/USGS_NED_",res,"_n", tileNorthing, "w", tileWesting, 
                "_ArcGrid.zip")
  expect_false(suppressWarnings(httr::http_error(url)))
  
  res <- "13"
  url <- paste0("https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/", res, "/ArcGrid/USGS_NED_",res,"_n", tileNorthing, "w", tileWesting, 
                "_ArcGrid.zip")
  expect_false(suppressWarnings(httr::http_error(url)))
  
  tileNorthing <- "350"
  url <- paste0("https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/", res, "/ArcGrid/USGS_NED_",res,"_n", tileNorthing, "w", tileWesting, 
                "_ArcGrid.zip")
  expect_true(suppressWarnings(httr::http_error(url)))
})

test_that("NED URLs fall back to legacy format", {
  res <- "1"
  tileNorthing <- "36"
  tileWesting <- "110"
  
  url <- paste0("https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/", res, "/ArcGrid/USGS_NED_",res,"_n", tileNorthing, "w", tileWesting, 
                "_ArcGrid.zip")
  
  # USGS is in the process of updating all of their URLs for the NED.
  # If the new URL doesn't exist, try the old one.
  if(httr::http_error(url))
    url <- paste0("https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/", res, "/ArcGrid/n", tileNorthing, "w", tileWesting, 
                  ".zip")
  
  expect_false(suppressWarnings(httr::http_error(url)))
  
  res <- "13"
  expect_false(suppressWarnings(httr::http_error(url)))
  
})

