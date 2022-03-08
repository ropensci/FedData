library(FedData)
library(httr)
context("National Elevation Dataset tests")

test_that("The NED tiles are available at the correct URL", {
  skip_on_cran()

  res <- "1"
  tileNorthing <- "35"
  tileWesting <- "100"
  url <- paste0(
    "https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/", res, "/TIFF/current/n", tileNorthing, "w", tileWesting, "/USGS_", res, "_n", tileNorthing, "w", tileWesting,
    ".tif"
  )
  expect_false(suppressWarnings(httr::http_error(url)))

  res <- "13"
  url <- paste0(
    "https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/", res, "/TIFF/current/n", tileNorthing, "w", tileWesting, "/USGS_", res, "_n", tileNorthing, "w", tileWesting,
    ".tif"
  )
  expect_false(suppressWarnings(httr::http_error(url)))

  tileNorthing <- "350"
  url <- paste0(
    "https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/", res, "/TIFF/current/n", tileNorthing, "w", tileWesting, "/USGS_", res, "_n", tileNorthing, "w", tileWesting,
    ".tif"
  )
  expect_true(suppressWarnings(httr::http_error(url)))
})
