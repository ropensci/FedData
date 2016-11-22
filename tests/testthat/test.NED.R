library(FedData)
library(httr)
context("National Elevation Dataset tests")

test_that("The NED tiles are available at the correct URL", {
  res <- "1"
  tileNorthing <- "35"
  tileWesting <- "100"
  url <- paste0("https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/",res,"/ArcGrid/n",tileNorthing,"w",tileWesting,".zip")
  expect_false(suppressWarnings(httr::http_error(url)))
  
  res <- "13"
  url <- paste0("https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/",res,"/ArcGrid/n",tileNorthing,"w",tileWesting,".zip")
  expect_false(suppressWarnings(httr::http_error(url)))
  
  tileNorthing <- "350"
  url <- paste0("https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/",res,"/ArcGrid/n",tileNorthing,"w",tileWesting,".zip")
  expect_true(suppressWarnings(httr::http_error(url)))
})
