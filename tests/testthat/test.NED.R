library(FedData)
library(httr)
context("National Elevation Dataset tests")

test_that("The NED tiles are available at the correct URL", {
  res <- "1"
  tileNorthing <- "35"
  tileWesting <- "100"
  url <- paste('ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/NED/',res,'/ArcGrid/n',tileNorthing,'w',tileWesting,'.zip',sep='')
  expect_false(suppressWarnings(httr::http_error(url)))
  
  res <- "13"
  url <- paste('ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/NED/',res,'/ArcGrid/n',tileNorthing,'w',tileWesting,'.zip',sep='')
  expect_false(suppressWarnings(httr::http_error(url)))
  
  tileNorthing <- "350"
  url <- paste('ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/NED/',res,'/ArcGrid/n',tileNorthing,'w',tileWesting,'.zip',sep='')
  expect_error(suppressWarnings(httr::http_error(url)))
})
