library(FedData)
library(httr)
context("Daymet tests")

test_that("The Daymet tiles are available at the correct URL", {
  year <- 1985
  element <- "prcp"
  tileID <- 11376
  url <- paste0('http://thredds.daac.ornl.gov/thredds/fileServer/ornldaac/1328/tiles/',year,'/',tileID,'_',year,'/',element,'.nc')
  expect_false(suppressWarnings(httr::http_error(url)))
  
  tileID <- 0
  url <- paste0('http://thredds.daac.ornl.gov/thredds/fileServer/ornldaac/1328/tiles/',year,'/',tileID,'_',year,'/',element,'.nc')
  expect_true(suppressWarnings(httr::http_error(url)))
})
