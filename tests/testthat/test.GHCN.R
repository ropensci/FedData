library(FedData)
library(httr)
context("Global Historical Climatology Network tests")

test_that("The GHCN inventory is available at the correct URL", {
  expect_false(suppressWarnings(httr::http_error("ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/ghcnd-inventory.txt")))
  expect_error(suppressWarnings(httr::http_error("ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/blah.txt")))
})

test_that("The GHCN daily weather stations are available at the correct URL", {
  expect_false(suppressWarnings(httr::http_error("ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/all/USC00051886.dly")))
  expect_error(suppressWarnings(httr::http_error("ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/all/xxxxxx.dly")))
})
