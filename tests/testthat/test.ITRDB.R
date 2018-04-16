library(FedData)
library(httr)
library(curl)
context("International Tree Ring Data Bank tests")

# test_that("The ITRDB is available at the correct URL", {
#   expect_false(suppressWarnings(httr::http_error("ftp://ftp.ncdc.noaa.gov/pub/data/paleo/treering/chronologies/")))
#   expect_error(suppressWarnings(httr::http_error("ftp://ftp.ncdc.noaa.gov/pub/data/paleo/treering/blah/")))
# })

test_that("ITRDB version files are available", {
  opts <- list(
    verbose = F, 
    noprogress = T,
    fresh_connect = T, 
    ftp_use_epsv = T, 
    forbid_reuse = T,
    dirlistonly = T)
  hand <- curl::new_handle()
  curl::handle_setopt(hand, .list = opts)
  
  url <- 'ftp://ftp.ncdc.noaa.gov/pub/data/paleo/treering/chronologies/'
  filenames = readLines(curl::curl(url, handle = hand))
  filenames = paste(url, filenames, sep = "")
  filenames <- filenames[grep("*.zip",filenames)]
  
  expect_true(length(filenames) > 0)
})
