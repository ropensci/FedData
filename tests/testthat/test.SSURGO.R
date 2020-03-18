library(FedData)
library(httr)
context("NRCS soils database (SSURGO) tests")

test_that("The SSURGO inventory dataset is available at the correct URL", {
  skip_on_cran()

  url <- "http://websoilsurvey.sc.egov.usda.gov/DataAvailability/SoilDataAvailabilityShapefile.zip"
  expect_false(suppressWarnings(httr::http_error(url)))

  url <- "http://websoilsurvey.sc.egov.usda.gov/DataAvailability/blah.zip"
  expect_true(suppressWarnings(httr::http_error(url)))
})

test_that("The SoilDB data queries work", {
  skip_on_cran()

  template <- "CO670"
  q <- paste0("SELECT areasymbol, saverest FROM sacatalog WHERE areasymbol IN (", paste(paste0("'", template, "'"), collapse = ","), ");")
  expect_is(FedData:::soils_query(q), "data.frame")

  template <- "blah"
  q <- paste0("SELECT areasymbol, saverest FROM sacatalog WHERE areasymbol IN (", paste(paste0("'", template, "'"), collapse = ","), ");")
  expect_error(FedData:::soils_query(q))
})

test_that("The SSURGO datasets are available at the correct URL", {
  skip_on_cran()

  template <- "CO670"
  q <- paste0("SELECT areasymbol, saverest FROM sacatalog WHERE areasymbol IN (", paste(paste0("'", template, "'"), collapse = ","), ");")
  q <- FedData:::soils_query(q)
  q$saverest <- as.Date(q$saverest, format = "%m/%d/%Y")

  url <- paste("https://websoilsurvey.sc.egov.usda.gov/DSD/Download/Cache/SSA/wss_SSA_", q$areasymbol, "_[", q$saverest, "].zip", sep = "")
  expect_error(suppressWarnings(curl::curl(url) %>% readLines(n = 1)), NA)

  url <- "http://websoilsurvey.sc.egov.usda.gov/DSD/Download/Cache/SSA/blah.zip"
  expect_error(suppressWarnings(curl::curl(url) %>% readLines(n = 1)))
})
