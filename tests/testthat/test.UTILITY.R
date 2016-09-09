library(FedData)
library(raster)
context("Utility tests")

test_that("substr_right returns correct strings", {
  expect_equal(substr_right("This shit is bananas",7), "bananas")
  expect_equal(substr_right("seize the day", 7), "the day")
  expect_equal(substr_right("reunion", 5), "union")
})