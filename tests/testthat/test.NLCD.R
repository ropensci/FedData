library(FedData)
library(httr)
library(magrittr)
context("National Land Cover Dataset tests")

test_that("The NLCD web coverage service is available at the correct URL", {
  skip_on_cran()

  year <- 2016
  dataset <- "Land_Cover"
  landmass <- "L48"

  coverage <- paste0("NLCD_", year, "_", dataset, "_", landmass)
  source <- paste0("https://www.mrlc.gov/geoserver/mrlc_download/", coverage, "/wcs")

  source %>%
    httr::modify_url(
      query = list(
        service = "WCS",
        version = "2.0.1",
        request = "DescribeCoverage",
        coverageid = coverage
      )
    )

  cat(
    source %>%
      httr::modify_url(
        query = list(
          service = "WCS",
          version = "2.0.1",
          request = "DescribeCoverage",
          coverageid = coverage
        )
      ) %>%
      httr::GET() %>%
      httr::status_code(),
    "\n"
  )

  expect_true(
    source %>%
      httr::modify_url(
        query = list(
          service = "WCS",
          version = "2.0.1",
          request = "DescribeCoverage",
          coverageid = coverage
        )
      ) %>%
      httr::GET() %>%
      httr::status_code() %>%
      identical(200L)
  )
})

test_that(
  "The NLCD provides the same data as a raw download",
  {
    skip_on_cran()

    ## Raw NLCD (Puerto Rico, 2001)
    raw_zip <- tempfile(fileext = ".zip")
    raw_tmp <- tempfile()
    "https://s3-us-west-2.amazonaws.com/mrlc/PR_landcover_wimperv_10-28-08_se5.zip" %>%
      httr::GET(httr::progress(), httr::write_disk(raw_zip, overwrite = TRUE))
    unzip(raw_zip, exdir = raw_tmp)
    raw <- terra::rast(paste0(raw_tmp, "/pr_landcover_wimperv_10-28-08_se5.img"))
    feddata <-
      FedData::get_nlcd(
        template = raw,
        label = "PR",
        year = 2001,
        landmass = "PR",
        force.redo = TRUE
      )

    raw %<>%
      raster::raster()

    expect_true(
      all((feddata - raw)[] == 0, na.rm = TRUE)
    )
  }
)
