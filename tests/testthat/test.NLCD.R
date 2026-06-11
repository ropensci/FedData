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
    "https://www.mrlc.gov/downloads/sciweb1/shared/mrlc/data-bundles/PR_landcover_wimperv_10-28-08_se5.zip" %>%
      httr::GET(httr::progress(), httr::write_disk(raw_zip, overwrite = TRUE))
    unzip(raw_zip, exdir = raw_tmp)
    raw <- terra::rast(file.path(raw_tmp, "pr_landcover_wimperv_10-28-08_se5.img"))
    feddata <-
      FedData::get_nlcd(
        template = raw,
        label = basename(tempfile()),
        year = 2001,
        landmass = "PR"
      )

    expect_true(
      all((feddata - raw)[] == 0, na.rm = TRUE)
    )
  }
)

test_that("The Annual NLCD web coverage services are available", {
  skip_on_cran()

  for (workspace in FedData:::nlcd_annual_wcs_coverages) {
    expect_identical(
      paste0("https://dmsdata.cr.usgs.gov/geoserver/", workspace, "/wcs") %>%
        httr::modify_url(
          query = list(
            service = "WCS",
            version = "2.0.1",
            request = "DescribeCoverage",
            coverageId = paste0(workspace, ":", sub("^mrlc_", "", workspace))
          )
        ) %>%
        httr::GET() %>%
        httr::status_code(),
      200L,
      info = workspace
    )
  }
})

test_that("get_nlcd_annual returns rasters on the native grid", {
  skip_on_cran()

  out <-
    get_nlcd_annual(
      template = FedData::meve,
      label = basename(tempfile()),
      year = 2024,
      product = "LndCov"
    )

  expect_s3_class(out, "tbl_df")
  expect_true(inherits(out$rast[[1]], "SpatRaster"))
  expect_identical(
    sf::st_crs(terra::crs(out$rast[[1]]))$epsg,
    5070L
  )
  expect_identical(unique(terra::res(out$rast[[1]])), 30)
})
