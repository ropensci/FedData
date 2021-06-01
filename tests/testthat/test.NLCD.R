library(FedData)
library(httr)
context("National Land Cover Dataset tests")

test_that("The NLCD web coverage service is available at the correct URL", {
  skip_on_cran()

  year <- 2016
  dataset <- "Land_Cover"
  landmass <- "L48"

  coverage <- paste0("NLCD_", year, "_", dataset, "_", landmass)
  source <- "https://www.mrlc.gov/geoserver/wcs"

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
  "The NLCD web coverage service provides the same data as a raw download",
  {
    skip_on_cran()

    year <- 2016
    dataset <- "Land_Cover"
    landmass <- "L48"

    coverage <- paste0("NLCD_", year, "_", dataset, "_", landmass)
    source <- "https://www.mrlc.gov/geoserver/wcs"

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
  }
)
