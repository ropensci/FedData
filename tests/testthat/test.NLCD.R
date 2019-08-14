library(FedData)
library(httr)
context("National Land Cover Dataset tests")

test_that("The NLCD web coverage service is available at the correct URL", {
  skip_on_cran()
  
  year = 2016
  dataset = "Land_Cover"
  landmass = "L48"
  
  coverage <- paste0("NLCD_",year,"_",dataset,"_",landmass)
  source <- paste0("https://www.mrlc.gov/geoserver/mrlc_display/",coverage,"/ows")
  
  cat(source %>% 
        httr::GET() %>% 
        httr::status_code(), "\n")
  
  expect_true(
    source %>% 
      httr::GET() %>% 
      httr::status_code() %>%
      identical(200L)
  )
  
})
