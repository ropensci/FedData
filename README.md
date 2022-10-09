
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN
version](https://www.r-pkg.org/badges/version/FedData)](https://cran.r-project.org/package=FedData)
[![CRAN downloads per
month](https://cranlogs.r-pkg.org/badges/FedData)](https://github.com/r-hub/cranlogs.app)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/grand-total/FedData)](https://github.com/r-hub/cranlogs.app)
[![R-CMD-check](https://github.com/ropensci/FedData/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ropensci/FedData/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/ropensci/FedData/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ropensci/FedData?branch=master)
[![Zenodo
DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.596344.svg)](https://doi.org/10.5281/zenodo.596344)
[![ROpenSci
Status](https://badges.ropensci.org/13_status.svg)](https://github.com/ropensci/software-review/issues/13)

**FedData version 3.0 has been released to CRAN! There are several
breaking changes in the FedData API from version 2.x. Please see
\[NEWS.md\] for a list of changes.**

`FedData` is an *R* package implementing functions to automate
downloading geospatial data available from several federated data
sources.

Currently, the package enables extraction from seven datasets:

- The [National Elevation Dataset (NED)](https://ned.usgs.gov) digital
  elevation models (1 and 1/3 arc-second; USGS)
- The [National Hydrography Dataset
  (NHD)](https://www.usgs.gov/national-hydrography/national-hydrography-dataset)
  (USGS)
- The [Soil Survey Geographic (SSURGO)
  database](https://websoilsurvey.sc.egov.usda.gov/) from the National
  Cooperative Soil Survey (NCSS), which is led by the Natural Resources
  Conservation Service (NRCS) under the USDA
- The [Global Historical Climatology Network
  (GHCN)](https://www.ncei.noaa.gov/products/land-based-station/global-historical-climatology-network-daily),
  coordinated by National Climatic Data Center at NOAA
- The [Daymet](https://daymet.ornl.gov/) gridded estimates of daily
  weather parameters for North America, version 4, available from the
  Oak Ridge National Laboratory’s Distributed Active Archive Center
  (DAAC)
- The [International Tree Ring Data Bank
  (ITRDB)](https://www.ncei.noaa.gov/products/paleoclimatology/tree-ring),
  coordinated by National Climatic Data Center at NOAA
- The [National Land Cover Database (NLCD)](https://www.mrlc.gov/)
- The [NASS Cropland Data
  Layer](https://www.nass.usda.gov/Research_and_Science/Cropland/SARS1a.php)
  from the National Agricultural Statistics Service

This package is designed with the large-scale geographic information
system (GIS) use-case in mind: cases where the use of dynamic
web-services is impractical due to the scale (spatial and/or temporal)
of analysis. It functions primarily as a means of downloading tiled or
otherwise spatially-defined datasets; additionally, it can preprocess
those datasets by extracting data within an area of interest (AoI),
defined spatially. It relies heavily on the
[**sf**](https://cran.r-project.org/package=sf) and
[**terra**](https://cran.r-project.org/package=terra) packages.

### Development

- [Kyle Bocinsky](https://www.bocinsky.io) - Montana Climate Office,
  Missoula, MT

### Contributors

- Dylan Beaudette - USDA-NRCS Soil Survey Office, Sonora, CA
- Jeffrey Hollister - US EPA Atlantic Ecology Division, Narragansett, RI
- Scott Chamberlain - ROpenSci and Museum of Paleontology at UC Berkeley

### Install `FedData`

- From CRAN:

``` r
install.packages("FedData")
```

- Development version from GitHub:

``` r
install.packages("devtools")
devtools::install_github("ropensci/FedData")
```

- Linux: Follow instructions for installing `sf` available at
  <https://r-spatial.github.io/sf/>.

### Demonstration

This demonstration script is available as an R Markdown document in the
GitHub repository: <https://github.com/ropensci/FedData>.

#### Load `FedData` and define a study area

``` r
# FedData Tester
library(FedData)
library(magrittr)

# FedData comes loaded with the boundary of Mesa Verde National Park, for testing
FedData::meve
```

#### Get and plot the National Elevation Dataset for the study area

``` r
# Get the NED (USA ONLY)
# Returns a raster
NED <- get_ned(
  template = FedData::meve,
  label = "meve"
)
# Plot with raster::plot
raster::plot(NED)
```

<img src="man/figures/README-NED-1.png" width="100%" />

#### Get and plot the Daymet dataset for the study area

``` r
# Get the DAYMET (North America only)
# Returns a raster
DAYMET <- get_daymet(
  template = FedData::meve,
  label = "meve",
  elements = c("prcp", "tmax"),
  years = 1980:1985
)
# Plot with raster::plot
raster::plot(DAYMET$tmax$X1985.10.23)
```

<img src="man/figures/README-DAYMET-1.png" width="100%" />

#### Get and plot the daily GHCN precipitation data for the study area

``` r
# Get the daily GHCN data (GLOBAL)
# Returns a list: the first element is the spatial locations of stations,
# and the second is a list of the stations and their daily data
GHCN.prcp <- get_ghcn_daily(
  template = FedData::meve,
  label = "meve",
  elements = c("prcp")
)
#> Warning: `select_()` was deprecated in dplyr 0.7.0.
#> ℹ Please use `select()` instead.
#> ℹ The deprecated feature was likely used in the dplyr package.
#>   Please report the issue at <]8;;https://github.com/tidyverse/dplyr/issueshttps://github.com/tidyverse/dplyr/issues]8;;>.
#> Warning in CPL_write_ogr(obj, dsn, layer, driver,
#> as.character(dataset_options), : GDAL Error 1: /private/var/folders/ys/
#> 7l0z3wlx7z14qxn9v0m9ckhw0000gq/T/Rtmp5iX7lU/FedData/extractions/ghcn/meve/
#> meve_GHCN_stations.shp does not appear to be a file or directory.
#> Warning: `filter_()` was deprecated in dplyr 0.7.0.
#> ℹ Please use `filter()` instead.
#> ℹ See vignette('programming') for more help
#> ℹ The deprecated feature was likely used in the dplyr package.
#>   Please report the issue at <]8;;https://github.com/tidyverse/dplyr/issueshttps://github.com/tidyverse/dplyr/issues]8;;>.
#> Warning: `funs_()` was deprecated in dplyr 0.7.0.
#> ℹ Please use `funs()` instead.
#> ℹ See vignette('programming') for more help
#> ℹ The deprecated feature was likely used in the dplyr package.
#>   Please report the issue at <]8;;https://github.com/tidyverse/dplyr/issueshttps://github.com/tidyverse/dplyr/issues]8;;>.
#> Warning: `funs()` was deprecated in dplyr 0.8.0.
#> ℹ Please use a list of either functions or lambdas:
#> 
#> # Simple named list: list(mean = mean, median = median)
#> 
#> # Auto named with `tibble::lst()`: tibble::lst(mean, median)
#> 
#> # Using lambdas list(~ mean(., trim = .2), ~ median(., na.rm = TRUE))
#> ℹ The deprecated feature was likely used in the dplyr package.
#>   Please report the issue at <]8;;https://github.com/tidyverse/dplyr/issueshttps://github.com/tidyverse/dplyr/issues]8;;>.
# Plot the NED again
raster::plot(NED)
# Plot the spatial locations
sp::plot(GHCN.prcp$spatial,
  pch = 1,
  add = TRUE
)
legend("bottomleft",
  pch = 1,
  legend = "GHCN Precipitation Records"
)
```

<img src="man/figures/README-GHCN-precipitation-1.png" width="100%" />

#### Get and plot the daily GHCN temperature data for the study area

``` r
# Elements for which you require the same data
# (i.e., minimum and maximum temperature for the same days)
# can be standardized using standardize==T
GHCN.temp <- get_ghcn_daily(
  template = FedData::meve,
  label = "meve",
  elements = c("tmin", "tmax"),
  years = 1980:1985,
  standardize = TRUE
)
#> Warning: `arrange_()` was deprecated in dplyr 0.7.0.
#> ℹ Please use `arrange()` instead.
#> ℹ See vignette('programming') for more help
#> ℹ The deprecated feature was likely used in the dplyr package.
#>   Please report the issue at <]8;;https://github.com/tidyverse/dplyr/issueshttps://github.com/tidyverse/dplyr/issues]8;;>.
# Plot the NED again
raster::plot(NED)
# Plot the spatial locations
sp::plot(GHCN.temp$spatial,
  add = TRUE,
  pch = 1
)
legend("bottomleft",
  pch = 1,
  legend = "GHCN Temperature Records"
)
```

<img src="man/figures/README-GHCN-temperature-1.png" width="100%" />

#### Get and plot the National Hydrography Dataset for the study area

``` r
# Get the NHD (USA ONLY)
get_nhd(
  template = FedData::meve,
  label = "meve"
) %>%
  plot_nhd(template = FedData::meve)
```

<img src="man/figures/README-NHD-1.png" width="100%" />

#### Get and plot the NRCS SSURGO data for the study area

``` r
# Get the NRCS SSURGO data (USA ONLY)
SSURGO.MEVE <- get_ssurgo(
  template = FedData::meve,
  label = "meve"
)
# Plot the NED again
raster::plot(NED)
# Plot the SSURGO mapunit polygons
plot(SSURGO.MEVE$spatial$geom,
  lwd = 0.1,
  add = TRUE
)
```

<img src="man/figures/README-SSURGO-1.png" width="100%" />

#### Get and plot the NRCS SSURGO data for particular soil survey areas

``` r
# Or, download by Soil Survey Area names
SSURGO.areas <- get_ssurgo(
  template = c("CO670", "CO075"),
  label = "CO_TEST"
)

# Let's just look at spatial data for CO675
SSURGO.areas.CO675 <-
  SSURGO.areas$spatial %>%
  dplyr::filter(AREASYMBOL == "CO075")

# And get the NED data under them for pretty plotting
NED.CO675 <- get_ned(
  template = SSURGO.areas.CO675,
  label = "SSURGO_CO675"
)

# Plot the SSURGO mapunit polygons, but only for CO675
raster::plot(NED.CO675)
plot(SSURGO.areas.CO675$geom,
  lwd = 0.1,
  add = TRUE
)
```

<img src="man/figures/README-SSURGO-area-1.png" width="100%" />

#### Get and plot the ITRDB chronology locations in the study area

``` r
# Get the ITRDB records
# Buffer MEVE, because there aren't any chronologies in the Park
ITRDB <- get_itrdb(
  template = FedData::meve %>%
    sf::st_buffer(50000),
  label = "meve",
  measurement.type = "Ring Width",
  chronology.type = "Standard"
)
#> Warning in eval(jsub, SDenv, parent.frame()): NAs introduced by coercion
#> Warning: attribute variables are assumed to be spatially constant throughout all
#> geometries

# Plot the MEVE buffer
plot(
  FedData::meve %>%
    sf::st_buffer(50000) %>%
    sf::st_transform(4326)
)
# Map the locations of the tree ring chronologies
plot(ITRDB$metadata$geometry,
  pch = 1,
  add = TRUE
)
legend("bottomleft",
  pch = 1,
  legend = "ITRDB chronologies"
)
```

<img src="man/figures/README-ITRDB-1.png" width="100%" />

#### Get and plot the National Land Cover Dataset for the study area

``` r
# Get the NLCD (USA ONLY)
# Returns a raster
NLCD <- get_nlcd(
  template = FedData::meve,
  year = 2011,
  label = "meve"
)

# Plot with raster::plot
raster::plot(NLCD)
```

<img src="man/figures/README-NLCD-1.png" width="100%" />

#### Get and plot the NASS Cropland Data Layer for the study area

``` r
# Get the NASS (USA ONLY)
# Returns a raster
NASS_CDL <- get_nass_cdl(
  template = FedData::meve,
  year = 2016,
  label = "meve"
)
# Plot with raster::plot
raster::plot(NASS_CDL)
```

<img src="man/figures/README-NASS-CDL-1.png" width="100%" />

``` r

# Get the NASS CDL classification table
raster::levels(NASS_CDL)[[1]]

# Also, a convenience function loading the NASS CDL categories and hex colors
cdl_colors()
```

------------------------------------------------------------------------

### Acknowledgements

This package is a product of SKOPE ([Synthesizing Knowledge of Past
Environments](https://www.openskope.org/)) and the [Village Ecodynamics
Project](https://crowcanyon.github.io/veparchaeology/) through grants
awarded to the [Crow Canyon Archaeological
Center](https://www.crowcanyon.org) and Washington State University by
the National Science Foundation. This software is licensed under the
[MIT license](https://opensource.org/licenses/MIT). Continuing
development is supported by the [Montana Climate
Office](https://climate.umt.edu).

FedData was reviewed for [rOpenSci](https://ropensci.org) by
[@jooolia](https://github.com/jooolia), and was greatly improved as a
result. [rOpenSci](https://ropensci.org) on-boarding was coordinated by
[@sckott](https://github.com/sckott).

<!-- [![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org) -->
