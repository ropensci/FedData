
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
[![R-CMD-check](https://github.com/ropensci/FedData/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/ropensci/FedData/actions/workflows/check-standard.yaml)
[![Codecov test
coverage](https://codecov.io/gh/ropensci/FedData/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ropensci/FedData?branch=master)
[![Zenodo
DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.596344.svg)](https://doi.org/10.5281/zenodo.596344)
[![ROpenSci
Status](https://badges.ropensci.org/13_status.svg)](https://github.com/ropensci/software-review/issues/13)

`FedData` is an *R* package implementing functions to automate
downloading geospatial data available from several federated data
sources.

Currently, the package enables extraction from nine datasets:

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
- The
  [PAD-US](https://www.usgs.gov/programs/gap-analysis-project/science/pad-us-data-overview)
  dataset of protected area boundaries from the USGS

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
  <https://r-spatial.github.io/sf/>, then install either from CRAN or
  GitHub.

### Getting Started

Check out our [Getting Started
article](https://docs.ropensci.org/FedData/articles/FedData.html).

### Acknowledgements

This package is a product of SKOPE ([Synthesizing Knowledge of Past
Environments](https://www.openskope.org/)) and the [Village Ecodynamics
Project](https://crowcanyon.github.io/veparchaeology/) through grants
awarded to the University of Montana, the [Crow Canyon Archaeological
Center](https://crowcanyon.org/), and Washington State University by the
National Science Foundation. This software is licensed under the [MIT
license](https://opensource.org/license/mit). Continuing development is
supported by the [Montana Climate Office](https://climate.umt.edu).

FedData was reviewed for [rOpenSci](https://ropensci.org) by
[@jooolia](https://github.com/jooolia), and was greatly improved as a
result. [rOpenSci](https://ropensci.org) on-boarding was coordinated by
[@sckott](https://github.com/sckott).

<!-- [![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org) -->
