---
name: FedData-release
layout: post
title: FedData - Getting assorted geospatial data into R
date: 2017-11-06
authors:
  - name: Kyle Bocinsky
categories:
  - technotes
tags:
  - R
  - package
  - data-access
  - packages
  - spatial
  - geospatial
  - review
  - onboarding
  - FedData
---

The package [FedData](https://github.com/ropensci/FedData) is now part of [rOpenSci](https://ropensci.org/). FedData includes functions to automate downloading geospatial data available from several federated data sources (mainly sources maintained by the US Federal government).

Currently, the package enables extraction from six datasets:

-   The [National Elevation Dataset (NED)](http://ned.usgs.gov) digital elevation models (1 and 1/3 arc-second; USGS)
-   The [National Hydrography Dataset (NHD)](http://nhd.usgs.gov) (USGS)
-   The [Soil Survey Geographic (SSURGO) database](http://websoilsurvey.sc.egov.usda.gov/) from the National Cooperative Soil Survey (NCSS), which is led by the Natural Resources Conservation Service (NRCS) under the USDA,
-   The [Global Historical Climatology Network (GHCN)](http://www.ncdc.noaa.gov/data-access/land-based-station-data/land-based-datasets/global-historical-climatology-network-ghcn), coordinated by National Climatic Data Center at NOAA,
-   The [Daymet](https://daymet.ornl.gov/) gridded estimates of daily weather parameters for North America, version 3, available from the Oak Ridge National Laboratory's Distributed Active Archive Center (DAAC), and
-   The [International Tree Ring Data Bank (ITRDB)](http://www.ncdc.noaa.gov/data-access/paleoclimatology-data/datasets/tree-ring), coordinated by National Climatic Data Center at NOAA.

FedData is designed with the large-scale geographic information system (GIS) use-case in mind: cases where the use of dynamic web-services is impractical due to the scale (spatial and/or temporal) of analysis. It functions primarily as a means of downloading tiled or otherwise spatially-defined datasets; additionally, it can preprocess those datasets by extracting data within an area of interest (AoI), defined spatially. It relies heavily on the [**sp**](https://cran.r-project.org/package=sp), [**raster**](https://cran.r-project.org/package=raster), and [**rgdal**](https://cran.r-project.org/package=rgdal) packages.

Acknowledgements
----------------

FedData is a product of SKOPE ([Synthesizing Knowledge of Past Environments](http://www.openskope.org)) and the [Village Ecodynamics Project](http://veparchaeology.org/).

FedData was reviewed for [rOpenSci](https://ropensci.org) by [@jooolia](https://github.com/jooolia), and was greatly improved as a result. [rOpenSci](https://ropensci.org) onboarding was coordinated by [@sckott](https://github.com/sckott).

TODO
----

**The current CRAN version of FedData, v2.4, will be the final minor CRAN release of FedData 2. FedData 3 will be released in the coming months, but some code built on FedData 2 will not be compatible with FedData 3.**

FedData was initially developed prior to widespread use of modern web mapping services and RESTful APIs by many Federal data-holders. Future releases of FedData will limit data transfer by utilizing server-side geospatial and data queries. We will also implement of data grammars from [dplyr](https://github.com/hadley/dplyr), tidy data structures, piping throughout ([magrittr](https://github.com/tidyverse/magrittr)), functional programming using [purrr](https://github.com/hadley/purrr), simple features for spatial data from [sf](https://github.com/edzer/sfr), and local data storage in OGC-compliant data formats (probably geojson and netCDF). I am also aiming for 100% testing coverage!

All that being said, much of the functionality of the FedData package could be spun off into more domain-specific packages. For example, ITRDB download functions could be part of the [dplR](https://r-forge.r-project.org/projects/dplr/) dendrochronology package; concepts/functions having to do with the GHCN data integrated into [rnoaa](https://github.com/ropensci/rnoaa); and Daymet concepts integrated into [daymetr](https://github.com/khufkens/daymetr). I welcome any and all suggestions about how to improve the utility of FedData; please [submit an issue](https://github.com/ropensci/FedData/issues).

Examples
--------

### Load `FedData` and define a study area

``` r
# FedData Tester
library(FedData)
library(magrittr)

# Extract data for the Village Ecodynamics Project "VEPIIN" study area:
# http://veparchaeology.org
vepPolygon <- polygon_from_extent(raster::extent(672800, 740000, 4102000, 4170000),
                                  proj4string = "+proj=utm +datum=NAD83 +zone=12")
```

### Get and plot the National Elevation Dataset for the study area

``` r
# Get the NED (USA ONLY)
# Returns a raster
NED <- get_ned(template = vepPolygon,
               label = "VEPIIN")
# Plot with raster::plot
raster::plot(NED)
```

![](https://github.com/ropensci/FedData/raw/master/inst/image/README-unnamed-chunk-6-1.png)

### Get and plot the Daymet dataset for the study area

``` r
# Get the DAYMET (North America only)
# Returns a raster
DAYMET <- get_daymet(template = vepPolygon,
               label = "VEPIIN",
               elements = c("prcp","tmax"),
               years = 1980:1985)
# Plot with raster::plot
raster::plot(DAYMET$tmax$X1985.10.23)
```

![](https://github.com/ropensci/FedData/raw/master/inst/image/README-unnamed-chunk-7-1.png)

### Get and plot the daily GHCN precipitation data for the study area

``` r
# Get the daily GHCN data (GLOBAL)
# Returns a list: the first element is the spatial locations of stations,
# and the second is a list of the stations and their daily data
GHCN.prcp <- get_ghcn_daily(template = vepPolygon, 
                            label = "VEPIIN", 
                            elements = c('prcp'))
# Plot the NED again
raster::plot(NED)
# Plot the spatial locations
sp::plot(GHCN.prcp$spatial,
         pch = 1,
         add = TRUE)
legend('bottomleft',
       pch = 1,
       legend="GHCN Precipitation Records")
```

![](https://github.com/ropensci/FedData/raw/master/inst/image/README-unnamed-chunk-8-1.png)

### Get and plot the daily GHCN temperature data for the study area

``` r
# Elements for which you require the same data
# (i.e., minimum and maximum temperature for the same days)
# can be standardized using standardize==T
GHCN.temp <- get_ghcn_daily(template = vepPolygon, 
                            label = "VEPIIN", 
                            elements = c('tmin','tmax'), 
                            years = 1980:1985,
                            standardize = TRUE)
# Plot the NED again
raster::plot(NED)
# Plot the spatial locations
sp::plot(GHCN.temp$spatial,
         add = TRUE,
         pch = 1)
legend('bottomleft',
       pch = 1,
       legend = "GHCN Temperature Records")
```

![](https://github.com/ropensci/FedData/raw/master/inst/image/README-unnamed-chunk-9-1.png)

### Get and plot the National Hydrography Dataset for the study area

``` r
# Get the NHD (USA ONLY)
NHD <- get_nhd(template = vepPolygon, 
               label = "VEPIIN")
# Plot the NED again
raster::plot(NED)
# Plot the NHD data
NHD %>%
  lapply(sp::plot,
         col = 'black',
         add = TRUE)
```

![](https://github.com/ropensci/FedData/raw/master/inst/image/README-unnamed-chunk-10-1.png)

### Get and plot the NRCS SSURGO data for the study area

``` r
# Get the NRCS SSURGO data (USA ONLY)
SSURGO.VEPIIN <- get_ssurgo(template = vepPolygon, 
                     label = "VEPIIN")
#> Warning: 1 parsing failure.
#> row # A tibble: 1 x 5 col     row     col               expected actual expected   <int>   <chr>                  <chr>  <chr> actual 1  1276 slope.r no trailing characters     .5 file # ... with 1 more variables: file <chr>
# Plot the NED again
raster::plot(NED)
# Plot the SSURGO mapunit polygons
plot(SSURGO.VEPIIN$spatial,
     lwd = 0.1,
     add = TRUE)
```

![](https://github.com/ropensci/FedData/raw/master/inst/image/README-unnamed-chunk-11-1.png)

### Get and plot the NRCS SSURGO data for particular soil survey areas

``` r
# Or, download by Soil Survey Area names
SSURGO.areas <- get_ssurgo(template = c("CO670","CO075"), 
                           label = "CO_TEST")

# Let's just look at spatial data for CO675
SSURGO.areas.CO675 <- SSURGO.areas$spatial[SSURGO.areas$spatial$AREASYMBOL=="CO075",]

# And get the NED data under them for pretty plotting
NED.CO675 <- get_ned(template = SSURGO.areas.CO675,
                            label = "SSURGO_CO675")
               
# Plot the SSURGO mapunit polygons, but only for CO675
plot(NED.CO675)
plot(SSURGO.areas.CO675,
     lwd = 0.1,
     add = TRUE)
```

![](https://github.com/ropensci/FedData/raw/master/inst/image/README-unnamed-chunk-12-1.png)

### Get and plot the ITRDB chronology locations in the study area

``` r
# Get the ITRDB records
ITRDB <- get_itrdb(template = vepPolygon,
                        label = "VEPIIN",
                        makeSpatial = TRUE)
# Plot the NED again
raster::plot(NED)
# Map the locations of the tree ring chronologies
plot(ITRDB$metadata,
     pch = 1,
     add = TRUE)
legend('bottomleft',
       pch = 1,
       legend = "ITRDB chronologies")
```

![](https://github.com/ropensci/FedData/raw/master/inst/image/README-unnamed-chunk-13-1.png)
