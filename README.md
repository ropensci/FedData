FedData
========

[![Build Status](https://api.travis-ci.org/bocinsky/FedData.png)](https://travis-ci.org/bocinsky/FedData)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/FedData)](https://github.com/metacran/cranlogs.app)
[![cran version](http://www.r-pkg.org/badges/version/FedData)](https://cran.r-project.org/package=FedData)

`FedData` is an *R* package implementing functions to automate downloading geospatial data available from several federated data sources (mainly sources maintained by the US Federal government). Currently, the package allows for retrieval of five datasets: 

* The [National Elevation Dataset (NED)](http://ned.usgs.gov) digital elevation models (1 and 1/3 arc-second; USGS)
* The [National Hydrography Dataset (NHD)](http://nhd.usgs.gov) (USGS)
* The [Soil Survey Geographic (SSURGO) database](http://websoilsurvey.sc.egov.usda.gov/) from the National Cooperative Soil Survey (NCSS), which is led by the Natural Resources Conservation Service (NRCS) under the USDA,
* The [Global Historical Climatology Network (GHCN)](http://www.ncdc.noaa.gov/data-access/land-based-station-data/land-based-datasets/global-historical-climatology-network-ghcn), coordinated by National Climatic Data Center at NOAA,
* The [Daymet](https://daymet.ornl.gov/) gridded estimates of daily weather parameters for North America, version 3, available from the Oak Ridge National Laboratory's Distributed Active Archive Center (DAAC), and
* The [International Tree Ring Data Bank (ITRDB)](http://www.ncdc.noaa.gov/data-access/paleoclimatology-data/datasets/tree-ring), coordinated by National Climatic Data Center at NOAA.

Additional data sources are in the works, including global DEM resources ([ETOPO1](https://www.ngdc.noaa.gov/mgg/global/global.html), [STRM](http://www2.jpl.nasa.gov/srtm/)), global soils ([HWSD](http://webarchive.iiasa.ac.at/Research/LUC/External-World-soil-database/HTML/)), [MODIS](https://modis.gsfc.nasa.gov/) satellite data products, the [National Atlas](http://nationalmap.gov/small_scale/) (US only), [Natural Earth](http://www.naturalearthdata.com/), and [WorldClim](http://www.worldclim.org/).

This package is designed with the large-scale geographic information system (GIS) use-case in mind: cases where the use of dynamic web-services is impractical due to the scale (spatial and/or temporal) of analysis. It functions primarily as a means of downloading tiled or otherwise spatially-defined datasets; additionally, it can preprocess those datasets by extracting data within an area of interest (AoI), defined spatially. It relies heavily on the [**sp**](https://cran.r-project.org/package=sp), [**raster**](https://cran.r-project.org/package=raster), and [**rgdal**](https://cran.r-project.org/package=rgdal) packages.

This package has been built and tested on a source (Homebrew) install of *R* on Mac OS 10.12 (Sierra), and has been successfully run on Ubuntu 16.04.1 LTS and binary installs of *R* on Mac OS 10.12 and Windows 10.

### Development
+ [Kyle Bocinsky](http://bocinsky.io) - Crow Canyon Archaeological Center, Cortez, CO

### Contributors
+ [Dylan Beaudette](http://casoilresource.lawr.ucdavis.edu/people/dylan-e-beaudette/) - USDA-NRCS Soil Survey Office, Sonora, CA
+ [Scott Chamberlain](http://scottchamberlain.info/) - ROpenSci and Museum of Paleontology at UC Berkeley

### Install `FedData`
+ CRAN:
```r
install.packages('FedData')
```

+ Development version from GitHub:
```r
install.packages("devtools")
library(devtools)
install_github("bocinsky/FedData")
library(FedData)
```

### Demonstration
This demo script is available in the `/inst` folder at the location of the installed package.

#### Load `FedData` and define a study area
```r
# FedData Tester
library(FedData)
library(magrittr)

# Set a directory for testing
testDir <- "~/FedData Test"
# and create it if necessary
dir.create(testDir, showWarnings=F, recursive=T)
setwd(testDir)

# Extract data for the Village Ecodynamics Project "VEPIIN" study area:
# http://village.anth.wsu.edu
vepPolygon <- polygon_from_extent(raster::extent(672800,740000,4102000,4170000),
                                  proj4string="+proj=utm +datum=NAD83 +zone=12")
```

#### Get and plot the National Elevation Dataset for the study area
```r
# Get the NED (USA ONLY)
# Returns a raster
NED <- get_ned(template=vepPolygon,
               label="VEPIIN")
# Plot with raster::plot
raster::plot(NED)
```
![thing](inst/img/NED.png)

#### Get and plot the Daymet dataset for the study area
```r
# Get the DAYMET (North America only)
# Returns a raster
DAYMET <- get_daymet(template=vepPolygon,
               label="VEPIIN",
               elements = c("prcp","tmax"),
               years = 1980:1985)
# Plot with raster::plot
raster::plot(DAYMET$tmax$X1985.10.23)
```
![thing](inst/img/DAYMET.png)

#### Get and plot the daily GHCN precipitation data for the study area
```r
# Get the daily GHCN data (GLOBAL)
# Returns a list: the first element is the spatial locations of stations,
# and the second is a list of the stations and their daily data
GHCN.prcp <- get_ghcn_daily(template=vepPolygon, 
                            label="VEPIIN", 
                            elements=c('prcp'))
# Plot the NED again
raster::plot(NED)
# Plot the spatial locations
sp::plot(GHCN.prcp$spatial, pch=1, add=T)
legend('bottomleft', pch=1, legend="GHCN Precipitation Records")
```
![thing](inst/img/GHCN_prcp.png)

#### Get and plot the daily GHCN temperature data for the study area
```r
# Elements for which you require the same data
# (i.e., minimum and maximum temperature for the same days)
# can be standardized using standardize==T
GHCN.temp <- get_ghcn_daily(template = vepPolygon, 
                            label = "VEPIIN", 
                            elements = c('tmin','tmax'), 
                            years = 1980:1985,
                            standardize = T)
# Plot the NED again
raster::plot(NED)
# Plot the spatial locations
sp::plot(GHCN.temp$spatial, add=T, pch=1)
legend('bottomleft', pch=1, legend="GHCN Temperature Records")
```
![thing](inst/img/GHCN_temp.png)

#### Get and plot the National Hydrography Dataset for the study area
```r
# Get the NHD (USA ONLY)
NHD <- get_nhd(template=vepPolygon, 
               label="VEPIIN")
# Plot the NED again
raster::plot(NED)
# Plot the NHD data
NHD %>%
  lapply(sp::plot, col='black', add=T)
```
![thing](inst/img/NHD.png)


#### Get and plot the NRCS SSURGO data for the study area
```r
# Get the NRCS SSURGO data (USA ONLY)
SSURGO.VEPIIN <- get_ssurgo(template=vepPolygon, 
                     label="VEPIIN")
# Plot the NED again
raster::plot(NED)
# Plot the SSURGO mapunit polygons
plot(SSURGO.VEPIIN$spatial,
     lwd=0.1,
     add=T)
```
![thing](inst/img/SSURGO_VEP.png)

#### Get and plot the NRCS SSURGO data for particular soil survey areas
```r
# Or, download by Soil Survey Area names
SSURGO.areas <- get_ssurgo(template=c("CO670","CO075"), 
                           label="CO_TEST")

# Let's just look at spatial data for CO675
SSURGO.areas.CO675 <- SSURGO.areas$spatial[SSURGO.areas$spatial$AREASYMBOL=="CO075",]

# And get the NED data under them for pretty plotting
NED.CO675 <- get_ned(template=SSURGO.areas.CO675,
                            label="SSURGO_CO675")
               
# Plot the SSURGO mapunit polygons, but only for CO675
plot(NED.CO675)
plot(SSURGO.areas.CO675,
     lwd=0.1,
     add=T)
```
![thing](inst/img/SSURGO_areas.png)

#### Get and plot the ITRDB chronology locations in the study area
```r
# Get the ITRDB records
ITRDB <- get_itrdb(template=vepPolygon,
                        label="VEPIIN",
                        makeSpatial=T)
# Plot the NED again
raster::plot(NED)
# Map the locations of the tree ring chronologies
plot(ITRDB$metadata, pch=1, add=T)
legend('bottomleft', pch=1, legend="ITRDB chronologies")
```
![thing](inst/img/ITRDB.png)
========

### Acknowledgements
This package is a product of SKOPE ([Synthesizing Knowledge of Past Environments](http://www.envirecon.org)) and the Village Ecodynamics Project. This software is licensed under the [MIT license](https://opensource.org/licenses/MIT).