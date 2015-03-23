# FedData
An *R* package implementing functions to automate downloading geospatial data available from several federated data sources (mainly sources maintained by the US Federal government). Currently, the package allows for retrieval of five datasets: 

* The [National Elevation Dataset (NED)](http://ned.usgs.gov) digital elevation models (1 and 1/3 arc-second; USGS)
* The [National Hydrography Dataset (NHD)](http://nhd.usgs.gov) (USGS)
* The [Soil Survey Geographic (SSURGO) database](http://websoilsurvey.sc.egov.usda.gov/) from the National Cooperative Soil Survey (NCSS), which is led by the Natural Resources Conservation Service (NRCS) under the USDA,
* The [Global Historical Climatology Network (GHCN)](http://www.ncdc.noaa.gov/data-access/land-based-station-data/land-based-datasets/global-historical-climatology-network-ghcn), coordinated by National Climatic Data Center at NOAA, and
* The [International Tree Ring Data Bank (ITRDB)](http://www.ncdc.noaa.gov/data-access/paleoclimatology-data/datasets/tree-ring), coordinated by National Climatic Data Center at NOAA.

Additional data sources are in the works, including global DEM resources ([ETOPO1](https://www.ngdc.noaa.gov/mgg/global/global.html), [STRM](http://www2.jpl.nasa.gov/srtm/)), global soils ([HWSD](http://webarchive.iiasa.ac.at/Research/LUC/External-World-soil-database/HTML/)), [MODIS](http://modis.gsfc.nasa.gov) satellite data products, the [National Atlas](http://nationalmap.gov/small_scale/) (US only), [Natural Earth](http://www.naturalearthdata.com), [PRISM](http://www.prism.oregonstate.edu), and [WorldClim](http://www.worldclim.org).

This package is designed with the large-scale geographic information system (GIS) use-case in mind: cases where the use of dynamic web-services is impractical due to the scale (spatial and/or temporal) of analysis. It functions primarily as a means of downloading tiled or otherwise spaticially-defined datasets; additionally, it can preprocess those datasets by extracting data within an area of interest (AoI), defined spatially. It relies heavily on the [**sp**](http://cran.r-project.org/package=sp), [**raster**](http://cran.r-project.org/package=raster), and [**rgdal**](http://cran.r-project.org/package=rgdal) packages, and requires three command line tools be installed by the user (and accesible through `system()` calls in *R*): [wget](https://www.gnu.org/software/wget/) (for downloading with timestamping), [GDAL](http://www.gdal.org) (for manipulating raster and vector spatial data), and [mdbtools](http://mdbtools.sourceforge.net) (for extracting data from access databases).

I **strongly** recommend [Homebrew](http://brew.sh) for installing *R* and each of these command-line tools. I suggest using these build parameters:
`brew install wget`
`brew install mdbtools`
`brew install gdal --complete --enable-mdb --enable-opencl --enable-unsupported --with-libkml --with-python`
`brew install r --with-openblas`

This package has been built and tested on Mac OS 10.10 (Yosemite), and has been successfully run on an Ubuntu Linux cluster.

### Installation
To install, use the following command in *R*:

`devtools::install_github("bocinsky/FedData")`

A demo script is available in the `/inst` folder at the location of the installed package.

### Acknowledgements
This package is a product of SKOPE ([Synthesized Knowledge of Past Environments](https://www.envirecon.org)) and the [Village Ecodynamics Project](http://village.anth.wsu.edu). This software is licensed under [GPLv3](https://www.gnu.org/copyleft/gpl.html).