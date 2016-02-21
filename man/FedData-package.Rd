\name{FedData-package}
\alias{FedData-package}
\alias{FedData}
\docType{package}
\title{
Scripts to automate downloading geospatial data available from the several federated data sources
}
\description{
This package contains scripts to automate downloading geospatial data available from the several federated data sources (mainly sources maintained by the US Federal government). Currently, the package allows for retrieval of five datasets: 
\itemize{
\item The \href{http://ned.usgs.gov}{National Elevation Dataset} digital elevation models (1 and 1/3 arc-second; USGS)
\item The \href{http://nhd.usgs.gov}{National Hydrography Dataset} (USGS)
\item The \href{http://websoilsurvey.sc.egov.usda.gov/}{Soil Survey Geographic} (SSURGO) database from the National Cooperative Soil Survey (NCSS), which is led by the Natural Resources Conservation Service (NRCS) under the USDA, and
\item The \href{http://www.ncdc.noaa.gov/data-access/land-based-station-data/land-based-datasets/global-historical-climatology-network-ghcn}{Global Historical Climatology Network} (GHCN), coordinated by National Climatic Data Center at NOAA.
\item The \href{http://www.ncdc.noaa.gov/data-access/paleoclimatology-data/datasets/tree-ring}{International Tree Ring Data Bank} (ITRDB), coordinated by National Climatic Data Center at NOAA.
}

Additional data sources are in the works, including global DEM resources (ETOPO1, ETOPO5, ETOPO30, SRTM), global soils (HWSD), MODIS satellite data products, the National Atlas (US), Natural Earth, PRISM, and WorldClim.
}
\details{
\tabular{ll}{
Package: \tab FedData\cr
Type: \tab FedData\cr
Version: \tab 2.0.6\cr
Date: \tab 2016-02-21\cr
License: \tab MIT\cr
}

This package is designed with the large-scale GIS use-case in mind: cases where the use of dynamic web-services is impractical due to the scale (spatial and/or temporal) of analysis. It functions primarily as a means of downloading tiled or otherwise spaticially-defined datasets; additionally, it can preprocess those datasets by extracting data within an area of interest (AoI), defined spatially. It relies heavily on the \pkg{sp}, \pkg{raster}, and \pkg{rgdal} packages.

There are three general types of methods available for each dataset (and others for particular datasets):

\tabular{ll}{
get...: \tab High-level function that allows the user to define an AoI ("template") and returns the dataset cropped/masked to the area of interest. Examples: get_ned(template, ...), get_nhd(template, ...), get_nrcs(template, ...).\cr

extract...: \tab Mid-level functions that automate extraction of tabular data from databases (such as the SSURGO soils tabular data) that may or may not be linked to spatial data. Often called by "get" functions.\cr

download...: \tab Low-level functions that automate downloading of raw tabular and spatial data from databases. Downloading is done at the smallest data-unit of each dataset; e.g., 1x1 degree tiles for the NED, HUC4 subregions for the NHD, and Study Areas for the SSURGO soils data. Downloading is sometimes (but not always) performed with timestamping, such that source files are only downloaded if the remote files are more recent than the local files. Often called by "get" functions. Returns a local path to the downloaded file. \cr

}

Additionally, most functions can be forced to "start fresh" in downloading or processing data by specifying \code{force.redo=TRUE} in the function call.

}

\author{
R. Kyle Bocinsky <bocinsky@gmail.com>
}

\keyword{ package }

\references{

Gesch, D.B. (2007) The National Elevation Dataset, in Maune, D., ed., \emph{Digital Elevation Model Technologies and Applications: The DEM Users Manual}. 2nd Edition. American Society for Photogrammetry and Remote Sensing, Bethesda, Maryland.

Gesch, D., Oimoen, M., Greenlee, S., Nelson, C., Steuck, M., and Tyler, D. (2002) The National Elevation Dataset. \emph{Photogrammetric Engineering and Remote Sensing} 68(1):5--11.

Grissino-Mayer HD, Fritts HC. (1997) The International Tree-Ring Data Bank: An enhanced global database serving the global scientific community. \emph{The Holocene} 7(2):235--238.

Menne, M., Durre, I., Korzeniewski, B., McNeal, S., Thomas, K., Yin, X., Anthony, S., Ray, R., Vose, R., B.E.Gleason, and Houston, T. (2012) \emph{Global Historical Climatology Network-Daily (GHCN-Daily), Version 3}. \url{http://doi.org/10.7289/V5D21VHZ}.

Soil Survey Staff, Natural Resources Conservation Service, United States Department of Agriculture. Soil Survey Geographic (SSURGO) Database. Available online at \url{http://sdmdataaccess.nrcs.usda.gov/}.

}

\examples{
\dontrun{
# FedData Tester
library(FedData)

# Set a directory for testing
testDir <- "~/FedData Test"

dir.create(testDir, showWarnings=F, recursive=T)
setwd(testDir)

# Extract data for the Village Ecodynamics Project "VEPIIN" study area:
# http://village.anth.wsu.edu
vepPolygon <- polygon_from_extent(raster::extent(672800,740000,4102000,4170000),
                                  proj4string="+proj=utm +datum=NAD83 +zone=12")

# Get the NED (USA ONLY)
# Returns a raster
NED <- get_ned(template=vepPolygon,
               label="VEPIIN")
# Plot with raster::plot
plot(NED)


# Get the daily GHCN data (GLOBAL)
# Returns a list: the first element is the spatial locations of stations,
# and the second is a list of the stations and their daily data
GHCN.prcp <- get_ghcn_daily(template=vepPolygon, 
                            label="VEPIIN", 
                            elements=c('prcp'))
# Plot the spatial locations
plot(GHCN.prcp$spatial, pch=1, add=T)
legend('bottomleft', pch=1, legend="GHCN Precipitation Records")

# Elements for which you require the same data
# (i.e., minimum and maximum temperature for the same days)
# can be standardized using standardize==T
GHCN.temp <- get_ghcn_daily(template=vepPolygon, 
                            label="VEPIIN", 
                            elements=c('tmin','tmax'), 
                            standardize=T)
# Plot the NED again
raster::plot(NED)
# Plot the spatial locations
plot(GHCN.temp$spatial, add=T, pch=1)
legend('bottomleft', pch=1, legend="GHCN Temperature Records")


# Get the NHD (USA ONLY)
NHD <- get_nhd(template=vepPolygon, 
               label="VEPIIN")
# Plot the NED again
raster::plot(NED)
# Plot the NHD data
plot(NHD$NHDFlowline, add=T)
plot(NHD$NHDLine, add=T)
plot(NHD$NHDArea, col='black', add=T)
plot(NHD$NHDWaterbody, col='black', add=T)

# Get the NRCS SSURGO data (USA ONLY)
SSURGO.VEPIIN <- get_ssurgo(template=vepPolygon, 
                     label="VEPIIN")
# Plot the NED again
raster::plot(NED)
# Plot the SSURGO mapunit polygons
plot(SSURGO.VEPIIN$spatial,
     lwd=0.1,
     add=T)

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


# Get the ITRDB records
ITRDB <- get_itrdb(template=vepPolygon,
                        label="VEPIIN",
                        makeSpatial=T)
# Plot the NED again
raster::plot(NED)
# Map the locations of the tree ring chronologies
plot(ITRDB$metadata, pch=1, add=T)
legend('bottomleft', pch=1, legend="ITRDB chronologies")
}
}