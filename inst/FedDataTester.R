# FedData Tester
library(FedData)
pkg_test("png")
pkg_test("RColorBrewer")

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
plot(GHCN.prcp$spatial)

# Elements for which you require the same data
# (i.e., minimum and maximum temperature for the same days)
# can be standardized using standardize==T
GHCN.temp <- get_ghcn_daily(template=vepPolygon, 
                            label="VEPIIN", 
                            elements=c('tmin','tmax'), 
                            standardize=T)
# Plot the spatial locations
plot(GHCN.temp$spatial)


# Get the NHD (USA ONLY)
NHD <- get_nhd(template=vepPolygon, 
               label="VEPIIN")
# Plot the NHD data
plot(NHD$NHDFlowline)
plot(NHD$NHDLine, add=T)
plot(NHD$NHDArea, col='black', add=T)
plot(NHD$NHDWaterbody, col='black', add=T)


# Get the NRCS SSURGO data (USA ONLY)
SSURGO.VEPIIN <- get_ssurgo(template=vepPolygon, 
                     label="VEPIIN")
# Plot the SSURGO mapunit polygons
plot(SSURGO.VEPIIN$spatial, lwd=0.1)

# Or, download by Soil Survey Area names
SSURGO.areas <- get_ssurgo(template=c("CO670","CO075"), 
                     label="CO_TEST")
# Plot the SSURGO mapunit polygons, but only for CO675
plot(SSURGO.areas$spatial[SSURGO.areas$spatial$AREASYMBOL=="CO075",], lwd=0.1)


# Get the ITRDB records
ITRDB <- get_itrdb(template=vepPolygon,
                        label="VEPIIN",
                        makeSpatial=T)
# Map the locations of the tree ring chronologies
plot(ITRDB$metadata)

