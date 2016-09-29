# FedData Tester
library(FedData)
library(magrittr)

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
raster::plot(NED)


# Get the DAYMET (North America only)
# Returns a raster
DAYMET <- get_daymet(template=vepPolygon,
               label="VEPIIN",
               elements = c("prcp","tmax"),
               years = 1980:1985)
# Plot with raster::plot
raster::plot(DAYMET$tmax$X1985.10.23)


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


# Get the NHD (USA ONLY)
NHD <- get_nhd(template=vepPolygon, 
               label="VEPIIN")
# Plot the NED again
raster::plot(NED)
# Plot the NHD data
NHD %>%
  lapply(sp::plot, col='black', add=T)

# Get the NRCS SSURGO data (USA ONLY)
SSURGO.VEPIIN <- get_ssurgo(template=vepPolygon, 
                     label="VEPIIN")
# Plot the NED again
raster::plot(NED)
# Plot the SSURGO mapunit polygons
plot(SSURGO.VEPIIN$spatial,
     lwd=0.1,
     add=T)

# # Or, download by Soil Survey Area names
# SSURGO.areas <- get_ssurgo(template=c("CO670","CO075"), 
#                      label="CO_TEST")

# # Let's just look at spatial data for CO675
# SSURGO.areas.CO675 <- SSURGO.areas$spatial[SSURGO.areas$spatial$AREASYMBOL=="CO075",]
# 
# # And get the NED data under them for pretty plotting
# NED.CO675 <- get_ned(template=SSURGO.areas.CO675,
#                      label="SSURGO_CO675")
# 
# # Plot the SSURGO mapunit polygons, but only for CO675
# plot(NED.CO675)
# plot(SSURGO.areas.CO675,
#      lwd=0.1,
#      add=T)


# Get the ITRDB records
ITRDB <- get_itrdb(template=vepPolygon,
                        label="VEPIIN",
                        makeSpatial=T)
# Plot the NED again
raster::plot(NED)
# Map the locations of the tree ring chronologies
plot(ITRDB$metadata, pch=1, add=T)
legend('bottomleft', pch=1, legend="ITRDB chronologies")

