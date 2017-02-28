# FedData Tester
library(FedData)
library(magrittr)

# Set a directory for testing
testDir <- "./FedData Test"

dir.create(testDir, showWarnings=F, recursive=T)
setwd(testDir)

# Extract data for the test study area:
# http://testarchaeology.org/
testPolygon <- polygon_from_extent(raster::extent(-110,-108,36.5,37.5),
                                   proj4string="+proj=longlat")

testPolygon %>%
  as("SpatialPolygonsDataFrame") %>%
writeOGR(dsn = "./FedData Test/",
         layer = "testPolygon",
         driver = "ESRI Shapefile",
         overwrite = T)

testPolygon %>%
  as("SpatialPolygonsDataFrame") %>%
  writeOGR(dsn = "./FedData Test/testPolygon.json",
           layer = "testPolygon",
           driver = "GeoJSON",
           overwrite = T)

unlink("./nc.shp")
nc = st_read(system.file("shape/nc.shp", package="sf"))
st_write(nc, "./nc.shp")
st_read("./nc.shp")

unlink("./test.geojson")
testPolygon %>%
  st_as_sf() %>%
  mutate(test = 1) %>%
  st_write("./test.geojson")

ogrDrivers() %>%
  dplyr::filter(write)

testPolygon %>%
  geojsonio::geojson_write(file = "./FedData Test/testPolygon")

test <- geojsonio::geojson_read("./FedData Test/testPolygon.geojson", method = "local", parse = T, what = "sp") 

test <- readOGR("./FedData Test/testPolygon.geojson")


# Get the NED (USA ONLY)
# Returns a raster
NED <- get_ned(template=testPolygon,
               label="test")
# Plot with raster::plot
raster::plot(NED)


# Get the DAYMET (North America only)
# Returns a raster
DAYMET <- get_daymet(template=testPolygon,
                     label="test",
                     elements = c("prcp","tmax"),
                     years = 1980:1985)
# Plot with raster::plot
raster::plot(DAYMET$tmax$X1985.10.23)


# Get the daily GHCN data (GLOBAL)
# Returns a list: the first element is the spatial locations of stations,
# and the second is a list of the stations and their daily data
GHCN.prcp <- get_ghcn_daily(template=testPolygon, 
                            label="test", 
                            elements=c('prcp'))
# Plot the NED again
raster::plot(NED)
# Plot the spatial locations
sp::plot(GHCN.prcp$spatial, pch=1, add=T)
legend('bottomleft', pch=1, legend="GHCN Precipitation Records")

# Elements for which you require the same data
# (i.e., minimum and maximum temperature for the same days)
# can be standardized using standardize==T
GHCN.temp <- get_ghcn_daily(template = testPolygon, 
                            label = "test", 
                            elements = c('tmin','tmax'), 
                            years = 1980:1985,
                            standardize = T)
# Plot the NED again
raster::plot(NED)
# Plot the spatial locations
sp::plot(GHCN.temp$spatial, add=T, pch=1)
legend('bottomleft', pch=1, legend="GHCN Temperature Records")


# Get the NHD (USA ONLY)
NHD <- get_nhd(template=testPolygon, 
               label="test")
# Plot the NED again
raster::plot(NED)
# Plot the NHD data
NHD %>%
  lapply(sp::plot, col='black', add=T)

# Get the NRCS SSURGO data (USA ONLY)
SSURGO.test <- get_ssurgo(template=testPolygon, 
                          label="test")
# Plot the NED again
raster::plot(NED)
# Plot the SSURGO mapunit polygons
plot(SSURGO.test$spatial,
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
ITRDB <- get_itrdb(template=testPolygon,
                   label="test",
                   makeSpatial=T)
# Plot the NED again
raster::plot(NED)
# Map the locations of the tree ring chronologies
plot(ITRDB$metadata, pch=1, add=T)
legend('bottomleft', pch=1, legend="ITRDB chronologies")

