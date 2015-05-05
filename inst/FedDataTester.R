# FedData Tester

# The development version is usually more up to date, but less stable
devtools::install_github("bocinsky/FedData")

# Install the CRAN version
# install.packages("FedData")

library(FedData)
pkgTest("raster")
pkgTest("png")
pkgTest("RColorBrewer")

setwd("~/Desktop/FedData Test")

# Get a random contiguous USA county for testing
curlDownload("http://dds.cr.usgs.gov/pub/data/nationalatlas/countyp010g.shp_nt00934.tar.gz",destdir=getwd())
untar("./countyp010g.shp_nt00934.tar.gz",exdir = "./countyp010g/")
county <- rgdal::readOGR("./countyp010g/countyp010g.shp","countyp010g")
county <- county[!(county$STATE %in% c("AK","VI","PR","HI")),]
county <- county[sample(1:length(county),1),]
# county <- county[which(county$NAME=='Napa'),]

# Get the NED (USA ONLY)
# Returns a raster
NED <- getNED(template=county, label=paste(county$STATE,'_',county$NAME, sep=''), res='1')

# Get the daily GHCN data (GLOBAL)
# Returns a list: the first element is the spatial locations of stations,
# and the second is a list of the stations and their daily data
GHCN.prcp <- getGHCNDaily(template=county, label=paste(county$STATE,'_',county$NAME, sep=''), elements=c('prcp'), standardize=F)
GHCN.temp <- getGHCNDaily(template=county, label=paste(county$STATE,'_',county$NAME, sep=''), elements=c('tmin','tmax'), standardize=T)

# Get the NHD (USA ONLY)
NHD <- getNHD(template=county, label=paste(county$STATE,'_',county$NAME, sep=''))

# Get the NRCS SSURGO data (USA ONLY)
SSURGO <- getSSURGO(template=county, label=paste(county$STATE,'_',county$NAME, sep=''))

# Get the ITRDB records
itrdb.data <- getITRDB(template=county, 
                       label=paste(county$STATE,'_',county$NAME, sep=''),
                       makeSpatial=T)



slope <- terrain(NED, opt='slope')
aspect <- terrain(NED, opt='aspect')
NED.hill <- hillShade(slope, aspect, 40, 230)

NED <- raster::mask(NED,county)
NED.hill <- raster::mask(NED.hill,county)

plot(NED)
plot(county, add=T)
plot(NHD$NHDFlowline, col='gray50', border='gray50', add=T)
plot(NHD$NHDWaterbody, col='gray50', border='gray50', add=T)
plot(NHD$NHDArea, col='gray50', border='gray50', add=T)
plot(county, add=T)
plot(GHCN.prcp[[1]], pch=17, add=T)
plot(GHCN.temp[[1]], pch=19, add=T)
