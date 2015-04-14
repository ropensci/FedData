# FedData Tester
devtools::install_github("bocinsky/FedData")
library(FedData)
library(raster)
library(png)
library(RColorBrewer)

setwd("~/Desktop/FedData Test")

# Get a random contiguous USA county for testing
curlDownload("http://dds.cr.usgs.gov/pub/data/nationalatlas/countyp010g.shp_nt00934.tar.gz",destdir=getwd())
untar("./countyp010g.shp_nt00934.tar.gz")
county <- rgdal::readOGR(".","countyp010g")
county <- county[!(county$STATE %in% c("AK","VI","PR","HI")),]
county <- county[sample(1:length(county),1),]
# county <- county[which(county$NAME=='Napa'),]

# Get the NED (USA ONLY)
# Returns a raster
NED <- getNED(template=county, label=paste(county$STATE,'_',county$NAME, sep=''), raw.dir="/Users/Bocinsky/Desktop/FedData Test/RAW/NED/",extraction.dir="/Users/Bocinsky/Desktop/FedData Test/EXTRACTIONS/NED/", res='1')

# Get the daily GHCN data (GLOBAL)
# Returns a list: the first element is the spatial locations of stations,
# and the second is a list of the stations and their daily data
GHCN.prcp <- getGHCNDaily(template=county, label=paste(county$STATE,'_',county$NAME, sep=''), elements=c('prcp'), standardize=F)
GHCN.temp <- getGHCNDaily(template=county, label=paste(county$STATE,'_',county$NAME, sep=''), elements=c('tmin','tmax'), standardize=T)

# Get the NHD (USA ONLY)
NHD <- getNHD(template=county, label=paste(county$STATE,'_',county$NAME, sep=''))

# Get the NRCS SSURGO data (USA ONLY)
SSURGO <- getSSURGO(template=county, label=paste(county$STATE,'_',county$NAME, sep=''))

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

# plot(SSURGO[['spatial']])

wine <- merge(SSURGO[['spatial']],SSURGO[['tabular']]$mucropyld[SSURGO[['tabular']]$mucropyld$cropname=='Wine grapes',],by.x='MUKEY',by.y="mukey",all=T)[,c("nonirryield.r","irryield.r")]
wine$yield.max <- pmax(wine$irryield.r,wine$nonirryield.r,na.rm=T)




# Plot NED
colors <- paste0(colorRampPalette(brewer.pal(9, "Greens"),bias=2)(1000),"CC")
plot(county)
plot(NED.hill, col=grey(30:100/100), maxpixels=2000000, useRaster=T, legend=FALSE,  xlab="", ylab="", axes=FALSE, main='', add=T)
plot(NED, col=colors, maxpixels=2000000, useRaster=T, legend=FALSE,  xlab="", ylab="", axes=FALSE, main='', add=T)
plot(county, add=T)

# Plot NHD
plot(NHD$NHDFlowline)
sp::plot(NHD$NHDWaterbody, col='black', add=T)
sp::plot(NHD$NHDArea, col='black', add=T)

# Plot average potential wine grape yield
colors <- paste0(colorRampPalette(brewer.pal(9, "Reds"))(1000),"CC")
wine$colors <- colors[wine$yield.max*100]
plot(county)
plot(NED.hill, col=grey(70:100/100), maxpixels=2000000, useRaster=T, legend=FALSE,  xlab="", ylab="", axes=FALSE, main='', add=T)
sp::plot(wine, col=wine$colors,border=NA, add=T)
plot(county, add=T)





itrdb.data <- getITRDB(recon.years=1:2000, 
                       calib.years=1924:1983, 
                       species=NULL, 
                       measurement.type="Ring Width", 
                       chronology.type="Standard", 
                       makeSpatial=T)


