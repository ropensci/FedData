pkgname <- "FedData"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('FedData')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("FedData-package")
### * FedData-package

flush(stderr()); flush(stdout())

### Name: FedData-package
### Title: Scripts to automate downloading geospatial data available from
###   the several federated data sources
### Aliases: FedData-package FedData
### Keywords: package

### ** Examples

## Not run: 
##D install.packages("FedData")
##D library(FedData)
##D 
##D # Get a random contiguous USA county for testing
##D wgetDownload(
##D   "http://dds.cr.usgs.gov/pub/data/nationalatlas/countyp010g.shp_nt00934.tar.gz"
##D   ,destdir=getwd())
##D 
##D untar("./countyp010g.shp_nt00934.tar.gz")
##D county <- rgdal::readOGR(".","countyp010g")
##D county <- county[!(county$STATE ##D 
##D county <- county[sample(1:length(county),1),]
##D 
##D # Get the NED (USA ONLY)
##D # Returns a raster
##D NED <- getNED(template=county,
##D   label=paste(county$STATE,'_',county$NAME, sep=''), res='1')
##D 
##D # Get the daily GHCN data (GLOBAL)
##D # Returns a list: the first element is the spatial locations of stations,
##D # and the second is a list of the stations and their daily data
##D GHCN.prcp <- getGHCNDaily(template=county, 
##D   label=paste(county$STATE,'_',county$NAME, sep=''), 
##D   elements=c('prcp'), 
##D   standardize=F)
##D   
##D GHCN.temp <- getGHCNDaily(template=county, 
##D   label=paste(county$STATE,'_',county$NAME, sep=''), 
##D   elements=c('tmin','tmax'), 
##D   standardize=T)
##D 
##D # Get the NHD (USA ONLY)
##D NHD <- getNHD(template=county, 
##D   label=paste(county$STATE,'_',county$NAME, sep=''))
##D 
##D # Get the NRCS SSURGO data (USA ONLY)
##D SSURGO <- getSSURGO(template=county, 
##D   label=paste(county$STATE,'_',county$NAME, sep=''))
##D   
##D # Get the ITRDB data
##D ITRDB <- getITRDB(template=county, 
##D   label=paste(county$STATE,'_',county$NAME, sep=''))
##D 
## End(Not run)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
