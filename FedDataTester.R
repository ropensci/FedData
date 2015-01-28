# FedData Tester
# devtools::install_github("bocinsky/FedData")
library(FedData)

test.raw.dir <- "/Users/Bocinsky/Desktop/FedDataTest/"

# Get a random contiguous USA county for testing
county <- readOGR("/Volumes/DATA/DATA/NATIONAL_ATLAS/countyp010","countyp010g")
county <- county[!(county$STATE %in% c("AK","HI","PR","VI")),]
county <- county[sample(1:length(county),1),]

# Get the NED
NED <- getNED(template=county, label=paste(county$STATE,'_',county$NAME, sep=''), res='1', raw.dir=test.raw.dir, force.redo=F)

# Get the daily GHCN data
GHCN.prcp <- getGHCNDaily(template=county, elements=c('prcp'), raw.dir=test.raw.dir, standardize=F)
GHCN.temp <- getGHCNDaily(template=county, elements=c('tmin','tmax'), raw.dir=test.raw.dir, standardize=T)

# Get the NHD
NHD <- getNHD(template=county, label=paste(county$STATE,'_',county$NAME, sep=''), raw.dir=test.raw.dir, force.redo=F)









plot(NED)
plot(county, add=T)
plot(GHCN.prcp[[1]], pch=17, add=T)
plot(GHCN.temp[[1]], pch=19, add=T)

