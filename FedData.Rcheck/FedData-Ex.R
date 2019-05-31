pkgname <- "FedData"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "FedData-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('FedData')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("get_daymet")
### * get_daymet

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_daymet
### Title: Download and crop the 1-km DAYMET daily weather dataset.
### Aliases: get_daymet

### ** Examples

## Not run: 
##D # Extract data for the Village Ecodynamics Project 'VEPIIN' study area:
##D # http://village.anth.wsu.edu
##D vepPolygon <- polygon_from_extent(raster::extent(672800,740000,4102000,4170000),
##D      proj4string='+proj=utm +datum=NAD83 +zone=12')
##D 
##D # Get the DAYMET (North America only)
##D # Returns a list of raster bricks
##D DAYMET <- get_daymet(template=vepPolygon,
##D                      label='VEPIIN',
##D                      elements = c('prcp','tmin','tmax'),
##D                      years = 1980:1985)
##D 
##D # Plot with raster::plot
##D plot(DAYMET$tmin$X1985.10.23)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_daymet", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get_ghcn_daily")
### * get_ghcn_daily

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_ghcn_daily
### Title: Download and crop the Global Historical Climate Network-Daily
###   data.
### Aliases: get_ghcn_daily

### ** Examples

## Not run: 
##D # Extract data for the Village Ecodynamics Project 'VEPIIN' study area:
##D # http://village.anth.wsu.edu
##D vepPolygon <- polygon_from_extent(raster::extent(672800,740000,4102000,4170000), 
##D      proj4string='+proj=utm +datum=NAD83 +zone=12')
##D 
##D # Get the daily GHCN data (GLOBAL)
##D # Returns a list: the first element is the spatial locations of stations,
##D # and the second is a list of the stations and their daily data
##D GHCN.prcp <- get_ghcn_daily(template=vepPolygon, label='VEPIIN', elements=c('prcp'))
##D 
##D # Plot the VEP polygon
##D plot(vepPolygon)
##D 
##D # Plot the spatial locations
##D plot(GHCN.prcp$spatial, pch=1, add=T)
##D legend('bottomleft', pch=1, legend='GHCN Precipitation Records')
##D 
##D # Elements for which you require the same data
##D # (i.e., minimum and maximum temperature for the same days)
##D # can be standardized using standardize==T
##D GHCN.temp <- get_ghcn_daily(template=vepPolygon, 
##D      label='VEPIIN', 
##D      elements=c('tmin','tmax'), 
##D      standardize=T)
##D 
##D # Plot the VEP polygon
##D plot(vepPolygon)
##D 
##D # Plot the spatial locations
##D plot(GHCN.temp$spatial, pch=1, add=T)
##D legend('bottomleft', pch=1, legend='GHCN Temperature Records')
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_ghcn_daily", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get_itrdb")
### * get_itrdb

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_itrdb
### Title: Download the latest version of the ITRDB, and extract given
###   parameters.
### Aliases: get_itrdb

### ** Examples

## Not run: 
##D # Extract data for the Village Ecodynamics Project 'VEPIIN' study area:
##D # http://village.anth.wsu.edu
##D vepPolygon <- polygon_from_extent(raster::extent(672800,740000,4102000,4170000), 
##D      proj4string='+proj=utm +datum=NAD83 +zone=12')
##D 
##D # Get the ITRDB records
##D ITRDB <- get_itrdb(template=vepPolygon, label='VEPIIN', makeSpatial=T)
##D 
##D # Plot the VEP polygon
##D plot(vepPolygon)
##D 
##D # Map the locations of the tree ring chronologies
##D plot(ITRDB$metadata, pch=1, add=T)
##D legend('bottomleft', pch=1, legend='ITRDB chronologies')
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_itrdb", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get_ned")
### * get_ned

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_ned
### Title: Download and crop the 1 (~30 meter) or 1/3 (~10 meter)
###   arc-second National Elevation Dataset.
### Aliases: get_ned

### ** Examples

## Not run: 
##D # Extract data for the Village Ecodynamics Project 'VEPIIN' study area:
##D # http://village.anth.wsu.edu
##D vepPolygon <- polygon_from_extent(raster::extent(672800,740000,4102000,4170000), 
##D      proj4string='+proj=utm +datum=NAD83 +zone=12')
##D 
##D # Get the NED (USA ONLY)
##D # Returns a raster
##D NED <- get_ned(template=vepPolygon, label='VEPIIN')
##D 
##D # Plot with raster::plot
##D plot(NED)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_ned", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get_nhd")
### * get_nhd

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_nhd
### Title: Download and crop the National Hydrography Dataset.
### Aliases: get_nhd

### ** Examples

## Not run: 
##D # Extract data for the Village Ecodynamics Project 'VEPIIN' study area:
##D # http://village.anth.wsu.edu
##D vepPolygon <- polygon_from_extent(raster::extent(672800,740000,4102000,4170000), 
##D      proj4string='+proj=utm +datum=NAD83 +zone=12')
##D 
##D # Get the NHD (USA ONLY)
##D NHD <- get_nhd(template=vepPolygon, label='VEPIIN')
##D 
##D # Plot the VEP polygon
##D plot(vepPolygon)
##D 
##D # Plot the NHD data
##D plot(NHD$NHDFlowline, add=T)
##D plot(NHD$NHDLine, add=T)
##D plot(NHD$NHDArea, col='black', add=T)
##D plot(NHD$NHDWaterbody, col='black', add=T)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_nhd", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get_nlcd")
### * get_nlcd

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_nlcd
### Title: Download and crop the National Land Cover Database.
### Aliases: get_nlcd

### ** Examples

## Not run: 
##D # Extract data for the Village Ecodynamics Project 'VEPIIN' study area:
##D # http://village.anth.wsu.edu
##D vepPolygon <- polygon_from_extent(raster::extent(672800,740000,4102000,4170000), 
##D      proj4string='+proj=utm +datum=NAD83 +zone=12')
##D 
##D # Get the NLCD (USA ONLY)
##D # Returns a raster
##D NLCD <- get_nlcd(template=vepPolygon, label='VEPIIN')
##D 
##D # Plot with raster::plot
##D plot(NLCD)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_nlcd", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get_ssurgo")
### * get_ssurgo

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_ssurgo
### Title: Download and crop data from the NRCS SSURGO soils database.
### Aliases: get_ssurgo

### ** Examples

## Not run: 
##D # Extract data for the Village Ecodynamics Project 'VEPIIN' study area:
##D # http://village.anth.wsu.edu
##D vepPolygon <- polygon_from_extent(raster::extent(672800,740000,4102000,4170000), 
##D      proj4string='+proj=utm +datum=NAD83 +zone=12')
##D 
##D # Get the NRCS SSURGO data (USA ONLY)
##D SSURGO.VEPIIN <- get_ssurgo(template=vepPolygon, label='VEPIIN')
##D 
##D # Plot the VEP polygon
##D plot(vepPolygon)
##D 
##D # Plot the SSURGO mapunit polygons
##D plot(SSURGO.VEPIIN$spatial, lwd=0.1, add=T)
##D 
##D # Or, download by Soil Survey Area names
##D SSURGO.areas <- get_ssurgo(template=c('CO670','CO075'), label='CO_TEST')
##D 
##D # Let's just look at spatial data for CO675
##D SSURGO.areas.CO675 <- SSURGO.areas$spatial[SSURGO.areas$spatial$AREASYMBOL=='CO075',]
##D 
##D # And get the NED data under them for pretty plotting
##D NED.CO675 <- get_ned(template=SSURGO.areas.CO675, label='SSURGO_CO675')
##D 
##D # Plot the SSURGO mapunit polygons, but only for CO675
##D plot(NED.CO675)
##D plot(SSURGO.areas.CO675, lwd=0.1, add=T)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_ssurgo", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pal_nlcd")
### * pal_nlcd

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pal_nlcd
### Title: NLCD colour map palettes
### Aliases: pal_nlcd

### ** Examples

## Not run: 
##D # Extract data for the Village Ecodynamics Project 'VEPIIN' study area:
##D # http://village.anth.wsu.edu
##D vepPolygon <- polygon_from_extent(raster::extent(672800,740000,4102000,4170000), 
##D      proj4string='+proj=utm +datum=NAD83 +zone=12')
##D 
##D NLCD <- get_nlcd(template=vepPolygon, label='VEPIIN')
##D NLCD <- as.matrix(table(raster::values(NLCD)))
##D cols <- dplyr::filter(pal_nlcd(), code %in% row.names(NLCD))
##D par(xpd = TRUE, mar = c(10, 3, 2, 1))
##D barplot(NLCD, beside = FALSE, col = cols$color) 
##D legend("bottom", legend = cols$description, fill = cols$color, 
##D        ncol = 2, inset = c(0, -0.6))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pal_nlcd", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
