## These are functions for downloading and creating the elevation model
## for the Village Ecodynamics Project simulation. UTM coordinates 
## representing the four sides of a simulated landscape are input in Master.R, 
## and a GeoTiff of the landscape is output in the NAD83 projection using the
## WGS84 ellipsoid.

## Author: R. Kyle Bocinsky
## Date: 09/14/2013

extractNED <- function(template, label, raw.dir, extraction.dir=NULL, res, drain=F, NHD.raw.dir=NULL, NRCS.raw.dir=NULL, force.redo=F){
  if(is.null(extraction.dir)){
    extraction.dir <- paste(raw.dir,"/EXTRACTIONS",sep='')
  }
  
  dir.create(paste(extraction.dir,"/",label,"/",sep=''), showWarnings=F, recursive=T)
  
  
  if(file.exists(paste(extraction.dir,"/",label,"/DEM_",res,".tif", sep='')) & !force.redo){
    extracted.DEM <- raster(paste(extraction.dir,"/",label,"/DEM_",res,".tif", sep=''))
    
  }else{
    extracted.DEM <- getNED(template, raw.dir=raw.dir, res=res)
    writeGDAL(as(extracted.DEM, "SpatialGridDataFrame"), paste(extraction.dir,"/",label,"/DEM_",res,".tif", sep=''), drivername="GTiff", type="Float32", options=c("INTERLEAVE=PIXEL", "COMPRESS=DEFLATE", "ZLEVEL=9"))
  }
  
  
  if(drain){
    if(file.exists(paste(extraction.dir,"/",label,"/DEM_",res,"_DRAINED.tif", sep='')) & !force.redo){
      dem.final <- raster(paste(extraction.dir,"/",label,"/DEM_",res,"_DRAINED.tif", sep=''))
      
    }else{
      if(is.null(NHD.raw.dir)){
        NHD.raw.dir <- readline("Please provide a path for the raw NHD data directory:")
      }
      if(is.null(NRCS.raw.dir)){
        NRCS.raw.dir <- readline("Please provide a path for the raw NRCS data directory:")
      }
      # Drain the DEM
      cat('Draining the DEM.\n\n')
      flush.console()
      dem.final <- drainDEM(template=template, dem=extracted.DEM, label=label, NED.raw.dir=raw.dir, NHD.raw.dir=NHD.raw.dir, NRCS.raw.dir=NRCS.raw.dir)
      writeGDAL(as(dem.final, "SpatialGridDataFrame"), paste(extraction.dir,"/",label,"/DEM_",res,"_DRAINED.tif", sep=''), drivername="GTiff", type="Float32", options=c("INTERLEAVE=PIXEL", "COMPRESS=DEFLATE", "ZLEVEL=9"))
    }
    
    return(dem.final)
  }else{
    return(extracted.DEM)
  }
}



# A function to automatically download and crop the National Elevation Dataset.
# x, a model raster, sp*, or extent object
# out.dir, a directory where the downloaded NED files should be stored
# res, the intended resolution as a character string. 
#     Options are "1" for the 1 arc-second NED (~30 meter), 
#     or "13" for the 1/3 arc-second NED (~10 meter)
# If x is a sp* or extent object, res must be provided!
getNED <- function(x, raw.dir="./Input/NED", res=NULL){  
  # Convert polygon to decimal degrees for data request
  extent.latlon <- extent(projectExtent(x, CRS("+proj=longlat")))
  
  if(is.null(res) & class(x) %in% c("RasterLayer","RasterStack","RasterBrick")){
    y.dim <- extent.latlon@ymax - extent.latlon@ymin
    y.res <- y.dim/nrow(x)
    if(y.res < (1/60)/60){
      res <- "13"
    }else{
      res <- "1"
    }
  }else if(is.null(res)){
    warning("If x is a sp* or extent object, res should be provided! \n Defaulting to 1 arc-second NED.\n", immediate.=T)
    res <- "1"
  }
  
  # Open USGS NED download service.
  # NED tiles are labeled by their northwest corner. Thus, coordinate 36.42N, -105.71W is in grid n37w106
  wests <- c(ceiling(abs(extent.latlon@xmax)),ceiling(abs(extent.latlon@xmin)))
  wests <- unique(append(wests,seq(wests[1],wests[2])))
  norths <- c(ceiling(abs(extent.latlon@ymin)),ceiling(abs(extent.latlon@ymax)))
  norths <- unique(append(norths,seq(norths[1],norths[2])))
  
  wests <- formatC(wests, width = 3, format = "d", flag = "0") 
  norths <- formatC(norths, width = 2, format = "d", flag = "0") 
  
  tiles <- vector("list", length(wests)*length(norths))
  
  cat("\nArea of interest includes",length(wests)*length(norths),"NED tiles.")
  
  # Automatically download 1x1 degree tiles from the USGS
  cat("\nChecking availability of 1x1 degree tiles from the USGS for the study area.\n")
  dir.create(paste(raw.dir,'/',res, sep=''), showWarnings = FALSE, recursive = TRUE)
  t <- 1
  for(w in wests){
    for(n in norths){
      if(!file.exists(paste(raw.dir,'/',res,'/n',n,'w',w,'.zip',sep=''))){
        if(url.exists(paste('ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/NED/',res,'/ArcGrid/n',n,'w',w,'.zip',sep=''))){
          cat(paste('\nDownloading: ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/NED/',res,'/ArcGrid/n',n,'w',w,'.zip\n',sep=''))
          flush.console()
          
          download.file(paste('ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/NED/',res,'/ArcGrid/n',n,'w',w,'.zip',sep=''), destfile=paste(raw.dir,'/',res,'/n',n,'w',w,'.zip',sep=''), mode="wb")
          
#           f = CFILE(paste(raw.dir,'/',res,'/n',n,'w',w,'.zip',sep=''), mode="wb")
#           curlPerform(url = paste('ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/NED/',res,'/ArcGrid/n',n,'w',w,'.zip',sep=''), writedata = f@ref)
#           close(f)
        }else{
          stop(paste('Please find a copy of the n',n,'w',w,'.zip NED grid in ESRI ARCGRID format, available from the USGS.',sep=''))
        }
      }
      unzip(paste(raw.dir,'/',res,'/n',n,'w',w,'.zip',sep=''),exdir=paste(raw.dir,'/',res,'/n',n,'w',w, sep=''))
      tiles[t] <- raster(readGDAL(paste(raw.dir,'/',res,'/n',n,'w',w,'/grdn',n,'w',w,'_',res,sep='')))
      unlink(paste(raw.dir,'/',res,'/n',n,'w',w,sep=''), recursive = TRUE)
      t <- t+1
    }
  }
  # Mosaic all tiles
  if(length(tiles)>1){
    cat('Mosaic-ing NED tiles.\n\n')
    flush.console()
    
    tiles$fun <- mean
    mosaic.all <- do.call(mosaic, tiles)
    
    rm(tiles)
    gc()
  }else{
    mosaic.all <- tiles[[1]]
  }
  
  # pre-crop
  mosaic.all <- crop(mosaic.all,extent.latlon, snap="out")
  return(mosaic.all)
}
