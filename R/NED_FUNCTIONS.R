# A function to automatically download and crop the National Elevation Dataset.
# x, a model raster, sp*, or extent object
# out.dir, a directory where the downloaded NED files should be stored
# res, the intended resolution as a character string. 
#     Options are "1" for the 1 arc-second NED (~30 meter), 
#     or "13" for the 1/3 arc-second NED (~10 meter)
# If template is a sp* or extent object, res must be provided!
# USES RCurl, raster, sp, rgdal
getNED <- function(template, label, res, raw.dir="./RAW/NED/", extraction.dir="./EXTRACTIONS/NED/", force.redo=F){  
  
  rasters.dir <- paste(extraction.dir,"/",label,"/rasters",sep='')
  
  dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(rasters.dir, showWarnings = FALSE, recursive = TRUE)
  
  if(file.exists(paste(rasters.dir,"/NED_",res,".tif", sep='')) & !force.redo){
    extracted.DEM <- raster::raster(paste(rasters.dir,"/NED_",res,".tif", sep=''))
    return(extracted.DEM)
  }
  
  # Convert polygon to decimal degrees for data request
  extent.latlon <- raster::extent(raster::projectExtent(template, sp::CRS("+proj=longlat +ellps=WGS84")))
  
  if(is.null(res) & (class(template) %in% c("RasterLayer","RasterStack","RasterBrick"))){
    y.dim <- extent.latlon@ymax - extent.latlon@ymin
    y.res <- y.dim/nrow(template)
    if(y.res < (1/60)/60){
      res <- "13"
    }else{
      res <- "1"
    }
  }else if(is.null(res)){
    warning("If template is a sp* or extent object, res should be provided! \n Defaulting to 1 arc-second NED.\n", immediate.=T)
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
      url <- paste('ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/NED/',res,'/ArcGrid/n',n,'w',w,'.zip',sep='')
      destdir <- paste(raw.dir,'/',res,'/',sep='')
      wgetDownload(url=url, destdir=destdir)

      unzip(paste(raw.dir,'/',res,'/n',n,'w',w,'.zip',sep=''),exdir=paste(raw.dir,'/',res,'/n',n,'w',w, sep=''))
      tiles[t] <- raster::raster(rgdal::readGDAL(paste(raw.dir,'/',res,'/n',n,'w',w,'/grdn',n,'w',w,'_',res,sep='')))
      unlink(paste(raw.dir,'/',res,'/n',n,'w',w,sep=''), recursive = TRUE)
      t <- t+1
    }
  }
  # Mosaic all tiles
  if(length(tiles)>1){
    cat('Mosaic-ing NED tiles.\n\n')
    flush.console()
    
    tiles$fun <- mean
    mosaic.all <- do.call(raster::mosaic, tiles)
    
    rm(tiles)
    gc()
  }else{
    mosaic.all <- tiles[[1]]
  }
  
  # Crop
  mosaic.all <- raster::crop(mosaic.all,sp::spTransform(template,sp::CRS(raster::projection(mosaic.all))), snap="out")
  
  rgdal::writeGDAL(as(mosaic.all, "SpatialGridDataFrame"), paste(rasters.dir,"/DEM_",res,".tif", sep=''), drivername="GTiff", type="Float32", options=c("INTERLEAVE=PIXEL", "COMPRESS=DEFLATE", "ZLEVEL=9"))
  
  return(mosaic.all)
}
