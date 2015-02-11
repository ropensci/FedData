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
  wests <- seq(ceiling(abs(extent.latlon@xmax)),ceiling(abs(extent.latlon@xmin)))
  norths <- seq(ceiling(abs(extent.latlon@ymin)),ceiling(abs(extent.latlon@ymax)))
  wests <- formatC(wests, width = 3, format = "d", flag = "0") 
  norths <- formatC(norths, width = 2, format = "d", flag = "0") 
  
  tilesLocations <- as.matrix(expand.grid(norths,wests,stringsAsFactors = FALSE))

  cat("\nArea of interest includes",nrow(tilesLocations),"NED tiles.")
  
  # Automatically download 1x1 degree tiles from the USGS
  cat("\nChecking availability of 1x1 degree tiles from the USGS for the study area.\n")
  dir.create(paste(raw.dir,'/',res, sep=''), showWarnings = FALSE, recursive = TRUE)
  
  # Download tiles
  tileFiles <- apply(tilesLocations,1,function(loc){
    outFile <- downloadNED(res=res, tileNorthing=loc[1], tileWesting=loc[2], raw.dir=raw.dir)
    return(outFile)
  })
  
  # Extract tiles from zip directories, and get a list of rasters
  # Force the rasters into memory
  tiles <- lapply(tileFiles,unzipNED)

  # Crop all tiles
  tiles <- lapply(tiles, function(tile){
    tryCatch(raster::crop(tile,sp::spTransform(template,sp::CRS(raster::projection(tile))), snap="out"),error=function(e){return(NULL)})
    })
  
  tiles <- tiles[!sapply(tiles,is.null)]
  
  # Mosaic all tiles
  if(length(tiles)>1){
    cat('Mosaic-ing NED tiles.\n\n')
    flush.console()
    
    tiles$fun <- mean
    tiles <- do.call(raster::mosaic, tiles)

    gc()
  }else{
    tiles <- tiles[[1]]
  }
  
  raster::writeRaster(tiles, paste(rasters.dir,"/DEM_",res,".tif", sep=''), datatype="FLT4S", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"),overwrite=T,setStatistics=FALSE)
  
  return(tiles)
}

downloadNED <- function(res, tileNorthing, tileWesting, raw.dir){
  
  url <- paste('ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/NED/',res,'/ArcGrid/n',tileNorthing,'w',tileWesting,'.zip',sep='')
  destdir <- paste(raw.dir,'/',res,'/',sep='')
  wgetDownload(url=url, destdir=destdir)
  
  return(normalizePath(paste(destdir,'n',tileNorthing,'w',tileWesting,'.zip',sep='')))
}

unzipNED <- function(file){
  unzip(file,exdir="./temp")
  
  dirs <- list.dirs("./temp",full.names = TRUE,recursive=F)
  dirs <- dirs[grepl("grdn",dirs)]
  
  tile <- raster::raster(rgdal::readGDAL(dirs))
  unlink("./temp", recursive = TRUE)
  return(tile)
}