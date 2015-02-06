# A function that loads the National Hydrography Dataset for a provided study area defined by "x."
getNHD <- function(template, label, raw.dir, extraction.dir=NULL, force.redo=FALSE){
  if(is.null(extraction.dir)){
    extraction.dir <- paste(raw.dir,"/EXTRACTIONS",sep='')
  }
  
  dir.create(paste(extraction.dir,"/",label,"/",sep=''), showWarnings = FALSE, recursive=T)
  dsn.vectors <- paste(extraction.dir,"/",label,"/vectors",sep='')
  dir.create(dsn.vectors, showWarnings = FALSE, recursive = TRUE)
  
  if(!force.redo & length(list.files(dsn.vectors))>0){
    files <- list.files(dsn.vectors)
    files <- files[grepl("shp",files)]
    files <- files[!grepl("template",files)]
    files <- gsub(".shp","",files)
    files <- files[order(files)]
    
    shapes <- lapply(files,function(file){
      readOGR(normalizePath(dsn.vectors),file, verbose=F)
    })
    names(shapes) <- files
    return(shapes)
  }
  
  if(class(template) %in% c("RasterLayer","RasterStack","RasterBrick")){
    template <- SPDFfromPolygon(spTransform(polygonFromExtent(template),CRS("+proj=longlat +ellps=GRS80")))
  }
    
  HUC4 <- getHUC4(template=template, raw.dir=raw.dir, force.redo=force.redo)
  
  area.list <- getAreaList(HUC4)
  
  shapes <- loadNHDSubregions(template=template, area.list=area.list, raw.dir=raw.dir, dsn.vectors=dsn.vectors, force.redo=force.redo)
  
  return(shapes)
}

# A method for downloading and processing the HUC4 shapefile from
# the National Hydrography Dataset.
getHUC4 <- function(template, raw.dir, force.redo=F){
  url <- 'ftp://ftp.igsb.uiowa.edu/gis_library/USA/huc_04.zip'
  destdir <- raw.dir
  wgetDownload(url=url, destdir=destdir)
  
  cat("Unzipping the NHD HUC4 dataset.\n")
  unzip(paste(raw.dir,"/huc_04.zip",sep=''),exdir=paste(raw.dir,"/huc_04",sep=''))
  
  cat("Loading the NHD HUC4 dataset.\n")
  HUC4 <- readOGR(normalizePath(paste(raw.dir,"/huc_04/",sep='')), layer="huc_04", verbose=FALSE)
  
  HUC4@proj4string <- CRS("+proj=utm +zone=15 +datum=NAD83 +ellps=WGS84")
  
  # Get a list of NHD subregions within the project study area
  HUC4 <- raster::crop(HUC4,spTransform(template,CRS(projection(HUC4))))
  
  unlink(paste(raw.dir,"/huc_04",sep=''), recursive = TRUE)
  
  return(HUC4)
}

# Get a list of HUC4 NHD subregions
getAreaList <- function(HUC4){
  area.list <- HUC4$HUC4
  area.list <- formatC(area.list, width = 4, format = "d", flag = "0")
  return(area.list)
}

# A method for loading subregions from the National Hydrography Dataset, given HUC4 areas.
loadNHDSubregions <- function(template, area.list, raw.dir, dsn.vectors, force.redo=FALSE){
  getNHDSubregions(area.list, raw.dir=raw.dir)
  
  shapes <- loadNHDLayers(template=template, area.list=area.list, raw.dir=raw.dir, dsn.vectors=dsn.vectors, force.redo=force.redo)
  
  return(shapes)
}

# A method for downloading medium resolution regions from the National
# Hydrography Dataset, given a list of regions.
getNHDSubregions <- function(area.list, raw.dir){
  # Download the corresponding NHD subregion files.
  for(area in area.list){
    url <- paste('ftp://nhdftp.usgs.gov/DataSets/Staged/SubRegions/FileGDB/MediumResolution/NHDM',area,'_92v200.zip',sep='')
    destdir <- raw.dir
    wgetDownload(url=url, destdir=destdir)
  }
}

# A method for loading layers from specified regions in the National Hydrology Dataset
loadNHDLayers <- function(template, area.list, raw.dir, dsn.vectors, force.redo=FALSE){
  # Get the spatial data for each area
  all.data <- lapply(area.list, function(area){
    unzip(paste(raw.dir,'/NHDM',area,'_92v200.zip',sep=''),exdir=paste(raw.dir, sep=''))
    
    # Get the path to the geodatabase
    dsn <- normalizePath(paste(raw.dir,"/NHDM",area,".gdb",sep=''))
    
    # List all layers in the geodatabase
    layers <- ogrListLayers(dsn)
    
    # Get each layer in the geodatabase
    shapes <- lapply(layers,function(layer){
      tryCatch(suppressWarnings(readOGR(dsn=dsn, layer=layer, verbose=F)),error=function(e) NULL)
    })
    names(shapes) <- layers
    
    # Rename the features to prepare for merging
    shapes <- lapply(shapes, function(shape){
      tryCatch(spChFIDs(shape, as.character(paste(area,"_",row.names(shape@data),sep=''))),error=function(e) NULL)
    })
    
    unlink(paste(raw.dir,"/NHDM",area,".gdb",sep=''), recursive = TRUE)
    
    return(shapes)
  })
  
  # Get all layer names
  layers <- unique(unlist(lapply(all.data,names)))
  
  # Merge like datasets
  merged.data <- lapply(layers,function(layer){
    shapes <- sapply(all.data, '[[', layer)
    null.shapes <- sapply(shapes,is.null)
    shapes <- do.call("rbind", shapes[!null.shapes])
    if(is.null(shapes)) return(shapes)
    shapes <- raster::crop(shapes,spTransform(template,CRS(projection(shapes))))
    if(is.null(shapes)) return(shapes)
#     shapes <- spTransform(shapes,CRS(projection(template)))
    suppressWarnings(writeOGR(shapes,dsn.vectors,layer,"ESRI Shapefile", overwrite_layer=TRUE))
    return(shapes)
  })
  names(merged.data) <- layers
  null.layers <- sapply(merged.data,is.null)
  merged.data <- merged.data[!null.layers]
  
  return(merged.data)
}