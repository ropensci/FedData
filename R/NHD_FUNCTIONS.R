#' Download and crop the National Hydrography Dataset.
#'
#' \code{get_nhd} returns a list of Spatial* objects extracted 
#' from the National Hydrography Dataset.
#' 
#' @param template A Raster* or Spatial* object to serve 
#' as a template for cropping.
#' @param label A character string naming the study area.
#' @param res A character string defining the resolution of the NHD to download.
#' Either "medium" (the default), or "high".
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing. Defaults to "./RAW/NHD/".
#' @param extraction.dir A character string indicating where the extracted and cropped NHD shapefiles should be put.
#' The directory will be created if missing. Defaults to "./EXTRACTIONS/NHD/".
#' @param force.redo If an extraction for this template and label already exists, should a new one be created?
#' @return A list of Spatial* objects extracted from the National Hydrography Dataset.
#' @export
get_nhd <- function(template, label, res='medium', raw.dir="./RAW/NHD/", extraction.dir="./EXTRACTIONS/NHD/", force.redo=FALSE){
  
  vectors.dir <- paste(extraction.dir,"/",label,"/spatial",sep='')
  
  dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(vectors.dir, showWarnings = FALSE, recursive = TRUE)
  
  if(!force.redo & length(list.files(vectors.dir))>0){
    files <- list.files(vectors.dir)
    files <- files[grepl("shp",files)]
    files <- files[!grepl("template",files)]
    files <- gsub(".shp","",files)
    files <- files[order(files)]
    
    shapes <- lapply(files,function(file){
      rgdal::readOGR(normalizePath(vectors.dir),file, verbose=F)
    })
    names(shapes) <- files
    return(shapes)
  }
  
  if(class(template) %in% c("RasterLayer","RasterStack","RasterBrick")){
    template <- spdf_from_polygon(sp::spTransform(polygon_from_extent(template),sp::CRS("+proj=longlat +ellps=GRS80")))
  }
  
  message("(Down)Loading the NHD HUC4 dataset.")
  HUC4 <- get_huc4(template=template, raw.dir=raw.dir)
  
  area.list <- formatC(HUC4$HUC4, width = 4, format = "d", flag = "0")
  
  # Get the spatial data for each area
  message("(Down)Loading the NHD subregion data.")
  subregionShapes <- lapply(area.list, function(area){
    return(get_nhd_subregion(template=template, area=area, res=res, raw.dir=raw.dir))
  })
  
  # Get all layer names
  layers <- unique(unlist(lapply(subregionShapes,names)))
  
  # Merge like datasets
  message("Merging all NHD data in study area.")
  allShapes <- lapply(layers,function(layer){
    shapes <- sapply(subregionShapes, '[[', layer)
    null.shapes <- sapply(shapes,is.null)
    shapes <- do.call("rbind", shapes[!null.shapes])
    if(is.null(shapes)) return(shapes)
    shapes <- raster::crop(shapes,sp::spTransform(template,sp::CRS(raster::projection(shapes))))
    if(is.null(shapes)) return(shapes)
    #     shapes <- spTransform(shapes,CRS(projection(template)))
    suppressWarnings(rgdal::writeOGR(shapes,vectors.dir,layer,"ESRI Shapefile", overwrite_layer=TRUE))
    return(shapes)
  })
  names(allShapes) <- layers
  
  # Remove null layers
  allShapes <- allShapes[!sapply(allShapes,is.null)]
    
  return(allShapes)
}

#' Download a zipped directory containing a shapefile of the HUC4 subregions of the NHD.
#'
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' @return A character string representing the full local path of the HUC4 zipped directory.
#' @export
download_huc4 <- function(raw.dir){
#   url <- "ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/WBD/FileGDB101/WBD_National.zip"
  url <- 'ftp://ftp.igsb.uiowa.edu/gis_library/USA/huc_04.zip'
  curl_download(url=url, destdir=raw.dir)
  return(normalizePath(paste(raw.dir,'huc_04.zip',sep='')))
#   return(normalizePath(paste(raw.dir,'WBD_National.zip',sep='')))
}


#' Download and crop a shapefile of the HUC4 
#' regions of the National Hydrography Dataset.
#'
#' \code{get_huc4} returns a \code{SpatialPolygonsDataFrame} of the HUC4 regions within
#' the specified \code{template}. If template is not provided, returns the entire HUC4 dataset.
#' 
#' @param template A Raster* or Spatial* object to serve 
#' as a template for cropping.
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing.
#' @return A \code{SpatialPolygonsDataFrame} of the HUC4 regions within
#' the specified \code{template}.
#' @export
get_huc4 <- function(template=NULL, raw.dir){
  tmpdir <- tempfile()
  if (!dir.create(tmpdir))
    stop("failed to create my temporary directory")
  
  huc4File <- download_huc4(raw.dir)
  
  utils::unzip(huc4File, exdir=tmpdir)
  
  HUC4 <- rgdal::readOGR(tmpdir, layer="huc_04", verbose=FALSE)
  
  HUC4@proj4string <- sp::CRS("+proj=utm +zone=15 +datum=NAD83 +ellps=WGS84")
  
  # Get a list of NHD subregions within the project study area
  if(!is.null(template)){
    if(class(template) %in% c("RasterLayer","RasterStack","RasterBrick")){
      template <- spdf_from_polygon(sp::spTransform(polygon_from_extent(template),sp::CRS("+proj=longlat +ellps=GRS80")))
    }
    
    HUC4 <- raster::crop(HUC4,sp::spTransform(template,sp::CRS(raster::projection(HUC4))))
  }
  
  unlink(tmpdir, recursive = TRUE)
  
  return(HUC4)
}



#' Download a zipped NHD HUC4 subregion.
#'
#' HUC4 subregion are specified by a unique character string, best obtained using the \code{\link{get_huc4}} function.
#' \code{download_nhd_subregion} returns the path to the downloaded zip file.
#'
#' @param area A 4-character string indicating the HUC4 NHD subregion to download.
#' @param res A character string defining the resolution of the NHD to download.
#' Either "medium" (the default), or "high".
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing.
#' @return A character string representing the full local path of the downloaded zip file.
#' @export
download_nhd_subregion <- function(area, res, raw.dir){
  if(res=="medium"){
    url <- paste('ftp://nhdftp.usgs.gov/DataSets/Staged/SubRegions/FileGDB/MediumResolution/NHDM',area,'_931v220.zip',sep='')
  }else{
    url <- paste('ftp://nhdftp.usgs.gov/DataSets/Staged/SubRegions/FileGDB/HighResolution/NHDH',area,'_931v220.zip',sep='')
  }
  
  destdir <- raw.dir
  curl_download(url=url, destdir=destdir)
  
  return(normalizePath(paste0(destdir,'/',basename(url))))
}

#' Download and crop data from a zipped HUC4 subregion
#' of the National Hydrography Dataset.
#'
#' \code{get_nhd_subregion} returns a list of \code{SpatialPolygonsDataFrame}s of the layers of the HUC4 subregion,
#'  within the specified \code{template}. If template is not provided, returns the entire HUC4 subregion.
#' 
#' @param template A Raster* or Spatial* object to serve 
#' as a template for cropping.
#' @param area A 4-character string indicating the HUC4 NHD subregion to download and crop.
#' @param res A character string defining the resolution of the NHD to download.
#' Either "medium" (the default), or "high".
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing.
#' @return A \code{SpatialPolygonsDataFrame} of the HUC4 regions within
#' the specified \code{template}.
#' @export
get_nhd_subregion <- function(template=NULL, area, res, raw.dir){
  tmpdir <- tempfile()
  if (!dir.create(tmpdir))
    stop("failed to create my temporary directory")
  
  file <- download_nhd_subregion(area=area, res=res, raw.dir=raw.dir)
  
  utils::unzip(file,exdir=tmpdir)
  
  # Get the path to the geodatabase
  dsn <- list.files(tmpdir, full.names=T)
  dsn <- dsn[grepl("gdb",dsn)]
  dsn <- normalizePath(dsn)
  
  # List all layers in the geodatabase
  layers <- rgdal::ogrListLayers(dsn)
  layers <- layers[grepl("NHD",layers)]
  
  # Get each layer in the geodatabase
  shapes <- lapply(layers,function(layer){
    tryCatch(suppressWarnings(rgdal::readOGR(dsn=dsn, layer=layer, verbose=F)),error=function(e) NULL)
  })
  names(shapes) <- layers
  
  # Rename the features to prepare for merging
  shapes <- lapply(shapes, function(shape){
    tryCatch(sp::spChFIDs(shape, paste(area,"_",shape$Permanent_Identifier,sep='')),error=function(e) NULL)
  })
  
  # Crop each feature to the template area
  if(!is.null(template)){
    if(class(template) %in% c("RasterLayer","RasterStack","RasterBrick")){
      template <- spdf_from_polygon(sp::spTransform(polygon_from_extent(template),sp::CRS("+proj=longlat +ellps=GRS80")))
    }
    
    shapes <- lapply(shapes,function(shape){
      tryCatch(raster::crop(shape,sp::spTransform(template,sp::CRS(raster::projection(shape)))),error=function(e) NULL)
    }) 
  }
  
  shapes <- shapes[!sapply(shapes,is.null)]
  
  unlink(tmpdir, recursive = TRUE)
  
  return(shapes)
}