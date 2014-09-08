# A function that loads the National Hydrography Dataset for a provided study area defined by "x."
extractNHD <- function(template, label, raw.dir, extraction.dir=NULL, remove.modern=TRUE, force.redo=FALSE){
  if(is.null(extraction.dir)){
    extraction.dir <- paste(raw.dir,"/EXTRACTIONS",sep='')
  }
  
  dir.create(paste(extraction.dir,"/",label,"/",sep=''), showWarnings = FALSE, recursive=T)
  dsn.vectors <- paste(extraction.dir,"/",label,"/vectors",sep='')
  dir.create(dsn.vectors, showWarnings = FALSE, recursive = TRUE)
  
  template.poly.latlon <- SPDFfromPolygon(spTransform(polygonFromExtent(template),CRS("+proj=longlat +ellps=GRS80")))
  
  # Write a polygon of the study area
  suppressWarnings(writeOGR(template.poly.latlon, dsn.vectors, "sim_poly","ESRI Shapefile", overwrite_layer=TRUE))  
  
  # Set the layers to be extracted
  layers <- c("NHDFlowline","NHDWaterbody","NHDArea")
  
  HUC4 <- loadHUC4(x=template, raw.dir=raw.dir, dsn.vectors=dsn.vectors, force.redo=force.redo)
  
  area.list <- getAreaList(HUC4)
  
  shapes <- loadNHDSubregions(x=template, layers=layers, area.list=area.list, raw.dir=raw.dir, dsn.vectors=dsn.vectors, force.redo=force.redo)
  
  splitAndExport(shapes=shapes, layers=layers, dsn.vectors=dsn.vectors)
  
  template.poly.latlon <- SPDFfromPolygon(spTransform(polygonFromExtent(template),CRS(projection(HUC4))))
  suppressWarnings(writeOGR(template.poly.latlon, dsn.vectors, "sim_poly","ESRI Shapefile", overwrite_layer=TRUE))
  
  streams <- readOGR(dsn.vectors,"Streams", verbose=F)
#   cat("\nNHD processed for study area defined by x")
  return(streams)
}

# A method for loading the HUC4 shapefile from the
# National Hydrography Dataset.
loadHUC4 <- function(x, raw.dir, dsn.vectors, force.redo=FALSE){
  if(!("HUC4" %in% ogrListLayers(dsn.vectors)) | force.redo){
    HUC4 <- getHUC4(x, raw.dir=raw.dir, force.redo=force.redo)
    writeOGR(HUC4, dsn.vectors, "HUC4","ESRI Shapefile", overwrite_layer=TRUE)
  }else{
    HUC4 <- readOGR(normalizePath(dsn.vectors),"HUC4", verbose=F)
  }
  return(HUC4)
}

# A method for downloading and processing the HUC4 shapefile from
# the National Hydrography Dataset.
getHUC4 <- function(x, raw.dir, force.redo=F){  
  if(url.exists('ftp://ftp.igsb.uiowa.edu/gis_library/basins/huc_04.zip')){
    system(paste("wget -np -nd -N ftp://ftp.igsb.uiowa.edu/gis_library/USA/huc_04.zip --directory-prefix=",raw.dir, sep=''))
  }else{
    stop("Unable to download HUC4 shapefile!")
  }
  
  cat("Unzipping the NHD HUC4 dataset.\n")
  unzip(paste(raw.dir,"/huc_04.zip",sep=''),exdir=paste(raw.dir,"/huc_04",sep=''))
  
  cat("Loading the NHD HUC4 dataset.\n")
  HUC4 <- readOGR(normalizePath(paste(raw.dir,"/huc_04/",sep='')), layer="huc_04", verbose=FALSE)
  
  HUC4@proj4string <- CRS("+proj=utm +zone=15 +datum=NAD83 +ellps=WGS84")
  
  # Get a list of NHD subregions within the project study area
  HUC4 <- crop.to.studyArea(HUC4,spTransform(x,CRS(projection(HUC4))))
  
  return(HUC4)
}

# A method for downloading medium resolution regions from the National
# Hydrography Dataset, given a list of regions.
getNHDSubregions <- function(area.list, raw.dir){
  # Download the corresponding NHD subregion files.
  unzip.new <- FALSE
  for(area in area.list){
    if(file.exists(paste(raw.dir,"/NHDM",area,".gdb",sep=''))) next
    
    if(url.exists(paste('ftp://nhdftp.usgs.gov/DataSets/Staged/SubRegions/FileGDB/MediumResolution/NHDM',area,'_92v200.zip',sep=''))){
      system(paste('wget -np -nd -N ftp://nhdftp.usgs.gov/DataSets/Staged/SubRegions/FileGDB/MediumResolution/NHDM',area,'_92v200.zip --directory-prefix=',raw.dir, sep=''))
      unzip(paste(raw.dir,'/NHDM',area,'_92v200.zip', sep=''),exdir=paste(raw.dir,'/',sep=''))
      unzip.new <- TRUE
    }else{
      stop("USGS NHD data unavailable!")
    }
  }
  
  if(unzip.new){
    # Ask the user to upgrade the file geodatabases
    cat("\nUnfortunately, the NHD data are in the ESRI ArcGIS 9.3.1 file geodatabase format, which cannot be read by R.
        \nIt must therefore be converted to the ArcGIS 10.0 or greater format. 
        \nTo do so, open up ArcCatalog10 and navigate to the",raw.dir,"directory under your working directory. 
        Select the",raw.dir,"directory. Your newly downloaded geodatabases will appear in the 'Contents' window. 
        Right-click the first geodatabase and select 'Properties' from the drop-down menu. 
        At the bottom of the 'General' tab there will be an 'Upgrade Geodatabase' button. 
        Click it, and follow the on-screen instructions. Repeat for all remaining geodatabases.")
    readline("\nPlease type  <Return>\t once all GDBs are converted.")
  }
}

# A method for loading specified layers fro specified regions in the National Hydrology Dataset
loadNHDLayers <- function(x, layers, area.list, raw.dir, dsn.vectors, force.redo=FALSE){
  layers.missing <- layers[!(layers %in% ogrListLayers(dsn.vectors))]
  
  if(length(layers.missing)>0){
    # Load USGS NHD data.
    # Load the Flowline, Point, and Waterbody data from the subregions specified.
#     cat("\nLoading the NHD layer data.\n")
    shapes <- vector("list",length(layers.missing))
    for(i in 1:length(layers.missing)){
      shapes[[i]] <- vector("list",length(area.list))
      for(j in 1:length(area.list)){
        
        dsn <- normalizePath(paste(raw.dir,"/",area.list[j],"_shapes",sep=''))
        shapes[[i]][[j]] <- suppressWarnings(readOGR(dsn, layer=layers[i], verbose=F))
        
        # Change all spatial IDs to prepare from merging
        if(class(shapes[[i]][[j]]) %in% c("SpatialLinesDataFrame","SpatialPolygonsDataFrame")){
          shapes[[i]][[j]] <- spChFIDs(shapes[[i]][[j]], as.character(paste(layers[i],shapes[[i]][[j]]$Permanent_,sep='')))
        }
      }
    }
    
    # Merge like datasets
    for(i in 1:length(layers.missing)){
      shapes[[i]] <- do.call("rbind", shapes[[i]])
      
      shapes[[i]] <- crop.to.studyArea(shapes[[i]],spTransform(x,CRS(projection(shapes[[i]]))))
      
      suppressWarnings(writeOGR(shapes[[i]],dsn.vectors,layers.missing[i],"ESRI Shapefile", overwrite_layer=TRUE))
    }  
  }
  
  shapes <- vector("list",length(layers))
  for(i in 1:length(layers)){
    shapes[[i]] <- readOGR(normalizePath(dsn.vectors),layers[i], verbose=FALSE)
  }
  
  # Transform to CRS of x
  #   cat("Transforming all shapefiles to the CRS of x\n")
  #   shapes <- lapply(shapes,function(y){ spTransform(y, CRS(projection(x))) })
  
  #   # Crop to area of x
  #   cat("Cropping all shapefiles to the extent of x\n")
  #   shapes <- lapply(shapes,function(y){ crop.to.studyArea(y,spTransform(x,CRS(projection(y)))) })
  
  return(shapes)
}

# A method for removing modern features from the National Hydrography Dataset.
removeModernFeatures <- function(shapes, layers){
#   cat("Removing modern or unwanted features\n")
  # Only keep flowlines that represent Streams/Rivers (46006/46003), or 
  # "Artificial Path," which are generally historic river courses under modern reservoirs
  if(!is.na(match("NHDFlowline",layers))){
    index <- match("NHDFlowline",layers)
    shapes[[index]] <- shapes[[index]][shapes[[index]]$FCode %in% c(46003,46006,55800),]
  }
  
  # Only keep waterbodies that are Lakes/ponds (39004/39001/39009), and reservoirs (39010)
  if(!is.na(match("NHDWaterbody",layers))){
    index <- match("NHDWaterbody",layers)
    if(any(shapes[[index]]$FCode %in% c(39004,39001,39009,39010))){
      shapes[[index]] <- shapes[[index]][shapes[[index]]$FCode %in% c(39004,39001,39009,39010),]
    }else{
      shapes[[index]] <- "MISSING"
    }
  }
  return(shapes)
}

# A method for splitting features in the National Hydrography Dataset by their feature type.
splitAndExport <- function(shapes, layers, dsn.vectors){
  # Split shapefiles by feature types.
#   cat("\nJoining shapefiles by feature types.\n")
  index <- match("NHDFlowline",layers)
  if(class(shapes[[index]])=="SpatialLinesDataFrame"){
    if(any(shapes[[index]]$FCode %in% c(55800,46006,46003))){
      Streams <- shapes[[index]][shapes[[index]]$FCode %in% c(55800,46006,46003),]
      suppressWarnings(writeOGR(Streams,dsn.vectors,"Streams","ESRI Shapefile", overwrite_layer=TRUE))
    }
  }else{
    Streams <- NULL
  }

  index <- match("NHDWaterbody",layers)
  if(class(shapes[[index]])=="SpatialPolygonsDataFrame"){
    if(any(shapes[[index]]$FCode %in% c(39004,39009,39010,39001))){
      Lakes <- shapes[[index]][shapes[[index]]$FCode %in% c(39004,39009,39010,39001),]
      Lakes$ReachCode <- NULL
    }else{
      Lakes <- NULL
    }
  }else{
    Lakes <- NULL
  }

  
  # "Areas" are large Bureau of Reclamation reservoirs
  index <- match("NHDArea",layers)
  if(class(shapes[[index]])=="SpatialLinesDataFrame"){
    if(any(shapes[[index]]$FCode %in% c(40309))){
      Areas <- shapes[[index]][shapes[[index]]$FCode %in% c(40309),]
    }else{
      Areas <- NULL
    }
  }else{
    Areas <- NULL
  }
  
  if(!is.null(Lakes) & !is.null(Areas)){
    Lakes <- spChFIDs(Lakes, as.character(paste("LAKE_",row.names(as(Lakes, "data.frame")),sep='')))
    Areas <- spChFIDs(Areas, as.character(paste("AREA_",row.names(as(Areas, "data.frame")),sep='')))
    
    # Join waterbodies and reservoir layers
    Reservoirs <- rbind(Lakes,Areas)
  }else{
    if(!is.null(Lakes)){
      Reservoirs <- Lakes
    }else if(!is.null(Areas)){
      Reservoirs <- Areas
    }else{
      Reservoirs <- NULL
    }
  }
  
  if(!is.null(Reservoirs)){
    # Remove holes, and small lakes (< 400x400 meters)
    Reservoirs <- Reservoirs[Reservoirs$AreaSqKm>0.16,]
    for(i in 1:length(Reservoirs)){
      Reservoirs@polygons[[i]] <- remove.holes(Reservoirs@polygons[[i]])
    }
    
    # Remove water bodies that touch no streams
    Reservoirs <- spTransform(Reservoirs,CRS(projection(Streams)))
    #   Reservoirs <- Reservoirs[!is.na(c(Reservoirs %over% Streams)),]
    Reservoirs <- Reservoirs[!is.na(c(as(Reservoirs,"SpatialPolygons") %over% as(Streams,"SpatialLines"))),]
    
    # Export final vector datasets for the study area.
#     cat("\nExporting final vector and raster data.\n")
    suppressWarnings(writeOGR(Reservoirs,dsn.vectors,"Reservoirs","ESRI Shapefile", overwrite_layer=TRUE))
  }
  
}

# A method for loading subregions from the National Hydrography Dataset, given HUC4 areas and the required layers.
loadNHDSubregions <- function(x, layers, area.list, raw.dir, dsn.vectors, remove.modern=TRUE, force.redo=FALSE){
  getNHDSubregions(area.list, raw.dir=raw.dir)
  
  shapes <- loadNHDLayers(x, layers=layers, area.list=area.list, raw.dir=raw.dir, dsn.vectors=dsn.vectors, force.redo=force.redo)
  
  # Remove modern features
  if(remove.modern){
    shapes <- removeModernFeatures(shapes, layers)
  }
  
  return(shapes)
}

# Get a list of HUC4 NHD subregions
getAreaList <- function(HUC4){
  area.list <- HUC4$HUC4
  area.list <- formatC(area.list, width = 4, format = "d", flag = "0")
  return(area.list)
}

# A wrapper for the RGEOS function "gLineMerge" specifically to merge
# Stream segments of the National Hydrography Dataset
mergeStreams <- function(x, proj4string){
  Streams.final <- vector("list")
  for(i in 1:length(x)){
    temp.merge <- gLineMerge(SpatialLines(list(x[i,]@lines[[1]]),proj4string=CRS(proj4string)))
    for(j in 1:length(temp.merge@lines)){
      for(k in 1:length(temp.merge@lines[[j]]@Lines)){
        Streams.final <- c(Streams.final,temp.merge@lines[[j]]@Lines[[k]])
      }
    }
  }
  Streams.lines.list <- lapply(Streams.final,function(...){ Lines(...,ID="lines") })
  for(i in 1:length(Streams.lines.list)){
    Streams.lines.list[[i]]@ID <- paste("Stream",i)
  }
  Streams <- SpatialLines(Streams.lines.list,proj4string=CRS(proj4string))
  return(Streams)
}

# Strong-arm the joining of stream segments that remain unmerged by gLineMerge.
forceMergeStreams <- function(streams, proj4string){
  # Do one more loop of gLineMerge to allow RGEOS to attempt to merge more segments
  done=FALSE
  while(!done){
    start.length <- length(streams)
    streams <- gLineMerge(streams,byid=TRUE)
    Streams.lines.list <- lapply(streams@lines[[1]]@Lines,function(...){ Lines(...,ID="lines") } )
    for(i in 1:length(Streams.lines.list)){
      Streams.lines.list[[i]]@ID <- paste("Stream",i)
    }
    streams <- SpatialLines(Streams.lines.list,proj4string=proj4string)
    if(length(streams) == start.length){
      done=TRUE
    }
  }
  
  # RGEOS doesn't seem to do a complete job for some reason, 
  # so this loop strong-arms the joining of stream segments that remain.
  done <- FALSE
  while(!done){
    start.length <- length(streams)
    gapped.list <- vector("list")
    joined <- vector("logical", length(streams))
    counter <- 1
    for(i in 1:length(streams)){
      for(j in 1:length(streams)){
        if(joined[i]) break
        if(joined[j]) next
        if(i==j) next
        if(all(tail(streams@lines[[i]]@Lines[[1]]@coords, n=1) == head(streams@lines[[j]]@Lines[[1]]@coords, n=1))){
          gapped.list <- c(gapped.list,Lines(Line(rbind(streams@lines[[i]]@Lines[[1]]@coords,streams@lines[[j]]@Lines[[1]]@coords)),ID=as.character(counter)))
          counter <- counter+1
          joined[c(i,j)] <- TRUE
          break
        }
      }
    }
    
    for(i in 1:length(streams)){
      if(!joined[i]) gapped.list <- c(gapped.list,streams@lines[[i]])
    }
    
    for(i in 1:length(gapped.list)){
      gapped.list[[i]]@ID <- paste("Stream ",i,sep='')
    }
    
    streams <- SpatialLines(gapped.list,proj4string=proj4string)
    
    if(length(streams)==start.length) done <- TRUE
  }
  
  return(streams)
}