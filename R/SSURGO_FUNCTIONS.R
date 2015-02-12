#' Download and crop data from the NRCS SSURGO soils database.
#'
#' This is an efficient method for spatially merging several different soil survey areas
#' as well as merging their tabular data.
#' 
#' \code{getSSURGO} returns a named list of length 2:
#' \enumerate{
#' \item "spatial": A \code{SpatialPolygonsDataFrame} of soil mapunits 
#' in the template, and 
#' \item "tabular": A named list of \code{\link{data.frame}s} with the SSURGO tabular data.
#' }
#' 
#' @param template A Raster* or Spatial* object to serve 
#' as a template for cropping.
#' @param label A character string naming the study area.
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing. Defaults to "./RAW/SSURGO/".
#' @param extraction.dir A character string indicating where the extracted and cropped SSURGO shapefiles should be put.
#' The directory will be created if missing. Defaults to "./EXTRACTIONS/SSURGO/".
#' @param force.redo If an extraction for this template and label already exists, should a new one be created? Defaults to FALSE.
#' @return A named list containing the "spatial" and "tabular" data.
getSSURGO <- function(template, label, raw.dir="./RAW/SSURGO/", extraction.dir="./EXTRACTIONS/SSURGO/", force.redo=FALSE){  
  vectors.dir <- paste(extraction.dir,"/",label,"/spatial",sep='')
  tables.dir <- paste(extraction.dir,"/",label,"/tabular",sep='')
  
  dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(vectors.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(tables.dir, showWarnings = FALSE, recursive = TRUE)
  
  if(!force.redo & length(list.files(vectors.dir))>0 & length(list.files(tables.dir))>0 & file.exists(paste(vectors.dir,"/SSURGOMapunits.shp",sep=''))){
    SSURGOMapunits <- rgdal::readOGR(normalizePath(vectors.dir),"SSURGOMapunits", verbose=F)
    
    files <- list.files(tables.dir)
    files <- files[grepl("csv",files)]
    files <- files[order(files)]
    
    tables <- lapply(files,function(file){
      read.csv(paste(normalizePath(tables.dir),'/',file,sep=''))
    })
    names(tables) <- files
    
    return(list(spatial=SSURGOMapunits,tabular=tables))
  }
  
  if(class(template) %in% c("RasterLayer","RasterStack","RasterBrick")){
    template <- SPDFfromPolygon(sp::spTransform(polygonFromExtent(template),sp::CRS("+proj=longlat +ellps=GRS80")))
  }
  
  # Get shapefile of SSURGO study areas in the template
  SSURGOAreas <- getSSURGOInventory(template=template, raw.dir=raw.dir)
  
  # Remove SSURGO study areas that are not available
  SSURGOAreas <- SSURGOAreas[SSURGOAreas@data$iscomplete != 0,]
  
  # Get data for each study area
  SSURGOData <- lapply(1:length(SSURGOAreas), function(i){
    cat("\n(Down)Loading SSURGO data for subregion",i,"of",length(SSURGOAreas))
    getSSURGOStudyArea(template=template, area=as.character(SSURGOAreas$areasymbol[i]), date=as.Date(SSURGOAreas$saverest[i]), raw.dir=raw.dir)
  })
  
  # Combine mapunits
  SSURGOPolys <- lapply(SSURGOData,"[[","spatial")
  
  # Merging all SSURGO Map Unit polygons
  cat("\nMerging all SSURGO Map Unit polygons")
  SSURGOPolys <- do.call("rbind", SSURGOPolys)
  
  # Crop to area of template
  cat("\nCropping all SSURGO Map Unit polygons to area of template")
  SSURGOPolys <- raster::crop(SSURGOPolys,sp::spTransform(template,sp::CRS(raster::projection(SSURGOPolys))))
  
  
  # Combine subregion data
  SSURGOTables <- lapply(SSURGOData,"[[","tabular")
  
  # Merging all SSURGO data tables
  cat("\nMerging all SSURGO data tables")
  tableNames <- unique(unlist(sapply(SSURGOTables,names)))
  tableNames <- tableNames[order(tableNames)]
  
  SSURGOTables <- lapply(tableNames,function(name){
    tables <- lapply(SSURGOTables,'[[',name)
    tables <- do.call("rbind", tables)
    tables <- unique(tables)
    return(tables)
  })
  
  names(SSURGOTables) <- tableNames
  
  
  # Extract only the mapunits in the study area, and iterate through the data structure
  SSURGOTables <- extractSSURGOData(tables=SSURGOTables, mapunits=SSURGOPolys)
  
  # Save the mapunit polygons
  suppressWarnings(rgdal::writeOGR(SSURGOPolys,vectors.dir,"SSURGOMapunits","ESRI Shapefile", overwrite_layer=TRUE))
  
  # Save the each data table as a csv
  junk <- lapply(names(SSURGOTables), function(tab){
    write.csv(SSURGOTables[[tab]],file=paste(tables.dir,'/',tab,'.csv',sep=''),row.names=F)
  })
  
  return(list(spatial=SSURGOPolys,tabular=SSURGOTables))
}

#' Download a zipped directory containing a shapefile of the SSURGO study areas.
#'
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' @return A character string representing the full local path of the SSURGO study areas zipped directory.
downloadSSURGOInventory <- function(raw.dir){
  # Import the shapefile of SSURGO study areas.
  # This is available at
  # http://soildatamart.sc.egov.usda.gov/download/StatusMaps/soilsa_a_SSURGO.zip
  url <- 'http://websoilsurvey.sc.egov.usda.gov/DataAvailability/SoilDataAvailabilityShapefile.zip'
  destdir <- raw.dir
  wgetDownload(url=url, destdir=destdir)
  return(normalizePath(paste(destdir,'SoilDataAvailabilityShapefile.zip',sep='')))
}

#' Download and crop a shapefile of the SSURGO study areas.
#'
#' \code{getSSURGOInventory} returns a \code{SpatialPolygonsDataFrame} of the SSURGO study areas within
#' the specified \code{template}. If template is not provided, returns the entire SSURGO inventory of study areas.
#' 
#' @param template A Raster* or Spatial* object to serve 
#' as a template for cropping.
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing.
#' @return A \code{SpatialPolygonsDataFrame} of the SSURGO study areas within
#' the specified \code{template}.
getSSURGOInventory <- function(template=NULL, raw.dir){
  tmpdir <- tempfile()
  if (!dir.create(tmpdir))
    stop("failed to create my temporary directory")
  
  file <- downloadSSURGOInventory(raw.dir=raw.dir)
  
  unzip(file,exdir=tmpdir)

  SSURGOAreas <- rgdal::readOGR(normalizePath(tmpdir), layer="soilsa_a_nrcs", verbose=FALSE)
  
  unlink(tmpdir, recursive = TRUE)
  
  if(!is.null(template)){
    
    if(class(template) %in% c("RasterLayer","RasterStack","RasterBrick")){
      template <- SPDFfromPolygon(sp::spTransform(polygonFromExtent(template),sp::CRS("+proj=longlat +ellps=GRS80")))
    }
    
    # Get a list of SSURGO study areas within the project study area
    SSURGOAreas <- raster::crop(SSURGOAreas,sp::spTransform(template,sp::CRS(raster::projection(SSURGOAreas))))
    
    # Check to see if all survey areas are available
    if(0 %in% SSURGOAreas@data$iscomplete){
      cat("WARNING! Some of the soil surveys in your area are unavailable.\n")
      cat("Soils and productivity data will have holes.\n")
      cat("Missing areas:\n")
      cat(as.vector(SSURGOAreas@data[SSURGOAreas@data$iscomplete==0,]$areasymbol))
      cat("\n\n")
      cat("Continuing with processing available soils.\n\n")
    }
  }
  
  return(SSURGOAreas)
}

#' Download a zipped directory containing the spatial and tabular data for a SSURGO study area.
#'
#' \code{downloadSSURGOStudyArea} first tries to download data including a state-specific Access
#' template, then the general US template.
#'
#' @param area A character string indicating the SSURGO study area to be downloaded.
#' @param date A character string indicating the date of the most recent update to the SSURGO 
#' area for these data. This information may be gleaned from the SSURGO Inventory (\code{\link{getSSURGOInventory}}).
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' @return A character string representing the full local path of the SSURGO study areas zipped directory.
downloadSSURGOStudyArea <- function(area, date, raw.dir){
  state <- substring(area,1,2)
  
  # Try to download with the state database, otherwise grab the US
  url <- paste("http://websoilsurvey.sc.egov.usda.gov/DSD/Download/Cache/SSA/wss_SSA_",area,"_soildb_",state,"_2003_[",date,"].zip",sep='')
  destdir <- raw.dir
  tryCatch(wgetDownload(url=url, destdir=destdir,nc=T,timestamping=F), warning = function(w) {
    url <- paste("http://websoilsurvey.sc.egov.usda.gov/DSD/Download/Cache/SSA/wss_SSA_",area,"_soildb_US_2003_[",date,"].zip",sep='')
    wgetDownload(url=url, destdir=destdir,nc=T,timestamping=F)
    state <<- "US"
  })
  
  return(normalizePath(paste(destdir,"wss_SSA_",area,"_soildb_",state,"_2003_[",date,"].zip",sep='')))
}

#' Download and crop the spatial and tabular data for a SSURGO study area.
#'
#' \code{getSSURGOStudyArea} returns a named list of length 2:
#' \enumerate{
#' \item "spatial": A \code{SpatialPolygonsDataFrame} of soil mapunits 
#' in the template, and 
#' \item "tabular": A named list of \code{\link{data.frame}s} with the SSURGO tabular data.
#' }
#' 
#' @param template A Raster* or Spatial* object to serve 
#' as a template for cropping. If missing, whose study area is returned
#' @param area A character string indicating the SSURGO study area to be downloaded.
#' @param date A character string indicating the date of the most recent update to the SSURGO 
#' area for these data. This information may be gleaned from the SSURGO Inventory (\code{\link{getSSURGOInventory}}).
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing.
#' @return A \code{SpatialPolygonsDataFrame} of the SSURGO study areas within
#' the specified \code{template}.
getSSURGOStudyArea <- function(template=NULL, area, date, raw.dir){
  tmpdir <- tempfile()
  if (!dir.create(tmpdir))
    stop("failed to create my temporary directory")
  
  file <- downloadSSURGOStudyArea(area=area, date=date, raw.dir=raw.dir)
  
  unzip(file,exdir=tmpdir)
  
  # Get spatial data
  mapunits <- rgdal::readOGR(paste(tmpdir,'/',area,'/spatial',sep=''), layer=paste("soilmu_a_",tolower(area),sep=''), verbose=F)
  
  # Crop to study area
  if(!is.null(template)){
    if(class(template) %in% c("RasterLayer","RasterStack","RasterBrick")){
      template <- SPDFfromPolygon(sp::spTransform(polygonFromExtent(template),sp::CRS("+proj=longlat +ellps=GRS80")))
    }
    
    mapunits <- raster::crop(mapunits,sp::spTransform(template,sp::CRS(raster::projection(mapunits))))    
  }

  # Change IDs, in case of merging later
  mapunits <- sp::spChFIDs(mapunits, as.character(paste(area,'_',row.names(mapunits@data),sep='')))
  
  # Read in all tables
  files <- list.files(paste(tmpdir,'/',area,'/tabular/',sep=''))
  tablesData <- lapply(files, function(file){
    tryCatch(return(read.delim(paste(tmpdir,'/',area,'/tabular/',file,sep=''), header=F,sep="|")), error = function(e){return(NULL)})
  })
  names(tablesData) <- files
  tablesData <- tablesData[!sapply(tablesData,is.null)]
  
  dbFile <- list.files(paste(tmpdir,'/',area,sep=''),full.names=T)
  dbFile <- dbFile[grepl("mdb",dbFile)]
  tablesHeaders <- Hmisc::mdb.get(dbFile)
  
  SSURGOTableMapping <- tablesData[["mstab.txt"]][,c(1,5)]
  names(SSURGOTableMapping) <- c("TABLE","FILE")
  SSURGOTableMapping[,"FILE"] <- paste(SSURGOTableMapping[,"FILE"],'.txt',sep='')
  
  tablesData <- tablesData[as.character(SSURGOTableMapping[,"FILE"])]
  tablesHeaders <- tablesHeaders[as.character(SSURGOTableMapping[,"TABLE"])]
  
  notNull <- (!sapply(tablesData,is.null) & !sapply(tablesHeaders,is.null))
  tablesData <- tablesData[notNull]
  tablesHeaders <- tablesHeaders[notNull]
  
  tables <- mapply(tablesData,tablesHeaders,FUN=function(theData,theHeader){
    names(theData) <- names(theHeader)
    return(theData)
  })
  
  names(tables) <- names(tablesHeaders)
  
  tables <- extractSSURGOData(tables=tables, mapunits=mapunits)
  
  unlink(tmpdir, recursive = TRUE)
  
  return(list(spatial=mapunits,tabular=tables))
}

#' Extract data from a SSURGO databse pertaining to a set of mapunits.
#'
#' \code{extractSSURGOData} creates a directed graph of the joins in a SSURGO tabular dataset,
#' and then iterates through the tables, only retaining data pertinant to a set of mapunits.
#' 
#' @param tables A list of SSURGO tabular data.
#' @param mapunits A \code{SpatialPolygonsDataFrame} of mapunits (likely dropped from SSURGO spatial data)
#' defining which mapunits to retain.
#' @return A list of extracted SSURGO tabular data.
extractSSURGOData <- function(tables,mapunits){
  mapunits <- as.character(unique(mapunits$MUKEY))
  
  mapping <- tables[['mdstatrshipdet']]
  mappingGraph <- igraph::graph.edgelist(as.matrix(mapping[,c("ltabphyname","rtabphyname")]))
  igraph::E(mappingGraph)$mapVar <- as.character(mapping[,'ltabcolphyname'])
  
  mappingGraph <- igraph::graph.neighborhood(mappingGraph,order=max(sapply(igraph::decompose.graph(mappingGraph),igraph::diameter))+1,nodes='mapunit', mode='out')[[1]]
  mapHierarchy <- igraph::shortest.paths(mappingGraph,'mapunit')
  mapHierarchy <- colnames(mapHierarchy)[order(mapHierarchy)]
  mapHierarchy <- mapHierarchy[-1]
  mapEdges <- cbind(igraph::get.edgelist(mappingGraph),igraph::E(mappingGraph)$mapVar)
  mapEdges <- mapEdges[match(mapHierarchy,mapEdges[,2]),]
  
  tables[['mapunit']] <- tables[['mapunit']][tables[['mapunit']][,'mukey'] %in% mapunits,]
  newTables <- apply(mapEdges,1,function(X){
    return(tables[[X[2]]][tables[[X[2]]][,X[3]] %in% tables[[X[1]]][,X[3]],])
  })
  names(newTables) <- mapEdges[,2]
  
  tables[names(newTables)] <- newTables
  
  tables <- tables[!sapply(tables,is.null)]
  
  return(tables)
}