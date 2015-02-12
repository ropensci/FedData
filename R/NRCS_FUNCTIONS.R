#' Download and crop data from the SSURGO SSURGO soils database.
#'
#' This is an efficient method for spatially merging several different soil survey areas
#' as well as merging their tabular data.
#' 
#' \code{getSSURGO} returns a named list of length 2:
#' \enumerate{
#' \item "spatial": A \code{\link{SpatialPolygonsDataFrame}} of soil mapunits 
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
  SSURGOAreas <- getSSURGOStudyAreas(template=template, raw.dir=raw.dir)
  
  # Remove SSURGO study areas that are not availables
  SSURGOAreas <- SSURGOAreas[SSURGOAreas@data$iscomplete != 0,]
  
  # Get all mapunit polygons for the study area
  SSURGOMapunits <- getSSURGOMapunits(template=template, areas=SSURGOAreas, raw.dir=raw.dir)
  
  # Save the mapunit polygons as a shapefile
  suppressWarnings(rgdal::writeOGR(SSURGOMapunits, vectors.dir, "SSURGOMapunits","ESRI Shapefile", overwrite_layer=TRUE))
  
  # Get all of the tabular data
  SSURGOData <- getSSURGOData(areas=SSURGOAreas, raw.dir=raw.dir)
  
  # Extract only the mapunits in the study area, and iterate through the data structure
  SSURGOData <- extractSSURGOData(tables=SSURGOData, mapunits=SSURGOMapunits)
  
  # Save the each data table as a csv
  junk <- lapply(names(SSURGOData), function(tab){
    write.csv(SSURGOData[[tab]],file=paste(tables.dir,'/',tab,'.csv',sep=''),row.names=F)
  })
  
  return(list(spatial=SSURGOMapunits,tabular=SSURGOData))
}

getSSURGOStudyAreas <- function(template=NULL, raw.dir){
  # Import the shapefile of SSURGO study areas.
  # This is available at
  # http://soildatamart.sc.egov.usda.gov/download/StatusMaps/soilsa_a_SSURGO.zip
  url <- 'http://websoilsurvey.sc.egov.usda.gov/DataAvailability/SoilDataAvailabilityShapefile.zip'
  destdir <- raw.dir
  wgetDownload(url=url, destdir=destdir)
  
  cat("Unzipping the SSURGO study areas.\n")
  unzip(paste(raw.dir,"/SoilDataAvailabilityShapefile.zip",sep=''),exdir=paste(raw.dir,"/SoilDataAvailabilityShapefile",sep=''))
  
  cat("Loading the SSURGO study areas.\n")
  SSURGOAreas <- rgdal::readOGR(normalizePath(paste(raw.dir,"/SoilDataAvailabilityShapefile/",sep='')), layer="soilsa_a_SSURGO", verbose=FALSE)
  
  if(is.null(template)){
    return(SSURGOAreas)
  }
  
  # Get a list of NHD subregions within the project study area
  SSURGOAreas <- raster::crop(SSURGOAreas,sp::spTransform(template,sp::CRS(projection(SSURGOAreas))))
  
  unlink(paste(raw.dir,"/SoilDataAvailabilityShapefile",sep=''), recursive = TRUE)
  
  # Check to see if all survey areas are available
  if(0 %in% SSURGOAreas@data$iscomplete){
    cat("WARNING! Some of the soil surveys in your area are unavailable.\n")
    cat("Soils and productivity data will have holes.\n")
    cat("Missing areas:\n")
    cat(as.vector(SSURGOAreas@data[SSURGOAreas@data$iscomplete==0,]$areasymbol))
    cat("\n\n")
    cat("Continuing with processing available soils.\n\n")
  }
  
  return(SSURGOAreas)
}

getSSURGOStudyAreaMapunits <- function(area,date,raw.dir){  
  cat("(Down)Loading soils polygons for",area,'\n')
  
  state <- substring(area,1,2)
  
  # Try to download with the state database, otherwise grab the US
  url <- paste("http://websoilsurvey.sc.egov.usda.gov/DSD/Download/Cache/SSA/wss_SSA_",area,"_soildb_",state,"_2003_[",date,"].zip",sep='')
  destdir <- raw.dir
  tryCatch(wgetDownload(url=url, destdir=destdir, nc=T), warning = function(w) {
    url <- paste("http://websoilsurvey.sc.egov.usda.gov/DSD/Download/Cache/SSA/wss_SSA_",area,"_soildb_US_2003_[",date,"].zip",sep='')
    wgetDownload(url=url, destdir=destdir, nc=T)
    state <<- "US"
  })
  
  unzip(paste(raw.dir,'/wss_SSA_',area,'_soildb_',state,'_2003_[',date,'].zip', sep=''),exdir=raw.dir)
  
  mapunits <- rgdal::readOGR(paste(raw.dir,'/',area,'/spatial',sep=''), layer=paste("soilmu_a_",tolower(area),sep=''), verbose=F)
  
  unlink(paste(raw.dir,'/',area,sep=''), recursive=T, force=T)
  
  return(mapunits)
}

getSSURGOMapunits <- function(template, areas, raw.dir){
  # Load raw SSURGO Map Unit polygons from the regions specified.
  cat("Loading the SSURGO soil survey polygons for each region.\n")
  SSURGOPolys <- vector("list", length(areas))
  
  template <- sp::spTransform(template,CRS(projection(areas)))
  
  for(i in 1:length(areas)){
    area <- as.character(areas$areasymbol[i])
    date <- as.Date(areas$saverest[i])
    poly <- getSSURGOStudyAreaMapunits(area=area, date=date, raw.dir=raw.dir)
    poly <- poly[!is.na(poly %over% template),]
    poly <- sp::spChFIDs(poly, as.character(paste(area,'_',row.names(poly@data),sep='')))
    SSURGOPolys[[i]] <- poly
  }
  
  # Merging all SSURGO Map Unit polygons
  cat("Merging all SSURGO Map Unit polygons\n")
  SSURGOPolys <- do.call("rbind", SSURGOPolys)
  
  # Crop to area of x
  cat("Cropping SSURGO Map Unit polygons to the extent of the template\n")
  SSURGOPolys <- raster::crop(SSURGOPolys,template)
  
  return(SSURGOPolys)
}

getSSURGOStudyAreaData <- function(area,date,raw.dir){  
  cat("(Down)Loading soils data for",area,'\n')
  state <- substring(area,1,2)
  
  # Try to download with the state database, otherwise grab the US
  url <- paste("http://websoilsurvey.sc.egov.usda.gov/DSD/Download/Cache/SSA/wss_SSA_",area,"_soildb_",state,"_2003_[",date,"].zip",sep='')
  destdir <- raw.dir
  tryCatch(wgetDownload(url=url, destdir=destdir, nc=T), warning = function(w) {
    url <- paste("http://websoilsurvey.sc.egov.usda.gov/DSD/Download/Cache/SSA/wss_SSA_",area,"_soildb_US_2003_[",date,"].zip",sep='')
    wgetDownload(url=url, destdir=destdir, nc=T)
    state <<- "US"
  })
  
  unzip(paste(raw.dir,'/wss_SSA_',area,'_soildb_',state,'_2003_[',date,'].zip', sep=''),exdir=raw.dir)
  
  # Read in all tables
  files <- list.files(paste(raw.dir,area,'/tabular/',sep=''))
  tablesData <- lapply(files, function(file){
    tryCatch(return(read.delim(paste(raw.dir,area,'/tabular/',file,sep=''), header=F,sep="|")), error = function(e){return(NULL)})
  })
  names(tablesData) <- files
  tablesData <- tablesData[!sapply(tablesData,is.null)]
  
  tablesHeaders <- Hmisc::mdb.get(paste(raw.dir,area,"/soildb_",state,"_2003.mdb",sep=''))
  
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
  
  unlink(paste(raw.dir,'/',area,sep=''), recursive=T, force=T)
  
  return(tables)
}

getSSURGOData <- function(areas, raw.dir){
  # Load raw SSURGO Map Unit polygons from the regions specified.
  cat("Loading the SSURGO soil survey data for each region.\n")
  SSURGOTables <- vector("list", length(areas))
  
  for(i in 1:length(areas)){
    area <- as.character(areas$areasymbol[i])
    date <- as.Date(areas$saverest[i])
    tables <- getSSURGOStudyAreaData(area=area, date=date, raw.dir=raw.dir)
    SSURGOTables[[i]] <- tables
  }
  
  # Merging all SSURGO data tables
  cat("Merging all SSURGO data tables\n")
  tableNames <- unique(unlist(sapply(SSURGOTables,names)))
  tableNames <- tableNames[order(tableNames)]
  
  SSURGOTables <- lapply(tableNames,function(name){
    tables <- lapply(SSURGOTables,'[[',name)
    tables <- do.call("rbind", tables)
    tables <- unique(tables)
    return(tables)
  })
  
  names(SSURGOTables) <- tableNames
  
  return(SSURGOTables)
}

extractSSURGOData <- function(tables,mapunits){
  mapunits <- as.character(unique(mapunits$MUKEY))
  
  mapping <- tables[['mdstatrshipdet']]
  mappingGraph <- igraph::graph.edgelist(as.matrix(mapping[,c("ltabphyname","rtabphyname")]))
  igraph::E(mappingGraph)$mapVar <- as.character(mapping[,'ltabcolphyname'])
  
  mappingGraph <- igraph::graph.neighborhood(mappingGraph,order=max(sapply(igraph::decompose.graph(mappingGraph),diameter))+1,nodes='mapunit', mode='out')[[1]]
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