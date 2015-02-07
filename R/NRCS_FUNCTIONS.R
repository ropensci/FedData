## All data is downloaded from the NRCS,
## available at http://websoilsurvey.sc.egov.usda.gov

## Author: R. Kyle Bocinsky
## Date: 02/14/2014

getNRCS <- function(template, label, raw.dir="./RAW/NRCS/", extraction.dir="./EXTRACTIONS/NRCS/", force.redo=FALSE){  
  vectors.dir <- paste(extraction.dir,"/",label,"/vectors",sep='')
  tables.dir <- paste(extraction.dir,"/",label,"/tables",sep='')
  
  dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(vectors.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(tables.dir, showWarnings = FALSE, recursive = TRUE)
  
  if(!force.redo & length(list.files(vectors.dir))>0 & length(list.files(tables.dir))>0){
    if(!file.exists(paste(vectors.dir,"/NRCSMapunits.shp",sep=''))) break
    NRCSMapunits <- rgdal::readOGR(normalizePath(vectors.dir),"NRCSMapunits", verbose=F)
    
    files <- list.files(tables.dir)
    files <- files[grepl("csv",files)]
    files <- files[order(files)]
    
    tables <- lapply(files,function(file){
      read.csv(paste(normalizePath(vectors.dir),file,'.csv',sep=''))
    })
    names(tables) <- files
    
    return(list(spatial=shapes,tabular=tables))
  }
  
  if(class(template) %in% c("RasterLayer","RasterStack","RasterBrick")){
    template <- SPDFfromPolygon(sp::spTransform(polygonFromExtent(template),sp::CRS("+proj=longlat +ellps=GRS80")))
  }
  
  # Get shapefile of NRCS study areas in the template
  NRCSAreas <- getNRCSStudyAreas(template=template, raw.dir=raw.dir)
  
  # Remove NRCS study areas that are not availables
  NRCSAreas <- NRCSAreas[NRCSAreas@data$iscomplete != 0,]
  
  # Get all mapunit polygons for the study area
  NRCSMapunits <- getNRCSMapunits(template=template, areas=NRCSAreas, raw.dir=raw.dir)
  
  # Save the mapunit polygons as a shapefile
  suppressWarnings(rgdal::writeOGR(NRCSMapunits, vectors.dir, "NRCSMapunits","ESRI Shapefile", overwrite_layer=TRUE))
  
  # Get all of the tabular data
  NRCSData <- getNRCSData(areas=NRCSAreas, raw.dir=raw.dir)
  
  # Extract only the mapunits in the study area, and iterate through the data structure
  NRCSData <- extractNRCSData(tables=NRCSData, mapunits=NRCSMapunits)
  
  # Save the each data table as a csv
  junk <- lapply(names(NRCSData), function(tab){
    write.csv(NRCSData[[tab]],file=paste(tables.dir,'/',tab,'.csv',sep=''),row.names=F)
  })
  
  return(list(spatial=NRCSMapunits,tabular=NRCSData))
}

getNRCSStudyAreas <- function(template=NULL, raw.dir){
  # Import the shapefile of NRCS study areas.
  # This is available at
  # http://soildatamart.sc.egov.usda.gov/download/StatusMaps/soilsa_a_nrcs.zip
  url <- 'http://websoilsurvey.sc.egov.usda.gov/DataAvailability/SoilDataAvailabilityShapefile.zip'
  destdir <- raw.dir
  wgetDownload(url=url, destdir=destdir)
  
  cat("Unzipping the NRCS study areas.\n")
  unzip(paste(raw.dir,"/SoilDataAvailabilityShapefile.zip",sep=''),exdir=paste(raw.dir,"/SoilDataAvailabilityShapefile",sep=''))
  
  cat("Loading the NRCS study areas.\n")
  NRCSAreas <- rgdal::readOGR(normalizePath(paste(raw.dir,"/SoilDataAvailabilityShapefile/",sep='')), layer="soilsa_a_nrcs", verbose=FALSE)
  
  if(is.null(template)){
    return(NRCSAreas)
  }
  
  # Get a list of NHD subregions within the project study area
  NRCSAreas <- raster::crop(NRCSAreas,sp::spTransform(template,sp::CRS(projection(NRCSAreas))))
  
  unlink(paste(raw.dir,"/SoilDataAvailabilityShapefile",sep=''), recursive = TRUE)
  
  # Check to see if all survey areas are available
  if(0 %in% NRCSAreas@data$iscomplete){
    cat("WARNING! Some of the soil surveys in your area are unavailable.\n")
    cat("Soils and productivity data will have holes.\n")
    cat("Missing areas:\n")
    cat(as.vector(NRCSAreas@data[NRCSAreas@data$iscomplete==0,]$areasymbol))
    cat("\n\n")
    cat("Continuing with processing available soils.\n\n")
  }
  
  return(NRCSAreas)
}

getNRCSStudyAreaMapunits <- function(area,date,raw.dir){  
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

getNRCSMapunits <- function(template, areas, raw.dir){
  # Load raw NRCS Map Unit polygons from the regions specified.
  cat("Loading the NRCS soil survey polygons for each region.\n")
  NRCSPolys <- vector("list", length(areas))
  
  template <- sp::spTransform(template,CRS(projection(areas)))
  
  for(i in 1:length(areas)){
    area <- as.character(areas$areasymbol[i])
    date <- as.Date(areas$saverest[i])
    poly <- getNRCSStudyAreaMapunits(area=area, date=date, raw.dir=raw.dir)
    poly <- poly[!is.na(poly %over% template),]
    poly <- sp::spChFIDs(poly, as.character(paste(area,'_',row.names(poly@data),sep='')))
    NRCSPolys[[i]] <- poly
  }
  
  # Merging all NRCS Map Unit polygons
  cat("Merging all NRCS Map Unit polygons\n")
  NRCSPolys <- do.call("rbind", NRCSPolys)
  
  # Crop to area of x
  cat("Cropping NRCS Map Unit polygons to the extent of the template\n")
  NRCSPolys <- raster::crop(NRCSPolys,template)
  
  return(NRCSPolys)
}

getNRCSStudyAreaData <- function(area,date,raw.dir){  
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
  
  NRCSTableMapping <- tablesData[["mstab.txt"]][,c(1,5)]
  names(NRCSTableMapping) <- c("TABLE","FILE")
  NRCSTableMapping[,"FILE"] <- paste(NRCSTableMapping[,"FILE"],'.txt',sep='')
  
  tablesData <- tablesData[as.character(NRCSTableMapping[,"FILE"])]
  tablesHeaders <- tablesHeaders[as.character(NRCSTableMapping[,"TABLE"])]
  
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

getNRCSData <- function(areas, raw.dir){
  # Load raw NRCS Map Unit polygons from the regions specified.
  cat("Loading the NRCS soil survey data for each region.\n")
  NRCSTables <- vector("list", length(areas))
  
  for(i in 1:length(areas)){
    area <- as.character(areas$areasymbol[i])
    date <- as.Date(areas$saverest[i])
    tables <- getNRCSStudyAreaData(area=area, date=date, raw.dir=raw.dir)
    NRCSTables[[i]] <- tables
  }
  
  # Merging all NRCS data tables
  cat("Merging all NRCS data tables\n")
  tableNames <- unique(unlist(sapply(NRCSTables,names)))
  tableNames <- tableNames[order(tableNames)]
  
  NRCSTables <- lapply(tableNames,function(name){
    tables <- lapply(NRCSTables,'[[',name)
    tables <- do.call("rbind", tables)
    tables <- unique(tables)
    return(tables)
  })
  
  names(NRCSTables) <- tableNames
  
  return(NRCSTables)
}

extractNRCSData <- function(tables,mapunits){
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