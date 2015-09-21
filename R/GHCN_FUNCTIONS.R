#' Download and crop the Global Historical Climate Network-Daily data.
#'
#' \code{get_ghcn_daily} returns a named list of length 2: 
#' \enumerate{
#' \item "spatial": A \code{SpatialPointsDataFrame} of the locations of GHCN weather stations 
#' in the template, and 
#' \item "tabular": A named list of \code{\link{data.frame}s} with the daily weather data for each station.
#' The name of each list item is the station ID.
#' }
#' 
#' @param template A Raster* or Spatial* object to serve 
#' as a template for cropping. If missing, all stations
#' will be downloaded!
#' @param label A character string naming the study area.
#' @param elements A character vector of elemets to extract.
#' Common elements include "tmin", "tmax", and "prcp".
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing. Defaults to "./RAW/GHCN/".
#' @param extraction.dir A character string indicating where the extracted and cropped GHCN shapefiles should be put.
#' The directory will be created if missing. Defaults to "./EXTRACTIONS/GHCN/".
#' @param standardize Select only common year/month/day? Defaults to FALSE.
#' @param force.redo If an extraction for this template and label already exists, should a new one be created? Defaults to FALSE.
#' @return A named list containing the "spatial" and "tabular" data.
#' @export
get_ghcn_daily <- function(template=NULL, label=NULL, elements=NULL, raw.dir="./RAW/GHCN/", extraction.dir="./EXTRACTIONS/GHCN/", standardize=F, force.redo=F){
  dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
  
  if(is.null(template) & is.null(label)){
    label <- "allStations"
  }
  
  if(!is.null(template) & is.null(label)){
    stop("Template provided but no label given.")
  }
  
  vectors.dir <- paste(extraction.dir,"/",label,"/spatial",sep='')
  tables.dir <- paste(extraction.dir,"/",label,"/tabular",sep='')
  
  dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(vectors.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(tables.dir, showWarnings = FALSE, recursive = TRUE)
  
  message("(Down)Loading GHCN station inventory.")
  if(!force.redo & file.exists(paste(vectors.dir,"/stations.shp",sep=''))){
    stations.sp <- rgdal::readOGR(dsn=vectors.dir,layer="stations",verbose=F)
  }else{
    stations.sp <- get_ghcn_inventory(template=template, raw.dir=raw.dir)
    suppressWarnings(rgdal::writeOGR(stations.sp, vectors.dir, "stations","ESRI Shapefile", overwrite_layer=TRUE))
  }
  
  # If the user didn't specify target elements, get them all.
  if(is.null(elements)){
    elements <- unique(stations.sp$ELEMENT)
  }
  
  stations.sp <- stations.sp[stations.sp@data[,"ELEMENT"] %in% toupper(elements),]
  
  if(standardize & !is.null(elements)){
    stations.sp.splits <- split(as.character(stations.sp$ELEMENT),f=stations.sp$ID, drop=T)
    stations.sp.splits.all <- sapply(stations.sp.splits,function(x){all(sapply(toupper(elements),function(y){y %in% x}))})
    stations.sp <- stations.sp[stations.sp$ID %in% names(stations.sp.splits.all)[stations.sp.splits.all],]
  }
  
  stations.out <- stations.sp[,c("ID","ELEMENT","YEAR_START","YEAR_END")]
  stations.sp <- stations.sp[!duplicated(stations.sp@data[,c("ID","LATITUDE","LONGITUDE")]),c("ID","LATITUDE","LONGITUDE")]
  
  if(!force.redo){
    daily <- tryCatch(lapply(elements,function(element){readRDS(paste(tables.dir,"/",element,".Rds",sep=''))}), warning = function(w){return(NULL)})
    if(!is.null(daily)){
      names(daily) <- elements
      
      daily <- lapply(as.character(stations.sp$ID),function(station){
        stationDaily <- tryCatch(lapply(daily,'[[',station),error=function(e){return(NULL)})
        stationDaily <- stationDaily[!sapply(stationDaily,is.null)]
        return(stationDaily)
      })
      names(daily) <- as.character(stations.sp$ID)
      daily <- daily[!sapply(daily,is.null)]
      # Make sure station names and elements are the same
      if(setequal(names(daily),stations.sp$ID) & all(sapply(daily,function(dat){setequal(names(dat),elements)}))){
        return(list(spatial=stations.out,tabular=daily))
      }
    }
  }
  
  daily <- lapply(stations.sp$ID,function(station){
    message("(Down)Loading GHCN station data for station ",as.character(station))
    return(get_ghcn_daily_station(ID=station, elements=elements, raw.dir=raw.dir, standardize=standardize, force.redo=force.redo))
  })
  names(daily) <- stations.sp$ID
  
  daily.split <- lapply(elements,function(element){
    lapply(daily,'[[',element)
  })
  names(daily.split) <- elements
  junk <- lapply(elements,function(element){
    saveRDS(daily.split[[element]],paste(tables.dir,"/",element,".Rds",sep=''),compress='xz')
  })
  
  
  return(list(spatial=stations.out,tabular=daily))
}

#' Download the daily data for a GHCN weather station.
#'
#' @param ID A character string giving the station ID.
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' @param force.redo If this weather station has been downloaded before, should it be updated? Defaults to FALSE.
#' @return A character string representing the full local path of the GHCN station data.
#' @export
download_ghcn_daily_station <- function(ID, raw.dir, force.redo=F){
  
  dir.create(raw.dir, recursive=T, showWarnings=F)
  
  url <- paste("ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/all/",ID,".dly",sep='')
  if(!force.redo){
    curl_download(url=url, destdir=raw.dir, timestamping=T)
  }else{
    curl_download(url=url, destdir=raw.dir, timestamping=F)
  }
  
  return(normalizePath(paste(raw.dir,ID,".dly",sep='')))
  
}

#' Download and extract the daily data for a GHCN weather station.
#'
#' \code{get_ghcn_daily_station} returns a named list of \code{\link{data.frame}s}, one for
#' each \code{elements}. If \code{elements} is undefined, it returns all available weather
#' tables for the station
#' 
#' @param ID A character string giving the station ID.
#' @param elements A character vector of elemets to extract.
#' Common elements include "tmin", "tmax", and "prcp".
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' @param standardize Select only common year/month/day? Defaults to FALSE.
#' @param force.redo If this weather station has been downloaded before, should it be updated? Defaults to FALSE.
#' @return A named list of \code{\link{data.frame}s}, one for each \code{elements}.
#' @export
get_ghcn_daily_station <- function(ID, elements=NULL, raw.dir, standardize=F, force.redo=F){
  
  file <- download_ghcn_daily_station(ID=ID, raw.dir=paste(raw.dir,"/DAILY/",sep=''), force.redo=force.redo)
  
  daily <- utils::read.fwf(file,c(11,4,2,4,rep(c(5,1,1,1),31)), stringsAsFactors=F)
  names(daily)[1:4] <- c("STATION","YEAR","MONTH","ELEMENT")
  
  if(is.null(elements)){
    elements <- unique(daily$ELEMENT)
  }
  
  daily <- daily[daily$ELEMENT %in% toupper(elements),c(2:4,seq(5,125,4))]
  daily[daily==-9999] <- NA
  names(daily) <- c("YEAR","MONTH","ELEMENT",paste("D",1:31,sep=''))
  
  ## Separate by element
  out.list <- lapply(elements, function(element){
    return(daily[daily$ELEMENT==toupper(element),-3])
  })
  
  ## If standardize, select only common year/month/day, and make NA if both not present
  if(standardize){
    yearMonths <- lapply(out.list, function(element){
      element <- element[order(element$YEAR,element$MONTH),]
      return(paste("Y",element[,c("YEAR")],"M",element[,c("MONTH")],sep=''))
    })
    
    all.yearMonths <- Reduce(intersect,yearMonths)
    
    out.list <- lapply(out.list, function(element){
      element.yearMonths <- paste("Y",element[,c("YEAR")],"M",element[,c("MONTH")],sep='')
      return(element[match(all.yearMonths,element.yearMonths),])
    })
    
    
  }
  
  names(out.list) <- elements
  
  return(out.list)
}

#' Download and crop the inventory of GHCN stations.
#'
#' \code{get_ghcn_inventory} returns a \code{SpatialPolygonsDataFrame} of the GHCN stations within
#' the specified \code{template}. If template is not provided, returns the entire GHCN inventory.
#' 
#' Stations with multiple elements will have multiple points. This allows for easy mapping of stations
#' by element availability.
#' 
#' @param template A Raster* or Spatial* object to serve 
#' as a template for cropping.
#' @param elements A character vector of elemets to extract.
#' Common elements include "tmin", "tmax", and "prcp".
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing.
#' @return A \code{SpatialPolygonsDataFrame} of the GHCN stations within
#' the specified \code{template}
#' @export
get_ghcn_inventory <- function(template=NULL, elements=NULL, raw.dir){
  if(!is.null(template) & (!(methods::is(template,"SpatialPolygonsDataFrame")) & !(methods::is(template,"SpatialPolygons")))){
    template <- polygon_from_extent(template)
  }
  
  template <- methods::as(template,"SpatialPolygons")
  
  url <- "ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/ghcnd-inventory.txt"
  destdir <- raw.dir
  curl_download(url=url, destdir=destdir)
  
  station.inventory <- utils::read.fwf(paste(raw.dir,"ghcnd-inventory.txt",sep=''),c(11,1,8,1,9,1,4,1,4,1,4), stringsAsFactors=F)[,seq(1,11,2)]
  names(station.inventory) <- c("ID","LATITUDE","LONGITUDE","ELEMENT","YEAR_START","YEAR_END")
  
  # Convert to SPDF
  stations.sp <- sp::SpatialPointsDataFrame(coords=station.inventory[,c("LONGITUDE","LATITUDE")],station.inventory,proj4string=sp::CRS("+proj=longlat +datum=WGS84"))

  if(!is.null(elements)){
    stations.sp <- stations.sp[stations.sp$ELEMENT %in% toupper(elements),]
  }

  if(!is.null(template)){
    stations.sp <- stations.sp[!is.na(sp::over(stations.sp,sp::spTransform(template,sp::CRS(raster::projection(stations.sp))))),]
  }
  
  return(stations.sp)
}

