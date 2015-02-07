getGHCNStations <- function(template=NULL, elements=NULL, standardize=F, raw.dir){
  if(!is.null(template) & !is(template,"SpatialPolygonsDataFrame")){
    template <- polygonFromExtent(template)
  }
  
  url <- "ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/ghcnd-inventory.txt"
  destdir <- raw.dir
  wgetDownload(url=url, destdir=destdir)
    
  system(paste("sed -i -E 's/#/ /' ",paste(raw.dir,"ghcnd-inventory.txt",sep=''),sep=''))
  system(paste("rm ",paste(raw.dir,"ghcnd-inventory.txt-E",sep=''),sep=''))
  
  station.inventory <- utils::read.fwf(paste(raw.dir,"ghcnd-inventory.txt",sep=''),c(11,1,8,1,9,1,4,1,4,1,4))[,seq(1,11,2)]
  names(station.inventory) <- c("ID","LATITUDE","LONGITUDE","ELEMENT","YEAR_START","YEAR_END")
  
  if(!is.null(elements)){
    station.inventory <- station.inventory[station.inventory[,"ELEMENT"] %in% toupper(elements),]
  }
  
  if(standardize & !is.null(elements)){
    station.inventory.splits <- split(as.character(station.inventory[,c("ELEMENT")]),f=station.inventory$ID, drop=T)
    station.inventory.splits.all <- sapply(station.inventory.splits,function(x){all(sapply(toupper(elements),function(y){y %in% x}))})
    station.inventory <- station.inventory[station.inventory$ID %in% names(station.inventory.splits.all)[station.inventory.splits.all],]
  }
  
#   station.inventory <- unique(station.inventory[,c("ID","LATITUDE","LONGITUDE")])
  
  # Convert to SPDF
  stations.sp <- SpatialPointsDataFrame(coords=station.inventory[,c("LONGITUDE","LATITUDE")],station.inventory,proj4string=sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  
  if(!is.null(template)){
    stations.sp <- stations.sp[!is.na(sp::over(stations.sp,sp::spTransform(template,sp::CRS(projection(stations.sp))))[,1]),]
  }

  return(stations.sp)
}

downloadGHCNDaily <- function(ID, raw.dir="../DATA/", force.redo=F){
  if(!file.exists(raw.dir)){
    dir.create(raw.dir, recursive=T)
  }
  
  urls <- paste("ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/all/",ID,".dly",sep='')
  out.files <- paste(raw.dir,ID,'.dly',sep='')
  
  if(!force.redo){
    existing.files <- sapply(out.files,file.exists)
    urls <- urls[!existing.files]
    out.files <- out.files[!existing.files]
  }

  out <- mapply(function(x,y) download.file(url=x,destfile=y,mode='wb'), urls, out.files)

}

getGHCNDaily <- function(template=NULL, elements=NULL, raw.dir="./RAW/GHCN/", standardize=F, force.redo=F){
  dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("\nGetting spatial data of GHCN stations")
  stations.sp <- getGHCNStations(template=template, raw.dir=raw.dir, elements=elements, standardize=standardize)
  
  # If the user didn't specify target elements, get them all.
  if(is.null(elements)){
    elements <- unique(stations.sp$ELEMENT)
  }
  
  stations.sp <- stations.sp[!duplicated(stations.sp@data[,c("ID","LATITUDE","LONGITUDE")]),c("ID","LATITUDE","LONGITUDE")]
  
  cat("\nDownloading daily GHCN data")
  downloadGHCNDaily(stations.sp$ID,raw.dir=raw.dir,force.redo=force.redo)
  
  daily <- lapply(1:length(stations.sp$ID),function(i){
    id <- stations.sp$ID[i]
    cat("\nProcessing daily GHCH station ",i,"of",length(stations.sp$ID))
    daily <- utils::read.fwf(paste(raw.dir,id,".dly",sep=''),c(11,4,2,4,rep(c(5,1,1,1),31)))
    names(daily)[1:4] <- c("STATION","YEAR","MONTH","ELEMENT")
    daily <- daily[daily$ELEMENT %in% toupper(elements),c(2:4,seq(5,125,4))]
    daily[daily==-9999] <- NA
    names(daily) <- c("YEAR","MONTH","ELEMENT",paste("D",1:31,sep=''))
    
    ## Seperate by element
    out.list <- lapply(elements, function(element){
      return(daily[daily$ELEMENT==toupper(element),-3])
    })
    
    ## If standardize, select only common year/month, and make NA if both not present
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
  })
  
  return(list(spatial=stations.sp,tabular=daily))
}