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

getGHCNDaily <- function(template=NULL, elements=NULL, label=NULL, raw.dir="./RAW/GHCN/", extraction.dir="./EXTRACTIONS/GHCN/", standardize=F, force.redo=F){
  dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
  
  if(is.null(template)){
    label <- "allStations"
  }
  
  vectors.dir <- paste(extraction.dir,"/",label,"/spatial",sep='')
  tables.dir <- paste(extraction.dir,"/",label,"/tabular",sep='')
  
  dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(vectors.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(tables.dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("\nGetting spatial data of GHCN stations")
  if(!force.redo & file.exists(paste(vectors.dir,"/stations.shp",sep=''))){
    stations.sp <- rgdal::readOGR(dsn=vectors.dir,layer="stations")
  }else{
    stations.sp <- getGHCNStations(template=template, raw.dir=raw.dir)
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
        return(tryCatch(lapply(daily,'[[',station),error=function(e){return(NULL)}))
      })
      names(daily) <- as.character(stations.sp$ID)
      daily <- daily[!sapply(daily,is.null)]
      # Make sure station names and elements are the same
      if(all(names(daily) %in% stations.sp$ID) & all(sapply(daily,function(dat){all(names(dat) %in% elements)}))){
        return(list(spatial=stations.out,tabular=daily))
      }
    }
  }
  
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