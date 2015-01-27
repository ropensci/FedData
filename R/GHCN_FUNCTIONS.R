downloadGHCN <- function(ID, data.dir="../DATA/"){
  if(!file.exists(data.dir)){
    dir.create(data.dir, recursive=T)
  }
  for(i in ID){
    if(!file.exists(paste(data.dir,i,".dly",sep=''))){      
      system(paste("wget -nd -N --directory-prefix=",data.dir," ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/all//",i,".dly",sep=''))
    }
  }
}

getGHCN <- function(ID, data.dir="../DATA/", element, months=c(1:12), fun){
  if(!file.exists(data.dir)){
    dir.create(data.dir, recursive=T)
  }
  
  for(i in ID){
    if(!file.exists(paste(data.dir,i,".dly",sep=''))){
      system(paste("wget -nd --output-document=",data.dir,i,".dly"," ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/all//",i,".dly",sep=''))
    }
  }
  
  series <- mclapply(ID,FUN=function(x,...){extractGHCNMonthly(x, data.dir=data.dir, element=element, months=months, fun=fun)}, mc.cores=detectCores())
  
  # convert to list of zoo objects
  series <- lapply(series, function(x) zoo(as.matrix(x[ ,-1, drop = FALSE]),as.character(x[,1])))
  
  # perform merge
  series <- do.call('merge', series)
  
  # to convert back to data frame
  series <- data.frame(YEAR = time(series), coredata(series))
  
  names(series) <- c("YEAR",ID)
  
  return(series)
  
}

extractGHCNDaily <- function(id, data.dir="../DATA/", element){
  
  daily <- read.fwf(paste(data.dir,id,".dly",sep=''),c(11,4,2,4,rep(c(5,1,1,1),31)))
  names(daily)[1:4] <- c("STATION","YEAR","MONTH","ELEMENT")
  daily <- daily[daily$ELEMENT==toupper(element),c(2:3,seq(5,125,4))]
  daily[daily==-9999] <- NA
  
  return(daily)
}


extractGHCNDailySpatial <- function(template, elements=c("tmin","tmax","prcp"), data.dir="../DATA/", standardize=F){
  #   system(paste("wget -nd -N --output-document=",paste(data.dir,"ghcnd-stations.txt",sep='')," ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily//ghcnd-stations.txt",sep=''))
  #   system(paste("sed -i -E 's/#/ /' ",paste(data.dir,"ghcnd-stations.txt",sep=''),sep=''))
  #   system(paste("rm ",paste(data.dir,"ghcnd-stations.txt-E",sep=''),sep=''))
  system(paste("wget -nd -N --output-document=",paste(data.dir,"ghcnd-inventory.txt",sep='')," ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily//ghcnd-inventory.txt",sep=''))
  system(paste("sed -i -E 's/#/ /' ",paste(data.dir,"ghcnd-inventory.txt",sep=''),sep=''))
  system(paste("rm ",paste(data.dir,"ghcnd-inventory.txt-E",sep=''),sep=''))
  
  #   station.meta <- read.fwf(paste(data.dir,"ghcnd-stations.txt",sep=''),c(11,1,8,1,9,1,6,1,2,1,30,1,3,1,3,1,5))[,seq(1,11,2)]
  #   names(station.meta) <- c("ID","LATITUDE","LONGITUDE","ELEVATION","STATE","NAME")
  #   station.meta$NAME <- gsub("^\\s+|\\s+$", "", station.meta$NAME)
  
  station.inventory <- read.fwf(paste(data.dir,"ghcnd-inventory.txt",sep=''),c(11,1,8,1,9,1,4,1,4,1,4))[,seq(1,11,2)]
  names(station.inventory) <- c("ID","LATITUDE","LONGITUDE","ELEMENT","YEAR_START","YEAR_END")
  
  station.inventory <- station.inventory[station.inventory[,"ELEMENT"] %in% toupper(elements),]
  
  station.inventory.splits <- split(as.character(station.inventory[,c("ELEMENT")]),f=station.inventory$ID, drop=T)
  station.inventory.splits.all <- sapply(station.inventory.splits,function(x){all(sapply(toupper(elements),function(y){y %in% x}))})

  station.inventory <- station.inventory[station.inventory$ID %in% names(station.inventory.splits.all)[station.inventory.splits.all],]
  station.inventory <- unique(station.inventory[,c("ID","LATITUDE","LONGITUDE")])
  
  # Convert to SPDF
  stations.sp <- SpatialPointsDataFrame(coords=station.inventory[,c("LONGITUDE","LATITUDE")],station.inventory[,c("ID"),drop=F],proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  
  stations.sp <- spTransform(stations.sp,CRS(projection(template)))
  
  stations.sp <- crop.to.studyArea(GHCN.stations,polygonFromExtent(template))
  
  downloadGHCN(stations.sp$ID,data.dir=data.dir)
  
  
  daily <- lapply(stations.sp$ID,function(id){
    daily <- read.fwf(paste(data.dir,id,".dly",sep=''),c(11,4,2,4,rep(c(5,1,1,1),31)))
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
  
  
  
  return(list(stations.sp,daily))
}

extractGHCNMonthly <- function(id, data.dir="../DATA/", element, months=c(1:12), fun=NULL){
  # If more months than a year, break
  if(length(months)>12){
    stop("ERROR! Too many months.")
  }
  
  # Process the months to get months from current, previous, and future years
  previous.year <- 12+months[months<1]
  current.year <- months[months>=1 & months<=12]
  next.year <- months[months>12]-12
  
  daily <- read.fwf(paste(data.dir,id,".dly",sep=''),c(11,4,2,4,rep(c(5,1,1,1),31)))
  names(daily)[1:4] <- c("STATION","YEAR","MONTH","ELEMENT")
  daily <- daily[daily$ELEMENT==toupper(element),c(2:3,seq(5,125,4))]
  daily[daily==-9999] <- NA
  
  if(fun=="sum"){
    daily <- cbind(daily[,1:2],AGG=rowSums(daily[,3:33], na.rm=T))
  }else if(fun=="mean"){
    daily <- cbind(daily[,1:2],AGG=rowMeans(daily[,3:33], na.rm=T))
  }
  daily$YEAR[daily$MONTH %in% previous.year] <- daily$YEAR[daily$MONTH %in% previous.year]+1
  daily <- daily[!(daily$MONTH %in% (1:12)[!(1:12 %in% c(previous.year,current.year,next.year))]),]
  daily <- daily[complete.cases(daily),]
  if(fun=="sum"){
    monthly <- aggregate(daily, by=list(daily$YEAR), FUN=sum)[,c(1,4)]
  }else if(fun=="mean"){
    monthly <- aggregate(daily, by=list(daily$YEAR), FUN=mean)[,c(1,4)]
  }
  names(monthly) <- c("YEAR",id)
  
  return(monthly)
}


getGHCNStationMetadataSPDF <- function(station.names=NULL, data.dir="../DATA/", write.output=T){
  if(!file.exists(paste(data.dir,"ghcnd-stations.txt",sep=''))){
    system(paste("wget -nd --output-document=",paste(data.dir,"ghcnd-stations.txt",sep='')," ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily//ghcnd-stations.txt",sep=''))
    system(paste("sed -i -E 's/#/ /' ",paste(data.dir,"ghcnd-stations.txt",sep=''),sep=''))
    system(paste("rm ",paste(data.dir,"ghcnd-stations.txt-E",sep=''),sep=''))
    system(paste("wget -nd --output-document=",paste(data.dir,"ghcnd-inventory.txt",sep='')," ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily//ghcnd-inventory.txt",sep=''))
    system(paste("sed -i -E 's/#/ /' ",paste(data.dir,"ghcnd-inventory.txt",sep=''),sep=''))
    system(paste("rm ",paste(data.dir,"ghcnd-inventory.txt-E",sep=''),sep=''))
  }
  
  # Read the USHCN station metadata
  station.meta <- read.fwf(paste(data.dir,"ghcnd-stations.txt",sep=''),c(11,1,8,1,9,1,6,1,2,1,30,1,3,1,3,1,5))[,seq(1,11,2)]
  names(station.meta) <- c("ID","LATITUDE","LONGITUDE","ELEVATION","STATE","NAME")
  station.meta$NAME <- gsub("^\\s+|\\s+$", "", station.meta$NAME)
  
  if(!is.null(station.names)){
    station.meta <- station.meta[station.meta$ID %in% station.names,]
  }
  
  # Convert to SPDF
  stations.sp <- SpatialPointsDataFrame(coords=station.meta[,c("LONGITUDE","LATITUDE")],station.meta,proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  
  stations.sp$NAME <- stations.sp$ID
  
  if(write.output){
    suppressWarnings(writeOGR(stations.sp,paste(data.dir,"ghcnd-stations.kml",sep=''),"ghcnd-stations",driver="KML", overwrite_layer=TRUE))
  }
  
  return(stations.sp)
}
