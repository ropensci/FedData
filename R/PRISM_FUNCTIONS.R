getPRISM_MONTHLYData <- function(template, type, label, out.dir, monthly.dir, force.redo=FALSE){
  if(!file.exists(paste(out.dir,label,"/", sep=''))){
    dir.create(paste(out.dir,label,"/", sep=''))
  }
  
  if(!force.redo & file.exists(paste(out.dir,label,"/",type,"_MONTHLY_FINAL_800.tif", sep=''))){
    all.data <- brick(paste(out.dir,label,"/",type,"_MONTHLY_FINAL_800.tif", sep=''))*1
    names(all.data) <- read.csv(paste(out.dir,label,"/",type,"_MONTHLY_FINAL_800_BANDS.csv", sep=''))$x
    #     names(all.data) <- as.Date(read.csv(paste(out.dir,label,"/",type,"_MONTHLY_FINAL_800_BANDS.csv", sep=''))$x, format="Y%Y.M%m")
    return(all.data)
  }
  
  # Convert template to decimal degrees for data request
  sim.poly.latlon <- spTransform(template, CRS("+proj=longlat +datum=NAD83"))
  
  monthly.files <- list.files(paste(monthly.dir,type,sep=''), recursive=T, full.names=T)
  
  if(!length(grep("*\\.bil$", monthly.files, value=TRUE))>0){
    system(paste("ls ",paste(monthly.dir,type,sep=''),'/* | parallel unzip -qqo {} -d ',paste(monthly.dir,type,sep=''),sep='')) 
    monthly.files <- list.files(paste(monthly.dir,type,sep=''), recursive=T, full.names=T)
  }
  
  monthly.files <- grep("*\\.bil$", monthly.files, value=TRUE)
  monthly.files <- grep("spqc", monthly.files, value=TRUE, invert=T)
  monthly.files <- grep("/cai", monthly.files, value=TRUE)
  
  all.data <- mclapply(as.list(monthly.files),function(monthly.file,...){ extractPRISM_MONTHLYMonth(template=sim.poly.latlon, file=monthly.file) }, mc.cores=detectCores())
  
  all.data <- brick(all.data)
  
  all.data <- all.data[[order(names(all.data))]]
  
  if(!file.exists(paste(out.dir,label,"/", sep=''))){
    dir.create(paste(out.dir,label,"/", sep=''))
  }
  
  writeRaster(all.data,paste(out.dir,label,"/",type,"_MONTHLY_FINAL_800.tif", sep=''), datatype="FLT4S", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"),overwrite=T,setStatistics=FALSE)
  #   writeGDAL(as(all.data, "SpatialGridDataFrame"), paste(out.dir,label,"/",type,"_MONTHLY_FINAL_800.nc", sep=''), drivername="netCDF", type="Float32", options=c("FORMAT=NC4","COMPRESS=DEFLATE","ZLEVEL=9"))
  write.csv(names(all.data),paste(out.dir,label,"/",type,"_MONTHLY_FINAL_800_BANDS.csv", sep=''), row.names=F, col.names=F)
  
#   unlink(paste(out.dir,'TEMP/', sep=''), recursive=T, force=T)
  
  gc()
  gc()
  
  return(NULL)
}

extractPRISM_MONTHLYMonth <- function(template=NULL, file){
  rast <- raster(file)
  projection(rast) <- CRS("+proj=longlat +datum=NAD83")
  
  # Crop mosaic to simulation extent
  if(!is.null(template)){
    rast <- crop(rast,template, snap="out") 
  }
  
  name <- names(rast)
  name <- substrRight(name, 6)
  name <- paste("Y",substr(name, 1, 4),".M",substr(name, 5, 6),sep="")
  names(rast) <- name
  # Return the cropped raster
  return(rast)
}



extractPRISM_MONTHLY <- function(template, out.dir, LT81.dir, label, force.redo=FALSE){
  
  types <- c("ppt","tmin","tmax")
  
  lapply(types,function(type){ getPRISM_MONTHLYData(type=type, label=label, template=template, out.dir=out.dir, monthly.dir=LT81.dir, force.redo=force.redo) })
  
  return(cat("PRISM Extracted\n"))
}


calcANNUALGDD_MONTHLY <- function(extraction.dir, months, t.base, t.cap=NULL){
  tmin.brick <- brick(paste(extraction.dir,"/tmin_MONTHLY_FINAL_800.tif", sep=''))
  names(tmin.brick) <- read.csv(paste(extraction.dir,"/tmin_MONTHLY_FINAL_800_BANDS.csv", sep=''), stringsAsFactors=F)[,1]
  
  tmax.brick <- brick(paste(extraction.dir,"/tmax_MONTHLY_FINAL_800.tif", sep=''))
  names(tmax.brick) <- read.csv(paste(extraction.dir,"/tmax_MONTHLY_FINAL_800_BANDS.csv", sep=''), stringsAsFactors=F)[,1]
  
  if(nlayers(tmin.brick)!=nlayers(tmax.brick) | any(names(tmin.brick)!=names(tmax.brick))){
    stop("tmin and tmax bricks must have same number of layers!")
  }
  
  ## Calculate an estimate of accumulated growing season GDD (Tbase=10, Tut=30)
  # Floor tmax and tmin at Tbase
  tmin.brick <- calc(tmin.brick,function(x) { x[x<t.base] <- t.base; return(x) })
  tmax.brick <- calc(tmax.brick,function(x) { x[x<t.base] <- t.base; return(x) })
  # Cap tmax and tmin at Tut
  if(!is.null(t.cap)){
    tmin.brick <- calc(tmin.brick,function(x) { x[x>t.cap] <- t.cap; return(x) })
    tmax.brick <- calc(tmax.brick,function(x) { x[x>t.cap] <- t.cap; return(x) })
  }
  
  GDD.brick <- ((tmin.brick+tmax.brick)/2)-t.base
  GDD.brick.months <- as.numeric(gsub("Y\\d{4}[.]M","",names(GDD.brick)))
  
  year.months <- 1:12
  days.per.month <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  
  GDD.brick.days <- as.numeric(mapply(gsub, year.months, days.per.month, GDD.brick.months))
  
  # Multiply by days per month, and convert to Fahrenheit GDD
  GDD.brick <- GDD.brick*GDD.brick.days*1.8
  
  # sum over months
  GDD.brick.total <- annualizePRISM_MONTHLY(GDD.brick, months=months, fun="sum")
  
  return(GDD.brick.total)
}

subsetPRISM_MONTHLY <- function(extraction.dir, element, months, fun){
  
  if(element=="prcp"){element <- "ppt"}
  
  all.brick <- brick(paste(extraction.dir,"/",element,"_MONTHLY_FINAL_800.tif", sep=''))
  names(all.brick) <- read.csv(paste(extraction.dir,"/",element,"_MONTHLY_FINAL_800_BANDS.csv", sep=''), stringsAsFactors=F)[,1]
  
  brick.names <- names(brick)
  brick.years <- gsub("([.].*)","",brick.names)
  brick.years <- as.numeric(gsub("Y","",brick.years))
  brick.months <- as.numeric(gsub("Y\\d{4}[.]M","",brick.names))
  
  annual.brick <- annualizePRISM_MONTHLY(brick=all.brick,months=months,fun=fun)
  
  return(annual.brick)
}

annualizePRISM_MONTHLY <- function(prism.brick, months=c(1:12), fun){
  brick.names <- names(prism.brick)
  brick.years <- gsub("([.].*)","",brick.names)
  brick.years <- as.numeric(gsub("Y","",brick.years))
  brick.months <- as.numeric(gsub("Y\\d{4}[.]M","",brick.names))
  
  # If more months than a year, break
  if(length(months)>12){
    stop("ERROR! Too many months.")
  }
  
  # Process the months to get months from current, previous, and future years
  previous.year <- 12+months[months<1]
  current.year <- months[months>=1 & months<=12]
  next.year <- months[months>12]-12
  no.year <- (1:12)[!(1:12 %in% c(previous.year,current.year,next.year))]
  
  brick.years[brick.months %in% previous.year] <- brick.years[brick.months %in% previous.year]+1
  brick.years[brick.months %in% next.year] <- brick.years[brick.months %in% next.year]-1
  brick.years[brick.months %in% no.year] <- 0
  
  signal.brick.temp <- setZ(prism.brick, brick.years, "year")
  
  if(fun=="sum"){
    signal.brick.temp <- zApply(signal.brick.temp, by=getZ(signal.brick.temp), fun=sum)
  }else if(fun=="mean"){
    signal.brick.temp <- zApply(signal.brick.temp, by=getZ(signal.brick.temp), fun=mean)
  }
  signal.brick.temp <- signal.brick.temp[[(1:nlayers(signal.brick.temp))[names(signal.brick.temp)!="X0"]]]
  #   signal.brick.temp <- round(signal.brick.temp, digits=3)
  
  return(signal.brick.temp)
}







getPRISM_MONTHLYPointsData <- function(sp.points, type, LT81.dir, year.range, months=c(1:12), fun='mean'){
  
  brick.years <- rep(year.range,each=12)
  brick.months <- rep(1:12,length(year.range))
  
  # If more months than a year, break
  if(length(months)>12){
    stop("ERROR! Too many months.")
  }
  
  # Process the months to get months from current, previous, and future years
  previous.year <- 12+months[months<1]
  current.year <- months[months>=1 & months<=12]
  next.year <- months[months>12]-12
  no.year <- (1:12)[!(1:12 %in% c(previous.year,current.year,next.year))]
  
  brick.years[brick.months %in% previous.year] <- brick.years[brick.months %in% previous.year]+1
  brick.years[brick.months %in% next.year] <- brick.years[brick.months %in% next.year]-1
  brick.years[brick.months %in% no.year] <- 0
  
  all.data <- lapply(year.range,function(year){ extractPRISM_MONTHLYPointsYear(type=type,sp.points=sp.points,LT81.dir=LT81.dir,year=year) })
  
  all.data <- do.call('rbind',all.data)
  
  if(fun=='mean'){
    all.data <- aggregate(all.data,by=list(brick.years), FUN=mean)
  }else if(fun=='sum'){
    all.data <- aggregate(all.data,by=list(brick.years), FUN=sum)
  }
  
  all.data <- all.data[all.data$Group.1 != 0,-1]
  
  if(fun=='mean'){
    all.data <- as.numeric(colMeans(all.data))
  }else if(fun=='sum'){
    all.data <- as.numeric(colSums(all.data))
  }
  
  return(all.data)
}

extractPRISM_MONTHLYPointsYear <- function(type, out.dir, LT81.dir, sp.points,year){
  unzip(paste(LT81.dir,type,'/',type,'_',year,'.zip',sep=''),exdir='TEMP/')
  
  months <- 1:12
  
  year.data <- mclapply(months,function(month){ extractPRISM_MONTHLYPointsMonth(type=type,sp.points=sp.points,year=year,month=month) }, mc.cores=detectCores())
  
  year.data <- do.call('rbind',year.data)
  
  unlink('TEMP/', recursive=T)
  
  return(year.data)
}

extractPRISM_MONTHLYPointsMonth <- function(type,sp.points,year,month){
  month.char <- formatC(month, width = 2, format = "d", flag = "0")
  
  file.path <- paste('TEMP/',year,'/cai_',type,'_us_us_30s_',year,month.char,'.bil', sep='')
  
  rast <- extractPRISMPointsFromFile(file.path=file.path,sp.points=sp.points)
  return(rast)
}

extractPRISM_MONTHLYPointsFromFile <- function(file.path,sp.points){
  temp <- raster(file.path)
  
  # Extract values at sp.points simulation extent
  temp <- extract(temp,sp.points) 
  
  # Return the cropped raster
  return(temp)
}









extractPRISM_DAILY <- function(template, types=c("ppt","tmin","tmax"), out.dir, daily.dir, label, year.range, force.redo=FALSE){
  poly <- polygonFromExtent(extent(template),projection(template))
  
  junk <- lapply(types,function(type){ getPRISM_DAILYData(type=type, label=label, template=template, out.dir=out.dir, daily.dir=daily.dir, year.range=year.range, force.redo=force.redo) })
  
  rm(junk)
  
  return(cat("PRISM Extracted\n"))
}

getPRISM_DAILYData <- function(type, label, out.dir, daily.dir, template, year.range, date.format="X%Y.%m.%d",force.redo=FALSE){
  if(!file.exists(paste(out.dir,label,"/", sep=''))){
    dir.create(paste(out.dir,label,"/", sep=''))
  }
  
  if(!force.redo & file.exists(paste(out.dir,label,"/",type,"_DAILY_FINAL_800.tif", sep=''))){
    all.data <- brick(paste(out.dir,label,"/",type,"_DAILY_FINAL_800.tif", sep=''))*1
    names(all.data) <- base::as.Date(read.csv(paste(out.dir,label,"/",type,"_DAILY_FINAL_800_BANDS.csv", sep=''),stringsAsFactors=F)$x, format=date.format)
    return(all.data)
  }
  
  # Convert UTMs to decimal degrees for data request
  sim.poly.latlon <- spTransform(template, CRS("+proj=longlat +datum=NAD83"))
  
  all.data <- mclapply(year.range,function(year){ extractPRISM_DAILYYear(type=type,template=sim.poly.latlon,out.dir=out.dir,daily.dir=daily.dir,year=year) },mc.cores=detectCores())
  
  all.data <- stack(unlist(all.data), quick=T)
  
  if(!file.exists(paste(out.dir,label,"/", sep=''))){
    dir.create(paste(out.dir,label,"/", sep=''))
  }
  
  writeRaster(all.data,paste(out.dir,label,"/",type,"_DAILY_FINAL_800.tif", sep=''), datatype="FLT4S", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"),overwrite=T,setStatistics=FALSE)
  
  #   writeRaster(all.data, paste(out.dir,label,"/",type,"_DAILY_FINAL_800.tif", sep=''), drivername="GTiff", type="Float32")
  #   writeGDAL(as(all.data, "SpatialGridDataFrame"), paste(out.dir,label,"/",type,"_DAILY_FINAL_800.tif", sep=''), drivername="GTiff", type="Float32")
  write.csv(names(all.data),paste(out.dir,label,"/",type,"_DAILY_FINAL_800_BANDS.csv", sep=''), row.names=F, col.names=F)
  
  gc()
  
  return(all.data)
}

extractPRISM_DAILYYear <- function(type,out.dir,daily.dir,template=NULL,year){
  unlink(paste(daily.dir,'/',type,'/',year,'/TEMP/',sep=''), recursive=T, force=T)
  
  daily.zips <- list.files(paste(daily.dir,'/',type,'/',year,'/',sep=''), full.names=T)
  junk <- sapply(daily.zips, function(i) unzip(i, exdir=paste(daily.dir,'/',type,'/',year,'/TEMP/',sep='')))
  rm(junk)          
  
  daily.files <- list.files(paste(daily.dir,'/',type,'/',year,'/TEMP',sep=''), full.names=T)
  daily.files <- grep("*\\.bil$", daily.files, value=TRUE)
  daily.files <- grep("early", daily.files, value=TRUE, invert=T)
  
  
  year.data <- suppressWarnings(stack(daily.files, quick=T, native=T))
  projection(year.data) <- CRS("+proj=longlat +datum=NAD83")
  year.data <- crop(year.data,template, snap='out') 
  
  names.data <- names(year.data)
  names.data <- gsub(paste("PRISM_",type,"_stable_4kmD1_", sep=''),"",names.data,)
  names.data <- gsub(paste("PRISM_",type,"_provisional_4kmD1_", sep=''),"",names.data,)
  names.data <- gsub("_bil","",names.data,)
  names.data <- as.Date(names.data,"%Y%m%d")
  names.data <- as.character(names.data,"Y%Y.M%m.D%d")
  names(year.data) <- names.data
  
  year.data <- year.data[[order(names(year.data))]]
  
  unlink(paste(daily.dir,'/',type,'/',year,'/TEMP/',sep=''), recursive=T)
  
  return(year.data)
}