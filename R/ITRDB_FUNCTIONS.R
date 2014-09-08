getITRDB <- function(raw.dir, output.dir, type='standard', download=F, force.redo=F){
  
  raw.dir <- paste(raw.dir,type,"/",sep='')
  dir.create(raw.dir, recursive=T, showWarnings=F)
  dir.create(output.dir, recursive=T, showWarnings=F)
  
  # First, get a list of all *numeric*.crn files in the treering section of the NCDC ftp site.
  # Only gets "standard" chronologies, i.e., files without a letter suffix.
  if(download){
    if(type=='standard'){
      system(paste('wget -ignore-case -np -nH -nd -nc -r -A "*[0-9][0-9].crn" --directory-prefix=',raw.dir,' ftp://ftp.ncdc.noaa.gov/pub/data/paleo/treering/chronologies/',sep=''))
    }else if(type=='whitened'){
      system(paste('wget -ignore-case -np -nH -nd -nc -r -A "*[0-9][0-9]r.crn" --directory-prefix=',raw.dir,' ftp://ftp.ncdc.noaa.gov/pub/data/paleo/treering/chronologies/',sep=''))
    }else if(type=='arstan'){
      system(paste('wget -ignore-case -np -nH -nd -nc -r -A "*[0-9][0-9]a.crn" --directory-prefix=',raw.dir,' ftp://ftp.ncdc.noaa.gov/pub/data/paleo/treering/chronologies/',sep=''))
    }else{
      stop("ERROR! Unknown chronology type. Please correct and retry.\n")
    }
  }
  
  ## A vector of the files in the output.dir
  dendros <- paste(raw.dir,list.files(raw.dir),sep='')
  
  if(file.exists(paste(output.dir,"/ITRDB_DATA.csv",sep='')) & !force.redo){
    existing.data <- read.csv(paste(output.dir,"/ITRDB_DATA.csv",sep=''))
    if((ncol(existing.data)-1)==length(dendros)){
      return(existing.data)
    }
  }
  
  # Some files are encoded in latin-1... change them to utf8
  #    blank <- lapply(dendros,FUN=function(x){system(paste("recode l1..utf8 ",x,sep='')); cat("Recoding ",x,"\n"); return()})
  
  # Get all series
  series.list <- lapply(dendros,function(x,...){ seriesToMatrix(x, raw.dir) })
  series.list.error <- dendros[which(grepl("ERROR",series.list))]
  if(length(series.list.error)>0){
    cat("ERROR! One or more improperly formatted chronologies. Please correct and retry.\n")
    print(series.list.error)
    stop()
  }
  
  # Extract and format metadata
  meta <- as.data.frame(do.call("rbind",lapply(series.list,FUN=function(x){x$meta})))
  names(meta) <- c("SERIES","NAME","SPECIES","ELEVATION","LAT","LON","START","END","CONTRIBUTOR")
  meta$ELEVATION <- gsub("M", "", meta$ELEVATION)
  meta$ELEVATION <- gsub("m", "", meta$ELEVATION)
  meta$CONTRIBUTOR <- gsub("  ", " ", meta$CONTRIBUTOR)
  SHemi <- grepl("-",meta$LAT)
  WHemi <- grepl("-",meta$LON)
  meta$LAT <- gsub("+", "", meta$LAT,fixed = TRUE)
  meta$LON <- gsub("-", "", meta$LON,fixed = TRUE)
  meta$CONTRIBUTOR <- toupper(meta$CONTRIBUTOR)
  meta$NAME <- toupper(meta$NAME)
  meta$SERIES <- toupper(meta$SERIES)
  meta$LAT <- format(as.numeric(meta$LAT)/100, nsmall=2)
  meta$LON <- format(as.numeric(meta$LON)/100,nsmall=2)
  
  locations <- meta[c("LAT","LON")]
  
  locations$LAT <- within(data.frame(lats=locations$LAT), {
    dms <- do.call(rbind, strsplit(as.character(lats), "\\."))
    LAT <- as.numeric(dms[,1]) + (as.numeric(dms[,2]))/60
    rm(dms)
    rm(lats)
  })
  
  locations$LON <- within(data.frame(lons=locations$LON), {
    dms <- do.call(rbind, strsplit(as.character(lons), "\\."))
    LON <- as.numeric(dms[,1]) + (as.numeric(dms[,2]))/60
    rm(dms)
    rm(lons)
  })
  
  meta$LAT <- locations$LAT
  meta$LON <- locations$LON
  
  if(any(SHemi)){
    meta[SHemi,]$LAT <- 0-meta[SHemi,]$LAT
  }
  
  if(any(WHemi)){
    meta[WHemi,]$LON <- 0-meta[WHemi,]$LON 
  }
  
  meta <- as.data.frame(lapply(meta,function(x) if(is.character(x)|is.factor(x)) gsub(", ","-",x) else x))
  
  # Merge all series into one DF
  data <- lapply(series.list,FUN=function(x){data.table::data.table(x$data, key=names(x$data))})
  series.merged <- as.data.frame(Reduce(function(...){merge(..., by="YEAR", all=T, sort=F)}, data))
  series.merged$YEAR <- as.numeric(series.merged$YEAR)
  series.merged <- series.merged[order(series.merged$YEAR),]
  
  # Sort in various ways
  YEAR <- series.merged$YEAR
  series.merged <- series.merged[,c(-1)]
  non.nas <- as.numeric(sapply(series.merged,FUN=function(x){match(FALSE,is.na(x))}))
  series.merged <- series.merged[,order(non.nas)]
  names(series.merged) <- toupper(names(series.merged))
  series.merged <- cbind(YEAR,series.merged)
  
  # Sort the metadata to match the order of the data
  meta$SERIES <- as.character(meta$SERIES)
  ordering <- unlist(lapply(2:length(names(series.merged)),FUN=function(x,...){which(meta$SERIES==names(series.merged)[x])}))
  meta <- meta[ordering,]
  
  # Write the metadata
  write.csv(meta,paste(output.dir,"/ITRDB_METADATA.csv",sep=''),row.names=F)
  
  # Write the final standardized series
  write.csv(series.merged,paste(output.dir,"/ITRDB_DATA.csv",sep=''),row.names=F)
  
  return(series.merged)
}

## A function to get a series and return it as a (year, z(width)) matrix
seriesToMatrix <- function(series, output.dir){
  cat("Processing ",series,"\n")
  
  if(suppressWarnings(tail(scan(series,what='character',nlines=1, quiet=T),n=1)!="RAW")){
    # Parse the header of the CRN file
    meta.1 <- as.character(read.fwf(series,c(6,3,52,4),skip=0,n=1,colClasses = "character", strip.white=T))
    meta.2 <- as.character(read.fwf(series,c(6,3,13,18,6,5,6,9,6,5),skip=1,n=1,colClasses = "character", strip.white=T))
    meta.3 <- as.character(read.fwf(series,c(6,3,60,12),skip=2,n=1,colClasses = "character", strip.white=T))
    
    # Get the name of the series
    name <- gsub(".crn","",series)
    name <- gsub(paste(output.dir,sep=''),"",name)
    
    meta <- c(name,meta.1[3:4],meta.2[c(5:7,9:10)],meta.3[3])
    
    digits.year <- max(nchar(meta[7]),nchar(meta[8]),4)
    
    fname <- series
    header=T
    encoding = getOption("encoding")
    ## Open the data file for reading
    con <- file(fname, encoding = encoding)
    on.exit(close(con))
    if(is.null(header)){
      ## Try to determine if the file has a header. This is failable.
      ## Find out if an ITRDB header (3 lines) in file
      hdr1 <- readLines(con, n=1)
      if(length(hdr1) == 0)
        stop("file is empty")
      if(nchar(hdr1) < 10)
        stop("first line in the crn file ends before col 10")
      yrcheck <- suppressWarnings(as.numeric(substr(hdr1, 7, 10)))
      if(is.null(yrcheck) || length(yrcheck)!=1 || is.na(yrcheck) |
           yrcheck < -1e04 || yrcheck > 1e04) {
        #         cat(gettext("There appears to be a header in the crn file\n",
        #                     domain="R-dplR"))
        is.head <- TRUE
      }
      else {
        #         cat(gettext("There does not appear to be a header in the crn file\n",
        #                     domain="R-dplR"))
        is.head <- FALSE # No header lines
      }
    } else if(!is.logical(header)){
      stop("'header' must be NULL or logical")
    } else{
      is.head <- header
    }
    if(is.head){
      ## Read 4th line - should be first data line
      dat1 <- readLines(con, n=4)
      if(length(dat1) < 4)
        stop("file has under 4 lines")
      dat1 <- dat1[4]
    } else{
      dat1 <- readLines(con, n=1)
      if(length(dat1) == 0)
        stop("file is empty")
    }
    if(nchar(dat1) < 10)
      stop("first data line ends before col 10")
    #     yrcheck <- as.numeric(substr(dat1, 7, 10))
    #     if(is.null(yrcheck) || length(yrcheck)!=1 || is.na(yrcheck) ||
    #          yrcheck < -1e04 || yrcheck > 1e04)
    #       stop(gettextf("cols %d-%d of first data line not a year", 7, 10,
    #                     domain="R-dplR"))
    ## Look at last line to determine if Chronology Statistics are present
    ## if nchar <=63 then there is a stats line
    nlines <- length(readLines(con, n=-1))
    ## Read file
    skip.lines <- ifelse(is.head, 3, 0)
    ## Do nothing. read.fwf closes (and destroys ?!?) the file connection
    on.exit()
    ## Get chron stats if needed
    suppressWarnings(chron.stats <- read.fwf(con, c(6, digits.year, 6, 6, 6, 7, 9, 9, 10),
                                             skip=nlines-1, strip.white=TRUE, colClasses="character"))
    ## Unintuitively, the connection object seems to have been destroyed
    ## by the previous read.fwf.  We need to create a new one.
    con <- file(fname, encoding = encoding)
    ## If columns 3 in chron.stats is an integer then there is no
    ## statistics line
    if((!is.int(chron.stats[[3]]) | (chron.stats[[3]]==0)) & !grepl(" ",chron.stats[[3]])){
      names(chron.stats) <-
        c("SiteID", "nYears", "AC[1]", "StdDev", "MeanSens",
          "MeanRWI", "IndicesSum", "IndicesSS", "MaxSeries")
      #       cat(gettext("Embedded chronology statistics\n", domain="R-dplR"))
      #       print(chron.stats)
      ## Really read file
      dat <- read.fwf(con, c(6, digits.year, rep(c(4, 3), 10)),
                      skip=skip.lines, n=nlines-skip.lines-1,
                      strip.white=TRUE, colClasses="character")
    } else {
      ## Really read file
      dat <- read.fwf(con, c(6, digits.year, rep(c(4, 3), 10)),
                      skip=skip.lines, n=nlines-skip.lines,
                      strip.white=TRUE, colClasses="character")
    }
    
    #     dat <- dat[complete.cases(dat),]
    
    dat[[1]] <- as.character(dat[[1]])
    for(i in 2:22){
      dat[[i]] <- as.numeric(as.character(dat[[i]]))
    }
    
    ## Remove any blank lines at the end of the file, for instance
    dat <- dat[!is.na(dat[[2]]), , drop=FALSE] # requires non-NA year
    
    series.label <- dat[[1]]
    series.ids <- unique(series.label)
    decade.yr <- as.numeric(dat[[2]])
    nseries <- length(series.ids)
    
    series.index <- match(series.label, series.ids)
    min.year <- (min(decade.yr) %/% 10) * 10
    max.year <- ((max(decade.yr)+10) %/% 10) * 10
    span <- max.year - min.year + 1
    ncol.crn.mat <- nseries + 1
    crn.mat <- matrix(NA_real_, ncol=ncol.crn.mat, nrow=span)
    colnames(crn.mat) <- c(as.character(series.ids), "samp.depth")
    rownames(crn.mat) <- min.year:max.year
    ## RWI
    x <- as.matrix(dat[seq(from=3, to=21, by=2)])
    ## All sample depths
    y <- as.matrix(dat[seq(from=4, to=22, by=2)])
    for(i in seq_len(nseries)){
      idx <- which(series.index == i)
      for(j in idx) {
        yr <- (decade.yr[j] %/% 10) * 10
        row.seq <- seq(from = yr - min.year + 1, by = 1, length.out = 10)
        crn.mat[row.seq, i] <- x[j, ]
        if(i == 1) {
          crn.mat[row.seq, ncol.crn.mat] <- y[j, ]
        }
      }
    }
    ## Clean up NAs
    crn.mat[which(crn.mat[, -ncol.crn.mat] == 9990)] <- NA # column-major order
    crn.mat <-
      crn.mat[!apply(is.na(crn.mat[, -ncol.crn.mat, drop=FALSE]), 1, all),
              ,
              drop=FALSE]
    ## If samp depth is all 1 then dump it
    sd.one <- all(crn.mat[, ncol.crn.mat] == 1)
    if(!is.na(sd.one) & sd.one) {
      save.names <- colnames(crn.mat)[-ncol.crn.mat]
      crn.mat <- crn.mat[, -ncol.crn.mat, drop=FALSE]
      crn.mat <- crn.mat / 1000
      crn.df <- as.data.frame(crn.mat)
      names(crn.df) <- save.names
      #       cat(gettext("All embedded sample depths are one...Dumping from matrix\n",
      #                   domain="R-dplR"))
    }else {
      seq.series <- seq_len(nseries)
      crn.mat[, seq.series] <- crn.mat[, seq.series] / 1000
      crn.df <- as.data.frame(crn.mat)
    }
    output <- crn.df
    
    output$YEAR <- row.names(output)
    if(ncol(output)==2){
      output <- output[,c(2,1)]
    }else if(ncol(output)==3){
      output <- output[,c(3,1)]
    }else{
      return("ERROR: IMPROPERLY FORMATTED CRN")
    }
    
    names(output) <- c("YEAR",name)
    
    out.mean <- mean(output[,2], na.rm=T)
    output[,2] <- output[,2]/out.mean
    
    final.out <- list(meta=meta,data=output)
    
    return(final.out)
    
  }else if(tail(scan(series,what='character',nlines=1, quiet=T),n=1)=="RAW"){
    raw.data <- scan(series,what='character', multi.line=F, fill=T, sep="\n", quiet=T)
    raw.data <- sub("s\x9fdl.","     ",raw.data)
    raw.data <- raw.data[grepl("*STD",raw.data)]
    
    
    meta <- raw.data[1:3]
    data <- raw.data[-c(1:3)]
    data <- sub("  STD","",data)
    
    
    meta.1 <- unlist(strsplit(meta[1], "\\s+"))[c(-1,-2)]
    place <- paste(meta.1[1:which(meta.1=="WIDTH_RING")-1], collapse=' ')
    type <- meta.1[which(meta.1=="WIDTH_RING")+1]
    
    meta.2 <- unlist(strsplit(meta[2], "\\s+"))[c(-1,-2)]
    meta.2 <- meta.2[-length(meta.2)]
    if(length(which(meta.2=="-"))>0){
      meta.2 <- meta.2[-which(meta.2=="-")]
    }
    title <- paste(meta.2[1:(length(meta.2)-4)], collapse=' ')
    elev <- meta.2[(length(meta.2)-3)]
    location <- meta.2[(length(meta.2)-2)]
    location.split <- unlist(strsplit(location,""))
    splits <- which(location.split=="-" | location.split==" ")-1
    if(splits[1]==0){
      splits <- splits[-1]
    }
    lat <- paste(location.split[1:splits[1]],collapse='')
    lon <- paste(location.split[(splits[1]+1):length(location.split)],collapse='')
    begin <- meta.2[(length(meta.2)-1)]
    end <- meta.2[(length(meta.2))]
    
    
    # Parse the header of the CRN file
    #     widths <- c(9, 21, 11, 4, 35)
    #     starts <- c(1,cumsum(widths)+1)
    #     stops <- cumsum(widths)
    #     starts <- head(starts,n=length(stops))
    #     meta.1 <- sapply(1:length(starts), function(ii){substr(meta[1], starts[ii], stops[ii])})[c(2,4)]
    #     meta.1 <- gsub("^\\s+|\\s+$", "", meta.1)
    
    #     widths <- c(9, 24, 15, 4, 6, 6, 5, 5, 3, 3)
    #     starts <- c(1,cumsum(widths)+1)
    #     stops <- cumsum(widths)
    #     starts <- head(starts,n=length(stops))
    #     meta.2 <- sapply(1:length(starts), function(ii){substr(meta[2], starts[ii], stops[ii])})[c(2,4:8)]
    #     meta.2 <- gsub("^\\s+|\\s+$", "", meta.2)
    
    widths <- c(9, 71)
    starts <- c(1,cumsum(widths)+1)
    stops <- cumsum(widths)
    starts <- head(starts,n=length(stops))
    meta.3 <- sapply(1:length(starts), function(ii){substr(meta[3], starts[ii], stops[ii])})[2]
    meta.3 <- sub("-","",meta.3)
    meta.3 <- sub("STD","",meta.3)
    meta.3 <- gsub("^\\s+|\\s+$", "", meta.3)
    
    # Get the name of the 
    name <- gsub(".crn","",series)
    name <- gsub(paste(output.dir,sep=''),"",name)
    
    meta <- c(name,place,type,elev,lat,lon,begin,end,meta.3)
    
    digits.year <- max(nchar(meta[7]),nchar(meta[8]),4)
    
    widths <- c(6, digits.year, rep(c(4, 3), 10))
    starts <- c(1,cumsum(widths)+1)
    stops <- cumsum(widths)
    starts <- head(starts,n=length(stops))
    dat <- as.data.frame(matrix(unlist(lapply(data,FUN=function(x){sapply(1:length(starts), function(ii){substr(x, starts[ii], stops[ii])})})),ncol=22,byrow=T))
    
    dat[[1]] <- as.character(dat[[1]])
    for(i in 2:22){
      dat[[i]] <- as.numeric(as.character(dat[[i]]))
    }
    
    series.name <- dat[[1]]
    series.ids <- unique(series.name)
    decade.yr <- dat[[2]]
    nseries <- length(series.ids)
    series.index <- match(series.name, series.ids)
    min.year <- (min(decade.yr) %/% 10) * 10
    max.year <- ((max(decade.yr)+10) %/% 10) * 10
    span <- max.year - min.year + 1
    ncol.crn.mat <- nseries + 1
    crn.mat <- matrix(NA_real_, ncol=ncol.crn.mat, nrow=span)
    colnames(crn.mat) <- c(as.character(series.ids), "samp.depth")
    rownames(crn.mat) <- min.year:max.year
    ## RWI
    x <- as.matrix(dat[seq(from=3, to=21, by=2)])
    ## All sample depths
    y <- as.matrix(dat[seq(from=4, to=22, by=2)])
    for(i in seq_len(nseries)){
      idx <- which(series.index == i)
      for(j in idx) {
        yr <- (decade.yr[j] %/% 10) * 10
        row.seq <- seq(from = yr - min.year + 1, by = 1, length.out = 10)
        crn.mat[row.seq, i] <- x[j, ]
        if(i == 1) {
          crn.mat[row.seq, ncol.crn.mat] <- y[j, ]
        }
      }
    }
    ## Clean up NAs
    crn.mat[which(crn.mat[, -ncol.crn.mat] == 9990)] <- NA # column-major order
    crn.mat <-
      crn.mat[!apply(is.na(crn.mat[, -ncol.crn.mat, drop=FALSE]), 1, all),
              ,
              drop=FALSE]
    
    seq.series <- seq_len(nseries)
    crn.mat[, seq.series] <- crn.mat[, seq.series] / 1000
    crn.df <- as.data.frame(crn.mat)
    
    output <- crn.df
    
    output$YEAR <- row.names(output)
    if(ncol(output)==2){
      output <- output[,c(2,1)]
    }else if(ncol(output)==3){
      output <- output[,c(3,1)]
    }else{
      return("ERROR: IMPROPERLY FORMATTED CRN")
    }
    
    names(output) <- c("YEAR",name)
    
    out.mean <- mean(output[,2], na.rm=T)
    output[,2] <- output[,2]/out.mean
    
    final.out <- list(meta=meta,data=output)
    
    return(final.out)
    
  }else{
    return(paste("ERROR! IMPROPERLY FORMATTED CRN: ",series,sep=''))
  }
  
}

### Function to check if x is equivalent to its integer
### representation. Note: Returns FALSE for values that fall outside
### the range of the integer type. The result has the same shape as x;
### at least vector and array x are supported.
is.int <- function(x) {
  suppressWarnings(y <- as.numeric(x) == as.integer(as.numeric(x)))
  y[is.na(y)] <- FALSE
  y
}

getITRDB_MANN08 <- function(data.dir){
  if(!file.exists(paste(data.dir,"mann2008datainfilled.txt",sep=''))){
    system(paste("wget -nd --output-document=",paste(data.dir,"mann2008datainfilled.txt",sep='')," ftp://ftp.ncdc.noaa.gov/pub/data/paleo/reconstructions/pcn/proxy/mann2008/mann2008datainfilled.txt",sep=''))
  }
  
  if(!file.exists(paste(data.dir,"1209proxynames.txt",sep=''))){
    system(paste("wget -nd --output-document=",paste(data.dir,"1209proxynames.txt",sep='')," ftp://ftp.ncdc.noaa.gov/pub/data/paleo/reconstructions/pcn/proxy/mann2008/1209proxynames.txt",sep=''))
  }
  
  mann2008.data <- read.csv(paste(data.dir,"mann2008datainfilled.txt",sep=''), header=T, skip=25)
  mann2008.meta <- read.delim(paste(data.dir,"1209proxynames.txt",sep=''),header=T)
  names(mann2008.data)[1] <- "YEAR"
  
  mann2008.meta <- mann2008.meta[mann2008.meta$Data.type.code %in% c(9000,3000),c(26,8,6,7,9)]
  names(mann2008.meta) <- c("NAME","ELEVATION","LATITUDE","LONGITUDE","SPECIES")
  
  mann2008.data <- mann2008.data[,names(mann2008.data) %in% c("YEAR",as.character(mann2008.meta$NAME))]
  mann2008.data <- data.matrix(mann2008.data)
  mann2008.data[mann2008.data=="NaN"] <- NA
  mann2008.data <- mann2008.data[c(-8:-1),]
  
  write.csv(mann2008.meta,paste(data.dir,"mann2008_meta.csv",sep=''), row.names=F)
  
  # Standardize columns
  years <- as.data.frame(mann2008.data[,1])
  data <- as.data.frame(scale(mann2008.data[,-1]))
  mann2008.data <- cbind(years,data)
  names(mann2008.data)[1] <- "YEAR"
  
  return(as.data.frame(mann2008.data))
  
}

getITRDBMetadataSPDF <- function(series.names, data.dir){
  dendros.study <- read.csv(paste(data.dir,"ITRDB_METADATA.csv",sep=''))
  dendros.study <- dendros.study[!is.na(dendros.study$LON) & (dendros.study$SERIES %in% series.names),]
  dendros.study <- dendros.study[order(dendros.study$SERIES),]
  
  dendros.study <- SpatialPointsDataFrame(dendros.study[,c("LON","LAT")],dendros.study, proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  suppressWarnings(writeOGR(dendros.study,paste(data.dir,"ITRDB_METADATA.kml",sep=''),"ITRDB_METADATA",driver="KML", overwrite_layer=TRUE))
  
  return(dendros.study)
}


getITRDBMetadataSPDF_MANN08 <- function(data.dir){
  dendros.study <- read.csv(paste(data.dir,"mann2008_meta.csv",sep=''))
  dendros.study <- dendros.study[!is.na(dendros.study$LONGITUDE),]
  dendros.study <- dendros.study[order(dendros.study$NAME),]
  
  dendros.study <- SpatialPointsDataFrame(dendros.study[,c("LONGITUDE","LATITUDE")],dendros.study, proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  suppressWarnings(writeOGR(dendros.study,paste(data.dir,"mann2008_meta.kml",sep=''),"mann2008_meta",driver="KML", overwrite_layer=TRUE))
  
  return(dendros.study)
}

sanitizeITRDBnames <- function(names){
  names <- gsub("NATIONAL MONUMENT","NM",names)
  names <- gsub("-.*","",names)
  names <- gsub("SITE.*","",names)
  names <- gsub("[(].*","",names)
  names <- gsub("[*].*","",names)
  names <- gsub("PIPO.*","",names)
  names <- gsub("  .*","",names)
  names <- gsub("[+].*","",names)
  names <- gsub("RES[.]","RESERVOIR",names)
  names <- gsub("[[:space:]]+$","",names)

  return(names)
}
