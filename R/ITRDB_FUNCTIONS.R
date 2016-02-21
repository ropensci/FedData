#' Download the latest version of the ITRDB, and extract given parameters.
#'
#' \code{get_itrdb} returns a named list of length 3: 
#' \enumerate{
#' \item "metadata": A data.table or \code{SpatialPointsDataFrame} (if \code{makeSpatial==TRUE}) of the locations 
#' and names of extracted ITRDB chrononlogies,
#' \item "widths": A matrix of tree-ring widths/densities given user selection, and
#' \item "depths": A matrix of tree-ring sample depths.
#' }
#' 
#' @param template A Raster* or Spatial* object to serve 
#' as a template for selecting chronologies. If missing, 
#' all available global chronologies are returned.
#' @param label A character string naming the study area.
#' @param recon.years A numeric vector of years over which reconstructions are needed; 
#' if missing, the union of all years in the available chronologies are given.
#' @param calib.years A numeric vector of all required years---chronologies without these years will be discarded; 
#' if missing, all available chronologies are given.
#' @param species A character vector of 4-letter tree species identifiers; 
#' if missing, all available chronologies are given.
#' @param measurement.type A character vector of measurement type identifiers. Options include:
#' \itemize{
#' \item "Total Ring Density"
#' \item "Earlywood Width"
#' \item "Earlywood Density"
#' \item "Latewood Width"
#' \item "Minimum Density"
#' \item "Ring Width"
#' \item "Latewood Density"
#' \item "Maximum Density"
#' \item "Latewood Percent"
#' }
#' if missing, all available chronologies are given.
#' @param chronology.type A character vector of chronology type identifiers. Options include:
#' \itemize{
#' \item "ARSTND"
#' \item "Low Pass Filter"
#' \item "Residual"
#' \item "Standard"
#' \item "Re-Whitened Residual"
#' \item "Measurements Only"
#' }
#' if missing, all available chronologies are given.
#' @param makeSpatial Should the metadata be presented as a SpatialPointsDataFrame? Defaults to FALSE.
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing. Defaults to "./RAW/ITRDB/".
#' @param extraction.dir A character string indicating where the extracted and cropped ITRDB dataset should be put.
#' The directory will be created if missing. Defaults to "./EXTRACTIONS/ITRDB/".
#' @param force.redo If an extraction already exists, should a new one be created? Defaults to FALSE.
#' @return A named list containing the "metadata", "widths", and "depths" data. 
#' @export
get_itrdb <- function(template=NULL, label=NULL, recon.years=NULL, calib.years=NULL, species=NULL, measurement.type=NULL, chronology.type=NULL, makeSpatial=F, raw.dir="./RAW/ITRDB/", extraction.dir="./EXTRACTIONS/ITRDB/", force.redo=FALSE){  
  dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
  
  if(is.null(template) & is.null(label)){
    label <- "allChronologies"
  }
  
  if(!is.null(template) & is.null(label)){
    stop("Template provided but no label given.")
  }
  
  dir.create(paste(extraction.dir,'/',label,'/',sep=''), showWarnings = FALSE, recursive = TRUE)
  
  if(!force.redo && !is.null(label) && file.exists(paste(extraction.dir,'/',label,'/',label,'.Rds',sep=''))){
    out <- readRDS(paste(extraction.dir,'/',label,'/',label,'.Rds',sep=''))
    return(out)
  }
  
  data <- download_itrdb(raw.dir=raw.dir,force.redo=force.redo)
  
  ## Nulling out to appease R CMD CHECK
  LAT <- LON <- START <- END <- SPECIES <- MEASUREMENT_TYPE <- CHRONOLOGY_TYPE <- NULL
  
  if(!is.null(calib.years)){
    data <- data[START < min(calib.years) & END > max(calib.years)]
  }
  
  if(!is.null(species)){
    data <- data[SPECIES %in% species]
  }
  
  if(!is.null(measurement.type)){
    data <- data[MEASUREMENT_TYPE %in% measurement.type]
  }
  
  if(!is.null(chronology.type)){
    data <- data[CHRONOLOGY_TYPE %in% chronology.type]
  }
  
  if(!is.null(template)){
    data <- data[LAT >= -90 & LAT <= 90 & LON >= -180 & LON <= 180]
    sp.data <- sp::SpatialPointsDataFrame(as.matrix(data[,c("LON","LAT"),with=F]),data=data[,"SERIES", with=F, drop=F], proj4string=sp::CRS("+proj=longlat +datum=WGS84"))
    data <- data[!is.na((sp.data %over% sp::spTransform(template,sp::CRS(raster::projection(sp.data)))))]
    rm(sp.data)
    if(dim(data)[[1]]==0) stop("No ITRDB chronologies within template polygon.")
  }
  
  if(is.null(recon.years)){
    if(dim(data)[[1]]==0) stop("No ITRDB chronologies within template polygon.")
    all.years <- range(as.numeric(unique(unlist(sapply(data[,data],rownames)))))
    recon.years <- all.years[[1]]:all.years[[2]]
  }
  
  all.data <- lapply(data[,data],function(df){
    df <- df[rownames(df) %in% recon.years,]
    df.years <- as.numeric(rownames(df))
    missing.years <- setdiff(recon.years,df.years)
    missing.mat <- matrix(nrow=length(missing.years),ncol=2)
    colnames(missing.mat) <- c("WIDTH","DEPTH")
    rownames(missing.mat) <- missing.years
    df.all <- rbind(df,missing.mat)
    df.all <- df.all[order(as.numeric(rownames(df.all))),]
  })
  names(all.data) <- data$SERIES
  
  widths <- sapply(all.data,function(df){data.matrix(df[,"WIDTH",drop=F])})
  depths <- sapply(all.data,function(df){data.matrix(df[,"DEPTH",drop=F])})
  
  colnames(widths) <- colnames(depths) <- names(all.data)
  rownames(widths) <- rownames(depths) <- recon.years
  
  data <- data[,data:=NULL]
  if(makeSpatial){
    data <- sp::SpatialPointsDataFrame(as.matrix(data[,c("LON","LAT"),with=F]),data=data, proj4string=sp::CRS("+proj=longlat +datum=WGS84"))
  }
  
  out <- list(metadata=data,widths=widths,depths=depths)
  if(!is.null(label)){
    saveRDS(out,file=paste(extraction.dir,'/',label,'/',label,'.Rds',sep=''), compress='xz')
  }
  
  return(out)
}

#' Download the latest version of the ITRDB.
#' 
#' Downloads and parses the latest zipped (numbered) version of the ITRDB.
#' This function includes improvements to the \code{\link{read_crn}} function from the 
#' \pkg{dplR} library. The principle changes are better parsing of metadata, and support
#' for the Schweingruber-type Tucson format. Chronologies that are unable to be read
#' are reported to the user.
#' 
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing. Defaults to "./RAW/ITRDB/".
#' @param force.redo If a download already exists, should a new one be created? Defaults to FALSE.
#' @return A data.table containing all of the ITRDB data. 
#' @export
download_itrdb <- function(raw.dir="./RAW/ITRDB/", force.redo=FALSE){
  dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
  
  opts <- list(
    verbose = F, 
    noprogress = T,
    fresh_connect = T, 
    ftp_use_epsv = T, 
    forbid_reuse = T,
    dirlistonly = T)
  hand <- curl::new_handle()
  curl::handle_setopt(hand, .list = opts)
  
  url <- 'ftp://ftp.ncdc.noaa.gov/pub/data/paleo/treering/chronologies/'
  filenames = readLines(curl::curl(url, handle = hand))
  filenames = paste(url, filenames, sep = "")
  filenames <- filenames[grep("*.zip",filenames)]
  
  for(file in filenames){
    download_data(url=file,destdir=raw.dir)
  }
  
  ## A vector of the files in the output.dir
  zips <- paste0(raw.dir,basename(filenames))
  
  version <- max(as.numeric(gsub("[^0-9]", "",zips)))
  zips <- zips[grepl(version,zips)]
  
  message("Extracting chronology data from ITRDB version ",version," dated ",as.character(as.Date(base::file.info(zips[[1]])$mtime)))
  
  if(!force.redo && 
       file.exists(paste(raw.dir,"ITRDB",version,".Rds",sep=''))){
    itrdb.metadata <- readRDS(paste(raw.dir,"/ITRDB",version,".Rds",sep=''))
  }else{
    all.data <- lapply(zips,function(file){
      
      tmpdir <- tempfile()
      if (!dir.create(tmpdir)) stop("failed to create my temporary directory")
      
      utils::unzip(file,exdir=tmpdir)
      
      crns <- list.files(tmpdir, full.names=T)
      crns <- crns[grepl("\\.crn",crns)]
      
      message("Extracting chronology data from ",length(crns)," files in ",file)
      
      records <- lapply(crns,function(this.crn){
        tryCatch(suppressWarnings(read_crn(this.crn)),error=function(e){return(NULL)})
      })
      
      unlink(tmpdir, recursive = TRUE)
      
      if(sum(sapply(records, is.null))>0){
        message(sum(sapply(records, is.null))," file(s) couldn't be extracted:")
        for(file in basename(crns)[sapply(records, is.null)])
          message(file)
        message("File(s) likely malformed!")
      }
      
      records <- unlist(records,recursive=F)
      
      return(records)
    })
    
    all.data <- unlist(all.data, recursive=F)
    
    ## Nulling out to appease R CMD CHECK
    LAT <- LON <- START <- END <- data <- NULL
    
    itrdb.metadata <- data.table::rbindlist(lapply(all.data,'[[','meta'))
    itrdb.data <- lapply(all.data,'[[','data')
    itrdb.metadata[,data:=itrdb.data]
    
    # Change latitude and longitude to numeric
    itrdb.metadata[,LAT:=as.numeric(as.character(LAT))]
    itrdb.metadata[,LON:=as.numeric(as.character(LON))]
    
    # remove records with unknown locations
    itrdb.metadata <- itrdb.metadata[!is.na(LAT) & !is.na(LON),]
    
    # or unknown years
    itrdb.metadata <- itrdb.metadata[!is.na(START) & !is.na(END),]
    
    # or nonsense years (There are ~50 weird ones)
    itrdb.metadata <- itrdb.metadata[START <= as.numeric(format(Sys.time(), "%Y")) & END <= as.numeric(format(Sys.time(), "%Y"))]
    
    saveRDS(itrdb.metadata,paste(raw.dir,"/ITRDB",version,".Rds",sep=''), compress="xz")
  }
  
  return(itrdb.metadata)
}

#' Read a Tucson-format chronology file.
#' 
#' This function includes improvements to the \code{read.crn} function from the 
#' \pkg{dplR} library. The principle changes are better parsing of metadata, and support
#' for the Schweingruber-type Tucson format. Chronologies that are unable to be read
#' are reported to the user. This function automatically recognizes Schweingruber-type files.
#' 
#' This wraps two other functions: \code{\link{read_crn_metadata}} \code{\link{read_crn_data}}.
#' 
#' @param file A character string path pointing to a \code{*.crn} file to be read.
#' @return A list containing the metadata and chronology.
#' @export
read_crn <- function(file){
  id <- toupper(gsub(".crn","",basename(file)))
  
  all.data <- scan(file,what='character', multi.line=F, fill=T, sep="\n", quiet=T)
  all.data <- iconv(all.data,from="latin1",to="")
  
  # This removes blank lines
  con <- file(file)
  writeLines(all.data,con)
  close(con)
  
  tails <- substr_right(all.data,3)
  if(any(tails=="RAW")){
    SCHWEINGRUBER <- T
  }else{
    SCHWEINGRUBER <- F
  }
  
  meta <- read_crn_metadata(file,SCHWEINGRUBER)
  data <- read_crn_data(file,SCHWEINGRUBER)
  
  if(!SCHWEINGRUBER){
    year.range <- range(as.numeric(rownames(data)))
  }else{
    year.range <- range(as.numeric(rownames(data[[1]])))
  }
  
  meta$START <- year.range[[1]]
  meta$END <- year.range[[2]]
  
  meta$FILE <- basename(file)
  
  if(SCHWEINGRUBER){
    out <- lapply(1:length(data), function(i){      
      if(names(data)[i]=='ARS') type.chronology <- "ARSTND"
      if(names(data)[i]=='RAW') type.chronology <- "Measurements Only"
      if(names(data)[i]=='RES') type.chronology <- "Residual"
      if(names(data)[i]=='STD') type.chronology <- "Standard"
      this.meta <- meta
      
      this.meta$CHRONOLOGY_TYPE <- type.chronology
      
      id <- regmatches(this.meta$SERIES,regexec(".*[0-9]", as.character(this.meta$SERIES)))
      typeID.chronology <- NA
      typeID.measurement <- NA
      
      if(this.meta$CHRONOLOGY_TYPE == "ARSTND") typeID.chronology <- "A"
      if(this.meta$CHRONOLOGY_TYPE == "Low Pass Filter") typeID.chronology <- "P"
      if(this.meta$CHRONOLOGY_TYPE == "Residual") typeID.chronology <- "R"
      if(this.meta$CHRONOLOGY_TYPE == "Standard") typeID.chronology <- "S"
      if(this.meta$CHRONOLOGY_TYPE == "Re-Whitened Residual") typeID.chronology <- "W"
      if(this.meta$CHRONOLOGY_TYPE == "Measurements Only") typeID.chronology <- "N"
      
      if(this.meta$MEASUREMENT_TYPE == "Total Ring Density") typeID.measurement <- "D"
      if(this.meta$MEASUREMENT_TYPE == "Earlywood Width") typeID.measurement <- "E"
      if(this.meta$MEASUREMENT_TYPE == "Earlywood Density") typeID.measurement <- "I"
      if(this.meta$MEASUREMENT_TYPE == "Latewood Width") typeID.measurement <- "L"
      if(this.meta$MEASUREMENT_TYPE == "Minimum Density") typeID.measurement <- "N"
      if(this.meta$MEASUREMENT_TYPE == "Ring Width") typeID.measurement <- "R"
      if(this.meta$MEASUREMENT_TYPE == "Latewood Density") typeID.measurement <- "T"
      if(this.meta$MEASUREMENT_TYPE == "Maximum Density") typeID.measurement <- "X"
      if(this.meta$MEASUREMENT_TYPE == "Latewood Percent") typeID.measurement <- "P"
      
      this.meta$SERIES <- paste(id,typeID.measurement,typeID.chronology,sep='')
      
      return(list(meta=this.meta,data=data[[i]]))
    })
    return(out)
  }
  
  return(list(list(meta=meta,data=data)))
}

#' Read metadata from a Tucson-format chronology file.
#' 
#' This function includes improvements to the \code{\link{read_crn}} function from the 
#' \pkg{dplR} library. The principle changes are better parsing of metadata, and support
#' for the Schweingruber-type Tucson format. Chronologies that are unable to be read
#' are reported to the user. The user (or \code{\link{read_crn}}) must tell the function whether
#' the file is a Schweingruber-type chronology.
#' 
#' Location information is converted to decimal degrees.
#' 
#' @param file A character string path pointing to a \code{*.crn} file to be read.
#' @param SCHWEINGRUBER Is the file in the Schweingruber-type Tucson format?
#' @return A data.frame containing the metadata.
#' @export
read_crn_metadata <- function(file,SCHWEINGRUBER){
  id <- toupper(gsub("\\.crn","",basename(file)))
  id <- gsub("_CRNS","",id)
  typeID <- regmatches(id,regexec("[0-9]+(.*)", id))[[1]][2]
  id <- regmatches(id,regexec(".*[0-9]", id))
  typeID.chronology <- NA
  typeID.measurement <- NA
  if(nchar(typeID)==2){
    typeID <- strsplit(typeID,split='')[[1]]
    typeID.measurement <- typeID[[1]]
    typeID.chronology <- typeID[[2]]
  }else if(SCHWEINGRUBER){
    typeID.measurement <- typeID[[1]]
  }else{
    typeID.chronology <- typeID[[1]]
  }
  if(is.na(typeID.chronology)) typeID.chronology <- ''
  if(is.na(typeID.measurement)) typeID.measurement <- ''
  
  type.chronology <- NA
  type.measurement <- NA
  
  if(typeID.chronology=="A") type.chronology <- "ARSTND"
  if(typeID.chronology=="P") type.chronology <- "Low Pass Filter"
  if(typeID.chronology=="R") type.chronology <- "Residual"
  if(typeID.chronology=="S") type.chronology <- "Standard"
  if(typeID.chronology=="W") type.chronology <- "Re-Whitened Residual"
  if(typeID.chronology=="N") type.chronology <- "Measurements Only"
  if(typeID.chronology==''){
    type.chronology <- "Standard"
    typeID.chronology <- "S"
  }
  if(is.na(type.chronology)){
    type.chronology <- "Standard"
    typeID.measurement <- typeID.chronology
    typeID.chronology <- "S"
  }
  
  if(typeID.measurement=="D") type.measurement <- "Total Ring Density"
  if(typeID.measurement=="E") type.measurement <- "Earlywood Width"
  if(typeID.measurement=="I") type.measurement <- "Earlywood Density"
  if(typeID.measurement=="L") type.measurement <- "Latewood Width"
  if(typeID.measurement=="N") type.measurement <- "Minimum Density"
  if(typeID.measurement=="R") type.measurement <- "Ring Width"
  if(typeID.measurement=="T") type.measurement <- "Latewood Density"
  if(typeID.measurement=="X") type.measurement <- "Maximum Density"
  if(typeID.measurement=="P") type.measurement <- "Latewood Percent"
  if(typeID.measurement=="W") type.measurement <- "Ring Width"
  if(typeID.measurement==''){
    type.measurement <- "Ring Width"
    typeID.measurement <- "R"
  }
  
  id <- paste(id,typeID.measurement,typeID.chronology,sep='')
  
  ## Nulling out to appease R CMD CHECK
  lats <- lons <- NULL
  
  if(!SCHWEINGRUBER){
    # Parse the header of the CRN file
    meta.1 <- as.character(utils::read.fwf(file,c(6,3,52,4), skip=0, n=1, colClasses = "character", strip.white=T, stringsAsFactors=F))
    meta.2 <- as.character(utils::read.fwf(file,c(6,3,13,18,6,5,6,9,6,5), skip=1,n=1, colClasses = "character", strip.white=T, stringsAsFactors=F))
    meta.3 <- as.character(utils::read.fwf(file,c(6,3,52,2,12), skip=2, n=1, colClasses = "character", strip.white=T, stringsAsFactors=F))
    
    meta <- c(id,meta.1[3:4],type.measurement,type.chronology,meta.2[c(5:7,9:10)],meta.3[3])
    
    names(meta) <- c("SERIES","NAME","SPECIES","MEASUREMENT_TYPE","CHRONOLOGY_TYPE","ELEVATION","LAT","LON","START","END","CONTRIBUTOR")
    
  }else{
    meta <- scan(file,what='character', n=3,multi.line=F, fill=T, sep="\n", quiet=T)
    meta <- sub("  RAW","",meta)
    
    meta[1] <- gsub(substr(meta[1],1,9),"",meta[1])
    meta[2] <- gsub(substr(meta[2],1,9),"",meta[2])
    meta[3] <- gsub(substr(meta[3],1,9),"",meta[3])
    
    
    
    meta.1 <- unlist(strsplit(meta[1], "\\s+"))
    place <- paste(meta.1[1:which(grepl("_",meta.1))-1], collapse=' ')
    
    species <- meta.1[which(grepl("_",meta.1))+1]
    
    meta.2 <- unlist(strsplit(meta[2], "\\s+\\s+\\s+"))
    meta.2 <- c(meta.2[1],substr(meta.2[2],1,as.numeric(regexpr('\\d',meta.2[2]))-1),substr(meta.2[2],as.numeric(regexpr('\\d',meta.2[2])),nchar(meta.2[2])))
    meta.2 <- c(meta.2[1],meta.2[2],unlist(strsplit(meta.2[3], "\\s+")))
    if(length(which(meta.2=="-"))>0){
      meta.2 <- meta.2[-which(meta.2=="-")]
    }
    title <- meta.2[1:2]
    elev <- meta.2[3]
    location <- meta.2[4]
    location.split <- unlist(strsplit(location,""))
    splits <- which(location.split=="-" | location.split==" ")-1
    if(length(splits)>0 && splits[1]==0){
      splits <- splits[-1]
    }
    lat <- tryCatch(paste(location.split[1:splits[1]],collapse=''),error=function(e){NA})
    lon <- tryCatch(paste(location.split[(splits[1]+1):length(location.split)],collapse=''),error=function(e){NA})
    begin <- meta.2[5]
    end <- meta.2[6]
    
    meta.3 <- meta[3]
    meta.3 <- sub("\\s-.*","",meta.3)
    
    meta <- c(id,place,species,type.measurement,type.chronology,elev,lat,lon,begin,end,meta.3)
    
    names(meta) <- c("SERIES","NAME","SPECIES","MEASUREMENT_TYPE","CHRONOLOGY_TYPE","ELEVATION","LAT","LON","START","END","CONTRIBUTOR")
    
  }
  
  meta[["ELEVATION"]] <- gsub("M", "", meta[["ELEVATION"]])
  meta[["ELEVATION"]] <- gsub("m", "", meta[["ELEVATION"]])
  meta[["CONTRIBUTOR"]] <- gsub("  ", " ", meta[["CONTRIBUTOR"]])
  meta[["LAT"]] <- gsub(" ", "0", meta[["LAT"]],fixed = TRUE)
  meta[["LON"]] <- gsub(" ", "0", meta[["LON"]],fixed = TRUE)
  meta[["LAT"]] <- gsub("+", "", meta[["LAT"]],fixed = TRUE)
  #   meta[["LON"]] <- gsub("-", "", meta[["LON"]],fixed = TRUE)
  meta[["CONTRIBUTOR"]] <- toupper(meta[["CONTRIBUTOR"]])
  meta[["NAME"]] <- toupper(meta[["NAME"]])
  meta[["SERIES"]] <- toupper(meta[["SERIES"]])
  
  meta[["LAT"]] <- format(as.numeric(meta[["LAT"]])/100, nsmall=2)
  meta[["LON"]] <- format(as.numeric(meta[["LON"]])/100, nsmall=2)
  
  if(!any(meta[["LAT"]]=="NA",meta[["LON"]]=="NA")){
    locations <- meta[c("LAT","LON")]
    
    locations[["LAT"]] <- within(data.frame(lats=locations[["LAT"]]), {
      dms <- do.call(rbind, strsplit(as.character(lats), "\\."))
      neg <- grepl("-",dms[,1])
      LAT <- abs(as.numeric(dms[,1])) + (as.numeric(dms[,2]))/60
      if(neg){
        LAT <- -1*LAT
      }
      rm(dms)
      rm(lats)
      rm(neg)
    })
    
    locations[["LON"]] <- within(data.frame(lons=locations[["LON"]]), {
      dms <- do.call(rbind, strsplit(as.character(lons), "\\."))
      neg <- grepl("-",dms[,1])
      LON <- abs(as.numeric(dms[,1])) + (as.numeric(dms[,2]))/60
      if(neg){
        LON <- -1*LON
      }
      rm(dms)
      rm(lons)
      rm(neg)
    })
    
    meta[["LAT"]] <- locations[["LAT"]]
    meta[["LON"]] <- locations[["LON"]]
    
  }
  
  meta <- data.frame(meta)
  if(nrow(meta)>ncol(meta)) meta <- data.frame(t(meta))
  
  return(meta)
  
}

#' Read chronology data from a Tucson-format chronology file.
#' 
#' This function includes improvements to the \code{\link{read_crn}} function from the 
#' \pkg{dplR} library. The principle changes are better parsing of metadata, and support
#' for the Schweingruber-type Tucson format. Chronologies that are unable to be read
#' are reported to the user. The user (or \code{\link{read_crn}}) must tell the function whether
#' the file is a Schweingruber-type chronology.
#' 
#' @param file A character string path pointing to a \code{*.crn} file to be read.
#' @param SCHWEINGRUBER Is the file in the Schweingruber-type Tucson format?
#' @return A data.frame containing the data, or if \code{SCHWEINGRUBER==T}, a list containing four types of data.
#' @export
read_crn_data <- function(file,SCHWEINGRUBER){  
  if(!SCHWEINGRUBER){
    years <- as.character(utils::read.fwf(file,c(6,3,13,18,6,5,6,9,6,5),skip=1,n=1,colClasses = "character", strip.white=T, stringsAsFactors=F))[9:10]
    
    if(any(grepl(" ",years))){
      years <- unlist(strsplit(years," "))
    }
    
    digits.year <- max(nchar(years),4)
    dig.dif <- digits.year-nchar(years[[1]])
    
    encoding = getOption("encoding")
    ## Open the data file for reading
    con <- file(file, encoding = encoding)
    on.exit(close(con))
    
    ## Read 4th line - should be first data line
    dat1 <- readLines(con, n=4)
    if(length(dat1) < 4)
      stop("file has under 4 lines")
    dat1 <- dat1[4]
    
    #     yearStart <- nchar(dat1)-70-digits.year
    
    yearStart <- as.numeric(regexpr(years[[1]],dat1))-1-dig.dif
    if(yearStart<0) yearStart <- 6
    #     if(yearStart>6) yearStart <- 6
    
    
    
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
    skip.lines <- 3
    ## Do nothing. read.fwf closes (and destroys ?!?) the file connection
    on.exit()
    ## Get chron stats if needed
    suppressWarnings(chron.stats <- utils::read.fwf(con, c(yearStart, digits.year, 6, 6, 6, 7, 9, 9, 10),
                                             skip=nlines-1, strip.white=TRUE, colClasses="character", stringsAsFactors=F))
    ## Unintuitively, the connection object seems to have been destroyed
    ## by the previous read.fwf.  We need to create a new one.
    con <- file(file, encoding = encoding)
    
    ## If columns 3 in chron.stats is an integer then there is no
    ## statistics line
    if((!is.integer(chron.stats[[3]]) | (chron.stats[[3]]==0)) & !grepl(" ",chron.stats[[3]])){
      names(chron.stats) <-
        c("SiteID", "nYears", "AC[1]", "StdDev", "MeanSens",
          "MeanRWI", "IndicesSum", "IndicesSS", "MaxSeries")

      ## Really read file
      dat <- utils::read.fwf(con, c(yearStart, digits.year, rep(c(4, 3), 10)),
                      skip=skip.lines, n=nlines-skip.lines-1,
                      strip.white=TRUE, colClasses="character", stringsAsFactors=F)
    } else {
      ## Really read file
      dat <- utils::read.fwf(con, c(yearStart, digits.year, rep(c(4, 3), 10)),
                      skip=skip.lines, n=nlines-skip.lines,
                      strip.white=TRUE, colClasses="character", stringsAsFactors=F)
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
    colnames(crn.mat) <- c("WIDTH", "DEPTH")
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
    crn.mat <- crn.mat[!apply(is.na(crn.mat[, -ncol.crn.mat, drop=FALSE]), 1, all),,drop=FALSE]
    
    seq.series <- seq_len(nseries)
    crn.mat[, "WIDTH"] <- crn.mat[, "WIDTH"] / 1000
    #     dt.widths <- data.table::data.table(t(crn.mat[,"WIDTH",drop=F]))
    #     dt.depths <- data.table::data.table(t(crn.mat[,"DEPTH",drop=F]))
    #     
    #     output <- list(WIDTH=dt.widths,DEPTH=dt.depths)
    
    return(data.frame(crn.mat))
    
  }else{
    raw.data <- scan(file,what='character', multi.line=F, fill=T, sep="\n", strip.white=T, quiet=T)
    tails <- toupper(substr_right(raw.data,3))
    raw.data <- split(raw.data,tails)
    
    raw.data <- lapply(raw.data,function(x){x[-c(1:3)]})
    
    meta <- scan(file,what='character', n=3,multi.line=F, fill=T, sep="\n", quiet=T)
    meta <- sub("  RAW","",meta)
    
    meta[2] <- gsub(substr(meta[2],1,9),"",meta[2])
    
    meta.2 <- unlist(strsplit(meta[2], "\\s+\\s+\\s+"))
    meta.2 <- c(meta.2[1],substr(meta.2[2],1,as.numeric(regexpr('\\d',meta.2[2]))-1),substr(meta.2[2],as.numeric(regexpr('\\d',meta.2[2])),nchar(meta.2[2])))
    meta.2 <- c(meta.2[1],meta.2[2],unlist(strsplit(meta.2[3], "\\s+")))
    if(length(which(meta.2=="-"))>0){
      meta.2 <- meta.2[-which(meta.2=="-")]
    }
    
    meta.2 <- meta.2[5:6]
    
    digits.year <- max(nchar(meta.2),4)
    
    widths <- c(6, digits.year, rep(c(4, 3), 10))
    starts <- c(1,cumsum(widths)+1)
    stops <- cumsum(widths)
    starts <- utils::head(starts,n=length(stops))
    dat <- lapply(raw.data,function(this.raw.data){
      dat <- as.data.frame(
        matrix(
          unlist(
            lapply(this.raw.data,FUN=function(x){
              sapply(1:length(starts), function(ii){substr(x, starts[ii], stops[ii])
              })
            })
          )
          ,ncol=22,byrow=T)
      )
      
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
      colnames(crn.mat) <- c("WIDTH", "DEPTH")
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
      crn.mat[, "WIDTH"] <- crn.mat[, "WIDTH"] / 1000
      #       dt.widths <- data.table::data.table(t(crn.mat[,"WIDTH",drop=F]))
      #       dt.depths <- data.table::data.table(t(crn.mat[,"DEPTH",drop=F]))
      #       
      #       output <- list(WIDTH=dt.widths,DEPTH=dt.depths)
      return(data.frame(crn.mat))
      
    })
    
    return(dat)
    
  }
}