#' Install and load a package.
#'
#'This is a convenience function that checks whether a package is installed, and if not, installs it.
#'
#' @param x A character string representing the name of a package.
pkgTest <- function(x){
  if(grepl("/",x)){
    pkgName <- basename(x)
  }else{
    pkgName <- x
  }
  if (!suppressWarnings(require(pkgName,character.only = TRUE)))
  {
    if(grepl("/",x)){
      suppressWarnings(devtools::install_github(x))
    }else{
      install.packages(x,dependencies=TRUE, repos="http://cran.rstudio.com")
    }
  }
  if(!suppressWarnings(require(pkgName,character.only = TRUE))) stop("Package not found")
}

#'Get the rightmost "n" characters of a character string.
#'
#' @param x A character string.
#' @param n The number of characters to retrieve.
#' @return A character string.
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#'Turn an extent object into a polygon
#'
#' @param x An \code{\link{extent}} object, or an object from which an extent object can be retrieved.
#' @param proj4string A PROJ.4 formatted string defining the required projection. If NULL, 
#' the function will attempt to get the projection from x using \code{\link{projection}}
#' @return A SpatialPolygons object.
polygonFromExtent <- function(x, proj4string=NULL){
  if(is.null(proj4string)){
    proj4string <- raster::projection(x)
  }
  
  if(class(x)!="extent"){
    x <- raster::extent(x)
  }
  
  extent.matrix <- rbind( c(x@xmin,x@ymin), c(x@xmin,x@ymax), c(x@xmax,x@ymax), c(x@xmax,x@ymin), c(x@xmin,x@ymin) ) # clockwise, 5 points to close it
  extent.SP <- sp::SpatialPolygons( list(sp::Polygons(list(sp::Polygon(extent.matrix)),"extent")), proj4string=sp::CRS(proj4string) )
  return(extent.SP)
}

#'Turn an SpatialPolygons object into a SpatialPolygonsDataFrame.
#'
#' @param x An SpatialPolygons object.
#' @return A SpatialPolygonsDataFrame object.
SPDFfromPolygon <- function(x){
  IDs <- sapply(slot(x, "polygons"), function(x) slot(x, "ID"))
  df <- data.frame(rep(0, length(IDs)), row.names=IDs)
  x <- sp::SpatialPolygonsDataFrame(x,df)
  return(x)
}

#'Get a logical vector of which elements in a vector are sequentially duplicated.
#'
#' @param x An vector of any type, or, if \code{rows}, a matrix.
#' @param rows Is x a matrix?
#' @return A logical vector of the same length as x.
sequential.duplicated <- function(x, rows=F){
  if(!rows){
    duplicates <- c(FALSE,unlist(lapply(1:(length(x)-1), function(i){duplicated(x[i:(i+1)])[2]})))
  }else{
    duplicates <- c(FALSE,unlist(lapply(1:(nrow(x)-1), function(i){duplicated(x[i:(i+1),])[2]})))
  }
  return(duplicates)
}

#'Use the \code{wget} command line tool to download a file.
#'
#' If both \code{timestamping} and \code{nc} are TRUE, timestamping behavior trumps nc.
#'
#' @param url The location of a file.
#' @param destdir Where the file should be downloaded to.
#' @param timestamping Should only newer files be downloaded?
#' @param nc Should files of the same type not be clobbered?
#' @return A logical vector of the same length as x.
wgetDownload <- function(url, destdir=getwd(), timestamping=T, nc=F){
  destdir <- gsub(" ","\\ ",destdir, fixed=T)
  if(!any(c(timestamping,nc))){
    status <- system(paste("wget -nd --quiet --directory-prefix=",destdir," ",url,sep=''))
  }else if(timestamping){
    status <- system(paste("wget -N -nd --quiet --directory-prefix=",destdir," ",url,sep=''))
  }else{
    status <- system(paste("wget -nc -nd --quiet --directory-prefix=",destdir," ",url,sep=''))
  }
  
  # If status is still not zero, report a warning
  if (status!=0)
    warning("Download of ",url," had nonzero exit status")
}

#' Use RCurl to download a file.
#'
#' This function makes it easy to implement timestamping and no-clobber of files.
#' Unlike \link{wgetDownload}, it doesn't require an external command-line tool to 
#' be installed.
#'
#' If both \code{timestamping} and \code{nc} are TRUE, nc behavior trumps timestamping.
#'
#' @param url The location of a file.
#' @param destdir Where the file should be downloaded to.
#' @param timestamping Should only newer files be downloaded?
#' @param nc Should files of the same type not be clobbered?
#' @param verbose Should cURL output be shown?
#' @param progress Should a progress bar be shown with cURL output?
#' @return A logical vector of the same length as x.
curlDownload <- function(url, destdir=getwd(), timestamping=T, nc=F, verbose=F, progress=F){
  
  destdir <- normalizePath(destdir)
  destfile <- paste0(destdir,'/',basename(url))
  #   destfile <- gsub(" ","\\ ",destfile, fixed=T)
  
  if(nc & file.exists(destfile)) return()
  
  if(timestamping & file.exists(destfile)){
    cat("\nDownloading file (if necessary):",url,"\n")
    temp.file <- paste0(tempdir(),"/",basename(url))
    f <- CFILE(temp.file, "wb")
    status <- curlPerform(url = url, 
                          writedata = f@ref,
                          verbose=verbose,
                          noprogress=!progress,
                          fresh.connect=T, 
                          ftp.use.epsv=F, 
                          forbid.reuse=T, 
                          timecondition=T, 
                          timevalue=base::file.info(destfile)$mtime)
    close(f)
    if(file.info(temp.file)$size > 0){
      file.copy(temp.file,destfile, overwrite=T)
    }
    
    #     status <- system(paste0("curl -R --globoff --create-dirs -z ",destfile," --url ",url," --output ",destfile))
  }else{
    cat("\nDownloading file:",url,"\n")
    temp.file <- paste0(tempdir(),"/",basename(url))
    f <- CFILE(temp.file, "wb")
    status <- curlPerform(url = url, 
                          writedata = f@ref,
                          verbose=verbose,
                          noprogress=!progress,
                          fresh.connect=T, 
                          ftp.use.epsv=F, 
                          forbid.reuse=T)
    close(f)
    file.copy(temp.file,destfile, overwrite=T)
    
    
    #     status <- system(paste0("curl -Rs --globoff --create-dirs --url ",url," --output ",destfile))
  }
  
  # If status is still not zero, report a warning
  if (!(status %in% c(0,3)))
    warning("Download of ",url," had nonzero exit status")
}