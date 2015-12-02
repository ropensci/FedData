#' Install and load a package.
#'
#'This is a convenience function that checks whether a package is installed, and if not, installs it.
#'
#' @param x A character string representing the name of a package.
#' @import data.table
#' @import sp
#' @export
pkg_test <- function(x){
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
      utils::install.packages(x,dependencies=TRUE, repos="http://cran.rstudio.com")
    }
  }
  if(!suppressWarnings(require(pkgName,character.only = TRUE))) stop("Package not found")
}

#'Get the rightmost "n" characters of a character string.
#'
#' @param x A character string.
#' @param n The number of characters to retrieve.
#' @return A character string.
#' @export
substr_right <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#'Turn an extent object into a polygon
#'
#' @param x An \code{\link{extent}} object, or an object from which an extent object can be retrieved.
#' @param proj4string A PROJ.4 formatted string defining the required projection. If NULL, 
#' the function will attempt to get the projection from x using \code{\link{projection}}
#' @return A SpatialPolygons object.
#' @export
polygon_from_extent <- function(x, proj4string=NULL){
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
#' @export
spdf_from_polygon <- function(x){
  IDs <- sapply((methods::slot(x, "polygons")), function(x){methods::slot(x, "ID")})
  df <- data.frame(rep(0, length(IDs)), row.names=IDs)
  x <- sp::SpatialPolygonsDataFrame(x,df)
  return(x)
}

#'Get a logical vector of which elements in a vector are sequentially duplicated.
#'
#' @param x An vector of any type, or, if \code{rows}, a matrix.
#' @param rows Is x a matrix?
#' @return A logical vector of the same length as x.
#' @export
sequential_duplicated <- function(x, rows=F){
  if(!rows){
    duplicates <- c(FALSE,unlist(lapply(1:(length(x)-1), function(i){duplicated(x[i:(i+1)])[2]})))
  }else{
    duplicates <- c(FALSE,unlist(lapply(1:(nrow(x)-1), function(i){duplicated(x[i:(i+1),])[2]})))
  }
  return(duplicates)
}

#'Unwraps a matrix and only keep the first n elements.
#'
#'A function that unwraps a matrix and only keeps the first n elements
#' n can be either a constant (in which case it will be repeated), or a vector
#' @param mat A matrix
#' @param n A numeric vector
#' @return A logical vector of the same length as x
#' @export
unwrap_rows <- function(mat,n){
  n <- rep_len(n,nrow(mat))
  i <- 0
  out <- lapply(1:nrow(mat),function(i){
    return(mat[i,1:n[i]])
  })
  return(as.numeric(do.call(c,out)))
}


#' Use RCurl to download a file.
#'
#' This function makes it easy to implement timestamping and no-clobber of files.
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
#' @export
curl_download <- function(url, destdir=getwd(), timestamping=T, nc=F, verbose=F, progress=F){
  
  destdir <- normalizePath(destdir)
  destfile <- paste0(destdir,'/',basename(url))
  
  if(nc & file.exists(destfile)) return()
  
  if(timestamping & file.exists(destfile)){
    message("Downloading file (if necessary): ",url)
    temp.file <- paste0(tempdir(),"/",basename(url))
    f <- RCurl::CFILE(temp.file, "wb")
    status <- RCurl::curlPerform(url = url, 
                          writedata = f@ref,
                          verbose=verbose,
                          noprogress=!progress,
                          fresh.connect=T, 
                          ftp.use.epsv=F, 
                          forbid.reuse=T, 
                          timecondition=T, 
                          timevalue=base::file.info(destfile)$mtime)
    RCurl::close(f)
    if(file.info(temp.file)$size > 0){
      file.copy(temp.file,destfile, overwrite=T)
    }
    
    #     status <- system(paste0("curl -R --globoff --create-dirs -z ",destfile," --url ",url," --output ",destfile))
  }else{
    message("Downloading file: ",url)
    temp.file <- paste0(tempdir(),"/",basename(url))
    f <- RCurl::CFILE(temp.file, "wb")
    status <- RCurl::curlPerform(url = url, 
                          writedata = f@ref,
                          verbose=verbose,
                          noprogress=!progress,
                          fresh.connect=T, 
                          ftp.use.epsv=F, 
                          forbid.reuse=T)
    RCurl::close(f)
    file.copy(temp.file,destfile, overwrite=T)
    
    
    #     status <- system(paste0("curl -Rs --globoff --create-dirs --url ",url," --output ",destfile))
  }
  
  # If status is still not zero, report a warning
  if (!(status %in% c(0,3)))
    warning("Download of ",url," had nonzero exit status")
}