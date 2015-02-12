#' Install and load a package.
#'
#'This is a convenience function that checks whether a package is installed, and if not, installs it.
#'
#' @param x A character string representing the name of a package.
pkgTest <- function(x){
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dependencies=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
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