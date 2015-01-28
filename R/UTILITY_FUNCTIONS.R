pkgTest <- function(x){
  if (!require(x,character.only = TRUE))
  {
    update.packages()
    install.packages(x,dep=TRUE, type="both")
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

createArea <- function(North,South,East,West,projection.string){
  # Create matrix of coordinates
  datainUTM<-matrix(c(East, West, West, East,East, North,North,South,South,North),nrow=5)
  
  # Set universal projection
  master.proj <- CRS(projection.string)
  
  # Create SpatialPolygon of simulation area
  sim.poly <- Polygons(list(Polygon(datainUTM, hole=FALSE)),ID='A')
  sim.poly <- SpatialPolygons(list(sim.poly), proj4string=master.proj)
  IDs <- sapply(slot(sim.poly, "polygons"), function(x) slot(x, "ID"))
  df <- data.frame(rep(0, length(IDs)), row.names=IDs)
  sim.poly <- SpatialPolygonsDataFrame(sim.poly,df)
  
  return(sim.poly)
}

polygonFromExtent <- function(x, proj4string=NULL){
  if(is.null(proj4string)){
    proj4string <- projection(x)
  }
  
  if(class(x)!="extent"){
    x <- extent(x)
  }
  
  extent.matrix <- rbind( c(x@xmin,x@ymin), c(x@xmin,x@ymax), c(x@xmax,x@ymax), c(x@xmax,x@ymin), c(x@xmin,x@ymin) ) # clockwise, 5 points to close it
  extent.SP <- SpatialPolygons( list(Polygons(list(Polygon(extent.matrix)),"extent")), proj4string=CRS(proj4string) )
  return(extent.SP)
}

SPDFfromPolygon <- function(x){
  IDs <- sapply(slot(x, "polygons"), function(x) slot(x, "ID"))
  df <- data.frame(rep(0, length(IDs)), row.names=IDs)
  x <- SpatialPolygonsDataFrame(x,df)
  return(x)
}

# A method to crop vector shapefiles to a given polygon.
cropToPoly <- function(x,y) {
  gI <- gIntersects(x,y,byid=TRUE)
  out <- vector(mode="list",length=length(which(gI)))
  ii <- 1
  if(length(out)==0) return(NULL)
  for (i in seq(along=gI)){
    if (gI[i]) {
      out[[ii]] <- gIntersection(x[i,], y)
      row.names(out[[ii]]) <- row.names(x)[i]
      ii <- ii+1
    }
  }
  out_class <- sapply(out, class)
  
  if(class(x)[1]=="SpatialLinesDataFrame" | class(x)[1]=="SpatialLines"){
    ri <- do.call("rbind", out[out_class == "SpatialLines"])
    if(class(x)[1]=="SpatialLines"){
      return(ri)
    }
    ri <- SpatialLinesDataFrame(ri,as.data.frame(x))
    return(ri)
  }else if(class(x)[1]=="SpatialPolygonsDataFrame"){
    ri <- do.call("rbind", out[out_class == "SpatialPolygons"])
    ri <- SpatialPolygonsDataFrame(ri,as.data.frame(x)[row.names(ri),])
    return(ri)
  }else if(class(x)[1]=="SpatialPointsDataFrame"){
    ri <- do.call("rbind", out[out_class == "SpatialPoints"])
    ri <- SpatialPointsDataFrame(ri,as.data.frame(x)[row.names(ri),],match.ID = TRUE)
    return(ri)
  }
}

sequential.duplicated <- function(x, rows=F){
  if(!rows){
    duplicates <- c(FALSE,unlist(lapply(1:(length(x)-1), function(i){duplicated(x[i:(i+1)])[2]})))
  }else{
    duplicates <- c(FALSE,unlist(lapply(1:(nrow(x)-1), function(i){duplicated(x[i:(i+1),])[2]})))
  }
  return(duplicates)
}

wgetDownload <- function(url, destdir, timestamping=T){
  if(timestamping){
    status <- system(paste("wget -N -nd --quiet --directory-prefix=",destdir," ",url,sep=''))
  }else{
    status <- system(paste("wget -nd --quiet --directory-prefix=",destdir," ",url,sep=''))
  }
  if (status!=0) 
    warning("Download of ",url," had nonzero exit status")
}