pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE, type="source")
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

# A method to crop vector shapefiles to a given study area.
crop.to.studyArea <- function(x,y) {
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

#plot.stacked makes a stacked plot where each y series is plotted on top
#of the each other using filled polygons
#
#Arguments include:
#'x' - a vector of values
#'y' - a matrix of data series (columns) corresponding to x
#'order.method' = c("as.is", "max", "first") 
#  "as.is" - plot in order of y column
#  "max" - plot in order of when each y series reaches maximum value
#  "first" - plot in order of when each y series first value > 0
#'col' - fill colors for polygons corresponding to y columns (will recycle)
#'border' - border colors for polygons corresponding to y columns (will recycle) (see ?polygon for details)
#'lwd' - border line width for polygons corresponding to y columns (will recycle)
#'...' - other plot arguments

plot.stacked <- function(
  x, y, 
  order.method = "as.is",
  ylab="", xlab="", 
  border = NULL, lwd=1, 
  col=rainbow(length(y[1,])),
  ylim=NULL,
  ...
){
  
  #   if(sum(y < 0) > 0) error("y cannot contain negative numbers")
  
  if(is.null(border)) border <- par("fg")
  border <- as.vector(matrix(border, nrow=ncol(y), ncol=1))
  col <- as.vector(matrix(col, nrow=ncol(y), ncol=1))
  lwd <- as.vector(matrix(lwd, nrow=ncol(y), ncol=1))
  
  if(order.method == "max") {
    ord <- order(apply(y, 2, which.max))
    y <- y[, ord]
    col <- col[ord]
    border <- border[ord]
  }
  
  if(order.method == "first") {
    ord <- order(apply(y, 2, function(x) min(which(r>0))))
    y <- y[, ord]
    col <- col[ord]
    border <- border[ord]
  }
  
  top.old <- x*0
  polys <- vector(mode="list", ncol(y))
  for(i in seq(polys)){
    top.new <- top.old + y[,i]
    polys[[i]] <- list(x=c(x, rev(x)), y=c(top.old, rev(top.new)))
    top.old <- top.new
  }
  
  if(is.null(ylim)) ylim <- range(sapply(polys, function(x) range(x$y, na.rm=TRUE)), na.rm=TRUE)
  plot(x,y[,1], ylab=ylab, xlab=xlab, ylim=ylim, t="n", ...)
  for(i in seq(polys)){
    polygon(polys[[i]], border=border[i], col=col[i], lwd=lwd[i])
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

scalebar.new <- function (d, xy = NULL, height = NULL, line.offset=c(0,0), side="right", lab.side='top', lonlat = NULL, label, adj = c(0.5, -0.5), lwd = 2, ...){
  pr <- par()
  if (is.null(lonlat)) {
    if (pr$usr[1] > -181 & pr$usr[2] < 181 & pr$yaxp[1] > 
          -200 & pr$yaxp[2] < 200) {
      lonlat <- TRUE
    }
    else {
      lonlat <- FALSE
    }
  }
  
  if (lonlat) {
    lat <- mean(pr$yaxp[1:2])
    if (missing(d)) {
      dx <- (pr$usr[2] - pr$usr[1])/10
      d <- pointDistance(cbind(0, lat), cbind(dx, lat), 
                         TRUE)
      d <- signif(d/1000, 2)
      label <- NULL
    }
    p <- cbind(0, lat)
    dd <- raster:::.destPoint(p, d * 1000)
    dd <- dd[1, 1]
  } else {
    if (missing(d)) {
      d <- round(10 * (pr$usr[2] - pr$usr[1])/10)/10
      label <- NULL
    }
    dd <- d
  }
  
  if (is.null(xy)) {
    padding = c(5, 5)/100
    parrange <- c(pr$usr[2] - pr$usr[1], pr$usr[4] - pr$usr[3])
    xy <- c(pr$usr[1] + (padding[1] * parrange[1]), pr$usr[3] + 
              (padding[2] * parrange[2]))
  }
  
  xy <- xy + line.offset
  
  if(side=='right'){
    xstart = xy[1]
    xend = xy[1] + dd
  }else{
    xstart = xy[1] - dd
    xend = xy[1]
  }
  
  if(is.null(height)){
    height <- dd * 0.1
  }
  
  rect(xleft=xstart, ybottom=xy[2], xright=xend, ytop=xy[2]+height, col='black',border=NA, lend=1)
  
#   lines(matrix(c(xstart, xy[2], xend, xy[2]), byrow = T, nrow = 2), lend=1, lwd = lwd, ...)
  
  if (missing(label)) {
    label <- paste(d)
  }
  if (is.null(label)) {
    label <- paste(d)
  }
  if (missing(adj)) {
    adj <- c(0.5, -0.2 - lwd/20)
  }
  
  if(lab.side=='top'){
    text(mean(c(xstart,xend)), xy[2], labels = label, adj = c(0.5,-0.5), 
         ...)
  }else if(lab.side=='right'){
    text(xend, xy[2], labels = label, adj = c(-0.1,0), 
         ...)
  }else if(lab.side=='left'){
    text(xstart, xy[2], labels = label, adj = c(1.1,0), 
         ...)
  }else if(lab.side=='bottom'){
    text(mean(c(xstart,xend)), xy[2], labels = label, adj = c(0.5,1.5), 
         ...)
  }
}