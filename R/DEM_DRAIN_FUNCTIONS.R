# The master function of the DEM draining.
# Given a digital elevation model (DEM), return
# a version with estimates for elevations under 
# reservoirs and dams.
# Save downloaded data along the way.
drainDEM <- function(template, dem, label, NED.raw.dir, NHD.raw.dir, NRCS.raw.dir){

  res.elevs.rast <- dem

  NRCS <- extractNRCS(template=template, label=label, raw.dir=NRCS.raw.dir, force.redo=F)
  NHD <- extractNHD(template=template, label=label, raw.dir=NHD.raw.dir, remove.modern=TRUE, force.redo=F)

  mapunits <- read.csv(paste(NRCS.raw.dir,"/EXTRACTIONS/",label,"/tables/mapunit.csv",sep=''))
  dam.mukeys <- mapunits[mapunits$muname=="Dam",]$mukey
  if(!any(NRCS$MUKEY %in% dam.mukeys)){
    return(res.elevs.rast)
  }
  dams <- NRCS[NRCS$MUKEY%in%dam.mukeys,]
  dams <- spTransform(dams,CRS(projection(res.elevs.rast)))
  
  writeOGR(dams, paste(NRCS.raw.dir,"/EXTRACTIONS/",label,"/vectors",sep=''), "Dams", "ESRI Shapefile", overwrite_layer=TRUE)
  
  Areas <- readOGR(paste(NHD.raw.dir,"/EXTRACTIONS/",label,"/vectors",sep=''),"Area")
  Streams <- readOGR(paste(NHD.raw.dir,"/EXTRACTIONS/",label,"/vectors",sep=''),"Streams")
  Reservoirs <- readOGR(paste(NHD.raw.dir,"/EXTRACTIONS/",label,"/vectors",sep=''),"Reservoirs")
  Dams <- readOGR(paste(NRCS.raw.dir,"/EXTRACTIONS/",label,"/vectors",sep=''),"Dams")
  
  Dams <- suppressWarnings(gBuffer(Dams,width=1*res(dem)[1]))
  Dams <- spTransform(Dams,CRS(projection(Reservoirs)))
  
  reservoirsWithDams <- suppressWarnings(gWithinDistance(Reservoirs,Dams,dist=10*res(dem)[1],byid=T))
  
  reservoirsWithDams.vector <- vector("logical",length(Reservoirs))
  for(i in 1:length(Reservoirs)){
    reservoirsWithDams.vector[i] <- any(reservoirsWithDams[,i])
  }
  
  Reservoirs <- Reservoirs[reservoirsWithDams.vector]
  
  damsWithReservoirs <- suppressWarnings(gWithinDistance(Dams,Reservoirs,dist=10*res(dem)[1],byid=T))
  
  damsWithReservoirs.vector <- vector("logical",length(Dams))
  for(i in 1:length(Dams)){
    damsWithReservoirs.vector[i] <- any(damsWithReservoirs[,i])
  }
  
  Dams <- Dams[damsWithReservoirs.vector,]
  
  Streams.named <- Streams[!is.na(Streams$GNIS_Name),]
  Streams.unnamed <- Streams[is.na(Streams$GNIS_Name),]
  Streams.named.merged <- gLineMerge(Streams.named,byid=T,id=Streams.named$GNIS_Name)
  Streams.unnamed.merged <- gLineMerge(Streams.unnamed)
  Streams <- rbind(Streams.named.merged,Streams.unnamed.merged)
  Streams <- mergeStreams(Streams, projection(dem))
  
  # Join the dam and water polygons, and rasterize
  Dams <- spTransform(Dams,CRS(projection(Reservoirs)))
  bad.data.vect <- gUnion(Dams,Reservoirs)
  bad.data.rast <- raster::extract(dem, bad.data.vect, cellnumbers=T)

  res.elevs.rast[bad.data.rast[[1]][,1]] <- NA
  
  
  # Find which streams have missing elevation data
  bad.data.vect <- spTransform(bad.data.vect,CRS(projection(dem)))
  Streams.gapped <- suppressWarnings(Streams[c(Streams %over% gBuffer(bad.data.vect, width=10*res(dem)[1]))==1])
  
  Streams.gapped <- forceMergeStreams(Streams.gapped,CRS(projection(dem)))
    
  res.elevs.rast <- bootstrapDrainDEM(gappedDEM=res.elevs.rast, streams=Streams.gapped, reservoirs=Reservoirs, orig.dem=dem)
  
  dem.final <- fillIDW(res.elevs.rast)
#   writeGDAL(as(dem.final, "SpatialGridDataFrame"), paste(NED.raw.dir,"/EXTRACTIONS/",label,"/DEM_",res,"_DRAINED.tif", sep=''), drivername="GTiff", type="Float32", options=c("INTERLEAVE=PIXEL", "COMPRESS=DEFLATE", "ZLEVEL=9"))
  
  # A 5x5 mean smooth
  dem.final <- focal(dem.final, w=matrix(1,nrow=5,ncol=5), fun=mean, na.rm=T,pad=T)
  
  return(dem.final)
}



# A function to get the nearest minimum value on a raster, given a cell location.
getNearestMinimum <- function(raster,cell){
  done <- FALSE
  radius <- 1
#   counter <- 1
  while(!done){
#     counter <- counter+1
    nRowCols <- (radius*2)+1
    neighbors <- matrix(nrow=nRowCols,ncol=nRowCols,data=1)
    neighbors[radius+1,radius+1] <- 0
    minimum <- suppressWarnings(min(raster[adjacent(raster,cell,directions=neighbors)[,2]], na.rm=T))
    radius <- radius+1
    if(is.finite(minimum)) done <- TRUE
#     if(counter>=100) done <- TRUE
  }
  return(minimum)
}

# A function that fills missing values in a raster using Inverse Distance Weighted interpolation
# with an IDW exponent of 0.5.
fillIDW <- function(gappedDEM){
  # Finally, interpolate the missing values using IDW interpolation!
  # This uses the known elevations to derive an Inverse Distance Weighted model
  # of spatial prediction, then predicts elevation values at the still-missing points.
  gappedDEM.df <- data.frame(coordinates(gappedDEM),values(gappedDEM))
  names(gappedDEM.df) <- c("x","y","ELEVATION")
  gappedDEM.df.known <- gappedDEM.df[!is.na(gappedDEM.df$ELEVATION),]
  gappedDEM.df.unknown <- gappedDEM.df[is.na(gappedDEM.df$ELEVATION),]
  gappedDEM.idw.model <- gstat(id = "ELEVATION", formula = ELEVATION~1, locations = ~x+y, data=gappedDEM.df.known, nmax=7, set=list(idp = 0.5))
  gappedDEM.df.new <- predict.gstat(gappedDEM.idw.model,newdata=gappedDEM.df.unknown)
  gappedDEM.df.new$CELL <- cellFromXY(gappedDEM,as.matrix(gappedDEM.df.new[,1:2]))
  gappedDEM.final <- gappedDEM
  gappedDEM.final[gappedDEM.df.new$CELL] <- gappedDEM.df.new$ELEVATION.pred
  return(gappedDEM.final)
}



# A function that loads the Natural Resources Conservation Service dam polygons
getNRCS_Dams <- function(x, aggFactor=1, raw.dir="~/IMPORTANT/DATA/NRCS", layers=NULL, out.dir=paste(output.vectors,"../",sep=''), dsn.vectors=output.vectors, force.redo=FALSE){  
  dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(dsn.vectors, showWarnings = FALSE, recursive = TRUE)
  suppressWarnings(writeOGR(SPDFfromPolygon(polygonFromExtent(extent(x),projection(x))), dsn.vectors, "sim_poly","ESRI Shapefile", overwrite_layer=TRUE))
  
  
  # Load the NRCS study areas
  cat("\nLoading map of NRCS survey areas\n")
  NRCS.areas <- loadNRCSStudyAreas(x, raw.dir=raw.dir, dsn.vectors=dsn.vectors, force.redo=force.redo)
  
  # Load the NRCS mapunit polygons
  cat("\nLoading NRCS mapunit polygons\n")
  NRCS.polys <- loadNRCSMapUnitPolygons(x=x, raw.dir=raw.dir, dsn.vectors=dsn.vectors, force.redo=force.redo)
  
  # Load the primary mapunit data
  if(!file.exists(paste(out.dir,"/mapunit.csv",sep='')) | force.redo){
    NRCS.mapunit <- loadAndAggregateSoilTable("mapunit",NRCS.areas, raw.dir=raw.dir)[,c(1:5,12:24)]
    names(NRCS.mapunit) <- c("musym","muname","mukind","mustatus","muacres","farmlndcl","muhelcl","muwathelcl","muwndhelcl",'interpfocus',"invesintens","iacornsr","nhiforsoigrp","nhspiagr","vtsepticsyscl","mucertstat","lkey","mukey")
    NRCS.mapunit <- NRCS.mapunit[NRCS.mapunit$mukey %in% unique(NRCS.polys$MUKEY),]
    write.csv(NRCS.mapunit,paste(out.dir,"/mapunit.csv",sep=''),row.names=F)
  }
  
  # Create a final NRCS soils polygon file
  NRCS.final <- NRCS.polys[,c(4)]
  
  soils <- suppressWarnings(readOGR(normalizePath(paste(dsn.vectors,sep='')), layer='soils'))
  mapunits <- read.csv(paste(out.dir,"/mapunit.csv",sep=''))
  dam.mukeys <- mapunits[mapunits$muname=="Dam",]$mukey
  dams <- soils[soils$MUKEY%in%dam.mukeys,]
  dams <- spTransform(dams,CRS(projection(x)))
  
  writeOGR(dams, dsn.vectors, "Dams","ESRI Shapefile", overwrite_layer=TRUE)
}

# This loop is the work-horse of the algorithm. It first calculates the
# elevations of cells along each stream. If a stream BOTH enters and exits
# a reservoir, then it linearly interpolates between the entry and exit 
# elevations for all elevations in between. It then assigns those elevations
# to the "holy" DEM along the stream, and also at a 1-cell buffer (this is dependent
# on the resolution of one's DEM).
#
# The algorithm then loops back and extracts the elevations along streams again.
# This time, streams that end at one of the streams that get interpolated in the 
# first loop will have "ending" elevations, and may be interpolated themselves, and
# interpolated values added to the DEM.
#
# This process is repeated until no more streams may be interpolated.
bootstrapDrainDEM <- function(gappedDEM,streams,reservoirs,orig.dem){
  done <- FALSE
  last <- FALSE
  extracts <- vector("list", length(streams))
  finished.lines <- vector("logical",length(extracts))
  counter <- 1
  while(!done | !last){
    
    cat("Finished lines: ",finished.lines,"\n")
    flush.console()
    if(done){
      last <- TRUE
    }
    cat("\nBootstrap-draining DEM, iteration",counter,"\n")
    flush.console()
    counter <- counter+1
    # A variable used to assess our progress
    start.length <- length(finished.lines[finished.lines==TRUE])
    
    # To make things faster, first calculate the elevations at the line-ends.
    # Only compute lines that have both ends.
    ends <- vector("list",length(streams))
    end.logic <- vector("logical",length(streams))
    for(i in 1:length(streams)){
      if(finished.lines[i]) next
      cat("Calculating elevation at line ends: ",i,"\n")
      flush.console()
      
      end.points <- rbind(head(streams[i]@lines[[1]]@Lines[[1]]@coords, n=1),tail(streams[i]@lines[[1]]@Lines[[1]]@coords, n=1))
      
      if(end.points[1,1]>xmax(gappedDEM)) end.points[1,1] <- xmax(gappedDEM)
      if(end.points[2,1]>xmax(gappedDEM)) end.points[2,1] <- xmax(gappedDEM)
      if(end.points[1,1]<xmin(gappedDEM)) end.points[1,1] <- xmin(gappedDEM)
      if(end.points[2,1]<xmin(gappedDEM)) end.points[2,1] <- xmin(gappedDEM)
      if(end.points[1,2]>ymax(gappedDEM)) end.points[1,2] <- ymax(gappedDEM)
      if(end.points[2,2]>ymax(gappedDEM)) end.points[2,2] <- ymax(gappedDEM)
      if(end.points[1,2]<ymin(gappedDEM)) end.points[1,2] <- ymin(gappedDEM)
      if(end.points[2,2]<ymin(gappedDEM)) end.points[2,2] <- ymin(gappedDEM)
      
      ends[[i]] <- raster::extract(gappedDEM,end.points,cellnumbers=T)
      
      if(is.na(ends[[i]][1,2])){
        minLocalValue <- suppressWarnings(min(gappedDEM[adjacent(gappedDEM,ends[[i]][1,1],directions=8,include=T)[,2]],na.rm=T))
        if(is.finite(minLocalValue)){
          gappedDEM[ends[[i]][1,1]] <- minLocalValue
          ends[[i]][1,2] <- minLocalValue
        }
      }
      
      if(is.na(ends[[i]][2,2])){
        minLocalValue <- suppressWarnings(min(gappedDEM[adjacent(gappedDEM,ends[[i]][2,1],directions=8,include=T)[,2]],na.rm=T))
        if(is.finite(minLocalValue)){
          gappedDEM[ends[[i]][2,1]] <- minLocalValue
          ends[[i]][2,2] <- minLocalValue
        }
      }
      
      end.logic[i] <- all(!is.na(ends[[i]][,2]))
      
    }
    
    # Calculate the elevations of the cells along each stream.
    # Don't re-calculate if a stream has already been interpolated.
    for(i in 1:length(streams)){
      if(finished.lines[i] | !end.logic[i]) next
      cat("Calculating elevations along line: ",i,"\n")
      flush.console()
      
      extracts[[i]] <- raster::extract(gappedDEM,streams[i], along=T, cellnumbers=T)[[1]]
    }
    
    if(last){
      ends <- vector("list",length(streams))
      for(i in 1:length(streams)){
        if(finished.lines[i]) next
        
        ends[[i]] <- raster::extract(gappedDEM,matrix(c(head(streams[i]@lines[[1]]@Lines[[1]]@coords, n=1),tail(streams[i]@lines[[1]]@Lines[[1]]@coords, n=1)),nrow=2, byrow=T),cellnumbers=TRUE)
        
        if(length((ends[[i]][,1])[is.na(ends[[i]][,2])])>0){
          minimum <- getNearestMinimum(gappedDEM,(ends[[i]][,1])[is.na(ends[[i]][,2])])
          gappedDEM[(ends[[i]][,1])[is.na(ends[[i]][,2])]] <- minimum
        }
        extracts[[i]] <- raster::extract(gappedDEM,streams[i], along=T, cellnumbers=T)[[1]]
      }
    }
    
    # For each extracted stream...
    for(i in 1:length(extracts)){
      # Don't re-process if already finished
      if(finished.lines[i] | is.null(extracts[[i]])){
        next
      }
      
      # Don't process if there are less than 4 total cells covered, or
      # if there are no NA cells. Consider these streams "finished"
      
      if(length(extracts[[i]][,2])<4 | all(!is.na(extracts[[i]][,2]))){
        finished.lines[i] <- TRUE
        next
      }
      
      cat("Processing line: ",i,"\n")
      flush.console()
      temp <- data.frame(1:length(extracts[[i]][,2]),extracts[[i]][,2],extracts[[i]][,1])
      names(temp) <- c("index","elevation","cell")
      if(all(is.na(head(temp$elevation,n=6))) | all(is.na(tail(temp$elevation,n=6)))) next
      
      temp$elevation <- make.monotonic(temp$elevation)
      
      gaps <- vector("list")
      temp$group <- 0
      group.counter <- 0
      for(j in 2:nrow(temp)){
        # Enters a batch of NAs
        if(!is.na(temp$elevation[j]) & is.na(temp$elevation[j+1])){
          if(!is.na(temp$elevation[j-1]) & j!=nrow(temp)){
            group.counter <- group.counter+1
          }
          temp$group[j] <- group.counter
        }else if(!is.na(temp$elevation[j]) & is.na(temp$elevation[j-1])){
          temp$group[j] <- group.counter
          if(!is.na(temp$elevation[j+1]) & !is.na(temp$elevation[j+2])){
            group.counter <- group.counter+1
          }
        }else{
          temp$group[j] <- group.counter
        }
      }
      
      temp.na <- temp[is.na(temp$elevation),]
      
      for(j in unique(temp$group)){
        if(any(is.na(temp[temp$group==j,]$elevation)) & length(temp[temp$group==j,]$elevation[!is.na(temp[temp$group==j,]$elevation)])>=2){
          temp[temp$group==j,]$elevation <- approx(temp[temp$group==j,]$index,temp[temp$group==j,]$elevation, xout=temp[temp$group==j,]$index)$y
        }
      }
      
      # Update the raster
      gappedDEM[temp$cell] <- temp$elevation
      
      local <- adjacent(gappedDEM,temp.na$cell, pairs=F, directions=16)
      
      gappedDEM[local] <-  gappedDEM[cellFromXY(gappedDEM,t(apply(xyFromCell(gappedDEM,local),1,function(x){ nearestPointOnLine(streams[i]@lines[[1]]@Lines[[1]]@coords, x) })))]
      
      gappedDEM <- setMinMax(gappedDEM)
      
      finished.lines[i] <- TRUE
    }
    if(length(finished.lines[finished.lines==TRUE])==start.length) done <- TRUE
  }
  
  if(any(!finished.lines)){
    # Lakes where no interpolated values could be determined (usually 
    # shallow stock-ponds) are re-filled with their original elevations.
    reservoirs <- spTransform(reservoirs,CRS(projection(gappedDEM)))
    unfinished.streams <- spTransform(Streams.gapped[!finished.lines],CRS(projection(gappedDEM)))
    if(any(is.na(c(reservoirs %over% unfinished.streams)))){
      reservoirs.no.fill <- reservoirs[is.na(c(reservoirs %over% unfinished.streams)),]
      reservoirs.no.fill.rast <- raster::extract(orig.dem,reservoirs.no.fill, cellnumbers=T)
      for(i in 1:length(reservoirs.no.fill.rast)){
        gappedDEM[reservoirs.no.fill.rast[[i]][,1]] <- reservoirs.no.fill.rast[[i]][,2]
      }
    }
  }
  
  return(gappedDEM)
}

# A function that takes a sometimes erratic stream profile
# and makes it monotonic (strictly downward or upward sloping) 
# for the purposes of interpolating across gaps. 
# The monotonic values do not replace the existing values in the DEM.
make.monotonic <- function(x){
  for(i in 2:length(x)){
    if(!is.na(x[i]) & !is.na(x[i-1])){
      if(!is.na(x[i-1]) & x[i]>x[i-1]){
        x[i] <- x[i-1]
      }
    }
  }
  return(x)
}