## All data is downloaded from the NRCS,
## available at http://websoilsurvey.sc.egov.usda.gov

## Author: R. Kyle Bocinsky
## Date: 02/14/2014

extractNRCS <- function(template, label, raw.dir, extraction.dir=NULL, SFNF.dir=NULL, NED.raw.dir=NULL, NHD.raw.dir=NULL, fillReservoirs=T, force.redo=F, return.rast=F){  
  if(is.null(extraction.dir)){
    extraction.dir <- paste(raw.dir,"/EXTRACTIONS",sep='')
  }
  dir.create(paste(extraction.dir,"/",label,"/",sep=''), showWarnings = FALSE, recursive=T)
  dsn.vectors <- paste(extraction.dir,"/",label,"/vectors",sep='')
  dir.create(dsn.vectors, showWarnings = FALSE, recursive = TRUE)
  
  tables.dir <- paste(extraction.dir,"/",label,"/tables",sep='')
  dir.create(tables.dir, showWarnings = FALSE, recursive = TRUE)
  
  template.poly.latlon <- SPDFfromPolygon(spTransform(polygonFromExtent(template),CRS("+proj=longlat +ellps=GRS80")))
  
  # Write a polygon of the study area
  suppressWarnings(writeOGR(template.poly.latlon, dsn.vectors, "sim_poly","ESRI Shapefile", overwrite_layer=TRUE))
  
  # Load the NRCS study areas
  #   cat("\nLoading map of NRCS survey areas\n")
  NRCS.areas <- loadNRCSStudyAreas(x=template, raw.dir=raw.dir, dsn.vectors=dsn.vectors, force.redo=force.redo)
  
  template.poly.latlon <- SPDFfromPolygon(spTransform(polygonFromExtent(template),CRS(projection(NRCS.areas))))
  suppressWarnings(writeOGR(template.poly.latlon, dsn.vectors, "sim_poly","ESRI Shapefile", overwrite_layer=TRUE))
  
  #   plot(NRCS.areas)
  
  # Load the NRCS mapunit polygons
  #     cat("\nLoading NRCS mapunit polygons\n")
  if(file.exists(paste(extraction.dir,"/",label,"/vectors/soils.shp", sep='')) & !force.redo){
    NRCS.polys <- readOGR(dsn.vectors, "soils", verbose=F)
  }else{
    NRCS.polys <- loadNRCSMapUnitPolygons(x=template, raw.dir=raw.dir, dsn.vectors=dsn.vectors, force.redo=force.redo)
    suppressWarnings(writeOGR(NRCS.polys, dsn.vectors, "soils","ESRI Shapefile", overwrite_layer=TRUE))
  }
  
  #   plot(NRCS.polys)
  
  # Load tabular soil data for study area
  #   cat("\nLoad tabular soil data for study area\n")
  getSoilData(x=template, areas=NRCS.areas[NRCS.areas@data$iscomplete != 0,], polys=NRCS.polys, raw.dir=raw.dir, dsn.vectors=dsn.vectors, tables.dir=tables.dir, force.redo=force.redo)
  
  #   texture <- calculateTexture(raw.dir=raw.dir, dsn.vectors=dsn.vectors, tables.dir=tables.dir)
  
  #   if("NM678" %in% NRCS.areas$areasymbol){
  #     if(!is.na(SFNF.dir)){
  #       if(NRCS.areas[NRCS.areas$areasymbol=="NM678",]$iscomplete == 0){
  #         if(is.null(SFNF.dir)){
  #           SFNF.dir <- readline("Please provide a path for the SFNF data:")
  #         }
  #         SFNF <- loadSFNF_Data(x=template, areas=NRCS.areas, polys=NRCS.polys, SFNF.dir=SFNF.dir, raw.dir=raw.dir, dsn.vectors=dsn.vectors, force.redo=force.redo)
  #         texture <- rbind(texture,SFNF)
  #       }
  #     }
  #   }
  
  #   names(texture) <- c("MUKEY","LEVEL","MUKEY_P","SAND","SILT","CLAY","TOP","BOTTOM")
  
  #   write.csv(texture,paste(tables.dir,"/level_texture.csv",sep=''),row.names=F)
  #   
  #   
  #   colors <- colorRampPalette(brewer.pal(9,"Greens"))(1000)
  #   props <- unique(texture[,c("mukey","mukey_p")])
  #   props$colors <- colors[round(props$mukey_p*1000, digits=0)]
  #   names(props) <- c("MUKEY","P","COLOR")
  #   NRCS.polys <- NRCS.polys[NRCS.polys$MUKEY %in% texture$mukey,]
  #   NRCS.polys@data <- merge(NRCS.polys@data,props, all.x=T)
  #   plot(NRCS.polys, border=NA, col=NRCS.polys$COLOR)
  #   plot(template.poly.latlon, add=T)
  
  #   TT.points.in.classes(tri.data= NRCS.chorizon[!is.na(NRCS.chorizon[,9]),6:8],class.sys= "USDA.TT", PiC.type='t')
  
  
  #   
  #   # Create a final NRCS soils polygon file
  #   NRCS.final <- NRCS.polys[,c(4)]
  
  
  # Create VEPII IDs for each MUKEY
  #   cat("\nCreate VEPII IDs for each MUKEY\n")
  #   NRCS.polys <- loadNRCSMapUnitPolygons(x=template, raw.dir=raw.dir, dsn.vectors=dsn.vectors, force.redo=F)
  #   NRCS.polys.mukeys <- data.frame(MUKEY=unique(NRCS.polys$MUKEY),ID=1:length(unique(NRCS.polys$MUKEY)))
  #   NRCS.polys <- sp::merge(NRCS.polys,NRCS.polys.mukeys)
  
  # Export final vector dataset for the study area.
  #     cat("Exporting final vector dataset for the study area.\n")
  
  if(fillReservoirs){
    cat("Filling reservoirs.\n")
    if(file.exists(paste(raw.dir,"EXTRACTIONS/",area.name,"/RASTERIZED_MUKEYS_1arcsec_filled.tif",sep='')) & !force.redo){
      return(NRCS.polys)
    }
    
    if(is.null(NED.raw.dir)){
      NED.raw.dir <- readline("Please provide a path for the raw NHD data directory:")
    }
    
    if(is.null(NHD.raw.dir)){
      NHD.raw.dir <- readline("Please provide a path for the raw NHD data directory:")
    }
    
    NED <- extractNED(template=sim.poly, label=area.name, raw.dir=NED.raw.dir, res="1", drain=T, force.redo=F)
    
    if(!file.exists(paste(raw.dir,"EXTRACTIONS/",area.name,"/RASTERIZED_MUKEYS_1arcsec.tif",sep=''))){
      NRCS.rast <- raster::rasterize(NRCS.polys,NED,field="MUKEY", na.rm=T)
      writeGDAL(as(NRCS.rast, "SpatialGridDataFrame"),paste(raw.dir,"EXTRACTIONS/",area.name,"/RASTERIZED_MUKEYS_1arcsec.tif",sep=''), drivername="GTiff", type="Int16", mvFlag=-32768, options=c("INTERLEAVE=PIXEL", "COMPRESS=DEFLATE", "ZLEVEL=9"))
    }else{
      NRCS.rast <- raster(paste(raw.dir,"EXTRACTIONS/",area.name,"/RASTERIZED_MUKEYS_1arcsec.tif",sep=''))
    }
    
    if(file.exists(paste(NHD.raw.dir,"EXTRACTIONS/",area.name,"/vectors/Reservoirs.shp", sep=''))){
      reservoirs <- readOGR(paste(NHD.raw.dir,"EXTRACTIONS/",area.name,"/vectors", sep=''),"Reservoirs", verbose=F)
    }else{
      NHD <- extractNED(template=sim.poly, label=area.name, raw.dir=NHD.raw.dir, force.redo=F)
    }
    
    if(file.exists(paste(raw.dir,"EXTRACTIONS/",area.name,"/vectors/Dams.shp", sep=''))){
      dams <- readOGR(paste(raw.dir,"EXTRACTIONS/",area.name,"/vectors", sep=''),"Dams", verbose=F)
    }else{
      mapunits <- read.csv(paste(raw.dir,"EXTRACTIONS/",label,"/tables/mapunit.csv",sep=''))
      dam.mukeys <- mapunits[mapunits$muname=="Dam",]$mukey
      
      dams <- NRCS.polys[NRCS.polys$MUKEY%in%dam.mukeys,]
      
      writeOGR(dams, paste(raw.dir,"EXTRACTIONS/",area.name,"/vectors",sep=''), "Dams", "ESRI Shapefile", overwrite_layer=TRUE)
      
    }
    #       areas <- readOGR(paste(NHD.dir,"EXTRACTIONS/",area.name,"/vectors", sep=''),"NHDArea", verbose=F)
    #       areas <- areas[areas$AreaSqKm>0.16,]
    #       for(i in 1:length(areas)){
    #         areas@polygons[[i]] <- remove.holes(areas@polygons[[i]])
    #       }
    
    # Get raster cells in reservoirs and dams, and set them to NA
    projection(dams) <- projection(reservoirs)
    bad.data.vect <- gUnion(dams,reservoirs)
    bad.data.vect <- spTransform(bad.data.vect,CRS(projection(NRCS.rast)))
    #       bad.data.vect <- gUnion(bad.data.vect,areas)
    bad.data.rast <- raster::extract(NRCS.rast, bad.data.vect, cellnumbers=T)
    NRCS.rast[bad.data.rast[[1]][,1]] <- NA
    
    # Fill the missing reservoir/dam soils using linear discriminant analysis
    NRCS.rast.filled <- fillReservoirSoils(gapped.soil.raster=NRCS.rast, dem.raster=NED,  label=area.name, raw.dir=paste(MASTER.DATA,"NRCS/",sep=''), force.redo=F)
    projection(NRCS.rast.filled) <- projection(NRCS.rast)
    NRCS.rast <- NRCS.rast.filled
    
    writeGDAL(as(NRCS.rast, "SpatialGridDataFrame"),paste(raw.dir,"EXTRACTIONS/",area.name,"/RASTERIZED_MUKEYS_1arcsec_filled.tif",sep=''), drivername="GTiff", type="Int16", mvFlag=-32768, options=c("INTERLEAVE=PIXEL", "COMPRESS=DEFLATE", "ZLEVEL=9"))
  }
  
  if(return.rast){
    if(!exists(NRCS.rast)){
      if(!file.exists(paste(raw.dir,"EXTRACTIONS/",area.name,"/RASTERIZED_MUKEYS_1arcsec.tif",sep=''))){
        NRCS.rast <- raster::rasterize(NRCS.polys,NED,field="MUKEY", na.rm=T)
        writeGDAL(as(NRCS.rast, "SpatialGridDataFrame"),paste(raw.dir,"EXTRACTIONS/",area.name,"/RASTERIZED_MUKEYS_1arcsec.tif",sep=''), drivername="GTiff", type="Int16", mvFlag=-32768, options=c("INTERLEAVE=PIXEL", "COMPRESS=DEFLATE", "ZLEVEL=9"))
      }else{
        NRCS.rast <- raster(paste(raw.dir,"EXTRACTIONS/",area.name,"/RASTERIZED_MUKEYS_1arcsec.tif",sep=''))
      }
    }
    return(NRCS.rast)
  }else{
    return(NRCS.polys)
  }
  
}


loadNRCSStudyAreas <- function(x, raw.dir="../Input/NRCS", dsn.vectors="Output/NRCS_vectors.gdb",force.redo=F){
  # Get NRCS study areas, and save them to disk
  if(!("NRCS_SurveyAreas" %in% ogrListLayers(dsn.vectors)) | force.redo){
    NRCS.areas <- getNRCSStudyAreas(x, raw.dir=raw.dir)
    suppressWarnings(writeOGR(NRCS.areas,dsn.vectors,"NRCS_SurveyAreas","ESRI Shapefile", overwrite_layer=TRUE))
  }else{
    NRCS.areas <- readOGR(dsn.vectors,"NRCS_SurveyAreas")
  }
  
  # Check to see if all survey areas are available
  if(0 %in% NRCS.areas@data$iscomplete){
    cat("WARNING! Some of the soil surveys in your area are unavailable.\n")
    cat("Soils and productivity data will have holes.\n")
    cat("Missing areas:\n")
    cat(as.vector(NRCS.areas@data[NRCS.areas@data$iscomplete==0,]$areasymbol))
    cat("\n\n")
    cat("Continuing with processing available soils.\n\n")
    
    #     NRCS.areas <- NRCS.areas[NRCS.areas@data$iscomplete != 0,]
  }
  
  return(NRCS.areas)
}


getNRCSStudyAreas <- function(x, raw.dir="../Input/NRCS"){
  # Import the shapefile of NRCS study areas.
  # This is available at
  # http://soildatamart.sc.egov.usda.gov/download/StatusMaps/soilsa_a_nrcs.zip
  if(url.exists('http://websoilsurvey.sc.egov.usda.gov/DataAvailability/SoilDataAvailabilityShapefile.zip')){
    system(paste("wget -np -nd -N http://websoilsurvey.sc.egov.usda.gov/DataAvailability/SoilDataAvailabilityShapefile.zip --directory-prefix=",raw.dir, sep=''))
  }else{
    stop("Unable to download NRCS study area shapefile! Please find a copy of the NRCS study areas shapefile, available from the USGS.")
  }
  unzip(paste(raw.dir,'/SoilDataAvailabilityShapefile.zip', sep=''),exdir=paste(raw.dir,'/SoilDataAvailabilityShapefile',sep=''))
  NRCS.areas <- suppressWarnings(readOGR(paste(raw.dir,'/SoilDataAvailabilityShapefile',sep=''), layer="soilsa_a_nrcs", verbose=FALSE))
  
  unlink(paste(raw.dir,'/SoilDataAvailabilityShapefile',sep=''), recursive=T, force=T)
  
  # Get a list of NHD subregions within the project study area
  NRCS.areas <- crop.to.studyArea(NRCS.areas,spTransform(x,CRS(projection(NRCS.areas))))
  
  return(NRCS.areas)
}

loadNRCSMapUnitPolygons <- function(x, raw.dir="../Input/NRCS", dsn.vectors="Output/NRCS_vectors.gdb", force.redo=F){
  if(!("soils" %in% ogrListLayers(dsn.vectors)) | force.redo){
    NRCS.areas <- loadNRCSStudyAreas(x, raw.dir=raw.dir, dsn.vectors=dsn.vectors,force.redo=F)
    
    NRCS.areas <- NRCS.areas[NRCS.areas@data$iscomplete != 0,]
    
    NRCS.polys <- getNRCSMapUnitPolygons(x, NRCS.areas, raw.dir=raw.dir)
    
    # Export final vector dataset for the study area.
    #     cat("Exporting final vector dataset for the study area.\n")
    suppressWarnings(writeOGR(NRCS.polys, dsn.vectors, "soils","ESRI Shapefile", overwrite_layer=TRUE))
  }else{
    NRCS.polys <- readOGR(dsn.vectors,"soils")
  }
  return(NRCS.polys)
}


getNRCSMapUnitPolygons <- function(x, NRCS.areas, raw.dir="../Input/NRCS"){
  # Load raw NRCS Map Unit polygons from the regions specified.
  cat("Loading the NRCS soil survey data for each region.\n")
  NRCS.polys <- vector("list", length(NRCS.areas))
  for(i in 1:length(NRCS.areas)){
    cat("\n(Down)Loading soils data for",as.character(NRCS.areas@data$areasymbol[i]),'\n')
    if(!file.exists(paste(raw.dir,'/wss_SSA_',NRCS.areas@data$areasymbol[i],'_[',as.Date(NRCS.areas@data$saverest[i]),'].zip', sep=''))){
      #       f = CFILE(paste(raw.dir,'/wss_SSA_',NRCS.areas@data$areasymbol[i],'_[',as.Date(NRCS.areas@data$saverest[i]),'].zip', sep=''), mode="wb")
      #       curlPerform(url = paste('http://websoilsurvey.sc.egov.usda.gov/DSD/Download/Cache/SSA/wss_SSA_',NRCS.areas@data$areasymbol[i],'_[',as.Date(NRCS.areas@data$saverest[i]),'].zip', sep=''), writedata = f@ref)
      #       close(f)
      
      system(paste("wget -np -N http://websoilsurvey.sc.egov.usda.gov/DSD/Download/Cache/SSA/wss_SSA_",NRCS.areas@data$areasymbol[i],"_[",as.Date(NRCS.areas@data$saverest[i]),"].zip --directory-prefix=",raw.dir, sep=''))
      
    }
    
    unzip(paste(raw.dir,'/wss_SSA_',NRCS.areas@data$areasymbol[i],'_[',as.Date(NRCS.areas@data$saverest[i]),'].zip', sep=''),exdir=raw.dir)
    
    NRCS.polys[[i]] <- readOGR(paste(raw.dir,'/wss_SSA_',NRCS.areas@data$areasymbol[i],'_[',as.Date(NRCS.areas@data$saverest[i]),']/spatial',sep=''), layer=paste("soilmu_a_",tolower(NRCS.areas@data$areasymbol[i]),sep=''), verbose=F)
    # Change all spatial IDs to prepare from merging
    NRCS.polys[[i]] <- spChFIDs(NRCS.polys[[i]], as.character(paste(NRCS.areas@data$areasymbol[i],'_',row.names(NRCS.polys[[i]]@data),sep='')))
    
    unlink(paste(raw.dir,'/wss_SSA_',NRCS.areas@data$areasymbol[i],'_[',as.Date(NRCS.areas@data$saverest[i]),']',sep=''), recursive=T, force=T)
    
  }
  
  # Merging all NRCS Map Unit polygons
  NRCS.polys <- do.call("rbind", NRCS.polys)
  
  # Crop to area of x
  cat("Cropping NRCS Map Unit polygons to the extent of x\n")
  NRCS.polys <- crop.to.studyArea(NRCS.polys,spTransform(x,CRS(projection(NRCS.polys))))
  
  # Make MUKEY column numeric, and not factor
  NRCS.polys$MUKEY <- as.numeric(as.character(NRCS.polys$MUKEY))
  
  NRCS.polys <- NRCS.polys[,c("MUKEY","AREASYMBOL")]
  
  return(NRCS.polys)
}


loadAndAggregateSoilTable <- function(x, NRCS.areas, raw.dir="../Input/NRCS"){
  tables <- vector("list",length(NRCS.areas))
  for(i in 1:length(NRCS.areas)){
    unzip(paste(raw.dir,'/wss_SSA_',NRCS.areas@data$areasymbol[i],'_[',as.Date(NRCS.areas@data$saverest[i]),'].zip', sep=''),exdir=raw.dir)
    
    if(length(readLines(paste(raw.dir,'/wss_SSA_',NRCS.areas@data$areasymbol[i],'_[',as.Date(NRCS.areas@data$saverest[i]),']/tabular/',x,'.txt',sep='')))>0){
      tables[[i]] <- read.delim(paste(raw.dir,'/wss_SSA_',NRCS.areas@data$areasymbol[i],'_[',as.Date(NRCS.areas@data$saverest[i]),']/tabular/',x,'.txt',sep=''), header=F,sep="|")
      tables[[i]]$area <- NRCS.areas@data$areasymbol[i]
    }
    
    unlink(paste(raw.dir,'/wss_SSA_',NRCS.areas@data$areasymbol[i],'_[',as.Date(NRCS.areas@data$saverest[i]),']',sep=''), recursive=T, force=T)
    
  }
  table <- do.call("rbind", tables)
}


loadHiResRaster <- function(x, dsn.rasters, aggFactor, force.redo=F){
  if(!file.exists(paste(dsn.rasters,'/soils_hiRes.tif',sep='')) | force.redo){
    NRCS.rast.hiRes <- createHiResRaster(x=x, aggFactor=aggFactor, dsn.rasters=dsn.rasters)
  }else{
    NRCS.rast.hiRes <- raster(paste(dsn.rasters,'/soils_hiRes.tif',sep=''))
  }
  return(NRCS.rast.hiRes)
}


createHiResRaster <- function(x, dsn.rasters, aggFactor){
  ############################################################################################################
  ### CREATE HI-RESOLUTION RASTER ###
  ## This is a method for generating a roughly weighted distribution of soil characteristics. The vector soils
  ## are rasterized onto a landscape that has a higher resolution than the final landscape (16 hi-res cells to 
  ## one model cell). All soil characteristics are first calculated at this higher resolution, then the final 
  ## rasters are aggregated to the model resolution, leading to a weighted distribution.
  cat("Creating high-resolution raster for rasterization of soil attributes.\n")
  
  NRCS.polys <- loadNRCSMapUnitPolygons(x=x, raw.dir=raw.dir, dsn.vectors=dsn.vectors, force.redo=force.redo)
  
  # Create a higher-resolution blank raster.
  sim.raster.hiRes <- raster(ext=extent(x),nrows=nrow(x)*aggFactor,ncols=ncol(x)*aggFactor,crs=CRS(projection(x)))
  
  # Create a raster version of the vector dataset.
  cat("Rasterizing soils by Map Unit keys.\n")
  NRCS.rast.hiRes <- raster::rasterize(NRCS.polys,sim.raster.hiRes,field="MUKEY")
  
  # Write rasters in GeoTiff format.
  cat("Writing high-resolution raster in GeoTiff format.\n")
  writeRaster(NRCS.rast.hiRes,paste(dsn.rasters,'/soils_hiRes',sep=''),format='GTiff',overwrite=TRUE)
  
  return(NRCS.rast.hiRes)
}


getSoilData <- function(x, areas, polys, raw.dir="../Input/NRCS", dsn.vectors="Output/vectors", tables.dir="Output/tables", force.redo=F){
  
  dbcon <- dbConnect(dbDriver("SQLite"), dbname = paste(raw.dir,"/soildb_US_2003.sqlite",sep=''))
  dbListTables(dbcon)
  
  if(!file.exists(paste(tables.dir,"/muaggatt.csv",sep='')) | force.redo){
    # Load the USGS survey data from the map units specified.
    NRCS.muaggatt <- loadAndAggregateSoilTable("muaggatt",areas, raw.dir=raw.dir)
    names(NRCS.muaggatt) <- c(dbListFields(dbcon, "muaggatt"),"area")
    #     NRCS.muaggatt <- NRCS.muaggatt[NRCS.muaggatt$mukey %in% unique(polys$MUKEY),]
    write.csv(NRCS.muaggatt,paste(tables.dir,"/muaggatt.csv",sep=''),row.names=F)
  }
  
  # Component data.
  if(!file.exists(paste(tables.dir,"/comp.csv",sep='')) | force.redo){
    NRCS.comp <- loadAndAggregateSoilTable("comp",areas, raw.dir=raw.dir)
    names(NRCS.comp)<- c(dbListFields(dbcon, "component"),"area")
    #     NRCS.comp <- NRCS.comp[NRCS.comp$mukey %in% unique(polys$MUKEY),]
    NRCS.comp <- correctSoilComponents(NRCS.comp)
    write.csv(NRCS.comp,paste(tables.dir,"/comp.csv",sep=''),row.names=F)
  }
  
  # Horizon data
  if(!file.exists(paste(tables.dir,"/chorizon.csv",sep='')) | force.redo){
    NRCS.chorizon <- loadAndAggregateSoilTable("chorizon",areas, raw.dir=raw.dir)
    names(NRCS.chorizon)<- c(dbListFields(dbcon, "chorizon"),"area")
    #     NRCS.chorizon <- NRCS.chorizon[NRCS.chorizon$cokey %in% unique(NRCS.comp$cokey),]
    write.csv(NRCS.chorizon,paste(tables.dir,"/chorizon.csv",sep=''),row.names=F)
  }  
  
  # Crop yield data
  if(!file.exists(paste(tables.dir,"/ccrpyd.csv",sep='')) | force.redo){
    NRCS.ccrpyd <- loadAndAggregateSoilTable("ccrpyd",areas, raw.dir=raw.dir)
    names(NRCS.ccrpyd)<- c(dbListFields(dbcon, "cocropyld"),"area")
    #     NRCS.chorizon <- NRCS.chorizon[NRCS.chorizon$cokey %in% unique(NRCS.comp$cokey),]
    write.csv(NRCS.ccrpyd,paste(tables.dir,"/ccrpyd.csv",sep=''),row.names=F)
  }  
  
  # Load the primary mapunit data
  if(!file.exists(paste(tables.dir,"/mapunit.csv",sep='')) | force.redo){
    NRCS.mapunit <- loadAndAggregateSoilTable("mapunit",areas, raw.dir=raw.dir)
    names(NRCS.mapunit) <- c(dbListFields(dbcon, "mapunit"),"area")
    #     NRCS.mapunit <- NRCS.mapunit[NRCS.mapunit$mukey %in% unique(polys$MUKEY),]
    write.csv(NRCS.mapunit,paste(tables.dir,"/mapunit.csv",sep=''),row.names=F)
  }
  
  # Horizon fragments data
  if(!file.exists(paste(tables.dir,"/chfrags.csv",sep='')) | force.redo){
    NRCS.chfrags <- loadAndAggregateSoilTable("chfrags",areas, raw.dir=raw.dir)
    names(NRCS.chfrags)<- c(dbListFields(dbcon, "chfrags"),"area")
    #     NRCS.chfrags <- NRCS.chfrags[NRCS.chfrags$chkey %in% unique(NRCS.chorizon$chkey),]
    write.csv(NRCS.chfrags,paste(tables.dir,"/chfrag.csv",sep=''),row.names=F)
  }  
}


correctSoilComponents <- function(x){
  out <- lapply(unique(x$mukey), function(i){    temp <- x[x$mukey==i,]; temp$comppct_r_corr <- temp$comppct_r/sum(temp$comppct_r); return(temp)})
  out <- do.call("rbind",out)
  return(out)
}


loadSoilData <- function(x, raw.dir="../Input/NRCS", dsn.vectors="Output/NRCS_vectors.gdb", tables.dir="Output/tables", force.redo=F){
  if(!file.exists(paste(tables.dir,"/",x,".csv",sep=''))){
    getSoilData(x, raw.dir=raw.dir, dsn.vectors=dsn.vectors, tables.dir=tables.dir, force.redo=F)
  }
  soil.data <- read.csv(paste(tables.dir,"/",x,".csv",sep=''))
  return(soil.data)
}

# This function splits the soil component at a 30cm depth, and calculates upper and lower soil size particle classes.
calculateTexture <- function(raw.dir="../Input/NRCS", dsn.vectors="Output/NRCS_vectors.gdb", tables.dir="Output/tables"){
  NRCS.chorizon <- loadSoilData("chorizon", raw.dir=raw.dir, dsn.vectors=dsn.vectors, tables.dir=tables.dir, force.redo=F)
  NRCS.comp <- loadSoilData("comp", raw.dir=raw.dir, dsn.vectors=dsn.vectors, tables.dir=tables.dir, force.redo=F)
  
  # Some soil particle size classes are only calculated for two of the three size classes; derive the third.
  # Or, if only calculated for one, make all NA. This is usually the case for water-logged clays in riverbeds.
  sanitize.textures <- function(x){
    out <- x
    nas <- names(x[which(is.na(x))])
    
    if(length(nas)==1){
      out[names(out)%in%nas] <- 100-sum(out[!(names(out) %in% nas)])
    }
    
    if(length(nas)==2){
      out[!(names(out) %in% nas)] <- NA
    }
    
    return(out)
    
  }
  
  NRCS.chorizon[,c("sandtotal_r","silttotal_r","claytotal_r")] <- t(apply(NRCS.chorizon[,c("sandtotal_r","silttotal_r","claytotal_r")],1,sanitize.textures))
  
  NRCS.chorizon <- NRCS.chorizon[order(NRCS.chorizon$hzdept_r),]
  NRCS.chorizon <- NRCS.chorizon[order(NRCS.chorizon$cokey),]
  
  
  
  # If a horizon spans the 30 cm depth, split it in two at 30 cm 
  for (i in 1:length(NRCS.chorizon$cokey)){
    if(NRCS.chorizon$hzdept_r[i]<30 & NRCS.chorizon$hzdepb_r[i]>30){
      temp <- NRCS.chorizon[i,]
      temp$hzdept_r <- 30
      NRCS.chorizon[i,]$hzdepb_r <- 30
      NRCS.chorizon <- rbind(NRCS.chorizon,temp)
    }
  }
  # Sort horizons by top depths and component keys
  NRCS.chorizon <- NRCS.chorizon[order(NRCS.chorizon$hzdept_r),]
  NRCS.chorizon <- NRCS.chorizon[order(NRCS.chorizon$cokey),]
  
  
  
  # Calculate final SAND, SILT, and CLAY values for each horizon
  NRCS.chorizon$sandtotal_r <- (NRCS.chorizon$hzdepb_r-NRCS.chorizon$hzdept_r)*NRCS.chorizon$sandtotal_r
  NRCS.chorizon$silttotal_r <- (NRCS.chorizon$hzdepb_r-NRCS.chorizon$hzdept_r)*NRCS.chorizon$silttotal_r
  NRCS.chorizon$claytotal_r <- (NRCS.chorizon$hzdepb_r-NRCS.chorizon$hzdept_r)*NRCS.chorizon$claytotal_r
  # Code horizons as being either 0-30 (level 1) or >30 (level 2) cm in depth
  NRCS.chorizon$level[NRCS.chorizon$hzdepb_r<=30] <- 1
  NRCS.chorizon$level[NRCS.chorizon$hzdepb_r>30] <- 2
  
  NRCS.chorizon[is.na(NRCS.chorizon$sandtotal_r) & NRCS.chorizon$level==1,]$level <- 0
  NRCS.chorizon[is.na(NRCS.chorizon$sandtotal_r) & NRCS.chorizon$level==2,]$level <- 3
  
  #   NRCS.chorizon <- NRCS.chorizon[!is.na(NRCS.chorizon$SAND),]
  
  # Sum counts by level
  NRCS.chorizon.final.counts <- aggregate(NRCS.chorizon[,c("sandtotal_r","silttotal_r","claytotal_r")],by=list(NRCS.chorizon$cokey,NRCS.chorizon$level),FUN=sum, na.pass=T)
  names(NRCS.chorizon.final.counts) <- c("cokey","level","sandtotal_r","silttotal_r","claytotal_r")
  NRCS.chorizon.final.counts <- NRCS.chorizon.final.counts[order(NRCS.chorizon.final.counts$level),]
  NRCS.chorizon.final.counts <- NRCS.chorizon.final.counts[order(NRCS.chorizon.final.counts$cokey),]
  
  NRCS.chorizon.final.depths.top <- aggregate(NRCS.chorizon[,"hzdept_r"],by=list(NRCS.chorizon$cokey,NRCS.chorizon$level),FUN=min, na.rm=TRUE)
  NRCS.chorizon.final.depths.bottom <- aggregate(NRCS.chorizon[,"hzdepb_r"],by=list(NRCS.chorizon$cokey,NRCS.chorizon$level),FUN=max, na.rm=TRUE)
  names(NRCS.chorizon.final.depths.top) <- c("cokey","level","TOP")
  names(NRCS.chorizon.final.depths.bottom) <- c("cokey","level","BOTTOM")
  
  NRCS.chorizon.final.depths.top <- NRCS.chorizon.final.depths.top[order(NRCS.chorizon.final.depths.top$level),]
  NRCS.chorizon.final.depths.top <- NRCS.chorizon.final.depths.top[order(NRCS.chorizon.final.depths.top$cokey),]
  
  NRCS.chorizon.final.depths.bottom <- NRCS.chorizon.final.depths.bottom[order(NRCS.chorizon.final.depths.bottom$level),]
  NRCS.chorizon.final.depths.bottom <- NRCS.chorizon.final.depths.bottom[order(NRCS.chorizon.final.depths.bottom$cokey),]
  
  NRCS.chorizon.final.counts <- merge(NRCS.chorizon.final.counts,NRCS.chorizon.final.depths.top, all=T)
  NRCS.chorizon.final.counts <- merge(NRCS.chorizon.final.counts,NRCS.chorizon.final.depths.bottom, all=T)
  
  
  COMP_TOTAL <- NRCS.chorizon.final.counts$sandtotal_r+NRCS.chorizon.final.counts$silttotal_r+NRCS.chorizon.final.counts$claytotal_r
  NRCS.chorizon.final.counts$SAND <- NRCS.chorizon.final.counts$sandtotal_r/COMP_TOTAL
  NRCS.chorizon.final.counts$SILT <- NRCS.chorizon.final.counts$silttotal_r/COMP_TOTAL
  NRCS.chorizon.final.counts$CLAY <- NRCS.chorizon.final.counts$claytotal_r/COMP_TOTAL
  
  NRCS.chorizon.final.counts <- NRCS.chorizon.final.counts[NRCS.chorizon.final.counts$level %in% c(1,2),]
  
  #   NRCS.chorizon.final.counts <- NRCS.chorizon.final.counts[NRCS.chorizon.final.counts$level==2 | (NRCS.chorizon.final.counts$level==1 & NRCS.chorizon.final.counts$TOP == 0)| (NRCS.chorizon.final.counts$level==1 & NRCS.chorizon.final.counts$BOTTOM == 30),]
  #   NRCS.chorizon.final.counts <- NRCS.chorizon.final.counts[NRCS.chorizon.final.counts$level==1 | (NRCS.chorizon.final.counts$level==2 & NRCS.chorizon.final.counts$BOTTOM >= 40),]
  #   
  NRCS.chorizon.final.counts[(NRCS.chorizon.final.counts$level==1),]$TOP <- 0
  dead.cokeys.1 <- NRCS.chorizon.final.counts[(NRCS.chorizon.final.counts$level==1 & NRCS.chorizon.final.counts$BOTTOM != 30),]$cokey
  dead.cokeys.2 <- NRCS.chorizon.final.counts[NRCS.chorizon.final.counts$level==2 & NRCS.chorizon.final.counts$BOTTOM < 40,]$cokey
  NRCS.chorizon.final.counts <- NRCS.chorizon.final.counts[!(NRCS.chorizon.final.counts$cokey %in% dead.cokeys.1),]
  NRCS.chorizon.final.counts <- NRCS.chorizon.final.counts[!(NRCS.chorizon.final.counts$cokey %in% dead.cokeys.2),]
  
  NRCS.comp <- merge(NRCS.comp,NRCS.chorizon.final.counts,all=T)
  NRCS.comp <- NRCS.comp[order(NRCS.comp$level),]
  NRCS.comp <- NRCS.comp[order(NRCS.comp$cokey),]
  NRCS.comp <- NRCS.comp[order(NRCS.comp$mukey),]
  
  # Weight each component texture classes by the component contribution to the map unit
  NRCS.comp$SAND_P <- NRCS.comp$comppct_r_corr*NRCS.comp$SAND
  NRCS.comp$SILT_P <- NRCS.comp$comppct_r_corr*NRCS.comp$SILT
  NRCS.comp$CLAY_P <- NRCS.comp$comppct_r_corr*NRCS.comp$CLAY
  
  # ... and by their layer thickness
  NRCS.comp$SAND_P <- (NRCS.comp$BOTTOM-NRCS.comp$TOP)*NRCS.comp$SAND_P
  NRCS.comp$SILT_P <- (NRCS.comp$BOTTOM-NRCS.comp$TOP)*NRCS.comp$SILT_P
  NRCS.comp$CLAY_P <- (NRCS.comp$BOTTOM-NRCS.comp$TOP)*NRCS.comp$CLAY_P
  
  NRCS.comp.final.depths.top <- aggregate(NRCS.comp$TOP,by=list(NRCS.comp$mukey,NRCS.comp$level),FUN=mean, na.rm=TRUE)
  NRCS.comp.final.depths.bottom <- aggregate(NRCS.comp$BOTTOM,by=list(NRCS.comp$mukey,NRCS.comp$level),FUN=mean, na.rm=TRUE)
  names(NRCS.comp.final.depths.top) <- c("mukey","level","TOP")
  names(NRCS.comp.final.depths.bottom) <- c("mukey","level","BOTTOM")
  NRCS.comp.final.depths <- merge(NRCS.comp.final.depths.top,NRCS.comp.final.depths.bottom, all=T)
  
  NRCS.comp <- NRCS.comp[!is.na(NRCS.comp$level),]
  
  NRCS.comp.final.counts <- aggregate(NRCS.comp[,c("comppct_r_corr","SAND_P","SILT_P","CLAY_P")],by=list(NRCS.comp$mukey,NRCS.comp$level),FUN=sum, na.rm=T)
  names(NRCS.comp.final.counts) <- c("mukey","level","mukey_p","SAND","SILT","CLAY")
  
  NRCS.comp.final <- merge(NRCS.comp.final.counts,NRCS.comp.final.depths, all=T)
  MU_TOTAL <- NRCS.comp.final$SAND+NRCS.comp.final$SILT+NRCS.comp.final$CLAY
  NRCS.comp.final$SAND <- NRCS.comp.final$SAND/MU_TOTAL
  NRCS.comp.final$SILT <- NRCS.comp.final$SILT/MU_TOTAL
  NRCS.comp.final$CLAY <- NRCS.comp.final$CLAY/MU_TOTAL
  
  NRCS.comp.final <- NRCS.comp.final[order(NRCS.comp.final$level),]
  NRCS.comp.final <- NRCS.comp.final[order(NRCS.comp.final$mukey,NRCS.comp.final$level),]
  
  return(NRCS.comp.final)
}

loadSFNF_Data <- function(x, areas, polys, SFNF.dir, raw.dir, dsn.vectors, force.redo){
  
  if(!("NM678" %in% polys$AREASYMBOL) | force.redo){
    SFNF.soils <- readOGR(paste(SFNF.dir,"/spatial", sep=''), layer='TEU_Land_Type')
    SFNF.soils <- spTransform(SFNF.soils,CRS(projection(polys)))
    
    SFNF.soils <- crop.to.studyArea(SFNF.soils,areas[areas$areasymbol=="NM678",])
    
    SFNF.soils$MAP_UNIT_S <- paste("SFNF_",SFNF.soils$MAP_UNIT_S,sep='')
    SFNF.soils <- SFNF.soils[,c('MAP_UNIT_S')]
    SFNF.soils$AREASYMBOL <- "NM678"
    
    names(SFNF.soils) <- c("MUKEY","AREASYMBOL")
    
    SFNF.soils <- spChFIDs(SFNF.soils, as.character(paste(SFNF.soils@data$AREASYMBOL,'_',row.names(SFNF.soils@data),sep='')))
    
    polys <- rbind(polys,SFNF.soils)
    suppressWarnings(writeOGR(polys, dsn.vectors, "soils","ESRI Shapefile", overwrite_layer=TRUE))
  }
  
  SFNF.soils.tabular <- read.csv(paste(SFNF.dir,"/SFmapLeg.csv",sep=''))
  SFNF.soils.tabular <- SFNF.soils.tabular[,c("MUNTNUMB","STEXTRCL","APCTCOMP")]
  names(SFNF.soils.tabular) <- c("MUKEY","TEXTURE","P_COMP")
  SFNF.soils.tabular$MUKEY <- paste("SFNF_",SFNF.soils.tabular$MUKEY,sep='')
  SFNF.soils.tabular <- SFNF.soils.tabular[SFNF.soils.tabular$MUKEY %in% polys$MUKEY,]
  SFNF.soils.tabular$TEXTURE <- as.character(SFNF.soils.tabular$TEXTURE)
  SFNF.soils.tabular$TEXTURE <- gsub("fine sandy loam","sandy loam",SFNF.soils.tabular$TEXTURE)
  SFNF.soils.tabular$TEXTURE <- gsub("silt loam","silty loam",SFNF.soils.tabular$TEXTURE)
  SFNF.soils.tabular$TEXTURE[SFNF.soils.tabular$TEXTURE=="---" | SFNF.soils.tabular$TEXTURE==""] <- NA
  
  # Get centroid CSS values on USDA texture triangle
  classes <- TT.classes.tbl(class.sys = "USDA.TT")
  vertices <- TT.css2xy(TT.vertices.tbl(class.sys = "USDA.TT")*100, geo=TT.geo.get(class.sys = "USDA.TT"))
  centroids <- t(apply(classes,1,FUN=function(x){sub.verts <- vertices[as.numeric(strsplit(as.character(x[3]),",")[[1]]),]; return(TT.polygon.centroids(sub.verts$xpos,sub.verts$ypos))}))
  colnames(centroids) <- c("xpos","ypos")
  centroids.css <- as.data.frame(TT.xy2css(centroids,geo=TT.geo.get(class.sys = "USDA.TT"))/100)
  centroids.css$TEXTURE <- classes[,2]
  
  SFNF.soils.tabular[!(SFNF.soils.tabular$P_COMP %in% centroids.css$TEXTURE),]
  
  SFNF.soils.tabular <- merge(SFNF.soils.tabular,centroids.css)[,-1]
  SFNF.soils.tabular <- SFNF.soils.tabular[order(SFNF.soils.tabular$MUKEY),]
  
  SFNF.soils.tabular$SAND <- SFNF.soils.tabular$SAND * SFNF.soils.tabular$P_COMP
  SFNF.soils.tabular$SILT <- SFNF.soils.tabular$SILT * SFNF.soils.tabular$P_COMP
  SFNF.soils.tabular$CLAY <- SFNF.soils.tabular$CLAY * SFNF.soils.tabular$P_COMP
  SFNF.soils.tabular <- aggregate(SFNF.soils.tabular[,c(2:5)],by=list(SFNF.soils.tabular$MUKEY),FUN=sum, na.rm=T)
  
  MU_TOTAL <- SFNF.soils.tabular$SAND+SFNF.soils.tabular$SILT+SFNF.soils.tabular$CLAY
  SFNF.soils.tabular$SAND <- SFNF.soils.tabular$SAND/MU_TOTAL
  SFNF.soils.tabular$SILT <- SFNF.soils.tabular$SILT/MU_TOTAL
  SFNF.soils.tabular$CLAY <- SFNF.soils.tabular$CLAY/MU_TOTAL
  names(SFNF.soils.tabular) <- c("mukey","mukey_p","CLAY","SILT","SAND")
  
  SFNF.soils.tabular.1 <- SFNF.soils.tabular
  SFNF.soils.tabular.1$level <- 1
  SFNF.soils.tabular.1$TOP <- 0
  SFNF.soils.tabular.1$BOTTOM <- 30
  
  SFNF.soils.tabular.2 <- SFNF.soils.tabular
  SFNF.soils.tabular.2$level <- 2
  SFNF.soils.tabular.2$TOP <- 30
  SFNF.soils.tabular.2$BOTTOM <- 152
  
  SFNF.soils.tabular <- rbind(SFNF.soils.tabular.1,SFNF.soils.tabular.2)
  SFNF.soils.tabular <- SFNF.soils.tabular[order(SFNF.soils.tabular$mukey,SFNF.soils.tabular$level),]
  
  return(SFNF.soils.tabular)
}

textureClass.rasters <- function(SAND.rast, SILT.rast, CLAY.rast, class.sys = "USDA.TT"){
  texture.css <- data.frame(SAND=getValues(SAND.rast)*100, SILT=getValues(SILT.rast)*100, CLAY=getValues(CLAY.rast)*100)
  texture.css$CELL <- row.names(texture.css)
  texture.css.incomplete <- texture.css[!complete.cases(texture.css),]
  texture.css.complete <- texture.css[complete.cases(texture.css),]
  texture.complete <- data.frame(CELL=as.numeric(texture.css.complete$CELL), TEXTURE=TT.points.in.classes(texture.css.complete, class.sys = class.sys, PiC.type="t"))
  texture.incomplete <- data.frame(CELL=as.numeric(texture.css.incomplete$CELL), TEXTURE=NA) 
  texture <- rbind(texture.complete,texture.incomplete)
  #   texture$CELL <- as.numeric(texture$CELL)
  texture <- texture[order(texture$CELL),]
  return(texture)
}

fillReservoirSoils <- function(gapped.soil.raster, dem.raster, label, raw.dir, force.redo=F){
  if(file.exists(paste(raw.dir,"/EXTRACTIONS/",label,"/",'RASTERIZED_MUKEYS_1arcsec_filled.tif',sep='')) & !force.redo){
    NRCS.rast.hiRes <- raster(readGDAL(paste(raw.dir,"/EXTRACTIONS/",label,"/",'RASTERIZED_MUKEYS_1arcsec_filled.tif',sep=''), silent=T))
    return(NRCS.rast.hiRes)
  }
  
  # Create a RAM disk for quick read/write in GRASS
  system("diskutil erasevolume HFS+ 'ramdisk' `hdiutil attach -nomount ram://16384000`")
  # Create GRASS things on RAM disk
  system("mkdir /Volumes/ramdisk/home")
  system("mkdir /Volumes/ramdisk/gisDbase")
  
  # Location of your GRASS installation:
  loc <- initGRASS("/usr/local/opt/grass-64/grass-6.4.4", home="/Volumes/ramdisk/home", gisDbase="/Volumes/ramdisk/gisDbase",mapset="PERMANENT",override=TRUE)
  gmeta6()
  gc()
  
  # # Compute the Topographic Wetness Index ("TWI") for the DEM
  # First, export the hiRes DEM as a GRASS raster dataset
  writeRAST6(as(dem.raster,"SpatialGridDataFrame"), flags=c("o","overwrite","quiet"), "DEM", useGDAL=TRUE)
  # Then set the GRASS mapping region to the DEM
  execGRASS("g.region", flags=c("quiet"), parameters=list(rast="DEM"))
  
  # Then, run the r.topidx function on the GRASS dataset
  execGRASS("r.topidx", flags=c("quiet","overwrite"), parameters=list(input="DEM",output="TWI"))
  twi.hiRes <- raster(readRAST6("TWI",ignore.stderr = TRUE))
  twi.hiRes[is.na(twi.hiRes)] <- 25
  
  # Discard GRASS
  system("hdiutil detach /volumes/ramdisk")
  
  # Calculate slope, aspect, and flow direction (in degrees) on the hiRes DEM
  slope.hiRes <- terrain(dem.raster, opt='slope', unit='degrees')
  aspect.hiRes <- terrain(dem.raster, opt='aspect', unit='degrees')
  slope.hiRes[is.na(slope.hiRes)] <- 0
  aspect.hiRes[is.na(aspect.hiRes)] <- 90
  
  # Aspect has to be split into sine and cosine components, because it is
  # circular.
  aspect.sine.hiRes <- sin_d(aspect.hiRes)
  aspect.cosine.hiRes <- cos_d(aspect.hiRes)
  
  # Coerce the rasters into attribute tables
  
  dem.hiRes.sp <- raster::as.data.frame(dem.raster, xy=T)
  nrcs.hiRes.sp <- raster::as.data.frame(gapped.soil.raster, xy=T)
  slope.hiRes.sp <- raster::as.data.frame(slope.hiRes, xy=T)
  aspect.sine.hiRes <- raster::as.data.frame(aspect.sine.hiRes, xy=T)
  aspect.cosine.hiRes <- raster::as.data.frame(aspect.cosine.hiRes, xy=T)
  twi.hiRes.sp <- raster::as.data.frame(twi.hiRes, xy=T)
  
  # Merge all tables together based on easting and northing
  merged.sp <- data.frame(x=dem.hiRes.sp$x, y=dem.hiRes.sp$y, elevation=dem.hiRes.sp$DEM_1_DRAINED, mukey=nrcs.hiRes.sp$layer, slope=slope.hiRes.sp$slope, aspect_sine=aspect.sine.hiRes$layer, aspect_cosine=aspect.cosine.hiRes$layer, twi=twi.hiRes.sp$TWI)
  
  # Rename merged variables to clarify
  names(merged.sp) <- c("easting","northing","elevation","mukey","slope","aspect_sin","aspect_cos","twi")
  
  merged.sp.nonNull <- merged.sp[!is.na(merged.sp$mukey),]
  merged.sp.null <- merged.sp[is.na(merged.sp$mukey),]
  
  # MUKEYs are a factor. Make them a factor
  merged.sp.nonNull$mukey <- as.character(merged.sp.nonNull$mukey)
  
  merged.sp.nonNull.prior <- merged.sp.nonNull
  
  # Perform the linear discriminant analysis
  merged.sp.nonNull.lda <- lda(mukey~easting+northing+elevation+slope+aspect_sin+aspect_cos+twi, data= merged.sp.nonNull)
  
  # Predict values for NULL soils, and append to soils.null
  merged.sp.null.predict <- predict(merged.sp.nonNull.lda,merged.sp.null[,c(1,2,3,5:8)])
  merged.sp.null$mukey <- merged.sp.null.predict$class
  
  # Convert MUKEYs to character for merging
  merged.sp.null$mukey2 <- as.character(merged.sp.null$mukey)
  merged.sp.nonNull$mukey2 <- as.character(merged.sp.nonNull$mukey)
  
  # Extract just the spatial and MUKEY data
  merged.sp.null <- merged.sp.null[,c(1,2,9)]
  merged.sp.nonNull <- merged.sp.nonNull[,c(1,2,9)]
  
  # Merge the null, nonNull, and edge datasets
  merged.sp.correct <- rbind(merged.sp.null, merged.sp.nonNull)
  
  # Change MUKEY to numeric for rasterization
  merged.sp.correct$mukey2 <- as.numeric(merged.sp.correct$mukey2)
  
  # Rasterize, first by making a SpatialPointsDataFrame, then by gridding it
  coordinates(merged.sp.correct) <- ~easting+northing
  gridded(merged.sp.correct) <- TRUE
  names(merged.sp.correct) <- "mukey"
  
  # Generate raster object, and write it
  NRCS.rast.hiRes <- raster(merged.sp.correct)
  writeGDAL(as(NRCS.rast.hiRes, "SpatialGridDataFrame"),paste(raw.dir,"/EXTRACTIONS/",label,"/",'RASTERIZED_MUKEYS_1arcsec_filled.tif',sep=''), drivername="GTiff", type="Int16", mvFlag=-32768, options=c("INTERLEAVE=PIXEL", "COMPRESS=DEFLATE", "ZLEVEL=9"))
  return(NRCS.rast.hiRes)
}

