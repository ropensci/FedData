#' Download and crop data from the NRCS SSURGO soils database.
#'
#' This is an efficient method for spatially merging several different soil survey areas
#' as well as merging their tabular data.
#' 
#' \code{get_ssurgo} returns a named list of length 2:
#' \enumerate{
#' \item 'spatial': A \code{SpatialPolygonsDataFrame} of soil mapunits 
#' in the template, and 
#' \item 'tabular': A named list of \code{\link{data.frame}s} with the SSURGO tabular data.
#' }
#' 
#' @param template A Raster* or Spatial* object to serve 
#' as a template for cropping; optionally, a vector of area names [e.g., c('IN087','IN088')] may be provided.
#' @param label A character string naming the study area.
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing. Defaults to './RAW/SSURGO/'.
#' @param extraction.dir A character string indicating where the extracted and cropped SSURGO shapefiles should be put.
#' The directory will be created if missing. Defaults to './EXTRACTIONS/SSURGO/'.
#' @param force.redo If an extraction for this template and label already exists, should a new one be created? Defaults to FALSE.
#' @return A named list containing the 'spatial' and 'tabular' data.
#' @export
#' @importFrom sp SpatialPointsDataFrame
#' @importFrom readr read_csv write_csv
#' @importFrom rgeos gIntersection
#' @examples
#' \dontrun{
#' # Extract data for the Village Ecodynamics Project 'VEPIIN' study area:
#' # http://village.anth.wsu.edu
#' vepPolygon <- polygon_from_extent(raster::extent(672800,740000,4102000,4170000), 
#'      proj4string='+proj=utm +datum=NAD83 +zone=12')
#' 
#' # Get the NRCS SSURGO data (USA ONLY)
#' SSURGO.VEPIIN <- get_ssurgo(template=vepPolygon, label='VEPIIN')
#' 
#' # Plot the VEP polygon
#' plot(vepPolygon)
#' 
#' # Plot the SSURGO mapunit polygons
#' plot(SSURGO.VEPIIN$spatial, lwd=0.1, add=T)
#' 
#' # Or, download by Soil Survey Area names
#' SSURGO.areas <- get_ssurgo(template=c('CO670','CO075'), label='CO_TEST')
#' 
#' # Let's just look at spatial data for CO675
#' SSURGO.areas.CO675 <- SSURGO.areas$spatial[SSURGO.areas$spatial$AREASYMBOL=='CO075',]
#' 
#' # And get the NED data under them for pretty plotting
#' NED.CO675 <- get_ned(template=SSURGO.areas.CO675, label='SSURGO_CO675')
#' 
#' # Plot the SSURGO mapunit polygons, but only for CO675
#' plot(NED.CO675)
#' plot(SSURGO.areas.CO675, lwd=0.1, add=T)
#' }
get_ssurgo <- function(template,
                       label,
                       raw.dir = "./RAW/SSURGO",
                       extraction.dir = paste0("./EXTRACTIONS/", label, "/SSURGO"), 
                       force.redo = FALSE) {
  
  raw.dir <- normalizePath(paste0(raw.dir,"/."), mustWork = FALSE)  
  extraction.dir <- normalizePath(paste0(extraction.dir,"/."), mustWork = FALSE)
  
  dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)
  
  files <- list.files(extraction.dir)
  files <- files[grepl("csv", files)]
  files <- files[order(files)]
  files <- files[grepl(label, files)]
  
  if (!force.redo & length(files) > 0 & file.exists(paste(extraction.dir, "/", label, "_SSURGO_", "Mapunits.shp", sep = ""))) {
    SSURGOMapunits <- rgdal::readOGR(dsn = normalizePath(extraction.dir), layer = paste0(label, "_SSURGO_", "Mapunits"), verbose = F)
    
    tables <- lapply(files, function(file) {
      suppressMessages(readr::read_csv(paste(normalizePath(extraction.dir), "/", file, sep = ""), progress = F))
    })
    names(tables) <- gsub(".csv", "", files)
    names(tables) <- gsub(paste0(label, "_SSURGO_"), "", names(tables))
    
    return(list(spatial = SSURGOMapunits, tabular = tables))
  }
  
  if (class(template) == "character") {
    q <- paste0("SELECT areasymbol, saverest FROM sacatalog WHERE areasymbol IN (", paste(paste0("'", template, "'"), collapse = ","), 
                ");")
    SSURGOAreas <- soilDB::SDA_query(q)
    
    template.poly <- template
    
  } else {
    
    if (class(template) %in% c("RasterLayer", "RasterStack", "RasterBrick")) {
      template.poly <- spdf_from_polygon(sp::spTransform(polygon_from_extent(template), sp::CRS("+proj=longlat +ellps=GRS80")))
    } else if (class(template) %in% c("SpatialPoints", "SpatialPointsDataFrame")) {
      suppressWarnings(template.poly <- raster::buffer(polygon_from_extent(template), width = 1e-06))
    } else {
      template.poly <- template
    }
    
    # Get shapefile of SSURGO study areas in the template
    SSURGOAreas <- get_ssurgo_inventory(template = template.poly, raw.dir = raw.dir)
    # Remove SSURGO study areas that are not available
    SSURGOAreas <- SSURGOAreas[SSURGOAreas$iscomplete != 0, ]
    
  }
  
  # Get data for each study area
  SSURGOData <- lapply(1:nrow(SSURGOAreas), function(i) {
    message("(Down)Loading SSURGO data for survey area ", i, " of ", nrow(SSURGOAreas), ": ", as.character(SSURGOAreas$areasymbol[i]))
    get_ssurgo_study_area(template = template.poly, area = as.character(SSURGOAreas$areasymbol[i]), date = as.Date(SSURGOAreas$saverest[i], 
                                                                                                                   format = "%m/%d/%Y"), raw.dir = raw.dir)
  })
  
  # Combine mapunits
  SSURGOPolys <- lapply(SSURGOData, "[[", "spatial")
  
  # Merging all SSURGO Map Unit polygons
  message("Merging all SSURGO Map Unit polygons")
  SSURGOPolys <- do.call("rbind", SSURGOPolys)
  
  # Crop to area of template
  if (!is.null(template) & !is.character(template)) {
    message("Cropping all SSURGO Map Unit polygons to template")
    if (class(template) %in% c("SpatialPoints", "SpatialPointsDataFrame")) {
      SSURGOPolys <- sp::spTransform(template, sp::CRS(raster::projection(SSURGOPolys))) %over% SSURGOPolys
      SSURGOPolys <- SpatialPointsDataFrame(template@coords, data = SSURGOPolys, proj4string = sp::CRS(raster::projection(template)))
    } else {
      SSURGOPolys <- raster::crop(SSURGOPolys, sp::spTransform(template.poly, sp::CRS(raster::projection(SSURGOPolys))))
    }
  }
  
  # Combine study area data
  SSURGOTables <- lapply(SSURGOData, "[[", "tabular")
  
  # Merging all SSURGO data tables
  message("Merging all SSURGO data tables")
  tableNames <- unique(unlist(sapply(SSURGOTables, names)))
  tableNames <- tableNames[order(tableNames)]
  
  # This function takes each table name, gets that table from each study area, and binds the rows of those tables.  Finally, it
  # removes any duplicate lines.
  SSURGOTables <- lapply(tableNames, function(name) {
    tables <- lapply(SSURGOTables, "[[", name)
    tables <- do.call("rbind", tables)
    tables <- unique(tables)
    return(tables)
  })
  
  names(SSURGOTables) <- tableNames
  
  # Extract only the mapunits in the study area, and iterate through the data structure
  SSURGOTables <- extract_ssurgo_data(tables = SSURGOTables, mapunits = as.character(unique(SSURGOPolys$MUKEY)))
  
  # Save the mapunit polygons
  suppressWarnings(rgdal::writeOGR(SSURGOPolys, 
                                   dsn = normalizePath(paste0(extraction.dir,"/.")),
                                   layer = paste0(label, "_SSURGO_Mapunits"),
                                   driver = "ESRI Shapefile",
                                   overwrite_layer = TRUE))
  
  # Save the each data table as a csv
  junk <- lapply(names(SSURGOTables), function(tab) {
    readr::write_csv(SSURGOTables[[tab]], path = paste(extraction.dir, "/", label, "_SSURGO_", tab, ".csv", sep = ""))
  })
  
  return(list(spatial = SSURGOPolys, tabular = SSURGOTables))
}

#' Download a zipped directory containing a shapefile of the SSURGO study areas.
#'
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' @return A character string representing the full local path of the SSURGO study areas zipped directory.
#' @export
#' @keywords internal
download_ssurgo_inventory <- function(raw.dir) {
  # Import the shapefile of SSURGO study areas.  This is available at
  # http://soildatamart.sc.egov.usda.gov/download/StatusMaps/soilsa_a_SSURGO.zip
  url <- "http://websoilsurvey.sc.egov.usda.gov/DataAvailability/SoilDataAvailabilityShapefile.zip"
  destdir <- raw.dir
  download_data(url = url, destdir = destdir)
  return(normalizePath(paste(destdir, "/SoilDataAvailabilityShapefile.zip", sep = "")))
}

#' Download and crop a shapefile of the SSURGO study areas.
#'
#' \code{get_ssurgo_inventory} returns a \code{SpatialPolygonsDataFrame} of the SSURGO study areas within
#' the specified \code{template}. If template is not provided, returns the entire SSURGO inventory of study areas.
#' 
#' @param template A Raster* or Spatial* object to serve 
#' as a template for cropping.
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing.
#' @return A \code{SpatialPolygonsDataFrame} of the SSURGO study areas within
#' the specified \code{template}.
#' @importFrom methods as
#' @export
#' @keywords internal
get_ssurgo_inventory <- function(template = NULL, raw.dir) {
  # If there is a template, only download the areas in the template Thanks to Dylan Beaudette for this method!
  if (!is.null(template)) {
    
    if (class(template) %in% c("RasterLayer", "RasterStack", "RasterBrick")) {
      template <- sp::spTransform(polygon_from_extent(template), sp::CRS("+proj=longlat +datum=WGS84"))
    } else {
      template <- sp::spTransform(template, sp::CRS("+proj=longlat +datum=WGS84"))
    }
    
    bounds <- polygon_from_extent(template)
    
    # Only download 1 square degree at a time to avoid oversized AOI error
    if((raster::xmax(bounds) - raster::xmin(bounds)) > 1 | (raster::ymax(bounds) - raster::ymin(bounds)) > 1){
      grid <- sp::GridTopology(cellcentre.offset=c(-179.5,-89.5),
                               cellsize=c(1,1),
                               cells.dim=c(360,180)) %>%
        sp::SpatialGrid(proj4string=CRS("+proj=longlat +datum=WGS84")) %>%
        as("SpatialPolygons")
      
      bounds <- rgeos::gIntersection( grid, bounds, byid = T )
      
    }
    
    n = NULL
    SSURGOAreas <- foreach::foreach(n = 1:length(bounds),
                                    .combine = rbind) %do% {
      
      bound <- sp::bbox(bounds[n])
      if (identical(bound[1, 1], bound[1, 2])) 
        bound[1, 2] <- bound[1, 2] + 1e-04
      if (identical(bound[2, 1], bound[2, 2])) 
        bound[2, 2] <- bound[2, 2] + 1e-04
      
      bbox.text <- paste(bound, collapse = ",")
      
      url <- paste("https://sdmdataaccess.nrcs.usda.gov/Spatial/SDMNAD83Geographic.wfs?Service=WFS&Version=1.0.0&Request=GetFeature&Typename=SurveyAreaPoly&BBOX=", 
                   bbox.text, sep = "")
      
      temp.file <- paste0(tempdir(), "/soils.gml")
      
      opts <- list(verbose = F, noprogress = T, fresh_connect = TRUE)
      hand <- curl::new_handle()
      curl::handle_setopt(hand, .list = opts)
      tryCatch(status <- curl::curl_download(url, destfile = temp.file, handle = hand), error = function(e) stop("Download of ", 
                                                                                                                 url, " failed!"))
      
      SSURGOAreas <- rgdal::readOGR(dsn = temp.file, 
                                    layer = "surveyareapoly", 
                                    disambiguateFIDs = TRUE, 
                                    stringsAsFactors = FALSE, 
                                    verbose = FALSE)
      raster::projection(SSURGOAreas) <- raster::projection(template)
      
      # Get a list of SSURGO study areas within the project study area
      if (class(template) == "SpatialPointsDataFrame" & length(template) == 1) {
        template <- polygon_from_extent(bounds, proj4string = raster::projection(template))
      }
      SSURGOAreas <- raster::crop(SSURGOAreas, sp::spTransform(template, sp::CRS(raster::projection(SSURGOAreas))))
      
      SSURGOAreas <- SSURGOAreas@data
      
      SSURGOAreas$saverest <- as.Date(SSURGOAreas$saverest, format = "%b %d %Y")
      
      return(SSURGOAreas)
    }
    
    SSURGOAreas <- unique(SSURGOAreas)
    
  } else {
    
    tmpdir <- tempfile()
    if (!dir.create(tmpdir)) 
      stop("failed to create my temporary directory")
    
    file <- download_ssurgo_inventory(raw.dir = raw.dir)
    
    utils::unzip(file, exdir = tmpdir)
    
    SSURGOAreas <- rgdal::readOGR(normalizePath(tmpdir), layer = "soilsa_a_nrcs", verbose = FALSE)@data
    
    unlink(tmpdir, recursive = TRUE)
  }
  
  # Check to see if all survey areas are available
  if (0 %in% SSURGOAreas$iscomplete) {
    warning("Some of the soil surveys in your area are unavailable.\n
            Soils and productivity data will have holes.\n
            Missing areas:\n", 
            as.vector(SSURGOAreas[SSURGOAreas$iscomplete == 0, ]$areasymbol), "\n\n
            Continuing with processing available soils.\n\n")
  }
  
  return(SSURGOAreas)
}

#' Download a zipped directory containing the spatial and tabular data for a SSURGO study area.
#'
#' \code{download_ssurgo_study_area} first tries to download data including a state-specific Access
#' template, then the general US template.
#'
#' @param area A character string indicating the SSURGO study area to be downloaded.
#' @param date A character string indicating the date of the most recent update to the SSURGO 
#' area for these data. This information may be gleaned from the SSURGO Inventory (\code{\link{get_ssurgo_inventory}}).
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' @return A character string representing the full local path of the SSURGO study areas zipped directory.
#' @export
#' @keywords internal
download_ssurgo_study_area <- function(area, date, raw.dir) {
  
  # Try to download with the state database, otherwise grab the US
  url <- paste("http://websoilsurvey.sc.egov.usda.gov/DSD/Download/Cache/SSA/wss_SSA_", area, "_[", date, "].zip", sep = "")
  destdir <- raw.dir
  download_data(url = url, destdir = destdir, nc = T)
  
  return(normalizePath(paste(destdir, "/wss_SSA_", area, "_[", date, "].zip", sep = "")))
}

#' Download and crop the spatial and tabular data for a SSURGO study area.
#'
#' \code{get_ssurgo_study_area} returns a named list of length 2:
#' \enumerate{
#' \item 'spatial': A \code{SpatialPolygonsDataFrame} of soil mapunits 
#' in the template, and 
#' \item 'tabular': A named list of \code{\link{data.frame}s} with the SSURGO tabular data.
#' }
#' 
#' @param template A Raster* or Spatial* object to serve 
#' as a template for cropping. If missing, whose study area is returned
#' @param area A character string indicating the SSURGO study area to be downloaded.
#' @param date A character string indicating the date of the most recent update to the SSURGO 
#' area for these data. This information may be gleaned from the SSURGO Inventory (\code{\link{get_ssurgo_inventory}}).
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing.
#' @return A \code{SpatialPolygonsDataFrame} of the SSURGO study areas within
#' the specified \code{template}.
#' @export
#' @keywords internal
get_ssurgo_study_area <- function(template = NULL, area, date, raw.dir) {
  tmpdir <- tempfile()
  if (!dir.create(tmpdir)) 
    stop("failed to create my temporary directory")
  
  file <- download_ssurgo_study_area(area = area, date = date, raw.dir = raw.dir)
  
  utils::unzip(file, exdir = tmpdir)
  
  # Get spatial data
  if (.Platform$OS.type == "windows") {
    suppressWarnings(mapunits <- rgdal::readOGR(paste0(tmpdir, "/", area, "/spatial"), layer = paste0("soilmu_a_", tolower(area)), 
                                                verbose = F))
  } else {
    suppressWarnings(mapunits <- rgdal::readOGR(tmpdir, layer = paste0("soilmu_a_", tolower(area)), verbose = F))
  }
  
  # Crop to study area
  if (!is.null(template) & !is.character(template)) {
    if (class(template) %in% c("RasterLayer", "RasterStack", "RasterBrick")) {
      template <- spdf_from_polygon(sp::spTransform(polygon_from_extent(template), sp::CRS("+proj=longlat +ellps=GRS80")))
    }
    
    mapunits <- raster::crop(mapunits, sp::spTransform(template, sp::CRS(raster::projection(mapunits))))
  }
  
  # Change IDs, in case of merging later
  mapunits <- sp::spChFIDs(mapunits, as.character(paste(area, "_", row.names(mapunits@data), sep = "")))
  
  # Read in all tables
  if (.Platform$OS.type == "windows") {
    files <- list.files(paste0(tmpdir, "/", area, "/tabular"), full.names = T)
    tablesData <- lapply(files, function(file) {
      tryCatch(return(utils::read.delim(file, header = F, sep = "|", stringsAsFactors = F)), error = function(e) {
        return(NULL)
      })
    })
    names(tablesData) <- basename(files)
    tablesData <- tablesData[!sapply(tablesData, is.null)]
  } else {
    files <- list.files(tmpdir, full.names = T, pattern = "tabular")
    tablesData <- lapply(files, function(file) {
      tryCatch(return(utils::read.delim(file, header = F, sep = "|", stringsAsFactors = F)), error = function(e) {
        return(NULL)
      })
    })
    names(tablesData) <- gsub(paste0(area, "\\\\tabular\\\\"), "", basename(files))
    tablesData <- tablesData[!sapply(tablesData, is.null)]
  }
  
  # tablesHeaders <- FedData::tablesHeaders
  
  SSURGOTableMapping <- tablesData[["mstab.txt"]][, c(1, 5)]
  names(SSURGOTableMapping) <- c("TABLE", "FILE")
  SSURGOTableMapping[, "FILE"] <- paste(SSURGOTableMapping[, "FILE"], ".txt", sep = "")
  
  tablesData <- tablesData[as.character(SSURGOTableMapping[, "FILE"])]
  tablesHeads <- tablesHeaders[as.character(SSURGOTableMapping[, "TABLE"])]
  
  notNull <- (!sapply(tablesData, is.null) & !sapply(tablesHeads, is.null))
  tablesData <- tablesData[notNull]
  tablesHeads <- tablesHeads[notNull]
  
  tables <- mapply(tablesData, tablesHeads, FUN = function(theData, theHeader) {
    names(theData) <- names(theHeader)
    return(theData)
  })
  
  names(tables) <- names(tablesHeads)
  
  tables <- extract_ssurgo_data(tables = tables, mapunits = as.character(unique(mapunits$MUKEY)))
  
  unlink(tmpdir, recursive = TRUE)
  
  return(list(spatial = mapunits, tabular = tables))
}

#' Extract data from a SSURGO database pertaining to a set of mapunits.
#'
#' \code{extract_ssurgo_data} creates a directed graph of the joins in a SSURGO tabular dataset,
#' and then iterates through the tables, only retaining data pertinent to a set of mapunits.
#' 
#' @param tables A list of SSURGO tabular data.
#' @param mapunits A character vector of mapunits (likely dropped from SSURGO spatial data)
#' defining which mapunits to retain.
#' @return A list of extracted SSURGO tabular data.
#' @export
#' @keywords internal
extract_ssurgo_data <- function(tables, mapunits) {
  
  mapping <- tables[["mdstatrshipdet"]]
  mappingGraph <- igraph::graph.edgelist(as.matrix(mapping[, c("ltabphyname", "rtabphyname")]))
  igraph::E(mappingGraph)$mapVar <- as.character(mapping[, "ltabcolphyname"])
  
  mappingGraph <- igraph::graph.neighborhood(mappingGraph, order = max(sapply(igraph::decompose.graph(mappingGraph), igraph::diameter)) + 
                                               1, nodes = "mapunit", mode = "out")[[1]]
  mapHierarchy <- igraph::shortest.paths(mappingGraph, "mapunit")
  mapHierarchy <- colnames(mapHierarchy)[order(mapHierarchy)]
  mapHierarchy <- mapHierarchy[-1]
  mapEdges <- cbind(igraph::get.edgelist(mappingGraph), igraph::E(mappingGraph)$mapVar)
  mapEdges <- mapEdges[match(mapHierarchy, mapEdges[, 2]), ]
  
  tables[["mapunit"]] <- tables[["mapunit"]][tables[["mapunit"]][, "mukey"] %in% mapunits, ]
  
  for (i in 1:nrow(mapEdges)) {
    X <- mapEdges[i, ]
    tables[[X[2]]] <- tables[[X[2]]][tables[[X[2]]][, X[3]] %in% tables[[X[1]]][, X[3]], ]
  }
  
  tables <- tables[!sapply(tables, is.null)]
  
  return(tables)
}
