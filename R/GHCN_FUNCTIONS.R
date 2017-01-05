# globalVariables(c(".", "ELEMENT","YEAR","MONTH","LONGITUDE","LATITUDE"))

#' Download and crop the Global Historical Climate Network-Daily data.
#'
#' \code{get_ghcn_daily} returns a named list of length 2: 
#' \enumerate{
#' \item "spatial": A \code{SpatialPointsDataFrame} of the locations of GHCN weather stations 
#' in the template, and 
#' \item "tabular": A named list of \code{\link{data.frame}s} with the daily weather data for each station.
#' The name of each list item is the station ID.
#' }
#' 
#' @param template A Raster* or Spatial* object to serve 
#' as a template for cropping. Alternatively, a character vector providing GHCN station IDs. If missing, all stations
#' will be downloaded!
#' @param label A character string naming the study area.
#' @param elements A character vector of elemets to extract.\cr
#' The five core elements are:\cr
#' PRCP = Precipitation (tenths of mm)\cr
#' SNOW = Snowfall (mm)\cr
#' SNWD = Snow depth (mm)\cr
#' TMAX = Maximum temperature (tenths of degrees C)\cr
#' TMIN = Minimum temperature (tenths of degrees C)\cr
#' \cr
#' The other elements are:\cr
#'   
#' ACMC = Average cloudiness midnight to midnight from 30-second 
#' ceilometer data (percent)\cr
#' ACMH = Average cloudiness midnight to midnight from 
#' manual observations (percent)\cr
#' ACSC = Average cloudiness sunrise to sunset from 30-second 
#' ceilometer data (percent)\cr
#' ACSH = Average cloudiness sunrise to sunset from manual 
#' observations (percent)\cr
#' AWDR = Average daily wind direction (degrees)\cr
#' AWND = Average daily wind speed (tenths of meters per second)\cr
#' DAEV = Number of days included in the multiday evaporation
#' total (MDEV)\cr
#' DAPR = Number of days included in the multiday precipiation 
#' total (MDPR)\cr
#' DASF = Number of days included in the multiday snowfall 
#' total (MDSF)\cr
#' DATN = Number of days included in the multiday minimum temperature 
#' (MDTN)\cr
#' DATX = Number of days included in the multiday maximum temperature 
#' (MDTX)\cr
#' DAWM = Number of days included in the multiday wind movement
#' (MDWM)\cr
#' DWPR = Number of days with non-zero precipitation included in 
#' multiday precipitation total (MDPR)\cr
#' EVAP = Evaporation of water from evaporation pan (tenths of mm)\cr
#' FMTM = Time of fastest mile or fastest 1-minute wind 
#' (hours and minutes, i.e., HHMM)\cr
#' FRGB = Base of frozen ground layer (cm)\cr
#' FRGT = Top of frozen ground layer (cm)\cr
#' FRTH = Thickness of frozen ground layer (cm)\cr
#' GAHT = Difference between river and gauge height (cm)\cr
#' MDEV = Multiday evaporation total (tenths of mm; use with DAEV)\cr
#' MDPR = Multiday precipitation total (tenths of mm; use with DAPR and 
#'                                      DWPR, if available)\cr
#' MDSF = Multiday snowfall total \cr
#' MDTN = Multiday minimum temperature (tenths of degrees C; use with DATN)\cr
#' MDTX = Multiday maximum temperature (tenths of degress C; use with DATX)\cr
#' MDWM = Multiday wind movement (km)\cr
#' MNPN = Daily minimum temperature of water in an evaporation pan 
#' (tenths of degrees C)\cr
#' MXPN = Daily maximum temperature of water in an evaporation pan 
#' (tenths of degrees C)\cr
#' PGTM = Peak gust time (hours and minutes, i.e., HHMM)\cr
#' PSUN = Daily percent of possible sunshine (percent)\cr
#' SN*# = Minimum soil temperature (tenths of degrees C)
#'   where * corresponds to a code
#' for ground cover and # corresponds to a code for soil 
#' depth.\cr
#' \cr
#' Ground cover codes include the following:\cr
#' 0 = unknown\cr
#' 1 = grass\cr
#' 2 = fallow\cr
#' 3 = bare ground\cr
#' 4 = brome grass\cr
#' 5 = sod\cr
#' 6 = straw multch\cr
#' 7 = grass muck\cr
#' 8 = bare muck\cr
#' \cr
#' Depth codes include the following:\cr
#' 1 = 5 cm\cr
#' 2 = 10 cm\cr
#' 3 = 20 cm\cr
#' 4 = 50 cm\cr
#' 5 = 100 cm\cr
#' 6 = 150 cm\cr
#' 7 = 180 cm\cr
#' \cr
#' SX*# = Maximum soil temperature (tenths of degrees C) 
#'   where * corresponds to a code for ground cover 
#' and # corresponds to a code for soil depth.\cr 
#' See SN*# for ground cover and depth codes. \cr
#' TAVG = Average temperature (tenths of degrees C)
#' [Note that TAVG from source 'S' corresponds
#' to an average for the period ending at
#' 2400 UTC rather than local midnight]\cr
#' THIC = Thickness of ice on water (tenths of mm)	\cr
#' TOBS = Temperature at the time of observation (tenths of degrees C)\cr
#' TSUN = Daily total sunshine (minutes)\cr
#' WDF1 = Direction of fastest 1-minute wind (degrees)\cr
#' WDF2 = Direction of fastest 2-minute wind (degrees)\cr
#' WDF5 = Direction of fastest 5-second wind (degrees)\cr
#' WDFG = Direction of peak wind gust (degrees)\cr
#' WDFI = Direction of highest instantaneous wind (degrees)\cr
#' WDFM = Fastest mile wind direction (degrees)\cr
#' WDMV = 24-hour wind movement (km)	   \cr
#' WESD = Water equivalent of snow on the ground (tenths of mm)\cr
#' WESF = Water equivalent of snowfall (tenths of mm)\cr
#' WSF1 = Fastest 1-minute wind speed (tenths of meters per second)\cr
#' WSF2 = Fastest 2-minute wind speed (tenths of meters per second)\cr
#' WSF5 = Fastest 5-second wind speed (tenths of meters per second)\cr
#' WSFG = Peak gust wind speed (tenths of meters per second)\cr
#' WSFI = Highest instantaneous wind speed (tenths of meters per second)\cr
#' WSFM = Fastest mile wind speed (tenths of meters per second)\cr
#' WT** = Weather Type where ** has one of the following values:\cr
#'   \cr
#' 01 = Fog, ice fog, or freezing fog (may include heavy fog)\cr
#' 02 = Heavy fog or heaving freezing fog (not always 
#'                                         distinquished from fog)\cr
#' 03 = Thunder\cr
#' 04 = Ice pellets, sleet, snow pellets, or small hail \cr
#' 05 = Hail (may include small hail)\cr
#' 06 = Glaze or rime \cr
#' 07 = Dust, volcanic ash, blowing dust, blowing sand, or 
#' blowing obstruction\cr
#' 08 = Smoke or haze \cr
#' 09 = Blowing or drifting snow\cr
#' 10 = Tornado, waterspout, or funnel cloud \cr
#' 11 = High or damaging winds\cr
#' 12 = Blowing spray\cr
#' 13 = Mist\cr
#' 14 = Drizzle\cr
#' 15 = Freezing drizzle \cr
#' 16 = Rain (may include freezing rain, drizzle, and freezing drizzle) \cr
#' 17 = Freezing rain \cr
#' 18 = Snow, snow pellets, snow grains, or ice crystals\cr
#' 19 = Unknown source of precipitation \cr
#' 21 = Ground fog \cr
#' 22 = Ice fog or freezing fog\cr
#' \cr
#' WV** = Weather in the Vicinity where ** has one of the following 
#' values:\cr
#' 01 = Fog, ice fog, or freezing fog (may include heavy fog)\cr
#' 03 = Thunder\cr
#' 07 = Ash, dust, sand, or other blowing obstruction\cr
#' 18 = Snow or ice crystals\cr
#' 20 = Rain or snow shower
#' @param years A numeric vector indicating which years to get.
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing. Defaults to "./RAW/GHCN/".
#' @param extraction.dir A character string indicating where the extracted and cropped GHCN shapefiles should be put.
#' The directory will be created if missing. Defaults to "./EXTRACTIONS/GHCN/".
#' @param standardize Select only common year/month/day? Defaults to FALSE.
#' @param force.redo If an extraction for this template and label already exists, should a new one be created? Defaults to FALSE.
#' @return A named list containing the "spatial" and "tabular" data.
#' @export
#' @importFrom readr read_fwf
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#' # Extract data for the Village Ecodynamics Project "VEPIIN" study area:
#' # http://village.anth.wsu.edu
#' vepPolygon <- polygon_from_extent(raster::extent(672800,740000,4102000,4170000), 
#'      proj4string="+proj=utm +datum=NAD83 +zone=12")
#' 
#' # Get the daily GHCN data (GLOBAL)
#' # Returns a list: the first element is the spatial locations of stations,
#' # and the second is a list of the stations and their daily data
#' GHCN.prcp <- get_ghcn_daily(template=vepPolygon, label="VEPIIN", elements=c('prcp'))
#' 
#' # Plot the VEP polygon
#' plot(vepPolygon)
#' 
#' # Plot the spatial locations
#' plot(GHCN.prcp$spatial, pch=1, add=T)
#' legend('bottomleft', pch=1, legend="GHCN Precipitation Records")
#' 
#' # Elements for which you require the same data
#' # (i.e., minimum and maximum temperature for the same days)
#' # can be standardized using standardize==T
#' GHCN.temp <- get_ghcn_daily(template=vepPolygon, 
#'      label="VEPIIN", 
#'      elements=c('tmin','tmax'), 
#'      standardize=T)
#' 
#' # Plot the VEP polygon
#' plot(vepPolygon)
#' 
#' # Plot the spatial locations
#' plot(GHCN.temp$spatial, pch=1, add=T)
#' legend('bottomleft', pch=1, legend="GHCN Temperature Records")
#' }
get_ghcn_daily <- function(template=NULL, label=NULL, elements=NULL, years = NULL, raw.dir="./RAW/GHCN/", extraction.dir=paste0("./EXTRACTIONS/",label,"/GHCN/"), standardize=F, force.redo=F){
  dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)
  
  if(is.null(template) & is.null(label)){
    label <- "allStations"
  }
  
  if(!is.null(template) & is.null(label)){
    stop("Template provided but no label given.")
  }
  
  message("(Down)Loading GHCN station inventory.")
  if(!force.redo & file.exists(paste0(extraction.dir,"/",label,"_GHCN_stations.shp"))){
    stations.sp <- rgdal::readOGR(dsn = extraction.dir,
                                  layer = paste0(label,"_GHCN_stations"),
                                  verbose = F)
  }else{
    stations.sp <- get_ghcn_inventory(template = template, raw.dir = raw.dir)
    suppressWarnings(rgdal::writeOGR(stations.sp,
                                     dsn = extraction.dir,
                                     layer = paste0(label,"_GHCN_stations"),
                                     driver = "ESRI Shapefile",
                                     overwrite_layer = TRUE)
    )
  }
  
  # If the user didn't specify target elements, get them all.
  if(!is.null(elements)){
    stations.sp <- stations.sp[stations.sp@data[,"ELEMENT"] %in% toupper(elements),]
    missing.elements <- setdiff(toupper(elements),unique(stations.sp$ELEMENT))
    if(length(missing.elements)==length(elements)) stop("No elements available from included stations during given years.")
    if(length(missing.elements)>0) warning("Elements not available: ",paste(missing.elements,collapse = ", "))
  }
  elements <- unique(stations.sp$ELEMENT)
  
  if(standardize){
    stations.sp.splits <- split(as.character(stations.sp$ELEMENT),f=stations.sp$ID, drop=T)
    stations.sp.splits.all <- sapply(stations.sp.splits,function(x){all(sapply(toupper(elements),function(y){y %in% x}))})
    stations.sp <- stations.sp[stations.sp$ID %in% names(stations.sp.splits.all)[stations.sp.splits.all],]
  }
  
  if(!is.null(years)){
    stations.sp <- stations.sp[stations.sp$YEAR_START <= min(years) & stations.sp$YEAR_END >= max(years),]
  }
  
  stations.sp <- stations.sp[!duplicated(stations.sp@data[,c("ID","LATITUDE","LONGITUDE")]),c("ID", "NAME")]
  
  if(!force.redo){
    daily <- tryCatch(lapply(elements,function(element){
      ifelse(!is.null(years),
             paste0(extraction.dir,"/",label,"_GHCN_",element,"_",min(years),"-",max(years),".Rds"),
             paste0(extraction.dir,"/",label,"_GHCN_",element,"_all_years.Rds")) %>%
        readr::read_rds()
    }), warning = function(w){return(NULL)})
    
    if(!is.null(daily)){
      names(daily) <- elements
      
      daily <- lapply(as.character(stations.sp$ID),function(station){
        stationDaily <- tryCatch(lapply(daily,'[[',station),error=function(e){return(NULL)})
        stationDaily <- stationDaily[!sapply(stationDaily,is.null)]
        return(stationDaily)
      })
      names(daily) <- as.character(stations.sp$ID)
      daily <- daily[!sapply(daily,is.null)]
      
      return(list(spatial=stations.sp,tabular=daily))
      
    }
  }
  
  daily <- lapply(stations.sp$ID,function(station){
    message("(Down)Loading GHCN station data for station ",as.character(station))
    tryCatch(station <- get_ghcn_daily_station(ID=station, elements=elements, years = years, raw.dir=raw.dir, standardize=standardize, force.redo=force.redo), error = function(e){message("Error (down)Loading GHCN station data for station ",as.character(station)); return(NULL)})
    return(station)
  })
  names(daily) <- stations.sp$ID
  
  daily <- daily[(lapply(daily,names) %>% sapply(length)) <= length(elements) & (lapply(daily,names) %>% sapply(length)) > 0]
  
  daily.split <- lapply(as.character(elements),function(element){
    lapply(daily,'[[',element)
  })
  names(daily.split) <- elements
  junk <- lapply(as.character(elements),function(element){
    
    saveRDS(daily.split[[element]],ifelse(!is.null(years),
                                          paste0(extraction.dir,"/",label,"_GHCN_",element,"_",min(years),"-",max(years),".Rds"),
                                          paste0(extraction.dir,"/",label,"_GHCN_",element,"_all_years.Rds")),compress='xz')
  })
  
  
  return(list(spatial=stations.sp,tabular=daily))
}

#' Download the daily data for a GHCN weather station.
#'
#' @param ID A character string giving the station ID.
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' @param force.redo If this weather station has been downloaded before, should it be updated? Defaults to FALSE.
#' @return A character string representing the full local path of the GHCN station data.
#' @export
#' @keywords internal
download_ghcn_daily_station <- function(ID, raw.dir, force.redo=F){
  
  dir.create(raw.dir, recursive=T, showWarnings=F)
  
  url <- paste("ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/all/",ID,".dly",sep='')
  if(!force.redo){
    download_data(url=url, destdir=raw.dir, timestamping=T)
  }else{
    download_data(url=url, destdir=raw.dir, timestamping=F)
  }
  
  return(normalizePath(paste(raw.dir,ID,".dly",sep='')))
  
}

#' Download and extract the daily data for a GHCN weather station.
#'
#' \code{get_ghcn_daily_station} returns a named list of \code{\link{data.frame}s}, one for
#' each \code{elements}. If \code{elements} is undefined, it returns all available weather
#' tables for the station
#' 
#' @param ID A character string giving the station ID.
#' @param elements A character vector of elemets to extract.\cr
#' The five core elements are:\cr
#' PRCP = Precipitation (tenths of mm)\cr
#' SNOW = Snowfall (mm)\cr
#' SNWD = Snow depth (mm)\cr
#' TMAX = Maximum temperature (tenths of degrees C)\cr
#' TMIN = Minimum temperature (tenths of degrees C)\cr
#' \cr
#' The other elements are:\cr
#'   
#' ACMC = Average cloudiness midnight to midnight from 30-second 
#' ceilometer data (percent)\cr
#' ACMH = Average cloudiness midnight to midnight from 
#' manual observations (percent)\cr
#' ACSC = Average cloudiness sunrise to sunset from 30-second 
#' ceilometer data (percent)\cr
#' ACSH = Average cloudiness sunrise to sunset from manual 
#' observations (percent)\cr
#' AWDR = Average daily wind direction (degrees)\cr
#' AWND = Average daily wind speed (tenths of meters per second)\cr
#' DAEV = Number of days included in the multiday evaporation
#' total (MDEV)\cr
#' DAPR = Number of days included in the multiday precipiation 
#' total (MDPR)\cr
#' DASF = Number of days included in the multiday snowfall 
#' total (MDSF)\cr
#' DATN = Number of days included in the multiday minimum temperature 
#' (MDTN)\cr
#' DATX = Number of days included in the multiday maximum temperature 
#' (MDTX)\cr
#' DAWM = Number of days included in the multiday wind movement
#' (MDWM)\cr
#' DWPR = Number of days with non-zero precipitation included in 
#' multiday precipitation total (MDPR)\cr
#' EVAP = Evaporation of water from evaporation pan (tenths of mm)\cr
#' FMTM = Time of fastest mile or fastest 1-minute wind 
#' (hours and minutes, i.e., HHMM)\cr
#' FRGB = Base of frozen ground layer (cm)\cr
#' FRGT = Top of frozen ground layer (cm)\cr
#' FRTH = Thickness of frozen ground layer (cm)\cr
#' GAHT = Difference between river and gauge height (cm)\cr
#' MDEV = Multiday evaporation total (tenths of mm; use with DAEV)\cr
#' MDPR = Multiday precipitation total (tenths of mm; use with DAPR and 
#'                                      DWPR, if available)\cr
#' MDSF = Multiday snowfall total \cr
#' MDTN = Multiday minimum temperature (tenths of degrees C; use with DATN)\cr
#' MDTX = Multiday maximum temperature (tenths of degress C; use with DATX)\cr
#' MDWM = Multiday wind movement (km)\cr
#' MNPN = Daily minimum temperature of water in an evaporation pan 
#' (tenths of degrees C)\cr
#' MXPN = Daily maximum temperature of water in an evaporation pan 
#' (tenths of degrees C)\cr
#' PGTM = Peak gust time (hours and minutes, i.e., HHMM)\cr
#' PSUN = Daily percent of possible sunshine (percent)\cr
#' SN*# = Minimum soil temperature (tenths of degrees C)
#'   where * corresponds to a code
#' for ground cover and # corresponds to a code for soil 
#' depth.\cr
#' \cr
#' Ground cover codes include the following:\cr
#' 0 = unknown\cr
#' 1 = grass\cr
#' 2 = fallow\cr
#' 3 = bare ground\cr
#' 4 = brome grass\cr
#' 5 = sod\cr
#' 6 = straw multch\cr
#' 7 = grass muck\cr
#' 8 = bare muck\cr
#' \cr
#' Depth codes include the following:\cr
#' 1 = 5 cm\cr
#' 2 = 10 cm\cr
#' 3 = 20 cm\cr
#' 4 = 50 cm\cr
#' 5 = 100 cm\cr
#' 6 = 150 cm\cr
#' 7 = 180 cm\cr
#' \cr
#' SX*# = Maximum soil temperature (tenths of degrees C) 
#'   where * corresponds to a code for ground cover 
#' and # corresponds to a code for soil depth.\cr 
#' See SN*# for ground cover and depth codes. \cr
#' TAVG = Average temperature (tenths of degrees C)
#' [Note that TAVG from source 'S' corresponds
#' to an average for the period ending at
#' 2400 UTC rather than local midnight]\cr
#' THIC = Thickness of ice on water (tenths of mm)	\cr
#' TOBS = Temperature at the time of observation (tenths of degrees C)\cr
#' TSUN = Daily total sunshine (minutes)\cr
#' WDF1 = Direction of fastest 1-minute wind (degrees)\cr
#' WDF2 = Direction of fastest 2-minute wind (degrees)\cr
#' WDF5 = Direction of fastest 5-second wind (degrees)\cr
#' WDFG = Direction of peak wind gust (degrees)\cr
#' WDFI = Direction of highest instantaneous wind (degrees)\cr
#' WDFM = Fastest mile wind direction (degrees)\cr
#' WDMV = 24-hour wind movement (km)	   \cr
#' WESD = Water equivalent of snow on the ground (tenths of mm)\cr
#' WESF = Water equivalent of snowfall (tenths of mm)\cr
#' WSF1 = Fastest 1-minute wind speed (tenths of meters per second)\cr
#' WSF2 = Fastest 2-minute wind speed (tenths of meters per second)\cr
#' WSF5 = Fastest 5-second wind speed (tenths of meters per second)\cr
#' WSFG = Peak gust wind speed (tenths of meters per second)\cr
#' WSFI = Highest instantaneous wind speed (tenths of meters per second)\cr
#' WSFM = Fastest mile wind speed (tenths of meters per second)\cr
#' WT** = Weather Type where ** has one of the following values:\cr
#'   \cr
#' 01 = Fog, ice fog, or freezing fog (may include heavy fog)\cr
#' 02 = Heavy fog or heaving freezing fog (not always 
#'                                         distinquished from fog)\cr
#' 03 = Thunder\cr
#' 04 = Ice pellets, sleet, snow pellets, or small hail \cr
#' 05 = Hail (may include small hail)\cr
#' 06 = Glaze or rime \cr
#' 07 = Dust, volcanic ash, blowing dust, blowing sand, or 
#' blowing obstruction\cr
#' 08 = Smoke or haze \cr
#' 09 = Blowing or drifting snow\cr
#' 10 = Tornado, waterspout, or funnel cloud \cr
#' 11 = High or damaging winds\cr
#' 12 = Blowing spray\cr
#' 13 = Mist\cr
#' 14 = Drizzle\cr
#' 15 = Freezing drizzle \cr
#' 16 = Rain (may include freezing rain, drizzle, and freezing drizzle) \cr
#' 17 = Freezing rain \cr
#' 18 = Snow, snow pellets, snow grains, or ice crystals\cr
#' 19 = Unknown source of precipitation \cr
#' 21 = Ground fog \cr
#' 22 = Ice fog or freezing fog\cr
#' \cr
#' WV** = Weather in the Vicinity where ** has one of the following 
#' values:\cr
#' 01 = Fog, ice fog, or freezing fog (may include heavy fog)\cr
#' 03 = Thunder\cr
#' 07 = Ash, dust, sand, or other blowing obstruction\cr
#' 18 = Snow or ice crystals\cr
#' 20 = Rain or snow shower
#' @param years A numeric vector indicating which years to get.
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' @param standardize Select only common year/month/day? Defaults to FALSE.
#' @param force.redo If this weather station has been downloaded before, should it be updated? Defaults to FALSE.
#' @return A named list of \code{\link{data.frame}s}, one for each \code{elements}.
#' @export
#' @keywords internal
get_ghcn_daily_station <- function(ID, elements=NULL, years = NULL, raw.dir, standardize=F, force.redo=F){
  
  file <- download_ghcn_daily_station(ID = ID, 
                                      raw.dir = paste(raw.dir,"/DAILY/",sep=''), 
                                      force.redo = force.redo)
  
  # GHCN files are fixed-width. The numbers here refer to those column widths.
  daily <- readr::read_fwf(file,
                           col_positions = readr::fwf_positions(start = c(1,12,16,18,seq(22,262,8)),
                                                                end = c(11,15,17,21,seq(26,266,8)),
                                                                col_names = c("STATION","YEAR","MONTH","ELEMENT",paste0("D",1:31))
                           ),
                           col_types = paste0(c("cicc",rep("i",31)),collapse = "")
  )
  
  # If the user didn't specify target elements, get them all.
  if(is.null(elements)){
    elements <- unique(daily$ELEMENT)
  }
  elements <- toupper(elements)
  daily %<>% dplyr::filter_(~ELEMENT %in% elements)
  missing.elements <- setdiff(toupper(elements),unique(daily$ELEMENT))
  if(length(missing.elements)>0) warning("Elements not available: ",paste(missing.elements,collapse = ", "))
  
  # If the user didn't specify target elements, get them all.
  if(!is.null(years)){
    daily %<>% dplyr::filter_(~YEAR %in% years)
  }
  
  daily %<>% dplyr::mutate_all(dplyr::funs_(quote(ifelse(. == -9999,NA,.))))
  
  ## Separate by element
  out.list <- lapply(elements, function(element){
    return(
      daily %>%
        dplyr::filter_(~ELEMENT == element) %>%
        dplyr::select_(quote(-ELEMENT))
    )
  })
  
  ## If standardize, select only common year/month/day, and make NA if both not present
  if(standardize){
    yearMonths <- lapply(out.list, function(element){
      element %<>%
        dplyr::arrange_(~YEAR,~MONTH)
      return(paste("Y",element[["YEAR"]],"M",element[["MONTH"]],sep=''))
    })
    
    all.yearMonths <- Reduce(intersect,yearMonths)
    
    out.list <- lapply(out.list, function(element){
      element.yearMonths <- paste("Y",element[["YEAR"]],"M",element[["MONTH"]],sep='')
      return(
        element %>%
          dplyr::filter(element.yearMonths %in% all.yearMonths)
      )
    })
  }
  
  names(out.list) <- elements
  
  return(out.list)
}

#' Download and crop the inventory of GHCN stations.
#'
#' \code{get_ghcn_inventory} returns a \code{SpatialPolygonsDataFrame} of the GHCN stations within
#' the specified \code{template}. If template is not provided, returns the entire GHCN inventory.
#' 
#' Stations with multiple elements will have multiple points. This allows for easy mapping of stations
#' by element availability.
#' 
#' @param template A Raster* or Spatial* object to serve 
#' as a template for cropping.
#' @param elements A character vector of elemets to extract.
#' Common elements include "tmin", "tmax", and "prcp".
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing.
#' @return A \code{SpatialPolygonsDataFrame} of the GHCN stations within
#' the specified \code{template}
#' @export
#' @keywords internal
get_ghcn_inventory <- function(template=NULL, elements=NULL, raw.dir){
  if(!is.null(template) & !(class(template) %in% c("SpatialPolygonsDataFrame","SpatialPolygons","character"))){
    template <- polygon_from_extent(template)
  }
  
  url <- "ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/ghcnd-inventory.txt"
  destdir <- raw.dir
  download_data(url=url, destdir=destdir)
  
  url <- "ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt"
  destdir <- raw.dir
  download_data(url=url, destdir=destdir)
  
  # GHCN files are fixed-width. The numbers here refer to those column widths.
  station.inventory <- readr::read_fwf(paste(raw.dir,"ghcnd-inventory.txt",sep=''),
                                       readr::fwf_positions(start = c(1,13,22,32,37,42),
                                                            end = c(11,20,30,35,40,45),
                                                            col_names = c("ID","LATITUDE","LONGITUDE","ELEMENT","YEAR_START","YEAR_END")),
                                       col_types = "cddcii")
  
  stations <- readr::read_fwf(paste(raw.dir,"ghcnd-stations.txt",sep=''),
                              readr::fwf_positions(start = c(1,13,22,32,42),
                                                   end = c(11,20,30,38,72),
                                                   col_names = c("ID","LATITUDE","LONGITUDE","ELEVATION","NAME")),
                              col_types = "cdddc")
  
  # Convert to SPDF
  stations.sp <- sp::SpatialPointsDataFrame(coords = station.inventory %>% 
                                              dplyr::select_(~LONGITUDE,~LATITUDE),
                                            data = station.inventory %>% 
                                              dplyr::left_join(stations %>% dplyr::select_("ID","NAME"), by="ID") %>%
                                              dplyr::select_("ID","NAME","LATITUDE","LONGITUDE","ELEMENT","YEAR_START","YEAR_END") %>%
                                              as.data.frame(),
                                            proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
  
  if(!is.null(elements)){
    stations.sp <- stations.sp[stations.sp$ELEMENT %in% toupper(elements),]
  }
  
  if(!is.null(template)){
    if(class(template) == "character"){
      missing.stations <- setdiff(template,unique(stations.sp$ID))
      if(length(missing.stations)>0) warning("Stations not available: ",paste(missing.stations,collapse = ", "))
      stations.sp <- stations.sp[stations.sp$ID %in% template,]
    }else{
      template <- methods::as(template,"SpatialPolygons")
      stations.sp <- stations.sp[!is.na(sp::over(stations.sp,sp::spTransform(template,sp::CRS(raster::projection(stations.sp))))),]
    }
  }
  
  return(stations.sp)
}

#' Convert a list of station data to a single data frame.
#'
#' \code{station_to_data_frame} returns a \code{data.frame} of the GHCN station data list.
#' 
#' This function unwraps the station data and merges all data into a single data frame, 
#' with the first column being in the \code{Date} class.
#' 
#' @param station.data A named list containing station data
#' @return A \code{data.frame} of the containing the unwrapped station data
#' @export
#' @keywords internal
station_to_data_frame <- function(station.data){
  data.list <- lapply(1:length(station.data),function(i){
    X <- station.data[[i]]
    
    # Get just the climate info
    annual.records <- as.matrix(X[,4:34])
    
    # Get the number of days per month in the records
    n.days <- Hmisc::monthDays(as.Date(paste(X$YEAR,X$MONTH,'01',sep='-')))
    
    ## Unnwrap each row, accounting for number of days in the month
    annual.records.unwrapped <- unwrap_rows(annual.records,n.days)
    
    dates <- as.Date(unlist(mapply(FUN = function(days,dates){paste(dates,days,sep="-")}, sapply(n.days,function(x){1:x}), paste(X$YEAR,X$MONTH,sep='-'))))
    
    annual.records.unwrapped <- data.table::data.table(dates,annual.records.unwrapped)
    names(annual.records.unwrapped) <- c("DATE",names(station.data)[i])
    data.table::setkey(annual.records.unwrapped,"DATE")
    
    return(annual.records.unwrapped)
  })
  
  return(Reduce(function(x,y) merge(x,y,all=TRUE),data.list))
}
