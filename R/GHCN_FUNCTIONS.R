# globalVariables(c('.', 'ELEMENT','YEAR','MONTH','LONGITUDE','LATITUDE'))

#' Download and crop the Global Historical Climate Network-Daily data.
#'
#' `get_ghcn_daily` returns a named list of length 2:
#' \enumerate{
#' \item 'spatial': A [`Simple Feature`][sf::sf] of the locations of GHCN weather stations
#' in the template, and
#' \item 'tabular': A named list of type [data.frame()] with the daily weather data for each station.
#' The name of each list item is the station ID.
#' }
#'
#' @param template An [`Simple Feature`][sf::sf]
#' or [`SpatRaster`][terra::SpatRaster] object to serve as a template for cropping.
#' Alternatively, a character vector providing GHCN station IDs. If missing, all stations
#' will be downloaded!
#' @param label A character string naming the study area.
#' @param elements A character vector of elements to extract.\cr
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
#' DAPR = Number of days included in the multiday precipitation
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
#' MDTX = Multiday maximum temperature (tenths of degrees C; use with DATX)\cr
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
#' (Note that TAVG from source 'S' corresponds
#' to an average for the period ending at
#' 2400 UTC rather than local midnight)\cr
#' THIC = Thickness of ice on water (tenths of mm)\cr
#' TOBS = Temperature at the time of observation (tenths of degrees C)\cr
#' TSUN = Daily total sunshine (minutes)\cr
#' WDF1 = Direction of fastest 1-minute wind (degrees)\cr
#' WDF2 = Direction of fastest 2-minute wind (degrees)\cr
#' WDF5 = Direction of fastest 5-second wind (degrees)\cr
#' WDFG = Direction of peak wind gust (degrees)\cr
#' WDFI = Direction of highest instantaneous wind (degrees)\cr
#' WDFM = Fastest mile wind direction (degrees)\cr
#' WDMV = 24-hour wind movement (km)\cr
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
#'                                         distinguished from fog)\cr
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
#' The directory will be created if missing. Defaults to './RAW/GHCN/'.
#' @param extraction.dir A character string indicating where the extracted and cropped GHCN shapefiles should be put.
#' The directory will be created if missing. Defaults to './EXTRACTIONS/GHCN/'.
#' @param standardize Select only common year/month/day? Defaults to FALSE.
#' @param force.redo If an extraction for this template and label already exists, should a new one be created? Defaults to FALSE.
#' @return A named list containing the 'spatial' and 'tabular' data.
#' @export
#' @importFrom readr read_fwf
#' @importFrom magrittr %>%
#' @examples
#' \donttest{
#' # Get the daily GHCN data (GLOBAL)
#' # Returns a list: the first element is the spatial locations of stations,
#' # and the second is a list of the stations and their daily data
#' GHCN.prcp <-
#'   get_ghcn_daily(
#'     template = FedData::meve,
#'     label = "meve",
#'     elements = c("prcp")
#'   )
#'
#' # Plot the VEP polygon
#' plot(meve)
#'
#' # Plot the spatial locations
#' plot(GHCN.prcp$spatial$geometry, pch = 1, add = TRUE)
#' legend("bottomleft", pch = 1, legend = "GHCN Precipitation Records")
#'
#' # Elements for which you require the same data
#' # (i.e., minimum and maximum temperature for the same days)
#' # can be standardized using `standardize = TRUE`
#' GHCN.temp <- get_ghcn_daily(
#'   template = FedData::meve,
#'   label = "meve",
#'   elements = c("tmin", "tmax"),
#'   standardize = TRUE
#' )
#'
#' # Plot the VEP polygon
#' plot(meve)
#'
#' # Plot the spatial locations
#' plot(GHCN.temp$spatial$geometry, pch = 1, add = TRUE)
#' legend("bottomleft", pch = 1, legend = "GHCN Temperature Records")
#' }
get_ghcn_daily <- function(template = NULL,
                           label = NULL,
                           elements = NULL,
                           years = NULL,
                           raw.dir =
                             file.path(
                               tempdir(),
                               "FedData",
                               "raw",
                               "ghcn"
                             ),
                           extraction.dir =
                             file.path(
                               tempdir(),
                               "FedData",
                               "extractions",
                               "ned",
                               label
                             ),
                           standardize = FALSE,
                           force.redo = FALSE) {
  # raw.dir <- normalizePath(paste0(raw.dir, "/."), mustWork = FALSE)
  # extraction.dir <- normalizePath(paste0(extraction.dir, "/."), mustWork = FALSE)

  dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)

  if (is.null(template) & is.null(label)) {
    label <- "allStations"
  }

  if (!is.null(template) & is.null(label)) {
    stop("Template provided but no label given.")
  }

  message("(Down)Loading GHCN station inventory.")
  raw_inventory <-
    file.path(
      extraction.dir,
      paste0(label, "_GHCN_stations.fgb")
    )
  if (force.redo | !file.exists(raw_inventory)) {
    unlink(raw_inventory, force = TRUE)

    stations.sf <-
      get_ghcn_inventory(template = template, raw.dir = raw.dir)

    stations.sf %>%
      sf::write_sf(
        dsn = raw_inventory
      )
  }

  stations.sf <-
    sf::read_sf(raw_inventory)

  # If the user didn't specify target elements, get them all.
  if (!is.null(elements)) {
    stations.sf <- stations.sf[stations.sf$ELEMENT %in% toupper(elements), ]
    missing.elements <- setdiff(toupper(elements), unique(stations.sf$ELEMENT))
    if (length(missing.elements) == length(elements)) {
      stop("No elements available from included stations during given years.")
    }
    if (length(missing.elements) > 0) {
      warning("Elements not available: ", paste(missing.elements, collapse = ", "))
    }
  }
  elements <- unique(stations.sf$ELEMENT)

  if (standardize) {
    stations.sf.splits <- split(as.character(stations.sf$ELEMENT),
      f = stations.sf$ID,
      drop = TRUE
    )
    stations.sf.splits.all <- sapply(stations.sf.splits, function(x) {
      all(sapply(toupper(elements), function(y) {
        y %in% x
      }))
    })
    stations.sf <- stations.sf[stations.sf$ID %in% names(stations.sf.splits.all)[stations.sf.splits.all], ]
  }

  if (!is.null(years)) {
    stations.sf <- stations.sf[stations.sf$YEAR_START <= min(years) & stations.sf$YEAR_END >= max(years), ]
  }

  stations.sf <- stations.sf[!duplicated(stations.sf[, c("ID", "LATITUDE", "LONGITUDE")]), c("ID", "NAME")]

  if (!force.redo) {
    daily <- tryCatch(lapply(elements, function(element) {
      ifelse(!is.null(years), paste0(extraction.dir, "/", label, "_GHCN_", element, "_", min(years), "-", max(years), ".Rds"),
        paste0(extraction.dir, "/", label, "_GHCN_", element, "_all_years.Rds")
      ) %>% readr::read_rds()
    }), warning = function(w) {
      return(NULL)
    })

    if (!is.null(daily)) {
      names(daily) <- elements

      daily <- lapply(as.character(stations.sf$ID), function(station) {
        stationDaily <- tryCatch(lapply(daily, "[[", station), error = function(e) {
          return(NULL)
        })
        stationDaily <- stationDaily[!sapply(stationDaily, is.null)]
        return(stationDaily)
      })
      names(daily) <- as.character(stations.sf$ID)
      daily <- daily[!sapply(daily, is.null)]

      return(list(spatial = stations.sf, tabular = daily))
    }
  }

  daily <- lapply(stations.sf$ID, function(station) {
    message("(Down)Loading GHCN station data for station ", as.character(station))
    tryCatch(station <- get_ghcn_daily_station(
      ID = station, elements = elements, years = years, raw.dir = raw.dir, standardize = standardize,
      force.redo = force.redo
    ), error = function(e) {
      message("Error (down)Loading GHCN station data for station ", as.character(station))
      return(NULL)
    })
    return(station)
  })
  names(daily) <- stations.sf$ID

  daily <- daily[(lapply(daily, names) %>% sapply(length)) <= length(elements) & (lapply(daily, names) %>% sapply(length)) >
    0]

  daily.split <- lapply(as.character(elements), function(element) {
    lapply(daily, "[[", element)
  })
  names(daily.split) <- elements
  junk <- lapply(as.character(elements), function(element) {
    saveRDS(daily.split[[element]], ifelse(!is.null(years), paste0(
      extraction.dir, "/", label, "_GHCN_", element, "_", min(years),
      "-", max(years), ".Rds"
    ), paste0(extraction.dir, "/", label, "_GHCN_", element, "_all_years.Rds")), compress = "xz")
  })


  return(list(spatial = stations.sf, tabular = daily))
}

#' Download the daily data for a GHCN weather station.
#'
#' @param ID A character string giving the station ID.
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' @param force.redo If this weather station has been downloaded before, should it be updated? Defaults to FALSE.
#' @return A character string representing the full local path of the GHCN station data.
#' @export
#' @keywords internal
download_ghcn_daily_station <- function(ID, raw.dir, force.redo = FALSE) {
  dir.create(raw.dir, recursive = TRUE, showWarnings = FALSE)

  url <- paste("https://www.ncei.noaa.gov/pub/data/ghcn/daily/all/", ID, ".dly", sep = "")
  if (!force.redo) {
    download_data(url = url, destdir = raw.dir, timestamping = TRUE)
  } else {
    download_data(url = url, destdir = raw.dir, timestamping = FALSE)
  }

  return(normalizePath(paste0(raw.dir, "/", ID, ".dly")))
}

#' Download and extract the daily data for a GHCN weather station.
#'
#' \code{get_ghcn_daily_station} returns a named list of \code{\link{data.frame}s}, one for
#' each \code{elements}. If \code{elements} is undefined, it returns all available weather
#' tables for the station
#'
#' @param ID A character string giving the station ID.
#' @param elements A character vector of elements to extract.\cr
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
#' DAPR = Number of days included in the multiday precipitation
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
#' MDTX = Multiday maximum temperature (tenths of degrees C; use with DATX)\cr
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
#' (Note that TAVG from source 'S' corresponds
#' to an average for the period ending at
#' 2400 UTC rather than local midnight)\cr
#' THIC = Thickness of ice on water (tenths of mm)\cr
#' TOBS = Temperature at the time of observation (tenths of degrees C)\cr
#' TSUN = Daily total sunshine (minutes)\cr
#' WDF1 = Direction of fastest 1-minute wind (degrees)\cr
#' WDF2 = Direction of fastest 2-minute wind (degrees)\cr
#' WDF5 = Direction of fastest 5-second wind (degrees)\cr
#' WDFG = Direction of peak wind gust (degrees)\cr
#' WDFI = Direction of highest instantaneous wind (degrees)\cr
#' WDFM = Fastest mile wind direction (degrees)\cr
#' WDMV = 24-hour wind movement (km)\cr
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
#'                                         distinguished from fog)\cr
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
get_ghcn_daily_station <- function(ID,
                                   elements = NULL,
                                   years = NULL,
                                   raw.dir,
                                   standardize = FALSE,
                                   force.redo = FALSE) {
  file <- download_ghcn_daily_station(ID = ID, raw.dir = paste(raw.dir, "/DAILY/", sep = ""), force.redo = force.redo)

  # GHCN files are fixed-width. The numbers here refer to those column widths.
  daily <- readr::read_fwf(file, col_positions = readr::fwf_positions(start = c(1, 12, 16, 18, seq(22, 262, 8)), end = c(
    11,
    15, 17, 21, seq(26, 266, 8)
  ), col_names = c("STATION", "YEAR", "MONTH", "ELEMENT", paste0("D", 1:31))), col_types = paste0(c(
    "cicc",
    rep("i", 31)
  ), collapse = ""))

  # If the user didn't specify target elements, get them all.
  if (is.null(elements)) {
    elements <- unique(daily$ELEMENT)
  }
  elements <- toupper(elements)
  daily %<>% dplyr::filter(ELEMENT %in% elements)
  missing.elements <- setdiff(toupper(elements), unique(daily$ELEMENT))
  if (length(missing.elements) > 0) {
    warning("Elements not available: ", paste(missing.elements, collapse = ", "))
  }

  # If the user didn't specify target elements, get them all.
  if (!is.null(years)) {
    daily %<>% dplyr::filter(YEAR %in% years)
  }

  daily %<>%
    dplyr::mutate(
      dplyr::across(
        dplyr::where(is.numeric),
        ~ dplyr::na_if(.x, -9999)
      )
    )

  ## Separate by element
  out.list <- lapply(elements, function(element) {
    return(daily %>% dplyr::filter(ELEMENT == element) %>% dplyr::select(-ELEMENT))
  })

  ## If standardize, select only common year/month/day, and make NA if both not present
  if (standardize) {
    yearMonths <- lapply(out.list, function(element) {
      element %<>% dplyr::arrange(YEAR, MONTH)
      return(paste("Y", element[["YEAR"]], "M", element[["MONTH"]], sep = ""))
    })

    all.yearMonths <- Reduce(intersect, yearMonths)

    out.list <- lapply(out.list, function(element) {
      element.yearMonths <- paste("Y", element[["YEAR"]], "M", element[["MONTH"]], sep = "")
      return(element %>% dplyr::filter(element.yearMonths %in% all.yearMonths))
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
#' @param template An [`Simple Feature`][sf::sf]
#' or [`SpatRaster`][terra::SpatRaster] object to serve as a template for cropping.
#' @param elements A character vector of elements to extract.
#' Common elements include 'tmin', 'tmax', and 'prcp'.
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing.
#' @return A [`Simple Feature`][sf::sf] of the GHCN stations within
#' the specified `template`
#' @export
#' @keywords internal
get_ghcn_inventory <- function(template = NULL, elements = NULL, raw.dir) {
  if (!is.null(template) & !inherits(template, c("SpatialPolygonsDataFrame", "SpatialPolygons", "character"))) {
    template %<>%
      template_to_sf() %>%
      sf::st_transform(4326) %>%
      polygon_from_extent()
  }

  # GHCN files are fixed-width. The numbers here refer to those column widths.
  station.inventory <-
    readr::read_fwf("https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-inventory.txt",
      readr::fwf_positions(
        start = c(1, 13, 22, 32, 37, 42),
        end = c(11, 20, 30, 35, 40, 45),
        col_names = c(
          "ID",
          "LATITUDE",
          "LONGITUDE",
          "ELEMENT",
          "YEAR_START",
          "YEAR_END"
        )
      ),
      col_types = "cddcii"
    )

  stations <-
    readr::read_fwf("https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt",
      readr::fwf_positions(
        start = c(1, 13, 22, 32, 42),
        end = c(11, 20, 30, 38, 72),
        col_names = c(
          "ID",
          "LATITUDE",
          "LONGITUDE",
          "ELEVATION",
          "NAME"
        )
      ),
      col_types = "cdddc"
    ) %>%
    sf::st_as_sf(
      coords = c("LONGITUDE", "LATITUDE"),
      crs = 4326
    )

  if (!is.null(template)) {
    if (inherits(template, "character")) {
      missing.stations <- setdiff(template, unique(stations$ID))
      if (length(missing.stations) > 0) {
        warning("Stations not available: ", paste(missing.stations, collapse = ", "))
      }

      stations %<>%
        dplyr::filter(ID %in% template)
    } else {
      suppressWarnings(
        stations %<>%
          sf::st_intersection(template)
      )
    }
  }

  # Convert to SF
  stations.sf <-
    stations %>%
    dplyr::select(ID, NAME) %>%
    dplyr::left_join(station.inventory,
      by = "ID"
    )

  if (!is.null(elements)) {
    stations.sf %<>%
      dplyr::filter(ELEMENT %in% toupper(elements))
  }

  return(stations.sf)
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
station_to_data_frame <- function(station.data) {
  station.data %>%
    dplyr::bind_rows(.id = "ELEMENT") %>%
    dplyr::mutate(ELEMENT = factor(ELEMENT,
      levels = names(station.data),
      ordered = TRUE
    )) %>%
    tidyr::pivot_longer(
      !c(
        ELEMENT,
        STATION,
        YEAR,
        MONTH
      ),
      names_to = "DAY"
    ) %>%
    stats::na.omit() %>%
    dplyr::mutate(DATE = lubridate::as_date(paste0(YEAR, MONTH, stringr::str_remove(DAY, "D")))) %>%
    dplyr::select(STATION, DATE, ELEMENT, value) %>%
    tidyr::pivot_wider(
      names_from = "ELEMENT",
      values_from = "value"
    ) %>%
    dplyr::select(!STATION)
}
