#' Download and crop the 1-km DAYMET daily weather dataset.
#'
#' \code{get_daymet} returns a \code{RasterBrick} of weather data cropped to a given
#' template study area.
#'
#' @param template A Raster* or Spatial* object to serve
#' as a template for cropping.
#' @param label A character string naming the study area.
#' @param elements A character vector of elements to extract.\cr
#' The available elements are:\cr
#' dayl = Duration of the daylight period in seconds per day. This calculation is based on the period of the day during which the sun is above a hypothetical flat horizon.\cr
#' prcp = Daily total precipitation in millimeters per day, sum of all forms converted to water-equivalent. Precipitation occurrence on any given day may be ascertained.\cr
#' srad = Incident shortwave radiation flux density in watts per square meter, taken as an average over the daylight period of the day. NOTE: Daily total radiation (MJ/m2/day) can be calculated as follows: ((srad (W/m2) * dayl (s/day)) / l,000,000)\cr
#' swe = Snow water equivalent in kilograms per square meter. The amount of water contained within the snowpack.\cr
#' tmax = Daily maximum 2-meter air temperature in degrees Celsius.\cr
#' tmin = Daily minimum 2-meter air temperature in degrees Celsius.\cr
#' vp = Water vapor pressure in pascals. Daily average partial pressure of water vapor.\cr
#' @param years A numeric vector of years to extract.
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing. Defaults to './RAW/DAYMET/'.
#' @param extraction.dir A character string indicating where the extracted and cropped DEM should be put.
#' The directory will be created if missing. Defaults to './EXTRACTIONS/DAYMET/'.
#' @param force.redo If an extraction for this template and label already exists, should a new one be created?
#' @return A named list of \code{RasterBrick}s of weather data cropped to the extent of the template.
#' @importFrom lubridate year
#' @importFrom magrittr %>% %$%
#' @importFrom readr read_rds write_rds
#' @importFrom raster brick projection mosaic writeRaster
#' @importFrom sp spTransform %over%
#' @importFrom foreach foreach %do%
#' @importFrom utils data
#' @export
#' @examples
#' \dontrun{
#' # Extract data for the Village Ecodynamics Project 'VEPIIN' study area:
#' # http://village.anth.wsu.edu
#' vepPolygon <- polygon_from_extent(raster::extent(672800,740000,4102000,4170000),
#'      proj4string='+proj=utm +datum=NAD83 +zone=12')
#'
#' # Get the DAYMET (North America only)
#' # Returns a list of raster bricks
#' DAYMET <- get_daymet(template=vepPolygon,
#'                      label='VEPIIN',
#'                      elements = c('prcp','tmin','tmax'),
#'                      years = 1980:1985)
#'
#' # Plot with raster::plot
#' plot(DAYMET$tmin$X1985.10.23)
#' }
get_daymet <- function(template,
                       label,
                       elements = NULL,
                       years = NULL,
                       raw.dir = "./RAW/DAYMET",
                       extraction.dir = paste0("./EXTRACTIONS/", label, "/DAYMET"),
                       force.redo = F) {
  
  raw.dir <- normalizePath(paste0(raw.dir,"/."), mustWork = FALSE)  
  extraction.dir <- normalizePath(paste0(extraction.dir,"/."), mustWork = FALSE) 
  
  dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)
  
  all.elements <- c("dayl", "prcp", "srad", "swe", "tmax", "tmin", "vp")
  
  if (is.null(elements))
    elements <- all.elements
  
  elements <- tolower(elements)
  
  missing.elements <- setdiff(elements, all.elements)
  if (length(missing.elements) > 0)
    warning("Elements not available: ", paste(missing.elements, collapse = ", "))
  elements <- setdiff(elements, missing.elements)
  if (length(elements) == 0)
    stop("No elements available")
  
  all.years <- 1980:(lubridate::year(Sys.time()) - 1)
  if (is.null(years))
    years <- all.years
  
  missing.years <- setdiff(years, all.years)
  if (length(missing.years) > 0)
    warning("Years not available: ", paste(missing.years, collapse = ", "))
  years <- setdiff(years, missing.years)
  if (length(years) == 0)
    stop("No years available")
  
  
  out.files <- paste0(extraction.dir, "/", label, "_DAYMET_", elements, "_", min(years), "-", max(years), ".tif")
  if (!force.redo & all(file.exists(out.files)) & file.exists(paste0(extraction.dir, "/", label, "_DAYMET_layer_names.Rds"))) {
    extracted.DAYMET <- out.files %>% lapply(FUN = function(x) {
      out <- raster::brick(x)
      names(out) <- readr::read_rds(paste0(extraction.dir, "/", label, "_DAYMET_layer_names.Rds"))
      return(out)
    })
    names(extracted.DAYMET) <- elements
    return(extracted.DAYMET)
  }
  
  data(tiles, envir = environment())
  
  template.latlon <- template %>% sp::spTransform(raster::projection(tiles))
  
  tile.ids <- tiles$TileID[!is.na(tiles %over% template.latlon)]
  
  tile.ids <- tile.ids[!is.na(tile.ids)]
  tile.ids <- unique(tile.ids)
  
  message("Area of interest includes ", length(tile.ids), " DAYMET tile(s).")
  
  # Download and crop tiles
  tiles <- lapply(tile.ids, function(tile) {
    return(FedData::get_daymet_tile(template = template,
                                    elements = elements,
                                    years = years,
                                    tileID = tile,
                                    raw.dir = raw.dir))
  })
  names(tiles) <- tile.ids
  
  # Mosaic all tiles
  if (length(tiles) > 1) {
    message("Mosaicking DAYMET tiles.")
    tiles <- foreach::foreach(element = elements) %do% {
      utils::flush.console()
      
      these.tiles <- lapply(tiles, "[[", element)
      these.tiles$fun <- mean
      names(these.tiles)[1:2] <- c("x", "y")
      out.tiles <- do.call(raster::mosaic, these.tiles)
      names(out.tiles) <- names(these.tiles$x)
      out.tiles
    }
  } else {
    tiles <- tiles[[1]]
  }
  names(tiles) <- elements
  
  tiles %>% mapply(x = ., y = names(tiles), FUN = function(x, y) {
    raster::writeRaster(x, paste0(extraction.dir, "/", label, "_DAYMET_", y, "_", min(years), "-", max(years), ".tif"), datatype = "FLT4S",
                        options = c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"), overwrite = T, setStatistics = FALSE)
  })
  
  tiles[[1]] %>% names() %>% readr::write_rds(paste0(extraction.dir, "/", label, "_DAYMET_layer_names.Rds"))
  
  
  return(tiles)
}

#' Download a netcdf tile from the 1-km DAYMET daily weather dataset.
#'
#' Tiles are specified by a tile ID, year, and element;
#' Tiles are downloaded in the NetCDF format. \code{download_daymet_tile} returns the path to the downloaded NetCDF tile files.
#'
#' @param tileID A numeric indicating the DAYMET tile ID number.
#' @param elements A character vector of elements to extract.\cr
#' The available elements are:\cr
#' dayl = Duration of the daylight period in seconds per day. This calculation is based on the period of the day during which the sun is above a hypothetical flat horizon.\cr
#' prcp = Daily total precipitation in millimeters per day, sum of all forms converted to water-equivalent. Precipitation occurrence on any given day may be ascertained.\cr
#' srad = Incident shortwave radiation flux density in watts per square meter, taken as an average over the daylight period of the day. NOTE: Daily total radiation (MJ/m2/day) can be calculated as follows: ((srad (W/m2) * dayl (s/day)) / l,000,000)\cr
#' swe = Snow water equivalent in kilograms per square meter. The amount of water contained within the snowpack.\cr
#' tmax = Daily maximum 2-meter air temperature in degrees Celsius.\cr
#' tmin = Daily minimum 2-meter air temperature in degrees Celsius.\cr
#' vp = Water vapor pressure in pascals. Daily average partial pressure of water vapor.\cr
#' @param years A numeric vector of years to extract.
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing. Defaults to './RAW/DAYMET/'.
#' @return A named list of character vectors, each representing the full local paths of the tile downloads.
#' @export
#' @keywords internal
#' @import ncdf4
#' @importFrom foreach foreach %do% %:%
download_daymet_tile <- function(tileID, elements, years, raw.dir) {
  
  # doParallel::registerDoParallel()
  out <- foreach::foreach(element = elements) %:%
    # foreach::foreach(year = sort(years), .combine = "c") %dopar% {
      foreach::foreach(year = sort(years), .combine = "c") %do% {
      destdir <- paste0(raw.dir, "/", tileID, "/", year)
      dir.create(destdir,
                 recursive = TRUE,
                 showWarnings = FALSE)
      url <- paste0("https://thredds.daac.ornl.gov/thredds/fileServer/ornldaac/1328/tiles/", year, "/", tileID, "_", year, "/",
                    element, ".nc")
      download_data(url = url,
                    destdir = destdir,
                    timestamping = FALSE,
                    nc = TRUE) %>%
        normalizePath(mustWork = T)
    }
  # doParallel::stopImplicitCluster()
  names(out) <- elements
  return(out)
}

#' Download and crop a netcdf tile from the 1-km DAYMET daily weather dataset.
#'
#' \code{get_daymet_tile} returns a list of \code{RasterBrick}s---one for each element---cropped within the specified \code{template}.
#' If template is not provided, returns the entire DAYMET tile.
#'
#' @param template A Raster* or Spatial* object to serve
#' as a template for cropping. If missing, entire tile is returned.
#' @param tileID A numeric indicating the DAYMET tile ID number.
#' @param elements A character vector of elements to extract.\cr
#' The available elements are:\cr
#' dayl = Duration of the daylight period in seconds per day. This calculation is based on the period of the day during which the sun is above a hypothetical flat horizon.\cr
#' prcp = Daily total precipitation in millimeters per day, sum of all forms converted to water-equivalent. Precipitation occurrence on any given day may be ascertained.\cr
#' srad = Incident shortwave radiation flux density in watts per square meter, taken as an average over the daylight period of the day. NOTE: Daily total radiation (MJ/m2/day) can be calculated as follows: ((srad (W/m2) * dayl (s/day)) / l,000,000)\cr
#' swe = Snow water equivalent in kilograms per square meter. The amount of water contained within the snowpack.\cr
#' tmax = Daily maximum 2-meter air temperature in degrees Celsius.\cr
#' tmin = Daily minimum 2-meter air temperature in degrees Celsius.\cr
#' vp = Water vapor pressure in pascals. Daily average partial pressure of water vapor.\cr
#' @param years A numeric vector of years to extract.
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing. Defaults to './RAW/NED/'.
#' @return A a list of \code{RasterBrick}s---one for each element---cropped within the specified \code{template}.
#' @export
#' @keywords internal
#' @importFrom foreach foreach %do%
#' @importFrom raster projection crop stack
#' @importFrom sp spTransform CRS
#' @importFrom lubridate year
#' @importFrom magrittr %<>% %>%
get_daymet_tile <- function(template, tileID, elements = NULL, years = NULL, raw.dir) {
  tmpdir <- tempfile()
  if (!dir.create(tmpdir))
    stop("failed to create my temporary directory")
  
  message("(Down)Loading DAYMET tile ", tileID)
  
  all.elements <- c("dayl", "prcp", "srad", "swe", "tmax", "tmin", "vp")
  elements <- tolower(elements)
  
  if (is.null(elements))
    elements <- all.elements
  
  missing.elements <- setdiff(elements, all.elements)
  if (length(missing.elements) > 0)
    warning("Elements not available: ", paste(missing.elements, collapse = ", "))
  elements <- setdiff(elements, missing.elements)
  if (length(elements) == 0)
    stop("No elements available")
  
  all.years <- 1980:(lubridate::year(Sys.time()) - 1)
  if (is.null(years))
    years <- all.years
  
  missing.years <- setdiff(years, all.years)
  if (length(missing.years) > 0)
    warning("Years not available: ", paste(missing.years, collapse = ", "))
  years <- setdiff(years, missing.years)
  if (length(years) == 0)
    stop("No years available")
  
  files <- download_daymet_tile(tileID = tileID, elements = elements, years = years, raw.dir = raw.dir)
  
  # doParallel::registerDoParallel()
  tiles <- foreach::foreach(element = files) %do% {
    # tiles <- foreach::foreach(element = files) %dopar% {
    tile <- foreach::foreach(file = element) %do% raster::brick(file)
    tile %<>% raster::stack(quick = TRUE)
    if (!is.null(template)) {
      tile <- tryCatch(tile %>% 
                         raster::crop(template %>% 
                                        sp::spTransform(tile %>% 
                                                          raster::projection() %>% 
                                                          sp::CRS()),
                                      snap = "out"),
                       error = function(e) {
                         tile %>% 
                           raster::crop(template %>% 
                                          sp::spTransform(tile %>% 
                                                            raster::projection() %>% 
                                                            sp::CRS()))
                       })
    }
  }
  # doParallel::stopImplicitCluster()
  names(tiles) <- elements
  
  unlink(tmpdir, recursive = TRUE)
  
  return(tiles)
}

#' The DAYMET tiles SpatialPolygonsDataFrame.
#'
#' A dataset containing the DAYMET tiles.
#'
#' @format A SpatialPolygonsDataFrame with 1060 features and 5 variables:
#' \describe{
#'   \item{TileID}{the numeric identifier of the tile}
#'   \item{XMin}{the minimum longitude of the tile}
#'   \item{XMax}{the maximum longitude of the tile}
#'   \item{YMin}{the minimum latitude of the tile}
#'   \item{YMax}{the maximum latitude of the tile}
#' }
#' @source \url{https://github.com/khufkens/daymetr/blob/master/data/DAYMET_grid.RData}
#' @keywords internal
"tiles"
