#' Download and crop the 1 (~30 meter) or 1/3 (~10 meter) arc-second National Elevation Dataset.
#'
#' \code{get_ned} returns a \code{RasterLayer} of elevation data cropped to a given
#' template study area.
#'
#' @param template A Raster* or Spatial* object to serve
#' as a template for cropping.
#' @param label A character string naming the study area.
#' @param res A character string representing the desired resolution of the NED. '1'
#' indicates the 1 arc-second NED (the default), while '13' indicates the 1/3 arc-second dataset.
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing. Defaults to './RAW/NED/'.
#' @param extraction.dir A character string indicating where the extracted and cropped DEM should be put.
#' The directory will be created if missing. Defaults to './EXTRACTIONS/NED/'.
#' @param raster.options a vector of options for raster::writeRaster.
#' @param force.redo If an extraction for this template and label already exists, should a new one be created?
#' @return A \code{RasterLayer} DEM cropped to the extent of the template.
#' @export
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#' # Get the NED (USA ONLY)
#' # Returns a raster
#' NED <- get_ned(template = FedData::meve, label = "meve")
#'
#' # Plot with raster::plot
#' plot(NED)
#' }
get_ned <- function(template,
                    label,
                    res = "1",
                    raw.dir = paste0(tempdir(), "/FedData/raw/ned"),
                    extraction.dir = paste0(tempdir(), "/FedData/extractions/ned/", label, "/"),
                    raster.options = c(
                      "COMPRESS=DEFLATE",
                      "ZLEVEL=9"
                    ),
                    force.redo = F) {
  raw.dir <- normalizePath(paste0(raw.dir, "/."), mustWork = FALSE)
  extraction.dir <- normalizePath(paste0(extraction.dir, "/."), mustWork = FALSE)

  dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)

  if (file.exists(paste0(extraction.dir, "/", label, "_NED_", res, ".tif")) & !force.redo) {
    extracted.DEM <- raster::raster(paste0(extraction.dir, "/", label, "_NED_", res, ".tif"))
    return(extracted.DEM)
  }

  template %<>%
    template_to_sf() %>%
    sf::st_transform(4326)

  extent.latlon <- raster::extent(template)

  # Open USGS NED download service.
  # NED tiles are labeled by their northwest corner.
  # Thus, coordinate 36.42N, -105.71W is in grid n37w106
  wests <- seq(ceiling(abs(extent.latlon@xmax)), ceiling(abs(extent.latlon@xmin)))
  norths <- seq(ceiling(abs(extent.latlon@ymin)), ceiling(abs(extent.latlon@ymax)))

  tilesLocations <- as.matrix(expand.grid(norths, wests, stringsAsFactors = FALSE))

  message("Area of interest includes ", nrow(tilesLocations), " NED tiles.")

  # Download and crop tiles
  loc <- NULL
  tiles <-
    purrr::map(
      1:nrow(tilesLocations),
      function(loc) {
        return(
          tryCatch(
            get_ned_tile(
              template = template,
              res = res,
              tileNorthing = tilesLocations[loc, 1],
              tileWesting = tilesLocations[loc, 2],
              raw.dir = raw.dir
            ),
            error = function(e) {
              message("WARNING: ", e$message)
              return(NULL)
            },
            warning = function(w) NULL
          )
        )
      }
    )

  if (all(sapply(tiles, is.null))) {
    stop("No NED tiles are available for your study area. Please check your input data and internet connection.")
  }
  tiles <- tiles[which(!sapply(tiles, is.null))]


  # Mosaic all tiles
  if (length(tiles) > 1) {
    message("Mosaicking NED tiles.")
    utils::flush.console()

    tiles$fun <- mean
    names(tiles)[1:2] <- c("x", "y")
    tiles <- do.call(raster::mosaic, tiles)

    gc()
  } else {
    tiles <- tiles[[1]]
  }

  tiles %>%
    terra::crop(.,
      sf::st_transform(template, sf::st_crs(raster::crs(.))),
      snap = "out",
      filename = paste(extraction.dir, "/", label, "_NED_", res, ".tif", sep = ""),
      datatype = "FLT4S",
      options = raster.options,
      overwrite = T,
      setStatistics = FALSE
    )
}

#' Load and crop tile from the 1 (~30 meter) or 1/3 (~10 meter) arc-second National Elevation Dataset.
#'
#' \code{get_ned_tile} returns a \code{RasterLayer} cropped within the specified \code{template}.
#' If template is not provided, returns the entire NED tile.
#'
#' @param template A Raster* or Spatial* object to serve
#' as a template for cropping. If missing, entire tile is returned.
#' @param res A character string representing the desired resolution of the NED. '1'
#' indicates the 1 arc-second NED (the default), while '13' indicates the 1/3 arc-second dataset.
#' @param tileNorthing An integer representing the northing (latitude, in degrees north of the equator) of the northwest corner of the tile to
#' be downloaded.
#' @param tileWesting An integer representing the westing (longitude, in degrees west of the prime meridian) of the northwest corner of the tile to
#' be downloaded.
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing. Defaults to './RAW/NED/'.
#' @return A \code{RasterLayer} cropped within the specified \code{template}.
#' @export
#' @importFrom magrittr %>%
#' @keywords internal
get_ned_tile <- function(template = NULL, res = "1", tileNorthing, tileWesting, raw.dir) {
  tmpdir <- tempfile()
  if (!dir.create(tmpdir)) {
    stop("failed to create my temporary directory")
  }

  message("(Down)Loading NED tile for ", tileNorthing, "N and ", tileWesting, "W.")

  tileWesting <- formatC(tileWesting, width = 3, format = "d", flag = "0")
  tileNorthing <- formatC(tileNorthing, width = 2, format = "d", flag = "0")

  url <- paste0(
    "/vsicurl/https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/", res, "/TIFF/current/n", tileNorthing, "w", tileWesting, "/USGS_", res, "_n", tileNorthing, "w", tileWesting,
    ".tif"
  )

  url %>%
    terra::rast() %>%
    terra::crop(.,
      sf::st_transform(template, sf::st_crs(terra::crs(.))),
      snap = "out"
    ) %>%
    raster::raster()
}
