#' Download and crop the 1 (~30 meter) or 1/3 (~10 meter) arc-second National Elevation Dataset.
#'
#' `get_ned` returns a `SpatRaster` of elevation data cropped to a given
#' template study area.
#'
#' @param template An [`Simple Feature`][sf::sf]
#' or [`SpatRaster`][terra::SpatRaster] object to serve as a template for cropping.
#' @param label A character string naming the study area.
#' @param res A character string representing the desired resolution of the NED. '1'
#' indicates the 1 arc-second NED (the default), while '13' indicates the 1/3 arc-second dataset.
#' @param extraction.dir A character string indicating where the extracted and cropped DEM should be put.
#' The directory will be created if missing.
#' @param raster.options a vector of GDAL options passed to [terra::writeRaster].
#' @param force.redo If an extraction for this template and label already exists, should a new one be created?
#' @return A `SpatRaster` DEM cropped to the extent of the template.
#' @export
#' @importFrom magrittr %>%
#' @examples
#' \donttest{
#' # Get the NED (USA ONLY)
#' # Returns a `SpatRaster`
#' NED <-
#'   get_ned(
#'     template = FedData::meve,
#'     label = "meve"
#'   )
#'
#' # Plot with terra::plot
#' terra::plot(NED)
#' }
get_ned <- function(template,
                    label,
                    res = "1",
                    extraction.dir = file.path(
                      tempdir(),
                      "FedData",
                      "extractions",
                      "ned",
                      label
                    ),
                    raster.options = c(
                      "COMPRESS=DEFLATE",
                      "ZLEVEL=9"
                    ),
                    force.redo = FALSE) {
  extraction.dir <- normalizePath(extraction.dir, mustWork = FALSE)

  dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)

  outfile <-
    paste0(extraction.dir, "/", label, "_NED_", res, ".tif")

  if (file.exists(outfile) & !force.redo) {
    return(terra::rast(outfile))
  }

  template %<>%
    template_to_sf() %>%
    sf::st_transform(4326)

  extent.latlon <- sf::st_bbox(template)

  # Open USGS NED download service.
  # NED tiles are labeled by their northwest corner.
  # Thus, coordinate 36.42N, -105.71W is in grid n37w106
  wests <- seq(ceiling(abs(extent.latlon["xmax"])), ceiling(abs(extent.latlon["xmin"])))
  norths <- seq(ceiling(abs(extent.latlon["ymin"])), ceiling(abs(extent.latlon["ymax"])))

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
              tileWesting = tilesLocations[loc, 2]
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
    stop("No NED tiles are available for your study area.
         Please check your input data and internet connection.")
  }
  tiles <- tiles[which(!sapply(tiles, is.null))]

  # Mosaic all tiles
  if (length(tiles) > 1) {
    message("Mosaicking NED tiles.")
    utils::flush.console()

    tiles %<>%
      terra::sprc() %>%
      terra::mosaic(fun = "mean")
  } else {
    tiles <- tiles[[1]]
  }

  tiles %>%
    terra::crop(.,
      sf::st_transform(template, sf::st_crs(terra::crs(.))),
      snap = "out",
      filename = outfile,
      datatype = "FLT4S",
      gdal = raster.options,
      overwrite = T
    )

  return(terra::rast(outfile))
}

#' Load and crop tile from the 1 (~30 meter) or 1/3 (~10 meter) arc-second National Elevation Dataset.
#'
#' `get_ned_tile` returns a`SpatRaster` cropped within the specified `template`.
#' If template is not provided, returns the entire NED tile.
#'
#' @param template An [`Simple Feature`][sf::sf]
#' or [`SpatRaster`][terra::SpatRaster] object to serve as a template for cropping.
#' If missing, entire tile is returned.
#' @param res A character string representing the desired resolution of the NED. '1'
#' indicates the 1 arc-second NED (the default), while '13' indicates the 1/3 arc-second dataset.
#' @param tileNorthing An integer representing the northing (latitude, in degrees north of the equator) of the northwest corner of the tile to
#' be downloaded.
#' @param tileWesting An integer representing the westing (longitude, in degrees west of the prime meridian) of the northwest corner of the tile to
#' be downloaded.
#' @return A `SpatRaster` cropped to the extent of the template.
#' @export
#' @importFrom magrittr %>%
#' @keywords internal
get_ned_tile <- function(template = NULL, res = "1", tileNorthing, tileWesting) {
  message("(Down)Loading NED tile for ", tileNorthing, "N and ", tileWesting, "W.")

  tileWesting <- formatC(tileWesting, width = 3, format = "d", flag = "0")
  tileNorthing <- formatC(tileNorthing, width = 2, format = "d", flag = "0")

  url <- paste0(
    "/vsicurl/https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/", res,
    "/TIFF/current/n", tileNorthing,
    "w", tileWesting,
    "/USGS_", res,
    "_n", tileNorthing,
    "w", tileWesting,
    ".tif"
  )

  url %>%
    terra::rast() %>%
    terra::crop(.,
      sf::st_transform(template, sf::st_crs(terra::crs(.))),
      snap = "out"
    )
}
