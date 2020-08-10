#' Download and crop the NASS Cropland Data Layer.
#'
#' \code{get_nass} returns a \code{RasterLayer} of NASS Cropland Data Layer cropped to a given
#' template study area.
#'
#' @param template A Raster* or Spatial* object to serve
#' as a template for cropping.
#' @param label A character string naming the study area.
#' @param year An integer representing the year of desired NASS Cropland Data Layer product.
#' Acceptable values are 2007--the last year.
#' @param extraction.dir A character string indicating where the extracted and cropped NASS data should be put.
#' The directory will be created if missing.
#' @param raster.options a vector of options for raster::writeRaster.
#' @param force.redo If an extraction for this template and label already exists, should a new one be created?
#' @param progress Draw a progress bar when downloading?
#' @return A \code{RasterLayer} cropped to the bounding box of the template.
#' @export
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#' # Extract data for the Mesa Verde National Park:
#'
#' # Get the NASS (USA ONLY)
#' # Returns a raster
#' NASS <- get_nass(
#'   template = paleocar::mvnp %>%
#'     sf::st_as_sf(),
#'   label = "MVNP",
#'   year = 2011
#' )
#'
#' # Plot with raster::plot
#' plot(NASS)
#' }
get_nass <- function(template,
                     label,
                     year = 2019,
                     extraction.dir = paste0(tempdir(), "/FedData/extractions/nass/", label, "/"),
                     raster.options = c(
                       "COMPRESS=DEFLATE",
                       "ZLEVEL=9",
                       "INTERLEAVE=BAND"
                     ),
                     force.redo = FALSE,
                     progress = TRUE) {
  extraction.dir <- normalizePath(paste0(extraction.dir, "/."), mustWork = FALSE)

  if (!("sf" %in% class(template))) {
    template %<>% sf::st_as_sf()
  }

  layer <- paste0("cdl_", year)
  source <- "https://nassgeodata.gmu.edu/CropScapeService/wms_cdlall"

  dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)

  outfile <- paste0(extraction.dir, "/", layer, ".tif")

  if (file.exists(outfile) & !force.redo) {
    return(raster::raster(outfile))
  }

  if (source %>%
    httr::GET() %>%
    httr::status_code() %>%
    identical(200L) %>%
    magrittr::not()) {
    stop("No web coverage service at ", source, ". See available services at https://nassgeodata.gmu.edu/CropScapeService/wms_cdlall?service=WCS&version=2.0.1&request=GetCapabilities")
  }



  template %<>%
    sf::st_transform("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m") %>%
    sf::st_bbox()

  # # Convert to floating grid
  # origin <- c(-2356095, 276915)
  # template["xmin"] <- floor((template["xmin"] - origin[[1]]) / 30)
  # template["xmax"] <- ceiling((template["xmax"] - origin[[1]]) / 30)
  # template["ymin"] <- floor((template["ymin"] - origin[[2]]) / -30)
  # template["ymax"] <- ceiling((template["ymax"] - origin[[2]]) / -30)
  #
  # template %<>%
  #   split_bbox(
  #     x = 2047 * 30,
  #     y = 2047 * 30
  #   )
  #
  # if (progress) {
  #   pb <- progress::progress_bar$new(total = length(template))
  # }
  #
  # out <-
  #   template %>%
  #   list() %>%
  #   purrr::map(function(x) {
  #     if (progress) {
  #       pb$tick()
  #     }
  #
  #     tf <- tempfile(fileext = ".tif")
  #
  #     response <-
  #       x %>%
  #       as.list() %$%
  #       httr::GET(source,
  #                 query = list(
  #                   service = "wcs",
  #                   version = "1.0.0",
  #                   request = "getcoverage",
  #                   coverage = layer,
  #                   crs = "epsg:102004",
  #                   resx = 30,
  #                   resy = 30,
  #                   bbox = paste0(c(xmin, ymin, xmax, ymax), collapse = ","),
  #                   format = "gtiff"
  #                 ),
  #                 httr::write_disk(
  #                   path = tf,
  #                   overwrite = TRUE
  #                 )
  #       )
  #
  #     if (httr::headers(response)$`content-type` != "image/tiff") {
  #       stop(response %>%
  #              httr::content(encoding = "UTF-8") %>%
  #              xml2::as_list() %$%
  #              ExceptionReport$Exception$ExceptionText[[1]])
  #     }
  #
  #     tf %>%
  #       raster::raster() %>%
  #       raster::readAll()
  #   }) %>%
  #   do.call(raster::merge, .) %>%
  #   raster::as.factor()

  tf <- tempfile(fileext = ".tif")

  response <-
    template %>%
    as.list() %$%
    httr::GET(source,
      query = list(
        service = "wcs",
        version = "1.0.0",
        request = "getcoverage",
        coverage = layer,
        crs = "epsg:102004",
        resx = 30,
        resy = 30,
        bbox = paste0(c(xmin, ymin, xmax, ymax), collapse = ","),
        format = "gtiff"
      ),
      httr::write_disk(
        path = tf,
        overwrite = TRUE
      )
    )

  if (httr::headers(response)$`content-type` != "image/tiff") {
    stop(response %>%
      httr::content(encoding = "UTF-8") %>%
      xml2::as_list() %$%
      ExceptionReport$Exception$ExceptionText[[1]])
  }

  out <-
    tf %>%
    raster::raster() %>%
    raster::readAll() %>%
    raster::as.factor()

  raster::colortable(out) <- nass$Color

  suppressWarnings(
    levels(out) <-
      nass %>%
      as.data.frame()
  )

  out %<>%
    raster::writeRaster(outfile,
      datatype = "INT1U",
      options = raster.options,
      overwrite = T,
      setStatistics = FALSE
    )

  return(out)
}
