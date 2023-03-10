#' Download and crop the NASS Cropland Data Layer.
#'
#' \code{get_nass_cdl} returns a \code{RasterLayer} of NASS Cropland Data Layer cropped to a given
#' template study area.
#'
#' @param template A Raster* or Spatial* object to serve
#' as a template for cropping.
#' @param label A character string naming the study area.
#' @param year An integer representing the year of desired NASS Cropland Data Layer product.
#' Acceptable values are 2007--the last year.
#' @param extraction.dir A character string indicating where the extracted and cropped NASS data should be put.
#' The directory will be created if missing.
#' @param raster.options a vector of options for terra::writeRaster.
#' @param force.redo If an extraction for this template and label already exists, should a new one be created?
#' @param progress Draw a progress bar when downloading?
#' @param ... Other parameters passed on to [get_nass_cdl].
#' @return A \code{RasterLayer} cropped to the bounding box of the template.
#' @export
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#' # Extract data for the Mesa Verde National Park:
#'
#' # Get the NASS CDL (USA ONLY)
#' # Returns a raster
#' NASS <-
#'   get_nass_cdl(
#'     template = FedData::meve,
#'     label = "meve",
#'     year = 2011
#'   )
#'
#' # Plot with raster::plot
#' plot(NASS)
#' }
get_nass_cdl <- function(template,
                         label,
                         year = 2019,
                         extraction.dir = paste0(tempdir(), "/FedData/"),
                         raster.options = c(
                           "COMPRESS=DEFLATE",
                           "ZLEVEL=9",
                           "INTERLEAVE=BAND"
                         ),
                         force.redo = FALSE,
                         progress = TRUE) {
  extraction.dir <- normalizePath(paste0(extraction.dir, "/."), mustWork = FALSE)

  template %<>% template_to_sf()

  layer <- paste0("cdl_", year)
  source <- "https://nassgeodata.gmu.edu/CropScapeService/wms_cdlall"

  dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)

  outfile <- paste0(extraction.dir, "/", label, "_", layer, "_nass.tif")

  if (file.exists(outfile) & !force.redo) {
    out <-
      terra::rast(outfile)

    outrast <-
      raster::raster(out)

    raster::colortable(outrast) <-
      nass$Color

    return(outrast)
  }

  if (source %>%
    httr::GET() %>%
    httr::status_code() %>%
    identical(200L) %>%
    magrittr::not()) {
    stop("No web coverage service at ", source, ". See available services at https://nassgeodata.gmu.edu/CropScapeService/wms_cdlall?service=WCS&version=2.0.1&request=GetCapabilities")
  }



  template %<>%
    sf::st_transform("+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs") %>%
    sf::st_bbox()

  tf <- tempfile(fileext = ".tif")

  response <-
    source %>%
    httr::GET(
      query = list(
        service = "WCS",
        version = "2.0.1",
        request = "GetCoverage",
        coverageid = layer,
        subset = paste0("x(", template["xmin"], ",", template["xmax"], ")"),
        subset = paste0("y(", template["ymin"], ",", template["ymax"], ")"),
        format = "image/tiff"
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
    terra::rast() %>%
    terra::as.factor()

  levels(out) <- dplyr::select(nass, ID, `Land Cover`)
  terra::coltab(out) <- nass$Color

  out %>%
    terra::writeRaster(
      filename = outfile,
      overwrite = TRUE,
      datatype = "INT1U",
      gdal = raster.options
    )


  out <-
    terra::rast(outfile)

  outrast <-
    raster::raster(out)

  raster::colortable(outrast) <-
    nass$Color

  return(outrast)
}

#' @export
#' @rdname get_nass_cdl
get_nass <- function(template, label, ...) {
  lifecycle::deprecate_warn("3.0.0", "get_nass()", "get_nass_cdl()",
    details = "`get_nass()` has become `get_nass_cdl()` to clarify the dataset provided. See `?get_nass_cdl`."
  )
  get_nass_cdl(template, label, ...)
}

#' @export
#' @rdname get_nass_cdl
get_cdl <- function(template, label, ...) {
  get_cdl(template, label, ...)
}

#' @export
#' @rdname get_nass_cdl
cdl_colors <- function() {
  nass
}
