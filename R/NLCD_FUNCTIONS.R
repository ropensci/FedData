#' Download and crop the National Land Cover Database.
#'
#' \code{get_nlcd} returns a \code{RasterLayer} of NLCD data cropped to a given
#' template study area. \code{nlcd_colors} and \code{pal_nlcd} return the NLCD
#' legend and color palette, as available through the
#' [MLRC website](https://www.mrlc.gov/data/legends/national-land-cover-database-2016-nlcd2016-legend).
#'
#' @param template A sf, Raster* or Spatial* object to serve
#' as a template for cropping.
#' @param label A character string naming the study area.
#' @param year An integer representing the year of desired NLCD product.
#' Acceptable values are 2019 (default), 2016, 2011, 2008, 2006, 2004, and 2001.
#' @param dataset A character string representing type of the NLCD product.
#' Acceptable values are 'landcover' (default), 'impervious', and
#' 'canopy' (2016 and 2011, L48 only).
#' @param landmass A character string representing the landmass to be extracted
#' Acceptable values are 'L48' (lower 48 US states, the default),
#' 'AK' (Alaska, 2011 and 2016 only), 'HI' (Hawaii, 2001 only), and
#' 'PR' (Puerto Rico, 2001 only).
#' @param extraction.dir A character string indicating where the extracted
#' and cropped NLCD data should be put. The directory will be created if missing.
#' @param raster.options a vector of options for terra::writeRaster.
#' @param force.redo If an extraction for this template and label already exists,
#' should a new one be created?
#' @return A \code{RasterLayer} cropped to the bounding box of the template.
#' @export
#' @importFrom magrittr %>% %<>%
#' @examples
#' \dontrun{
#' # Extract data for the Mesa Verde National Park:
#'
#' # Get the NLCD (USA ONLY)
#' # Returns a raster
#' NLCD <-
#'   get_nlcd(
#'     template = FedData::meve,
#'     label = "meve",
#'     year = 2016
#'   )
#'
#' # Plot with raster::plot
#' plot(NLCD)
#' }
get_nlcd <- function(template,
                     label,
                     year = 2019,
                     dataset = c("landcover", "impervious", "canopy"),
                     landmass = "L48",
                     extraction.dir = paste0(
                       tempdir(),
                       "/FedData/"
                     ),
                     raster.options = c(
                       "COMPRESS=DEFLATE",
                       "ZLEVEL=9"
                     ),
                     force.redo = FALSE) {
  extraction.dir <-
    normalizePath(paste0(extraction.dir, "/."), mustWork = FALSE)

  template %<>% template_to_sf()

  dataset <- match.arg(dataset)
  dataset <- switch(dataset,
    landcover = "Land_Cover",
    impervious = "Impervious",
    canopy = "Tree_Canopy"
  )

  dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)

  outfile <-
    paste0(extraction.dir, "/", label, "_NLCD_", dataset, "_", year, ".tif")

  if (file.exists(outfile) & !force.redo) {
    return(raster::raster(outfile))
  }

  src <- "wcs"

  if (src == "wcs") {
    coverage <- paste0("NLCD_", year, "_", dataset, "_", landmass)
    source <- paste0("https://www.mrlc.gov/geoserver/mrlc_download/", coverage, "/wcs")


    # This code uses the (oft-changing) MRLC web services.
    if (source %>%
      httr::GET(query = list(
        service = "WCS",
        version = "2.0.1",
        request = "DescribeCoverage",
        coverageid = coverage
      )) %>%
      httr::status_code() %>%
      identical(200L) %>%
      magrittr::not()) {
      stop("No web coverage service at ", source, ". See available services at https://www.mrlc.gov/geoserver/ows?service=WCS&version=2.0.1&request=GetCapabilities")
    }

    template %<>%
      sf::st_transform("+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>%
      # sf::st_transform(3857) %>%
      sf::st_bbox()

    axis_labels <-
      source %>%
      httr::GET(
        query = list(
          service = "WCS",
          version = "2.0.1",
          request = "DescribeCoverage",
          coverageid = coverage
        )
      ) %>%
      httr::content(encoding = "UTF-8") %>%
      xml2::as_list() %$%
      CoverageDescriptions %$%
      CoverageDescription$boundedBy$Envelope %>%
      attr("axisLabels") %>%
      stringr::str_split(" ") %>%
      unlist()

    source %>%
      httr::GET(
        query = list(
          service = "WCS",
          version = "2.0.1",
          request = "GetCoverage",
          coverageid = coverage,
          subset = paste0(axis_labels[[1]], "(", template["xmin"], ",", template["xmax"], ")"),
          subset = paste0(axis_labels[[2]], "(", template["ymin"], ",", template["ymax"], ")")
        ),
        httr::write_disk(
          path = outfile,
          overwrite = TRUE
        )
      )
  }


  if (src == "cog") {
    source <- "https://storage.googleapis.com/feddata-r/nlcd/"
    file <- paste0(year, "_", dataset, "_", landmass, ".tif")

    path <- paste0(source, file)

    if (path %>%
      httr::HEAD() %>%
      httr::status_code() %>%
      identical(200L) %>%
      magrittr::not()) {
      stop(
        "NLCD data are not available for dataset '", dataset, "', year '", year,
        "', and landmass '", landmass,
        "'. Please see available datasets at https://www.mrlc.gov/data."
      )
    }

    template %<>%
      template_to_sf()

    out <-
      paste0("/vsicurl/", path) %>%
      terra::rast() %>%
      terra::crop(.,
        sf::st_transform(template, sf::st_crs(terra::crs(.))),
        snap = "out",
        filename = outfile,
        datatype = "INT1U",
        gdal = raster.options,
        overwrite = TRUE
      )
  }

  if (dataset == "Land_Cover") {
    out <-
      outfile %>%
      raster::raster() %>%
      raster::readAll() %>%
      raster::as.factor()

    raster::colortable(out) <- nlcd$Color

    suppressWarnings(
      levels(out) <-
        nlcd %>%
        as.data.frame()
    )

    out %>%
      terra::rast() %>%
      terra::writeRaster(outfile,
        datatype = "INT1U",
        gdal = raster.options,
        overwrite = TRUE
      )
  }

  return(raster::raster(outfile))
}

#' @export
#' @rdname get_nlcd
nlcd_colors <- function() {
  stats::na.omit(nlcd)
}

#' @export
#' @rdname get_nlcd
pal_nlcd <- function() {
  stats::na.omit(nlcd)
}
