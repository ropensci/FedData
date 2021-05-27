#' Download and crop the National Land Cover Database.
#'
#' \code{get_nlcd} returns a \code{RasterLayer} of NLCD data cropped to a given
#' template study area.
#'
#' NOTE: Prior to FedData version 3.0.0, the `get_nlcd` function returned
#' data in the web Mercator coordinate reference system available through
#' the [MRLC web mapping services](https://www.mrlc.gov/geoserver/web/), rather
#' than data in the NLCD's native projection (a flavor of North American Albers).
#' This function now returns data in the original CRS.
#'
#' @param template A sf, Raster* or Spatial* object to serve
#' as a template for cropping.
#' @param label A character string naming the study area.
#' @param year An integer representing the year of desired NLCD product.
#' Acceptable values are 2016 (default), 2011, 2008 (landcover only), 2006, 2004, and 2001.
#' @param dataset A character string representing type of the NLCD product.
#' Acceptable values are 'Impervious', 'Land_Cover', 'Tree_Canopy' (2011 and 2016 only),
#' @param landmass A character string representing the landmass to be extracted
#' Acceptable values are 'L48' (lower 48 US states, the default), 'AK' (Alaska), 'HI' (Hawaii), and 'PR' (Puerto Rico).
#' @param extraction.dir A character string indicating where the extracted and cropped NLCD data should be put.
#' The directory will be created if missing.
#' @param force.redo If an extraction for this template and label already exists, should a new one be created?
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
#'     year = 2011,
#'     landmass = "L48"
#'   )
#'
#' # Plot with raster::plot
#' plot(NLCD)
#' }
get_nlcd <- function(template,
                     label,
                     year = 2016,
                     dataset = "Land_Cover",
                     landmass = "L48",
                     extraction.dir = paste0(
                       tempdir(),
                       "/FedData/"
                     ),
                     force.redo = F) {
  extraction.dir <- normalizePath(paste0(extraction.dir, "/."), mustWork = FALSE)

  template %<>% template_to_sf()

  coverage <- paste0("NLCD_", year, "_", dataset, "_", landmass)
  source <- "https://www.mrlc.gov/geoserver/wcs"

  dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)

  outfile <- paste0(extraction.dir, "/", label, "_", coverage, "_nlcd.tif")

  if (file.exists(outfile) & !force.redo) {
    return(raster::raster(outfile))
  }

  # if (source %>%
  #   httr::GET() %>%
  #   httr::status_code() %>%
  #   identical(200L) %>%
  #   magrittr::not()) {
  #   stop("No web coverage service at ", source, ". See available services at https://www.mrlc.gov/geoserver/ows?service=WCS&version=2.0.1&request=GetCapabilities")
  # }

  template %<>%
    sf::st_transform("+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>%
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
    httr::content() %>%
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

  outfile %>%
    raster::raster()
}
