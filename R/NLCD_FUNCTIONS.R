#' Download and crop the National Land Cover Database.
#'
#' \code{get_nlcd} returns a \code{RasterLayer} of NLCD data cropped to a given
#' template study area.
#'
#' @param template A Raster* or Spatial* object to serve
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
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#' # Extract data for the Mesa Verde National Park:
#'
#' # Get the NLCD (USA ONLY)
#' # Returns a raster
#' NLCD <- get_nlcd(paleocar::mvnp %>%
#'   sf::st_as_sf(),
#' label = "MVNP",
#' year = 2011,
#' landmass = "L48"
#' )
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
  source <- paste0("https://www.mrlc.gov/geoserver/mrlc_display/", coverage, "/ows")

  dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)

  outfile <- paste0(extraction.dir, "/", label, "_", coverage, "_nlcd.tif")

  if (file.exists(outfile) & !force.redo) {
    return(raster::raster(outfile))
  }

  if (source %>%
    httr::GET() %>%
    httr::status_code() %>%
    identical(200L) %>%
    magrittr::not()) {
    stop("No web coverage service at ", source, ". See available services at https://www.mrlc.gov/geoserver/ows?service=WCS&version=2.0.1&request=GetCapabilities")
  }

  template %<>%
    sf::st_transform(3857) %>%
    sf::st_bbox()

  source %>%
    httr::GET(
      query = list(
        service = "WCS",
        version = "2.0.1",
        request = "GetCoverage",
        coverageid = coverage,
        subset = paste0('X("', template["xmin"], '","', template["xmax"], '")'),
        subset = paste0('Y("', template["ymin"], '","', template["ymax"], '")')
      ),
      httr::write_disk(
        path = outfile,
        overwrite = TRUE
      )
    )

  outfile %>%
    raster::raster()
}
