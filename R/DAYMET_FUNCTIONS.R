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
#' @param region The name of a region. The available regions are:\cr
#' na = North America\cr
#' hi = Hawaii\cr
#' pr = Puerto Rico\cr
#' @param tempo The frequency of the data. The available tempos are:\cr
#' day = Daily data\cr
#' mon = Monthly summary data\cr
#' ann = Annual summary data\cr
#' @param extraction.dir A character string indicating where the extracted and cropped DEM should be put.
#' Defaults to a temporary directory.
#' @param raster.options a vector of options for raster::writeRaster.
#' @param force.redo If an extraction for this template and label already exists in extraction.dir,
#' should a new one be created?
#' @param progress Draw a progress bar when downloading?
#' @return A named list of \code{RasterBrick}s of weather data cropped to the extent of the template.
#' @importFrom magrittr %>% %$% %<>% %T>%
#' @export
#' @examples
#' \dontrun{
#' # Extract data for the Village Ecodynamics Project 'VEPIIN' study area:
#' # http://village.anth.wsu.edu
#' vepPolygon <- polygon_from_extent(raster::extent(672800, 740000, 4102000, 4170000),
#'   proj4string = "+proj=utm +datum=NAD83 +zone=12"
#' )
#'
#' # Get the DAYMET (North America only)
#' # Returns a list of raster bricks
#' DAYMET <- get_daymet(
#'   template = vepPolygon,
#'   label = "VEPIIN",
#'   elements = c("prcp", "tmin", "tmax"),
#'   years = 1980:1985
#' )
#'
#' # Plot with raster::plot
#' plot(DAYMET$tmin$X1985.10.23)
#' }
get_daymet <- function(template,
                       label,
                       elements = c("dayl", "prcp", "srad", "swe", "tmax", "tmin", "vp"),
                       years = 1980:(lubridate::year(Sys.time()) - 1),
                       region = "na",
                       tempo = "day",
                       extraction.dir = paste0(tempdir(), "/FedData/extractions/daymet/", label, "/"),
                       raster.options = c(
                         "COMPRESS=DEFLATE",
                         "ZLEVEL=9",
                         "INTERLEAVE=BAND"
                       ),
                       force.redo = F,
                       progress = TRUE) {
  extraction.dir %<>%
    paste0("/") %>%
    normalizePath(mustWork = FALSE) %T>%
    dir.create(
      showWarnings = FALSE,
      recursive = TRUE
    )

  all.regions <- c("na", "hi", "pr")

  if (length(region) > 1) {
    stop("Please select only one region.")
  }
  if (!(region %in% all.regions)) {
    stop("`region` must be one of c('na', 'hi', 'pr').")
  }



  all.tempos <- c("day", "mon", "ann")

  if (length(tempo) > 1) {
    stop("Please select only one tempo.")
  }
  if (!(tempo %in% all.tempos)) {
    stop("`tempo` must be one of c('day', 'mon', 'ann').")
  }



  all.elements <- c("dayl", "prcp", "srad", "swe", "tmax", "tmin", "vp")

  elements %<>% tolower()

  missing.elements <- setdiff(elements, all.elements)
  if (length(missing.elements) > 0) {
    stop("Elements not available: ", paste(missing.elements, collapse = ", "), ".\n
         Please select among c('dayl', 'prcp', 'srad', 'swe', 'tmax', 'tmin', 'vp').")
  }

  if (tempo != "day" & length(base::setdiff(elements, c("prcp", "tmax", "tmin", "vp"))) > 0) {
    warning("Only elements in c('prcp', 'tmax', 'tmin', 'vp') 
            are available for monthly or annual data.")

    elements <- base::intersect(elements, c("prcp", "tmax", "tmin", "vp"))
  }


  all.years <- 1980:(lubridate::year(Sys.time()) - 1)

  missing.years <- setdiff(years, all.years)
  if (length(missing.years) > 0) {
    stop("Years not available: ", paste(missing.years, collapse = ", "), ".\n
         Please select among 1980--", lubridate::year(Sys.time()) - 1, ".")
  }

  years <- setdiff(years, missing.years)
  if (length(years) == 0) {
    stop("No years available")
  }



  template_bbox <-
    template %>%
    sf::st_bbox() %>%
    sf::st_as_sfc() %>%
    sf::st_transform(4326) %>%
    sf::st_bbox()


  out.files <- paste0(extraction.dir, "/", label, "_", elements, "_", tempo, ".tif")

  if (!force.redo & all(file.exists(out.files))) {
    out.files %>%
      purrr::map(function(x) {
        x %>%
          raster::brick() %>%
          raster::readAll()
      }) %>%
      magrittr::set_names(elements)
  }

  if (progress) {
    pb <- progress::progress_bar$new(total = length(elements) * length(years))
  }

  out <-
    elements %>%
    magrittr::set_names(., .) %>%
    purrr::map(function(element) {
      years %>%
        sort() %>%
        magrittr::set_names(., .) %>%
        purrr::map(function(year) {
          if (progress) {
            pb$tick()
          }

          download_daymet_thredds(
            bbox = template_bbox %>%
              paste0(collapse = ","),
            element = element,
            year = year,
            region = region,
            tempo = tempo
          )
        }) %>%
        raster::stack()
    }) %>%
    purrr::map(function(x) {
      t_bb <-
        template_bbox %>%
        sf::st_as_sfc() %>%
        sf::st_transform(
          raster::projection(x)
        ) %>%
        sf::st_as_sf()

      x %>%
        raster::crop(t_bb, snap = "out")
    })

  out %>%
    purrr::iwalk(function(x, i) {
      raster::writeRaster(x,
        paste0(extraction.dir, "/", label, "_", i, "_", tempo, ".tif"),
        options = raster.options,
        overwrite = T,
        setStatistics = FALSE
      )
    })

  return(out)
}

#' Download the 1-km DAYMET daily weather dataset for a region as a netcdf.
#'
#' Data are downloaded in the NetCDF format. \code{download_daymet_thredds} returns the path to the downloaded NetCDF file.
#'
#' @param bbox the bounding box in WGS84 coordinates as a comma-separated character vector
#' "xmin,ymin,xmax,ymax"
#' @param element An element to extract.\cr
#' The available elements are:\cr
#' dayl = Duration of the daylight period in seconds per day. This calculation is based on the period of the day during which the sun is above a hypothetical flat horizon.\cr
#' prcp = Daily total precipitation in millimeters per day, sum of all forms converted to water-equivalent. Precipitation occurrence on any given day may be ascertained.\cr
#' srad = Incident shortwave radiation flux density in watts per square meter, taken as an average over the daylight period of the day. NOTE: Daily total radiation (MJ/m2/day) can be calculated as follows: ((srad (W/m2) * dayl (s/day)) / l,000,000)\cr
#' swe = Snow water equivalent in kilograms per square meter. The amount of water contained within the snowpack.\cr
#' tmax = Daily maximum 2-meter air temperature in degrees Celsius.\cr
#' tmin = Daily minimum 2-meter air temperature in degrees Celsius.\cr
#' vp = Water vapor pressure in pascals. Daily average partial pressure of water vapor.\cr
#' @param year An integer year to extract.
#' @param region The name of a region. The available regions are:\cr
#' na = North America\cr
#' hi = Hawaii\cr
#' pr = Puerto Rico\cr
#' @param tempo The frequency of the data. The available tempos are:\cr
#' day = Daily data\cr
#' mon = Monthly summary data\cr
#' ann = Annual summary data\cr
#' @return A named list of character vectors, each representing the full local paths of the tile downloads.
#' @importFrom magrittr %>% %$% %<>% %T>%
#' @keywords internal
download_daymet_thredds <-
  function(bbox,
           element,
           year,
           region,
           tempo) {
    tf <- tempfile(fileext = ".nc")

    element <-
      c(
        "dayl" = NA,
        "prcp" = "ttl",
        "srad" = NA,
        "swe" = NA,
        "tmax" = "avg",
        "tmin" = "avg",
        "vp" = "avg"
      )[element]

    region <-
      c(
        "na" = "na",
        "hi" = "hawaii",
        "pr" = "puertorico"
      )[region]

    tempo <-
      c(
        "ann" = 1343,
        "mon" = 1345,
        "day" = 1328
      )[tempo]



    if (tempo == 1328) {
      url <- paste0(
        "https://thredds.daac.ornl.gov/thredds/wcs/ornldaac/",
        tempo,
        "/",
        year,
        "/daymet_v3_",
        names(element),
        "_",
        year,
        "_",
        region,
        ".nc4"
      )
    } else {
      url <- paste0(
        "https://thredds.daac.ornl.gov/thredds/wcs/ornldaac/",
        tempo,
        "/daymet_v3_",
        names(element),
        "_",
        names(tempo),
        element,
        "_",
        year,
        "_",
        names(region),
        ".nc4"
      )
    }

    response <-
      httr::GET(url,
        query = list(
          service = "WCS",
          version = "1.0.0",
          request = "GetCoverage",
          format = "NetCDF3",
          coverage = names(element),
          bbox = bbox
        ),
        httr::write_disk(
          path = tf,
          overwrite = TRUE
        )
      )

    if (httr::headers(response)$`content-type` != "application/x-netcdf") {
      stop(response %>%
        httr::content(type = "text/xml", encoding = "UTF-8") %>%
        xml2::as_list() %$%
        ServiceExceptionReport$ServiceException[[1]])
    }

    suppressWarnings(
      out <-
        tf %>%
        raster::stack()
    )

    raster::projection(out) <-
      "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

    out %<>%
      raster::setExtent((sf::st_bbox(out) * 1000)[c("xmin", "xmax", "ymin", "ymax")])

    out %>%
      raster::readAll()
  }
