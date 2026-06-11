#' Download and crop the National Land Cover Database.
#'
#' \code{get_nlcd} returns a [`SpatRaster`][terra::SpatRaster] of NLCD data cropped to a given
#' template study area. \code{nlcd_colors} and \code{pal_nlcd} return the NLCD
#' legend and color palette, as available through the
#' [MLRC website](https://www.mrlc.gov/data/legends/national-land-cover-database-class-legend-and-description).
#'
#' @param template An [`Simple Feature`][sf::sf]
#' or [`terra`][terra::SpatRaster] object to serve as a template for cropping.
#' @param label A character string naming the study area.
#' @param year An integer representing the year of desired NLCD product.
#' Acceptable values are 2019 (default), 2016, 2011, 2008, 2006, 2004, and 2001.
#' The L48 data set for 2021 is corrupted on the NLCD Mapserver, and is thus
#' not available through FedData.
#' @param dataset A character string representing type of the NLCD product.
#' Acceptable values are 'landcover' (default), 'impervious', and
#' 'canopy'.
#' @param landmass A character string representing the landmass to be extracted
#' Acceptable values are 'L48' (lower 48 US states, the default),
#' 'AK' (Alaska, 2001, 2011 and 2016 only), 'HI' (Hawaii, 2001 only), and
#' 'PR' (Puerto Rico, 2001 only).
#' @param extraction.dir A character string indicating where the extracted
#' and cropped NLCD data should be put. The directory will be created if missing.
#' @param raster.options a vector of GDAL options passed to [terra::writeRaster].
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
#' # Plot with terra::plot
#' terra::plot(NLCD)
#' }
get_nlcd <- function(template,
                     label,
                     year = 2021,
                     dataset = "landcover",
                     landmass = "L48",
                     extraction.dir = file.path(
                       tempdir(),
                       "FedData",
                       "extractions",
                       "nlcd",
                       label
                     ),
                     raster.options = c(
                       "COMPRESS=DEFLATE",
                       "ZLEVEL=9"
                     ),
                     force.redo = FALSE) {
  extraction.dir <-
    normalizePath(paste0(extraction.dir, "/."), mustWork = FALSE)

  template %<>%
    template_to_sf()

  dataset <-
    match.arg(dataset,
      choices = c(
        "landcover",
        "impervious",
        "canopy"
      )
    )

  dataset <-
    switch(dataset,
      landcover = "Land_Cover",
      impervious = "Impervious",
      canopy = "Tree_Canopy"
    )

  dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)

  outfile <-
    paste0(extraction.dir, "/", label, "_NLCD_", dataset, "_", year, ".tif")

  if (file.exists(outfile) & !force.redo) {
    return(terra::rast(outfile))
  }

  nlcd_base_url <-
    "https://www.mrlc.gov/geoserver/ows"

  if (dataset == "Tree_Canopy" & landmass == "L48") {
    # Because MRLC did an update, and of course they didn't use the same naming convention
    coverage <- paste0("mrlc_download__nlcd_tcc_conus_", year, "_v2021-4")
  } else {
    coverage <- paste0("mrlc_download__NLCD_", year, "_", dataset, "_", landmass)
  }



  source <-
    httr::modify_url(
      nlcd_base_url,
      query =
        list(
          version = "2.0.1",
          coverageid = coverage
        )
    )

  # This code uses the (oft-changing) MRLC web services.
  if (nlcd_base_url %>%
    httr::GET(query = list(
      service = "WCS",
      version = "2.0.1",
      request = "DescribeCoverage",
      coverageid = coverage
    )) %>%
    httr::status_code() %>%
    identical(200L) %>%
    magrittr::not()) {
    stop(
      "No web coverage service at ",
      source, ". See available services at ",
      nlcd_base_url, "?service=WCS&version=2.0.1&request=GetCapabilities"
    )
  }

  out <-
    tryCatch(
      paste0("WCS:", source) |>
        terra::rast(),
      error = function(e) {
        # Some coverages (e.g., NLCD_2016_Land_Cover_AK) advertise a
        # CRS code that GDAL cannot interpret, which makes the WCS driver
        # error out. Fall back to requesting the coverage directly.
        tryCatch(
          get_nlcd_wcs_direct(
            nlcd_base_url = nlcd_base_url,
            coverage = coverage,
            template = template
          ),
          error = function(e) {
            stop(
              "Web coverage service at ",
              source, " is corrupted. You may have better luck downloading the file at https://www.mrlc.gov/data."
            )
          }
        )
      }
    )

  out <-
    out |>
    terra::crop(
      sf::st_transform(
        template,
        terra::crs(out)
      ) |>
        terra::vect(),
      snap = "out",
      mask = TRUE
    ) |>
    magrittr::set_names(coverage)

  if (dataset == "Land_Cover") {
    out <-
      out |>
      terra::as.factor()

    levels(out) <-
      nlcd_colors() %>%
      as.data.frame()

    terra::coltab(out) <-
      nlcd_colors() |>
      dplyr::transmute(
        value = as.integer(ID),
        color = Color
      ) |>
      as.data.frame()
  }

  terra::writeRaster(
    x = out,
    filename = outfile,
    datatype = "INT1U",
    gdal = raster.options,
    overwrite = TRUE
  )

  return(terra::rast(outfile))
}

# Request an NLCD coverage directly with a WCS GetCoverage request,
# subset to the (snapped) bounding box of the template.
#
# This is a fallback for coverages whose DescribeCoverage response
# advertises a CRS that GDAL cannot interpret (e.g., the Alaska NLCD
# coverages, which advertise the nonexistent code EPSG:675225); the
# GeoTIFFs returned by GetCoverage embed the true CRS. Because the
# advertised CRS cannot be used to project the template into coverage
# coordinates, a small probe request at the corner of the coverage
# envelope is made first to learn the true CRS.
get_nlcd_wcs_direct <- function(nlcd_base_url, coverage, template) {
  describe <-
    httr::GET(
      nlcd_base_url,
      query = list(
        service = "WCS",
        version = "2.0.1",
        request = "DescribeCoverage",
        coverageid = coverage
      )
    ) %>%
    httr::content(encoding = "UTF-8")

  envelope <-
    xml2::xml_find_first(describe, "//*[local-name()='Envelope']")

  axes <-
    strsplit(xml2::xml_attr(envelope, "axisLabels"), " ")[[1]]

  lower <-
    envelope |>
    xml2::xml_find_first("./*[local-name()='lowerCorner']") |>
    xml2::xml_text() |>
    strsplit(" ") |>
    unlist() |>
    as.numeric()

  upper <-
    envelope |>
    xml2::xml_find_first("./*[local-name()='upperCorner']") |>
    xml2::xml_text() |>
    strsplit(" ") |>
    unlist() |>
    as.numeric()

  cells <-
    describe |>
    xml2::xml_find_first(
      "//*[local-name()='GridEnvelope']/*[local-name()='high']"
    ) |>
    xml2::xml_text() |>
    strsplit(" ") |>
    unlist() |>
    as.numeric() |>
    magrittr::add(1)

  res <- (upper - lower) / cells

  get_coverage <- function(xmin, xmax, ymin, ymax) {
    tmp <- tempfile(fileext = ".tif")

    resp <-
      httr::GET(
        nlcd_base_url,
        query = list(
          service = "WCS",
          version = "2.0.1",
          request = "GetCoverage",
          coverageid = coverage,
          format = "image/geotiff",
          subset = paste0(axes[[1]], "(", xmin, ",", xmax, ")"),
          subset = paste0(axes[[2]], "(", ymin, ",", ymax, ")")
        ),
        httr::write_disk(tmp, overwrite = TRUE)
      )

    if (
      httr::http_error(resp) ||
        !grepl("tiff", httr::headers(resp)[["content-type"]])
    ) {
      stop("GetCoverage request for ", coverage, " failed.")
    }

    terra::rast(tmp)
  }

  # Probe the corner of the coverage to learn its true CRS
  probe <-
    get_coverage(
      xmin = lower[[1]],
      xmax = lower[[1]] + 3 * res[[1]],
      ymin = lower[[2]],
      ymax = lower[[2]] + 3 * res[[2]]
    )

  bbox <-
    template |>
    sf::st_transform(terra::crs(probe)) |>
    sf::st_bbox()

  # Snap the template bounding box outward to the coverage grid,
  # and clamp it to the coverage envelope
  bbox <-
    c(
      xmin = lower[[1]] +
        floor((bbox[["xmin"]] - lower[[1]]) / res[[1]]) * res[[1]],
      ymin = lower[[2]] +
        floor((bbox[["ymin"]] - lower[[2]]) / res[[2]]) * res[[2]],
      xmax = lower[[1]] +
        ceiling((bbox[["xmax"]] - lower[[1]]) / res[[1]]) * res[[1]],
      ymax = lower[[2]] +
        ceiling((bbox[["ymax"]] - lower[[2]]) / res[[2]]) * res[[2]]
    )

  bbox <-
    c(
      xmin = max(bbox[["xmin"]], lower[[1]]),
      ymin = max(bbox[["ymin"]], lower[[2]]),
      xmax = min(bbox[["xmax"]], upper[[1]]),
      ymax = min(bbox[["ymax"]], upper[[2]])
    )

  if (
    bbox[["xmin"]] >= bbox[["xmax"]] ||
      bbox[["ymin"]] >= bbox[["ymax"]]
  ) {
    stop("The provided template does not overlap the coverage ", coverage, ".")
  }

  get_coverage(
    xmin = bbox[["xmin"]],
    xmax = bbox[["xmax"]],
    ymin = bbox[["ymin"]],
    ymax = bbox[["ymax"]]
  )
}

#' Download and crop the Annual National Land Cover Database.
#'
#' \code{get_nlcd_annual} returns a [`SpatRaster`][terra::SpatRaster] of NLCD data cropped to a given
#' template study area. The Annual NLCD is currently only available for the conterminous United States.
#' More information about the Annual NLCD product is available on the
#' [Annual NLCD web page](https://www.mrlc.gov/data/project/annual-nlcd).
#'
#' Data are downloaded using the
#' [USGS Web Coverage Service](https://dmsdata.cr.usgs.gov/geoserver/web/)
#' for the Annual NLCD, in the native coordinate reference system and resolution
#' (CONUS Albers, EPSG:5070, at 30 m) and snapped to the native grid.
#'
#' @param template An [`Simple Feature`][sf::sf]
#' or [`terra`][terra::SpatRaster] object to serve as a template for cropping.
#' @param label A character string naming the study area.
#' @param year An integer vector representing the year of desired NLCD product.
#' Acceptable values are currently 1985 through 2024 (defaults to 2024).
#' @param product A character vector representing type of the NLCD product.
#' Defaults to 'LndCov' (Land Cover).\cr
#' LndCov = Land Cover\cr
#' LndChg = Land Cover Change\cr
#' LndCnf = Land Cover Confidence\cr
#' FctImp = Fractional Impervious Surface\cr
#' ImpDsc = Impervious Descriptor\cr
#' SpcChg = Spectral Change Day of Year\cr
#' @param region A character string representing the region to be extracted
#' Acceptable values are 'CU' (Conterminous US, the default),
#' 'AK' (Alaska), and 'HI' (Hawaii). **Currently, only 'CU' is available.**
#' @param collection An integer representing the collection number.
#' **Currently, only '1' is available.**
#' @param version An integer representing the version number.
#' **Currently, only '1' is available.**
#' @param extraction.dir A character string indicating where the extracted
#' and cropped NLCD data should be put. The directory will be created if missing.
#' @param raster.options a vector of GDAL options passed to [terra::writeRaster].
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
#' NLCD_ANNUAL <-
#'   get_nlcd_annual(
#'     template = FedData::meve,
#'     label = "meve",
#'     year = 2020,
#'     product =
#'       c(
#'         "LndCov",
#'         "LndChg",
#'         "LndCnf",
#'         "FctImp",
#'         "ImpDsc",
#'         "SpcChg"
#'       )
#'   )
#'
#' NLCD_ANNUAL
#' }
get_nlcd_annual <-
  function(template,
           label,
           year = 2024,
           product = "LndCov",
           region = "CU",
           collection = 1,
           version = 1,
           extraction.dir = file.path(
             tempdir(),
             "FedData",
             "extractions",
             "nlcd_annual",
             label
           ),
           raster.options = c(
             "COMPRESS=DEFLATE",
             "ZLEVEL=9"
           ),
           force.redo = FALSE) {
    if (any(!(year %in% 1985:2024))) {
      stop("'year' must only contain integers between 1985 and 2024.")
    }

    if (
      any(
        !(
          product %in%
            c(
              "LndCov",
              "LndChg",
              "LndCnf",
              "FctImp",
              "ImpDsc",
              "SpcChg"
            )
        )
      )
    ) {
      stop(
        "'product' must be one or more of 'LndCov', 'LndChg', 'LndCnf', 'FctImp', 'ImpDsc', 'SpcChg'."
      )
    }

    if (
      any(
        !(
          region %in%
            c("CU")
        )
      )
    ) {
      stop(
        "'region' currently must be 'CU' (the default)."
      )
    }

    if (
      any(
        !(
          collection %in%
            c(1)
        )
      )
    ) {
      stop(
        "'collection' currently must be '1' (the default)."
      )
    }

    if (
      any(
        !(
          version %in%
            c(1)
        )
      )
    ) {
      stop(
        "'version' currently must be '1' (the default)."
      )
    }

    extraction.dir <-
      normalizePath(extraction.dir, mustWork = FALSE)

    dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)

    template %<>% template_to_sf()

    files <-
      tidyr::expand_grid(
        product =
          product |>
            unique() |>
            factor(
              levels = c(
                "LndCov",
                "LndChg",
                "LndCnf",
                "FctImp",
                "ImpDsc",
                "SpcChg"
              ),
              ordered = TRUE
            ),
        year =
          year |>
            unique() |>
            sort() |>
            as.integer(),
        region = region |>
          unique() |>
          factor(
            levels = c(
              "CU",
              "AK",
              "HI"
            ),
            ordered = TRUE
          ),
        collection = collection |>
          unique() |>
          sort() |>
          as.integer(),
        version = version |>
          unique() |>
          sort() |>
          as.integer()
      ) %>%
      dplyr::mutate(
        file =
          glue::glue("Annual_NLCD_{product}_{year}_{region}_C{collection}V{version}.tif"),
        outfile =
          glue::glue("{extraction.dir}/{label}_{file}")
      )

    # The native grid of the Annual NLCD Collection 1 products:
    # CONUS Albers (EPSG:5070), 30m resolution, with cell edges
    # on multiples of 30m
    nlcd_annual_crs <- "EPSG:5070"
    nlcd_annual_res <- 30
    nlcd_annual_extent <-
      c(
        xmin = -2415600,
        ymin = 164820,
        xmax = 2384400,
        ymax = 3314820
      )

    template %<>%
      sf::st_transform(nlcd_annual_crs)

    # Snap the template bounding box outward to the native NLCD grid
    bbox <- sf::st_bbox(template)
    bbox <-
      c(
        xmin = floor(bbox[["xmin"]] / nlcd_annual_res) * nlcd_annual_res,
        ymin = floor(bbox[["ymin"]] / nlcd_annual_res) * nlcd_annual_res,
        xmax = ceiling(bbox[["xmax"]] / nlcd_annual_res) * nlcd_annual_res,
        ymax = ceiling(bbox[["ymax"]] / nlcd_annual_res) * nlcd_annual_res
      )

    if (
      bbox[["xmin"]] >= nlcd_annual_extent[["xmax"]] ||
        bbox[["xmax"]] <= nlcd_annual_extent[["xmin"]] ||
        bbox[["ymin"]] >= nlcd_annual_extent[["ymax"]] ||
        bbox[["ymax"]] <= nlcd_annual_extent[["ymin"]]
    ) {
      stop("The provided template is not within the specified region.")
    }

    bbox <-
      c(
        xmin = max(bbox[["xmin"]], nlcd_annual_extent[["xmin"]]),
        ymin = max(bbox[["ymin"]], nlcd_annual_extent[["ymin"]]),
        xmax = min(bbox[["xmax"]], nlcd_annual_extent[["xmax"]]),
        ymax = min(bbox[["ymax"]], nlcd_annual_extent[["ymax"]])
      )

    template_nlcd <-
      terra::rast(
        xmin = bbox[["xmin"]],
        xmax = bbox[["xmax"]],
        ymin = bbox[["ymin"]],
        ymax = bbox[["ymax"]],
        resolution = nlcd_annual_res,
        crs = nlcd_annual_crs
      )

    read_nlcd_annual <-
      function(x, prod, yr, outfile) {
        if (!force.redo) {
          if (
            file.exists(outfile)
          ) {
            if (
              compare_rast_dims(terra::rast(outfile), template_nlcd)
            ) {
              return(outfile)
            }
          }
        }

        tmp <- tempfile(fileext = ".tif")

        workspace <- nlcd_annual_wcs_coverages[[as.character(prod)]]

        wcs_url <-
          httr::modify_url(
            "https://dmsdata.cr.usgs.gov/geoserver/",
            path = c("geoserver", workspace, "wcs"),
            query = list(
              service = "WCS",
              version = "1.0.0",
              request = "GetCoverage",
              coverage = paste0(
                workspace, ":", sub("^mrlc_", "", workspace)
              ),
              CRS = nlcd_annual_crs,
              BBOX = paste(
                bbox[["xmin"]], bbox[["ymin"]],
                bbox[["xmax"]], bbox[["ymax"]],
                sep = ","
              ),
              time = paste0(yr, "-01-01T00:00:00.000Z"),
              format = "image/geotiff",
              resx = nlcd_annual_res,
              resy = nlcd_annual_res
            )
          )

        resp <-
          httr::GET(
            wcs_url,
            httr::write_disk(tmp, overwrite = TRUE)
          )

        if (
          httr::http_error(resp) ||
            !grepl("tiff", httr::headers(resp)[["content-type"]])
        ) {
          stop(
            "Download of ", x, " from the Annual NLCD web coverage service at ",
            "https://dmsdata.cr.usgs.gov/geoserver/", workspace, "/wcs failed! ",
            "You may have better luck downloading the data manually from ",
            "https://www.mrlc.gov/data."
          )
        }

        suppressWarnings(
          out <-
            terra::rast(tmp)
        )

        # The web coverage service introduces sub-millimeter floating-point
        # jitter in the returned geotransform; snap back to the native grid
        if (
          identical(terra::ncol(out), terra::ncol(template_nlcd)) &&
            identical(terra::nrow(out), terra::nrow(template_nlcd))
        ) {
          terra::ext(out) <- terra::ext(template_nlcd)
        } else {
          out %<>%
            terra::resample(template_nlcd, method = "near")
        }

        if (stringr::str_detect(x, "LndCov")) {
          levels(out) <-
            nlcd_colors() %>%
            as.data.frame()

          terra::coltab(out) <-
            nlcd_colors() |>
            dplyr::select(ID, Color) |>
            as.data.frame()

          terra::writeRaster(
            x = out,
            filename = outfile,
            datatype = "INT1U",
            gdal = raster.options,
            overwrite = TRUE
          )
        }

        if (stringr::str_detect(x, "LndChg")) {
          terra::writeRaster(
            x = terra::as.factor(out),
            filename = outfile,
            datatype = "INT2U",
            gdal = raster.options,
            overwrite = TRUE
          )
        }

        if (stringr::str_detect(x, "LndCnf|FctImp")) {
          terra::writeRaster(
            x = out,
            filename = outfile,
            datatype = "INT1U",
            gdal = raster.options,
            overwrite = TRUE
          )
        }

        if (stringr::str_detect(x, "ImpDsc")) {
          impdsc_coltab <- terra::coltab(out)[[1]]

          levels(out) <-
            tibble::tibble(
              ID = 0:2,
              Class = c(
                "Non-Urban",
                "Roads",
                "Urban"
              )
            ) %>%
            as.data.frame()

          if (is.null(impdsc_coltab)) {
            impdsc_coltab <-
              tibble::tibble(
                value = 0:2,
                red = c(0L, 230L, 38L),
                green = c(0L, 0L, 115L),
                blue = c(0L, 0L, 0L),
                alpha = c(0L, 255L, 255L)
              )
          }

          terra::coltab(out) <-
            impdsc_coltab |>
            dplyr::filter(value %in% 0:2) |>
            as.data.frame()

          terra::writeRaster(
            x = out,
            filename = outfile,
            datatype = "INT1U",
            gdal = raster.options,
            overwrite = TRUE
          )
        }

        if (stringr::str_detect(x, "SpcChg")) {
          terra::writeRaster(
            x = out,
            filename = outfile,
            datatype = "INT2U",
            gdal = raster.options,
            overwrite = TRUE
          )
        }

        return(outfile)
      }

    return(
      files |>
        dplyr::rowwise() |>
        dplyr::mutate(
          rast = read_nlcd_annual(file, product, year, outfile) |>
            terra::rast() |>
            list()
        ) |>
        dplyr::select(!file)
    )
  }

# The USGS web coverage service workspaces serving each
# Annual NLCD product in its native CRS and resolution
nlcd_annual_wcs_coverages <-
  c(
    LndCov = "mrlc_Land-Cover-Native_conus_year_data",
    LndChg = "mrlc_Land-Cover-Change-Native_conus_year_data",
    LndCnf = "mrlc_Land-Cover-Confidence-Native_conus_year_data",
    FctImp = "mrlc_Fractional-Impervious-Surface-Native_conus_year_data",
    ImpDsc = "mrlc_Impervious-Descriptor-Native_conus_year_data",
    SpcChg = "mrlc_Spectral-Change-Day-of-Year-Native_conus_year_data"
  )

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

list_nlcd <-
  function() {
    (httr::GET(
      "https://www.mrlc.gov/geoserver/ows",
      query = list(
        service = "WCS",
        acceptversions = "2.0.1",
        request = "GetCapabilities"
      )
    ) |>
      httr::content() |>
      xml2::as_list())$Capabilities$Contents |>
      purrr::map(\(x)
      list(
        Title = x$Title[[1]],
        id = x$CoverageId[[1]]
      )) |>
      magrittr::set_names(NULL) |>
      purrr::transpose() |>
      purrr::map(unlist) |>
      tibble::as_tibble() |>
      dplyr::filter(
        stringr::str_starts(id, "mrlc_download"),
        stringr::str_detect(id, "NLCD|nlcd"),
        stringr::str_detect(id, "Science_Product", negate = TRUE),
        stringr::str_detect(id, "Disturbance_Date", negate = TRUE),
        stringr::str_detect(id, "Fractional_Component", negate = TRUE),
        stringr::str_detect(id, "descriptor", negate = TRUE),
        stringr::str_detect(id, "Descriptor", negate = TRUE),
        stringr::str_detect(id, "Pixels", negate = TRUE),
        stringr::str_detect(id, "Annual_NLCD", negate = TRUE),
        stringr::str_detect(id, "Count", negate = TRUE),
        stringr::str_detect(id, "Index", negate = TRUE)
      ) |>
      dplyr::mutate(
        Title = stringr::str_replace(
          Title,
          "nlcd_tcc_conus",
          "NLCD"
        ),
        Title = stringr::str_replace(
          Title,
          "v2021-4",
          "Tree_Canopy_L48"
        ),
        Title = stringr::str_replace(
          Title,
          "Land_Cover",
          "landcover"
        ),
        Title = stringr::str_replace(
          Title,
          "Impervious",
          "impervious"
        ),
        Title = stringr::str_replace(
          Title,
          "Tree_Canopy",
          "canopy"
        )
      ) |>
      tidyr::separate_wider_delim(
        cols = Title,
        delim = "_",
        names = c(
          "NLCD",
          "year",
          "dataset",
          "landmass"
        )
      ) |>
      dplyr::select(!NLCD) |>
      dplyr::arrange(landmass, dataset, year) |>
      print(n = 200)
  }

