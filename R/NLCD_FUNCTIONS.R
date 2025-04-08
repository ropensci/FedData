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

  source <-
    httr::modify_url(
      nlcd_base_url,
      query =
        list(
          version = "2.0.1",
          coverageid = coverage
        )
    )

  out <-
    tryCatch(
      paste0("WCS:", source) |>
        terra::rast(),
      error = function(e) {
        stop(
          "Web coverage service at ",
          source, " is corrupted. You may have better luck downloading the file at https://www.mrlc.gov/data."
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

#' Download and crop the Annual National Land Cover Database.
#'
#' \code{get_nlcd_annual} returns a [`SpatRaster`][terra::SpatRaster] of NLCD data cropped to a given
#' template study area. The Annual NLCD is currently only available for the conterminous United States.
#' More information about the Annual NLCD product is available on the
#' [Annual NLCD web page](https://www.mrlc.gov/data/project/annual-nlcd).
#'
#' @param template An [`Simple Feature`][sf::sf]
#' or [`terra`][terra::SpatRaster] object to serve as a template for cropping.
#' @param label A character string naming the study area.
#' @param year An integer vector representing the year of desired NLCD product.
#' Acceptable values are currently 1985 through 2023 (defaults to 2023).
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
#' **Currently, only '0' is available.**
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
           year = 2023,
           product = "LndCov",
           region = "CU",
           collection = 1,
           version = 0,
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
    nlcd_annual_bucket <- "https://s3-us-west-2.amazonaws.com/mrlc"

    if (any(!(year %in% 1985:2023))) {
      stop("'year' must only contain integers between 1985 and 2023.")
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
            c(0)
        )
      )
    ) {
      stop(
        "'version' currently must be '0' (the default)."
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

    suppressWarnings(
      template_nlcd <-
        file.path(
          "/vsicurl",
          nlcd_annual_bucket,
          files$file[[1]]
        ) |>
        terra::rast() |>
        terra::rast()
    )

    template %<>%
      sf::st_transform(
        terra::crs(template_nlcd)
      )

    tryCatch(
      template_nlcd %<>%
        terra::crop(template,
          snap = "out"
        ),
      error =
        function(e) {
          stop("The provided template is not within the specified region.")
        }
    )

    read_nlcd_annual <-
      function(x, outfile) {
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

        suppressWarnings(
          out <-
            terra::rast(
              file.path(
                "/vsicurl",
                nlcd_annual_bucket,
                x
              )
            ) |>
            terra::crop(template_nlcd)
        )

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

          terra::coltab(out) <-
            terra::coltab(out)[[1]] |>
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
          rast = read_nlcd_annual(file, outfile) |>
            terra::rast() |>
            list()
        ) |>
        dplyr::select(!file)
    )
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
