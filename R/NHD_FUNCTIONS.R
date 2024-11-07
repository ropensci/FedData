#' Download and crop the National Hydrography Dataset.
#'
#' \code{get_nhd} returns a list of [`Simple Feature`][sf::sf] objects extracted
#' from the National Hydrography Dataset.
#'
#' @param template An [`Simple Feature`][sf::sf]
#' or [`SpatRaster`][terra::SpatRaster] object to serve as a template for cropping.
#' @param label A character string naming the study area.
#' @param nhdplus Extract data from the USGS NHDPlus High Resolution service (experimental)
#' @param extraction.dir A character string indicating where the extracted and cropped NHD data should be put.
#' @param force.redo If an extraction for this template and label already exists, should a new one be created?
#' @return A list of `sf` collections extracted from the National Hydrography Dataset.
#' @importFrom magrittr %>% %<>%
#' @export
#' @examples
#' \donttest{
#' # Get the NHD (USA ONLY)
#' NHD <- get_nhd(
#'   template = FedData::meve,
#'   label = "meve"
#' )
#' NHD
#' NHD %>%
#'   plot_nhd(template = FedData::meve)
#' }
get_nhd <-
  function(template,
           label,
           nhdplus = FALSE,
           extraction.dir = file.path(
             tempdir(),
             "FedData",
             "extractions",
             "nhd",
             label
           ),
           force.redo = FALSE) {
    extraction.dir <-
      normalizePath(extraction.dir,
        mustWork = FALSE
      )

    dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)
    out_file <- paste0(label, "_nhd.gpkg")

    out_dsn <- file.path(extraction.dir, out_file)

    if (!force.redo & file.exists(out_dsn)) {
      return(read_sf_all(out_dsn))
    }

    template %<>%
      template_to_sf() %>%
      sf::st_as_sfc() %>%
      sf::st_union() %>%
      sf::st_cast("POLYGON")

    if (nhdplus) {
      layers <-
        c(
          "NHDPoint",
          "NetworkNHDFlowline",
          "NonNetworkNHDFlowline",
          "NHDLine",
          "NHDArea",
          "NHDWaterbody"
        )

      nhd_out <-
        "https://hydro.nationalmap.gov/arcgis/rest/services/NHDPlus_HR/MapServer" %>%
        arcgislayers::arc_open() %>%
        arcgislayers::get_layers(
          name = layers
        ) %>%
        purrr::map(
          ~ tryCatch(
            arcgislayers::arc_select(
              .x,
              filter_geom =
                template
            ),
            error = function(e) {
              NULL
            }
          )
        ) %>%
        magrittr::set_names(layers) %$%
        list(
          Point = NHDPoint,
          Flowline = list(
            NetworkNHDFlowline,
            NonNetworkNHDFlowline
          ) %>%
            dplyr::bind_rows(),
          Line = NHDLine,
          Area = NHDArea,
          Waterbody = NHDWaterbody
        )
    } else {
      layers <-
        c(
          "Point",
          "Flowline - Large Scale",
          "Line - Large Scale ",
          "Area - Large Scale",
          "Waterbody - Large Scale"
        )

      nhd_out <-
        "https://hydro.nationalmap.gov/arcgis/rest/services/nhd/MapServer" %>%
        arcgislayers::arc_open() %>%
        arcgislayers::get_layers(
          name = layers
        ) %>%
        purrr::map(
          ~ tryCatch(
            arcgislayers::arc_select(
              .x,
              filter_geom =
                template
            ),
            error = function(e) {
              NULL
            }
          )
        ) %>%
        magrittr::set_names(layers) %$%
        list(
          Point = Point,
          Flowline = `Flowline - Large Scale`,
          Line = `Line - Large Scale `,
          Area = `Area - Large Scale`,
          Waterbody = `Waterbody - Large Scale`
        )
    }

    null_elements <- purrr::map_lgl(nhd_out, is.null)
    if (all(null_elements)) stop("No NHD data present within template.")
    nhd_for_crs <- nhd_out[!null_elements][[1]]

    suppressWarnings({
      suppressMessages({
        nhd_out %>%
          purrr::compact() %>%
          purrr::map(sf::st_make_valid) %>%
          purrr::map(
            sf::st_intersection,
            template %>%
              sf::st_geometry() %>%
              sf::st_transform(sf::st_crs(nhd_for_crs))
          ) %>%
          write_sf_all(dsn = out_dsn)
      })
    })

    return(read_sf_all(out_dsn))
  }


#' A basic plotting function for NHD data.
#'
#' This is more of an example than anything
#'
#' @param x The result of [get_nhd].
#' @param template An [`Simple Feature`][sf::sf]
#' or [`SpatRaster`][terra::SpatRaster] object to serve as a template for cropping.
#' @return A `ggplot2` panel of plots
#' @importFrom magrittr %>% %<>%
#' @export
#' @examples
#' \donttest{
#' # Get the NHD (USA ONLY)
#' NHD <- get_nhd(
#'   template = FedData::meve,
#'   label = "meve"
#' )
#' NHD
#' NHD %>%
#'   plot_nhd(template = FedData::meve)
#' }
plot_nhd <-
  function(x,
           template = NULL) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package \"ggplot2\" needed for this function to work. Please install it.",
        call. = FALSE
      )
    }

    template %<>%
      template_to_sf() %>%
      sf::st_transform(sf::st_crs(x[[1]]))

    g <-
      x %>%
      purrr::map(~ dplyr::select(.x, geom)) %>%
      purrr::imap(~ dplyr::mutate(.x, Dataset = .y)) %>%
      do.call("rbind", .) %>%
      ggplot2::ggplot() +
      ggplot2::geom_sf() +
      ggplot2::facet_wrap(
        facets = "Dataset",
        ncol = 2
      ) +
      ggplot2::theme_minimal()

    if (!is.null(template)) {
      g <-
        g +
        ggplot2::geom_sf(
          data = template,
          color = "red",
          fill = NA
        )
    }

    g
  }


#' Download and crop the Watershed Boundary Dataset.
#'
#' \code{get_wbd} returns an [`Simple Feature`][sf::sf] collection of the HUC 12 regions within
#' the specified \code{template}.
#'
#' @param template An [`Simple Feature`][sf::sf]
#' or [`SpatRaster`][terra::SpatRaster] object to serve as a template for cropping.
#' @param label A character string naming the study area.
#' @param extraction.dir A character string indicating where the extracted and cropped NHD data should be put.
#' @param force.redo If an extraction for this template and label already exists, should a new one be created?
#' @return An `sf` collection of the HUC 12 regions within
#' the specified \code{template}.
#' @export
get_wbd <- function(template,
                    label,
                    extraction.dir = file.path(
                      tempdir(),
                      "FedData",
                      "extractions",
                      "nhd",
                      label
                    ),
                    force.redo = FALSE) {
  extraction.dir <-
    normalizePath(extraction.dir,
      mustWork = FALSE
    )

  dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)

  out_file <- paste0(label, "_nhd_wbd.gpkg")
  out_dsn <- file.path(extraction.dir, out_file)

  if (!force.redo & file.exists(out_dsn)) {
    return(read_sf_all(out_dsn))
  }

  "https://hydro.nationalmap.gov/arcgis/rest/services/NHDPlus_HR/MapServer/" %>%
    arcgislayers::arc_open() %>%
    arcgislayers::get_layer(
      name = "WBDHU12"
    ) %>%
    arcgislayers::arc_select(
      filter_geom =
        template %>%
          template_to_sf() %>%
          sf::st_as_sfc()
    ) %>%
    sf::write_sf(dsn = out_dsn)

  return(
    sf::read_sf(out_dsn)
  )
}
