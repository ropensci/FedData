#' Download and crop data from the NRCS SSURGO soils database.
#'
#' This is an efficient method for spatially merging several different soil survey areas
#' as well as merging their tabular data.
#'
#' \code{get_ssurgo} returns a named list of length 2:
#' \enumerate{
#' \item 'spatial': A [`Simple Feature`][sf::sf] of soil mapunits
#' in the template, and
#' \item 'tabular': A named list of \code{\link{data.frame}s} with the SSURGO tabular data.
#' }
#'
#' @param template An [`Simple Feature`][sf::sf]
#' or [`SpatRaster`][terra::SpatRaster] object to serve as a template for cropping.
#' Optionally, a vector of area names, e.g., `c('IN087','IN088')` may be provided.
#' @param label A character string naming the study area.
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing. Defaults to './RAW/SSURGO/'.
#' @param extraction.dir A character string indicating where the extracted and cropped SSURGO shapefiles should be put.
#' The directory will be created if missing. Defaults to './EXTRACTIONS/SSURGO/'.
#' @param force.redo If an extraction for this template and label already exists, should a new one be created? Defaults to FALSE.
#' @return A named list containing the 'spatial' and 'tabular' data.
#' @export
#' @importFrom readr read_csv write_csv
#' @examples
#' \donttest{
#' # Get the NRCS SSURGO data (USA ONLY)
#' SSURGO.MEVE <-
#'   get_ssurgo(
#'     template = FedData::meve,
#'     label = "meve"
#'   )
#'
#' # Plot the VEP polygon
#' plot(meve)
#'
#' # Plot the SSURGO mapunit polygons
#' plot(SSURGO.MEVE$spatial["MUKEY"],
#'   lwd = 0.1,
#'   add = TRUE
#' )
#'
#' # Or, download by Soil Survey Area names
#' SSURGO.areas <-
#'   get_ssurgo(
#'     template = c("CO670", "CO075"),
#'     label = "CO_TEST"
#'   )
#'
#' # Let's just look at spatial data for CO675
#' SSURGO.areas.CO675 <-
#'   SSURGO.areas$spatial[SSURGO.areas$spatial$AREASYMBOL == "CO075", ]
#'
#' # And get the NED data under them for pretty plotting
#' NED.CO675 <-
#'   get_ned(
#'     template = SSURGO.areas.CO675,
#'     label = "SSURGO_CO675"
#'   )
#'
#' # Plot the SSURGO mapunit polygons, but only for CO675
#' terra::plot(NED.CO675)
#' plot(
#'   SSURGO.areas.CO675$geom,
#'   lwd = 0.1,
#'   add = TRUE
#' )
#' }
get_ssurgo <- function(template,
                       label,
                       raw.dir = paste0(tempdir(), "/FedData/raw/ssurgo"),
                       extraction.dir = paste0(tempdir(), "/FedData/"),
                       force.redo = FALSE) {
  raw.dir <- normalizePath(paste0(raw.dir, "/."), mustWork = FALSE)
  extraction.dir <- normalizePath(paste0(extraction.dir, "/."), mustWork = FALSE)

  dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)

  outfile <- paste0(extraction.dir, "/", label, "_ssurgo.gpkg")

  if (!force.redo & file.exists(outfile)) {
    SSURGOData <-
      outfile %>%
      sf::st_layers() %$%
      name %>%
      magrittr::set_names(., .) %>%
      purrr::map(~ sf::read_sf(outfile, layer = .x))

    return(list(
      spatial = SSURGOData$geometry,
      tabular = purrr::list_modify(SSURGOData, geometry = NULL)
    ))
  }

  if (inherits(template, "character")) {
    q <- paste0(
      "SELECT areasymbol, saverest FROM sacatalog WHERE areasymbol IN (", paste(paste0("'", template, "'"), collapse = ","),
      ");"
    )
    SSURGOAreas <- soils_query(q)
  } else {
    template %<>% template_to_sf()

    # Get shapefile of SSURGO study areas in the template
    SSURGOAreas <-
      get_ssurgo_inventory(template = template, raw.dir = raw.dir)

    if (!any(SSURGOAreas$iscomplete == 1)) {
      stop("There are no complete soil surveys in your study area.")
    }


    # Remove SSURGO study areas that are not available
    SSURGOAreas <- SSURGOAreas[SSURGOAreas$iscomplete != 0, ]
  }

  # Get data for each study area
  SSURGOData <-
    purrr::map2(
      .x = SSURGOAreas$areasymbol,
      .y = SSURGOAreas$saverest,
      .f = function(area, date) {
        message("(Down)Loading SSURGO data for survey area ", as.character(area))
        get_ssurgo_study_area(
          template = template,
          area = as.character(area),
          date = as.Date(date, format = "%m/%d/%Y"),
          raw.dir = raw.dir
        )
      }
    ) %>%
    purrr::transpose()


  # Combine mapunits
  SSURGOData$spatial %<>%
    do.call("rbind", .)

  # Crop to template
  if (!identical(class(template), "character")) {
    suppressWarnings(
      SSURGOData$spatial %<>%
        dplyr::ungroup() %>%
        sf::st_cast("MULTIPOLYGON") %>%
        sf::st_cast("POLYGON")
    )

    suppressMessages(
      suppressWarnings(
        SSURGOData$spatial %<>%
          dplyr::filter(
            SSURGOData$spatial %>%
              sf::st_intersects(
                template %>%
                  sf::st_transform(
                    sf::st_crs(SSURGOData$spatial)
                  )
              ) %>%
              purrr::map_lgl(~ (length(.x) > 0))
          )
      )
    )
  }

  SSURGOData$tabular %<>%
    purrr::transpose() %>%
    purrr::map(function(x) {
      x %>%
        purrr::compact() %>%
        purrr::discard(~ nrow(.x) == 0) %>%
        purrr::map(function(y) {
          if (any(names(y) == "musym")) {
            y %<>% dplyr::mutate(musym = as.character(musym))
          }
          y
        }) %>%
        dplyr::bind_rows() %>%
        readr::type_convert(col_types = readr::cols())
    })

  # Extract only the mapunits in the study area,
  # and iterate through the data structure
  SSURGOData$tabular %<>%
    extract_ssurgo_data(as.character(SSURGOData$spatial$MUKEY))

  # SSURGOData$tabular$mapunit %<>%
  #   dplyr::mutate(mukey = as.character(mukey)) %>%
  #   dplyr::left_join(SSURGOData$spatial %>%
  #                      dplyr::select(MUKEY),
  #                    by = c("mukey" = "MUKEY")) %>%
  #   sf::st_as_sf()

  unlink(outfile)

  SSURGOData %$%
    c(geometry = list(spatial), tabular) %>%
    purrr::iwalk(function(x, n) {
      sf::write_sf(x, dsn = outfile, layer = n)
    })

  SSURGOData <-
    outfile %>%
    sf::st_layers() %$%
    name %>%
    magrittr::set_names(., .) %>%
    purrr::map(~ sf::read_sf(outfile, layer = .x))

  return(list(spatial = SSURGOData$geometry, tabular = purrr::list_modify(SSURGOData, geometry = NULL)))
}

#' Download a zipped directory containing a shapefile of the SSURGO study areas.
#'
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' @return A character string representing the full local path of the SSURGO study areas zipped directory.
#' @export
#' @keywords internal
download_ssurgo_inventory <- function(raw.dir, ...) {
  # Import the shapefile of SSURGO study areas.  This is available at
  # http://soildatamart.sc.egov.usda.gov/download/StatusMaps/soilsa_a_SSURGO.zip
  url <- "https://websoilsurvey.sc.egov.usda.gov/DataAvailability/SoilDataAvailabilityShapefile.zip"
  destdir <- raw.dir
  download_data(url = url, destdir = destdir, ...)
  return(normalizePath(paste(destdir, "/SoilDataAvailabilityShapefile.zip", sep = "")))
}

#' Download and crop a shapefile of the SSURGO study areas.
#'
#' \code{get_ssurgo_inventory} returns a \code{SpatialPolygonsDataFrame} of the SSURGO study areas within
#' the specified \code{template}. If template is not provided, returns the entire SSURGO inventory of study areas.
#'
#' @param template An [`Simple Feature`][sf::sf]
#' or [`SpatRaster`][terra::SpatRaster] object to serve as a template for cropping.
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing.
#' @return A \code{SpatialPolygonsDataFrame} of the SSURGO study areas within
#' the specified \code{template}.
#' @export
#' @keywords internal
get_ssurgo_inventory <- function(template = NULL, raw.dir) {
  if (!is.null(template)) {
    template %<>%
      template_to_sf() %>%
      sf::st_transform(4326)
  }

  # If there is a template, only download the areas in the template. Thanks to Dylan Beaudette for this method!
  if (
    !is.null(template) &&
      httr::status_code(
        httr::RETRY(
          verb = "GET",
          url = "https://sdmdataaccess.nrcs.usda.gov/Spatial/SDMWGS84Geographic.wfs"
        )
      ) == 200
  ) {
    bounds <-
      template %>%
      sf::st_bbox() %>%
      sf::st_as_sfc()

    # Only download 1 square degree at a time to avoid oversized AOI error
    if ((sf::st_bbox(template)[["xmax"]] - sf::st_bbox(template)[["xmin"]]) > 1 |
      (sf::st_bbox(template)[["ymax"]] - sf::st_bbox(template)[["ymin"]]) > 1) {
      bounds %<>%
        sf::st_intersection(grid)
    }

    SSURGOAreas <-
      bounds %>%
      purrr::map_dfr(function(x) {
        bound <-
          x %>%
          sf::st_bbox()

        if (identical(bound["xmin"], bound["xmax"])) {
          bound["xmax"] <- bound["xmax"] + 1e-04
        }
        if (identical(bound["ymin"], bound["ymax"])) {
          bound["ymax"] <- bound["ymax"] + 1e-04
        }

        bbox.text <- paste(bound, collapse = ",")

        temp.file <- paste0(tempdir(), "/soils.gml")

        httr::RETRY(
          verb = "GET",
          url = "https://sdmdataaccess.nrcs.usda.gov/Spatial/SDMWGS84Geographic.wfs",
          query = list(
            Service = "WFS",
            Version = "1.1.0",
            Request = "GetFeature",
            Typename = "SurveyAreaPoly",
            BBOX = bbox.text,
            SRSNAME = "EPSG:4326",
            OUTPUTFORMAT = "GML3"
          ),
          httr::write_disk(temp.file,
            overwrite = TRUE
          )
        )

        tryCatch(
          suppressMessages(
            suppressWarnings(
              sf::read_sf(temp.file, drivers = "GML") %>%
                dplyr::mutate(saverest = as.Date(lubridate::parse_date_time(saverest, orders = "b d Y HMOp", locale = "en_US"))) %>%
                # sf::st_intersection(template) %>%
                sf::st_drop_geometry()
            )
          ),
          error = function(e) {
            return(NULL)
          }
        )
      }) %>%
      dplyr::distinct() %>%
      dplyr::arrange(areasymbol)
  } else {
    tmpdir <- tempfile()
    if (!dir.create(tmpdir)) {
      stop("failed to create my temporary directory")
    }

    file <- download_ssurgo_inventory(raw.dir = raw.dir)

    utils::unzip(file, exdir = tmpdir)

    SSURGOAreas <- sf::read_sf(normalizePath(tmpdir), layer = "soilsa_a_nrcs")

    if (!is.null(template)) {
      SSURGOAreas %<>%
        sf::st_make_valid() %>%
        sf::st_intersection(sf::st_transform(template, sf::st_crs(SSURGOAreas)))
    }

    unlink(tmpdir, recursive = TRUE)
  }

  # Check to see if all survey areas are available
  if (0 %in% SSURGOAreas$iscomplete) {
    warning(
      "Some of the soil surveys in your area are unavailable.\n
            Soils and productivity data will have holes.\n
            Missing areas:\n",
      paste0(as.vector(SSURGOAreas[SSURGOAreas$iscomplete == 0, ]$areasymbol), collapse = "\n"), "\n\n
            Continuing with processing available soils.\n\n"
    )
  }

  return(SSURGOAreas)
}

#' Download a zipped directory containing the spatial and tabular data for a SSURGO study area.
#'
#' \code{download_ssurgo_study_area} first tries to download data including a state-specific Access
#' template, then the general US template.
#'
#' @param area A character string indicating the SSURGO study area to be downloaded.
#' @param date A character string indicating the date of the most recent update to the SSURGO
#' area for these data. This information may be gleaned from the SSURGO Inventory (\code{\link{get_ssurgo_inventory}}).
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' @return A character string representing the full local path of the SSURGO study areas zipped directory.
#' @export
#' @keywords internal
download_ssurgo_study_area <- function(area, date, raw.dir) {
  # Try to download with the state database, otherwise grab the US
  url <- paste("http://websoilsurvey.sc.egov.usda.gov/DSD/Download/Cache/SSA/wss_SSA_", area, "_[", date, "].zip", sep = "")
  destdir <- raw.dir
  download_data(url = url, destdir = destdir, nc = TRUE)

  return(normalizePath(paste(destdir, "/wss_SSA_", area, "_[", date, "].zip", sep = "")))
}

#' Download and crop the spatial and tabular data for a SSURGO study area.
#'
#' \code{get_ssurgo_study_area} returns a named list of length 2:
#' \enumerate{
#' \item 'spatial': A [`Simple Feature`][sf::sf] of soil mapunits
#' in the template, and
#' \item 'tabular': A named list of \code{\link{data.frame}s} with the SSURGO tabular data.
#' }
#'
#' @param template An [`Simple Feature`][sf::sf]
#' or [`SpatRaster`][terra::SpatRaster] object to serve as a template for cropping.
#' If missing, whose study area is returned
#' @param area A character string indicating the SSURGO study area to be downloaded.
#' @param date A character string indicating the date of the most recent update to the SSURGO
#' area for these data. This information may be gleaned from the SSURGO Inventory (\code{\link{get_ssurgo_inventory}}).
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing.
#' @return A \code{SpatialPolygonsDataFrame} of the SSURGO study areas within
#' the specified \code{template}.
#' @export
#' @keywords internal
get_ssurgo_study_area <- function(template = NULL, area, date, raw.dir) {
  tmpdir <- tempfile()
  if (!dir.create(tmpdir)) {
    stop("failed to create my temporary directory")
  }

  file <- download_ssurgo_study_area(area = area, date = date, raw.dir = raw.dir)

  utils::unzip(file, exdir = tmpdir)
  suppressMessages({
    mapunits <-
      sf::read_sf(paste0(tmpdir, "/", area, "/spatial"),
        layer = paste0("soilmu_a_", tolower(area))
      ) %>%
      sf::st_make_valid()
  })

  # Read in all tables
  tablesData <-
    paste0(tmpdir, "/", area, "/tabular") %>%
    list.files(full.names = TRUE) %>%
    magrittr::set_names(
      .,
      tools::file_path_sans_ext(basename(.))
    ) %>%
    purrr::map(
      function(file) {
        # Hack to bypass one line file bug in readr::read_delim
        if (length(readLines(file)) == 1) {
          write("\n", file, append = TRUE)
        }

        tryCatch(
          return(
            readr::read_delim(file,
              col_names = FALSE,
              col_types = readr::cols(.default = readr::col_character()),
              delim = "|"
            )
          ),
          error = function(e) {
            return(NULL)
          }
        )
      }
    ) %>%
    purrr::compact()


  SSURGOTableMapping <-
    tablesData$mstab %>%
    dplyr::select(1, 5) %>%
    magrittr::set_names(c("TABLE", "FILE"))

  tablesData <- tablesData[SSURGOTableMapping$FILE]
  tablesHeads <- tablesHeaders[SSURGOTableMapping$TABLE]

  notNull <- (!sapply(tablesData, is.null) & !sapply(tablesHeads, is.null))
  tablesData <- tablesData[notNull]
  tablesHeads <- tablesHeads[notNull]

  tables <-
    mapply(tablesData,
      tablesHeads,
      FUN = function(theData, theHeader) {
        names(theData) <- names(theHeader)
        return(theData)
      }
    )

  names(tables) <- names(tablesHeads)

  tables %<>%
    extract_ssurgo_data(as.character(mapunits$MUKEY))

  unlink(tmpdir, recursive = TRUE)

  return(list(spatial = mapunits, tabular = tables))
}

#' Extract data from a SSURGO database pertaining to a set of mapunits.
#'
#' \code{extract_ssurgo_data} creates a directed graph of the joins in a SSURGO tabular dataset,
#' and then iterates through the tables, only retaining data pertinent to a set of mapunits.
#'
#' @param tables A list of SSURGO tabular data.
#' @param mapunits A character vector of mapunits (likely dropped from SSURGO spatial data)
#' defining which mapunits to retain.
#' @return A list of extracted SSURGO tabular data.
#' @export
#' @keywords internal
extract_ssurgo_data <- function(tables, mapunits) {
  mapping <- tables$mdstatrshipdet
  mappingGraph <- igraph::graph.edgelist(as.matrix(mapping[, c("ltabphyname", "rtabphyname")]))
  igraph::E(mappingGraph)$mapVar <- as.character(mapping$ltabcolphyname)

  mappingGraph <- igraph::graph.neighborhood(mappingGraph, order = max(sapply(igraph::decompose.graph(mappingGraph), igraph::diameter)) +
    1, nodes = "mapunit", mode = "out")[[1]]
  mapHierarchy <- igraph::shortest.paths(mappingGraph, "mapunit")
  mapHierarchy <- colnames(mapHierarchy)[order(mapHierarchy)]
  mapHierarchy <- mapHierarchy[-1]
  mapEdges <- cbind(igraph::get.edgelist(mappingGraph), igraph::E(mappingGraph)$mapVar)
  mapEdges <- mapEdges[match(mapHierarchy, mapEdges[, 2]), ]

  tables$mapunit %<>%
    dplyr::filter(mukey %in% mapunits)

  for (i in 1:nrow(mapEdges)) {
    X <- mapEdges[i, ]
    tables[[X[2]]] <- tables[[X[2]]][tables[[X[2]]][[X[3]]] %in% tables[[X[1]]][[X[3]]], ]
  }

  return(tables)
}

#' Submit a Soil Data Access (SDA) Query
#'
#' \code{soils_query} submit an SQL query to retrieve data from the Soil Data Mart.
#' Please see https://sdmdataaccess.sc.egov.usda.gov/Query.aspx for guidelines
#'
#' @param q A character string representing a SQL query to the SDA service
#' @return A tibble returned from the SDA service
#' @keywords internal
soils_query <- function(q) {
  tryCatch(
    httr::POST(
      url = "https://sdmdataaccess.sc.egov.usda.gov/tabular/post.rest",
      body = list(
        query = q,
        format = "JSON+COLUMNNAME"
      ),
      encode = "form"
    ) %>%
      httr::content(as = "parse", encoding = "UTF-8") %$%
      Table %>%
      purrr::transpose() %>%
      magrittr::set_names(., purrr::map_chr(., ~ (.x[[1]]))) %>%
      purrr::map(~ (unlist(.x[-1]))) %>%
      tibble::as_tibble(),
    error = function(e) {
      stop("Improperly formatted SDA SQL request")
    }
  )
}
