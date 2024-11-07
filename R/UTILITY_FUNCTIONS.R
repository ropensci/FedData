# Make CRAN check not complain about "." and package data
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".",
    "element",
    "nlcd_tiles",
    "nlcd_landcover_pam",
    "nlcd_canopy_pam",
    "nlcd_impervious_pam",
    "daymet_tiles",
    "NewDataSet",
    "Table",
    "saverest",
    "areasymbol",
    "tablesHeaders",
    "xmin",
    "xmax",
    "ymin",
    "ymax",
    "xsize",
    "ExceptionReport",
    "name",
    "ServiceExceptionReport",
    "year",
    "AREASYMBOL",
    "Area - Large Scale",
    "Flowline - Large Scale",
    "Line - Large Scale",
    "MUKEY",
    "MUSYM",
    "NHDArea",
    "NHDLine",
    "NHDPoint",
    "NHDWaterbody",
    "Point",
    "SPATIALVER",
    "Waterbody - Large Scale",
    "Line - Large Scale ",
    "NetworkNHDFlowline",
    "NonNetworkNHDFlowline",
    "mukey",
    "musym",
    "spatial",
    "tabular",
    "CoverageDescriptions",
    "CoverageDescription",
    "value",
    "Land Cover",
    "ID",
    "YEAR",
    "ELEMENT",
    "MONTH",
    "NAME",
    "properties",
    "objectIds",
    "Color",
    "cats",
    "STATION",
    "DAY",
    "DATE",
    "outfile"
  ))
}

#' Get the rightmost 'n' characters of a character string.
#'
#' @param x A character string.
#' @param n The number of characters to retrieve.
#' @return A character string.
#' @export
#' @keywords internal
substr_right <- function(x, n) {
  substr(x, nchar(x) - n + 1, nchar(x))
}

#' Turn an extent object into a polygon
#'
#' @param x An object from which an bounding box object can be retrieved.
#' @param proj4string A PROJ.4 formatted string defining the required projection.
#' @return A [`Simple Feature`][sf::sf] object.
#' @export
#' @keywords internal
polygon_from_extent <- function(x, proj4string = NULL) {
  x %<>%
    template_to_sf() %>%
    sf::st_bbox() %>%
    sf::st_as_sfc()

  if (!is.null(proj4string)) {
    x %<>%
      sf::st_transform(proj4string)
  }
  return(x)
}

template_to_sf <-
  function(template) {
    if (inherits(template, c(
      "RasterLayer",
      "RasterStack",
      "RasterBrick",
      "Extent",
      "SpatRaster",
      "SpatVector"
    ))) {
      template %<>%
        sf::st_bbox() %>%
        sf::st_as_sfc()
    }

    template %<>%
      sf::st_as_sf()

    return(template)
  }

read_sf_all <- function(dsn) {
  dsn %>%
    sf::st_layers() %$%
    name %>%
    magrittr::set_names(., .) %>%
    purrr::map(~ sf::read_sf(
      dsn = dsn,
      layer = .x
    ))
}

write_sf_all <-
  function(x, dsn) {
    if (is.null(names(x))) {
      stop("'x' must be a named list.")
    }

    unlink(dsn,
      recursive = TRUE,
      force = TRUE
    )

    x %>%
      purrr::iwalk(
        ~ sf::write_sf(.x,
          dsn = dsn,
          layer = .y,
          delete_layer = TRUE
        )
      )
  }

#' Get a logical vector of which elements in a vector are sequentially duplicated.
#'
#' @param x An vector of any type, or, if \code{rows}, a matrix.
#' @param rows Is x a matrix?
#' @return A logical vector of the same length as x.
#' @export
#' @keywords internal
sequential_duplicated <- function(x, rows = FALSE) {
  if (!rows) {
    duplicates <- c(FALSE, unlist(lapply(1:(length(x) - 1), function(i) {
      duplicated(x[i:(i + 1)])[2]
    })))
  } else {
    duplicates <- c(FALSE, unlist(lapply(1:(nrow(x) - 1), function(i) {
      duplicated(x[i:(i + 1), ])[2]
    })))
  }
  return(duplicates)
}

#' Unwraps a matrix and only keep the first n elements.
#'
#' A function that unwraps a matrix and only keeps the first n elements
#' n can be either a constant (in which case it will be repeated), or a vector
#' @param mat A matrix
#' @param n A numeric vector
#' @return A logical vector of the same length as x
#' @export
#' @keywords internal
unwrap_rows <- function(mat, n) {
  n <- rep_len(n, nrow(mat))
  i <- 0
  out <- lapply(1:nrow(mat), function(i) {
    return(mat[i, 1:n[i]])
  })
  return(as.numeric(do.call(c, out)))
}


#' Splits a bbox into a list of bboxes less than a certain size
#'
#' @param x The maximum x size of the resulting bounding boxes
#' @param y The maximum y size of the resulting bounding boxes; defaults to x
#' @return A list of bbox objects
#' @export
#' @keywords internal
split_bbox <- function(bbox, x, y = x) {
  if (bbox[["xmin"]] > bbox[["xmax"]]) {
    x <- -1 * x
  }
  if (bbox[["ymin"]] > bbox[["ymax"]]) {
    y <- -1 * y
  }

  xs <- c(
    seq(
      bbox[["xmin"]],
      bbox[["xmax"]],
      x
    ),
    bbox["xmax"]
  )
  xs <-
    tibble::tibble(
      xmin = xs[1:(length(xs) - 1)],
      xmax = xs[2:length(xs)]
    )

  ys <- c(
    seq(
      bbox[["ymin"]],
      bbox[["ymax"]],
      y
    ),
    bbox[["ymax"]]
  )

  ys <-
    tibble::tibble(
      ymin = ys[1:(length(ys) - 1)],
      ymax = ys[2:length(ys)]
    )

  tidyr::crossing(xs, ys) %>%
    dplyr::rowwise() %>%
    dplyr::group_split() %>%
    purrr::map(as.list) %>%
    purrr::map(unlist) %>%
    purrr::map(
      magrittr::set_names,
      c("xmin", "xmax", "ymin", "ymax")
    ) %>%
    purrr::map(sf::st_bbox,
      crs = sf::st_crs(bbox)
    )
}


#' Use curl to download a file.
#'
#' This function makes it easy to implement timestamping and no-clobber of files.
#'
#' If both \code{timestamping} and \code{nc} are TRUE, nc behavior trumps timestamping.
#'
#' @param url The location of a file.
#' @param destdir Where the file should be downloaded to.
#' @param timestamping Should only newer files be downloaded?
#' @param nc Should files of the same type not be clobbered?
#' @param verbose Should cURL output be shown?
#' @param progress Should a progress bar be shown with cURL output?
#' @return A character string of the file path to the downloaded file.
#' @export
#' @keywords internal
download_data <-
  function(url,
           destdir = getwd(),
           timestamping = TRUE,
           nc = FALSE,
           verbose = FALSE,
           progress = FALSE) {
    destdir <- normalizePath(paste0(destdir, "/."))
    destfile <- paste0(destdir, "/", basename(url))
    temp.file <- paste0(tempdir(), "/", basename(url))

    if (nc & file.exists(destfile)) {
      message("Local file exists. Returning.")
      return(destfile)
    } else if (timestamping & file.exists(destfile)) {
      message("Downloading file (if necessary): ", url)
      opts <- list(
        verbose = verbose, noprogress = !progress,
        fresh_connect = TRUE, ftp_use_epsv = FALSE, forbid_reuse = TRUE,
        timecondition = TRUE, timevalue = base::file.info(destfile)$mtime
      )
      hand <- curl::new_handle()
      curl::handle_setopt(hand, .list = opts)
      tryCatch(status <- curl::curl_fetch_disk(url, path = temp.file, handle = hand),
        error = function(e) {
          message(
            "Download of ",
            url, " failed. Reverting to already cached file."
          )
          return(destfile)
        }
      )

      if (file.info(temp.file)$size > 0) {
        file.copy(temp.file, destfile, overwrite = TRUE)
      }
      return(destfile)
    } else {
      message("Downloading file: ", url)
      opts <- list(
        verbose = verbose,
        noprogress = !progress,
        fresh_connect = TRUE,
        ftp_use_epsv = FALSE,
        forbid_reuse = TRUE
      )
      hand <- curl::new_handle()
      curl::handle_setopt(hand, .list = opts)
      tryCatch(
        status <- curl::curl_fetch_disk(url,
          path = destfile,
          handle = hand
        ),
        error = function(e) stop("Download of ", url, " failed!")
      )
      return(destfile)
    }
    return(destfile)
  }

#' Check whether a web service is unavailable, and stop function if necessary.
#'
#' @param x The path to the web service.
#' @return Error if service unavailable.
#' @export
#' @keywords internal
check_service <- function(x) {
  if (x %>%
    httr::GET() %>%
    httr::status_code() %>%
    identical(200L) %>%
    magrittr::not()) {
    stop("Web service currently unavailable: ", source)
  }
}

#' Strip query parameters from a URL
#'
#' @param x The URL to be modified
#' @return The URL without parameters
#' @export
#' @keywords internal
url_base <- function(x) {
  x %<>% httr::parse_url()
  x$query <- list()
  x %<>% httr::build_url()
}

#' Replace NULLs
#'
#' @description Replace all the empty values in a list
#' @param x A list
#' @returns A list with NULLs replaced by NA
#' @examples
#' list(a = NULL, b = 1, c = list(foo = NULL, bar = NULL)) %>% replace_null()
#' @export

replace_null <- function(x) {
  is.na(x) <- x == "NULL"
  x
}

list_to_tibble <-
  function(x) {
    nms <- x %>%
      purrr::map(names) %>%
      purrr::reduce(union)

    test <- x %>%
      purrr::transpose(.names = nms) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ replace_null(.x))) %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ purrr::flatten(.x)))
  }

split_n <- function(x, n) {
  split(x, ceiling(seq_along(x) / n))
}

compare_rast_dims <-
  function(x, y) {
    x_dims <-
      c(
        terra::ncol(x),
        terra::nrow(x),
        terra::res(x)
      )

    y_dims <-
      c(
        terra::ncol(y),
        terra::nrow(y),
        terra::res(y)
      )

    all(x_dims == y_dims)
  }
