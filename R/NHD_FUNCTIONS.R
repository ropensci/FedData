#' Download and crop the National Hydrography Dataset.
#'
#' \code{get_nhd} returns a list of Spatial* objects extracted 
#' from the National Hydrography Dataset.
#' 
#' @param template A Raster* or Spatial* object to serve 
#' as a template for cropping.
#' @param label A character string naming the study area.
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing. Defaults to './RAW/NHD/'.
#' @param extraction.dir A character string indicating where the extracted and cropped NHD shapefiles should be put.
#' The directory will be created if missing. Defaults to './EXTRACTIONS/NHD/'.
#' @param force.redo If an extraction for this template and label already exists, should a new one be created?
#' @return A list of Spatial* objects extracted from the National Hydrography Dataset.
#' @importFrom magrittr %>% %<>%
#' @export
#' @examples
#' \dontrun{
#' # Extract data for the Village Ecodynamics Project 'VEPIIN' study area:
#' # http://village.anth.wsu.edu
#' vepPolygon <- polygon_from_extent(raster::extent(672800,740000,4102000,4170000), 
#'      proj4string='+proj=utm +datum=NAD83 +zone=12')
#' 
#' # Get the NHD (USA ONLY)
#' NHD <- get_nhd(template=vepPolygon, label='VEPIIN')
#' 
#' # Plot the VEP polygon
#' plot(vepPolygon)
#' 
#' # Plot the NHD data
#' plot(NHD$Flowline, add=T)
#' plot(NHD$Line, add=T)
#' plot(NHD$Area, col='black', add=T)
#' plot(NHD$Waterbody, col='black', add=T)
#' }
get_nhd <- function(template,
                    label,
                    raw.dir = paste0(tempdir(),"/FedData/raw/nhd"),
                    extraction.dir = paste0(tempdir(),"/FedData/extractions/nhd/",label,"/"),
                    force.redo = FALSE) {
  
  raw.dir <- normalizePath(paste0(raw.dir,"/"), mustWork = FALSE)  
  extraction.dir <- normalizePath(paste0(extraction.dir,"/"), mustWork = FALSE)
  
  dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)
  
  out_dsn <- paste0(extraction.dir, label, "_nhd.gpkg")
  
  if (!force.redo & file.exists(out_dsn)) {
    return(read_sf_all(out_dsn))
  }
  
  unlink(out_dsn, 
         recursive = TRUE, 
         force = TRUE)
  
  template %<>%
    template_to_sf()
  
  message("(Down)Loading the NHD HUC4 dataset.")
  HUC4 <- 
    get_huc4(template = template, 
             raw.dir = raw.dir)
  
  area_list <- 
    HUC4$HUC4 %>% 
    as.character()
  
  # Get the spatial data for each area
  message("(Down)Loading the NHD subregion data.")
  subregion_shapes <- 
    area_list %>%
    purrr::map(
      ~get_nhd_subregion(template = template,
                         area = .x,
                         raw.dir = raw.dir)
    )
  
  shapes <-
    subregion_shapes %>%
    purrr::transpose() %>%
    purrr::map(~do.call(rbind, .x))
  
  shapes %>%
    purrr::iwalk(
      ~sf::write_sf(.x,
                    dsn = paste0(extraction.dir, label, "_nhd.gpkg"),
                    layer = .y,
                    delete_layer = TRUE)
    )
  
  return(read_sf_all(out_dsn))
}

#' Download a zipped directory containing a shapefile of the HUC4 subregions of the NHD.
#'
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' @return A character string representing the full local path of the HUC4 zipped directory.
#' @export
#' @keywords internal
download_huc4 <- function(raw.dir) {
  # url <- 'ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/WBD/FileGDB101/WBD_National.zip'
  # url <- "ftp://ftp.igsb.uiowa.edu/gis_library/USA/huc_04.zip"
  url <- "https://github.com/ropensci/FedData/raw/master/data-raw/nhd_huc4.gpkg.zip"
  tryCatch(download_data(url = url, destdir = raw.dir), error = function(e){message("HUC4 download not available. Using local version.")})
  return(normalizePath(paste(raw.dir, "/nhd_huc4.gpkg.zip", sep = "")))
  # return(normalizePath(paste(raw.dir,'WBD_National.zip',sep='')))
}


#' Download and crop a shapefile of the HUC4 
#' regions of the National Hydrography Dataset.
#'
#' \code{get_huc4} returns a \code{SpatialPolygonsDataFrame} of the HUC4 regions within
#' the specified \code{template}. If template is not provided, returns the entire HUC4 dataset.
#' 
#' @param template A Raster* or Spatial* object to serve 
#' as a template for cropping.
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing.
#' @return A \code{SpatialPolygonsDataFrame} of the HUC4 regions within
#' the specified \code{template}.
#' @export
#' @keywords internal
get_huc4 <- function(template = NULL, raw.dir) {
  tmpdir <- tempfile()
  if (!dir.create(tmpdir)) 
    stop("failed to create my temporary directory")
  
  huc4File <- download_huc4(raw.dir)
  
  utils::unzip(huc4File, exdir = tmpdir)
  
  HUC4 <- sf::st_read(stringr::str_c(tmpdir,"/data-raw/nhd_huc4.gpkg"),
                      quiet = TRUE)
  
  # Get a list of NHD subregions within the project study area
  if (!is.null(template)) {
    
    template %<>%
      template_to_sf()
    
    suppressMessages(
      suppressWarnings(
        HUC4 %<>%
          sf::st_intersection(template %>% 
                                sf::st_transform(
                                  sf::st_crs(HUC4)
                                )
          )
      )
    )
  }
  
  unlink(tmpdir, recursive = TRUE)
  
  return(HUC4)
}



#' Download a zipped NHD HUC4 subregion.
#'
#' HUC4 subregion are specified by a unique character string, best obtained using the \code{\link{get_huc4}} function.
#' \code{download_nhd_subregion} returns the path to the downloaded zip file.
#'
#' @param area A 4-character string indicating the HUC4 NHD subregion to download.
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing.
#' @return A character string representing the full local path of the downloaded zip file.
#' @export
#' @keywords internal
download_nhd_subregion <- function(area, raw.dir) {
  url <- paste0("https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU4/HighResolution/GDB/NHD_H_",area,"_HU4_GDB.zip")
  
  destdir <- raw.dir
  download_data(url = url, destdir = destdir)
  
  return(normalizePath(paste0(destdir, "/", basename(url))))
}

#' Download and crop data from a zipped HUC4 subregion
#' of the National Hydrography Dataset.
#'
#' \code{get_nhd_subregion} returns a list of \code{SpatialPolygonsDataFrame}s of the layers of the HUC4 subregion,
#'  within the specified \code{template}. If template is not provided, returns the entire HUC4 subregion.
#' 
#' @param template A Raster* or Spatial* object to serve 
#' as a template for cropping.
#' @param area A 4-character string indicating the HUC4 NHD subregion to download and crop.
#' @param raw.dir A character string indicating where raw downloaded files should be put.
#' The directory will be created if missing.
#' @return A \code{SpatialPolygonsDataFrame} of the HUC4 regions within
#' the specified \code{template}.
#' @export
#' @keywords internal
get_nhd_subregion <- function(template = NULL, area, raw.dir) {
  tmpdir <- tempfile()
  if (!dir.create(tmpdir)) 
    stop("failed to create my temporary directory")
  
  file <- 
    download_nhd_subregion(area = area, 
                           raw.dir = raw.dir)
  
  utils::unzip(file, exdir = tmpdir)
  
  # Get the path to the geodatabase
  dsn <- 
    list.files(tmpdir, 
               full.names = T, 
               pattern = "*.gdb") %>%
    normalizePath()
  
  layers <- 
    dsn %>%
    sf::st_layers()
  
  suppressWarnings(
    suppressMessages(
      shapes <- 
        layers[['name']][!is.na(unlist(layers['geomtype']))] %>%
        stringr::str_subset("NHD") %>%
        stringr::str_subset("EventFC", negate = TRUE) %>%
        magrittr::set_names(.,.) %>%
        purrr::map(function(x){
          shape <-
            sf::read_sf(dsn = dsn, layer = x)
          
          shape %>%
            sf::st_zm() %>%
            sf::st_intersection(template %>%
                                  sf::st_transform(
                                    sf::st_crs(shape)
                                  )
            )
          
        }) %>%
        purrr::compact()
    )
  )
  
  names(shapes) %<>%
    stringr::str_remove("NHD")
  
  unlink(tmpdir, recursive = TRUE)
  
  return(shapes)
}