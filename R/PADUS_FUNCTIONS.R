#' Download and crop the PAD-US Dataset.
#'
#' `get_padus` returns a list of `sf` objects extracted
#' from the PAD-US Dataset. Data are retrieved directly from
#' [PAD-US ArcGIS Web Services](https://www.usgs.gov/programs/gap-analysis-project/science/pad-us-web-services).
#'
#' [PAD-US](https://www.usgs.gov/programs/gap-analysis-project/science/pad-us-data-overview) is America’s official national inventory of U.S. terrestrial and
#' marine protected areas that are dedicated to the preservation of biological
#' diversity and to other natural, recreation and cultural uses, managed for
#' these purposes through legal or other effective means. PAD-US also includes
#' the best available aggregation of federal land and marine areas provided
#' directly by managing agencies, coordinated through the Federal Geographic
#' Data Committee Federal Lands Working Group.
#'
#' @param template An [`Simple Feature`][sf::sf]
#' or [`SpatRaster`][terra::SpatRaster] object to serve as a template for cropping.
#' Optionally, a vector of unit names, e.g., `c('Mesa Verde National Park','Ute Mountain Reservation')` may be provided.
#' @param label A character string naming the study area.
#' @param layer A character vector containing one or more PAD-US Layers.
#' By default, the **Manager_Name** layer is downloaded.
#' \itemize{
#' \item **Protection_Status_by_GAP_Status_Code**: [PAD-US Protection Status by GAP Status Code](https://usgs.maps.arcgis.com/home/item.html?id=98fce3fb0c8241ce8847e9f7d0d212e9) — Service representing a measure of management intent to permanently protect biodiversity. GAP 1&2 areas are primarily managed for biodiversity, GAP 3 are managed for multiple uses including conservation and extraction, GAP 4 no known mandate for biodiversity protection. GAP Status Codes 1-3 are displayed, GAP 4 areas included but not displayed.
#' \item **Public_Access**: [PAD-US Public Access](https://usgs.maps.arcgis.com/home/item.html?id=c91a5655a1be428daeb778888e60db24) — Service representing general level of public access permitted in the area - Open, Restricted (permit, seasonal), Closed. Public Access Unknown areas not displayed. Use to show general categories of public access (however, not all areas have been locally reviewed).
#' \item **Fee_Managers**: [PAD-US Fee Managers](https://usgs.maps.arcgis.com/home/item.html?id=a516f406610b42d8a755c25665324201) — Service providing manager or administrative agency names standardized nationally. Use for categorization by manager name, with detailed federal managers and generic state/local/other managers. Where available this layer includes fee simple parcels from the Fee feature class plus DOD and Tribal areas from the Proclamation feature class.
#' \item **Manager_Name**: [PAD-US Manager Name](https://usgs.maps.arcgis.com/home/item.html?id=d1c464bac63d44bc88a4a059dc69d259) — Service representing fine level manager or administrative agency name standardized for the Nation (USFS, BLM, State Fish and Wildlife, State Parks and Rec, City, NGO, etc). This map is based on the PAD-US Combined Proclamation, Marine, Fee, Designation, Easement feature class. DOD and Tribal areas shown with 50% transparency. Use for categorization by manager name, with detailed federal managers and generic state/local/other managers.
#' \item **Manager_Type**: [PAD-US Manager Type](https://usgs.maps.arcgis.com/home/item.html?id=f71b5910fbcf4349b384b298f1247eb7) — Service representing coarse level land manager description from "Agency Type" Domain, "Manager Type" Field (for example, Federal, Tribal, State, Local Gov, Private). Use for broad categorization of manager levels, for general depictions of who manages what areas.
#' \item **Federal_Fee_Managers_Authoritative**: [PAD-US Federal Fee Managers Authoritative](https://usgs.maps.arcgis.com/home/item.html?id=995e04a677374976931c502c2b095402) — Service describing authoritative fee data for federal managers or administrative agencies by name. U.S. Department of Defense and Tribal areas shown from the Proclamation feature class. Use to depict authoritative fee data for individual federal management agencies (no state, local or private lands). This service does not include designations that often overlap state, private or other inholdings. U.S. Department of Defense internal land ownership is not represented but is implied Federal. See the Federal Management Agencies service for a combined view of fee ownership, designations, and easements.
#' \item **Federal_Management_Agencies**: [PAD-US Federal Management Agencies](https://usgs.maps.arcgis.com/home/item.html?id=d85301944bc048fa8d07b6040b4c9f69) — Service providing Federal managers or administrative agencies by name. Use to depict individual federal management agencies (no state, local or private lands). This map is based on the Combined Proclamation, Marine, Fee, Designation, Easement feature class.
#' \item **Protection_Mechanism_Category**: [PAD-US Protection Mechanism Category](https://usgs.maps.arcgis.com/home/item.html?id=80c7a138989c4cf8937c4502371d6e99) — Service representing the protection mechanism category including fee simple, internal management designations, easements, leases and agreements, and Marine Areas. Use to show categories of land tenure for all protected areas, including marine areas.
#' \item **Proclamation_and_Other_Planning_Boundaries**: [PAD-US Proclamation and Other Planning Boundaries](https://usgs.maps.arcgis.com/home/item.html?id=d478e6d7848a498ba3691b54a2f7cb95) — Service representing boundaries that provide additional context. Administrative agency name standardized for the nation (DOD, FWS, NPS, USFS, Tribal). Boundaries shown with outline only, as proclamation data do not depict actual ownership or management. Use to show outline of agency proclamation, approved acquisition or other planning boundaries where internal ownership is not depicted.
#' }
#' @param extraction.dir A character string indicating where the extracted and cropped PAD-US data should be put.
#' @param force.redo If an extraction for this template and label already exists, should a new one be created?
#' @return A list of [sf::sf] collections extracted from the PAD-US Dataset.
#' @importFrom magrittr %>% %<>%
#' @export
#' @examples
#' \dontrun{
#' # Get the PAD-US (USA ONLY)
#' PADUS <- get_padus(
#'   template = FedData::meve,
#'   label = "meve"
#' )
#' PADUS
#' }
get_padus <-
  function(template,
           label,
           layer = c(
             # "Protection_Status_by_GAP_Status_Code",
             # "Public_Access",
             # "Fee_Manager",
             "Manager_Name" # ,
             # "Manager_Type",
             # "Federal_Fee_Managers_Authoritative",
             # "Federal_Management_Agencies",
             # "Protection_Mechanism_Category",
             # "Proclamation_and_Other_Planning_Boundaries",
             # "Fee_Topology_Overlaps"
           ),
           extraction.dir = file.path(
             tempdir(),
             "FedData",
             "extractions",
             "padus",
             label
           ),
           force.redo = FALSE) {
    padus_base_url <-
      "https://services.arcgis.com/v01gqwM5QqNysAAi/arcgis/rest/services"

    padus_services <-
      c(
        "Protection_Status_by_GAP_Status_Code" = "PADUS_Protection_Status_by_GAP_Status_Code",
        "Public_Access" = "PADUS_Public_Access",
        "Fee_Managers" = "Fee_Managers_PADUS",
        "Manager_Name" = "Manager_Name_PADUS",
        "Manager_Type" = "Manager_Type_PADUS",
        "Federal_Fee_Managers_Authoritative" = "Federal_Fee_Managers_Authoritative_PADUS",
        "Federal_Management_Agencies" = "Federal_Management_Agencies_PADUS",
        "Protection_Mechanism_Category" = "Protection_Mechanism_Category_PADUS",
        "Proclamation_and_Other_Planning_Boundaries" = "Proclamation_and_Other_Planning_Boundaries_PADUS"
      )

    if (any(!(layer %in% names(padus_services)))) {
      # layer[which(!(layer %in% names(padus_services)))]
      stop("Requested layers must be one or more of the available layers. Please see `?get_padus`.")
    }

    extraction.dir <-
      normalizePath(extraction.dir,
        mustWork = FALSE
      )

    dir.create(extraction.dir, showWarnings = FALSE, recursive = TRUE)
    out_file <- paste0(label, "_padus.gpkg")

    out_dsn <- file.path(extraction.dir, out_file)

    if (!force.redo & file.exists(out_dsn)) {
      return(read_sf_all(out_dsn))
    }

    if (inherits(template, "character")) {
      padus_out <-
        padus_services[layer] %>%
        purrr::map(
          function(x) {
            file.path(padus_base_url, x, "FeatureServer") %>%
              httr::parse_url() %>%
              httr::modify_url(.,
                path = c(.$path, 0, "query"),
                query = list(
                  outFields = "*",
                  f = "pgeojson",
                  where =
                    paste0(
                      "Unit_Nm IN (",
                      paste(paste0("'", template, "'"), collapse = ","),
                      ")"
                    )
                )
              ) %>%
              sf::read_sf()
          }
        )
    } else {
      padus_out <-
        padus_services[layer] %>%
        purrr::map(
          function(x) {
            agol_filter_httr(
              url = file.path(padus_base_url, x, "FeatureServer"),
              layer_name = 0,
              geom = template
            )
          }
        )
    }

    null_elements <- purrr::map_lgl(padus_out, is.null)
    if (all(null_elements)) stop("No PAD-US data present within template.")
    padus_for_crs <- padus_out[!null_elements][[1]]

    padus_out %<>%
      purrr::compact() %>%
      purrr::map(sf::st_make_valid)

    if (!inherits(template, "character")) {
      sf_state <- sf::sf_use_s2()
      suppressMessages(sf::sf_use_s2(FALSE))
      suppressWarnings({
        suppressMessages({
          padus_out %<>%
            purrr::map(
              sf::st_intersection,
              template %>%
                sf::st_geometry() %>%
                sf::st_transform(sf::st_crs(padus_for_crs))
            )
        })
      })
      suppressMessages(sf::sf_use_s2(sf_state))
    }

    write_sf_all(padus_out, dsn = out_dsn)

    return(read_sf_all(out_dsn))
  }
