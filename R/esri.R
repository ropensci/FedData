esri_describe <-
  function(url) {
    description <-
      url %>%
      httr::GET(query = list(f = "json")) %>%
      httr::stop_for_status(
        task =
          paste0("describe the NHD dataset at:\n", url)
      ) %>%
      httr::content(
        as = "text",
        encoding = "UTF-8"
      ) %>%
      jsonlite::fromJSON()

    description
  }

esri_query <-
  function(url,
           layers = NULL,
           geom = NULL) {
    url %<>%
      url_base()

    all_layers <- esri_describe(url)$layers

    layer_ids <-
      all_layers$id %>%
      magrittr::set_names(stringr::str_trim(all_layers$name))

    if (!is.null(layers)) {
      layer_ids %<>%
        magrittr::extract(layers)
    }

    query <-
      list(
        where = "1=1",
        outFields = "*",
        f = "geoJSON"
      )

    if (!is.null(geom)) {
      query$geometry <-
        geom %>%
        sf::st_transform(4326) %>%
        sf::st_bbox() %>%
        paste0(collapse = ",")

      query$geometryType <- "esriGeometryEnvelope"

      query$inSR <- 4326
    }

    layer_ids %>%
      purrr::map(function(x) {
        max_count <-
          esri_describe(url)$maxRecordCount / 2

        ids <-
          httr::POST(
            url = paste0(url, x, "/query"),
            body = c(query,
              returnIdsOnly = TRUE
            )
          ) %>%
          httr::content(
            as = "text",
            encoding = "UTF-8"
          ) %>%
          jsonlite::fromJSON() %>%
          magrittr::extract2("objectIds")

        if (is.null(ids)) {
          return(NULL)
        }

        ids %<>%
          split_n(max_count)

        ids %>%
          purrr::map(function(i) {
            httr::POST(
              url = paste0(url, x, "/query"),
              body = list(
                where = "1=1",
                outFields = "*",
                f = "geoJSON",
                objectIds = paste0(i, collapse = ",")
              )
            ) %>%
              httr::content(
                as = "text",
                encoding = "UTF-8"
              ) %>%
              sf::read_sf()
          }) %>%
          do.call("rbind", .)
      })
  }
