# Scaffolds the common pattern of selecting a layer and filter a geometry from an ArcGIS feature service.

This function does **not** use the `arcgislayers` package, which has has
had compatibility issues for several commonly used platforms.

## Usage

``` r
agol_filter_httr(url, layer_name = NULL, geom, simplify = TRUE)
```

## Arguments

- url:

  the url of the remote resource. Must be of length one.

- layer_name:

  the name(s) associated with the layer you want to retrieve. Can be a
  character vector. If `NULL` (the default), iterates through all
  layers.

- geom:

  an object of class `bbox`, `sfc` or `sfg` used to filter query results
  based on a predicate function.

- simplify:

  when only one layer exists, just return the `sf` object or
  `data.frame`, otherwise return a list of these objects.

## Value

An `sf` object, or a `data.frame`, or a list of these objects if
`layer_name == NULL` or if `length(layer_name) > 1`. Missing layers
return "NULL".

## Examples

``` r
if (FALSE) { # \dontrun{

# Get a single layer
agol_filter_httr(
  url = "https://hydro.nationalmap.gov/arcgis/rest/services/NHDPlus_HR/MapServer/",
  layer_name = "WBDHU12",
  geom = FedData::meve
)

# Can be returned as a list
agol_filter_httr(
  url = "https://hydro.nationalmap.gov/arcgis/rest/services/NHDPlus_HR/MapServer/",
  layer_name = "WBDHU12",
  geom = FedData::meve,
  simplify = FALSE
)

# Get a list with all layers
agol_filter_httr(
  url = "https://hydro.nationalmap.gov/arcgis/rest/services/NHDPlus_HR/MapServer/",
  geom = FedData::meve
)

# Or include a vector of layer names
# Note that missing layers are returned as `NULL` values
agol_filter_httr(
  url = "https://hydro.nationalmap.gov/arcgis/rest/services/NHDPlus_HR/MapServer/",
  layer_name = c(
    "NHDPoint",
    "NetworkNHDFlowline",
    "NonNetworkNHDFlowline",
    "NHDLine",
    "NHDArea",
    "NHDWaterbody"
  ),
  geom = FedData::meve
)
} # }
```
