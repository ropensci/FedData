# Download and crop the inventory of GHCN stations.

`get_ghcn_inventory` returns a `SpatialPolygonsDataFrame` of the GHCN
stations within the specified `template`. If template is not provided,
returns the entire GHCN inventory.

## Usage

``` r
get_ghcn_inventory(template = NULL, elements = NULL, raw.dir)
```

## Arguments

- template:

  An
  [`Simple Feature`](https://r-spatial.github.io/sf/reference/sf.html)
  or
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object to serve as a template for cropping.

- elements:

  A character vector of elements to extract. Common elements include
  'tmin', 'tmax', and 'prcp'.

- raw.dir:

  A character string indicating where raw downloaded files should be
  put. The directory will be created if missing.

## Value

A [`Simple Feature`](https://r-spatial.github.io/sf/reference/sf.html)
of the GHCN stations within the specified `template`

## Details

Stations with multiple elements will have multiple points. This allows
for easy mapping of stations by element availability.
