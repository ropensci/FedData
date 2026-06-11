# Download and crop a shapefile of the SSURGO study areas.

`get_ssurgo_inventory` returns a `SpatialPolygonsDataFrame` of the
SSURGO study areas within the specified `template`. If template is not
provided, returns the entire SSURGO inventory of study areas.

## Usage

``` r
get_ssurgo_inventory(template = NULL, raw.dir)
```

## Arguments

- template:

  An
  [`Simple Feature`](https://r-spatial.github.io/sf/reference/sf.html)
  or
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object to serve as a template for cropping.

- raw.dir:

  A character string indicating where raw downloaded files should be
  put. The directory will be created if missing.

## Value

A `SpatialPolygonsDataFrame` of the SSURGO study areas within the
specified `template`.
