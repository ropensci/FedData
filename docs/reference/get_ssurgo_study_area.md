# Download and crop the spatial and tabular data for a SSURGO study area.

`get_ssurgo_study_area` returns a named list of length 2:

1.  'spatial': A
    [`Simple Feature`](https://r-spatial.github.io/sf/reference/sf.html)
    of soil mapunits in the template, and

2.  'tabular': A named list of
    [`data.frame`](https://rdrr.io/r/base/data.frame.html)`s` with the
    SSURGO tabular data.

## Usage

``` r
get_ssurgo_study_area(template = NULL, area, date, raw.dir)
```

## Arguments

- template:

  An
  [`Simple Feature`](https://r-spatial.github.io/sf/reference/sf.html)
  or
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object to serve as a template for cropping. If missing, whose study
  area is returned

- area:

  A character string indicating the SSURGO study area to be downloaded.

- date:

  A character string indicating the date of the most recent update to
  the SSURGO area for these data. This information may be gleaned from
  the SSURGO Inventory
  ([`get_ssurgo_inventory`](https://docs.ropensci.org/FedData/reference/get_ssurgo_inventory.md)).

- raw.dir:

  A character string indicating where raw downloaded files should be
  put. The directory will be created if missing.

## Value

A `SpatialPolygonsDataFrame` of the SSURGO study areas within the
specified `template`.
