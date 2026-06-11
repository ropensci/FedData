# Load and crop tile from the 1 (~30 meter) or 1/3 (~10 meter) arc-second National Elevation Dataset.

`get_ned_tile` returns a`SpatRaster` cropped within the specified
`template`. If template is not provided, returns the entire NED tile.

## Usage

``` r
get_ned_tile(template = NULL, res = "1", tileNorthing, tileWesting)
```

## Arguments

- template:

  An
  [`Simple Feature`](https://r-spatial.github.io/sf/reference/sf.html)
  or
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object to serve as a template for cropping. If missing, entire tile is
  returned.

- res:

  A character string representing the desired resolution of the NED. '1'
  indicates the 1 arc-second NED (the default), while '13' indicates the
  1/3 arc-second dataset.

- tileNorthing:

  An integer representing the northing (latitude, in degrees north of
  the equator) of the northwest corner of the tile to be downloaded.

- tileWesting:

  An integer representing the westing (longitude, in degrees west of the
  prime meridian) of the northwest corner of the tile to be downloaded.

## Value

A `SpatRaster` cropped to the extent of the template.
