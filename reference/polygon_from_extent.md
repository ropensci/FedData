# Turn an extent object into a polygon

Turn an extent object into a polygon

## Usage

``` r
polygon_from_extent(x, proj4string = NULL)
```

## Arguments

- x:

  An object from which an bounding box object can be retrieved.

- proj4string:

  A PROJ.4 formatted string defining the required projection.

## Value

A [`Simple Feature`](https://r-spatial.github.io/sf/reference/sf.html)
object.
