# A basic plotting function for NHD data.

This is more of an example than anything

## Usage

``` r
plot_nhd(x, template = NULL)
```

## Arguments

- x:

  The result of
  [get_nhd](https://docs.ropensci.org/FedData/reference/get_nhd.md).

- template:

  An
  [`Simple Feature`](https://r-spatial.github.io/sf/reference/sf.html)
  or
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object to serve as a template for cropping.

## Value

A `ggplot2` panel of plots

## Examples

``` r
if (FALSE) { # \dontrun{
# Get the NHD (USA ONLY)
NHD <- get_nhd(
  template = FedData::meve,
  label = "meve"
)
NHD
NHD %>%
  plot_nhd(template = FedData::meve)
} # }
```
