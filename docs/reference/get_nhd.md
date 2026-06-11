# Download and crop the National Hydrography Dataset.

`get_nhd` returns a list of
[`Simple Feature`](https://r-spatial.github.io/sf/reference/sf.html)
objects extracted from the National Hydrography Dataset.

## Usage

``` r
get_nhd(
  template,
  label,
  nhdplus = FALSE,
  extraction.dir = file.path(tempdir(), "FedData", "extractions", "nhd", label),
  force.redo = FALSE
)
```

## Arguments

- template:

  An
  [`Simple Feature`](https://r-spatial.github.io/sf/reference/sf.html)
  or
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object to serve as a template for cropping.

- label:

  A character string naming the study area.

- nhdplus:

  Extract data from the USGS NHDPlus High Resolution service
  (experimental)

- extraction.dir:

  A character string indicating where the extracted and cropped NHD data
  should be put.

- force.redo:

  If an extraction for this template and label already exists, should a
  new one be created?

## Value

A list of `sf` collections extracted from the National Hydrography
Dataset.

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
