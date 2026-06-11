# Download and crop the Watershed Boundary Dataset.

`get_wbd` returns an
[`Simple Feature`](https://r-spatial.github.io/sf/reference/sf.html)
collection of the HUC 12 regions within the specified `template`.

## Usage

``` r
get_wbd(
  template,
  label,
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

- extraction.dir:

  A character string indicating where the extracted and cropped NHD data
  should be put.

- force.redo:

  If an extraction for this template and label already exists, should a
  new one be created?

## Value

An `sf` collection of the HUC 12 regions within the specified
`template`.
