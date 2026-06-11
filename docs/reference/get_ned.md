# Download and crop the 1 (~30 meter) or 1/3 (~10 meter) arc-second National Elevation Dataset.

`get_ned` returns a `SpatRaster` of elevation data cropped to a given
template study area.

## Usage

``` r
get_ned(
  template,
  label,
  res = "1",
  extraction.dir = file.path(tempdir(), "FedData", "extractions", "ned", label),
  raster.options = c("COMPRESS=DEFLATE", "ZLEVEL=9"),
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

- res:

  A character string representing the desired resolution of the NED. '1'
  indicates the 1 arc-second NED (the default), while '13' indicates the
  1/3 arc-second dataset.

- extraction.dir:

  A character string indicating where the extracted and cropped DEM
  should be put. The directory will be created if missing.

- raster.options:

  a vector of GDAL options passed to
  [terra::writeRaster](https://rspatial.github.io/terra/reference/writeRaster.html).

- force.redo:

  If an extraction for this template and label already exists, should a
  new one be created?

## Value

A `SpatRaster` DEM cropped to the extent of the template.

## Examples

``` r
if (FALSE) { # \dontrun{
# Get the NED (USA ONLY)
# Returns a `SpatRaster`
NED <-
  get_ned(
    template = FedData::meve,
    label = "meve"
  )

# Plot with terra::plot
terra::plot(NED)
} # }
```
