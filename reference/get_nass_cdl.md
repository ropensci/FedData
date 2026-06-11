# Download and crop the NASS Cropland Data Layer.

`get_nass_cdl` returns a
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
of NASS Cropland Data Layer cropped to a given template study area.

## Usage

``` r
get_nass_cdl(
  template,
  label,
  year = 2019,
  extraction.dir = paste0(tempdir(), "/FedData/"),
  raster.options = c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"),
  force.redo = FALSE,
  progress = TRUE
)

get_nass(template, label, ...)

get_cdl(template, label, ...)

cdl_colors()
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

- year:

  An integer representing the year of desired NASS Cropland Data Layer
  product. Acceptable values are 2007–the last year.

- extraction.dir:

  A character string indicating where the extracted and cropped NASS
  data should be put. The directory will be created if missing.

- raster.options:

  a vector of options for terra::writeRaster.

- force.redo:

  If an extraction for this template and label already exists, should a
  new one be created?

- progress:

  Draw a progress bar when downloading?

- ...:

  Other parameters passed on to get_nass_cdl.

## Value

A
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
cropped to the bounding box of the template.

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract data for the Mesa Verde National Park:

# Get the NASS CDL (USA ONLY)
# Returns a raster
NASS <-
  get_nass_cdl(
    template = FedData::meve,
    label = "meve",
    year = 2011
  )

# Plot with terra::plot
terra::plot(NASS)
} # }
```
