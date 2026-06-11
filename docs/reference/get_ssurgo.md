# Download and crop data from the NRCS SSURGO soils database.

This is an efficient method for spatially merging several different soil
survey areas as well as merging their tabular data.

## Usage

``` r
get_ssurgo(
  template,
  label,
  raw.dir = paste0(tempdir(), "/FedData/raw/ssurgo"),
  extraction.dir = paste0(tempdir(), "/FedData/"),
  force.redo = FALSE
)
```

## Arguments

- template:

  An
  [`Simple Feature`](https://r-spatial.github.io/sf/reference/sf.html)
  or
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object to serve as a template for cropping. Optionally, a vector of
  area names, e.g., `c('IN087','IN088')` may be provided.

- label:

  A character string naming the study area.

- raw.dir:

  A character string indicating where raw downloaded files should be
  put. The directory will be created if missing. Defaults to
  './RAW/SSURGO/'.

- extraction.dir:

  A character string indicating where the extracted and cropped SSURGO
  shapefiles should be put. The directory will be created if missing.
  Defaults to './EXTRACTIONS/SSURGO/'.

- force.redo:

  If an extraction for this template and label already exists, should a
  new one be created? Defaults to FALSE.

## Value

A named list containing the 'spatial' and 'tabular' data.

## Details

`get_ssurgo` returns a named list of length 2:

1.  'spatial': A
    [`Simple Feature`](https://r-spatial.github.io/sf/reference/sf.html)
    of soil mapunits in the template, and

2.  'tabular': A named list of
    [`data.frame`](https://rdrr.io/r/base/data.frame.html)`s` with the
    SSURGO tabular data.

## Examples

``` r
if (FALSE) { # \dontrun{
# Get the NRCS SSURGO data (USA ONLY)
SSURGO.MEVE <-
  get_ssurgo(
    template = FedData::meve,
    label = "meve"
  )

# Plot the VEP polygon
plot(meve)

# Plot the SSURGO mapunit polygons
plot(SSURGO.MEVE$spatial["MUKEY"],
  lwd = 0.1,
  add = TRUE
)

# Or, download by Soil Survey Area names
SSURGO.areas <-
  get_ssurgo(
    template = c("CO670", "CO075"),
    label = "CO_TEST"
  )

# Let's just look at spatial data for CO675
SSURGO.areas.CO675 <-
  SSURGO.areas$spatial[SSURGO.areas$spatial$AREASYMBOL == "CO075", ]

# And get the NED data under them for pretty plotting
NED.CO675 <-
  get_ned(
    template = SSURGO.areas.CO675,
    label = "SSURGO_CO675"
  )

# Plot the SSURGO mapunit polygons, but only for CO675
terra::plot(NED.CO675)
plot(
  SSURGO.areas.CO675$geom,
  lwd = 0.1,
  add = TRUE
)
} # }
```
