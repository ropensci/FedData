# Download and crop the 1-km DAYMET v4 daily weather dataset.

`get_daymet()` is deprecated, and currently non-functional. In 2025, the
ORNL Distributed Active Archive Center retired the public THREDDS data
server that FedData used to subset and download gridded Daymet data.
Gridded Daymet data are now distributed through the NASA Earthdata
Cloud, which requires (free) authentication with a [NASA Earthdata
Login](https://urs.earthdata.nasa.gov), and are also available from
[Google Earth
Engine](https://developers.google.com/earth-engine/datasets/catalog/NASA_ORNL_DAYMET_V4).
Point ("single-pixel") extractions remain freely available via the
[Daymet single-pixel extraction tool](https://daymet.ornl.gov/getdata),
accessible from R with the
[daymetr](https://github.com/bluegreen-labs/daymetr) package. Support
for gridded Daymet downloads with Earthdata authentication may return in
a future release of FedData.

`get_daymet` returned a
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
of weather data cropped to a given template study area.

## Usage

``` r
get_daymet(
  template,
  label,
  elements = c("dayl", "prcp", "srad", "swe", "tmax", "tmin", "vp"),
  years = 1980:(lubridate::year(Sys.time()) - 1),
  region = "na",
  tempo = "day",
  extraction.dir = file.path(tempdir(), "FedData", "extractions", "daymet", label),
  raster.options = c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"),
  force.redo = FALSE,
  progress = TRUE
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

- elements:

  A character vector of elements to extract.\
  The available elements are:\
  dayl = Duration of the daylight period in seconds per day. This
  calculation is based on the period of the day during which the sun is
  above a hypothetical flat horizon.\
  prcp = Daily total precipitation in millimeters per day, sum of all
  forms converted to water-equivalent. Precipitation occurrence on any
  given day may be ascertained.\
  srad = Incident shortwave radiation flux density in watts per square
  meter, taken as an average over the daylight period of the day. NOTE:
  Daily total radiation (MJ/m2/day) can be calculated as follows: ((srad
  (W/m2) \* dayl (s/day)) / l,000,000)\
  swe = Snow water equivalent in kilograms per square meter. The amount
  of water contained within the snowpack.\
  tmax = Daily maximum 2-meter air temperature in degrees Celsius.\
  tmin = Daily minimum 2-meter air temperature in degrees Celsius.\
  vp = Water vapor pressure in pascals. Daily average partial pressure
  of water vapor.\

- years:

  A numeric vector of years to extract.

- region:

  The name of a region. The available regions are:\
  na = North America\
  hi = Hawaii\
  pr = Puerto Rico\

- tempo:

  The frequency of the data. The available tempos are:\
  day = Daily data\
  mon = Monthly summary data\
  ann = Annual summary data\

- extraction.dir:

  A character string indicating where the extracted and cropped DEM
  should be put. Defaults to a temporary directory.

- raster.options:

  a vector of GDAL options passed to
  [terra::writeRaster](https://rspatial.github.io/terra/reference/writeRaster.html).

- force.redo:

  If an extraction for this template and label already exists in
  extraction.dir, should a new one be created?

- progress:

  Draw a progress bar when downloading?

## Value

A named list of `SpatRaster`s of weather data cropped to the extent of
the template.

## Examples

``` r
if (FALSE) { # \dontrun{
library(terra)

# Get the DAYMET (North America only)
# Returns a list of raster bricks
DAYMET <- get_daymet(
  template = FedData::meve,
  label = "meve",
  elements = c("prcp", "tmin", "tmax"),
  years = 1985
)

# Plot with terra::plot
plot(DAYMET$tmin$`1985-10-23`)
} # }
```
