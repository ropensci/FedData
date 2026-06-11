# Download and crop the National Land Cover Database.

`get_nlcd` returns a
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
of NLCD data cropped to a given template study area. `nlcd_colors` and
`pal_nlcd` return the NLCD legend and color palette, as available
through the [MLRC
website](https://www.mrlc.gov/data/legends/national-land-cover-database-class-legend-and-description).

## Usage

``` r
get_nlcd(
  template,
  label,
  year = 2021,
  dataset = "landcover",
  landmass = "L48",
  extraction.dir = file.path(tempdir(), "FedData", "extractions", "nlcd", label),
  raster.options = c("COMPRESS=DEFLATE", "ZLEVEL=9"),
  force.redo = FALSE
)

nlcd_colors()

pal_nlcd()
```

## Arguments

- template:

  An
  [`Simple Feature`](https://r-spatial.github.io/sf/reference/sf.html)
  or
  [`terra`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object to serve as a template for cropping.

- label:

  A character string naming the study area.

- year:

  An integer representing the year of desired NLCD product. Acceptable
  values are 2019 (default), 2016, 2011, 2008, 2006, 2004, and 2001. The
  L48 data set for 2021 is corrupted on the NLCD Mapserver, and is thus
  not available through FedData.

- dataset:

  A character string representing type of the NLCD product. Acceptable
  values are 'landcover' (default), 'impervious', and 'canopy'.

- landmass:

  A character string representing the landmass to be extracted
  Acceptable values are 'L48' (lower 48 US states, the default), 'AK'
  (Alaska, 2001, 2011 and 2016 only), 'HI' (Hawaii, 2001 only), and 'PR'
  (Puerto Rico, 2001 only).

- extraction.dir:

  A character string indicating where the extracted and cropped NLCD
  data should be put. The directory will be created if missing.

- raster.options:

  a vector of GDAL options passed to
  [terra::writeRaster](https://rspatial.github.io/terra/reference/writeRaster.html).

- force.redo:

  If an extraction for this template and label already exists, should a
  new one be created?

## Value

A `RasterLayer` cropped to the bounding box of the template.

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract data for the Mesa Verde National Park:

# Get the NLCD (USA ONLY)
# Returns a raster
NLCD <-
  get_nlcd(
    template = FedData::meve,
    label = "meve",
    year = 2016
  )

# Plot with terra::plot
terra::plot(NLCD)
} # }
```
