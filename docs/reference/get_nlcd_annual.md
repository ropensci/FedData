# Download and crop the Annual National Land Cover Database.

`get_nlcd_annual` returns a
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
of NLCD data cropped to a given template study area. The Annual NLCD is
currently only available for the conterminous United States. More
information about the Annual NLCD product is available on the [Annual
NLCD web page](https://www.mrlc.gov/data/project/annual-nlcd).

## Usage

``` r
get_nlcd_annual(
  template,
  label,
  year = 2024,
  product = "LndCov",
  region = "CU",
  collection = 1,
  version = 1,
  extraction.dir = file.path(tempdir(), "FedData", "extractions", "nlcd_annual", label),
  raster.options = c("COMPRESS=DEFLATE", "ZLEVEL=9"),
  force.redo = FALSE
)
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

  An integer vector representing the year of desired NLCD product.
  Acceptable values are currently 1985 through 2024 (defaults to 2024).

- product:

  A character vector representing type of the NLCD product. Defaults to
  'LndCov' (Land Cover).\
  LndCov = Land Cover\
  LndChg = Land Cover Change\
  LndCnf = Land Cover Confidence\
  FctImp = Fractional Impervious Surface\
  ImpDsc = Impervious Descriptor\
  SpcChg = Spectral Change Day of Year\

- region:

  A character string representing the region to be extracted Acceptable
  values are 'CU' (Conterminous US, the default), 'AK' (Alaska), and
  'HI' (Hawaii). **Currently, only 'CU' is available.**

- collection:

  An integer representing the collection number. **Currently, only '1'
  is available.**

- version:

  An integer representing the version number. **Currently, only '1' is
  available.**

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

## Details

Data are downloaded using the [USGS Web Coverage
Service](https://dmsdata.cr.usgs.gov/geoserver/web/) for the Annual
NLCD, in the native coordinate reference system and resolution (CONUS
Albers, EPSG:5070, at 30 m) and snapped to the native grid.

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract data for the Mesa Verde National Park:

# Get the NLCD (USA ONLY)
# Returns a raster
NLCD_ANNUAL <-
  get_nlcd_annual(
    template = FedData::meve,
    label = "meve",
    year = 2020,
    product =
      c(
        "LndCov",
        "LndChg",
        "LndCnf",
        "FctImp",
        "ImpDsc",
        "SpcChg"
      )
  )

NLCD_ANNUAL
} # }
```
