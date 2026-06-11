# Download the latest version of the ITRDB, and extract given parameters.

`get_itrdb` returns a named list of length 3:

1.  'metadata': A data frame or
    [`Simple Feature`](https://r-spatial.github.io/sf/reference/sf.html)
    (if `makeSpatial==TRUE`) of the locations and names of extracted
    ITRDB chronologies,

2.  'widths': A matrix of tree-ring widths/densities given user
    selection, and

3.  'depths': A matrix of tree-ring sample depths.

## Usage

``` r
get_itrdb(
  template = NULL,
  label = NULL,
  recon.years = NULL,
  calib.years = NULL,
  species = NULL,
  measurement.type = NULL,
  chronology.type = NULL,
  raw.dir = paste0(tempdir(), "/FedData/raw/itrdb"),
  extraction.dir = ifelse(!is.null(label), paste0(tempdir(),
    "/FedData/extractions/itrdb/", label, "/"), paste0(tempdir(),
    "/FedData/extractions/itrdb")),
  force.redo = FALSE
)
```

## Arguments

- template:

  An
  [`Simple Feature`](https://r-spatial.github.io/sf/reference/sf.html)
  or
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object to serve as a template for cropping. If missing, all available
  global chronologies are returned.

- label:

  A character string naming the study area.

- recon.years:

  A numeric vector of years over which reconstructions are needed; if
  missing, the union of all years in the available chronologies are
  given.

- calib.years:

  A numeric vector of all required years—chronologies without these
  years will be discarded; if missing, all available chronologies are
  given.

- species:

  A character vector of 4-letter tree species identifiers; if missing,
  all available chronologies are given.

- measurement.type:

  A character vector of measurement type identifiers. Options include:

  - 'Total Ring Density'

  - 'Earlywood Width'

  - 'Earlywood Density'

  - 'Latewood Width'

  - 'Minimum Density'

  - 'Ring Width'

  - 'Latewood Density'

  - 'Maximum Density'

  - 'Latewood Percent'

  if missing, all available chronologies are given.

- chronology.type:

  A character vector of chronology type identifiers. Options include:

  - 'ARSTND'

  - 'Low Pass Filter'

  - 'Residual'

  - 'Standard'

  - 'Re-Whitened Residual'

  - 'Measurements Only'

  if missing, all available chronologies are given.

- raw.dir:

  A character string indicating where raw downloaded files should be
  put. The directory will be created if missing.

- extraction.dir:

  A character string indicating where the extracted and cropped ITRDB
  dataset should be put. The directory will be created if missing.

- force.redo:

  If an extraction already exists, should a new one be created? Defaults
  to FALSE.

## Value

A named list containing the 'metadata', 'widths', and 'depths' data.

## Examples

``` r
if (FALSE) { # \dontrun{
# Get the ITRDB records
ITRDB <- get_itrdb(
  template = FedData::meve,
  label = "meve"
)

# Plot the VEP polygon
plot(meve)

# Map the locations of the tree ring chronologies
plot(ITRDB$metadata$geometry, pch = 1, add = TRUE)
legend("bottomleft", pch = 1, legend = "ITRDB chronologies")
} # }
```
