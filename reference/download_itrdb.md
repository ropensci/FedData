# Download the latest version of the ITRDB.

Downloads and parses the latest zipped (numbered) version of the ITRDB.
This function includes improvements to the
[`read_crn`](https://docs.ropensci.org/FedData/reference/read_crn.md)
function from the dplR library. The principle changes are better parsing
of metadata, and support for the Schweingruber-type Tucson format.
Chronologies that are unable to be read are reported to the user.

## Usage

``` r
download_itrdb(
  raw.dir = paste0(tempdir(), "/FedData/raw/itrdb"),
  force.redo = FALSE
)
```

## Arguments

- raw.dir:

  A character string indicating where raw downloaded files should be
  put. The directory will be created if missing. Defaults to
  './RAW/ITRDB/'.

- force.redo:

  If a download already exists, should a new one be created? Defaults to
  FALSE.

## Value

A data frame containing all of the ITRDB data.
