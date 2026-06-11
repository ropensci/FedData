# Read a Tucson-format chronology file.

This function includes improvements to the `read.crn` function from the
dplR library. The principle changes are better parsing of metadata, and
support for the Schweingruber-type Tucson format. Chronologies that are
unable to be read are reported to the user. This function automatically
recognizes Schweingruber-type files.

## Usage

``` r
read_crn(file)
```

## Arguments

- file:

  A character string path pointing to a `*.crn` file to be read.

## Value

A list containing the metadata and chronology.

## Details

This wraps two other functions:
[`read_crn_metadata`](https://docs.ropensci.org/FedData/reference/read_crn_metadata.md)
[`read_crn_data`](https://docs.ropensci.org/FedData/reference/read_crn_data.md).
