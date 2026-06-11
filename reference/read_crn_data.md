# Read chronology data from a Tucson-format chronology file.

This function includes improvements to the
[`read_crn`](https://docs.ropensci.org/FedData/reference/read_crn.md)
function from the dplR library. The principle changes are better parsing
of metadata, and support for the Schweingruber-type Tucson format.
Chronologies that are unable to be read are reported to the user. The
user (or
[`read_crn`](https://docs.ropensci.org/FedData/reference/read_crn.md))
must tell the function whether the file is a Schweingruber-type
chronology.

## Usage

``` r
read_crn_data(file, SCHWEINGRUBER)
```

## Arguments

- file:

  A character string path pointing to a `*.crn` file to be read.

- SCHWEINGRUBER:

  Is the file in the Schweingruber-type Tucson format?

## Value

A data.frame containing the data, or if `SCHWEINGRUBER==T`, a list
containing four types of data.
