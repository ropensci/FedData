# Extract data from a SSURGO database pertaining to a set of mapunits.

`extract_ssurgo_data` creates a directed graph of the joins in a SSURGO
tabular dataset, and then iterates through the tables, only retaining
data pertinent to a set of mapunits.

## Usage

``` r
extract_ssurgo_data(tables, mapunits)
```

## Arguments

- tables:

  A list of SSURGO tabular data.

- mapunits:

  A character vector of mapunits (likely dropped from SSURGO spatial
  data) defining which mapunits to retain.

## Value

A list of extracted SSURGO tabular data.
