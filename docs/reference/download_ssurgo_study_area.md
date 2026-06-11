# Download a zipped directory containing the spatial and tabular data for a SSURGO study area.

`download_ssurgo_study_area` first tries to download data including a
state-specific Access template, then the general US template.

## Usage

``` r
download_ssurgo_study_area(area, date, raw.dir)
```

## Arguments

- area:

  A character string indicating the SSURGO study area to be downloaded.

- date:

  A character string indicating the date of the most recent update to
  the SSURGO area for these data. This information may be gleaned from
  the SSURGO Inventory
  ([`get_ssurgo_inventory`](https://docs.ropensci.org/FedData/reference/get_ssurgo_inventory.md)).

- raw.dir:

  A character string indicating where raw downloaded files should be
  put.

## Value

A character string representing the full local path of the SSURGO study
areas zipped directory.
