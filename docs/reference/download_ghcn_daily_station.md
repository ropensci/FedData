# Download the daily data for a GHCN weather station.

Download the daily data for a GHCN weather station.

## Usage

``` r
download_ghcn_daily_station(ID, raw.dir, force.redo = FALSE)
```

## Arguments

- ID:

  A character string giving the station ID.

- raw.dir:

  A character string indicating where raw downloaded files should be
  put.

- force.redo:

  If this weather station has been downloaded before, should it be
  updated? Defaults to FALSE.

## Value

A character string representing the full local path of the GHCN station
data.
