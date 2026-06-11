# Convert a list of station data to a single data frame.

`station_to_data_frame` returns a `data.frame` of the GHCN station data
list.

## Usage

``` r
station_to_data_frame(station.data)
```

## Arguments

- station.data:

  A named list containing station data

## Value

A `data.frame` of the containing the unwrapped station data

## Details

This function unwraps the station data and merges all data into a single
data frame, with the first column being in the `Date` class.
