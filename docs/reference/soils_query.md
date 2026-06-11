# Submit a Soil Data Access (SDA) Query

`soils_query` submit an SQL query to retrieve data from the Soil Data
Mart. Please see https://sdmdataaccess.sc.egov.usda.gov/Query.aspx for
guidelines

## Usage

``` r
soils_query(q)
```

## Arguments

- q:

  A character string representing a SQL query to the SDA service

## Value

A tibble returned from the SDA service
