# Download the 1-km DAYMET daily weather dataset for a region as a netcdf.

Data are downloaded in the NetCDF format. `download_daymet_thredds`
returns the path to the downloaded NetCDF file.

## Usage

``` r
download_daymet_thredds(bbox, element, year, region, tempo)
```

## Arguments

- bbox:

  the bounding box in WGS84 coordinates as a comma-separated character
  vector "xmin,ymin,xmax,ymax"

- element:

  An element to extract.  
  The available elements are:  
  dayl = Duration of the daylight period in seconds per day. This
  calculation is based on the period of the day during which the sun is
  above a hypothetical flat horizon.  
  prcp = Daily total precipitation in millimeters per day, sum of all
  forms converted to water-equivalent. Precipitation occurrence on any
  given day may be ascertained.  
  srad = Incident shortwave radiation flux density in watts per square
  meter, taken as an average over the daylight period of the day. NOTE:
  Daily total radiation (MJ/m2/day) can be calculated as follows: ((srad
  (W/m2) \* dayl (s/day)) / l,000,000)  
  swe = Snow water equivalent in kilograms per square meter. The amount
  of water contained within the snowpack.  
  tmax = Daily maximum 2-meter air temperature in degrees Celsius.  
  tmin = Daily minimum 2-meter air temperature in degrees Celsius.  
  vp = Water vapor pressure in pascals. Daily average partial pressure
  of water vapor.  

- year:

  An integer year to extract.

- region:

  The name of a region. The available regions are:  
  na = North America  
  hi = Hawaii  
  pr = Puerto Rico  

- tempo:

  The frequency of the data. The available tempos are:  
  day = Daily data  
  mon = Monthly summary data  
  ann = Annual summary data  

## Value

A named list of character vectors, each representing the full local
paths of the tile downloads.
