# Download and extract the daily data for a GHCN weather station.

`get_ghcn_daily_station` returns a named list of
[`data.frame`](https://rdrr.io/r/base/data.frame.html)`s`, one for each
`elements`. If `elements` is undefined, it returns all available weather
tables for the station

## Usage

``` r
get_ghcn_daily_station(
  ID,
  elements = NULL,
  years = NULL,
  raw.dir,
  standardize = FALSE,
  force.redo = FALSE
)
```

## Arguments

- ID:

  A character string giving the station ID.

- elements:

  A character vector of elements to extract.  
  The five core elements are:  
  PRCP = Precipitation (tenths of mm)  
  SNOW = Snowfall (mm)  
  SNWD = Snow depth (mm)  
  TMAX = Maximum temperature (tenths of degrees C)  
  TMIN = Minimum temperature (tenths of degrees C)  
    
  The other elements are:  

  ACMC = Average cloudiness midnight to midnight from 30-second
  ceilometer data (percent)  
  ACMH = Average cloudiness midnight to midnight from manual
  observations (percent)  
  ACSC = Average cloudiness sunrise to sunset from 30-second ceilometer
  data (percent)  
  ACSH = Average cloudiness sunrise to sunset from manual observations
  (percent)  
  AWDR = Average daily wind direction (degrees)  
  AWND = Average daily wind speed (tenths of meters per second)  
  DAEV = Number of days included in the multiday evaporation total
  (MDEV)  
  DAPR = Number of days included in the multiday precipitation total
  (MDPR)  
  DASF = Number of days included in the multiday snowfall total (MDSF)  
  DATN = Number of days included in the multiday minimum temperature
  (MDTN)  
  DATX = Number of days included in the multiday maximum temperature
  (MDTX)  
  DAWM = Number of days included in the multiday wind movement (MDWM)  
  DWPR = Number of days with non-zero precipitation included in multiday
  precipitation total (MDPR)  
  EVAP = Evaporation of water from evaporation pan (tenths of mm)  
  FMTM = Time of fastest mile or fastest 1-minute wind (hours and
  minutes, i.e., HHMM)  
  FRGB = Base of frozen ground layer (cm)  
  FRGT = Top of frozen ground layer (cm)  
  FRTH = Thickness of frozen ground layer (cm)  
  GAHT = Difference between river and gauge height (cm)  
  MDEV = Multiday evaporation total (tenths of mm; use with DAEV)  
  MDPR = Multiday precipitation total (tenths of mm; use with DAPR and
  DWPR, if available)  
  MDSF = Multiday snowfall total  
  MDTN = Multiday minimum temperature (tenths of degrees C; use with
  DATN)  
  MDTX = Multiday maximum temperature (tenths of degrees C; use with
  DATX)  
  MDWM = Multiday wind movement (km)  
  MNPN = Daily minimum temperature of water in an evaporation pan
  (tenths of degrees C)  
  MXPN = Daily maximum temperature of water in an evaporation pan
  (tenths of degrees C)  
  PGTM = Peak gust time (hours and minutes, i.e., HHMM)  
  PSUN = Daily percent of possible sunshine (percent)  
  SN\*# = Minimum soil temperature (tenths of degrees C) where \*
  corresponds to a code for ground cover and \# corresponds to a code
  for soil depth.  
    
  Ground cover codes include the following:  
  0 = unknown  
  1 = grass  
  2 = fallow  
  3 = bare ground  
  4 = brome grass  
  5 = sod  
  6 = straw multch  
  7 = grass muck  
  8 = bare muck  
    
  Depth codes include the following:  
  1 = 5 cm  
  2 = 10 cm  
  3 = 20 cm  
  4 = 50 cm  
  5 = 100 cm  
  6 = 150 cm  
  7 = 180 cm  
    
  SX\*# = Maximum soil temperature (tenths of degrees C) where \*
  corresponds to a code for ground cover and \# corresponds to a code
  for soil depth.  
  See SN\*# for ground cover and depth codes.  
  TAVG = Average temperature (tenths of degrees C) (Note that TAVG from
  source 'S' corresponds to an average for the period ending at 2400 UTC
  rather than local midnight)  
  THIC = Thickness of ice on water (tenths of mm)  
  TOBS = Temperature at the time of observation (tenths of degrees C)  
  TSUN = Daily total sunshine (minutes)  
  WDF1 = Direction of fastest 1-minute wind (degrees)  
  WDF2 = Direction of fastest 2-minute wind (degrees)  
  WDF5 = Direction of fastest 5-second wind (degrees)  
  WDFG = Direction of peak wind gust (degrees)  
  WDFI = Direction of highest instantaneous wind (degrees)  
  WDFM = Fastest mile wind direction (degrees)  
  WDMV = 24-hour wind movement (km)  
  WESD = Water equivalent of snow on the ground (tenths of mm)  
  WESF = Water equivalent of snowfall (tenths of mm)  
  WSF1 = Fastest 1-minute wind speed (tenths of meters per second)  
  WSF2 = Fastest 2-minute wind speed (tenths of meters per second)  
  WSF5 = Fastest 5-second wind speed (tenths of meters per second)  
  WSFG = Peak gust wind speed (tenths of meters per second)  
  WSFI = Highest instantaneous wind speed (tenths of meters per
  second)  
  WSFM = Fastest mile wind speed (tenths of meters per second)  
  WT\*\* = Weather Type where \*\* has one of the following values:  
    
  01 = Fog, ice fog, or freezing fog (may include heavy fog)  
  02 = Heavy fog or heaving freezing fog (not always distinguished from
  fog)  
  03 = Thunder  
  04 = Ice pellets, sleet, snow pellets, or small hail  
  05 = Hail (may include small hail)  
  06 = Glaze or rime  
  07 = Dust, volcanic ash, blowing dust, blowing sand, or blowing
  obstruction  
  08 = Smoke or haze  
  09 = Blowing or drifting snow  
  10 = Tornado, waterspout, or funnel cloud  
  11 = High or damaging winds  
  12 = Blowing spray  
  13 = Mist  
  14 = Drizzle  
  15 = Freezing drizzle  
  16 = Rain (may include freezing rain, drizzle, and freezing drizzle)  
  17 = Freezing rain  
  18 = Snow, snow pellets, snow grains, or ice crystals  
  19 = Unknown source of precipitation  
  21 = Ground fog  
  22 = Ice fog or freezing fog  
    
  WV\*\* = Weather in the Vicinity where \*\* has one of the following
  values:  
  01 = Fog, ice fog, or freezing fog (may include heavy fog)  
  03 = Thunder  
  07 = Ash, dust, sand, or other blowing obstruction  
  18 = Snow or ice crystals  
  20 = Rain or snow shower

- years:

  A numeric vector indicating which years to get.

- raw.dir:

  A character string indicating where raw downloaded files should be
  put.

- standardize:

  Select only common year/month/day? Defaults to FALSE.

- force.redo:

  If this weather station has been downloaded before, should it be
  updated? Defaults to FALSE.

## Value

A named list of
[`data.frame`](https://rdrr.io/r/base/data.frame.html)`s`, one for each
`elements`.
