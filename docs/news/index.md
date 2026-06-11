# Changelog

## FedData 4.0.1

CRAN release: 2024-03-16

- Updated the getting started article to not include a degree symbol in
  the Daymet graph, which caused compilation errors in Windows. Fixes
  [Issue](https://github.com/ropensci/FedData/issues/106)
  [\#106](https://github.com/ropensci/FedData/issues/106).
- Updated
  [`get_nlcd()`](https://docs.ropensci.org/FedData/reference/get_nlcd.md)
  to include the 2021 NLCD as the default, in response to
  [Issue](https://github.com/ropensci/FedData/issues/105)
  [\#105](https://github.com/ropensci/FedData/issues/105).
- Updated outdated package description
- Bumped GDAL version req to \>= 3.1.0 to accommodate storing spatial
  vectors as FlatGeoBufs
- Added `arcgislayers` dependency and retired self-written esri
  functions. Closes
  [Issue](https://github.com/ropensci/FedData/issues/109)
  [\#109](https://github.com/ropensci/FedData/issues/109).

## FedData 4.0.0

CRAN release: 2023-10-03

- Updated the [README](https://docs.ropensci.org/FedData/news/README.md)
  and moved examples to an article
- Added the [PAD-US
  dataset](https://www.usgs.gov/programs/gap-analysis-project/science/pad-us-data-overview)
  in response to [Issue](https://github.com/ropensci/FedData/issues/100)
  [\#100](https://github.com/ropensci/FedData/issues/100).
- Updated
  [`get_nlcd()`](https://docs.ropensci.org/FedData/reference/get_nlcd.md)
  to fix bug [Issue](https://github.com/ropensci/FedData/issues/101)
  [\#101](https://github.com/ropensci/FedData/issues/101) in tree canopy
  data, whose URL changed.
- Updated
  [`get_ghcn_daily()`](https://docs.ropensci.org/FedData/reference/get_ghcn_daily.md)
  to respond to slightly different file formatting on https server
  [Issue](https://github.com/ropensci/FedData/issues/99)
  [\#99](https://github.com/ropensci/FedData/issues/99).
- Removed dependencies on `sp` and `raster`. All functions now return
  either `sf` or `terra` objects.

## FedData 3.0.4

CRAN release: 2023-05-25

- Updated
  [`get_ghcn_daily()`](https://docs.ropensci.org/FedData/reference/get_ghcn_daily.md)
  and related functions to use https addresses instead of the legacy ftp
  server [Issue](https://github.com/ropensci/FedData/issues/99)
  [\#99](https://github.com/ropensci/FedData/issues/99).
- Updated
  [`get_nhd()`](https://docs.ropensci.org/FedData/reference/get_nhd.md)
  to gracefully handle situations where no point data (or any NHD data)
  are available within an area of interest, using solution offered in
  [Issue](https://github.com/ropensci/FedData/issues/98)
  [\#98](https://github.com/ropensci/FedData/issues/98).
- Updated
  [`get_nhd()`](https://docs.ropensci.org/FedData/reference/get_nhd.md)
  to handle some strange geometries present in the NHD
  [Issue](https://github.com/ropensci/FedData/issues/98)
  [\#98](https://github.com/ropensci/FedData/issues/98).

## FedData 3.0.3

CRAN release: 2023-03-10

- Updated
  [`get_nhd()`](https://docs.ropensci.org/FedData/reference/get_nhd.md)
  to correctly access the extraction.dir if it already exists
  [Issue](https://github.com/ropensci/FedData/issues/95)
  [\#95](https://github.com/ropensci/FedData/issues/95).
  [`get_nhd()`](https://docs.ropensci.org/FedData/reference/get_nhd.md)
  and
  [`get_ned()`](https://docs.ropensci.org/FedData/reference/get_ned.md)
  now both use [`file.path()`](https://rdrr.io/r/base/file.path.html) in
  constructing file outputs.
- Updated
  [`get_daymet()`](https://docs.ropensci.org/FedData/reference/get_daymet.md)
  to use `terra` for writing raster data, in response to
  [Issue](https://github.com/ropensci/FedData/issues/95)
  [\#96](https://github.com/ropensci/FedData/issues/96). All functions
  in `FedData` now use `terra` for writing.
- Cleaned up documentation.

## FedData 3.0.2

CRAN release: 2023-02-26

- Removed `rgdal` dependency throughout and instead require `raster` \>=
  3.6.3.
- Switched back to using the Web Coverave Service for the NLCD, as it is
  now providing data in the native CRS.

## FedData 3.0.1

CRAN release: 2022-11-28

- Updated
  [`get_nlcd()`](https://docs.ropensci.org/FedData/reference/get_nlcd.md)
  to read CRS data correctly
  [Issue](https://github.com/ropensci/FedData/issues/91)
  [\#91](https://github.com/ropensci/FedData/issues/91) by requiring
  `rgdal` be installed.
- Changed Daymet tempo codes
  [Issue](https://github.com/ropensci/FedData/issues/92)
  [\#92](https://github.com/ropensci/FedData/issues/92)

## FedData 3.0.0

CRAN release: 2022-10-11

- All `get_*()` functions now return `sf` or `raster` objects.
- Changed `SDA_query()` to
  [`soils_query()`](https://docs.ropensci.org/FedData/reference/soils_query.md)
  to avoid namespace masking with SoilDB.
- Added
  [`get_nass_cdl()`](https://docs.ropensci.org/FedData/reference/get_nass_cdl.md)
  to retrieve the NASS Cropland Data Layer
- Updated
  [`get_daymet()`](https://docs.ropensci.org/FedData/reference/get_daymet.md)
  to pull from ORNL WCS, and fixed bug in
  [Issue](https://github.com/ropensci/FedData/issues/49)
  [\#49](https://github.com/ropensci/FedData/issues/49)
- Fixed issue where
  [`soils_query()`](https://docs.ropensci.org/FedData/reference/soils_query.md)
  was only returning first SSURGO study area
  ([Issue](https://github.com/ropensci/FedData/issues/51)
  [\#51](https://github.com/ropensci/FedData/issues/51))
- Fixed issue where date parsing was USA-specific by leveraging
  [`lubridate::parse_date_time`](https://lubridate.tidyverse.org/reference/parse_date_time.html)
  ([Issue](https://github.com/ropensci/FedData/issues/61)
  [\#61](https://github.com/ropensci/FedData/issues/61))
- [`get_ssurgo()`](https://docs.ropensci.org/FedData/reference/get_ssurgo.md)
  now saves in the [GeoPackage file format](http://www.geopackage.org)
- [`get_nhd()`](https://docs.ropensci.org/FedData/reference/get_nhd.md)
  now access ESRI web services. Users can optionally get data from
  NHDPlus.
- [`get_nlcd()`](https://docs.ropensci.org/FedData/reference/get_nlcd.md)
  now provides data in native CRS (CONUS Albers), rather than
  web-mercator ([Issue](https://github.com/ropensci/FedData/issues/77)
  [\#77](https://github.com/ropensci/FedData/issues/77)) by pulling from
  self-hosted Cloud-Optimized GeoTiffs.
- [Issue](https://github.com/ropensci/FedData/issues/88)
  [\#88](https://github.com/ropensci/FedData/issues/88) identified an
  issue with GDAL versions prior GDAL 3 reading COGs using the vsicurl
  functionality in GDAL. Accordingly, FedData3 requires GDAL \>= 3.0.0.
- Added Mesa Verde National Park as exemplar region, and removed use of
  `paleocar` package.
- [`get_ned()`](https://docs.ropensci.org/FedData/reference/get_ned.md)
  now pulls from USGS NED Cloud-Optimized GeoTiffs available at
  <https://prd-tnm.s3.amazonaws.com/index.html?prefix=StagedProducts/Elevation/>.

## FedData 2.5.7

CRAN release: 2019-04-22

- Removing many internet resource tests from CRAN, to satisfy: ‘Packages
  which use Internet resources should fail gracefully with an
  informative message if the resource is not available (and not give a
  check warning nor error).’

## FedData 2.5.6

CRAN release: 2019-01-11

- Built-in access to the Soils Data Analysis query service to remove
  dependency on soilDB package.

## FedData 2.5.5

CRAN release: 2018-08-09

- Fixed issue ([\#41](https://github.com/ropensci/FedData/issues/41))
  that occurs when mosaicking NLCD tiles that are not cropped. When they
  aren’t cropped, the NLCD data is never read into memory, and the
  temporary file that the raster was created from gets destroyed.
  Solution: Force NLCD data into memory prior to mosaicking.
- Added (non-CRAN) test for issue
  [\#41](https://github.com/ropensci/FedData/issues/41)

## FedData 2.5.4

CRAN release: 2018-07-05

- Fixed issue in downloading NED tiles.

## FedData 2.5.3

CRAN release: 2018-05-21

- Added httr to package imports.

## FedData 2.5.2

CRAN release: 2018-03-12

- Updated NHD HUC4 to copy stored on Github.
- Fixed bug in ITRDB that caused some chronologies not to be read.

## FedData 2.5.1

CRAN release: 2018-01-24

- Switch to laze-loading data.
- Updated NHD paths to new National Map directory structure.

## FedData 2.5.0

CRAN release: 2017-12-15

- Added functions for the National Land Cover Database.

## FedData 2.4.7

CRAN release: 2017-11-20

- SSURGO fixed test where supplying an unavailable survey area now
  returns NULL instead of an error.
- SSURGO zip directory encoding changes as of late October 2017 forced
  changes in the FedData:::get_ssurgo_study_area function.
- Fixed issue where NHD template wouldn’t load because they added a jpeg
  preview to the directory.

## FedData 2.4.6

CRAN release: 2017-08-18

- DAYMET functions now do *not* operate in parallel. This was breaking
  the download functions.
- Final update for version 2 of FedData.
- Accepted to ROpenSci! Migrating to the ROpenSci organization on
  GitHub.

## FedData 2.4.3

- writeOGR for SSURGO and NHD were failing on Windows when the
  `extraction.dir` included a trailing slash. Paths are now normalized
  to remove the trailing slash.

## FedData 2.4.2

- Updated the `get_ned` function to provide more useful errors and
  warnings when downloads are unsuccessful.

## FedData 2.4.1

- Added pkgdown site.
- SSURGO functions (e.g., `get_ssurgo`) now doesn’t bomb on large (\> 1
  billion sq meter) requests. Now, the area of interest is broken into
  smaller chunks to build the download list.

## FedData 2.4.0

CRAN release: 2017-01-20

- Added a `NEWS.md` file to track changes to the package.
- Updated DAYMET functions to fix a bug that downloaded only one tile at
  a time.
- Linted all code.
