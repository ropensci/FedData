# FedData (development version)
- Updated the getting started article to not include a degree symbol in the Daymet graph, which caused compilation errors in Windows. Fixes [Issue #106](https://github.com/ropensci/FedData/issues/106).
- Updated `get_nlcd()` to include the 2021 NLCD as the default, in response to [Issue #105](https://github.com/ropensci/FedData/issues/105).
- Updated outdated package description
- Bumped GDAL version req to >= 3.1.0 to accommodate storing spatial vectors as FlatGeoBufs
- Added `arcgislayers` dependency and retired self-written esri functions. Closes [Issue #109](https://github.com/ropensci/FedData/issues/109).

# FedData 4.0.0
-   Updated the [README](README.md) and moved examples to an article
-   Added the [PAD-US dataset](https://www.usgs.gov/programs/gap-analysis-project/science/pad-us-data-overview)
in response to [Issue #100](https://github.com/ropensci/FedData/issues/100).
-   Updated `get_nlcd()` to fix bug [Issue #101](https://github.com/ropensci/FedData/issues/101) in tree canopy data, whose URL changed.
-   Updated `get_ghcn_daily()` to respond to slightly different file formatting on https server [Issue #99](https://github.com/ropensci/FedData/issues/99).
-   Removed dependencies on `sp` and `raster`. All functions now return either `sf` or `terra` objects.


# FedData 3.0.4

-   Updated `get_ghcn_daily()` and related functions to use https addresses 
instead of the legacy ftp server [Issue #99](https://github.com/ropensci/FedData/issues/99).
-   Updated `get_nhd()` to gracefully handle situations where no point data (or any NHD data) are available within an area of interest, using solution offered in [Issue #98](https://github.com/ropensci/FedData/issues/98).
-   Updated `get_nhd()` to handle some strange geometries present in the NHD [Issue #98](https://github.com/ropensci/FedData/issues/98).

# FedData 3.0.3

-   Updated `get_nhd()` to correctly access the extraction.dir if it already exists [Issue #95](https://github.com/ropensci/FedData/issues/95). `get_nhd()` and `get_ned()`
now both use `file.path()` in constructing file outputs.
-   Updated `get_daymet()` to use `terra` for writing raster data, in response to
[Issue #96](https://github.com/ropensci/FedData/issues/95). All functions in 
`FedData` now use `terra` for writing.
-   Cleaned up documentation.

# FedData 3.0.2

-   Removed `rgdal` dependency throughout and instead require `raster` >= 3.6.3.
-   Switched back to using the Web Coverave Service for the NLCD, as it is now providing data in the native CRS.

# FedData 3.0.1

-   Updated `get_nlcd()` to read CRS data correctly [Issue #91](https://github.com/ropensci/FedData/issues/91) by requiring `rgdal` be installed.
-   Changed Daymet tempo codes [Issue #92](https://github.com/ropensci/FedData/issues/92)

# FedData 3.0.0

-   All `get_*()` functions now return `sf` or `raster` objects.
-   Changed `SDA_query()` to `soils_query()` to avoid namespace masking with SoilDB.
-   Added `get_nass_cdl()` to retrieve the NASS Cropland Data Layer
-   Updated `get_daymet()` to pull from ORNL WCS, and fixed bug in [Issue #49](https://github.com/ropensci/FedData/issues/49)
-   Fixed issue where `soils_query()` was only returning first SSURGO study area ([Issue #51](https://github.com/ropensci/FedData/issues/51))
-   Fixed issue where date parsing was USA-specific by leveraging `lubridate::parse_date_time` ([Issue #61](https://github.com/ropensci/FedData/issues/61))
-   `get_ssurgo()` now saves in the [GeoPackage file format](http://www.geopackage.org)
-   `get_nhd()` now access ESRI web services. Users can optionally get data from NHDPlus.
-   `get_nlcd()` now provides data in native CRS (CONUS Albers), rather than web-mercator ([Issue #77](https://github.com/ropensci/FedData/issues/77)) by pulling from self-hosted Cloud-Optimized GeoTiffs.
-   [Issue #88](https://github.com/ropensci/FedData/issues/88) identified an issue with GDAL versions prior GDAL 3 reading COGs using the vsicurl functionality in GDAL. Accordingly, FedData3 requires GDAL \>= 3.0.0.
-   Added Mesa Verde National Park as exemplar region, and removed use of `paleocar` package.
-   `get_ned()` now pulls from USGS NED Cloud-Optimized GeoTiffs available at <https://prd-tnm.s3.amazonaws.com/index.html?prefix=StagedProducts/Elevation/>.

# FedData 2.5.7

-   Removing many internet resource tests from CRAN, to satisfy: 'Packages which use Internet resources should fail gracefully with an informative message if the resource is not available (and not give a check warning nor error).'

# FedData 2.5.6

-   Built-in access to the Soils Data Analysis query service to remove dependency on soilDB package.

# FedData 2.5.5

-   Fixed issue (#41) that occurs when mosaicking NLCD tiles that are not cropped. When they aren't cropped, the NLCD data is never read into memory, and the temporary file that the raster was created from gets destroyed. Solution: Force NLCD data into memory prior to mosaicking.
-   Added (non-CRAN) test for issue #41

# FedData 2.5.4

-   Fixed issue in downloading NED tiles.

# FedData 2.5.3

-   Added httr to package imports.

# FedData 2.5.2

-   Updated NHD HUC4 to copy stored on Github.
-   Fixed bug in ITRDB that caused some chronologies not to be read.

# FedData 2.5.1

-   Switch to laze-loading data.
-   Updated NHD paths to new National Map directory structure.

# FedData 2.5.0

-   Added functions for the National Land Cover Database.

# FedData 2.4.7

-   SSURGO fixed test where supplying an unavailable survey area now returns NULL instead of an error.
-   SSURGO zip directory encoding changes as of late October 2017 forced changes in the FedData:::get_ssurgo_study_area function.
-   Fixed issue where NHD template wouldn't load because they added a jpeg preview to the directory.

# FedData 2.4.6

-   DAYMET functions now do *not* operate in parallel. This was breaking the download functions.
-   Final update for version 2 of FedData.
-   Accepted to ROpenSci! Migrating to the ROpenSci organization on GitHub.

# FedData 2.4.3

-   writeOGR for SSURGO and NHD were failing on Windows when the `extraction.dir` included a trailing slash. Paths are now normalized to remove the trailing slash.

# FedData 2.4.2

-   Updated the `get_ned` function to provide more useful errors and warnings when downloads are unsuccessful.

# FedData 2.4.1

-   Added pkgdown site.
-   SSURGO functions (e.g., `get_ssurgo`) now doesn't bomb on large (\> 1 billion sq meter) requests. Now, the area of interest is broken into smaller chunks to build the download list.

# FedData 2.4.0

-   Added a `NEWS.md` file to track changes to the package.
-   Updated DAYMET functions to fix a bug that downloaded only one tile at a time.
-   Linted all code.
