# FedData 2.5.2
* Updated NHD HUC4 to copy stored on Github.
* Fixed bug in ITRDB that caused some chronologies not to be read.

# FedData 2.5.1
* Switch to laze-loading data.
* Updated NHD paths to new National Map directory structure.

# FedData 2.5.0
* Added functions for the National Land Cover Database.

# FedData 2.4.7
* SSURGO fixed test where supplying an unavailable survey area now returns NULL instead of an error.
* SSURGO zip directory encoding changes as of late October 2017 forced changes in the FedData:::get_ssurgo_study_area function.
* Fixed issue where NHD template wouldn't load because they added a jpeg preview to the directory.

# FedData 2.4.6
* DAYMET functions now do *not* operate in parallel. This was breaking the download functions.
* Final update for version 2 of FedData.
* Accepted to ROpenSci! Migrating to the ROpenSci organization on GitHub.

# FedData 2.4.3
* writeOGR for SSURGO and NHD were failing on Windows when the `extraction.dir` included a trailing slash. Paths are now normalized to remove the trailing slash.

# FedData 2.4.2
* Updated the `get_ned` function to provide more useful errors and warnings when downloads are unsuccessful.

# FedData 2.4.1
* Added pkgdown site.
* SSURGO functions (e.g., `get_ssurgo`) now doesn't bomb on large (> 1 billion sq meter) requests. Now, the area of interest is broken into smaller chunks to build the download list.

# FedData 2.4.0
* Added a `NEWS.md` file to track changes to the package.
* Updated DAYMET functions to fix a bug that downloaded only one tile at a time.
* Linted all code.



