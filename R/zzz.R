.onAttach <- function(libname, pkgname) {
  msg <-
    "You have loaded FedData v4.
As of FedData v4 we have retired
dependencies on the `sp` and `raster` packages.
All functions in FedData v4 return `terra` (raster)
or `sf` (vector) objects by default, and there may be
other breaking changes."

  base::packageStartupMessage(msg)
}
