.onAttach <- function(libname, pkgname) {
  msg <-
    "You have loaded FedData v4.0.0,
the latest major version of FedData. We have retired
FedData dependencies on the `sp` and `raster` packages.
All functions in FedData v4 return `terra` (raster)
or `sf` (vector) objects by default, and there may be
other breaking changes."

  base::packageStartupMessage(msg)
}
