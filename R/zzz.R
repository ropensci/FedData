.onAttach <- function(libname, pkgname) {
  msg <-
    paste0(
      "\nThis is version ", utils::packageVersion(pkgname),
      " of ", pkgname, ".\n",
      "On August 1, 2023, we will release FedData v4 on CRAN,
the next major version of FedData. We will be retiring
FedData dependencies on the `sp` and `raster` packages.
All functions in FedData v4 will return `terra` (raster)
or `sf` (vector) objects by default, and there may be
other breaking changes."
    )

  base::packageStartupMessage(msg)
}
