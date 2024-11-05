---
editor_options: 
  markdown: 
    wrap: 72
---

# CRAN comments

This is a minor release from FedData v4.0.1 to FedData v4.1.0. Please
see NEWS.md for release details.

FedData was removed from CRAN due to the removal of an dependency. That
package has since been restored.

Please note that in recent CRAN submissions, CRAN maintainers have been 
seeing URL check issues in checks by servers hosted outside the USA. This is
likely due to US Government servers not "liking" automated checking from
servers outside the USA.

------------------------------------------------------------------------

## R CMD check results

`devtools::check()` result:

**Test environment:** aarch64-apple-darwin20 (64-bit), R 4.4.2 (2024-10-31)

**0 errors ✔ \| 0 warnings ✔ \| 0 notes ✔**

------------------------------------------------------------------------

## GitHub Actions results:

**Test environments:**

-   macOS-latest, R release
-   windows-latest, R release
-   ubuntu-latest, R devel
-   ubuntu-latest, R release
-   ubuntu-latest, R oldrel-1

**0 errors ✔ \| 0 warnings ✔ \| 0 notes ✔**

------------------------------------------------------------------------

## R-hub check results:

**Test environments:**

-   ubuntu-release
-   linux (R-devel)
-   macos (R-devel)
-   macos-arm64 (R-devel)
-   windows (R-devel)

**0 errors ✔ \| 0 warnings ✔ \| 0 notes ✔**

------------------------------------------------------------------------

## win-builder results

`devtools::check_win_devel()` result:

**Test environment:** Windows Server 2022, R-devel, 64 bit

**0 errors ✔ \| 0 warnings ✔ \| 1 note ✖**

- NOTE: Found the following (possibly) invalid URLs:
  - URL: https://www.ncei.noaa.gov/products/land-based-station/global-historical-climatology-network-daily
    - From: README.md
    - Status: Error
    - Message: Empty reply from server
  - URL: https://www.ncei.noaa.gov/products/paleoclimatology/tree-ring
    - From: README.md
    - Status: Error
    - Message: Empty reply from server
- Explanation: Manually accessed each URL, 
and confirmed a http 200 status for each. Possibly due to check 
server outside USA trying to access USA .gov URLs.

------------------------------------------------------------------------

## revdepcheck results

We checked 0 reverse dependencies, comparing R CMD check results 
across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

