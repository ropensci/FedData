---
editor_options: 
  markdown: 
    wrap: 72
---

# CRAN comments

This is a patch release from FedData v3.0.1 to FedData v3.0.2. Please
see NEWS.md for release details. 

Resubmission from 2023-02-23. Fixed by updating URLs:
Found the following (possibly) invalid URLs:
   URL: https://opensource.org/licenses/MIT (moved to https://opensource.org/license/mit/)
     From: README.md
     Status: 301
     Message: Moved Permanently

------------------------------------------------------------------------

## R CMD check results

`devtools::check()` result:

**Test environment:** local MacOS Version 13.2.1 install, R 4.2.2

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

`rhub::check_for_cran()` result:

**Test environment:** Windows Server 2022, R-devel, 64 bit

**0 errors ✔ \| 0 warnings ✔ \| 1 note ✖**

-   checking for detritus in the temp directory ...
    -   NOTE: Found the following files/directories:
        'lastMiKTeXException'
    -   Explanation: As noted in [R-hub issue
        #503](https://github.com/r-hub/rhub/issues/503), this could be
        due to a bug/crash in MiKTeX and can likely be ignored.

------------------------------------------------------------------------

`rhub::check_for_cran()` result:

**Test environment:** Ubuntu Linux 20.04.1 LTS, R-release, GCC

**0 errors ✔ \| 0 warnings ✔ \| 0 notes ✖**

------------------------------------------------------------------------

`rhub::check_for_cran()` result:

**Test environment:** Fedora Linux, R-devel, clang, gfortran

**0 errors ✔ \| 0 warnings ✔ \| 1 note ✖**

-   checking HTML version of manual ...
    -   NOTE: Skipping checking HTML validation: no command 'tidy' found
    -   Explanation: The note about HTML validation only occurs on the
        Fedora dev install and does not seem critical. The HTML version
        of the manual is able to be validated on other platforms.

------------------------------------------------------------------------

## win-builder results

`devtools::check_win_devel()` result:

**Test environment:** x86_64-w64-mingw32 (64-bit)

**0 errors ✔ \| 0 warnings ✔ \| 1 note ✖**

- NOTE: Found the following (possibly) invalid URLs: URL:
    <https://www.ncei.noaa.gov/products/land-based-station/global-historical-climatology-network-daily>
    From: README.md Status: Error Message: Empty reply from server URL:
    <https://www.ncei.noaa.gov/products/paleoclimatology/tree-ring>
    From: README.md Status: Error Message: Empty reply from server
- Explanation: Manually accessed each URL, 
and confirmed a http 200 status for each.

------------------------------------------------------------------------

## revdepcheck results

We checked one reverse dependencies, comparing R CMD check results
across CRAN and dev versions of this package.

-   We saw 0 new problems
-   We failed to check 0 packages
