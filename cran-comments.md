---
editor_options: 
  markdown: 
    wrap: 72
---

# CRAN comments

This is a major release from FedData v3.0.4 to FedData v4.0.0. Please
see NEWS.md for release details. 

Please note that in recent CRAN submissions, CRAN maintainers have been 
seeing URL check issues in checks by servers hosted outside the USA. This is
likely due to US Government servers not "liking" automated checking from
servers outside the USA.

------------------------------------------------------------------------

## R CMD check results

`devtools::check()` result:

**Test environment:** aarch64-apple-darwin20 (64-bit), R 4.3.1 (2023-06-16)

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

**Test environment:** x86_64-w64-mingw32, R Under development (unstable) (2023-07-21 r84722 ucrt)

**0 errors ✔ \| 0 warnings ✔ \| 2 notes ✖**

-   checking for non-standard things in the check directory ...
    -   NOTE: Found the following files/directories:
        ''NULL''
    -   Explanation: I was unable to locate a file called "NULL" 
        in the check directory, and only R-hub checks are throwing this
        note. Related to [R-hub issue #560](https://github.com/r-hub/rhub/issues/560) 
        and can probably be ignored.

-   checking for detritus in the temp directory ...
    -   NOTE: Found the following files/directories:
        'lastMiKTeXException'
    -   Explanation: As noted in [R-hub issue
        #503](https://github.com/r-hub/rhub/issues/503), this could be
        due to a bug/crash in MiKTeX and can likely be ignored.

------------------------------------------------------------------------

`rhub::check_for_cran()` result:

**Test environment:** Ubuntu Linux 20.04.1 LTS, R-release, GCC

**0 errors ✔ \| 0 warnings ✔ \| 1 notes ✖**

-   checking HTML version of manual ...
    -   NOTE: Skipping checking HTML validation: no command 'tidy' found
    -   Explanation: The note about HTML validation only occurs on the
        Fedora dev install and does not seem critical. The HTML version
        of the manual is able to be validated on other platforms.
        
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

We checked one reverse dependency, comparing R CMD check results
across CRAN and dev versions of this package.

-   We saw 0 new problems
-   We failed to check 0 packages
