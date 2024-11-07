---
editor_options: 
  markdown: 
    wrap: 72
---

# CRAN comments

This is a minor release from FedData v4.0.1 to FedData v4.1.0. Please
see NEWS.md for release details.

FedData was Archived on 2024-08-19 as requires archived package 'arcgislayers'.
The 'arcgislayers' package has since been restored.

My recent re-submission was rejected for several issues that were not
detailed in the checks below. All issues have been addressed. Note that checks
that use --run-donttest will take a long-ish time.

These are the issues we've addressed:

- Please omit the redundant "Functions to" from the title and description.
- Please do not start the description with the title, "This package", package name, or similar.
- Please provide a link to the used webservices to the description field of your DESCRIPTION file in the form
<http:...> or <https:...>
with angle brackets for auto-linking and no space after 'http:' and 'https:'.
- Please write TRUE and FALSE instead of T and F. Please don't use "T" or "F" as vector names. -> Warning: 'T' and 'F' instead of TRUE and FALSE:
  - man/download_ghcn_daily_station.Rd:
    - `download_ghcn_daily_station(ID, raw.dir, force.redo = F)`
  - man/get_ghcn_daily_station.Rd:
  ```
    get_ghcn_daily_station(
      ID,
      elements = NULL,
      years = NULL,
      raw.dir,
      standardize = F,
      force.redo = F
    )
    ```
  - man/get_ghcn_daily.Rd:
  ```
    get_ghcn_daily(
      template = NULL,
      label = NULL,
      elements = NULL,
      years = NULL,
      raw.dir = file.path(tempdir(), "FedData", "raw", "ghcn"),
      extraction.dir = file.path(tempdir(), "FedData", "extractions", "ned", label),
      standardize = F,
      force.redo = F
    )
    ```
  - man/sequential_duplicated.Rd:
    `sequential_duplicated(x, rows = F)`

- Please add `\value` to .Rd files regarding exported methods and explain the functions results in the documentation. Please write about the structure of the output (class) and also what the output means. (If a function does not return a value, please document that too, e.g. `\value{No return value, called for side effects}` or similar) -> Missing Rd-tags:
     - pipe.Rd: `\arguments`,  `\value`
     - replace_null.Rd: `\value`

- `\dontrun{}` should only be used if the example really cannot be executed (e.g. because of missing additional software, missing API keys, ...) by the user. That's why wrapping examples in `\dontrun{}` adds the comment ("# Not run:") as a warning for the user. Does not seem necessary. Please replace `\dontrun` with `\donttest`.

- Please unwrap the examples if they are executable in < 5 sec, or replace `\dontrun{}` with `\donttest{}`.

- Please put functions which download data in `\donttest{}`.

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

**Test environment:** Windows Server 2022 x64 (build 20348), R-devel, 64 bit

**0 errors ✔ \| 0 warnings ✔ \| 1 note ✖**

- CRAN repository db overrides:
  - X-CRAN-Comment: Archived on 2024-08-19 as requires archived package
    'arcgislayers'.
- Explanation: This is a resubmission and update; 
'arcgislayers' is back on CRAN.

------------------------------------------------------------------------

## revdepcheck results

We checked 0 reverse dependencies, comparing R CMD check results 
across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

