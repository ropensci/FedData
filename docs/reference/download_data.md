# Use curl to download a file.

This function makes it easy to implement timestamping and no-clobber of
files.

## Usage

``` r
download_data(
  url,
  destdir = getwd(),
  timestamping = TRUE,
  nc = FALSE,
  verbose = FALSE,
  progress = FALSE
)
```

## Arguments

- url:

  The location of a file.

- destdir:

  Where the file should be downloaded to.

- timestamping:

  Should only newer files be downloaded?

- nc:

  Should files of the same type not be clobbered?

- verbose:

  Should cURL output be shown?

- progress:

  Should a progress bar be shown with cURL output?

## Value

A character string of the file path to the downloaded file.

## Details

If both `timestamping` and `nc` are TRUE, nc behavior trumps
timestamping.
