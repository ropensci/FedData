# Skip a test when a remote web service is unreachable or failing.
#
# Distinguishes service outages (no response at the connection level,
# throttling, or server errors) from genuine test failures (the service
# responds, but a specific resource is missing or has moved): only the
# former are skipped. This keeps transient outages of federal data
# services from failing continuous integration runs.
skip_if_service_unavailable <- function(url) {
  status <-
    tryCatch(
      httr::status_code(
        httr::GET(
          url,
          # Only request the first byte, so that resource URLs
          # can be probed without downloading them
          httr::add_headers(Range = "bytes=0-0"),
          httr::timeout(60)
        )
      ),
      error = function(e) NA_integer_
    )

  testthat::skip_if(
    is.na(status) || status == 429 || status >= 500,
    message = paste0("Service unavailable: ", url)
  )

  invisible(TRUE)
}

# Expect that a web resource is available, retrying to ride out
# flapping services. Passes as soon as any attempt succeeds. If every
# attempt fails, a persistent client error (4xx other than 429) fails
# the test -- it indicates a moved or broken URL -- while anything
# outage-shaped (connection failures, throttling, or server errors)
# skips instead.
expect_resource_available <- function(url, tries = 3, wait = 20) {
  status <- NA_integer_

  for (i in seq_len(tries)) {
    status <-
      tryCatch(
        httr::status_code(
          httr::GET(
            url,
            httr::add_headers(Range = "bytes=0-0"),
            httr::timeout(60)
          )
        ),
        error = function(e) NA_integer_
      )

    if (!is.na(status) && status < 400) {
      return(testthat::succeed())
    }

    if (i < tries) {
      Sys.sleep(wait)
    }
  }

  if (!is.na(status) && status < 500 && status != 429) {
    testthat::fail(
      paste0("Resource unavailable (HTTP ", status, "): ", url)
    )
  } else {
    testthat::skip(
      paste0("Service unavailable: ", url)
    )
  }
}
