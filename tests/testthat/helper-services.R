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
        httr::GET(url, httr::timeout(60))
      ),
      error = function(e) NA_integer_
    )

  testthat::skip_if(
    is.na(status) || status == 429 || status >= 500,
    message = paste0("Service unavailable: ", url)
  )

  invisible(TRUE)
}
