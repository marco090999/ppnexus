# Shared fixtures for the ppnexus test suite.

unit_square <- function() {
  pp_window(lower = c(x = 0, y = 0), upper = c(x = 1, y = 1))
}

# A rectangle whose cubature weights have non-integer reciprocals. This
# exercises the non-integer Berman--Turner responses of fit_pp().
odd_rectangle <- function() {
  pp_window(lower = c(x = 0, y = 0), upper = c(x = 1.37, y = 0.83))
}

space_time_window <- function(tmax = 10) {
  pp_window(
    lower = c(x = 0, y = 0, time = 0),
    upper = c(x = 1, y = 1, time = tmax),
    axis_roles = c(x = "space", y = "space", time = "time")
  )
}

n_points <- function(pattern) {
  nrow(pattern$coords)
}

# Reference implementation, in plain R, of the triggered ETAS intensity
# computed by the Fortran routine src/etas_trigger.f90.
etas_trigger_reference <- function(tflag, k, c, p, g, d, q,
                                   x, y, t, m, predictor) {
  n <- length(t)
  out <- numeric(n)

  for (i in seq_len(n)[-1L]) {
    j <- seq_len(i - 1L)
    dt <- t[i] - t[j]
    keep <- dt > 0
    j <- j[keep]
    dt <- dt[keep]

    temporal <- (dt + c)^(-p) * exp(predictor[j])

    if (tflag > 0L) {
      contribution <- temporal
    } else {
      ds <- (x[i] - x[j])^2 + (y[i] - y[j])^2
      contribution <- temporal * (ds / exp(g * m[j]) + d)^(-q)
    }

    out[i] <- k * sum(contribution)
  }

  out
}

# Number of pages drawn on a PDF device while evaluating `expr`.
count_pdf_pages <- function(expr) {
  file <- tempfile(fileext = ".pdf")
  on.exit(unlink(file), add = TRUE)
  grDevices::pdf(file, compress = FALSE)
  device_open <- TRUE
  on.exit(if (device_open) grDevices::dev.off(), add = TRUE, after = FALSE)
  force(expr)
  grDevices::dev.off()
  device_open <- FALSE

  content <- readLines(file, warn = FALSE, encoding = "bytes")
  pages <- grepl("/Type /Page\\b", content, useBytes = TRUE) &
    !grepl("/Type /Pages", content, useBytes = TRUE)
  as.integer(sum(pages))
}
