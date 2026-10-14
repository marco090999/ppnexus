fitted_loglinear <- function() {
  x <- simulate_pp(
    unit_square(),
    log_intensity = function(coords) log(400) + coords$x - coords$y,
    bound = "user",
    lambda_upper = 400 * exp(1),
    seed = 11
  )
  list(pattern = x, fit = fit_pp(x, ~ x + y, nc = 10))
}

test_that("raw residuals are observed minus expected counts per cell", {
  obj <- fitted_loglinear()
  res <- residuals_pp(obj$fit, type = "raw", nc = 5)
  d <- res$data

  expect_equal(d$.pp_raw, d$.pp_observed - d$.pp_expected)
  expect_equal(sum(d$.pp_observed), n_points(obj$pattern))
  expect_identical(nrow(d), 25L)
})

test_that("Pearson and deviance residuals follow their textbook definitions", {
  obj <- fitted_loglinear()
  d <- residuals_pp(obj$fit, type = "pearson", nc = 5)$data

  o <- d$.pp_observed
  e <- d$.pp_expected
  expect_equal(d$.pp_pearson, (o - e) / sqrt(e))

  dev <- sign(o - e) * sqrt(2 * (ifelse(o > 0, o * log(o / e), 0) - (o - e)))
  expect_equal(d$.pp_deviance, dev)
})

test_that("cell expectations add up to the integrated fitted intensity", {
  obj <- fitted_loglinear()
  b <- unname(coef(obj$fit$gam))
  analytic <- exp(b[1]) * (exp(b[2]) - 1) / b[2] * (exp(b[3]) - 1) / b[3]

  res <- residuals_pp(obj$fit, type = "raw", nc = 20)
  expect_equal(sum(res$data$.pp_expected), analytic, tolerance = 2e-3)
})

test_that("residual-profile envelopes are drawn once per variable", {
  x <- simulate_pp(unit_square(), expected_n = 80, seed = 27)
  fit <- fit_pp(x, ~ x + y, nc = 6)
  diagnostic <- diagnose_first_order_pp(fit, nsim = 19, alpha = 0.1, nc = 5,
                                        seed = 28, verbose = FALSE)

  expect_identical(count_pdf_pages(plot(diagnostic)), 2L)
  expect_identical(count_pdf_pages(plot(diagnostic, variable = "y")), 1L)
})
