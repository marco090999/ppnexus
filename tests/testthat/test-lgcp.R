simulate_reference_lgcp <- function(output = "pattern") {
  simulate_lgcp_pp(
    window = unit_square(),
    expected_n = 400,
    cov_model = "exponential",
    covariance_params = c(sigma2 = 0.5, alpha = 0.1),
    field_n = 25,
    output = output,
    seed = 2
  )
}

test_that("LGCP simulation is reproducible and centres the latent field", {
  a <- suppressWarnings(simulate_reference_lgcp("diagnostic"))
  b <- suppressWarnings(simulate_reference_lgcp("diagnostic"))
  expect_identical(a$pattern$coords, b$pattern$coords)

  field <- a$lgcp$field
  # exp(Z - sigma2 / 2) has unit mean, so the expected count is preserved.
  expect_equal(field$centred_value, field$value - 0.5 / 2)
  expect_equal(a$lgcp$marginal_variance, 0.5)
})

test_that("the fitted LGCP K-function solves the exponential-model integral", {
  skip_on_cran()
  x <- suppressWarnings(simulate_reference_lgcp())
  r <- seq(0.01, 0.25, length.out = 15)
  fit <- fit_lgcp_pp(pattern = x, formula = ~1, cov_model = "exponential",
                     r = r, nc = 8, integration_n = 150, seed = 3)
  th <- fit$theta

  # K(r) = 2 pi int_0^r s exp(C(s)) ds with C(s) = sigma2 exp(-s / alpha).
  reference <- vapply(r, function(rr) {
    2 * pi * stats::integrate(
      function(s) s * exp(th[["sigma2"]] * exp(-s / th[["alpha"]])), 0, rr
    )$value
  }, numeric(1))

  expect_equal(fit$second_order$theoretical_K, reference, tolerance = 1e-2)
  expect_identical(names(th), c("sigma2", "alpha"))
  expect_identical(fit$second_order$optimisation$convergence, 0L)
})

test_that("minimum contrast recovers the covariance parameters", {
  skip_on_cran()
  x <- suppressWarnings(simulate_reference_lgcp())
  fit <- fit_lgcp_pp(pattern = x, formula = ~1, cov_model = "exponential",
                     r = seq(0.01, 0.25, length.out = 15), nc = 8,
                     integration_n = 150, seed = 3)

  # Loose bounds: a single realisation of a 400-point LGCP.
  expect_gt(fit$theta[["sigma2"]], 0.2)
  expect_lt(fit$theta[["sigma2"]], 1.0)
  expect_gt(fit$theta[["alpha"]], 0.04)
  expect_lt(fit$theta[["alpha"]], 0.25)
})
