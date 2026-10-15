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
  # Regression test: the non-actionable grid-bound warning of simulate_pp()
  # is no longer emitted for LGCPs.
  expect_no_warning(a <- simulate_reference_lgcp("diagnostic"))
  b <- simulate_reference_lgcp("diagnostic")
  expect_identical(a$pattern$coords, b$pattern$coords)

  field <- a$lgcp$field
  # exp(Z - sigma2 / 2) has unit mean, so the expected count is preserved.
  expect_equal(field$centred_value, field$value - 0.5 / 2)
  expect_equal(a$lgcp$marginal_variance, 0.5)
})

test_that("the fitted LGCP K-function solves the exponential-model integral", {
  skip_on_cran()
  x <- simulate_reference_lgcp()
  r <- seq(0.01, 0.25, length.out = 15)
  fit <- fit_lgcp_pp(pattern = x, formula = ~1, cov_model = "exponential",
                     r = r, nc = 8,
                     control = lgcp_control(integration_n = 150))
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
  x <- simulate_reference_lgcp()
  fit <- fit_lgcp_pp(pattern = x, formula = ~1, cov_model = "exponential",
                     r = seq(0.01, 0.25, length.out = 15), nc = 8,
                     control = lgcp_control(integration_n = 150))

  # Loose bounds: a single realisation of a 400-point LGCP.
  expect_gt(fit$theta[["sigma2"]], 0.2)
  expect_lt(fit$theta[["sigma2"]], 1.0)
  expect_gt(fit$theta[["alpha"]], 0.04)
  expect_lt(fit$theta[["alpha"]], 0.25)
})

test_that("fit_lgcp_pp forwards first-order settings and validates control", {
  x <- simulate_reference_lgcp()
  r <- seq(0.01, 0.25, length.out = 10)

  fit <- fit_lgcp_pp(x, ~1, cov_model = "exponential", r = r, nc = 6,
                     control = list(integration_n = 100))
  expect_s3_class(fit, "pp_lgcp_fit")
  expect_identical(fit$settings$integration_n, 100)

  expect_error(
    fit_lgcp_pp(x, ~1, cov_model = "exponential", r = r, bogus = 1),
    "Unknown argument"
  )
  expect_error(
    fit_lgcp_pp(x, ~1, cov_model = "exponential", r = r,
                control = list(bogus = 1)),
    "Unknown `control` setting"
  )

  # Previous interface: still accepted, with a warning.
  expect_warning(
    old <- fit_lgcp_pp(x, ~1, cov_model = "exponential", r = r, nc = 6,
                       integration_n = 100),
    "previous interface"
  )
  expect_equal(old$theta, fit$theta)
})

test_that("a custom covariance function gives the built-in estimates", {
  # Regression test: with unbounded parameters L-BFGS-B stopped at the
  # starting values; Nelder-Mead is now used when no bound is finite.
  skip_on_cran()
  x <- simulate_reference_lgcp()
  r <- seq(0.01, 0.25, length.out = 12)
  builtin <- fit_lgcp_pp(x, ~1, cov_model = "exponential", r = r, nc = 8,
                         control = lgcp_control(integration_n = 150))
  custom <- fit_lgcp_pp(
    x, ~1,
    cov_model = function(h, params) params[["sigma2"]] * exp(-h / params[["alpha"]]),
    r = r, nc = 8,
    control = lgcp_control(start = c(sigma2 = 1, alpha = 0.2), integration_n = 150)
  )

  expect_false(isTRUE(all.equal(unname(custom$theta), c(1, 0.2))))
  expect_equal(custom$theta, builtin$theta, tolerance = 0.02)
})

test_that("a spatio-temporal LGCP is fitted with the default bounds", {
  # Regression test: the default start and bounds of `alpha_time` carried the
  # name of the time axis, so they became NA and the fit stopped with an error.
  skip_on_cran()
  w <- space_time_window(tmax = 5)
  x <- simulate_lgcp_pp(
    window = w,
    log_intensity = function(coords) rep(log(60), nrow(coords)),
    cov_model = "separable_exponential",
    covariance_params = c(sigma2 = 0.5, alpha_space = 0.1, alpha_time = 0.5),
    field_n = c(x = 10, y = 10, time = 10),
    output = "pattern",
    seed = 4
  )
  fit <- fit_lgcp_pp(x, ~1, cov_model = "separable_exponential",
                     r = seq(0.02, 0.25, length.out = 6),
                     u = seq(0.2, 1.5, length.out = 5), nc = c(x = 5, y = 5, time = 5),
                     control = lgcp_control(integration_n = 60))
  expect_identical(names(fit$theta), c("sigma2", "alpha_space", "alpha_time"))
  expect_true(all(is.finite(fit$theta)))
})

test_that("multi-start minimum contrast never worsens the criterion", {
  skip_on_cran()
  x <- simulate_reference_lgcp()
  r <- seq(0.01, 0.25, length.out = 12)
  single <- fit_lgcp_pp(x, ~1, cov_model = "exponential", r = r, nc = 8,
                        control = lgcp_control(integration_n = 150, multistart = FALSE))
  multi <- fit_lgcp_pp(x, ~1, cov_model = "exponential", r = r, nc = 8,
                       control = lgcp_control(integration_n = 150))
  expect_lte(multi$second_order$contrast_value,
             single$second_order$contrast_value * (1 + 1e-8))
  expect_error(lgcp_control(multistart = NA), "multistart")
})
