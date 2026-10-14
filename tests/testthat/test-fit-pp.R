# For an intercept-only Poisson model the cubature likelihood has the closed
# form maximiser exp(beta0) = n / |W|, because the cubature weights sum to the
# hypervolume of the window.

test_that("intercept-only fits recover n / |W| exactly in 2, 3 and space-time dimensions", {
  w2 <- pp_window(lower = c(x = 0, y = 0), upper = c(x = 2, y = 1))
  x2 <- simulate_pp(w2, expected_n = 200, seed = 1)
  f2 <- fit_pp(x2, ~1, nc = c(x = 6, y = 3))
  expect_equal(unname(exp(coef(f2$gam))), n_points(x2) / 2, tolerance = 1e-6)
  expect_equal(f2$cubature_summary$sum_weights, 2)

  w3 <- pp_window(lower = c(x = 0, y = 0, z = 0), upper = c(x = 1, y = 2, z = 0.5))
  x3 <- simulate_pp(w3, expected_n = 150, seed = 2)
  f3 <- fit_pp(x3, ~1, nc = 4)
  expect_equal(unname(exp(coef(f3$gam))), n_points(x3) / 1, tolerance = 1e-6)

  wst <- space_time_window(tmax = 10)
  xst <- simulate_pp(wst, expected_n = 120, seed = 3)
  fst <- fit_pp(xst, ~1, nc = c(x = 4, y = 4, time = 5))
  expect_equal(unname(exp(coef(fst$gam))), n_points(xst) / 10, tolerance = 1e-6)
})

test_that("categorical marks give the per-level counts in closed form", {
  set.seed(2)
  type <- factor(rep(c("A", "B", "C"), c(10, 30, 50)))
  x <- pp_pattern(
    data.frame(x = runif(90), y = runif(90), type = type),
    coords = c("x", "y"),
    marks = "type",
    window = unit_square()
  )
  fit <- fit_pp(x, ~type, categorical_marks = "type", nc = 5)

  expect_equal(
    unname(exp(coef(fit$gam))),
    c(10, 30 / 10, 50 / 10),
    tolerance = 1e-6
  )
  expect_equal(
    predict_pp_fit(fit, at = data.frame(x = 0.5, y = 0.5), type = "ground"),
    90,
    tolerance = 1e-6
  )
})

test_that("joint continuous marks integrate to the observed count", {
  set.seed(5)
  x <- pp_pattern(
    data.frame(x = runif(80), y = runif(80), m = runif(80, 0, 2)),
    coords = c("x", "y"),
    marks = "m",
    window = unit_square()
  )

  # Regression test: no spurious warning about `method`/`args`.
  expect_no_warning(
    fit <- fit_pp(
      x, ~m,
      continuous_marks = list(m = list(mode = "joint", bounds = c(0, 2))),
      nc = c(x = 4, y = 4, m = 4)
    )
  )
  expect_equal(fit$cubature_summary$continuous_hypervolume, 2)

  # Ground intensity = integral over the mark of exp(b0 + b1 m).
  b <- unname(coef(fit$gam))
  ground <- exp(b[1]) * (exp(2 * b[2]) - 1) / b[2]
  expect_equal(
    predict_pp_fit(fit, at = data.frame(x = 0.3, y = 0.7), type = "ground"),
    ground,
    tolerance = 1e-3
  )
  expect_equal(integrate_pp_fit(fit, nc = 4), 80, tolerance = 1e-3)
})

test_that("log-linear intensities are recovered within sampling error", {
  truth <- c(log(400), 1, -1)
  x <- simulate_pp(
    unit_square(),
    log_intensity = function(coords) truth[1] + truth[2] * coords$x + truth[3] * coords$y,
    bound = "user",
    lambda_upper = 400 * exp(1),
    seed = 11
  )
  fit <- fit_pp(x, ~ x + y, nc = 10)

  est <- unname(coef(fit$gam))
  se <- sqrt(diag(vcov(fit$gam)))
  expect_true(all(abs(est - truth) < 4 * se))
})

test_that("predictions and integrals agree with the closed-form fitted intensity", {
  x <- simulate_pp(
    unit_square(),
    log_intensity = function(coords) log(400) + coords$x - coords$y,
    bound = "user",
    lambda_upper = 400 * exp(1),
    seed = 11
  )
  fit <- fit_pp(x, ~ x + y, nc = 10)
  b <- unname(coef(fit$gam))

  at <- data.frame(x = c(0.2, 0.7), y = c(0.3, 0.9))
  expect_equal(
    predict_pp_fit(fit, at = at, type = "ground"),
    exp(b[1] + b[2] * at$x + b[3] * at$y),
    tolerance = 1e-8
  )

  analytic <- exp(b[1]) * (exp(b[2]) - 1) / b[2] * (exp(b[3]) - 1) / b[3]
  expect_equal(integrate_pp_fit(fit, nc = 40), analytic, tolerance = 1e-3)

  left_half <- rbind(lower = c(x = 0, y = 0), upper = c(x = 0.5, y = 1))
  analytic_left <- exp(b[1]) * (exp(b[2] / 2) - 1) / b[2] * (exp(b[3]) - 1) / b[3]
  expect_equal(
    integrate_pp_fit(fit, region = left_half, nc = c(x = 20, y = 40)),
    analytic_left,
    tolerance = 1e-3
  )
})

test_that("covariate effects are recovered from an external field", {
  support <- expand.grid(x = seq(0, 1, by = 0.25), y = seq(0, 1, by = 0.25))
  support$v <- 3 + 2 * support$x - support$y
  v <- pp_covariate(support, coords = c("x", "y"), value = "v",
                    window = unit_square(), name = "v")

  x <- simulate_pp(
    unit_square(),
    log_intensity = function(coords) log(100) + 0.5 * (3 + 2 * coords$x - coords$y),
    bound = "user",
    lambda_upper = 100 * exp(2.5),
    seed = 12
  )
  fit <- fit_pp(x, ~v, covariates = list(v = v), nc = 10)

  est <- unname(coef(fit$gam))[2]
  se <- unname(sqrt(diag(vcov(fit$gam))))[2]
  expect_lt(abs(est - 0.5), 4 * se)
})

test_that("non-integer cubature responses do not produce warnings", {
  # Regression test: stats::poisson()$aic used dpois() on non-integer
  # Berman--Turner responses and emitted one warning per quadrature row.
  x <- simulate_pp(odd_rectangle(), expected_n = 100, seed = 1)
  expect_no_warning(fit <- fit_pp(x, ~x, nc = 7))
  expect_true(is.finite(fit$aic))
  expect_equal(
    unname(exp(coef(fit_pp(x, ~1, nc = 7)$gam))),
    n_points(x) / odd_rectangle()$hypervolume,
    tolerance = 1e-6
  )
})

test_that("the continuous Poisson family matches stats::poisson for integer data", {
  fam <- .pp_poisson_family()
  y <- c(0, 1, 3, 7)
  mu <- c(0.5, 2, 3, 6)
  wt <- c(1, 2, 1, 0.5)

  expect_equal(fam$aic(y, 1, mu, wt, 0), stats::poisson()$aic(y, 1, mu, wt, 0))
  expect_equal(fam$linkfun(mu), log(mu))
})

test_that("refit_pp reproduces the original fit on the original pattern", {
  x <- simulate_pp(unit_square(), expected_n = 60, seed = 13)
  fit <- fit_pp(x, ~ x + y, nc = 6)
  refit <- refit_pp(fit, pattern = x)

  expect_equal(coef(refit$gam), coef(fit$gam), tolerance = 1e-8)
})

test_that("simulate_pp_fit draws counts consistent with the fitted integral", {
  x <- simulate_pp(unit_square(), expected_n = 200, seed = 14)
  fit <- fit_pp(x, ~1, nc = 6)
  sims <- simulate_pp_fit(fit, nsim = 50, bound_n = 8, seed = 15)
  counts <- vapply(sims, n_points, integer(1))

  lambda <- n_points(x)
  expect_lt(abs(mean(counts) - lambda), 4 * sqrt(lambda / 50))
})
