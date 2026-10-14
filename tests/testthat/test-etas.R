etas_params <- function(k0 = 0.0005) {
  c(mu = 0.5, k0 = k0, c = 0.05, p = 1.2, gamma = 0, d = 0.05, q = 2)
}

simulate_reference_etas <- function(k0 = 0.0005, max_events = 5000, t.lag = 500) {
  simulate_etas_pp(
    params = etas_params(k0),
    m0 = 2.5, b = 1, tmin = 0, t.lag = t.lag,
    long.range = c(0, 1), lat.range = c(0, 1), longlat.to.km = FALSE,
    n_support_obs = 50, mult_support = 1,
    trig_lp_coefs = c(m_rel = 0.3),
    seed = 19, max_events = max_events
  )
}

test_that("the compiled triggering routine matches a plain R implementation", {
  set.seed(1)
  n <- 150L
  x <- runif(n)
  y <- runif(n)
  t <- sort(runif(n, 0, 100))
  t[c(20, 21)] <- t[20]          # ties must not trigger each other
  m <- rexp(n)
  predictor <- rnorm(n, 0, 0.3)
  par <- list(k = 0.02, c = 0.01, p = 1.1, g = 0.5, d = 1, q = 1.5)

  for (tflag in c(0L, 1L)) {
    native <- .Fortran(
      C_pp_etas_trigger,
      tflag = tflag, n = n, mu = 1,
      k = par$k, c = par$c, p = par$p, g = par$g, d = par$d, q = par$q,
      x = x, y = y, t = t, m = m, predictor = predictor, l = double(n)
    )$l
    reference <- etas_trigger_reference(
      tflag, par$k, par$c, par$p, par$g, par$d, par$q, x, y, t, m, predictor
    )

    expect_equal(native, reference, tolerance = 1e-12)
    expect_identical(native[1], 0)
  }
})

test_that("ETAS simulation is reproducible and internally consistent", {
  a <- simulate_reference_etas()
  b <- simulate_reference_etas()
  expect_identical(a$cat.sim, b$cat.sim)

  cat <- a$cat.sim
  expect_false(a$exploded)
  expect_identical(a$n_events, nrow(cat))
  expect_identical(a$n_background + a$n_triggered, a$n_events)
  expect_identical(cat$event_id, seq_len(nrow(cat)))
  expect_true(all(diff(cat$time) >= 0))
  expect_true(all(cat$magn1 >= 2.5))
  expect_true(all(cat$long >= 0 & cat$long <= 1 & cat$lat >= 0 & cat$lat <= 1))

  # Every triggered event occurs after its parent.
  children <- which(cat$father_id > 0L)
  expect_true(all(cat$time[children] > cat$time[cat$father_id[children]]))
})

test_that("background counts and magnitudes follow their generating laws", {
  sim <- simulate_reference_etas()

  # Background events: Poisson with mean mu * T = 0.5 * 500.
  expect_lt(abs(sim$n_background - 250), 4 * sqrt(250))

  # Gutenberg-Richter: m - m0 ~ Exp(b log 10), mean 1 / log(10).
  excess <- sim$cat.sim$magn1 - 2.5
  expect_lt(
    abs(mean(excess) - 1 / log(10)),
    4 / (log(10) * sqrt(length(excess)))
  )
})

test_that("a supercritical simulation stops cleanly at max_events", {
  # Regression test: the explosion guard previously left new offspring
  # without identifiers and failed internal validation.
  expect_warning(
    sim <- simulate_reference_etas(k0 = 0.003, max_events = 2000, t.lag = 200),
    "Maximum number of events exceeded"
  )
  expect_true(sim$exploded)
  expect_identical(sim$cat.sim$event_id, seq_len(nrow(sim$cat.sim)))
  expect_gt(nrow(sim$cat.sim), 2000)
})

test_that("etas_catalog converts degrees to kilometres with the etasFLP factor", {
  set.seed(42)
  n <- 20
  w <- pp_window(
    lower = c(long = 10, lat = 40, time = 0),
    upper = c(long = 12, lat = 42, time = 100),
    axis_roles = c(long = "space", lat = "space", time = "time")
  )
  eq <- pp_pattern(
    data.frame(long = runif(n, 10, 12), lat = runif(n, 40, 42),
               time = sort(runif(n, 0, 100)), z = runif(n, 1, 20),
               magn1 = 2.5 + rexp(n, 2)),
    coords = c("long", "lat", "time"), marks = c("z", "magn1"), window = w
  )

  km <- etas_catalog(eq, longlat.to.km = TRUE)
  deg <- etas_catalog(eq, longlat.to.km = FALSE)

  expect_equal(km$spatial_scale, 6371.3 * pi / 180)
  expect_equal(deg$spatial_scale, 1)
  expect_identical(nrow(km$data), as.integer(n))
})

test_that("fitted ETAS probabilities and intensities are coherent", {
  skip_on_cran()
  sim <- simulate_reference_etas()
  fit <- suppressWarnings(fit_etas_pp(
    sim$catalog, magn.threshold = 2.5, formula.bg = ~1, ndeclust = 2,
    iterlim = 30, verbose = FALSE, trace = FALSE, verbose.convergence = FALSE
  ))

  expect_equal(fit$background_probability + fit$triggered_probability,
               rep(1, sim$n_events))
  expect_equal(fit$triggered_probability,
               fit$triggering_intensity / fit$total_intensity)
  expect_true(all(fit$total_intensity > 0))
  expect_equal(fit$loglik, -fit$negative_loglik)

  declustered <- decluster_etas_pp(fit)
  expect_identical(nrow(declustered), sim$n_events)
})

test_that("ETAS fitting recovers the simulation parameters", {
  skip_on_cran()
  truth <- c(mu = 0.5, k0 = 1.5e-5, c = 0.05, p = 1.2, gamma = 0, d = 0.001, q = 2)
  sim <- simulate_etas_pp(
    params = truth, m0 = 2.5, b = 1, tmin = 0, t.lag = 500,
    long.range = c(0, 1), lat.range = c(0, 1), longlat.to.km = FALSE,
    n_support_obs = 50, mult_support = 1, trig_lp_coefs = c(m_rel = 0.3),
    seed = 19
  )
  fit <- suppressWarnings(fit_etas_pp(
    sim$catalog, magn.threshold = 2.5, formula.bg = ~1, ndeclust = 5,
    iterlim = 100, verbose = FALSE, trace = FALSE, verbose.convergence = FALSE
  ))

  expect_true(fit$convergence$converged)
  shown <- c("mu", "c", "p", "d", "q")
  z <- (fit$parameters[shown] - truth[shown]) / fit$standard_errors[shown]
  expect_true(all(abs(z) < 4))
})
