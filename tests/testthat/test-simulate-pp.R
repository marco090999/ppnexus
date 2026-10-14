test_that("simulate_pp is reproducible for a fixed seed", {
  a <- simulate_pp(unit_square(), expected_n = 30, seed = 9)
  b <- simulate_pp(unit_square(), expected_n = 30, seed = 9)
  c <- simulate_pp(unit_square(), expected_n = 30, seed = 10)

  expect_identical(a$coords, b$coords)
  expect_false(identical(a$coords, c$coords))
})

test_that("simulated points lie inside the window in every dimension", {
  w3 <- pp_window(lower = c(x = 0, y = 0, z = 0), upper = c(x = 1, y = 2, z = 0.5))
  x <- simulate_pp(w3, expected_n = 300, seed = 5)

  ranges <- apply(x$coords, 2, range)
  expect_true(all(ranges[1, ] >= w3$bounds["lower", ]))
  expect_true(all(ranges[2, ] <= w3$bounds["upper", ]))
})

test_that("homogeneous counts have Poisson mean and variance", {
  sims <- simulate_pp(unit_square(), expected_n = 50, nsim = 200, seed = 4)
  counts <- vapply(sims, n_points, integer(1))

  # Mean of 200 Poisson(50) counts has standard error 0.5.
  expect_lt(abs(mean(counts) - 50), 4 * 0.5)
  # Dispersion index of a Poisson sample is close to one.
  expect_gt(var(counts) / mean(counts), 0.7)
  expect_lt(var(counts) / mean(counts), 1.4)
})

test_that("inhomogeneous thinning follows the requested intensity", {
  # lambda(x, y) = 2000 * exp(2 x) on the unit square.
  x <- simulate_pp(
    unit_square(),
    log_intensity = function(coords) log(2000) + 2 * coords$x,
    bound = "user",
    lambda_upper = 2000 * exp(2),
    seed = 21
  )
  n <- n_points(x)
  expected_n <- 2000 * (exp(2) - 1) / 2
  expect_lt(abs(n - expected_n), 4 * sqrt(expected_n))

  # The x-coordinates have density proportional to exp(2 x), whose mean is
  # (e^2 + 1) / (2 (e^2 - 1)).
  mean_x <- (exp(2) + 1) / (2 * (exp(2) - 1))
  expect_lt(abs(mean(x$coords[, "x"]) - mean_x), 0.02)
})

test_that("a user bound below the intensity is detected", {
  expect_error(
    simulate_pp(
      unit_square(),
      log_intensity = function(coords) log(100) + 2 * coords$x,
      bound = "user",
      lambda_upper = 100,
      on_bound_violation = "error",
      seed = 1
    )
  )
})
