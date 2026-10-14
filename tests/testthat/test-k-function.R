test_that("the translation-corrected estimator matches a hand calculation", {
  pts <- data.frame(x = c(0.1, 0.2, 0.5), y = c(0.1, 0.1, 0.5))
  pattern <- pp_pattern(pts, coords = c("x", "y"), window = unit_square())
  r <- c(0.05, 0.15, 1)
  lambda <- rep(3, 3)

  k <- K_pp(pattern, lambda = lambda, r = r)

  P <- as.matrix(pts)
  hand <- vapply(r, function(rr) {
    total <- 0
    for (i in 1:3) for (j in 1:3) {
      if (i != j && sqrt(sum((P[i, ] - P[j, ])^2)) <= rr) {
        overlap <- prod(1 - abs(P[i, ] - P[j, ]))
        total <- total + 1 / (lambda[i] * lambda[j] * overlap)
      }
    }
    total
  }, numeric(1))

  expect_equal(unname(k$Khat), hand)
  expect_equal(k$correction, "translate")
})

test_that("K scales with the inverse square of the supplied intensity", {
  x <- simulate_pp(unit_square(), expected_n = 300, seed = 6)
  n <- n_points(x)
  r <- c(0.05, 0.1)

  k_default <- K_pp(x, r = r)
  k_const <- K_pp(x, lambda = rep(n, n), r = r)
  k_double <- K_pp(x, lambda = rep(2 * n, n), r = r)

  expect_equal(k_default$Khat, k_const$Khat)
  expect_equal(k_double$Khat / k_const$Khat, c(0.25, 0.25), ignore_attr = TRUE)
})

test_that("theoretical Poisson K is the volume of a ball (times 2u in space-time)", {
  r <- c(0, 0.05, 0.1)

  x2 <- simulate_pp(unit_square(), expected_n = 50, seed = 1)
  expect_equal(unname(K_pp(x2, r = r)$Ktheo), pi * r^2)

  w3 <- pp_window(lower = c(x = 0, y = 0, z = 0), upper = c(x = 1, y = 1, z = 1))
  x3 <- simulate_pp(w3, expected_n = 50, seed = 2)
  expect_equal(unname(K_pp(x3, r = r)$Ktheo), 4 / 3 * pi * r^3)

  xst <- simulate_pp(space_time_window(), expected_n = 50, seed = 3)
  kst <- K_pp(xst, r = c(0.05, 0.1), u = c(0.5, 1))
  expect_equal(
    unname(kst$Ktheo),
    outer(pi * c(0.05, 0.1)^2, 2 * c(0.5, 1)),
    ignore_attr = TRUE
  )
})

test_that("the estimator is close to the Poisson value for CSR", {
  x <- simulate_pp(unit_square(), expected_n = 2000, seed = 6)
  k <- K_pp(x, r = seq(0.02, 0.1, by = 0.02))
  expect_true(all(abs(k$Khat / k$Ktheo - 1) < 0.05))

  w3 <- pp_window(lower = c(x = 0, y = 0, z = 0), upper = c(x = 1, y = 1, z = 1))
  x3 <- simulate_pp(w3, expected_n = 3000, seed = 7)
  k3 <- K_pp(x3, r = c(0.05, 0.1))
  expect_true(all(abs(k3$Khat / k3$Ktheo - 1) < 0.05))
})

test_that("local contributions average to the global estimate", {
  x <- simulate_pp(unit_square(), expected_n = 300, seed = 6)
  r <- c(0.05, 0.1)
  k <- K_pp(x, r = r, local = TRUE)

  expect_identical(dim(k$Klocal), c(length(r), n_points(x)))
  expect_equal(unname(rowMeans(k$Klocal)), unname(k$Khat))

  kl <- Klocal_pp(x, r = r)
  expect_equal(kl$Klocal, k$Klocal)
})

test_that("K_fit_pp with an intercept-only model equals K_pp with constant intensity", {
  x <- simulate_pp(unit_square(), expected_n = 150, seed = 16)
  fit <- fit_pp(x, ~1, nc = 5)
  r <- c(0.05, 0.1, 0.15)

  expect_equal(
    unname(K_fit_pp(fit, r = r)$Khat),
    unname(K_pp(x, r = r)$Khat),
    tolerance = 1e-6
  )
})

test_that("plotting a spatial K envelope draws the GET figure", {
  # Regression test: plot.pp_envelope() did not print the ggplot returned by
  # GET, so nothing was drawn.
  x <- simulate_pp(unit_square(), expected_n = 60, seed = 17)
  fit <- fit_pp(x, ~1, nc = 5)
  env <- envelope_pp(fit, nsim = 19, alpha = 0.1, r = seq(0, 0.2, length.out = 6),
                     simulation_args = list(bound_n = 6), seed = 18,
                     verbose = FALSE)

  drawn <- FALSE
  local_mocked_bindings(
    print.ggplot = function(x, ...) {
      drawn <<- TRUE
      invisible(x)
    },
    .package = "ggplot2"
  )
  pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot(env)
  expect_true(drawn)
})
