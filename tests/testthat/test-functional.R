rank_two_curves <- function(n = 40, seed = 8) {
  set.seed(seed)
  u <- seq(0, 1, length.out = 31)
  support <- data.frame(x = runif(n), y = runif(n))
  curves <- outer(rnorm(n, sd = 2), sqrt(2) * sin(pi * u)) +
    outer(rnorm(n), sqrt(2) * cos(2 * pi * u)) + 5
  list(u = u, support = support, curves = curves)
}

test_that("FPCA eigenfunctions are orthonormal under the quadrature weights", {
  d <- rank_two_curves()
  fc <- pp_functional_covariate(d$support, d$curves, d$u,
                                coords = c("x", "y"), name = "f", ncomp = 2)
  E <- fc$eigenfunctions
  gram <- t(E) %*% (E * fc$weights)

  expect_equal(unname(gram), diag(2), tolerance = 1e-10)
  expect_equal(sum(fc$weights), diff(range(d$u)))
})

test_that("rank-two curves are reconstructed exactly from two components", {
  d <- rank_two_curves()
  fc <- pp_functional_covariate(d$support, d$curves, d$u,
                                coords = c("x", "y"), name = "f", ncomp = 2)

  rebuilt <- sweep(fc$scores %*% t(fc$eigenfunctions), 2, fc$mean_function, "+")
  expect_equal(unname(rebuilt), unname(d$curves), tolerance = 1e-10)
  expect_equal(unname(fc$cumulative_pve[2]), 1, tolerance = 1e-10)
})

test_that("scores are centred, uncorrelated, and have the eigenvalues as variances", {
  d <- rank_two_curves()
  fc <- pp_functional_covariate(d$support, d$curves, d$u,
                                coords = c("x", "y"), name = "f", ncomp = 2)
  s <- functional_scores(fc)

  expect_equal(unname(colMeans(s)), c(0, 0), tolerance = 1e-10)
  expect_equal(cor(s)[1, 2], 0, tolerance = 1e-10)
  expect_equal(unname(apply(s, 2, var)), unname(fc$eigenvalues), tolerance = 1e-10)
  expect_identical(colnames(s), c("f_pc1", "f_pc2"))
})

test_that("scores are the weighted inner products with the eigenfunctions", {
  d <- rank_two_curves()
  fc <- pp_functional_covariate(d$support, d$curves, d$u,
                                coords = c("x", "y"), name = "f", ncomp = 2)
  centred <- sweep(d$curves, 2, fc$mean_function)

  expect_equal(
    unname(centred %*% (fc$eigenfunctions * fc$weights)),
    unname(fc$scores),
    tolerance = 1e-10
  )
  expect_equal(
    unname(as.matrix(project_functional_scores(fc, d$curves))),
    unname(as.matrix(functional_scores(fc))),
    tolerance = 1e-10
  )
})

test_that("functional_coefficient is the eigenfunction expansion of the score effects", {
  d <- rank_two_curves()
  pattern <- pp_pattern(d$support, coords = c("x", "y"), window = unit_square())
  fc <- pp_functional_covariate(support = pattern, curves = d$curves,
                                grid = d$u, name = "f", ncomp = 2)
  fields <- as_pp_covariates(fc)
  fit <- fit_pp(pattern, stats::reformulate(names(fields)),
                covariates = fields, nc = 5, q = 2)

  beta <- functional_coefficient(fit, fc)
  E <- fc$eigenfunctions
  b <- coef(fit$gam)[colnames(E)]
  V <- vcov(fit$gam)[colnames(E), colnames(E)]
  se <- sqrt(rowSums((E %*% V) * E))

  expect_equal(beta$estimate, unname(drop(E %*% b)), tolerance = 1e-10)
  expect_equal(beta$standard_error, unname(se), tolerance = 1e-10)
  expect_equal(beta$lower, beta$estimate - stats::qnorm(0.975) * beta$standard_error,
               tolerance = 1e-10)
})
