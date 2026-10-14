linear_field <- function() {
  support <- expand.grid(x = seq(0, 1, by = 0.25), y = seq(0, 1, by = 0.25))
  support$v <- 3 + 2 * support$x - support$y
  pp_covariate(
    support,
    coords = c("x", "y"),
    value = "v",
    window = unit_square(),
    name = "v"
  )
}

test_that("interpolation reproduces the field exactly at support points", {
  field <- linear_field()
  at <- data.frame(x = c(0.25, 0.5, 0.75), y = c(0.5, 0.75, 0))
  truth <- 3 + 2 * at$x - at$y

  expect_equal(evaluate_covariate(field, at, method = "idw"), truth)
  expect_equal(evaluate_covariate(field, at, method = "nearest"), truth)
})

test_that("nearest-neighbour evaluation returns the closest support value", {
  field <- linear_field()
  at <- data.frame(x = c(0.30, 0.61), y = c(0.26, 0.90))
  # Nearest support points are (0.25, 0.25) and (0.5, 1).
  expect_equal(
    evaluate_covariate(field, at, method = "nearest"),
    c(3 + 2 * 0.25 - 0.25, 3 + 2 * 0.5 - 1)
  )
})

test_that("inverse-distance weighting of a constant field is constant", {
  support <- expand.grid(x = seq(0, 1, by = 0.25), y = seq(0, 1, by = 0.25))
  support$value <- 7
  field <- pp_covariate(
    support,
    coords = c("x", "y"),
    value = "value",
    window = unit_square(),
    name = "constant"
  )
  set.seed(1)
  at <- data.frame(x = runif(10), y = runif(10))
  expect_equal(evaluate_covariate(field, at, method = "idw"), rep(7, 10))
})

test_that("geometric supports return Euclidean distances to the nearest source", {
  sources <- pp_covariate(
    data.frame(x = c(0.2, 0.8), y = c(0.2, 0.7)),
    coords = c("x", "y"),
    value = NULL,
    window = unit_square(),
    name = "distance"
  )
  at <- data.frame(x = c(0.15, 0.50, 0.85), y = c(0.20, 0.50, 0.80))

  src <- cbind(c(0.2, 0.8), c(0.2, 0.7))
  expected <- apply(as.matrix(at), 1, function(p) {
    min(sqrt((src[, 1] - p[1])^2 + (src[, 2] - p[2])^2))
  })

  expect_equal(evaluate_covariate(sources, at), unname(expected))
})

test_that("extract_covariates agrees with evaluate_covariate column by column", {
  field <- linear_field()
  sources <- pp_covariate(
    data.frame(x = c(0.2, 0.8), y = c(0.25, 0.75)),
    coords = c("x", "y"),
    window = unit_square(),
    name = "distance"
  )
  g <- pp_grid(unit_square(), n = 4)
  values <- extract_covariates(
    covariates = list(v = field, distance = sources),
    at = g,
    methods = c(v = "idw", distance = "distance_min")
  )

  expect_equal(values$v, evaluate_covariate(field, g, method = "idw"))
  expect_equal(values$distance, evaluate_covariate(sources, g))
})

test_that("unsupported interpolation methods are rejected", {
  expect_error(
    evaluate_covariate(linear_field(), data.frame(x = 0.5, y = 0.5), method = "linear"),
    "not available"
  )
})
