test_that("pp_window stores bounds, dimension and hypervolume", {
  w2 <- pp_window(lower = c(x = 0, y = 0), upper = c(x = 2, y = 0.5))
  expect_s3_class(w2, "pp_window")
  expect_identical(w2$dimension, 2L)
  expect_equal(w2$hypervolume, 1)
  expect_equal(unname(w2$bounds["upper", ]), c(2, 0.5))

  w3 <- pp_window(lower = c(x = 0, y = 0, z = -1), upper = c(x = 1, y = 2, z = 2))
  expect_identical(w3$dimension, 3L)
  expect_equal(w3$hypervolume, 6)

  wst <- space_time_window(tmax = 365)
  expect_equal(wst$hypervolume, 365)
  expect_identical(wst$axis_info$role, c("space", "space", "time"))
})

test_that("pp_window rejects degenerate or inconsistent bounds", {
  expect_error(pp_window(lower = c(x = 0, y = 0), upper = c(x = 0, y = 1)))
  expect_error(pp_window(lower = c(x = 1, y = 0), upper = c(x = 0, y = 1)))
  expect_error(pp_window(lower = c(x = 0, y = 0), upper = c(x = 1)))
})

test_that("pp_pattern keeps coordinates and marks in input order", {
  dat <- data.frame(
    x = c(0.1, 0.5, 0.9),
    y = c(0.2, 0.4, 0.6),
    magnitude = c(1.5, 2.5, 3.5),
    type = c("A", "B", "A")
  )
  pat <- pp_pattern(
    dat,
    coords = c("x", "y"),
    marks = c("magnitude", "type"),
    window = unit_square()
  )

  expect_s3_class(pat, "pp")
  expect_equal(unname(pat$coords[, "x"]), dat$x)
  expect_equal(unname(pat$coords[, "y"]), dat$y)
  expect_equal(pat$marks$magnitude, dat$magnitude)
  expect_equal(as.character(pat$marks$type), dat$type)
})

test_that("pp_pattern rejects points outside the observation window", {
  expect_error(
    pp_pattern(
      data.frame(x = c(0.5, 1.5), y = c(0.5, 0.5)),
      coords = c("x", "y"),
      window = unit_square()
    )
  )
})

test_that("pp_grid places cell centres on a regular lattice", {
  w <- pp_window(lower = c(x = 0, y = 0), upper = c(x = 1, y = 2))
  g <- pp_grid(w, n = c(x = 5, y = 4))

  expect_s3_class(g, "pp_grid")
  expect_identical(nrow(g$coords), 20L)
  expect_equal(sort(unique(g$coords[, "x"])), seq(0, 1, length.out = 5))
  expect_equal(sort(unique(g$coords[, "y"])), seq(0, 2, length.out = 4))
})

test_that("pp_window has a compact print method", {
  w <- pp_window(
    lower = c(x = 0, y = 0, time = 0),
    upper = c(x = 10, y = 5, time = 30),
    axis_roles = c(x = "space", y = "space", time = "time"),
    axis_units = c(x = "km", y = "km", time = "days")
  )
  out <- capture.output(print(w))
  expect_identical(out[1], "<pp_window>")
  expect_true(any(grepl("time \\[days\\] \\(time\\): \\[0, 30\\]", out)))
  expect_true(any(grepl("Hypervolume: 1500", out)))
  expect_invisible(print(w))
})
