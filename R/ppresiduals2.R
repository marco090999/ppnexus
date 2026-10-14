# VERSION MARKER --------------------------------------------------------------
# This is ppresiduals.R version 0.2.0 (visual diagnostics update).
# Expected checks after source():
#   .pp_res_default_nc(fit)  # coarsens the fitting grid
#   names(formals(plot.pp_residuals))
#     # includes zlim, robust, show_colorbar, show_cell_boundaries, show_points
.pp_residuals_version <- "0.2.0"

# First-order diagnostics for Poisson point-process models -------------------
# Version: 2026-07-03 (coarser default grids and improved residual graphics)
#
# This module builds on ppdata.R, ppfit.R, and ppdiagnostic.R. It provides
# regular diagnostic grids, integrated fitted ground intensity, cell residuals,
# residual profiles with respect to coordinates and fitted covariates, and
# plug-in refit global-envelope diagnostics for first-order model checking.
#
# The implementation assumes the hyperrectangular pp_window objects used by
# the rest of the package. Residuals are calculated for the ground process:
# for joint marked models, the fitted joint intensity is marginalised over the
# mark domain by predict_pp_fit(type = "ground").


# Dependencies and shared validation -----------------------------------------

.pp_res_require <- function(functions, context) {
  lookup_env <- environment()

  missing <- functions[
    !vapply(
      functions,
      exists,
      logical(1),
      envir = lookup_env,
      mode = "function",
      inherits = TRUE
    )
  ]

  if (length(missing) > 0L) {
    stop(
      sprintf(
        paste0(
          "%s requires function(s) not currently available: %s. ",
          "Source ppdata.R, ppfit.R, and ppdiagnostic.R before ppresiduals.R."
        ),
        context,
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.pp_res_validate_fit <- function(fit) {
  .pp_res_require(
    c(
      ".pp_diag_validate_fit", ".pp_diag_model_spec",
      ".pp_diag_validate_positive_integer", ".pp_diag_validate_seed",
      ".pp_diag_validate_probability", ".pp_diag_evaluate_sources",
      "predict_pp_fit"
    ),
    context = "First-order residual diagnostics"
  )

  .pp_diag_validate_fit(fit, require_model_spec = TRUE)
  invisible(fit)
}


.pp_res_validate_logical <- function(x, arg_name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", arg_name), call. = FALSE)
  }

  isTRUE(x)
}


.pp_res_axis_names <- function(fit) {
  spec <- .pp_diag_model_spec(fit)
  as.character(spec$axis_info$name)
}


.pp_res_validate_max_cells <- function(max_cells) {
  .pp_diag_validate_positive_integer(max_cells, "max_cells", min_value = 1L)
}


.pp_res_validate_region <- function(fit, region = NULL) {
  axis_names <- .pp_res_axis_names(fit)
  window_bounds <- fit$model_spec$window$bounds[, axis_names, drop = FALSE]
  rownames(window_bounds) <- c("lower", "upper")

  if (is.null(region)) {
    return(window_bounds)
  }

  if (!is.matrix(region) && !is.data.frame(region)) {
    stop(
      "`region` must be NULL or a two-row numeric matrix/data frame of bounds.",
      call. = FALSE
    )
  }

  region <- as.matrix(region)

  if (!is.numeric(region) || nrow(region) != 2L ||
      is.null(colnames(region)) ||
      !setequal(colnames(region), axis_names)) {
    stop(
      paste0(
        "`region` must be a two-row numeric matrix with columns named exactly ",
        "as the fitted coordinate axes."
      ),
      call. = FALSE
    )
  }

  region <- region[, axis_names, drop = FALSE]
  rownames(region) <- c("lower", "upper")

  if (any(!is.finite(region))) {
    stop("`region` must contain only finite bounds.", call. = FALSE)
  }

  if (any(region["lower", ] >= region["upper", ])) {
    stop("Every lower region bound must be strictly smaller than its upper bound.", call. = FALSE)
  }

  tolerance <- .Machine$double.eps^0.5 * pmax(
    1,
    abs(window_bounds["upper", ] - window_bounds["lower", ])
  )

  if (any(region["lower", ] < window_bounds["lower", ] - tolerance) ||
      any(region["upper", ] > window_bounds["upper", ] + tolerance)) {
    stop("`region` must lie entirely inside the fitted observation window.", call. = FALSE)
  }

  region
}


.pp_res_default_nc <- function(fit) {
  axis_names <- .pp_res_axis_names(fit)
  spec <- .pp_diag_model_spec(fit)

  # Diagnostic residual maps should be materially coarser than the fitting
  # cubature. Otherwise, most cells contain too few expected events and the
  # displayed Pearson/deviance residuals are dominated by Poisson noise.
  target_expected_per_cell <- 5
  max_cells_per_axis <- 10L

  candidate <- spec$quadrature$nc_resolved

  if (is.null(candidate)) {
    candidate <- spec$quadrature$nc_requested
  }

  if (is.null(candidate) && !is.null(fit$settings$nc)) {
    candidate <- fit$settings$nc
  }

  if (is.null(candidate)) {
    candidate <- rep.int(max_cells_per_axis, length(axis_names))
    names(candidate) <- axis_names
  }

  if (!is.numeric(candidate) || any(!is.finite(candidate)) ||
      any(candidate != as.integer(candidate)) || any(candidate < 1L)) {
    stop(
      "The retained cubature resolution in `fit` is invalid.",
      call. = FALSE
    )
  }

  if (!is.null(names(candidate))) {
    missing_axes <- setdiff(axis_names, names(candidate))

    if (length(missing_axes) > 0L) {
      stop(
        paste0(
          "The retained cubature resolution does not contain all coordinate axes: ",
          paste(missing_axes, collapse = ", "), "."
        ),
        call. = FALSE
      )
    }

    candidate <- as.integer(candidate[axis_names])
  } else if (length(candidate) == 1L) {
    candidate <- rep.int(as.integer(candidate), length(axis_names))
  } else if (length(candidate) == length(axis_names)) {
    candidate <- as.integer(candidate)
  } else {
    stop(
      paste0(
        "The retained cubature resolution cannot be aligned with the fitted ",
        "coordinate axes. Supply `nc` explicitly."
      ),
      call. = FALSE
    )
  }

  n_events <- nrow(fit$pattern$coords)
  d <- length(axis_names)
  target_cells <- max(1L, floor(n_events / target_expected_per_cell))
  per_axis <- max(1L, floor(target_cells^(1 / d)))

  out <- pmin(
    candidate,
    rep.int(per_axis, d),
    rep.int(max_cells_per_axis, d)
  )
  out <- pmax(1L, as.integer(out))
  names(out) <- axis_names
  out
}

.pp_res_normalize_nc <- function(fit, nc = NULL, max_cells = 1e6) {
  axis_names <- .pp_res_axis_names(fit)
  max_cells <- .pp_res_validate_max_cells(max_cells)

  if (is.null(nc)) {
    nc <- .pp_res_default_nc(fit)
  }

  if (!is.numeric(nc) || any(!is.finite(nc)) ||
      any(nc != as.integer(nc)) || any(nc < 1L)) {
    stop("`nc` must contain positive finite integer cell counts.", call. = FALSE)
  }

  if (!is.null(names(nc))) {
    if (anyDuplicated(names(nc))) {
      stop("`nc` cannot have duplicated axis names.", call. = FALSE)
    }

    missing_axes <- setdiff(axis_names, names(nc))
    extra_axes <- setdiff(names(nc), axis_names)

    if (length(missing_axes) > 0L || length(extra_axes) > 0L) {
      stop(
        sprintf(
          "Named `nc` must match the fitted coordinate axes exactly: %s.",
          paste(axis_names, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    out <- as.integer(nc[axis_names])
  } else if (length(nc) == 1L) {
    out <- rep.int(as.integer(nc), length(axis_names))
  } else if (length(nc) == length(axis_names)) {
    out <- as.integer(nc)
  } else {
    stop(
      sprintf(
        "`nc` must be scalar, named by axes, or have length %d.",
        length(axis_names)
      ),
      call. = FALSE
    )
  }

  names(out) <- axis_names
  n_cells <- prod(as.double(out))

  if (!is.finite(n_cells) || n_cells > max_cells) {
    stop(
      sprintf(
        "`nc` defines %.0f cells, exceeding `max_cells = %.0f`.",
        n_cells,
        as.double(max_cells)
      ),
      call. = FALSE
    )
  }

  out
}


.pp_res_index_grid <- function(nc) {
  parts <- lapply(nc, seq_len)
  names(parts) <- names(nc)

  out <- do.call(
    expand.grid,
    c(parts, list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  )

  out <- as.matrix(out)
  storage.mode(out) <- "integer"
  colnames(out) <- names(nc)
  out
}


.pp_res_new_grid <- function(
    fit,
    coords,
    cell_lower,
    cell_upper,
    tile_index,
    breaks,
    source
) {
  axis_names <- .pp_res_axis_names(fit)
  coords <- as.matrix(coords)
  cell_lower <- as.matrix(cell_lower)
  cell_upper <- as.matrix(cell_upper)
  tile_index <- as.matrix(tile_index)

  if (!identical(colnames(coords), axis_names) ||
      !identical(colnames(cell_lower), axis_names) ||
      !identical(colnames(cell_upper), axis_names) ||
      !identical(colnames(tile_index), axis_names)) {
    stop("Internal error: residual-grid axes are inconsistent.", call. = FALSE)
  }

  cell_volume <- apply(cell_upper - cell_lower, 1L, prod)

  if (any(!is.finite(cell_volume)) || any(cell_volume <= 0)) {
    stop("Internal error: residual-grid cell volumes are invalid.", call. = FALSE)
  }

  out <- list(
    coords = coords,
    cell_lower = cell_lower,
    cell_upper = cell_upper,
    cell_volume = as.numeric(cell_volume),
    tile_index = tile_index,
    breaks = breaks,
    resolution = stats::setNames(
      vapply(breaks, function(z) length(z) - 1L, integer(1)),
      axis_names
    ),
    bounds = rbind(
      lower = vapply(breaks, min, numeric(1)),
      upper = vapply(breaks, max, numeric(1))
    ),
    window = fit$model_spec$window,
    axis_info = fit$model_spec$axis_info,
    time_info = fit$model_spec$time_info,
    crs = fit$model_spec$crs,
    source = source,
    n_cells = nrow(coords)
  )

  class(out) <- "pp_residual_grid"
  out
}


.pp_res_regular_grid <- function(fit, nc = NULL, region = NULL, max_cells = 1e6) {
  nc <- .pp_res_normalize_nc(fit, nc = nc, max_cells = max_cells)
  axis_names <- names(nc)
  region <- .pp_res_validate_region(fit, region)

  index <- .pp_res_index_grid(nc)
  lower <- region["lower", axis_names]
  upper <- region["upper", axis_names]
  widths <- upper - lower

  coords <- matrix(
    NA_real_,
    nrow = nrow(index),
    ncol = length(axis_names),
    dimnames = list(NULL, axis_names)
  )
  cell_lower <- coords
  cell_upper <- coords
  breaks <- vector("list", length(axis_names))
  names(breaks) <- axis_names

  for (j in seq_along(axis_names)) {
    axis_name <- axis_names[j]
    step <- widths[j] / nc[j]
    edge <- lower[j] + seq.int(0L, nc[j]) * step
    breaks[[axis_name]] <- as.numeric(edge)

    cell_lower[, j] <- edge[index[, j]]
    cell_upper[, j] <- edge[index[, j] + 1L]
    coords[, j] <- (cell_lower[, j] + cell_upper[, j]) / 2
  }

  .pp_res_new_grid(
    fit = fit,
    coords = coords,
    cell_lower = cell_lower,
    cell_upper = cell_upper,
    tile_index = index,
    breaks = breaks,
    source = "regular_cells"
  )
}


.pp_res_grid_from_pp_grid <- function(fit, grid, max_cells = 1e6) {
  if (!inherits(grid, "pp_grid")) {
    stop("`grid` must inherit from class 'pp_grid'.", call. = FALSE)
  }

  axis_names <- .pp_res_axis_names(fit)

  if (is.null(grid$coords) || !is.matrix(grid$coords) ||
      !identical(colnames(grid$coords), axis_names)) {
    stop(
      "The supplied `pp_grid` does not have the fitted coordinate axes.",
      call. = FALSE
    )
  }

  if (!is.null(grid$window) &&
      !isTRUE(all.equal(grid$window$bounds[, axis_names, drop = FALSE],
                         fit$model_spec$window$bounds[, axis_names, drop = FALSE]))) {
    stop("The supplied `pp_grid` must use the same window as `fit`.", call. = FALSE)
  }

  grid_coords <- grid$coords[, axis_names, drop = FALSE]
  axis_values <- lapply(
    seq_along(axis_names),
    function(j) sort(unique(as.numeric(grid_coords[, j])))
  )
  names(axis_values) <- axis_names
  resolution <- vapply(axis_values, length, integer(1))

  if (any(resolution < 1L)) {
    stop("The supplied `pp_grid` has an empty axis.", call. = FALSE)
  }

  n_cells <- prod(as.double(resolution))
  max_cells <- .pp_res_validate_max_cells(max_cells)

  if (n_cells > max_cells) {
    stop(
      sprintf(
        "The supplied `pp_grid` implies %.0f cells, exceeding `max_cells = %.0f`.",
        n_cells,
        as.double(max_cells)
      ),
      call. = FALSE
    )
  }

  row_key <- do.call(
    paste,
    c(as.data.frame(grid_coords, check.names = FALSE), sep = "\r")
  )

  if (nrow(grid_coords) != n_cells || anyDuplicated(row_key)) {
    stop(
      paste0(
        "`pp_grid` must contain one unique point for every Cartesian combination ",
        "of its axis values."
      ),
      call. = FALSE
    )
  }

  bounds <- fit$model_spec$window$bounds[, axis_names, drop = FALSE]
  index <- .pp_res_index_grid(stats::setNames(resolution, axis_names))
  coords <- matrix(
    NA_real_,
    nrow = nrow(index),
    ncol = length(axis_names),
    dimnames = list(NULL, axis_names)
  )
  cell_lower <- coords
  cell_upper <- coords
  breaks <- vector("list", length(axis_names))
  names(breaks) <- axis_names

  for (j in seq_along(axis_names)) {
    axis_name <- axis_names[j]
    values <- axis_values[[axis_name]]
    lower <- bounds["lower", axis_name]
    upper <- bounds["upper", axis_name]

    tolerance <- .Machine$double.eps^0.5 * max(1, abs(c(lower, upper)))

    if (values[1L] < lower - tolerance || values[length(values)] > upper + tolerance) {
      stop(
        sprintf(
          "The supplied `pp_grid` has values outside the fitted window on axis `%s`.",
          axis_name
        ),
        call. = FALSE
      )
    }

    edge <- c(lower, (values[-1L] + values[-length(values)]) / 2, upper)
    edge <- as.numeric(edge)

    if (any(diff(edge) <= 0)) {
      stop(
        sprintf("Could not derive non-overlapping cells from `pp_grid` axis `%s`.", axis_name),
        call. = FALSE
      )
    }

    breaks[[axis_name]] <- edge
    coords[, j] <- values[index[, j]]
    cell_lower[, j] <- edge[index[, j]]
    cell_upper[, j] <- edge[index[, j] + 1L]
  }

  .pp_res_new_grid(
    fit = fit,
    coords = coords,
    cell_lower = cell_lower,
    cell_upper = cell_upper,
    tile_index = index,
    breaks = breaks,
    source = "pp_grid"
  )
}


.pp_res_validate_grid <- function(fit, grid) {
  if (!inherits(grid, "pp_residual_grid")) {
    stop("`grid` must inherit from class 'pp_residual_grid'.", call. = FALSE)
  }

  axis_names <- .pp_res_axis_names(fit)
  required <- c(
    "coords", "cell_lower", "cell_upper", "cell_volume", "tile_index",
    "breaks", "resolution", "bounds"
  )
  missing <- setdiff(required, names(grid))

  if (length(missing) > 0L) {
    stop(
      sprintf(
        "The supplied residual grid is missing: %s.",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!identical(colnames(grid$coords), axis_names) ||
      !identical(colnames(grid$cell_lower), axis_names) ||
      !identical(colnames(grid$cell_upper), axis_names) ||
      !identical(colnames(grid$tile_index), axis_names)) {
    stop("The residual grid axes are incompatible with `fit`.", call. = FALSE)
  }

  if (length(grid$cell_volume) != nrow(grid$coords) ||
      any(!is.finite(grid$cell_volume)) || any(grid$cell_volume <= 0)) {
    stop("The residual grid contains invalid cell volumes.", call. = FALSE)
  }

  fitted_bounds <- fit$model_spec$window$bounds[, axis_names, drop = FALSE]
  lower_matrix <- matrix(
    rep(fitted_bounds["lower", ], each = nrow(grid$coords)),
    nrow = nrow(grid$coords),
    ncol = length(axis_names),
    dimnames = list(NULL, axis_names)
  )
  upper_matrix <- matrix(
    rep(fitted_bounds["upper", ], each = nrow(grid$coords)),
    nrow = nrow(grid$coords),
    ncol = length(axis_names),
    dimnames = list(NULL, axis_names)
  )
  tolerance <- matrix(
    rep(.Machine$double.eps^0.5 * pmax(1, abs(fitted_bounds["upper", ] - fitted_bounds["lower", ])),
        each = nrow(grid$coords)),
    nrow = nrow(grid$coords),
    ncol = length(axis_names),
    dimnames = list(NULL, axis_names)
  )

  if (any(grid$cell_lower < lower_matrix - tolerance) ||
      any(grid$cell_upper > upper_matrix + tolerance)) {
    stop("The residual grid contains cells outside the fitted observation window.", call. = FALSE)
  }

  invisible(grid)
}


.pp_res_prepare_grid <- function(
    fit,
    grid = NULL,
    nc = NULL,
    region = NULL,
    max_cells = 1e6
) {
  if (!is.null(grid) && !is.null(nc)) {
    stop("Supply either `grid` or `nc`, not both.", call. = FALSE)
  }

  if (!is.null(grid) && !is.null(region)) {
    stop("`region` can be supplied only when `grid` is NULL.", call. = FALSE)
  }

  if (is.null(grid)) {
    return(
      .pp_res_regular_grid(
        fit = fit,
        nc = nc,
        region = region,
        max_cells = max_cells
      )
    )
  }

  if (inherits(grid, "pp_residual_grid")) {
    .pp_res_validate_grid(fit, grid)
    return(grid)
  }

  if (inherits(grid, "pp_grid")) {
    return(.pp_res_grid_from_pp_grid(fit, grid, max_cells = max_cells))
  }

  stop(
    paste0(
      "`grid` must be NULL, a `pp_residual_grid` created by `residual_grid_pp()`, ",
      "or a regular `pp_grid`."
    ),
    call. = FALSE
  )
}


.pp_res_cell_ids <- function(pattern, grid) {
  coords <- pattern$coords[, colnames(grid$coords), drop = FALSE]
  axis_names <- colnames(coords)
  n <- nrow(coords)

  if (n == 0L) {
    return(integer(0))
  }

  index <- matrix(NA_integer_, nrow = n, ncol = length(axis_names))

  for (j in seq_along(axis_names)) {
    axis_name <- axis_names[j]
    edge <- grid$breaks[[axis_name]]

    index[, j] <- findInterval(
      x = coords[, axis_name],
      vec = edge,
      rightmost.closed = TRUE,
      all.inside = TRUE
    )
  }

  multipliers <- c(1L, cumprod(as.integer(grid$resolution[-length(axis_names)])))
  as.integer(1L + (index - 1L) %*% multipliers)
}


.pp_res_poisson_deviance <- function(observed, expected) {
  out <- numeric(length(observed))
  positive <- observed > 0

  contribution <- expected
  contribution[positive] <-
    observed[positive] * log(observed[positive] / expected[positive]) -
    (observed[positive] - expected[positive])

  out <- sign(observed - expected) * sqrt(pmax(0, 2 * contribution))
  out
}


.pp_res_weighted_tabulate <- function(index, weights, nbins) {
  out <- numeric(nbins)

  if (length(index) == 0L) {
    return(out)
  }

  if (length(index) != length(weights)) {
    stop("`index` and `weights` must have the same length.", call. = FALSE)
  }

  valid <- !is.na(index) & is.finite(index) &
    index >= 1L & index <= nbins & is.finite(weights)

  if (!any(valid)) {
    return(out)
  }

  sums <- tapply(
    as.numeric(weights[valid]),
    as.integer(index[valid]),
    sum
  )

  out[as.integer(names(sums))] <- as.numeric(sums)
  out
}


.pp_res_residual_column <- function(type) {
  switch(
    type,
    raw = ".pp_raw",
    pearson = ".pp_pearson",
    deviance = ".pp_deviance",
    stop("Internal error: unknown residual type.", call. = FALSE)
  )
}


.pp_res_expected_counts <- function(
    fit,
    grid,
    mark_n = NULL,
    chunk_size = 100000L,
    max_mark_combinations = 100000L
) {
  lambda <- predict_pp_fit(
    object = fit,
    at = grid$coords,
    type = "ground",
    mark_n = mark_n,
    chunk_size = chunk_size,
    max_mark_combinations = max_mark_combinations
  )

  expected <- lambda * grid$cell_volume

  if (any(!is.finite(expected)) || any(expected <= 0)) {
    stop("Fitted expected counts are non-finite or non-positive.", call. = FALSE)
  }

  list(lambda = lambda, expected = expected)
}


# Diagnostic grids and integrated intensity ----------------------------------

#' Construct a diagnostic grid for first-order residual analysis
#'
#' Creates a regular hyperrectangular partition of the observation window
#' retained by a fitted Poisson point-process model. The grid is used to
#' integrate the fitted intensity and to calculate cell residuals and residual
#' profiles.
#'
#' @details
#' First-order diagnostics compare the observed number of events in each cell
#' with its fitted expected count. Diagnostic grids should normally be coarser
#' than the cubature used for estimation: when cells are too small, expected
#' counts are close to zero and Pearson or deviance residual maps may be
#' dominated by sampling noise.
#'
#' If `nc = NULL`, the fitting-cubature resolution is automatically coarsened.
#' The default targets approximately five observed events per cell, uses at
#' most ten cells along each axis, and never exceeds the retained fitting
#' resolution. This is an exploratory default rather than a model-selection
#' rule, so sensitivity to alternative resolutions should be considered.
#'
#' A scalar `nc` is repeated over all coordinate axes. Alternatively, provide
#' one value per axis or a vector named exactly by the fitted axes. The current
#' implementation assumes a finite hyperrectangular `pp_window`.
#'
#' @param fit A fitted model inheriting from `pp_fit`, as returned by
#'   [fit_pp()]. The fit must retain a finite observation window and its current
#'   `model_spec`.
#' @param nc Optional positive integer grid resolution. It may be scalar, have
#'   one entry per fitted coordinate axis, or be named exactly by the axes. If
#'   `NULL`, an automatically coarsened resolution is used.
#' @param max_cells Positive integer giving the maximum permitted total number
#'   of cells. The product of the axis-specific resolutions must not exceed this
#'   value.
#'
#' @return An object of class `pp_residual_grid` containing cell midpoints in
#'   `coords`, axis-specific `breaks`, `resolution`, `cell_volume`, `n_cells`,
#'   covered `region`, coordinate metadata, and the grid-construction `source`.
#'
#' @examples
#' w <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' x <- simulate_pp(
#'   window = w,
#'   expected_n = 60,
#'   output = "pattern",
#'   seed = 21
#' )
#' fit <- fit_pp(x, formula = ~ x + y, nc = 6)
#'
#' g_default <- residual_grid_pp(fit)
#' print(g_default)
#'
#' g_custom <- residual_grid_pp(fit, nc = c(x = 5, y = 4))
#' g_custom$resolution
#'
#' @references
#' Baddeley, A., Rubak, E., and Turner, R. (2015).
#' *Spatial Point Patterns: Methodology and Applications with R*.
#' Chapman and Hall/CRC.
#'
#' @seealso [fit_pp()], [integrate_pp_fit()], [residuals_pp()],
#'   [residual_profile_pp()].
#' @family first-order point-process diagnostics
#' @export
residual_grid_pp <- function(
    fit,
    nc = NULL,
    max_cells = 1e6
) {
  .pp_res_validate_fit(fit)

  .pp_res_regular_grid(
    fit = fit,
    nc = nc,
    region = NULL,
    max_cells = max_cells
  )
}


#' @param x An object of class `pp_residual_grid`.
#' @param ... Additional arguments. They are currently ignored.
#'
#' @return Invisibly returns `x`.
#'
#' @rdname residual_grid_pp
#' @method print pp_residual_grid
#' @export
print.pp_residual_grid <- function(x, ...) {
  cat("<pp_residual_grid>\n")
  cat("First-order diagnostic grid\n\n")
  cat("Cells: ", x$n_cells, "\n", sep = "")
  cat(
    "Resolution: ",
    paste(paste0(names(x$resolution), " = ", x$resolution), collapse = ", "),
    "\n",
    sep = ""
  )
  cat("Source: ", x$source, "\n", sep = "")
  invisible(x)
}


#' Integrate a fitted ground intensity
#'
#' Numerically approximates the fitted mean measure over the complete
#' observation window or over a user-specified hyperrectangular subregion.
#'
#' @details
#' For a region \eqn{B} and fitted ground intensity
#' \eqn{\hat\lambda_g(u)}, the expected event count is
#' \deqn{\hat\Lambda(B)=\int_B \hat\lambda_g(u)\,du.}
#' The function applies the midpoint rule on a regular grid,
#' \deqn{\hat\Lambda(B)\approx\sum_j \hat\lambda_g(c_j)|B_j|,}
#' where \eqn{c_j} and \eqn{|B_j|} are the midpoint and hypervolume of cell
#' \eqn{B_j}.
#'
#' The **ground-process intensity** is always integrated. For a joint marked
#' model, [predict_pp_fit()] first marginalises the fitted joint intensity over
#' the declared continuous and categorical mark domains. Increasing `nc` or
#' `mark_n` may improve numerical accuracy but increases computation and memory
#' use, especially in high dimension.
#'
#' A supplied `grid` is used directly and cannot be combined with `region` or
#' `nc`. Otherwise a grid is constructed over `region`; if `region = NULL`, the
#' complete fitted window is used.
#'
#' @param fit A fitted model inheriting from `pp_fit`.
#' @param region Optional two-row numeric matrix or data frame. Columns must be
#'   named exactly as the fitted coordinate axes; the first and second rows give
#'   lower and upper bounds. The region must lie inside the fitted window.
#' @param grid Optional `pp_residual_grid` or compatible `pp_grid`. It cannot be
#'   combined with `region` or `nc`.
#' @param nc Optional positive integer grid resolution used when `grid = NULL`.
#'   It may be scalar, have one entry per coordinate axis, or be named by axes.
#' @param mark_n Optional positive integer integration resolution for each
#'   continuous joint mark.
#' @param chunk_size Positive integer maximum number of joint-intensity
#'   evaluations performed in one block.
#' @param max_mark_combinations Positive integer maximum number of combined
#'   mark states allowed during ground-intensity marginalisation.
#' @param max_cells Positive integer maximum number of coordinate-domain cells.
#' @param return_grid Logical. If `TRUE`, also return the grid, midpoint
#'   intensities, and cellwise expected counts.
#'
#' @return If `return_grid = FALSE`, an unnamed numeric scalar. Otherwise, a
#'   list with components `integral`, `grid`, `lambda`, and `expected`.
#'
#' @examples
#' w <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' x <- simulate_pp(
#'   window = w,
#'   expected_n = 70,
#'   output = "pattern",
#'   seed = 22
#' )
#' fit <- fit_pp(x, formula = ~ x + y, nc = 6)
#'
#' integrate_pp_fit(fit, nc = 8)
#'
#' left_half <- rbind(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 0.5, y = 1)
#' )
#' integrate_pp_fit(fit, region = left_half, nc = c(x = 4, y = 8))
#'
#' integration <- integrate_pp_fit(fit, nc = 6, return_grid = TRUE)
#' head(integration$expected)
#'
#' @references
#' Baddeley, A., Rubak, E., and Turner, R. (2015).
#' *Spatial Point Patterns: Methodology and Applications with R*.
#' Chapman and Hall/CRC.
#'
#' @seealso [predict_pp_fit()], [residual_grid_pp()], [residuals_pp()].
#' @family first-order point-process diagnostics
#' @export
integrate_pp_fit <- function(
    fit,
    region = NULL,
    grid = NULL,
    nc = NULL,
    mark_n = NULL,
    chunk_size = 100000L,
    max_mark_combinations = 100000L,
    max_cells = 1e6,
    return_grid = FALSE
) {
  .pp_res_validate_fit(fit)
  return_grid <- .pp_res_validate_logical(return_grid, "return_grid")

  diagnostic_grid <- .pp_res_prepare_grid(
    fit = fit,
    grid = grid,
    nc = nc,
    region = region,
    max_cells = max_cells
  )

  predicted <- .pp_res_expected_counts(
    fit = fit,
    grid = diagnostic_grid,
    mark_n = mark_n,
    chunk_size = chunk_size,
    max_mark_combinations = max_mark_combinations
  )

  integral <- sum(predicted$expected)

  if (!isTRUE(return_grid)) {
    return(unname(integral))
  }

  list(
    integral = unname(integral),
    grid = diagnostic_grid,
    lambda = predicted$lambda,
    expected = predicted$expected
  )
}


# Cell residuals --------------------------------------------------------------

#' Compute first-order cell residuals
#'
#' Partitions the fitted coordinate domain into hyperrectangular cells and
#' compares observed counts with counts expected under the fitted ground
#' intensity.
#'
#' @details
#' Let \eqn{B_j} be diagnostic cell \eqn{j}, let
#' \eqn{N_j=N(B_j)} be its observed count, and let
#' \deqn{\hat\mu_j=\int_{B_j}\hat\lambda_g(u)\,du}
#' be its fitted expected count. Three residual definitions are calculated:
#'
#' - raw innovation,
#'   \deqn{R_j=N_j-\hat\mu_j;}
#' - Pearson residual,
#'   \deqn{R_j^P=(N_j-\hat\mu_j)/\sqrt{\hat\mu_j};}
#' - signed Poisson-deviance residual,
#'   \deqn{R_j^D=\operatorname{sign}(N_j-\hat\mu_j)
#'   \sqrt{2\{N_j\log(N_j/\hat\mu_j)-(N_j-\hat\mu_j)\}},}
#'   with \eqn{N_j\log(N_j/\hat\mu_j)=0} when \eqn{N_j=0}.
#'
#' Expected counts are obtained by midpoint integration. The fitted intensity is
#' always the ground intensity; joint marked models are therefore marginalised
#' over their declared mark domains before expectations are calculated.
#'
#' Cell residuals are exploratory and depend on grid resolution. Very fine
#' grids may contain many cells with small expected counts and highly variable
#' Pearson or deviance residuals. The default grid is intentionally coarser
#' than the fitting cubature. Formal simulation-based assessment is available
#' through [residual_profile_pp()] and [envelope_residual_pp()].
#'
#' @param fit A fitted model inheriting from `pp_fit`.
#' @param type Residual copied to the generic `.pp_residual` column:
#'   `"pearson"`, `"raw"`, or `"deviance"`. All three definitions are retained
#'   in the returned object.
#' @param grid Optional `pp_residual_grid` or compatible `pp_grid`. If `NULL`, a
#'   regular diagnostic grid is constructed.
#' @param nc Optional positive integer resolution used when `grid = NULL`. It
#'   may be scalar, have one value per axis, or be named exactly by the axes.
#' @param mark_n Optional integration resolution for continuous joint marks.
#' @param chunk_size Positive integer maximum number of joint-intensity
#'   evaluations per block.
#' @param max_mark_combinations Positive integer maximum number of combined
#'   mark states used during marginalisation.
#' @param max_cells Positive integer maximum number of diagnostic cells.
#'
#' @return An object of class `pp_residuals`. Component `data` contains cell
#'   midpoints, observed and expected counts, fitted midpoint intensities, and
#'   all residual definitions. The object also retains `variable_values`,
#'   `grid`, the observed `pattern`, coordinate metadata, summary `diagnostics`,
#'   and numerical `settings`.
#'
#' @examples
#' w <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' x <- simulate_pp(
#'   window = w,
#'   expected_n = 80,
#'   output = "pattern",
#'   seed = 23
#' )
#' fit <- fit_pp(x, formula = ~ x + y, nc = 7)
#'
#' r_pearson <- residuals_pp(fit, type = "pearson")
#' print(r_pearson)
#' summary(r_pearson)
#'
#' head(r_pearson$data[c(
#'   ".pp_observed", ".pp_expected", ".pp_raw",
#'   ".pp_pearson", ".pp_deviance"
#' )])
#'
#' r_raw <- residuals_pp(
#'   fit,
#'   type = "raw",
#'   nc = c(x = 5, y = 5)
#' )
#' plot(r_raw, what = "raw", show_points = TRUE)
#'
#' @references
#' Baddeley, A., Rubak, E., and Møller, J. (2005). Residual analysis for
#' spatial point processes (with discussion). *Journal of the Royal Statistical
#' Society: Series B*, 67, 617--666.
#'
#' Baddeley, A., Rubak, E., and Turner, R. (2015).
#' *Spatial Point Patterns: Methodology and Applications with R*.
#' Chapman and Hall/CRC.
#'
#' @seealso [residual_grid_pp()], [integrate_pp_fit()],
#'   [residual_profile_pp()], [envelope_residual_pp()].
#' @family first-order point-process diagnostics
#' @export
residuals_pp <- function(
    fit,
    type = c("pearson", "raw", "deviance"),
    grid = NULL,
    nc = NULL,
    mark_n = NULL,
    chunk_size = 100000L,
    max_mark_combinations = 100000L,
    max_cells = 1e6
) {
  .pp_res_validate_fit(fit)
  type <- match.arg(type)

  diagnostic_grid <- .pp_res_prepare_grid(
    fit = fit,
    grid = grid,
    nc = nc,
    region = NULL,
    max_cells = max_cells
  )

  predicted <- .pp_res_expected_counts(
    fit = fit,
    grid = diagnostic_grid,
    mark_n = mark_n,
    chunk_size = chunk_size,
    max_mark_combinations = max_mark_combinations
  )

  cell_id <- .pp_res_cell_ids(fit$pattern, diagnostic_grid)
  observed <- tabulate(cell_id, nbins = diagnostic_grid$n_cells)
  expected <- predicted$expected
  raw <- observed - expected
  pearson <- raw / sqrt(expected)
  deviance <- .pp_res_poisson_deviance(observed, expected)

  data <- as.data.frame(diagnostic_grid$coords, check.names = FALSE)
  data$.pp_cell_id <- seq_len(diagnostic_grid$n_cells)
  data$.pp_observed <- as.numeric(observed)
  data$.pp_expected <- as.numeric(expected)
  data$.pp_lambda <- as.numeric(predicted$lambda)
  data$.pp_raw <- as.numeric(raw)
  data$.pp_pearson <- as.numeric(pearson)
  data$.pp_deviance <- as.numeric(deviance)
  data$.pp_residual <- data[[.pp_res_residual_column(type)]]

  source_values <- .pp_diag_evaluate_sources(fit, diagnostic_grid$coords)
  variable_values <- as.data.frame(diagnostic_grid$coords, check.names = FALSE)

  if (ncol(source_values) > 0L) {
    overlap <- intersect(names(variable_values), names(source_values))

    if (length(overlap) > 0L) {
      stop(
        sprintf(
          "Fitted covariate source names conflict with coordinate axes: %s.",
          paste(overlap, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    variable_values <- cbind(variable_values, source_values)
  }

  diagnostics <- list(
    n_cells = diagnostic_grid$n_cells,
    observed_total = sum(observed),
    expected_total = sum(expected),
    raw_total = sum(raw),
    pearson_mean = mean(pearson),
    pearson_sd = stats::sd(pearson),
    deviance_sum = sum(deviance^2),
    maximum_abs_pearson = max(abs(pearson))
  )

  out <- list(
    call = match.call(),
    type = type,
    residual_column = ".pp_residual",
    data = data,
    variable_values = variable_values,
    grid = diagnostic_grid,
    pattern = fit$pattern,
    axis_info = fit$model_spec$axis_info,
    time_info = fit$model_spec$time_info,
    crs = fit$model_spec$crs,
    formula = fit$formula,
    diagnostics = diagnostics,
    settings = list(
      mark_n = mark_n,
      chunk_size = chunk_size,
      max_mark_combinations = max_mark_combinations
    )
  )

  class(out) <- "pp_residuals"
  out
}


#' Methods for cell-residual objects
#'
#' Print, summarise, and plot objects returned by [residuals_pp()].
#'
#' @details
#' The print and summary methods report grid size, observed and expected totals,
#' the total raw innovation, and descriptive summaries of the selected
#' residual. The plotting method displays residuals, counts, or fitted midpoint
#' intensities.
#'
#' One-dimensional domains are shown as stem plots. In two dimensions the full
#' map is displayed. In higher dimensions, two axes are shown and all remaining
#' axes are fixed at median grid values or user-specified values matched to the
#' nearest grid cells.
#'
#' Residual quantities use a divergent scale centred at zero. By default,
#' colour limits are robust to a single extreme cell; set `robust = FALSE` to
#' use the full range or supply `zlim` explicitly.
#'
#' @param x A `pp_residuals` object for `print()` and `plot()`, or a
#'   `summary_pp_residuals` object for `print.summary_pp_residuals()`.
#' @param object A `pp_residuals` object to summarise.
#' @param digits Number of significant digits used when printing.
#' @param axes Optional names of two distinct fitted coordinate axes to display.
#' @param fixed Optional named numeric vector for axes not displayed. Values are
#'   matched to the nearest grid midpoint.
#' @param what Quantity to plot: `"residual"`, `"raw"`, `"pearson"`,
#'   `"deviance"`, `"observed"`, `"expected"`, or `"lambda"`.
#' @param col Optional vector of colours.
#' @param zlim Optional numeric vector of length two giving colour-scale limits.
#' @param robust Logical; use robust automatic limits.
#' @param robust_prob Probability strictly between `0.5` and `1` used to derive
#'   robust limits.
#' @param show_colorbar Logical; display a labelled colour bar.
#' @param show_cell_boundaries Logical; draw diagnostic-cell boundaries.
#' @param show_points Logical; overlay observed events. In a
#'   higher-dimensional slice, only events in the selected fixed cells appear.
#' @param point_pch Plotting symbol for overlaid events.
#' @param point_cex Point-size multiplier for overlaid events.
#' @param main Optional plot title.
#' @param ... Additional arguments. Plotting arguments are passed to
#'   [graphics::plot()] or [graphics::image()]; printing methods ignore them.
#'
#' @return Printing methods invisibly return `x`.
#'   `summary.pp_residuals()` returns a `summary_pp_residuals` object.
#'   `plot.pp_residuals()` invisibly returns the displayed slice data.
#'
#' @examples
#' w <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' x <- simulate_pp(
#'   window = w,
#'   expected_n = 70,
#'   output = "pattern",
#'   seed = 24
#' )
#' fit <- fit_pp(x, formula = ~ x + y, nc = 6)
#' res <- residuals_pp(fit, type = "pearson", nc = 5)
#'
#' print(res)
#' print(summary(res))
#' plot(res, what = "residual")
#' plot(res, what = "expected", robust = FALSE, show_points = TRUE)
#'
#' @seealso [residuals_pp()], [residual_profile_pp()].
#' @family first-order point-process diagnostics
#' @rdname pp_residuals-methods
#' @method print pp_residuals
#' @export
print.pp_residuals <- function(x, digits = 4, ...) {
  cat("<pp_residuals>\n")
  cat("First-order Poisson cell residuals\n\n")
  cat("Type: ", x$type, "\n", sep = "")
  cat("Cells: ", x$diagnostics$n_cells, "\n", sep = "")
  cat("Observed count: ", x$diagnostics$observed_total, "\n", sep = "")
  cat("Expected count: ", formatC(x$diagnostics$expected_total, digits = digits, format = "fg"), "\n", sep = "")
  cat("Total raw residual: ", formatC(x$diagnostics$raw_total, digits = digits, format = "fg"), "\n", sep = "")
  invisible(x)
}


#' @rdname pp_residuals-methods
#' @method summary pp_residuals
#' @export
summary.pp_residuals <- function(object, ...) {
  values <- object$data$.pp_residual

  out <- list(
    type = object$type,
    diagnostics = object$diagnostics,
    residual_summary = summary(values),
    quantiles = stats::quantile(values, probs = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)),
    resolution = object$grid$resolution
  )

  class(out) <- "summary_pp_residuals"
  out
}


#' @rdname pp_residuals-methods
#' @method print summary_pp_residuals
#' @export
print.summary_pp_residuals <- function(x, digits = 4, ...) {
  cat("First-order residual summary\n\n")
  cat("Type: ", x$type, "\n", sep = "")
  cat("Observed total: ", x$diagnostics$observed_total, "\n", sep = "")
  cat("Expected total: ", formatC(x$diagnostics$expected_total, digits = digits, format = "fg"), "\n", sep = "")
  cat("Total raw residual: ", formatC(x$diagnostics$raw_total, digits = digits, format = "fg"), "\n", sep = "")
  cat("Pearson residual SD: ", formatC(x$diagnostics$pearson_sd, digits = digits, format = "fg"), "\n", sep = "")
  cat("Maximum |Pearson residual|: ", formatC(x$diagnostics$maximum_abs_pearson, digits = digits, format = "fg"), "\n\n", sep = "")
  print(x$residual_summary)
  invisible(x)
}


.pp_res_plot_value <- function(x, what) {
  if (identical(what, "residual")) {
    return(list(
      values = x$data$.pp_residual,
      label = paste0(x$type, " residual")
    ))
  }

  column <- switch(
    what,
    raw = ".pp_raw",
    pearson = ".pp_pearson",
    deviance = ".pp_deviance",
    observed = ".pp_observed",
    expected = ".pp_expected",
    lambda = ".pp_lambda",
    stop("Internal error: unknown plot quantity.", call. = FALSE)
  )

  list(values = x$data[[column]], label = what)
}


.pp_res_slice_indices <- function(x, axes, fixed = NULL) {
  axis_names <- x$axis_info$name
  remaining <- setdiff(axis_names, axes)
  data <- x$data

  if (is.null(fixed)) {
    fixed <- numeric(0)
  }

  if (!is.numeric(fixed)) {
    stop("`fixed` must be NULL or a numeric vector.", call. = FALSE)
  }

  if (length(fixed) > 0L && is.null(names(fixed))) {
    stop("`fixed` must be named by fixed coordinate axes.", call. = FALSE)
  }

  if (length(fixed) > 0L) {
    unknown <- setdiff(names(fixed), remaining)

    if (length(unknown) > 0L) {
      stop(
        sprintf(
          "`fixed` contains axes not available for slicing: %s.",
          paste(unknown, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  keep <- rep(TRUE, nrow(data))
  selected <- numeric(length(remaining))
  names(selected) <- remaining

  for (axis_name in remaining) {
    values <- sort(unique(data[[axis_name]]))
    target <- if (axis_name %in% names(fixed)) {
      fixed[[axis_name]]
    } else {
      stats::median(values)
    }

    closest <- values[which.min(abs(values - target))]
    selected[[axis_name]] <- closest
    keep <- keep & abs(data[[axis_name]] - closest) <= .Machine$double.eps^0.5 * max(1, abs(closest))
  }

  list(keep = keep, selected = selected)
}


.pp_res_plot_limits <- function(values, residual_scale, zlim = NULL,
                                robust = TRUE, robust_prob = 0.98) {
  finite_values <- as.numeric(values[is.finite(values)])

  if (length(finite_values) == 0L) {
    stop("The selected slice contains no finite values to plot.", call. = FALSE)
  }

  if (!is.null(zlim)) {
    if (!is.numeric(zlim) || length(zlim) != 2L || any(!is.finite(zlim)) ||
        zlim[1L] >= zlim[2L]) {
      stop("`zlim` must contain two finite values with lower < upper.", call. = FALSE)
    }

    return(as.numeric(zlim))
  }

  if (!is.logical(robust) || length(robust) != 1L || is.na(robust)) {
    stop("`robust` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.numeric(robust_prob) || length(robust_prob) != 1L ||
      !is.finite(robust_prob) || robust_prob <= 0.5 || robust_prob >= 1) {
    stop("`robust_prob` must be a single number strictly between 0.5 and 1.", call. = FALSE)
  }

  if (isTRUE(residual_scale)) {
    limit <- if (isTRUE(robust)) {
      as.numeric(stats::quantile(abs(finite_values), probs = robust_prob, names = FALSE))
    } else {
      max(abs(finite_values))
    }

    if (!is.finite(limit) || limit <= 0) {
      limit <- 1
    }

    return(c(-limit, limit))
  }

  limits <- if (isTRUE(robust)) {
    as.numeric(stats::quantile(
      finite_values,
      probs = c(1 - robust_prob, robust_prob),
      names = FALSE
    ))
  } else {
    range(finite_values)
  }

  if (!all(is.finite(limits)) || limits[1L] == limits[2L]) {
    centre <- if (all(is.finite(limits))) mean(limits) else 0
    spread <- max(0.5, abs(centre) * 0.25)
    limits <- centre + c(-spread, spread)
  }

  limits
}


.pp_res_draw_colorbar <- function(zlim, col, label) {
  n_col <- length(col)
  y <- seq(zlim[1L], zlim[2L], length.out = n_col)

  # The colour bar occupies a narrow layout panel: use small side margins so
  # that it fits on ordinary device sizes (the default margins alone are
  # wider than the panel on a 6-inch device).
  main_mar <- graphics::par("mar")
  old_par <- graphics::par(mar = c(main_mar[1L], 0.4, main_mar[3L], 2.6))
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::image(
    x = 1,
    y = y,
    z = matrix(seq_len(n_col), nrow = 1L, ncol = n_col),
    col = col,
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  graphics::axis(4L, at = pretty(zlim))
  graphics::box()
  graphics::mtext(label, side = 3L, line = 0.25, cex = 0.85)
}


.pp_res_points_in_slice <- function(x, axes, selected) {
  if (is.null(x$pattern) || is.null(x$pattern$coords)) {
    return(NULL)
  }

  points <- x$pattern$coords[, x$axis_info$name, drop = FALSE]
  remaining <- names(selected)

  for (axis_name in remaining) {
    centre <- selected[[axis_name]]
    edge <- x$grid$breaks[[axis_name]]
    cell_index <- which.min(abs((edge[-1L] + edge[-length(edge)]) / 2 - centre))
    lower <- edge[cell_index]
    upper <- edge[cell_index + 1L]

    if (cell_index == length(edge) - 1L) {
      keep <- points[, axis_name] >= lower & points[, axis_name] <= upper
    } else {
      keep <- points[, axis_name] >= lower & points[, axis_name] < upper
    }

    points <- points[keep, , drop = FALSE]
  }

  points[, axes, drop = FALSE]
}


#' @rdname pp_residuals-methods
#' @method plot pp_residuals
#' @export
plot.pp_residuals <- function(
    x,
    axes = NULL,
    fixed = NULL,
    what = c("residual", "raw", "pearson", "deviance", "observed", "expected", "lambda"),
    col = NULL,
    zlim = NULL,
    robust = TRUE,
    robust_prob = 0.98,
    show_colorbar = TRUE,
    show_cell_boundaries = TRUE,
    show_points = FALSE,
    point_pch = 16,
    point_cex = 0.55,
    main = NULL,
    ...
) {
  if (!inherits(x, "pp_residuals")) {
    stop("`x` must inherit from class 'pp_residuals'.", call. = FALSE)
  }

  what <- match.arg(what)
  axis_names <- x$axis_info$name
  d <- length(axis_names)
  quantity <- .pp_res_plot_value(x, what)
  residual_scale <- what %in% c("residual", "raw", "pearson", "deviance")

  if (!is.logical(show_colorbar) || length(show_colorbar) != 1L || is.na(show_colorbar)) {
    stop("`show_colorbar` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(show_cell_boundaries) || length(show_cell_boundaries) != 1L || is.na(show_cell_boundaries)) {
    stop("`show_cell_boundaries` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(show_points) || length(show_points) != 1L || is.na(show_points)) {
    stop("`show_points` must be TRUE or FALSE.", call. = FALSE)
  }

  if (d == 1L) {
    graphics::plot(
      x$data[[axis_names]],
      quantity$values,
      type = "h",
      xlab = axis_names,
      ylab = quantity$label,
      main = if (is.null(main)) "First-order residuals" else main,
      ...
    )
    if (isTRUE(residual_scale)) {
      graphics::abline(h = 0, lty = 2)
    }
    return(invisible(list(data = x$data, axes = axis_names, values = quantity$values)))
  }

  if (is.null(axes)) {
    axes <- axis_names[seq_len(2L)]
  }

  if (!is.character(axes) || length(axes) != 2L || anyDuplicated(axes) ||
      !all(axes %in% axis_names)) {
    stop("`axes` must name two distinct fitted coordinate axes.", call. = FALSE)
  }

  slice <- .pp_res_slice_indices(x, axes = axes, fixed = fixed)
  data <- x$data[slice$keep, , drop = FALSE]
  values <- quantity$values[slice$keep]

  x_values <- sort(unique(data[[axes[1L]]]))
  y_values <- sort(unique(data[[axes[2L]]]))
  z <- matrix(NA_real_, nrow = length(x_values), ncol = length(y_values))
  ix <- match(data[[axes[1L]]], x_values)
  iy <- match(data[[axes[2L]]], y_values)
  z[cbind(ix, iy)] <- values

  finite_z <- z[is.finite(z)]
  if (length(finite_z) == 0L) {
    stop("The selected slice contains no finite values to plot.", call. = FALSE)
  }

  if (is.null(col)) {
    if (isTRUE(residual_scale)) {
      col <- grDevices::hcl.colors(101L, palette = "Blue-Red 3", rev = TRUE)
    } else {
      col <- grDevices::hcl.colors(101L, palette = "YlOrRd", rev = FALSE)
    }
  }

  zlim <- .pp_res_plot_limits(
    values = finite_z,
    residual_scale = residual_scale,
    zlim = zlim,
    robust = robust,
    robust_prob = robust_prob
  )

  slice_label <- if (length(slice$selected) == 0L) {
    ""
  } else {
    paste(
      paste0(names(slice$selected), " = ", formatC(slice$selected, digits = 4, format = "fg")),
      collapse = ", "
    )
  }

  default_main <- paste0("First-order ", quantity$label, " map")
  if (nzchar(slice_label)) {
    default_main <- paste0(default_main, " (", slice_label, ")")
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  if (isTRUE(show_colorbar)) {
    graphics::layout(matrix(c(1L, 2L), nrow = 1L), widths = c(4.7, 1.15))
    on.exit(graphics::layout(matrix(1L, nrow = 1L, ncol = 1L)), add = TRUE)
  }

  x_breaks <- x$grid$breaks[[axes[1L]]]
  y_breaks <- x$grid$breaks[[axes[2L]]]

  graphics::image(
    x = x_values,
    y = y_values,
    z = z,
    col = col,
    zlim = zlim,
    xlim = range(x_breaks),
    ylim = range(y_breaks),
    xlab = axes[1L],
    ylab = axes[2L],
    main = if (is.null(main)) default_main else main,
    useRaster = FALSE,
    ...
  )

  if (isTRUE(show_cell_boundaries)) {
    graphics::abline(v = x_breaks, h = y_breaks, col = grDevices::adjustcolor("grey20", alpha.f = 0.35), lwd = 0.4)
  }

  if (isTRUE(show_points)) {
    points_slice <- .pp_res_points_in_slice(x, axes = axes, selected = slice$selected)
    if (!is.null(points_slice) && nrow(points_slice) > 0L) {
      graphics::points(
        points_slice[, 1L],
        points_slice[, 2L],
        pch = point_pch,
        cex = point_cex,
        col = grDevices::adjustcolor("black", alpha.f = 0.70)
      )
    }
  }

  if (isTRUE(show_colorbar)) {
    .pp_res_draw_colorbar(
      zlim = zlim,
      col = col,
      label = if (isTRUE(residual_scale)) quantity$label else quantity$label
    )
  }

  invisible(list(
    axes = axes,
    fixed = slice$selected,
    x = x_values,
    y = y_values,
    z = z,
    values = values,
    zlim = zlim
  ))
}


# Residual profiles -----------------------------------------------------------

.pp_res_profile_variables <- function(residuals, variable) {
  available <- names(residuals$variable_values)

  if (is.character(variable) && length(variable) == 1L && identical(variable, "all")) {
    return(available)
  }

  if (!is.character(variable) || length(variable) == 0L ||
      any(!nzchar(variable)) || anyDuplicated(variable)) {
    stop("`variable` must be 'all' or a non-empty character vector.", call. = FALSE)
  }

  missing <- setdiff(variable, available)

  if (length(missing) > 0L) {
    stop(
      sprintf(
        "Unknown profile variable(s): %s.",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  variable
}


.pp_res_validate_profile_breaks <- function(values, breaks, variable) {
  if (!is.numeric(breaks) || length(breaks) < 2L ||
      any(!is.finite(breaks))) {
    stop(
      sprintf("Breaks for profile `%s` must be finite numeric values.", variable),
      call. = FALSE
    )
  }

  breaks <- sort(unique(as.numeric(breaks)))

  if (length(breaks) < 2L || any(diff(breaks) <= 0)) {
    stop(
      sprintf("Breaks for profile `%s` must define at least one non-empty interval.", variable),
      call. = FALSE
    )
  }

  tolerance <- .Machine$double.eps^0.5 * max(1, abs(c(values, breaks)))

  if (min(values) < breaks[1L] - tolerance ||
      max(values) > breaks[length(breaks)] + tolerance) {
    stop(
      sprintf(
        "Breaks for profile `%s` do not cover all diagnostic-grid values.",
        variable
      ),
      call. = FALSE
    )
  }

  breaks
}


.pp_res_default_profile_breaks <- function(values, bins = NULL) {
  values <- as.numeric(values)
  unique_values <- sort(unique(values))

  if (length(unique_values) == 1L) {
    delta <- max(1, abs(unique_values)) * 0.5
    return(c(unique_values - delta, unique_values + delta))
  }

  if (is.null(bins)) {
    bins <- min(20L, length(unique_values))
  }

  if (!is.numeric(bins) || length(bins) != 1L || !is.finite(bins) ||
      bins != as.integer(bins) || bins < 1L) {
    stop("`bins` must be NULL or one positive integer.", call. = FALSE)
  }

  bins <- min(as.integer(bins), length(unique_values))
  probs <- seq(0, 1, length.out = bins + 1L)
  breaks <- unique(as.numeric(stats::quantile(values, probs = probs, names = FALSE, type = 8)))

  if (length(breaks) < 2L) {
    breaks <- c(min(values), max(values))
  }

  breaks[1L] <- min(values)
  breaks[length(breaks)] <- max(values)
  breaks
}


.pp_res_resolve_breaks <- function(variable, values, breaks, bins, n_variables) {
  if (is.null(breaks)) {
    return(.pp_res_default_profile_breaks(values, bins = bins))
  }

  if (is.list(breaks)) {
    if (is.null(names(breaks)) || !variable %in% names(breaks)) {
      stop(
        paste0(
          "When `breaks` is a list, it must be named and contain each requested ",
          "profile variable."
        ),
        call. = FALSE
      )
    }

    return(.pp_res_validate_profile_breaks(values, breaks[[variable]], variable))
  }

  if (n_variables > 1L) {
    stop("With multiple profile variables, `breaks` must be a named list.", call. = FALSE)
  }

  .pp_res_validate_profile_breaks(values, breaks, variable)
}


.pp_res_make_profile <- function(
    residuals,
    variable,
    form,
    residual_type,
    breaks
) {
  values <- residuals$variable_values[[variable]]

  if (!is.numeric(values) || any(!is.finite(values))) {
    stop(
      sprintf("Profile variable `%s` must be finite and numeric.", variable),
      call. = FALSE
    )
  }

  index <- findInterval(
    x = values,
    vec = breaks,
    rightmost.closed = TRUE,
    all.inside = TRUE
  )
  n_bins <- length(breaks) - 1L
  data <- residuals$data

  observed <- .pp_res_weighted_tabulate(index, data$.pp_observed, n_bins)
  expected <- .pp_res_weighted_tabulate(index, data$.pp_expected, n_bins)
  raw <- .pp_res_weighted_tabulate(index, data$.pp_raw, n_bins)
  pearson <- .pp_res_weighted_tabulate(index, data$.pp_pearson, n_bins)
  deviance <- .pp_res_weighted_tabulate(index, data$.pp_deviance, n_bins)
  n_cells <- tabulate(index, nbins = n_bins)

  if (identical(form, "cumulative")) {
    observed <- cumsum(observed)
    expected <- cumsum(expected)
    raw <- cumsum(raw)
    pearson <- cumsum(pearson)
    deviance <- cumsum(deviance)
    n_cells <- cumsum(n_cells)
    profile_index <- breaks[-1L]
  } else {
    profile_index <- vapply(
      seq_len(n_bins),
      function(j) {
        in_bin <- index == j
        if (any(in_bin)) mean(values[in_bin]) else mean(breaks[j:(j + 1L)])
      },
      numeric(1)
    )
  }

  selected <- switch(
    residual_type,
    raw = raw,
    pearson = pearson,
    deviance = deviance
  )

  out <- list(
    call = match.call(),
    variable = variable,
    form = form,
    residual_type = residual_type,
    breaks = breaks,
    data = data.frame(
      index = profile_index,
      n_cells = as.integer(n_cells),
      observed = observed,
      expected = expected,
      raw = raw,
      pearson = pearson,
      deviance = deviance,
      residual = selected,
      check.names = FALSE
    )
  )

  class(out) <- "pp_residual_profile"
  out
}


#' Construct first-order residual profiles
#'
#' Aggregates cellwise residual contributions along one or more coordinate axes
#' or fitted numeric covariates, producing cumulative or binned one-dimensional
#' diagnostic profiles.
#'
#' @details
#' Residual profiles reduce a multidimensional residual field to curves indexed
#' by a coordinate or covariate. For `form = "cumulative"`, cell contributions
#' are accumulated up to successive upper breakpoints. Persistent positive or
#' negative excursions can reveal ranges where the model systematically
#' underpredicts or overpredicts events. For `form = "binned"`, contributions
#' are summed separately within each interval.
#'
#' If `breaks = NULL`, intervals are based on empirical quantiles of the
#' cell-level variable. Explicit breaks must be finite, strictly increasing,
#' and span the available values. With several variables, explicit breaks must
#' be supplied as a named list.
#'
#' If `x` is a `pp_fit`, [residuals_pp()] is called first. If `x` is already a
#' `pp_residuals` object, its grid and cell contributions are reused.
#' `variable = "all"` includes all fitted coordinate axes and all available
#' finite numeric fitted-covariate sources.
#'
#' Profiles are exploratory unless calibrated against simulations. Use
#' [envelope_residual_pp()] or [diagnose_first_order_pp()] for plug-in refit
#' global-envelope diagnostics.
#'
#' @param x A `pp_fit` object or an object returned by [residuals_pp()].
#' @param variable Character vector naming coordinate axes or fitted numeric
#'   covariates. The special value `"all"` requests every available variable.
#' @param form `"cumulative"` or `"binned"`.
#' @param residual Residual contribution to aggregate: `"raw"`, `"pearson"`, or
#'   `"deviance"`.
#' @param bins Optional positive integer number of quantile-based intervals used
#'   when `breaks = NULL`.
#' @param breaks Optional numeric break vector for one variable, or a named list
#'   of numeric break vectors when several variables are requested.
#' @param grid Optional `pp_residual_grid` or compatible `pp_grid`, used only
#'   when `x` is a `pp_fit`.
#' @param nc Optional diagnostic-grid resolution used only when `x` is a fit and
#'   `grid = NULL`.
#' @param mark_n Optional continuous joint-mark integration resolution passed to
#'   [residuals_pp()] when `x` is a fit.
#' @param chunk_size Positive integer block size for joint-intensity evaluation.
#' @param max_mark_combinations Positive integer limit on combined mark states.
#' @param max_cells Positive integer maximum number of diagnostic cells.
#'
#' @return For one variable, an object of class `pp_residual_profile` containing
#'   the variable, profile form, residual type, breaks, and a data frame with
#'   profile indices, observed and expected counts, all residual aggregates,
#'   and the selected curve. For several variables, an object of class
#'   `pp_residual_profile_set` containing a named list of profiles.
#'
#' @examples
#' w <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' x <- simulate_pp(
#'   window = w,
#'   expected_n = 90,
#'   output = "pattern",
#'   seed = 25
#' )
#' fit <- fit_pp(x, formula = ~ x + y, nc = 7)
#'
#' px <- residual_profile_pp(
#'   fit,
#'   variable = "x",
#'   form = "cumulative",
#'   residual = "raw",
#'   nc = 6
#' )
#' print(px)
#' print(summary(px))
#' plot(px)
#'
#' res <- residuals_pp(fit, type = "pearson", nc = 5)
#' profiles <- residual_profile_pp(
#'   res,
#'   variable = "all",
#'   form = "binned",
#'   residual = "pearson",
#'   bins = 5
#' )
#' print(profiles)
#'
#' @references
#' Baddeley, A., Rubak, E., and Møller, J. (2005). Residual analysis for
#' spatial point processes (with discussion). *Journal of the Royal Statistical
#' Society: Series B*, 67, 617--666.
#'
#' Baddeley, A., Rubak, E., and Turner, R. (2015).
#' *Spatial Point Patterns: Methodology and Applications with R*.
#' Chapman and Hall/CRC.
#'
#' @seealso [residuals_pp()], [envelope_residual_pp()],
#'   [diagnose_first_order_pp()].
#' @family first-order point-process diagnostics
#' @export
residual_profile_pp <- function(
    x,
    variable = "all",
    form = c("cumulative", "binned"),
    residual = c("raw", "pearson", "deviance"),
    bins = NULL,
    breaks = NULL,
    grid = NULL,
    nc = NULL,
    mark_n = NULL,
    chunk_size = 100000L,
    max_mark_combinations = 100000L,
    max_cells = 1e6
) {
  form <- match.arg(form)
  residual <- match.arg(residual)

  residuals <- if (inherits(x, "pp_fit")) {
    residuals_pp(
      fit = x,
      type = residual,
      grid = grid,
      nc = nc,
      mark_n = mark_n,
      chunk_size = chunk_size,
      max_mark_combinations = max_mark_combinations,
      max_cells = max_cells
    )
  } else if (inherits(x, "pp_residuals")) {
    x
  } else {
    stop("`x` must be a `pp_fit` or `pp_residuals` object.", call. = FALSE)
  }

  variables <- .pp_res_profile_variables(residuals, variable)
  profiles <- vector("list", length(variables))
  names(profiles) <- variables

  for (i in seq_along(variables)) {
    variable_i <- variables[i]
    values_i <- residuals$variable_values[[variable_i]]
    breaks_i <- .pp_res_resolve_breaks(
      variable = variable_i,
      values = values_i,
      breaks = breaks,
      bins = bins,
      n_variables = length(variables)
    )

    profiles[[i]] <- .pp_res_make_profile(
      residuals = residuals,
      variable = variable_i,
      form = form,
      residual_type = residual,
      breaks = breaks_i
    )
  }

  if (length(profiles) == 1L) {
    return(profiles[[1L]])
  }

  out <- list(
    call = match.call(),
    profiles = profiles,
    variables = variables,
    form = form,
    residual_type = residual,
    residuals = residuals
  )

  class(out) <- "pp_residual_profile_set"
  out
}


#' Methods for first-order residual profiles
#'
#' Print, summarise, and plot objects returned by [residual_profile_pp()].
#'
#' @details
#' The methods for a single profile report its variable, form, residual type,
#' terminal values, and largest absolute excursion, or draw the residual curve
#' with a zero reference line. Methods for `pp_residual_profile_set` report or
#' plot all stored variables in a multi-panel display.
#'
#' @param x A `pp_residual_profile` or `pp_residual_profile_set` object for
#'   printing and plotting, or a `summary_pp_residual_profile` object for
#'   `print.summary_pp_residual_profile()`.
#' @param object A `pp_residual_profile` object to summarise.
#' @param digits Number of significant digits used when printing.
#' @param main Optional title for an individual profile plot.
#' @param ... Additional arguments. For plotting, they are passed to
#'   [graphics::plot()] and, for a profile set, forwarded to each panel.
#'
#' @return Printing and plotting methods invisibly return `x`.
#'   `summary.pp_residual_profile()` returns a
#'   `summary_pp_residual_profile` object.
#'
#' @examples
#' w <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' x <- simulate_pp(
#'   window = w,
#'   expected_n = 70,
#'   output = "pattern",
#'   seed = 26
#' )
#' fit <- fit_pp(x, formula = ~ x + y, nc = 6)
#'
#' profile_x <- residual_profile_pp(fit, variable = "x", nc = 5)
#' print(profile_x)
#' print(summary(profile_x))
#' plot(profile_x)
#'
#' profile_set <- residual_profile_pp(
#'   fit,
#'   variable = c("x", "y"),
#'   form = "binned",
#'   bins = 4,
#'   nc = 5
#' )
#' print(profile_set)
#' plot(profile_set)
#'
#' @seealso [residual_profile_pp()], [envelope_residual_pp()].
#' @family first-order point-process diagnostics
#' @rdname pp_residual_profile-methods
#' @method print pp_residual_profile
#' @export
print.pp_residual_profile <- function(x, ...) {
  cat("<pp_residual_profile>\n")
  cat("Variable: ", x$variable, "\n", sep = "")
  cat("Form: ", x$form, "\n", sep = "")
  cat("Residual contribution: ", x$residual_type, "\n", sep = "")
  cat("Profile points: ", nrow(x$data), "\n", sep = "")
  invisible(x)
}


#' @rdname pp_residual_profile-methods
#' @method summary pp_residual_profile
#' @export
summary.pp_residual_profile <- function(object, ...) {
  out <- list(
    variable = object$variable,
    form = object$form,
    residual_type = object$residual_type,
    final_observed = utils::tail(object$data$observed, 1L),
    final_expected = utils::tail(object$data$expected, 1L),
    final_residual = utils::tail(object$data$residual, 1L),
    maximum_abs_residual = max(abs(object$data$residual))
  )

  class(out) <- "summary_pp_residual_profile"
  out
}


#' @rdname pp_residual_profile-methods
#' @method print summary_pp_residual_profile
#' @export
print.summary_pp_residual_profile <- function(x, digits = 4, ...) {
  cat("Residual profile summary\n\n")
  cat("Variable: ", x$variable, "\n", sep = "")
  cat("Form: ", x$form, "\n", sep = "")
  cat("Final observed count: ", formatC(x$final_observed, digits = digits, format = "fg"), "\n", sep = "")
  cat("Final expected count: ", formatC(x$final_expected, digits = digits, format = "fg"), "\n", sep = "")
  cat("Final residual: ", formatC(x$final_residual, digits = digits, format = "fg"), "\n", sep = "")
  cat("Maximum absolute residual: ", formatC(x$maximum_abs_residual, digits = digits, format = "fg"), "\n", sep = "")
  invisible(x)
}


#' @rdname pp_residual_profile-methods
#' @method plot pp_residual_profile
#' @export
plot.pp_residual_profile <- function(x, main = NULL, ...) {
  if (!inherits(x, "pp_residual_profile")) {
    stop("`x` must inherit from class 'pp_residual_profile'.", call. = FALSE)
  }

  graphics::plot(
    x$data$index,
    x$data$residual,
    type = "l",
    xlab = x$variable,
    ylab = paste0(x$residual_type, " residual"),
    main = if (is.null(main)) {
      paste0(x$form, " residual profile: ", x$variable)
    } else {
      main
    },
    ...
  )
  graphics::abline(h = 0, lty = 2)
  invisible(x)
}


#' @rdname pp_residual_profile-methods
#' @method print pp_residual_profile_set
#' @export
print.pp_residual_profile_set <- function(x, ...) {
  cat("<pp_residual_profile_set>\n")
  cat("Variables: ", paste(x$variables, collapse = ", "), "\n", sep = "")
  cat("Form: ", x$form, "\n", sep = "")
  cat("Residual contribution: ", x$residual_type, "\n", sep = "")
  invisible(x)
}


#' @rdname pp_residual_profile-methods
#' @method plot pp_residual_profile_set
#' @export
plot.pp_residual_profile_set <- function(x, ...) {
  if (!inherits(x, "pp_residual_profile_set")) {
    stop("`x` must inherit from class 'pp_residual_profile_set'.", call. = FALSE)
  }

  n_profiles <- length(x$profiles)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  ncol <- ceiling(sqrt(n_profiles))
  nrow <- ceiling(n_profiles / ncol)
  graphics::par(mfrow = c(nrow, ncol))

  for (profile in x$profiles) {
    plot(profile, ...)
  }

  invisible(x)
}


# Global envelopes for first-order residual profiles -------------------------

.pp_res_profile_list <- function(x) {
  if (inherits(x, "pp_residual_profile")) {
    return(stats::setNames(list(x), x$variable))
  }

  if (inherits(x, "pp_residual_profile_set")) {
    return(x$profiles)
  }

  stop("Internal error: invalid residual-profile object.", call. = FALSE)
}


.pp_res_get_p_interval <- function(envelope) {
  p_interval <- attr(envelope, "p_interval", exact = TRUE)

  if (is.null(p_interval)) {
    p_interval <- envelope$p_interval
  }

  if (!is.null(p_interval)) {
    p_interval <- as.numeric(p_interval)
    if (length(p_interval) == 2L && all(is.finite(p_interval))) {
      return(stats::setNames(sort(p_interval), c("lower", "upper")))
    }
  }

  p <- attr(envelope, "p", exact = TRUE)

  if (is.null(p)) {
    p <- envelope$p
  }

  if (is.null(p)) {
    p <- envelope$p.value
  }

  if (!is.null(p)) {
    p <- as.numeric(p)
    if (length(p) == 1L && is.finite(p)) {
      return(stats::setNames(rep(p, 2L), c("lower", "upper")))
    }
  }

  stats::setNames(c(NA_real_, NA_real_), c("lower", "upper"))
}


.pp_res_validate_simulation_args <- function(simulation_args) {
  if (!is.list(simulation_args) ||
      (length(simulation_args) > 0L &&
       (is.null(names(simulation_args)) || any(!nzchar(names(simulation_args)))))) {
    stop("`simulation_args` must be a named list or an empty list.", call. = FALSE)
  }

  prohibited <- intersect(names(simulation_args), c("fit", "nsim", "seed", "output"))

  if (length(prohibited) > 0L) {
    stop(
      sprintf(
        "`simulation_args` cannot contain: %s.",
        paste(prohibited, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  simulation_args
}


#' Global-envelope tests for first-order residual profiles
#'
#' Compares observed residual profiles with profiles generated by repeated
#' simulation and refitting under a fitted Poisson point-process model.
#'
#' @details
#' The function implements a plug-in parametric-bootstrap diagnostic:
#'
#' 1. calculate observed cell residuals and profiles;
#' 2. simulate a pattern from the fitted model with [simulate_pp_fit()];
#' 3. refit the retained model specification with [refit_pp()];
#' 4. calculate simulated profiles on the observed grid and breakpoints;
#' 5. repeat `nsim` times and apply [GET::global_envelope_test()].
#'
#' Refitting each simulation accounts for the main effect of parameter
#' estimation. The implementation is nevertheless a **plug-in refit**
#' procedure: it does not perform the additional nested simulation required for
#' exact composite-hypothesis adjustment. Envelopes and p-values should
#' therefore be interpreted as practical model diagnostics.
#'
#' With `variable = "all"`, a separate test is returned for every fitted
#' coordinate and available finite numeric fitted covariate. All simulated
#' curves use the observed diagnostic grid and profile breaks.
#'
#' The procedure can be computationally expensive. Its cost grows with `nsim`,
#' cubature size, number of variables, and mark marginalisation. It is not
#' suitable for models with substantial probability of an empty simulated
#' pattern, because the current [fit_pp()] cannot refit an empty pattern.
#'
#' @param fit A fitted model inheriting from `pp_fit` and retaining a reusable
#'   `model_spec`.
#' @param variable Character vector naming profile variables, or `"all"`.
#' @param nsim Integer at least two giving the number of primary simulations and
#'   refits. Larger values give finer p-value resolution.
#' @param grid Optional `pp_residual_grid` or compatible `pp_grid` used for the
#'   observed and all simulated residuals.
#' @param nc Optional diagnostic-grid resolution used when `grid = NULL`.
#' @param form `"cumulative"` or `"binned"`.
#' @param residual Residual contribution: `"raw"`, `"pearson"`, or
#'   `"deviance"`.
#' @param bins Optional number of quantile-based profile intervals.
#' @param breaks Optional break vector for one variable or named list for
#'   several variables, passed to [residual_profile_pp()].
#' @param alpha Significance level strictly between zero and one.
#' @param type Global-envelope ordering passed to
#'   [GET::global_envelope_test()], for example `"erl"`.
#' @param alternative Alternative passed to [GET::global_envelope_test()]:
#'   `"two.sided"`, `"less"`, or `"greater"`.
#' @param mark_n Optional integration resolution for continuous joint marks.
#' @param chunk_size Positive integer maximum number of joint-intensity
#'   evaluations per block.
#' @param max_mark_combinations Positive integer maximum number of combined mark
#'   states.
#' @param max_cells Positive integer maximum number of diagnostic cells.
#' @param simulation_args Optional named list of additional arguments passed to
#'   [simulate_pp_fit()]. It cannot contain `fit`, `nsim`, `seed`, or `output`.
#' @param keep_simulations Logical; retain simulated `pp` patterns.
#' @param keep_refits Logical; retain all simulated `pp_fit` refits.
#' @param seed Optional finite integer random seed for the full simulation
#'   sequence.
#' @param verbose Logical; print simulation progress.
#' @param ... Additional named arguments passed to
#'   [GET::global_envelope_test()].
#'
#' @return An object of class `pp_residual_envelope`. It contains observed
#'   residuals and profiles, a named `envelopes` list with `GET` curve sets,
#'   envelope results and simulated curves, retained settings, and optionally
#'   simulated patterns and refits.
#'
#' @examples
#' w <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' x <- simulate_pp(
#'   window = w,
#'   expected_n = 80,
#'   output = "pattern",
#'   seed = 27
#' )
#' fit <- fit_pp(x, formula = ~ x + y, nc = 6)
#'
#' # A substantive analysis should use substantially more simulations.
#' \donttest{
#' env <- envelope_residual_pp(
#'   fit,
#'   variable = "x",
#'   nsim = 19,
#'   nc = 5,
#'   form = "cumulative",
#'   residual = "raw",
#'   seed = 28,
#'   verbose = FALSE
#' )
#' print(env)
#' print(summary(env))
#' plot(env, variable = "x")
#' }
#'
#' @references
#' Baddeley, A., Rubak, E., and Møller, J. (2005). Residual analysis for
#' spatial point processes (with discussion). *Journal of the Royal Statistical
#' Society: Series B*, 67, 617--666.
#'
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H., and Hahn, U. (2017).
#' Global envelope tests for spatial processes. *Journal of the Royal
#' Statistical Society: Series B*, 79, 381--404.
#'
#' @seealso [residual_profile_pp()], [simulate_pp_fit()], [refit_pp()],
#'   [diagnose_first_order_pp()], [GET::global_envelope_test()].
#' @family first-order point-process diagnostics
#' @export
envelope_residual_pp <- function(
    fit,
    variable = "all",
    nsim = 199L,
    grid = NULL,
    nc = NULL,
    form = c("cumulative", "binned"),
    residual = c("raw", "pearson", "deviance"),
    bins = NULL,
    breaks = NULL,
    alpha = 0.05,
    type = "erl",
    alternative = c("two.sided", "less", "greater"),
    mark_n = NULL,
    chunk_size = 100000L,
    max_mark_combinations = 100000L,
    max_cells = 1e6,
    simulation_args = list(),
    keep_simulations = FALSE,
    keep_refits = FALSE,
    seed = NULL,
    verbose = TRUE,
    ...
) {
  .pp_res_require(c("simulate_pp_fit", "refit_pp"), "Residual-profile envelopes")
  .pp_res_validate_fit(fit)

  if (!requireNamespace("GET", quietly = TRUE)) {
    stop(
      "Package 'GET' is required for residual-profile envelopes; add it to Imports and install it.",
      call. = FALSE
    )
  }

  nsim <- .pp_diag_validate_positive_integer(nsim, "nsim", min_value = 2L)
  alpha <- .pp_diag_validate_probability(alpha, "alpha")
  form <- match.arg(form)
  residual <- match.arg(residual)
  alternative <- match.arg(alternative)
  seed <- .pp_diag_validate_seed(seed)
  keep_simulations <- .pp_res_validate_logical(keep_simulations, "keep_simulations")
  keep_refits <- .pp_res_validate_logical(keep_refits, "keep_refits")
  verbose <- .pp_res_validate_logical(verbose, "verbose")
  simulation_args <- .pp_res_validate_simulation_args(simulation_args)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  observed_residuals <- residuals_pp(
    fit = fit,
    type = residual,
    grid = grid,
    nc = nc,
    mark_n = mark_n,
    chunk_size = chunk_size,
    max_mark_combinations = max_mark_combinations,
    max_cells = max_cells
  )

  observed_profiles_object <- residual_profile_pp(
    x = observed_residuals,
    variable = variable,
    form = form,
    residual = residual,
    bins = bins,
    breaks = breaks
  )
  observed_profiles <- .pp_res_profile_list(observed_profiles_object)
  profile_names <- names(observed_profiles)
  profile_breaks <- lapply(observed_profiles, function(z) z$breaks)

  simulated_curves <- lapply(
    observed_profiles,
    function(z) matrix(NA_real_, nrow = nrow(z$data), ncol = nsim)
  )
  curve_indices <- lapply(observed_profiles, function(z) z$data$index)
  simulations <- if (isTRUE(keep_simulations)) vector("list", nsim) else NULL
  refits <- if (isTRUE(keep_refits)) vector("list", nsim) else NULL

  for (i in seq_len(nsim)) {
    if (isTRUE(verbose)) {
      message(sprintf("Residual diagnostic simulation %d of %d", i, nsim))
    }

    simulation <- tryCatch(
      do.call(
        simulate_pp_fit,
        c(list(fit = fit, nsim = 1L, output = "pattern"), simulation_args)
      ),
      error = function(error) {
        stop(
          sprintf("Simulation %d failed: %s", i, conditionMessage(error)),
          call. = FALSE
        )
      }
    )

    if (nrow(simulation$coords) == 0L) {
      stop(
        paste0(
          "A simulated pattern was empty, so it cannot be refitted by the current `fit_pp()`. ",
          "This diagnostic is unsuitable for models with substantial probability of zero events."
        ),
        call. = FALSE
      )
    }

    refit_i <- tryCatch(
      refit_pp(fit = fit, pattern = simulation),
      error = function(error) {
        stop(
          sprintf("Refit %d failed: %s", i, conditionMessage(error)),
          call. = FALSE
        )
      }
    )

    residuals_i <- tryCatch(
      residuals_pp(
        fit = refit_i,
        type = residual,
        grid = observed_residuals$grid,
        mark_n = mark_n,
        chunk_size = chunk_size,
        max_mark_combinations = max_mark_combinations,
        max_cells = max_cells
      ),
      error = function(error) {
        stop(
          sprintf("Residual calculation for simulation %d failed: %s", i, conditionMessage(error)),
          call. = FALSE
        )
      }
    )

    profiles_i_object <- residual_profile_pp(
      x = residuals_i,
      variable = profile_names,
      form = form,
      residual = residual,
      bins = bins,
      breaks = profile_breaks
    )
    profiles_i <- .pp_res_profile_list(profiles_i_object)

    for (name_i in profile_names) {
      observed_index <- curve_indices[[name_i]]
      simulated_index <- profiles_i[[name_i]]$data$index

      if (length(observed_index) != length(simulated_index) ||
          !isTRUE(all.equal(observed_index, simulated_index))) {
        stop(
          sprintf(
            "Simulation %d produced a profile grid incompatible with the observed `%s` profile.",
            i,
            name_i
          ),
          call. = FALSE
        )
      }

      simulated_curves[[name_i]][, i] <- profiles_i[[name_i]]$data$residual
    }

    if (isTRUE(keep_simulations)) {
      simulations[[i]] <- simulation
    }

    if (isTRUE(keep_refits)) {
      refits[[i]] <- refit_i
    }
  }

  envelopes <- vector("list", length(profile_names))
  names(envelopes) <- profile_names

  for (name_i in profile_names) {
    profile_i <- observed_profiles[[name_i]]
    curve_set_i <- GET::curve_set(
      r = profile_i$data$index,
      obs = profile_i$data$residual,
      sim = simulated_curves[[name_i]]
    )

    envelope_i <- do.call(
      GET::global_envelope_test,
      c(
        list(
          curve_set_i,
          type = type,
          alpha = alpha,
          alternative = alternative
        ),
        list(...)
      )
    )

    envelopes[[name_i]] <- list(
      curve_set = curve_set_i,
      envelope = envelope_i,
      index = profile_i$data$index,
      observed = profile_i$data$residual,
      simulated = simulated_curves[[name_i]],
      profile = profile_i
    )
  }

  out <- list(
    call = match.call(),
    method = "plug_in_refit",
    observed_residuals = observed_residuals,
    observed_profiles = observed_profiles,
    envelopes = envelopes,
    settings = list(
      variable = profile_names,
      nsim = nsim,
      form = form,
      residual = residual,
      bins = bins,
      alpha = alpha,
      type = type,
      alternative = alternative,
      mark_n = mark_n,
      chunk_size = chunk_size,
      max_mark_combinations = max_mark_combinations,
      max_cells = max_cells,
      simulation_args = simulation_args,
      seed = seed
    ),
    simulations = simulations,
    refits = refits
  )

  class(out) <- "pp_residual_envelope"
  out
}


#' Diagnose the first-order component of a Poisson point-process model
#'
#' High-level wrapper for [envelope_residual_pp()] that labels the resulting
#' plug-in residual-envelope analysis as a first-order model diagnostic.
#'
#' @details
#' All arguments are forwarded unchanged to [envelope_residual_pp()]. The
#' returned object inherits from both `pp_first_order_diagnostic` and
#' `pp_residual_envelope`, so it supports the same print, summary, and plot
#' methods. See [envelope_residual_pp()] for residual definitions, the
#' simulation-refitting algorithm, interpretation, computational considerations,
#' and the limitation of plug-in rather than nested calibration.
#'
#' @inheritParams envelope_residual_pp
#'
#' @return An object inheriting from `pp_first_order_diagnostic` and
#'   `pp_residual_envelope`.
#'
#' @examples
#' w <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' x <- simulate_pp(
#'   window = w,
#'   expected_n = 80,
#'   output = "pattern",
#'   seed = 29
#' )
#' fit <- fit_pp(x, formula = ~ x + y, nc = 6)
#'
#' \donttest{
#' diagnostic <- diagnose_first_order_pp(
#'   fit,
#'   variable = c("x", "y"),
#'   nsim = 19,
#'   nc = 5,
#'   seed = 30,
#'   verbose = FALSE
#' )
#' print(diagnostic)
#' print(summary(diagnostic))
#' plot(diagnostic, variable = "x")
#' }
#'
#' @seealso [envelope_residual_pp()], [residuals_pp()],
#'   [residual_profile_pp()], [diagnose_pp()].
#' @family first-order point-process diagnostics
#' @export
diagnose_first_order_pp <- function(...) {
  out <- envelope_residual_pp(...)
  class(out) <- c("pp_first_order_diagnostic", "pp_residual_envelope")
  out
}


#' Methods for residual-envelope diagnostics
#'
#' Print, summarise, and plot objects returned by [envelope_residual_pp()] and
#' [diagnose_first_order_pp()].
#'
#' @details
#' Printing reports the number of simulations, profile form, residual
#' contribution, tested variables, significance level, and plug-in calibration.
#' The summary method extracts the p-value interval supplied by `GET` for each
#' variable. If only one p-value is available, it is used as both endpoints.
#'
#' The plotting method displays either the global envelope for one variable or
#' its observed residual profile. Methods for `pp_first_order_diagnostic`
#' delegate to the same implementations.
#'
#' @param x A `pp_residual_envelope` or `pp_first_order_diagnostic` object for
#'   printing and plotting, or a `summary_pp_residual_envelope` object for
#'   `print.summary_pp_residual_envelope()`.
#' @param object A residual-envelope or first-order diagnostic object to
#'   summarise.
#' @param digits Number of significant digits used when printing.
#' @param variable Optional name of the profile variable to display. It defaults
#'   to the first available variable.
#' @param what `"envelope"` or `"profile"`.
#' @param ... Additional arguments. Plotting arguments are passed to the
#'   selected `GET` plot or residual-profile plot.
#'
#' @return Printing and plotting methods invisibly return `x`. Summary methods
#'   return a `summary_pp_residual_envelope` object containing diagnostic
#'   settings and p-value intervals by variable.
#'
#' @examples
#' w <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' x <- simulate_pp(
#'   window = w,
#'   expected_n = 70,
#'   output = "pattern",
#'   seed = 31
#' )
#' fit <- fit_pp(x, formula = ~ x + y, nc = 6)
#'
#' \donttest{
#' diagnostic <- diagnose_first_order_pp(
#'   fit,
#'   variable = "x",
#'   nsim = 19,
#'   nc = 5,
#'   seed = 32,
#'   verbose = FALSE
#' )
#' print(diagnostic)
#' print(summary(diagnostic))
#' plot(diagnostic, variable = "x", what = "envelope")
#' plot(diagnostic, variable = "x", what = "profile")
#' }
#'
#' @seealso [envelope_residual_pp()], [diagnose_first_order_pp()],
#'   [residual_profile_pp()].
#' @family first-order point-process diagnostics
#' @rdname pp_residual_envelope-methods
#' @method print pp_residual_envelope
#' @export
print.pp_residual_envelope <- function(x, ...) {
  cat("<pp_residual_envelope>\n")
  cat("Plug-in refit global-envelope diagnostic for first-order residual profiles\n\n")
  cat("Primary simulations: ", x$settings$nsim, "\n", sep = "")
  cat("Profile form: ", x$settings$form, "\n", sep = "")
  cat("Residual contribution: ", x$settings$residual, "\n", sep = "")
  cat("Variables: ", paste(names(x$envelopes), collapse = ", "), "\n", sep = "")
  cat("Alpha: ", x$settings$alpha, "\n", sep = "")
  cat("Note: p-values use plug-in refit calibration, not nested composite adjustment.\n")
  invisible(x)
}


#' @rdname pp_residual_envelope-methods
#' @method summary pp_residual_envelope
#' @export
summary.pp_residual_envelope <- function(object, ...) {
  rows <- lapply(
    names(object$envelopes),
    function(name_i) {
      p_interval <- .pp_res_get_p_interval(object$envelopes[[name_i]]$envelope)

      data.frame(
        variable = name_i,
        p_lower = unname(p_interval[["lower"]]),
        p_upper = unname(p_interval[["upper"]]),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  )

  out <- list(
    method = object$method,
    nsim = object$settings$nsim,
    alpha = object$settings$alpha,
    form = object$settings$form,
    residual = object$settings$residual,
    p_values = do.call(rbind, rows)
  )

  class(out) <- "summary_pp_residual_envelope"
  out
}


#' @rdname pp_residual_envelope-methods
#' @method print summary_pp_residual_envelope
#' @export
print.summary_pp_residual_envelope <- function(x, digits = 4, ...) {
  cat("First-order residual-envelope summary\n\n")
  cat("Method: ", x$method, "\n", sep = "")
  cat("Primary simulations: ", x$nsim, "\n", sep = "")
  cat("Alpha: ", x$alpha, "\n", sep = "")
  cat("Profile form: ", x$form, "\n", sep = "")
  cat("Residual contribution: ", x$residual, "\n\n", sep = "")

  print(x$p_values, row.names = FALSE, digits = digits)
  invisible(x)
}


#' @rdname pp_residual_envelope-methods
#' @method plot pp_residual_envelope
#' @export
plot.pp_residual_envelope <- function(
    x,
    variable = NULL,
    what = c("envelope", "profile"),
    ...
) {
  if (!inherits(x, "pp_residual_envelope")) {
    stop("`x` must inherit from class 'pp_residual_envelope'.", call. = FALSE)
  }

  what <- match.arg(what)
  variables <- names(x$envelopes)

  if (is.null(variable)) {
    variable <- variables[1L]
  }

  if (!is.character(variable) || length(variable) != 1L ||
      !variable %in% variables) {
    stop(
      sprintf(
        "`variable` must name one available residual profile: %s.",
        paste(variables, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  object <- x$envelopes[[variable]]

  if (identical(what, "profile")) {
    plot(object$profile, ...)
    return(invisible(x))
  }

  plotted <- plot(object$envelope, ...)

  # GET may return a ggplot object. Inside an S3 wrapper it must be printed
  # explicitly; otherwise nothing is drawn in non-interactive method calls.
  if (inherits(plotted, "ggplot")) {
    print(plotted)
  }

  invisible(x)
}


#' @rdname pp_residual_envelope-methods
#' @method print pp_first_order_diagnostic
#' @export
print.pp_first_order_diagnostic <- function(x, ...) {
  print.pp_residual_envelope(x, ...)
}


#' @rdname pp_residual_envelope-methods
#' @method summary pp_first_order_diagnostic
#' @export
summary.pp_first_order_diagnostic <- function(object, ...) {
  summary.pp_residual_envelope(object, ...)
}


#' @rdname pp_residual_envelope-methods
#' @method plot pp_first_order_diagnostic
#' @export
plot.pp_first_order_diagnostic <- function(x, ...) {
  plot.pp_residual_envelope(x, ...)
}
