# Poisson point-process diagnostics ------------------------------------------
#
# This module builds on ppdata.R, ppfit.R, and ppsimul.R. It provides
# inhomogeneous K-functions for spatial and spatio-temporal point patterns,
# numerical marginalisation of fitted joint marked intensities, simulation and
# refitting from a pp_fit object, and plug-in global-envelope diagnostics.
#
# The implementation currently assumes hyperrectangular pp_window objects and
# uses translation edge correction.


# Validation and shared helpers ----------------------------------------------

.pp_diag_validate_positive_integer <- function(x, arg_name, min_value = 1L) {
  if (!is.numeric(x) ||
      length(x) != 1L ||
      !is.finite(x) ||
      x != as.integer(x) ||
      x < min_value) {
    stop(
      sprintf(
        "`%s` must be one finite integer greater than or equal to %d.",
        arg_name,
        min_value
      ),
      call. = FALSE
    )
  }

  as.integer(x)
}


.pp_diag_validate_probability <- function(x, arg_name) {
  if (!is.numeric(x) ||
      length(x) != 1L ||
      !is.finite(x) ||
      x <= 0 ||
      x >= 1) {
    stop(
      sprintf("`%s` must be one number strictly between zero and one.", arg_name),
      call. = FALSE
    )
  }

  as.numeric(x)
}


.pp_diag_validate_seed <- function(seed) {
  if (is.null(seed)) {
    return(NULL)
  }

  if (!is.numeric(seed) ||
      length(seed) != 1L ||
      !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("`seed` must be NULL or one finite integer.", call. = FALSE)
  }

  as.integer(seed)
}


.pp_diag_validate_pattern <- function(pattern, require_window = TRUE) {
  if (!inherits(pattern, "pp")) {
    stop("`pattern` must inherit from class 'pp'.", call. = FALSE)
  }

  if (is.null(pattern$coords) || !is.matrix(pattern$coords)) {
    stop("`pattern` must contain a numeric coordinate matrix.", call. = FALSE)
  }

  if (!is.numeric(pattern$coords) || any(!is.finite(pattern$coords))) {
    stop("`pattern$coords` must contain only finite numeric values.", call. = FALSE)
  }

  if (is.null(pattern$axis_info) || !is.data.frame(pattern$axis_info)) {
    stop("`pattern` is missing valid axis information.", call. = FALSE)
  }

  axis_names <- pattern$axis_info$name

  if (!identical(colnames(pattern$coords), axis_names)) {
    stop(
      "The coordinate matrix and axis information of `pattern` are inconsistent.",
      call. = FALSE
    )
  }

  if (isTRUE(require_window)) {
    if (is.null(pattern$window) || !inherits(pattern$window, "pp_window")) {
      stop(
        "`pattern` must carry a finite `pp_window`; diagnostics require an observation window.",
        call. = FALSE
      )
    }

    if (is.null(pattern$window$bounds) ||
        !all(axis_names %in% colnames(pattern$window$bounds))) {
      stop("The pattern window is inconsistent with its coordinate axes.", call. = FALSE)
    }
  }

  invisible(pattern)
}


.pp_diag_validate_fit <- function(fit, require_model_spec = TRUE) {
  if (!inherits(fit, "pp_fit")) {
    stop("`fit` must inherit from class 'pp_fit'.", call. = FALSE)
  }

  if (is.null(fit$gam) || !inherits(fit$gam, "gam")) {
    stop("`fit` does not contain a valid `mgcv::gam` object.", call. = FALSE)
  }

  if (is.null(fit$pattern)) {
    stop("`fit` does not retain the fitted point pattern.", call. = FALSE)
  }

  .pp_diag_validate_pattern(fit$pattern)

  if (isTRUE(require_model_spec) &&
      (is.null(fit$model_spec) || !is.list(fit$model_spec))) {
    stop(
      paste0(
        "`fit` has no `model_spec`. Refit the model with the current version of ",
        "`fit_pp()` before using prediction, simulation, refitting, or diagnostics."
      ),
      call. = FALSE
    )
  }

  invisible(fit)
}


.pp_diag_window_bounds <- function(window, axis_names) {
  bounds <- window$bounds[, axis_names, drop = FALSE]

  if (!is.numeric(bounds) ||
      nrow(bounds) != 2L ||
      any(!is.finite(bounds))) {
    stop("Window bounds must be a finite numeric two-row matrix.", call. = FALSE)
  }

  rownames(bounds) <- c("lower", "upper")

  widths <- bounds["upper", ] - bounds["lower", ]

  if (any(widths <= 0)) {
    stop("Every observation-window axis must have strictly positive width.", call. = FALSE)
  }

  bounds
}


.pp_diag_window_volume <- function(window, axis_names = window$axis_info$name) {
  bounds <- .pp_diag_window_bounds(window, axis_names)
  prod(as.numeric(bounds["upper", ] - bounds["lower", ]))
}


.pp_diag_axis_sets <- function(pattern) {
  spatial <- pattern$axis_info$name[pattern$axis_info$role == "space"]
  time <- pattern$axis_info$name[pattern$axis_info$role == "time"]

  if (length(spatial) == 0L) {
    stop("A K-function requires at least one spatial axis.", call. = FALSE)
  }

  if (length(time) > 1L) {
    stop("At most one time axis is supported.", call. = FALSE)
  }

  list(
    spatial = spatial,
    time = time,
    is_st = length(time) == 1L
  )
}


.pp_diag_unit_ball_volume <- function(d) {
  if (!is.numeric(d) || length(d) != 1L || d < 1L || d != as.integer(d)) {
    stop("`d` must be a strictly positive integer.", call. = FALSE)
  }

  pi^(d / 2) / gamma(d / 2 + 1)
}


.pp_diag_validate_lags <- function(x, arg_name) {
  if (!is.numeric(x) ||
      length(x) == 0L ||
      any(!is.finite(x)) ||
      any(x < 0)) {
    stop(
      sprintf("`%s` must be a non-empty numeric vector of non-negative finite values.", arg_name),
      call. = FALSE
    )
  }

  x <- sort(unique(as.numeric(x)))

  if (length(x) == 0L) {
    stop(sprintf("`%s` must contain at least one value.", arg_name), call. = FALSE)
  }

  x
}


.pp_diag_default_r <- function(pattern, n = 20L) {
  n <- .pp_diag_validate_positive_integer(n, "n", min_value = 2L)
  axes <- .pp_diag_axis_sets(pattern)$spatial
  bounds <- .pp_diag_window_bounds(pattern$window, axes)
  max_r <- min(bounds["upper", ] - bounds["lower", ]) / 4

  if (!is.finite(max_r) || max_r <= 0) {
    stop("Could not determine a positive default maximum spatial lag.", call. = FALSE)
  }

  seq(0, max_r, length.out = n)
}


.pp_diag_default_u <- function(pattern, n = 20L) {
  n <- .pp_diag_validate_positive_integer(n, "n", min_value = 2L)
  axes <- .pp_diag_axis_sets(pattern)

  if (!axes$is_st) {
    return(NULL)
  }

  bounds <- .pp_diag_window_bounds(pattern$window, axes$time)
  max_u <- (bounds["upper", 1L] - bounds["lower", 1L]) / 4

  if (!is.finite(max_u) || max_u <= 0) {
    stop("Could not determine a positive default temporal lag.", call. = FALSE)
  }

  seq(0, max_u, length.out = n)
}


.pp_diag_validate_lambda <- function(lambda, n, window_volume) {
  if (n == 0L) {
    if (is.null(lambda)) {
      return(numeric(0))
    }

    if (!is.numeric(lambda) || length(lambda) != 0L) {
      stop("For an empty pattern, `lambda` must be NULL or numeric(0).", call. = FALSE)
    }

    return(numeric(0))
  }

  if (is.null(lambda)) {
    lambda <- rep(n / window_volume, n)
  } else if (!is.numeric(lambda) || !(length(lambda) %in% c(1L, n))) {
    stop(
      "`lambda` must be NULL, one positive number, or a vector with one value per observed event.",
      call. = FALSE
    )
  }

  if (length(lambda) == 1L) {
    lambda <- rep(as.numeric(lambda), n)
  }

  lambda <- as.numeric(lambda)

  if (any(!is.finite(lambda)) || any(lambda <= 0)) {
    stop("Every intensity in `lambda` must be finite and strictly positive.", call. = FALSE)
  }

  lambda
}


.pp_diag_first_ge <- function(x, grid) {
  # First index k for which grid[k] >= x. Values above max(grid) become NA.
  idx <- findInterval(x, grid)
  out <- idx + 1L

  exact <- idx >= 1L
  if (any(exact)) {
    tol <- .Machine$double.eps^0.5 * pmax(1, abs(x[exact]), abs(grid[idx[exact]]))
    equal <- abs(x[exact] - grid[idx[exact]]) <= tol
    out[which(exact)[equal]] <- idx[exact][equal]
  }

  out[idx == 0L] <- 1L
  out[out > length(grid)] <- NA_integer_
  out
}


.pp_diag_cumsum2d <- function(x) {
  out <- apply(x, 2L, cumsum)
  out <- t(apply(out, 1L, cumsum))
  dimnames(out) <- dimnames(x)
  out
}


.pp_diag_model_spec <- function(fit) {
  .pp_diag_validate_fit(fit, require_model_spec = TRUE)
  spec <- fit$model_spec

  required <- c(
    "formula", "window", "axis_info", "covariate_sources",
    "continuous_marks", "categorical_levels", "mark_structure",
    "quadrature", "gam_args"
  )

  missing <- setdiff(required, names(spec))

  if (length(missing) > 0L) {
    stop(
      sprintf(
        "`fit$model_spec` is incomplete; missing component(s): %s.",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  spec
}


.pp_diag_target_coords <- function(fit, at) {
  spec <- .pp_diag_model_spec(fit)
  axis_names <- spec$axis_info$name

  if (inherits(at, "pp")) {
    .pp_diag_validate_pattern(at)

    if (!setequal(colnames(at$coords), axis_names)) {
      stop("The target pattern does not have the fitted model axes.", call. = FALSE)
    }

    out <- at$coords[, axis_names, drop = FALSE]
  } else {
    if (!is.data.frame(at) && !is.matrix(at)) {
      stop("`at` must be a pp object, data frame, or matrix.", call. = FALSE)
    }

    at <- as.data.frame(at, check.names = FALSE, stringsAsFactors = FALSE)

    missing_axes <- setdiff(axis_names, names(at))
    if (length(missing_axes) > 0L) {
      stop(
        sprintf(
          "`at` is missing fitted coordinate axis/axes: %s.",
          paste(missing_axes, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    out <- as.matrix(at[, axis_names, drop = FALSE])
  }

  storage.mode(out) <- "double"

  if (any(!is.finite(out))) {
    stop("Target coordinates must be finite numeric values.", call. = FALSE)
  }

  bounds <- .pp_diag_window_bounds(spec$window, axis_names)
  outside <- sweep(out, 2L, bounds["lower", ], "<") |
    sweep(out, 2L, bounds["upper", ], ">")

  if (any(outside)) {
    stop("Some target coordinates lie outside the fitted observation window.", call. = FALSE)
  }

  colnames(out) <- axis_names
  out
}


.pp_diag_evaluate_sources <- function(fit, coords) {
  spec <- .pp_diag_model_spec(fit)
  sources <- spec$covariate_sources

  if (length(sources) == 0L) {
    return(data.frame(row.names = seq_len(nrow(coords))))
  }

  out <- data.frame(row.names = seq_len(nrow(coords)))
  target <- as.data.frame(coords, check.names = FALSE)

  for (name_i in names(sources)) {
    source_i <- sources[[name_i]]

    if (is.null(source_i$object) || !inherits(source_i$object, "pp_covariate")) {
      stop(
        sprintf(
          "`fit$model_spec$covariate_sources$%s` lacks a valid `pp_covariate` object.",
          name_i
        ),
        call. = FALSE
      )
    }

    method_i <- source_i$method
    args_i <- source_i$args

    if (is.null(args_i)) {
      args_i <- list()
    }

    value_i <- do.call(
      evaluate_covariate,
      c(
        list(
          covariate = source_i$object,
          at = target,
          method = method_i
        ),
        args_i
      )
    )

    if (!is.numeric(value_i) ||
        length(value_i) != nrow(coords) ||
        any(!is.finite(value_i))) {
      stop(
        sprintf(
          "Evaluation of fitted covariate `%s` did not return finite numeric values.",
          name_i
        ),
        call. = FALSE
      )
    }

    out[[name_i]] <- as.numeric(value_i)
  }

  rownames(out) <- NULL
  out
}


.pp_diag_joint_mark_names <- function(fit) {
  spec <- .pp_diag_model_spec(fit)
  out <- spec$mark_structure$joint_continuous

  if (is.null(out)) {
    out <- character(0)
  }

  as.character(out)
}


.pp_diag_field_mark_names <- function(fit) {
  spec <- .pp_diag_model_spec(fit)
  out <- spec$mark_structure$field_continuous

  if (is.null(out)) {
    out <- character(0)
  }

  as.character(out)
}


.pp_diag_categorical_mark_names <- function(fit) {
  spec <- .pp_diag_model_spec(fit)
  out <- spec$mark_structure$categorical

  if (is.null(out)) {
    out <- names(spec$categorical_levels)
  }

  as.character(out)
}


.pp_diag_prepare_marks <- function(fit, marks, n) {
  spec <- .pp_diag_model_spec(fit)
  joint_names <- .pp_diag_joint_mark_names(fit)
  categorical_names <- .pp_diag_categorical_mark_names(fit)
  needed <- c(joint_names, categorical_names)

  if (length(needed) == 0L) {
    if (!is.null(marks) && ncol(as.data.frame(marks)) > 0L) {
      warning("`marks` are ignored because the fitted model has no joint marks.", call. = FALSE)
    }

    return(data.frame(row.names = seq_len(n)))
  }

  if (is.null(marks)) {
    stop(
      sprintf(
        "Joint intensity prediction requires mark(s): %s.",
        paste(needed, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!is.data.frame(marks) && !is.matrix(marks)) {
    stop("`marks` must be a data frame or matrix.", call. = FALSE)
  }

  marks <- as.data.frame(marks, check.names = FALSE, stringsAsFactors = FALSE)

  if (nrow(marks) != n) {
    stop("`marks` must have one row per target coordinate.", call. = FALSE)
  }

  missing_marks <- setdiff(needed, names(marks))
  if (length(missing_marks) > 0L) {
    stop(
      sprintf(
        "`marks` is missing required mark(s): %s.",
        paste(missing_marks, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  out <- data.frame(row.names = seq_len(n))

  for (name_i in joint_names) {
    values <- marks[[name_i]]

    if (!is.numeric(values) || any(!is.finite(values))) {
      stop(
        sprintf("Joint continuous mark `%s` must be finite and numeric.", name_i),
        call. = FALSE
      )
    }

    bounds_i <- spec$continuous_marks[[name_i]]$bounds
    tol <- .Machine$double.eps^0.5 * max(1, abs(bounds_i))

    if (any(values < bounds_i[1L] - tol | values > bounds_i[2L] + tol)) {
      stop(
        sprintf(
          "Values of joint continuous mark `%s` fall outside its fitted bounds.",
          name_i
        ),
        call. = FALSE
      )
    }

    out[[name_i]] <- as.numeric(values)
  }

  for (name_i in categorical_names) {
    values <- as.character(marks[[name_i]])
    levels_i <- as.character(spec$categorical_levels[[name_i]])

    if (anyNA(values)) {
      stop(sprintf("Categorical mark `%s` contains missing values.", name_i), call. = FALSE)
    }

    unknown <- setdiff(unique(values), levels_i)
    if (length(unknown) > 0L) {
      stop(
        sprintf(
          "Categorical mark `%s` has value(s) absent from the fitted levels: %s.",
          name_i,
          paste(unknown, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    out[[name_i]] <- factor(values, levels = levels_i)
  }

  rownames(out) <- NULL
  out
}


.pp_diag_build_newdata <- function(fit, coords, marks = NULL) {
  spec <- .pp_diag_model_spec(fit)
  coords <- as.matrix(coords)
  axis_names <- spec$axis_info$name
  n <- nrow(coords)

  out <- as.data.frame(coords[, axis_names, drop = FALSE], check.names = FALSE)
  source_values <- .pp_diag_evaluate_sources(fit, coords)

  if (ncol(source_values) > 0L) {
    out <- cbind(out, source_values)
  }

  mark_values <- .pp_diag_prepare_marks(fit, marks = marks, n = n)

  if (ncol(mark_values) > 0L) {
    out <- cbind(out, mark_values)
  }

  formula_variables <- .pp_fit_formula_variables(spec$formula)
  missing_variables <- setdiff(formula_variables, names(out))

  if (length(missing_variables) > 0L) {
    stop(
      sprintf(
        "The fitted formula cannot be evaluated because new data lack: %s.",
        paste(missing_variables, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  out
}


.pp_diag_predict_joint <- function(fit, coords, marks = NULL) {
  newdata <- .pp_diag_build_newdata(fit, coords = coords, marks = marks)

  eta <- tryCatch(
    as.numeric(stats::predict(fit$gam, newdata = newdata, type = "link")),
    error = function(error) {
      stop(
        sprintf(
          "Prediction from the fitted GAM failed: %s",
          conditionMessage(error)
        ),
        call. = FALSE
      )
    }
  )

  if (anyNA(eta) || any(is.nan(eta)) || any(eta > log(.Machine$double.xmax))) {
    stop("Fitted log-intensity could not be evaluated safely.", call. = FALSE)
  }

  out <- numeric(length(eta))
  finite <- is.finite(eta)
  out[finite] <- exp(eta[finite])
  out[is.infinite(eta) & eta < 0] <- 0

  if (any(!is.finite(out)) || any(out < 0)) {
    stop("Fitted intensity contains invalid values.", call. = FALSE)
  }

  out
}


.pp_diag_default_mark_n <- function(fit, mark_n) {
  joint_names <- .pp_diag_joint_mark_names(fit)

  if (length(joint_names) == 0L) {
    return(integer(0))
  }

  if (is.null(mark_n)) {
    nc <- fit$model_spec$quadrature$nc_resolved

    if (!is.null(nc) && !is.null(names(nc)) && all(joint_names %in% names(nc))) {
      out <- as.integer(nc[joint_names])
    } else {
      out <- rep.int(15L, length(joint_names))
    }
  } else {
    if (!is.numeric(mark_n) ||
        any(!is.finite(mark_n)) ||
        any(mark_n != as.integer(mark_n)) ||
        any(mark_n < 2L)) {
      stop(
        "`mark_n` must contain integers greater than or equal to two.",
        call. = FALSE
      )
    }

    if (!is.null(names(mark_n))) {
      if (!setequal(names(mark_n), joint_names) || anyDuplicated(names(mark_n))) {
        stop("Named `mark_n` must match the joint continuous-mark names exactly.", call. = FALSE)
      }

      out <- as.integer(mark_n[joint_names])
    } else if (length(mark_n) == 1L) {
      out <- rep.int(as.integer(mark_n), length(joint_names))
    } else if (length(mark_n) == length(joint_names)) {
      out <- as.integer(mark_n)
    } else {
      stop(
        "`mark_n` must be scalar, named by joint marks, or have one entry per joint mark.",
        call. = FALSE
      )
    }
  }

  names(out) <- joint_names
  out
}


.pp_diag_mark_integration_grid <- function(fit, mark_n = NULL, midpoint = TRUE) {
  spec <- .pp_diag_model_spec(fit)
  joint_names <- .pp_diag_joint_mark_names(fit)
  categorical_names <- .pp_diag_categorical_mark_names(fit)
  mark_n <- .pp_diag_default_mark_n(fit, mark_n)

  continuous_parts <- list()
  continuous_weights <- list()

  for (name_i in joint_names) {
    bounds_i <- as.numeric(spec$continuous_marks[[name_i]]$bounds)
    n_i <- mark_n[[name_i]]
    width_i <- diff(bounds_i)

    values_i <- if (isTRUE(midpoint)) {
      bounds_i[1L] + ((seq_len(n_i) - 0.5) / n_i) * width_i
    } else {
      seq(bounds_i[1L], bounds_i[2L], length.out = n_i)
    }

    continuous_parts[[name_i]] <- values_i
    continuous_weights[[name_i]] <- rep(width_i / n_i, n_i)
  }

  categorical_parts <- lapply(
    categorical_names,
    function(name_i) as.character(spec$categorical_levels[[name_i]])
  )
  names(categorical_parts) <- categorical_names

  parts <- c(continuous_parts, categorical_parts)

  if (length(parts) == 0L) {
    return(data.frame(.pp_mark_weight = 1, row.names = 1L))
  }

  grid <- do.call(
    expand.grid,
    c(parts, list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  )

  weight <- rep.int(1, nrow(grid))

  for (name_i in joint_names) {
    lookup <- stats::setNames(continuous_weights[[name_i]], as.character(continuous_parts[[name_i]]))
    weight <- weight * unname(lookup[as.character(grid[[name_i]])])
  }

  for (name_i in categorical_names) {
    grid[[name_i]] <- factor(
      as.character(grid[[name_i]]),
      levels = as.character(spec$categorical_levels[[name_i]])
    )
  }

  grid$.pp_mark_weight <- weight
  rownames(grid) <- NULL
  grid
}


.pp_diag_validate_mark_grid_size <- function(grid, max_mark_combinations) {
  max_mark_combinations <- .pp_diag_validate_positive_integer(
    max_mark_combinations,
    "max_mark_combinations",
    min_value = 1L
  )

  if (nrow(grid) > max_mark_combinations) {
    stop(
      sprintf(
        paste0(
          "The mark-integration grid contains %d combinations, exceeding ",
          "`max_mark_combinations = %d`. Reduce `mark_n` or increase the limit."
        ),
        nrow(grid),
        max_mark_combinations
      ),
      call. = FALSE
    )
  }

  invisible(grid)
}


.pp_diag_ground_intensity <- function(
    fit,
    coords,
    mark_n = NULL,
    chunk_size = 100000L,
    max_mark_combinations = 100000L
) {
  .pp_diag_validate_fit(fit, require_model_spec = TRUE)
  coords <- .pp_diag_target_coords(fit, coords)
  n <- nrow(coords)

  if (n == 0L) {
    return(numeric(0))
  }

  chunk_size <- .pp_diag_validate_positive_integer(
    chunk_size,
    "chunk_size",
    min_value = 1L
  )

  grid <- .pp_diag_mark_integration_grid(fit, mark_n = mark_n)
  .pp_diag_validate_mark_grid_size(grid, max_mark_combinations)

  n_grid <- nrow(grid)
  weight <- grid$.pp_mark_weight
  mark_grid <- grid[, setdiff(names(grid), ".pp_mark_weight"), drop = FALSE]
  block_n <- max(1L, floor(chunk_size / n_grid))

  out <- numeric(n)
  starts <- seq.int(1L, n, by = block_n)

  for (start in starts) {
    end <- min(n, start + block_n - 1L)
    ids <- start:end
    n_block <- length(ids)

    coords_block <- coords[ids, , drop = FALSE]
    expanded_coords <- coords_block[rep(seq_len(n_block), each = n_grid), , drop = FALSE]
    expanded_marks <- mark_grid[rep(seq_len(n_grid), times = n_block), , drop = FALSE]

    lambda_joint <- .pp_diag_predict_joint(
      fit = fit,
      coords = expanded_coords,
      marks = expanded_marks
    )

    lambda_matrix <- matrix(lambda_joint, nrow = n_block, ncol = n_grid, byrow = TRUE)
    out[ids] <- as.vector(lambda_matrix %*% weight)
  }

  if (any(!is.finite(out)) || any(out <= 0)) {
    stop(
      paste0(
        "Numerical marginalisation produced a non-positive or non-finite ground intensity. ",
        "Check mark bounds, the fitted model, and `mark_n`."
      ),
      call. = FALSE
    )
  }

  out
}


#' Predict fitted joint or ground-process intensity
#'
#' Evaluates the intensity implied by a fitted Poisson point-process model at
#' new coordinate locations. Predictions may be returned for the complete
#' joint coordinate--mark model or for the ground process obtained after
#' marginalising over all joint marks.
#'
#' @details
#' Let \eqn{\lambda(\mathbf{s},\mathbf{m},\mathbf{c})} denote the fitted joint
#' intensity at coordinates \eqn{\mathbf{s}}, continuous joint marks
#' \eqn{\mathbf{m}}, and categorical marks \eqn{\mathbf{c}}. The corresponding
#' ground-process intensity is
#'
#' \deqn{
#' \lambda_g(\mathbf{s}) =
#' \sum_{\mathbf{c}}\int_{\mathcal{M}}
#' \lambda(\mathbf{s},\mathbf{m},\mathbf{c})\,d\mathbf{m}.
#' }
#'
#' With `type = "joint"`, the function evaluates the fitted intensity at the
#' supplied coordinates and mark values. Continuous and categorical joint marks
#' required by the fitted model must be supplied through `marks`.
#'
#' With `type = "ground"`, continuous joint marks are integrated numerically
#' over the bounds retained in the fitted model using a tensor-product midpoint
#' rule, while categorical marks are summed over their complete fitted level
#' combinations. Marks treated as covariate fields are evaluated at `at` and
#' are not integrated out.
#'
#' The convenience wrapper `ground_intensity_pp()` always requests ground
#' intensity. The S3 method `predict.pp_fit()` provides the same interface
#' through [stats::predict()].
#'
#' @param object A fitted object of class `"pp_fit"` created by [fit_pp()].
#' @param fit A fitted object of class `"pp_fit"`. This is the first argument of
#'   `ground_intensity_pp()` and is equivalent to `object`.
#' @param at,newdata A `pp` object, data frame, or numeric matrix containing the
#'   fitted coordinate axes. Columns may appear in any order when they are
#'   named, but all coordinate axes retained by the fitted model must be
#'   present. Additional columns are ignored.
#' @param type Character string specifying the target intensity. `"joint"`
#'   evaluates the complete coordinate--mark intensity; `"ground"` integrates
#'   and sums over joint marks.
#' @param marks Mark values used when `type = "joint"`. Supply `NULL` for an
#'   unmarked model, or a data frame, matrix, or named list containing every
#'   continuous joint mark and categorical mark required by the fitted model.
#'   A one-row object is recycled over all requested coordinates.
#' @param mark_n Optional numerical integration resolution for continuous joint
#'   marks when `type = "ground"`. It may be one integer recycled over all joint
#'   marks, one integer per joint mark, or a named integer vector matching the
#'   joint-mark names. Every value must be at least two. By default, the
#'   resolved fitting-cubature resolution is reused when available; otherwise
#'   15 midpoint intervals are used per joint mark.
#' @param chunk_size Positive integer giving the maximum approximate number of
#'   joint-intensity evaluations processed in one block during ground-intensity
#'   marginalisation. Smaller values reduce memory use.
#' @param max_mark_combinations Positive integer limiting the Cartesian product
#'   of continuous midpoint states and categorical level combinations used for
#'   marginalisation.
#' @param ... Additional arguments required by the S3 generic. They are
#'   currently ignored.
#'
#' @return A numeric vector with one fitted intensity for each row of `at` or
#'   `newdata`. Ground-intensity predictions are strictly positive and finite;
#'   joint-intensity predictions are non-negative and finite.
#'
#' @examples
#' set.seed(1001)
#' window <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' pattern <- pp_pattern(
#'   data = data.frame(x = runif(50), y = runif(50)),
#'   coords = c("x", "y"),
#'   window = window
#' )
#' fit <- fit_pp(
#'   pattern = pattern,
#'   formula = ~ 1,
#'   nc = c(x = 6, y = 6)
#' )
#'
#' locations <- data.frame(
#'   x = c(0.25, 0.75),
#'   y = c(0.25, 0.75)
#' )
#'
#' predict_pp_fit(fit, at = locations, type = "ground")
#' ground_intensity_pp(fit, at = locations)
#' predict(fit, newdata = locations, type = "ground")
#'
#' @references
#' Daley, D. J. and Vere-Jones, D. (2003). *An Introduction to the
#' Theory of Point Processes, Volume I: Elementary Theory and Methods*,
#' second edition. Springer.
#'
#' Baddeley, A., Rubak, E., and Turner, R. (2015). *Spatial Point Patterns:
#' Methodology and Applications with R*. Chapman and Hall/CRC.
#'
#' @seealso
#' [fit_pp()] for fitting Poisson point-process models, [K_fit_pp()] for
#' fitted-intensity K-functions, and [integrate_pp_fit()] for numerical
#' integration of the ground intensity.
#'
#' @family fitted Poisson model operations
#' @export
predict_pp_fit <- function(
    object,
    at,
    type = c("ground", "joint"),
    marks = NULL,
    mark_n = NULL,
    chunk_size = 100000L,
    max_mark_combinations = 100000L,
    ...
) {
  .pp_diag_validate_fit(object, require_model_spec = TRUE)
  type <- match.arg(type)
  coords <- .pp_diag_target_coords(object, at)

  if (identical(type, "joint")) {
    return(.pp_diag_predict_joint(object, coords = coords, marks = marks))
  }

  .pp_diag_ground_intensity(
    fit = object,
    coords = coords,
    mark_n = mark_n,
    chunk_size = chunk_size,
    max_mark_combinations = max_mark_combinations
  )
}


#' @rdname predict_pp_fit
#' @export
ground_intensity_pp <- function(
    fit,
    at,
    mark_n = NULL,
    chunk_size = 100000L,
    max_mark_combinations = 100000L
) {
  predict_pp_fit(
    object = fit,
    at = at,
    type = "ground",
    mark_n = mark_n,
    chunk_size = chunk_size,
    max_mark_combinations = max_mark_combinations
  )
}


#' @rdname predict_pp_fit
#' @method predict pp_fit
#' @export
predict.pp_fit <- function(
    object,
    newdata,
    type = c("ground", "joint"),
    marks = NULL,
    mark_n = NULL,
    chunk_size = 100000L,
    max_mark_combinations = 100000L,
    ...
) {
  predict_pp_fit(
    object = object,
    at = newdata,
    type = match.arg(type),
    marks = marks,
    mark_n = mark_n,
    chunk_size = chunk_size,
    max_mark_combinations = max_mark_combinations
  )
}


# Inhomogeneous K-functions --------------------------------------------------

.pp_diag_pair_weights <- function(
    coords_i,
    coords_all,
    bounds,
    lambda_i,
    lambda_all,
    correction
) {
  n_i <- nrow(coords_i)
  n_all <- nrow(coords_all)

  denominator <- outer(lambda_i, lambda_all, "*")

  if (identical(correction, "none")) {
    return(1 / denominator)
  }

  widths <- bounds["upper", ] - bounds["lower", ]
  overlap <- matrix(1, nrow = n_i, ncol = n_all)

  for (j in seq_len(ncol(coords_i))) {
    shift <- abs(outer(coords_i[, j], coords_all[, j], "-"))
    overlap <- overlap * pmax(0, widths[j] - shift)
  }

  volume <- prod(as.numeric(widths))
  out <- matrix(0, nrow = n_i, ncol = n_all)
  valid <- overlap > 0
  out[valid] <- volume / overlap[valid] / denominator[valid]
  out
}


.pp_diag_spatial_distances <- function(coords_i, coords_all) {
  out <- matrix(0, nrow = nrow(coords_i), ncol = nrow(coords_all))

  for (j in seq_len(ncol(coords_i))) {
    delta <- outer(coords_i[, j], coords_all[, j], "-")
    out <- out + delta^2
  }

  sqrt(out)
}

.pp_diag_weighted_tabulate <- function(index, weights, nbins) {
  out <- numeric(nbins)

  if (length(index) == 0L) {
    return(out)
  }

  if (length(index) != length(weights)) {
    stop(
      "`index` and `weights` must have the same length.",
      call. = FALSE
    )
  }

  valid <- !is.na(index) &
    is.finite(index) &
    index >= 1L &
    index <= nbins &
    is.finite(weights)

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


.pp_diag_spatial_k <- function(
    coords,
    lambda,
    r,
    bounds,
    correction,
    local,
    chunk_size
) {
  n <- nrow(coords)
  nr <- length(r)
  volume <- prod(as.numeric(bounds["upper", ] - bounds["lower", ]))
  global_bins <- numeric(nr)
  local_out <- if (isTRUE(local)) {
    matrix(0, nrow = nr, ncol = n, dimnames = list(as.character(r), NULL))
  } else {
    NULL
  }

  if (n < 2L) {
    return(list(
      Khat = numeric(nr),
      Klocal = local_out
    ))
  }

  starts <- seq.int(1L, n, by = chunk_size)

  for (start in starts) {
    end <- min(n, start + chunk_size - 1L)
    ids <- start:end
    coords_i <- coords[ids, , drop = FALSE]

    distances <- .pp_diag_spatial_distances(coords_i, coords)
    weights <- .pp_diag_pair_weights(
      coords_i = coords_i,
      coords_all = coords,
      bounds = bounds,
      lambda_i = lambda[ids],
      lambda_all = lambda,
      correction = correction
    )

    weights[cbind(seq_along(ids), ids)] <- 0
    first_r <- .pp_diag_first_ge(as.vector(distances), r)
    valid <- !is.na(first_r) & as.vector(weights) > 0

    if (any(valid)) {
      bins <- .pp_diag_weighted_tabulate(
        index = first_r[valid],
        weights = as.vector(weights)[valid],
        nbins = nr
      )
      global_bins <- global_bins + bins
    }

    if (isTRUE(local)) {
      for (ii in seq_along(ids)) {
        first_i <- .pp_diag_first_ge(distances[ii, ], r)
        weights_i <- weights[ii, ]
        valid_i <- !is.na(first_i) & weights_i > 0
        bins_i <- numeric(nr)

        if (any(valid_i)) {
          bins_i <- .pp_diag_weighted_tabulate(
            index = first_i[valid_i],
            weights = weights_i[valid_i],
            nbins = nr
          )
        }

        local_out[, ids[ii]] <- (n / volume) * cumsum(bins_i)
      }
    }
  }

  list(
    Khat = cumsum(global_bins) / volume,
    Klocal = local_out
  )
}


.pp_diag_st_k <- function(
    coords_space,
    coords_time,
    lambda,
    r,
    u,
    bounds_space,
    bounds_time,
    correction,
    local,
    chunk_size
) {
  n <- nrow(coords_space)
  nr <- length(r)
  nu <- length(u)
  all_bounds <- cbind(bounds_space, bounds_time)
  volume <- prod(as.numeric(all_bounds["upper", ] - all_bounds["lower", ]))
  global_bins <- matrix(0, nrow = nr, ncol = nu)
  local_out <- if (isTRUE(local)) {
    array(
      0,
      dim = c(nr, nu, n),
      dimnames = list(as.character(r), as.character(u), NULL)
    )
  } else {
    NULL
  }

  if (n < 2L) {
    return(list(
      Khat = global_bins,
      Klocal = local_out
    ))
  }

  coords_all <- cbind(coords_space, coords_time)
  starts <- seq.int(1L, n, by = chunk_size)

  for (start in starts) {
    end <- min(n, start + chunk_size - 1L)
    ids <- start:end
    space_i <- coords_space[ids, , drop = FALSE]
    time_i <- coords_time[ids]

    spatial_distance <- .pp_diag_spatial_distances(space_i, coords_space)
    temporal_distance <- abs(outer(time_i, coords_time, "-"))

    coords_i_full <- cbind(space_i, time_i)
    weights <- .pp_diag_pair_weights(
      coords_i = coords_i_full,
      coords_all = coords_all,
      bounds = all_bounds,
      lambda_i = lambda[ids],
      lambda_all = lambda,
      correction = correction
    )

    weights[cbind(seq_along(ids), ids)] <- 0
    first_r <- .pp_diag_first_ge(as.vector(spatial_distance), r)
    first_u <- .pp_diag_first_ge(as.vector(temporal_distance), u)
    valid <- !is.na(first_r) & !is.na(first_u) & as.vector(weights) > 0

    if (any(valid)) {
      bin_index <- first_r[valid] + (first_u[valid] - 1L) * nr
      bins <- .pp_diag_weighted_tabulate(
        index = bin_index,
        weights = as.vector(weights)[valid],
        nbins = nr * nu
      )
      global_bins <- global_bins + matrix(bins, nrow = nr, ncol = nu)
    }

    if (isTRUE(local)) {
      for (ii in seq_along(ids)) {
        first_r_i <- .pp_diag_first_ge(spatial_distance[ii, ], r)
        first_u_i <- .pp_diag_first_ge(temporal_distance[ii, ], u)
        weights_i <- weights[ii, ]
        valid_i <- !is.na(first_r_i) & !is.na(first_u_i) & weights_i > 0
        bins_i <- matrix(0, nrow = nr, ncol = nu)

        if (any(valid_i)) {
          bin_index_i <- first_r_i[valid_i] + (first_u_i[valid_i] - 1L) * nr
          bins_i <- matrix(
            .pp_diag_weighted_tabulate(
              index = bin_index_i,
              weights = weights_i[valid_i],
              nbins = nr * nu
            ),
            nrow = nr,
            ncol = nu
          )
        }

        local_out[, , ids[ii]] <- (n / volume) * .pp_diag_cumsum2d(bins_i)
      }
    }
  }

  list(
    Khat = .pp_diag_cumsum2d(global_bins) / volume,
    Klocal = local_out
  )
}


#' Estimate an inhomogeneous K-function
#'
#' Estimates an intensity-reweighted K-function for a spatial point pattern in
#' any number of spatial dimensions or an inhomogeneous space-time K-function
#' for a pattern containing one temporal axis.
#'
#' @details
#' For a \eqn{d}-dimensional spatial point pattern observed in a
#' hyperrectangular window \eqn{W}, the implemented estimator has the form
#'
#' \deqn{
#' \widehat K_{\mathrm{inhom}}(r) =
#' \frac{1}{|W|}
#' \sum_i\sum_{j\ne i}
#' \frac{
#' \mathbf{1}\{\|\mathbf{s}_i-\mathbf{s}_j\|\le r\}
#' e_{ij}}
#' {\lambda_i\lambda_j},
#' }
#'
#' where \eqn{\lambda_i} is the supplied ground-process intensity at event
#' \eqn{i}. With translation correction,
#' \eqn{e_{ij}=|W|/|W\cap(W+\mathbf{s}_i-\mathbf{s}_j)|}; for
#' `correction = "none"`, \eqn{e_{ij}=1}. Because `pp_window` objects are
#' hyperrectangles, the translated overlap is calculated exactly as the product
#' of the overlapping side lengths.
#'
#' Under a homogeneous Poisson reference process, the theoretical spatial
#' curve is
#'
#' \deqn{K_0(r)=v_d r^d,\qquad
#' v_d=\frac{\pi^{d/2}}{\Gamma(d/2+1)}.}
#'
#' For a spatio-temporal pattern, pairs are counted when their spatial distance
#' is at most \eqn{r} and their absolute temporal separation is at most
#' \eqn{u}. The Poisson reference surface is
#'
#' \deqn{K_0(r,u)=2v_d r^d u.}
#'
#' The estimator always uses the intensity of the *ground process*. Therefore,
#' when a fitted model contains joint marks, use [K_fit_pp()] so that the joint
#' intensity is first integrated over continuous marks and summed over
#' categorical levels.
#'
#' If `local = TRUE`, eventwise contributions are retained. They are scaled so
#' that their arithmetic mean at each lag equals the corresponding global
#' estimate. Since pairwise distances are evaluated in chunks, `chunk_size`
#' controls memory use without changing the estimator.
#'
#' @param pattern A point-pattern object created by [pp_pattern()] and carrying
#'   a finite hyperrectangular [pp_window()].
#' @param lambda `NULL` or a strictly positive finite numeric vector with one
#'   ground-process intensity value per observed event. If `NULL`, the
#'   homogeneous estimate \eqn{n/|W|} is used at every event.
#' @param r Optional numeric vector of non-negative spatial lags. Values are
#'   sorted and duplicated values are removed. By default, 20 lags are placed
#'   between zero and one quarter of the shortest spatial side of the window.
#' @param u Optional numeric vector of non-negative temporal lags, used only for
#'   a spatio-temporal pattern. Values are sorted and duplicates removed. By
#'   default, 20 lags are placed between zero and one quarter of the temporal
#'   window length. Supplying `u` for a purely spatial pattern is an error.
#' @param correction Edge correction. `"translate"` uses exact translation
#'   correction for the hyperrectangular window; `"none"` applies no edge
#'   correction.
#' @param local Logical value indicating whether eventwise local K contributions
#'   should be retained.
#' @param chunk_size Positive integer giving the number of focal events whose
#'   pairwise distances and weights are processed together. Lower values reduce
#'   peak memory use.
#'
#' @return An object of class `"pp_K"` and either `"pp_K_spatial"` or
#'   `"pp_K_st"`. Its main components are:
#'
#' - `Khat`: estimated K curve or space-time K surface;
#' - `Ktheo`: homogeneous Poisson reference;
#' - `residual`: `Khat - Ktheo`;
#' - `r` and, for space-time patterns, `u`: evaluated lags;
#' - `Klocal`: local contributions when requested, otherwise `NULL`;
#' - `lambda`: eventwise ground intensities used in the estimator;
#' - `correction`, `n_points`, and axis/window metadata.
#'
#' @examples
#' set.seed(1002)
#' window <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' pattern <- pp_pattern(
#'   data = data.frame(x = runif(80), y = runif(80)),
#'   coords = c("x", "y"),
#'   window = window
#' )
#'
#' k <- K_pp(
#'   pattern,
#'   r = seq(0, 0.25, length.out = 12)
#' )
#' k
#' summary(k)
#' plot(k, what = "residual")
#'
#' # Retain contributions associated with each observed event
#' k_local <- K_pp(
#'   pattern,
#'   r = seq(0, 0.25, length.out = 12),
#'   local = TRUE
#' )
#' plot(k_local, what = "local", local_id = 1)
#'
#' # A small spatio-temporal example
#' window_st <- pp_window(
#'   lower = c(x = 0, y = 0, time = 0),
#'   upper = c(x = 1, y = 1, time = 10),
#'   axis_roles = c(x = "space", y = "space", time = "time")
#' )
#' pattern_st <- pp_pattern(
#'   data = data.frame(
#'     x = runif(40),
#'     y = runif(40),
#'     time = runif(40, 0, 10)
#'   ),
#'   coords = c("x", "y", "time"),
#'   axis_roles = c("space", "space", "time"),
#'   window = window_st
#' )
#' k_st <- K_pp(
#'   pattern_st,
#'   r = seq(0, 0.25, length.out = 6),
#'   u = seq(0, 2.5, length.out = 6)
#' )
#'
#' @references
#' Ripley, B. D. (1976). The second-order analysis of stationary point
#' processes. *Journal of Applied Probability*, 13, 255--266.
#'
#' Baddeley, A. J., Møller, J., and Waagepetersen, R. (2000).
#' Non- and semi-parametric estimation of interaction in inhomogeneous point
#' patterns. *Statistica Neerlandica*, 54, 329--350.
#'
#' Gabriel, E. and Diggle, P. J. (2009). Second-order analysis of
#' inhomogeneous spatio-temporal point process data. *Statistica Neerlandica*,
#' 63, 43--51.
#'
#' @seealso
#' [Klocal_pp()] for the local wrapper, [K_fit_pp()] for K-functions based on a
#' fitted intensity, and [envelope_pp()] for simulation-based global envelopes.
#'
#' @family second-order diagnostics
#' @export
K_pp <- function(
    pattern,
    lambda = NULL,
    r = NULL,
    u = NULL,
    correction = c("translate", "none"),
    local = FALSE,
    chunk_size = 2000L
) {
  .pp_diag_validate_pattern(pattern, require_window = TRUE)
  correction <- match.arg(correction)
  chunk_size <- .pp_diag_validate_positive_integer(
    chunk_size,
    "chunk_size",
    min_value = 1L
  )

  if (!is.logical(local) || length(local) != 1L || is.na(local)) {
    stop("`local` must be TRUE or FALSE.", call. = FALSE)
  }

  axes <- .pp_diag_axis_sets(pattern)
  n <- nrow(pattern$coords)
  all_axes <- pattern$axis_info$name
  volume <- .pp_diag_window_volume(pattern$window, all_axes)
  lambda <- .pp_diag_validate_lambda(lambda, n = n, window_volume = volume)

  if (is.null(r)) {
    r <- .pp_diag_default_r(pattern)
  }
  r <- .pp_diag_validate_lags(r, "r")

  if (!axes$is_st) {
    if (!is.null(u)) {
      stop("`u` can be supplied only for spatio-temporal point patterns.", call. = FALSE)
    }

    bounds_space <- .pp_diag_window_bounds(pattern$window, axes$spatial)
    result <- .pp_diag_spatial_k(
      coords = pattern$coords[, axes$spatial, drop = FALSE],
      lambda = lambda,
      r = r,
      bounds = bounds_space,
      correction = correction,
      local = local,
      chunk_size = chunk_size
    )

    d <- length(axes$spatial)
    Ktheo <- .pp_diag_unit_ball_volume(d) * r^d

    out <- list(
      call = match.call(),
      type = "spatial",
      Khat = result$Khat,
      Ktheo = Ktheo,
      residual = result$Khat - Ktheo,
      r = r,
      u = NULL,
      Klocal = result$Klocal,
      lambda = lambda,
      correction = correction,
      local = local,
      n_points = n,
      spatial_dimension = d,
      spatial_axes = axes$spatial,
      time_axis = NULL,
      window = pattern$window,
      pattern = pattern
    )
  } else {
    if (is.null(u)) {
      u <- .pp_diag_default_u(pattern)
    }
    u <- .pp_diag_validate_lags(u, "u")

    bounds_space <- .pp_diag_window_bounds(pattern$window, axes$spatial)
    bounds_time <- .pp_diag_window_bounds(pattern$window, axes$time)

    result <- .pp_diag_st_k(
      coords_space = pattern$coords[, axes$spatial, drop = FALSE],
      coords_time = pattern$coords[, axes$time],
      lambda = lambda,
      r = r,
      u = u,
      bounds_space = bounds_space,
      bounds_time = bounds_time,
      correction = correction,
      local = local,
      chunk_size = chunk_size
    )

    d <- length(axes$spatial)
    Ktheo <- outer(
      .pp_diag_unit_ball_volume(d) * r^d,
      u,
      function(a, b) 2 * a * b
    )
    dimnames(Ktheo) <- list(as.character(r), as.character(u))
    dimnames(result$Khat) <- list(as.character(r), as.character(u))

    out <- list(
      call = match.call(),
      type = "spatio_temporal",
      Khat = result$Khat,
      Ktheo = Ktheo,
      residual = result$Khat - Ktheo,
      r = r,
      u = u,
      Klocal = result$Klocal,
      lambda = lambda,
      correction = correction,
      local = local,
      n_points = n,
      spatial_dimension = d,
      spatial_axes = axes$spatial,
      time_axis = axes$time,
      window = pattern$window,
      pattern = pattern
    )
  }

  class(out) <- c("pp_K", if (axes$is_st) "pp_K_st" else "pp_K_spatial")
  out
}


#' Calculate local inhomogeneous K contributions
#'
#' Convenience wrapper for `K_pp(..., local = TRUE)`. It returns the global
#' inhomogeneous K estimate together with the contribution associated with each
#' observed event.
#'
#' @param pattern A point-pattern object created by [pp_pattern()] and carrying
#'   a finite hyperrectangular [pp_window()].
#' @param lambda `NULL` or a strictly positive finite numeric vector with one
#'   ground-process intensity value per observed event. If `NULL`, the
#'   homogeneous estimate \eqn{n/|W|} is used.
#' @param r Optional vector of non-negative spatial lags. See [K_pp()].
#' @param u Optional vector of non-negative temporal lags for a
#'   spatio-temporal pattern. See [K_pp()].
#' @param correction Edge correction passed to [K_pp()]: `"translate"` or
#'   `"none"`.
#' @param chunk_size Positive integer giving the number of focal events
#'   processed together.
#'
#' @return An object of class `"pp_K"` with a non-`NULL` `Klocal` component.
#' For a spatial pattern, `Klocal` is a lag-by-event matrix. For a
#' spatio-temporal pattern, it is a spatial-lag by temporal-lag by event array.
#' At every lag, the arithmetic mean over events equals the global `Khat`.
#'
#' @examples
#' set.seed(1003)
#' window <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' pattern <- pp_pattern(
#'   data = data.frame(x = runif(40), y = runif(40)),
#'   coords = c("x", "y"),
#'   window = window
#' )
#'
#' local_k <- Klocal_pp(
#'   pattern,
#'   r = seq(0, 0.25, length.out = 8)
#' )
#' dim(local_k$Klocal)
#' plot(local_k, what = "local", local_id = 1)
#'
#' @seealso
#' [K_pp()] for the underlying estimator and [Klocal_fit_pp()] for local
#' contributions using a fitted ground intensity.
#'
#' @family second-order diagnostics
#' @export
Klocal_pp <- function(
    pattern,
    lambda = NULL,
    r = NULL,
    u = NULL,
    correction = c("translate", "none"),
    chunk_size = 2000L
) {
  K_pp(
    pattern = pattern,
    lambda = lambda,
    r = r,
    u = u,
    correction = match.arg(correction),
    local = TRUE,
    chunk_size = chunk_size
  )
}


#' Calculate an inhomogeneous K-function from a fitted Poisson model
#'
#' Evaluates the fitted ground-process intensity at the events of a point
#' pattern and uses it in [K_pp()]. The companion `Klocal_fit_pp()` fixes
#' `local = TRUE`.
#'
#' @details
#' For unmarked models and models containing only marks treated as covariate
#' fields, the fitted intensity can be evaluated directly at the event
#' coordinates. For joint marked models, the ground intensity is obtained from
#'
#' \deqn{
#' \lambda_g(\mathbf{s}) =
#' \sum_{\mathbf{c}}\int_{\mathcal{M}}
#' \lambda(\mathbf{s},\mathbf{m},\mathbf{c})\,d\mathbf{m},
#' }
#'
#' using the midpoint marginalisation implemented by [predict_pp_fit()].
#' Consequently, the resulting K-function diagnoses second-order structure in
#' the ground process rather than dependence in the mark distribution.
#'
#' When `pattern` differs from the pattern used for fitting, it must retain the
#' same coordinate axes and exactly the same observation-window bounds. This is
#' the route used internally after simulation and refitting in [envelope_pp()].
#'
#' @param fit A fitted object of class `"pp_fit"` created by [fit_pp()].
#' @param pattern A `pp` object on which to evaluate the statistic. It defaults
#'   to the observed pattern stored in `fit` and must have the same coordinate
#'   axes and observation-window bounds.
#' @param r Optional vector of non-negative spatial lags. See [K_pp()].
#' @param u Optional vector of non-negative temporal lags for a
#'   spatio-temporal pattern. See [K_pp()].
#' @param correction Edge correction passed to [K_pp()]: `"translate"` or
#'   `"none"`.
#' @param local Logical value indicating whether local eventwise contributions
#'   should be retained. `Klocal_fit_pp()` always uses `TRUE`.
#' @param chunk_size Positive integer controlling both the block size used for
#'   pairwise K calculations and the minimum prediction block size.
#' @param mark_n Optional midpoint integration resolution for continuous joint
#'   marks. See [predict_pp_fit()].
#' @param max_mark_combinations Positive integer limiting the continuous and
#'   categorical mark-state combinations used when calculating ground
#'   intensity.
#'
#' @return An object of class `"pp_K"` as described in [K_pp()]. The `lambda`
#'   component contains fitted ground intensities evaluated at the events in
#'   `pattern`.
#'
#' @examples
#' set.seed(1004)
#' window <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' pattern <- pp_pattern(
#'   data = data.frame(x = runif(60), y = runif(60)),
#'   coords = c("x", "y"),
#'   window = window
#' )
#' fit <- fit_pp(
#'   pattern = pattern,
#'   formula = ~ x + y,
#'   nc = c(x = 7, y = 7)
#' )
#'
#' fitted_k <- K_fit_pp(
#'   fit,
#'   r = seq(0, 0.25, length.out = 10)
#' )
#' plot(fitted_k, what = "residual")
#'
#' fitted_local_k <- Klocal_fit_pp(
#'   fit,
#'   r = seq(0, 0.25, length.out = 10)
#' )
#' dim(fitted_local_k$Klocal)
#'
#' @references
#' Baddeley, A. J., Møller, J., and Waagepetersen, R. (2000).
#' Non- and semi-parametric estimation of interaction in inhomogeneous point
#' patterns. *Statistica Neerlandica*, 54, 329--350.
#'
#' Gabriel, E. and Diggle, P. J. (2009). Second-order analysis of
#' inhomogeneous spatio-temporal point process data. *Statistica Neerlandica*,
#' 63, 43--51.
#'
#' @seealso
#' [predict_pp_fit()] for the fitted ground intensity, [K_pp()] for direct
#' calculation from supplied intensities, and [envelope_pp()] for global
#' simulation envelopes.
#'
#' @family second-order diagnostics
#' @export
K_fit_pp <- function(
    fit,
    pattern = fit$pattern,
    r = NULL,
    u = NULL,
    correction = c("translate", "none"),
    local = FALSE,
    chunk_size = 2000L,
    mark_n = NULL,
    max_mark_combinations = 100000L
) {
  .pp_diag_validate_fit(fit, require_model_spec = TRUE)
  .pp_diag_validate_pattern(pattern, require_window = TRUE)

  if (!isTRUE(all.equal(pattern$window$bounds, fit$model_spec$window$bounds)) ||
      !identical(pattern$axis_info$name, fit$model_spec$axis_info$name)) {
    stop(
      "`pattern` must use the same observation window and coordinate axes as `fit`.",
      call. = FALSE
    )
  }

  lambda <- predict_pp_fit(
    object = fit,
    at = pattern,
    type = "ground",
    mark_n = mark_n,
    chunk_size = max(chunk_size, 1L),
    max_mark_combinations = max_mark_combinations
  )

  K_pp(
    pattern = pattern,
    lambda = lambda,
    r = r,
    u = u,
    correction = match.arg(correction),
    local = local,
    chunk_size = chunk_size
  )
}


#' @rdname K_fit_pp
#' @export
Klocal_fit_pp <- function(
    fit,
    pattern = fit$pattern,
    r = NULL,
    u = NULL,
    correction = c("translate", "none"),
    chunk_size = 2000L,
    mark_n = NULL,
    max_mark_combinations = 100000L
) {
  K_fit_pp(
    fit = fit,
    pattern = pattern,
    r = r,
    u = u,
    correction = match.arg(correction),
    local = TRUE,
    chunk_size = chunk_size,
    mark_n = mark_n,
    max_mark_combinations = max_mark_combinations
  )
}


#' Methods for inhomogeneous K-function objects
#'
#' Print, summarise, and plot objects returned by [K_pp()] or [K_fit_pp()].
#'
#' @param x,object An object inheriting from `"pp_K"`. For
#'   `print.summary_pp_K()`, `x` is an object returned by `summary.pp_K()`.
#' @param digits Positive integer controlling the number of significant digits
#'   printed for numerical summaries.
#' @param what Component to display with `plot.pp_K()`: `"residual"` plots
#'   `Khat - Ktheo`; `"K"` plots the estimate and, for spatial objects, its
#'   Poisson reference; `"local"` plots the contribution of one event and
#'   requires that the object was computed with `local = TRUE`.
#' @param local_id Integer identifying the observed event to display when
#'   `what = "local"`.
#' @param ... Additional graphical arguments passed to [graphics::plot()] for a
#'   spatial curve or [graphics::image()] for a spatio-temporal surface. They
#'   are ignored by the print and summary methods.
#'
#' @details
#' Spatial objects are plotted as one-dimensional curves over `r`.
#' Spatio-temporal objects are displayed as images over the `(r, u)` lag grid.
#' The summary reports the maximum absolute departure of the estimated K
#' function from its homogeneous Poisson reference over the evaluated grid.
#'
#' @return
#' `print.pp_K()`, `print.summary_pp_K()`, and `plot.pp_K()` return their input
#' invisibly. `summary.pp_K()` returns an object of class
#' `"summary_pp_K"` containing dimensions, correction, lag counts, and the
#' maximum absolute K residual.
#'
#' @examples
#' set.seed(1005)
#' window <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' pattern <- pp_pattern(
#'   data = data.frame(x = runif(50), y = runif(50)),
#'   coords = c("x", "y"),
#'   window = window
#' )
#' k <- K_pp(
#'   pattern,
#'   r = seq(0, 0.25, length.out = 10),
#'   local = TRUE
#' )
#'
#' print(k)
#' summary(k)
#' plot(k, what = "K")
#' plot(k, what = "local", local_id = 1)
#'
#' @seealso
#' [K_pp()], [K_fit_pp()], and [envelope_pp()].
#'
#' @family second-order diagnostics
#' @rdname pp_K-methods
#' @method print pp_K
#' @export
print.pp_K <- function(x, digits = 4, ...) {
  cat("<pp_K>\n")
  cat(
    if (identical(x$type, "spatial")) {
      "Inhomogeneous spatial K-function\n\n"
    } else {
      "Inhomogeneous spatio-temporal K-function\n\n"
    }
  )
  cat("Points: ", x$n_points, "\n", sep = "")
  cat("Spatial dimension: ", x$spatial_dimension, "\n", sep = "")
  cat("Correction: ", x$correction, "\n", sep = "")
  cat("Local contributions: ", if (isTRUE(x$local)) "yes" else "no", "\n", sep = "")
  cat("Spatial lags: ", length(x$r), "\n", sep = "")

  if (!is.null(x$u)) {
    cat("Temporal lags: ", length(x$u), "\n", sep = "")
  }

  invisible(x)
}


#' @rdname pp_K-methods
#' @method summary pp_K
#' @export
summary.pp_K <- function(object, ...) {
  max_abs_residual <- max(abs(object$residual), na.rm = TRUE)

  out <- list(
    type = object$type,
    n_points = object$n_points,
    spatial_dimension = object$spatial_dimension,
    correction = object$correction,
    n_r = length(object$r),
    n_u = if (is.null(object$u)) 0L else length(object$u),
    max_abs_residual = max_abs_residual,
    local = object$local
  )

  class(out) <- "summary_pp_K"
  out
}


#' @rdname pp_K-methods
#' @method print summary_pp_K
#' @export
print.summary_pp_K <- function(x, digits = 4, ...) {
  cat("Inhomogeneous K-function summary\n\n")
  cat("Type: ", x$type, "\n", sep = "")
  cat("Points: ", x$n_points, "\n", sep = "")
  cat("Spatial dimension: ", x$spatial_dimension, "\n", sep = "")
  cat("Correction: ", x$correction, "\n", sep = "")
  cat("Spatial lags: ", x$n_r, "\n", sep = "")

  if (x$n_u > 0L) {
    cat("Temporal lags: ", x$n_u, "\n", sep = "")
  }

  cat(
    "Maximum absolute departure from Poisson reference: ",
    formatC(x$max_abs_residual, format = "fg", digits = digits),
    "\n",
    sep = ""
  )

  invisible(x)
}


#' @rdname pp_K-methods
#' @method plot pp_K
#' @export
plot.pp_K <- function(
    x,
    what = c("residual", "K", "local"),
    local_id = NULL,
    ...
) {
  what <- match.arg(what)

  if (identical(x$type, "spatial")) {
    y <- switch(
      what,
      residual = x$residual,
      K = x$Khat,
      local = {
        if (is.null(x$Klocal)) {
          stop("Local contributions were not requested; use `local = TRUE`.", call. = FALSE)
        }
        if (is.null(local_id) ||
            !is.numeric(local_id) ||
            length(local_id) != 1L ||
            local_id < 1L ||
            local_id > ncol(x$Klocal)) {
          stop("`local_id` must identify one observed point.", call. = FALSE)
        }
        x$Klocal[, as.integer(local_id)]
      }
    )

    graphics::plot(
      x$r,
      y,
      type = "l",
      xlab = "Spatial lag r",
      ylab = if (identical(what, "residual")) "Khat(r) - K0(r)" else "K(r)",
      ...
    )

    if (identical(what, "K")) {
      graphics::lines(x$r, x$Ktheo, lty = 2)
      graphics::legend(
        "topleft",
        legend = c("Estimated", "Poisson reference"),
        lty = c(1, 2),
        bty = "n"
      )
    } else if (identical(what, "residual")) {
      graphics::abline(h = 0, lty = 2)
    }

    return(invisible(x))
  }

  z <- switch(
    what,
    residual = x$residual,
    K = x$Khat,
    local = {
      if (is.null(x$Klocal)) {
        stop("Local contributions were not requested; use `local = TRUE`.", call. = FALSE)
      }
      if (is.null(local_id) ||
          !is.numeric(local_id) ||
          length(local_id) != 1L ||
          local_id < 1L ||
          local_id > dim(x$Klocal)[3L]) {
        stop("`local_id` must identify one observed point.", call. = FALSE)
      }
      x$Klocal[, , as.integer(local_id)]
    }
  )

  graphics::image(
    x = x$r,
    y = x$u,
    z = z,
    xlab = "Spatial lag r",
    ylab = "Temporal lag u",
    main = if (identical(what, "residual")) {
      "Khat(r, u) - K0(r, u)"
    } else {
      "Inhomogeneous K(r, u)"
    },
    ...
  )

  invisible(x)
}


# Simulation and refitting ---------------------------------------------------

.pp_diag_joint_domain <- function(fit) {
  spec <- .pp_diag_model_spec(fit)
  axis_names <- spec$axis_info$name
  bounds_coords <- .pp_diag_window_bounds(spec$window, axis_names)
  joint_names <- .pp_diag_joint_mark_names(fit)
  categorical_names <- .pp_diag_categorical_mark_names(fit)

  bounds <- bounds_coords

  if (length(joint_names) > 0L) {
    mark_bounds <- vapply(
      joint_names,
      function(name_i) as.numeric(spec$continuous_marks[[name_i]]$bounds),
      numeric(2L)
    )
    mark_bounds <- as.matrix(mark_bounds)
    colnames(mark_bounds) <- joint_names
    rownames(mark_bounds) <- c("lower", "upper")
    bounds <- cbind(bounds, mark_bounds)
  }

  categorical_grid <- if (length(categorical_names) == 0L) {
    data.frame(.pp_category_id = 1L)
  } else {
    parts <- lapply(
      categorical_names,
      function(name_i) as.character(spec$categorical_levels[[name_i]])
    )
    names(parts) <- categorical_names
    do.call(
      expand.grid,
      c(parts, list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
    )
  }

  for (name_i in categorical_names) {
    categorical_grid[[name_i]] <- factor(
      as.character(categorical_grid[[name_i]]),
      levels = as.character(spec$categorical_levels[[name_i]])
    )
  }

  list(
    bounds = bounds,
    coordinate_axes = axis_names,
    joint_mark_names = joint_names,
    categorical_mark_names = categorical_names,
    categorical_grid = categorical_grid,
    continuous_volume = prod(as.numeric(bounds["upper", ] - bounds["lower", ])),
    n_categories = nrow(categorical_grid),
    augmented_volume = prod(as.numeric(bounds["upper", ] - bounds["lower", ])) * nrow(categorical_grid)
  )
}


.pp_diag_resolve_bound_n <- function(bounds, n_categories, bound_n, bound_max_points) {
  d <- ncol(bounds)
  axis_names <- colnames(bounds)
  bound_max_points <- .pp_diag_validate_positive_integer(
    bound_max_points,
    "bound_max_points",
    min_value = 2L
  )

  max_continuous <- floor(bound_max_points / n_categories)

  if (max_continuous < 2L^d) {
    stop(
      paste0(
        "`bound_max_points` is too small for at least two grid points on every ",
        "continuous joint-domain axis and every categorical combination."
      ),
      call. = FALSE
    )
  }

  if (is.null(bound_n)) {
    per_axis <- max(2L, floor(max_continuous^(1 / d)))
    out <- rep.int(as.integer(per_axis), d)
    names(out) <- axis_names

    while (prod(as.double(out)) > max_continuous) {
      j <- which.max(out)
      out[j] <- out[j] - 1L
    }

    return(out)
  }

  if (!is.numeric(bound_n) ||
      any(!is.finite(bound_n)) ||
      any(bound_n != as.integer(bound_n)) ||
      any(bound_n < 2L)) {
    stop("`bound_n` must contain integers greater than or equal to two.", call. = FALSE)
  }

  if (!is.null(names(bound_n))) {
    if (!setequal(names(bound_n), axis_names) || anyDuplicated(names(bound_n))) {
      stop("Named `bound_n` must match all joint-domain axis names exactly.", call. = FALSE)
    }
    out <- as.integer(bound_n[axis_names])
  } else if (length(bound_n) == 1L) {
    out <- rep.int(as.integer(bound_n), d)
  } else if (length(bound_n) == d) {
    out <- as.integer(bound_n)
  } else {
    stop(
      "`bound_n` must be scalar, named by joint-domain axes, or have one entry per axis.",
      call. = FALSE
    )
  }

  names(out) <- axis_names

  if (prod(as.double(out)) * n_categories > bound_max_points) {
    stop(
      sprintf(
        "The requested joint-domain bound grid has %.0f points, exceeding `bound_max_points = %.0f`.",
        prod(as.double(out)) * n_categories,
        bound_max_points
      ),
      call. = FALSE
    )
  }

  out
}


.pp_diag_joint_bound_grid <- function(fit, bound_n, bound_max_points, bound_safety) {
  domain <- .pp_diag_joint_domain(fit)
  resolution <- .pp_diag_resolve_bound_n(
    bounds = domain$bounds,
    n_categories = domain$n_categories,
    bound_n = bound_n,
    bound_max_points = bound_max_points
  )

  sequences <- lapply(
    seq_along(resolution),
    function(j) {
      seq(
        domain$bounds["lower", j],
        domain$bounds["upper", j],
        length.out = resolution[j]
      )
    }
  )
  names(sequences) <- names(resolution)

  continuous_grid <- do.call(
    expand.grid,
    c(sequences, list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  )

  coordinate_grid <- continuous_grid[, domain$coordinate_axes, drop = FALSE]
  mark_grid <- continuous_grid[, domain$joint_mark_names, drop = FALSE]

  if (domain$n_categories > 1L ||
      length(domain$categorical_mark_names) > 0L) {
    coordinate_grid <- coordinate_grid[
      rep(seq_len(nrow(coordinate_grid)), each = domain$n_categories),
      ,
      drop = FALSE
    ]
    mark_grid <- mark_grid[
      rep(seq_len(nrow(mark_grid)), each = domain$n_categories),
      ,
      drop = FALSE
    ]
    categorical_grid <- domain$categorical_grid[
      rep(seq_len(domain$n_categories), times = nrow(continuous_grid)),
      ,
      drop = FALSE
    ]
  } else {
    categorical_grid <- data.frame(row.names = seq_len(nrow(continuous_grid)))
  }

  marks <- cbind(mark_grid, categorical_grid)
  marks <- marks[, setdiff(names(marks), ".pp_category_id"), drop = FALSE]

  lambda <- .pp_diag_predict_joint(
    fit = fit,
    coords = as.matrix(coordinate_grid),
    marks = marks
  )

  max_lambda <- max(lambda)

  if (!is.finite(max_lambda) || max_lambda < 0) {
    stop("Could not obtain a finite joint-intensity bound from the fitted model.", call. = FALSE)
  }

  list(
    lambda_upper = max_lambda * bound_safety,
    max_lambda_on_grid = max_lambda,
    resolution = resolution,
    n_grid = length(lambda)
  )
}


.pp_diag_sample_joint_candidates <- function(domain, n) {
  bounds <- domain$bounds
  coord_names <- domain$coordinate_axes
  joint_names <- domain$joint_mark_names
  categorical_names <- domain$categorical_mark_names

  continuous <- matrix(
    NA_real_,
    nrow = n,
    ncol = ncol(bounds),
    dimnames = list(NULL, colnames(bounds))
  )

  if (n > 0L) {
    for (j in seq_len(ncol(bounds))) {
      continuous[, j] <- stats::runif(
        n,
        min = bounds["lower", j],
        max = bounds["upper", j]
      )
    }
  }

  coords <- as.data.frame(continuous[, coord_names, drop = FALSE], check.names = FALSE)
  marks <- as.data.frame(continuous[, joint_names, drop = FALSE], check.names = FALSE)

  if (length(categorical_names) > 0L) {
    category_id <- sample.int(domain$n_categories, size = n, replace = TRUE)
    cats <- domain$categorical_grid[category_id, categorical_names, drop = FALSE]
    rownames(cats) <- NULL
    marks <- cbind(marks, cats)
  }

  list(coords = coords, marks = marks)
}


.pp_diag_build_simulated_pattern <- function(fit, coords, marks, simulation_metadata) {
  spec <- .pp_diag_model_spec(fit)
  coords <- as.data.frame(coords, check.names = FALSE)
  rownames(coords) <- NULL
  n <- nrow(coords)

  source_values <- .pp_diag_evaluate_sources(fit, as.matrix(coords))
  field_names <- .pp_diag_field_mark_names(fit)

  if (length(field_names) > 0L) {
    missing_fields <- setdiff(field_names, names(source_values))
    if (length(missing_fields) > 0L) {
      stop(
        sprintf(
          "Fitted field mark(s) cannot be reconstructed: %s.",
          paste(missing_fields, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    field_marks <- source_values[, field_names, drop = FALSE]
    marks <- cbind(marks, field_marks)
  }

  if (is.null(marks)) {
    marks <- data.frame(row.names = seq_len(n))
  }

  marks <- as.data.frame(marks, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(marks) <- NULL

  categorical_names <- .pp_diag_categorical_mark_names(fit)
  for (name_i in categorical_names) {
    if (name_i %in% names(marks)) {
      marks[[name_i]] <- factor(
        as.character(marks[[name_i]]),
        levels = as.character(spec$categorical_levels[[name_i]])
      )
    }
  }

  data <- cbind(coords, marks)
  mark_names <- names(marks)

  pattern <- pp_pattern(
    data = data,
    coords = spec$axis_info$name,
    marks = if (length(mark_names) == 0L) NULL else mark_names,
    window = spec$window
  )

  # Covariate values are intentionally not attached to `pattern$covariates`.
  # `refit_pp()` passes all retained external and pattern-derived fields as
  # explicit pp_covariate objects; attaching duplicate per-event columns here
  # would make `fit_pp()` detect two sources with the same name. For audit
  # purposes, the values are retained only in the simulation metadata.
  simulation_metadata$evaluated_covariates <- source_values[, setdiff(names(source_values), field_names), drop = FALSE]

  pattern$metadata$simulation <- simulation_metadata
  pattern
}


#' Simulate from a fitted Poisson point-process model
#'
#' Generates independent point patterns from a fitted `"pp_fit"` model by
#' thinning a homogeneous Poisson process on the fitted model's complete joint
#' coordinate--mark domain.
#'
#' @details
#' For a joint intensity bounded by \eqn{\lambda_{\max}} on an augmented domain
#' \eqn{\mathcal{D}}, the algorithm first generates
#'
#' \deqn{N^\star\sim\mathrm{Poisson}(\lambda_{\max}|\mathcal{D}|)}
#'
#' uniform candidate states and accepts a candidate \eqn{z} with probability
#'
#' \deqn{\lambda(z)/\lambda_{\max}.}
#'
#' The augmented domain contains all coordinate axes, continuous marks declared
#' with `mode = "joint"`, and every fitted combination of categorical-mark
#' levels. Its reference measure is therefore the product of continuous volume
#' and the counting measure over categorical combinations. Continuous marks
#' treated as fields are reconstructed deterministically from their retained
#' covariate source after thinning.
#'
#' With `bound = "grid"`, the joint intensity is evaluated over a regular
#' tensor-product grid and multiplied by `bound_safety`. A numerical grid does
#' not mathematically guarantee a global upper bound. If a subsequently sampled
#' candidate exceeds it, `on_bound_violation = "restart"` discards the complete
#' thinning attempt, increases the bound, and restarts. This preserves exact
#' thinning conditional on eventually obtaining a valid bound. Alternatively,
#' `"error"` stops immediately.
#'
#' A user-supplied `lambda_upper` must bound the *joint* intensity, not the
#' ground intensity. A loose bound remains valid but can make simulation
#' inefficient.
#'
#' @param fit A fitted object of class `"pp_fit"` created by [fit_pp()]. It must
#'   retain a current `model_spec`.
#' @param nsim Positive integer giving the number of independent point patterns
#'   to generate.
#' @param bound Method used to obtain the thinning upper bound. `"grid"`
#'   estimates it numerically; `"user"` uses `lambda_upper`.
#' @param lambda_upper A verified finite non-negative upper bound for the fitted
#'   joint intensity. It is required and used only when `bound = "user"`.
#' @param bound_n Optional grid resolution used when `bound = "grid"`. It may be
#'   one integer recycled over all continuous joint-domain axes, one integer per
#'   axis, or a named vector matching all coordinate and joint continuous-mark
#'   axes. Every value must be at least two. If `NULL`, a balanced resolution is
#'   selected subject to `bound_max_points`.
#' @param bound_max_points Positive integer limiting the total number of bound
#'   grid evaluations after multiplying the continuous grid size by the number
#'   of categorical-mark combinations.
#' @param bound_safety Finite multiplier greater than or equal to one applied to
#'   the maximum intensity found on a numerical bound grid.
#' @param on_bound_violation Action taken if a candidate intensity exceeds the
#'   current upper bound: `"restart"` increases the bound and repeats the entire
#'   thinning attempt; `"error"` stops.
#' @param restart_safety Finite multiplier strictly greater than one used to
#'   enlarge the bound after a violation.
#' @param max_restarts Positive integer limiting full thinning restarts within
#'   each simulation.
#' @param max_candidates Positive integer limiting both the expected and realised
#'   number of candidates in one thinning attempt.
#' @param output Output format. `"pattern"` returns simulated `pp` objects.
#'   `"diagnostic"` returns objects of class `"pp_fit_simulation"` containing
#'   the pattern together with candidate counts, acceptance rate, restarts, and
#'   bound information.
#' @param seed `NULL` or one finite integer used to initialise R's random-number
#'   generator before all simulations.
#'
#' @return If `nsim = 1`, a single `pp` object for `output = "pattern"` or a
#'   `"pp_fit_simulation"` object for `output = "diagnostic"`. If `nsim > 1`, a
#'   named list `sim_1`, ..., `sim_nsim` containing the requested output type.
#'   Simulated patterns retain fitting-compatible axes, window, marks, and
#'   simulation metadata.
#'
#' @examples
#' set.seed(1006)
#' window <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' pattern <- pp_pattern(
#'   data = data.frame(x = runif(50), y = runif(50)),
#'   coords = c("x", "y"),
#'   window = window
#' )
#' fit <- fit_pp(
#'   pattern = pattern,
#'   formula = ~ 1,
#'   nc = c(x = 6, y = 6)
#' )
#'
#' simulated <- simulate_pp_fit(
#'   fit,
#'   bound_n = c(x = 8, y = 8),
#'   seed = 1007
#' )
#' simulated
#'
#' simulation_info <- simulate_pp_fit(
#'   fit,
#'   bound_n = 8,
#'   output = "diagnostic",
#'   seed = 1007
#' )
#' simulation_info$acceptance_rate
#'
#' simulations <- simulate_pp_fit(
#'   fit,
#'   nsim = 2,
#'   bound_n = 8,
#'   seed = 1008
#' )
#' names(simulations)
#'
#' @references
#' Lewis, P. A. W. and Shedler, G. S. (1979). Simulation of nonhomogeneous
#' Poisson processes by thinning. *Naval Research Logistics Quarterly*, 26,
#' 403--413.
#'
#' Daley, D. J. and Vere-Jones, D. (2003). *An Introduction to the
#' Theory of Point Processes, Volume I: Elementary Theory and Methods*,
#' second edition. Springer.
#'
#' @seealso
#' [simulate_pp()] for direct simulation from a user-specified intensity,
#' [refit_pp()] for refitting the retained model specification, and
#' [envelope_pp()] for repeated simulation and refitting.
#'
#' @family fitted Poisson model operations
#' @export
simulate_pp_fit <- function(
    fit,
    nsim = 1L,
    bound = c("grid", "user"),
    lambda_upper = NULL,
    bound_n = NULL,
    bound_max_points = 50000L,
    bound_safety = 1.20,
    on_bound_violation = c("restart", "error"),
    restart_safety = 1.20,
    max_restarts = 10L,
    max_candidates = 1e6,
    output = c("pattern", "diagnostic"),
    seed = NULL
) {
  .pp_diag_validate_fit(fit, require_model_spec = TRUE)
  nsim <- .pp_diag_validate_positive_integer(nsim, "nsim", min_value = 1L)
  bound <- match.arg(bound)
  on_bound_violation <- match.arg(on_bound_violation)
  output <- match.arg(output)
  seed <- .pp_diag_validate_seed(seed)
  bound_max_points <- .pp_diag_validate_positive_integer(
    bound_max_points,
    "bound_max_points",
    min_value = 2L
  )
  max_restarts <- .pp_diag_validate_positive_integer(
    max_restarts,
    "max_restarts",
    min_value = 1L
  )
  max_candidates <- .pp_diag_validate_positive_integer(
    max_candidates,
    "max_candidates",
    min_value = 1L
  )

  if (!is.numeric(bound_safety) || length(bound_safety) != 1L ||
      !is.finite(bound_safety) || bound_safety < 1) {
    stop("`bound_safety` must be one finite number greater than or equal to one.", call. = FALSE)
  }

  if (!is.numeric(restart_safety) || length(restart_safety) != 1L ||
      !is.finite(restart_safety) || restart_safety <= 1) {
    stop("`restart_safety` must be one finite number strictly greater than one.", call. = FALSE)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  domain <- .pp_diag_joint_domain(fit)

  if (identical(bound, "user")) {
    if (!is.numeric(lambda_upper) ||
        length(lambda_upper) != 1L ||
        !is.finite(lambda_upper) ||
        lambda_upper < 0) {
      stop(
        "`lambda_upper` must be one finite non-negative number when `bound = 'user'`.",
        call. = FALSE
      )
    }

    initial_bound <- as.numeric(lambda_upper)
    bound_details <- list(
      method = "user",
      lambda_upper_initial = initial_bound,
      max_lambda_on_grid = NA_real_,
      resolution = NULL,
      n_grid = NULL
    )
  } else {
    grid_bound <- .pp_diag_joint_bound_grid(
      fit = fit,
      bound_n = bound_n,
      bound_max_points = bound_max_points,
      bound_safety = bound_safety
    )

    initial_bound <- grid_bound$lambda_upper
    bound_details <- c(list(method = "grid", lambda_upper_initial = initial_bound), grid_bound[-1L])
  }

  one_simulation <- function() {
    current_bound <- initial_bound
    restarts <- 0L

    repeat {
      expected_candidates <- current_bound * domain$augmented_volume

      if (!is.finite(expected_candidates) || expected_candidates > max_candidates) {
        stop(
          sprintf(
            paste0(
              "The joint thinning process would require an expected %.0f candidates, ",
              "exceeding `max_candidates = %.0f`."
            ),
            expected_candidates,
            max_candidates
          ),
          call. = FALSE
        )
      }

      n_candidates <- stats::rpois(1L, lambda = expected_candidates)

      if (n_candidates > max_candidates) {
        stop(
          sprintf(
            "The joint thinning process generated %d candidates, exceeding `max_candidates = %d`.",
            n_candidates,
            max_candidates
          ),
          call. = FALSE
        )
      }

      candidates <- .pp_diag_sample_joint_candidates(domain, n_candidates)
      lambda_candidate <- .pp_diag_predict_joint(
        fit = fit,
        coords = as.matrix(candidates$coords),
        marks = candidates$marks
      )

      tolerance <- .Machine$double.eps^0.5 * max(1, current_bound)
      violation <- lambda_candidate > current_bound + tolerance

      if (any(violation)) {
        observed_max <- max(lambda_candidate[violation])

        if (identical(on_bound_violation, "error")) {
          stop(
            sprintf(
              paste0(
                "The joint thinning bound was violated: %.8g exceeds lambda_upper = %.8g. ",
                "Supply a verified bound or use `on_bound_violation = 'restart'`."
              ),
              observed_max,
              current_bound
            ),
            call. = FALSE
          )
        }

        restarts <- restarts + 1L
        if (restarts > max_restarts) {
          stop(
            "The joint thinning bound was violated more often than `max_restarts`.",
            call. = FALSE
          )
        }

        current_bound <- max(
          current_bound * restart_safety,
          observed_max * restart_safety
        )
        next
      }

      accepted <- if (current_bound == 0) {
        rep(FALSE, n_candidates)
      } else {
        stats::runif(n_candidates) <= lambda_candidate / current_bound
      }

      metadata <- list(
        mode = "joint_thinning_from_pp_fit",
        lambda_upper = current_bound,
        lambda_upper_initial = initial_bound,
        n_candidates = n_candidates,
        n_accepted = sum(accepted),
        acceptance_rate = if (n_candidates == 0L) 0 else mean(accepted),
        thinning_restarts = restarts,
        bound = bound_details
      )

      pattern <- .pp_diag_build_simulated_pattern(
        fit = fit,
        coords = candidates$coords[accepted, , drop = FALSE],
        marks = candidates$marks[accepted, , drop = FALSE],
        simulation_metadata = metadata
      )

      if (identical(output, "pattern")) {
        return(pattern)
      }

      out <- c(metadata, list(pattern = pattern))
      class(out) <- "pp_fit_simulation"
      return(out)
    }
  }

  results <- lapply(seq_len(nsim), function(i) one_simulation())

  if (nsim == 1L) {
    return(results[[1L]])
  }

  names(results) <- paste0("sim_", seq_len(nsim))
  results
}


#' Refit a Poisson point-process model to a new pattern
#'
#' Reconstructs the public model specification retained by a `"pp_fit"` object
#' and applies the same first-order model to another compatible point pattern.
#' This function is intended primarily for parametric bootstrap and Monte Carlo
#' diagnostics.
#'
#' @details
#' The refit reuses:
#'
#' - the one-sided model formula and `mgcv::gam()` arguments;
#' - external covariate objects and their evaluation settings;
#' - continuous joint-mark bounds and continuous mark-field specifications;
#' - categorical marks and their complete fitted levels;
#' - cubature layout, resolved resolution, and computational limits.
#'
#' The new pattern must have exactly the same coordinate-axis names and
#' observation-window bounds as the original fit. It must also contain all
#' marks needed by the retained model. Covariate values are re-evaluated from
#' the retained sources rather than copied from the original observed events.
#'
#' Event-specific weights are deliberately not transferred to the new pattern.
#' Thus, if the original model was fitted with fractional `event_weights`, the
#' refit uses unit weights for the new events. This is usually appropriate for
#' simulated catalogues but should be considered when interpreting bootstrap
#' diagnostics.
#'
#' @param fit A fitted object of class `"pp_fit"` produced by [fit_pp()].
#' @param pattern A new `pp` object with exactly the same coordinate axes and
#'   observation-window bounds as `fit`, and with all marks required by the
#'   retained formula.
#' @param seed `NULL` or one finite integer. It affects the reconstructed
#'   cubature only when the retained dummy layout is `"random"`.
#' @param keep_cubature Optional storage level passed to [fit_pp()]:
#'   `"summary"`, `"full"`, or `"none"`. If `NULL`, the value retained by the
#'   original fit is reused.
#'
#' @return A new object of class `"pp_fit"` fitted to `pattern`.
#'
#' @examples
#' set.seed(1009)
#' window <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' pattern <- pp_pattern(
#'   data = data.frame(x = runif(50), y = runif(50)),
#'   coords = c("x", "y"),
#'   window = window
#' )
#' fit <- fit_pp(
#'   pattern = pattern,
#'   formula = ~ x + y,
#'   nc = c(x = 6, y = 6)
#' )
#'
#' simulated <- simulate_pp_fit(
#'   fit,
#'   bound_n = 8,
#'   seed = 1010
#' )
#' refitted <- refit_pp(fit, pattern = simulated)
#' summary(refitted)
#'
#' @seealso
#' [fit_pp()] for the original fitting interface, [simulate_pp_fit()] for
#' compatible simulations, and [envelope_pp()] for automated simulation,
#' refitting, and K-function comparison.
#'
#' @family fitted Poisson model operations
#' @export
refit_pp <- function(
    fit,
    pattern,
    seed = NULL,
    keep_cubature = NULL
) {
  .pp_diag_validate_fit(fit, require_model_spec = TRUE)
  .pp_diag_validate_pattern(pattern, require_window = TRUE)

  seed <- .pp_diag_validate_seed(seed)
  spec <- .pp_diag_model_spec(fit)

  if (!identical(pattern$axis_info$name, spec$axis_info$name) ||
      !isTRUE(all.equal(pattern$window$bounds, spec$window$bounds))) {
    stop(
      "`pattern` must have exactly the same coordinate axes and window bounds as `fit`.",
      call. = FALSE
    )
  }

  sources <- spec$covariate_sources

  external_names <- names(sources)[
    vapply(
      sources,
      function(z) !identical(z$source, "continuous_mark_field"),
      logical(1)
    )
  ]

  covariates <- lapply(
    sources[external_names],
    function(z) z$object
  )

  methods <- vapply(
    sources[external_names],
    function(z) z$method,
    character(1)
  )

  args <- lapply(
    sources[external_names],
    function(z) z$args
  )

  if (length(covariates) == 0L) {
    covariates <- NULL
    methods <- NULL
    args <- NULL
  } else {
    names(covariates) <- external_names
    names(methods) <- external_names
    names(args) <- external_names
  }

  nc <- spec$quadrature$nc_resolved

  if (is.null(nc)) {
    nc <- spec$quadrature$nc_requested
  }

  if (is.null(keep_cubature)) {
    keep_cubature <- spec$keep_cubature
  }

  categorical_levels_refit <- spec$categorical_levels
  categorical_marks_refit <- names(categorical_levels_refit)

  if (is.null(categorical_levels_refit) ||
      length(categorical_levels_refit) == 0L) {
    categorical_levels_refit <- NULL
    categorical_marks_refit <- NULL
  }

  continuous_marks_refit <- spec$continuous_marks

  if (is.null(continuous_marks_refit) ||
      length(continuous_marks_refit) == 0L) {
    continuous_marks_refit <- NULL
  }

  fit_pp(
    pattern = pattern,
    formula = spec$formula,
    covariates = covariates,
    covariate_methods = methods,
    covariate_args = args,
    continuous_marks = continuous_marks_refit,
    categorical_marks = categorical_marks_refit,
    categorical_levels = categorical_levels_refit,
    q = spec$quadrature$q,
    dummy_layout = spec$quadrature$dummy_layout,
    nc = nc,
    max_tiles = spec$quadrature$max_tiles,
    max_quadrature_points = spec$quadrature$max_quadrature_points,
    max_mark_combinations = spec$quadrature$max_mark_combinations,
    seed = seed,
    gam_args = spec$gam_args,
    keep_cubature = keep_cubature
  )
}


# Global-envelope diagnostics ------------------------------------------------

.pp_diag_curve_set <- function(K_observed, K_simulated) {
  if (identical(K_observed$type, "spatial")) {
    obs <- as.numeric(K_observed$Khat)
    sim <- vapply(K_simulated, function(z) as.numeric(z$Khat), numeric(length(obs)))
    r_index <- K_observed$r
    grid_index <- data.frame(r = K_observed$r)
  } else {
    obs <- as.numeric(K_observed$Khat)
    sim <- vapply(K_simulated, function(z) as.numeric(z$Khat), numeric(length(obs)))
    r_index <- seq_along(obs)
    grid_index <- expand.grid(
      r = K_observed$r,
      u = K_observed$u,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  }

  list(
    obs = obs,
    sim = sim,
    r_index = r_index,
    grid_index = grid_index
  )
}


#' Global-envelope diagnostic for a fitted Poisson point process
#'
#' Performs a plug-in parametric-bootstrap diagnostic based on the
#' intensity-reweighted ground-process K-function and a global envelope test
#' from the `GET` package.
#'
#' @details
#' The procedure consists of four steps:
#'
#' 1. calculate [K_fit_pp()] for the observed pattern;
#' 2. simulate `nsim` independent patterns from the fitted model with
#'    [simulate_pp_fit()];
#' 3. refit the same model specification to every simulated pattern with
#'    [refit_pp()] and recalculate its K-function;
#' 4. combine the observed and simulated curves in [GET::curve_set()] and call
#'    [GET::global_envelope_test()].
#'
#' Refitting in every primary simulation accounts for the most immediate effect
#' of parameter estimation. The implementation nevertheless remains a
#' *plug-in* test: it does not perform the second simulation level required for
#' a fully adjusted Monte Carlo test of a composite hypothesis. The reported
#' p-value should therefore be interpreted as an approximate calibration under
#' the fitted model.
#'
#' For a spatial pattern, the K curves are passed to `GET` on their natural
#' `r` grid. For a spatio-temporal pattern, each K surface is flattened in a
#' fixed `(r, u)` order; `curve_index` retains the mapping back to the original
#' lag grid. The plotting method displays the GET envelope directly for spatial
#' curves. For a spatio-temporal object, it displays the observed residual
#' surface because the GET object is stored on the flattened index.
#'
#' Values such as `nsim = 199` may be sufficient for exploratory work, whereas
#' 499 or more simulations provide finer Monte Carlo resolution for final
#' analyses. The required number depends on the selected global-envelope
#' ordering and inferential purpose.
#'
#' @param fit A fitted object of class `"pp_fit"` created by [fit_pp()].
#' @param nsim Integer greater than or equal to two giving the number of primary
#'   simulations and refits.
#' @param r Optional non-negative spatial lags used for both observed and
#'   simulated K-functions. Defaults are determined by [K_pp()].
#' @param u Optional non-negative temporal lags for a spatio-temporal fit.
#'   Defaults are determined by [K_pp()].
#' @param correction Edge correction passed to [K_fit_pp()]: `"translate"` or
#'   `"none"`.
#' @param alpha Number strictly between zero and one giving the significance
#'   level of the global envelope.
#' @param type Character string selecting an ordering supported by
#'   [GET::global_envelope_test()]. The default `"erl"` uses extreme rank
#'   length.
#' @param alternative Alternative hypothesis passed to
#'   [GET::global_envelope_test()]: `"two.sided"`, `"less"`, or `"greater"`.
#' @param mark_n Optional midpoint integration resolution used to obtain ground
#'   intensity from continuous joint-mark models. See [predict_pp_fit()].
#' @param chunk_size Positive integer controlling pairwise K-function block
#'   calculations.
#' @param max_mark_combinations Positive integer limiting the mark-state grid
#'   used during fitted ground-intensity marginalisation.
#' @param simulation_args Named list of additional arguments forwarded to
#'   [simulate_pp_fit()] for every primary simulation. It may contain options
#'   such as `bound`, `lambda_upper`, `bound_n`, or bound-safety controls, but
#'   cannot contain `fit`, `nsim`, `seed`, or `output`.
#' @param keep_simulations Logical value indicating whether all simulated
#'   `pp` objects should be retained in the returned object. The default saves
#'   memory.
#' @param keep_refits Logical value indicating whether all simulated-data
#'   `"pp_fit"` objects should be retained. This can require substantial memory.
#' @param seed `NULL` or one finite integer used to initialise the complete
#'   simulation sequence.
#' @param verbose Logical value indicating whether progress messages should be
#'   emitted for the primary simulations.
#' @param ... Additional named arguments passed to
#'   [GET::global_envelope_test()].
#'
#' @return An object of class `"pp_envelope"` with components:
#'
#' - `observed`: the observed `"pp_K"` object;
#' - `simulated`: a list of K-functions from simulated and refitted patterns;
#' - `curve_set`: the `"curve_set"` supplied to `GET`;
#' - `envelope`: the result returned by [GET::global_envelope_test()];
#' - `curve_index`: the correspondence between flattened curve positions and
#'   spatial or spatio-temporal lags;
#' - `settings`: simulation, K-function, and envelope settings;
#' - `simulations` and `refits`: retained objects when requested.
#'
#' @examples
#' \donttest{
#' set.seed(1011)
#' window <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' pattern <- pp_pattern(
#'   data = data.frame(x = runif(60), y = runif(60)),
#'   coords = c("x", "y"),
#'   window = window
#' )
#' fit <- fit_pp(
#'   pattern = pattern,
#'   formula = ~ x + y,
#'   nc = c(x = 6, y = 6)
#' )
#'
#' # nsim = 9 is used only to keep this example short.
#' # Use substantially more simulations for scientific inference.
#' envelope <- envelope_pp(
#'   fit,
#'   nsim = 9,
#'   alpha = 0.10,
#'   r = seq(0, 0.25, length.out = 10),
#'   simulation_args = list(bound_n = 8),
#'   seed = 1012,
#'   verbose = FALSE
#' )
#' envelope
#' summary(envelope)
#' plot(envelope)
#' }
#'
#' @references
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H., and Hahn, U. (2017).
#' Global envelope tests for spatial processes. *Journal of the Royal
#' Statistical Society: Series B*, 79, 381--404.
#'
#' Baddeley, A. J., Møller, J., and Waagepetersen, R. (2000).
#' Non- and semi-parametric estimation of interaction in inhomogeneous point
#' patterns. *Statistica Neerlandica*, 54, 329--350.
#'
#' Gabriel, E. and Diggle, P. J. (2009). Second-order analysis of
#' inhomogeneous spatio-temporal point process data. *Statistica Neerlandica*,
#' 63, 43--51.
#'
#' @seealso
#' [diagnose_pp()] for the high-level class alias, [K_fit_pp()] for the observed
#' statistic, [simulate_pp_fit()] for fitted-model simulation, and
#' [diagnose_first_order_pp()] for residual-based first-order diagnostics.
#'
#' @family second-order diagnostics
#' @export
envelope_pp <- function(
    fit,
    nsim = 199L,
    r = NULL,
    u = NULL,
    correction = c("translate", "none"),
    alpha = 0.05,
    type = "erl",
    alternative = c("two.sided", "less", "greater"),
    mark_n = NULL,
    chunk_size = 2000L,
    max_mark_combinations = 100000L,
    simulation_args = list(),
    keep_simulations = FALSE,
    keep_refits = FALSE,
    seed = NULL,
    verbose = TRUE,
    ...
) {
  .pp_diag_validate_fit(fit, require_model_spec = TRUE)

  if (!requireNamespace("GET", quietly = TRUE)) {
    stop(
      "Package 'GET' is required for global-envelope diagnostics; add it to Imports and install it.",
      call. = FALSE
    )
  }

  nsim <- .pp_diag_validate_positive_integer(nsim, "nsim", min_value = 2L)
  alpha <- .pp_diag_validate_probability(alpha, "alpha")
  correction <- match.arg(correction)
  alternative <- match.arg(alternative)
  seed <- .pp_diag_validate_seed(seed)

  if (!is.list(simulation_args) ||
      (length(simulation_args) > 0L && (is.null(names(simulation_args)) || any(!nzchar(names(simulation_args)))))) {
    stop("`simulation_args` must be a named list or an empty list.", call. = FALSE)
  }

  prohibited <- intersect(
    names(simulation_args),
    c("fit", "nsim", "seed", "output")
  )

  if (length(prohibited) > 0L) {
    stop(
      sprintf(
        "`simulation_args` cannot contain: %s.",
        paste(prohibited, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!is.logical(keep_simulations) || length(keep_simulations) != 1L || is.na(keep_simulations)) {
    stop("`keep_simulations` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(keep_refits) || length(keep_refits) != 1L || is.na(keep_refits)) {
    stop("`keep_refits` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  observed <- K_fit_pp(
    fit = fit,
    r = r,
    u = u,
    correction = correction,
    local = FALSE,
    chunk_size = chunk_size,
    mark_n = mark_n,
    max_mark_combinations = max_mark_combinations
  )

  simulations <- if (isTRUE(keep_simulations)) vector("list", nsim) else NULL
  refits <- if (isTRUE(keep_refits)) vector("list", nsim) else NULL
  K_simulated <- vector("list", nsim)

  for (i in seq_len(nsim)) {
    if (isTRUE(verbose)) {
      message(sprintf("Diagnostic simulation %d of %d", i, nsim))
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
          "A simulated pattern was empty, so the model cannot be refitted by the current `fit_pp()`. ",
          "This diagnostic is not suitable for models with substantial probability of zero events."
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

    K_i <- tryCatch(
      K_fit_pp(
        fit = refit_i,
        r = observed$r,
        u = observed$u,
        correction = correction,
        local = FALSE,
        chunk_size = chunk_size,
        mark_n = mark_n,
        max_mark_combinations = max_mark_combinations
      ),
      error = function(error) {
        stop(
          sprintf("K-function calculation for simulation %d failed: %s", i, conditionMessage(error)),
          call. = FALSE
        )
      }
    )

    K_simulated[[i]] <- K_i

    if (isTRUE(keep_simulations)) {
      simulations[[i]] <- simulation
    }

    if (isTRUE(keep_refits)) {
      refits[[i]] <- refit_i
    }
  }

  curves <- .pp_diag_curve_set(observed, K_simulated)
  curve_set <- GET::curve_set(
    r = curves$r_index,
    obs = curves$obs,
    sim = curves$sim
  )

  envelope <- do.call(
    GET::global_envelope_test,
    c(
      list(
        curve_set,
        type = type,
        alpha = alpha,
        alternative = alternative
      ),
      list(...)
    )
  )

  out <- list(
    call = match.call(),
    method = "plug_in_refit",
    observed = observed,
    simulated = K_simulated,
    curve_set = curve_set,
    envelope = envelope,
    curve_index = curves$grid_index,
    settings = list(
      nsim = nsim,
      alpha = alpha,
      type = type,
      alternative = alternative,
      correction = correction,
      mark_n = mark_n,
      chunk_size = chunk_size,
      max_mark_combinations = max_mark_combinations,
      simulation_args = simulation_args,
      seed = seed
    ),
    simulations = simulations,
    refits = refits
  )

  class(out) <- "pp_envelope"
  out
}


#' Diagnose a fitted Poisson point-process model
#'
#' High-level wrapper around [envelope_pp()]. It performs the same plug-in
#' simulation, refitting, K-function, and global-envelope procedure but adds the
#' more descriptive class `"pp_diagnostic"` to the result.
#'
#' @param ... Arguments passed unchanged to [envelope_pp()], including `fit`,
#'   `nsim`, lag grids, simulation controls, and `GET` options.
#'
#' @return An object inheriting from both `"pp_diagnostic"` and
#' `"pp_envelope"`.
#'
#' @examples
#' \donttest{
#' set.seed(1013)
#' window <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' pattern <- pp_pattern(
#'   data = data.frame(x = runif(50), y = runif(50)),
#'   coords = c("x", "y"),
#'   window = window
#' )
#' fit <- fit_pp(
#'   pattern = pattern,
#'   formula = ~ 1,
#'   nc = c(x = 6, y = 6)
#' )
#'
#' diagnostic <- diagnose_pp(
#'   fit = fit,
#'   nsim = 9,
#'   alpha = 0.10,
#'   r = seq(0, 0.25, length.out = 8),
#'   simulation_args = list(bound_n = 8),
#'   seed = 1014,
#'   verbose = FALSE
#' )
#' diagnostic
#' }
#'
#' @seealso
#' [envelope_pp()] for the complete method description and
#' [diagnose_first_order_pp()] for residual-based first-order diagnostics.
#'
#' @family second-order diagnostics
#' @export
diagnose_pp <- function(...) {
  out <- envelope_pp(...)
  class(out) <- c("pp_diagnostic", "pp_envelope")
  out
}


#' Methods for global-envelope diagnostic objects
#'
#' Print, summarise, and plot objects returned by [envelope_pp()] and
#' [diagnose_pp()].
#'
#' @param x,object An object inheriting from `"pp_envelope"` or
#'   `"pp_diagnostic"`. For `print.summary_pp_envelope()`, `x` is the result of
#'   `summary.pp_envelope()`.
#' @param what Component displayed by `plot.pp_envelope()`. `"envelope"` plots
#'   the `GET` global envelope for a spatial K curve. `"residual"` plots the
#'   observed `Khat - Ktheo`. For a spatio-temporal object, `"envelope"` also
#'   displays the observed residual surface because the GET curves are stored
#'   on a flattened index.
#' @param ... Additional arguments passed to the underlying `GET` plot method
#'   or to [plot.pp_K()]. They are ignored by the print and summary methods.
#'
#' @return
#' The print and plot methods return their input invisibly.
#' `summary.pp_envelope()` returns an object of class
#' `"summary_pp_envelope"` containing the diagnostic settings and, when
#' available, the summary of the stored `GET` envelope. Methods for
#' `"pp_diagnostic"` delegate to the corresponding `"pp_envelope"` methods.
#'
#' @examples
#' \donttest{
#' set.seed(1015)
#' window <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' pattern <- pp_pattern(
#'   data = data.frame(x = runif(50), y = runif(50)),
#'   coords = c("x", "y"),
#'   window = window
#' )
#' fit <- fit_pp(
#'   pattern = pattern,
#'   formula = ~ 1,
#'   nc = c(x = 6, y = 6)
#' )
#' envelope <- envelope_pp(
#'   fit,
#'   nsim = 9,
#'   alpha = 0.10,
#'   r = seq(0, 0.25, length.out = 8),
#'   simulation_args = list(bound_n = 8),
#'   seed = 1016,
#'   verbose = FALSE
#' )
#'
#' print(envelope)
#' summary(envelope)
#' plot(envelope, what = "envelope")
#' plot(envelope, what = "residual")
#' }
#'
#' @seealso
#' [envelope_pp()], [diagnose_pp()], and [plot.pp_K()].
#'
#' @family second-order diagnostics
#' @rdname pp_envelope-methods
#' @method print pp_envelope
#' @export
print.pp_envelope <- function(x, ...) {
  cat("<pp_envelope>\n")
  cat("Plug-in refit global-envelope diagnostic\n\n")
  cat("Primary simulations: ", x$settings$nsim, "\n", sep = "")
  cat("Envelope ordering: ", x$settings$type, "\n", sep = "")
  cat("Alternative: ", x$settings$alternative, "\n", sep = "")
  cat("Alpha: ", x$settings$alpha, "\n", sep = "")
  cat("Correction: ", x$settings$correction, "\n", sep = "")
  cat("K type: ", x$observed$type, "\n", sep = "")
  cat("Note: p-values use plug-in refit calibration, not double-simulation adjustment.\n")
  invisible(x)
}


#' @rdname pp_envelope-methods
#' @method summary pp_envelope
#' @export
summary.pp_envelope <- function(object, ...) {
  envelope_summary <- tryCatch(summary(object$envelope), error = function(e) NULL)

  out <- list(
    method = object$method,
    nsim = object$settings$nsim,
    alpha = object$settings$alpha,
    type = object$settings$type,
    alternative = object$settings$alternative,
    correction = object$settings$correction,
    K_type = object$observed$type,
    GET_summary = envelope_summary
  )

  class(out) <- "summary_pp_envelope"
  out
}


#' @rdname pp_envelope-methods
#' @method print summary_pp_envelope
#' @export
print.summary_pp_envelope <- function(x, ...) {
  cat("Global-envelope diagnostic summary\n\n")
  cat("Method: ", x$method, "\n", sep = "")
  cat("Primary simulations: ", x$nsim, "\n", sep = "")
  cat("Alpha: ", x$alpha, "\n", sep = "")
  cat("Ordering: ", x$type, "\n", sep = "")
  cat("Alternative: ", x$alternative, "\n", sep = "")
  cat("K type: ", x$K_type, "\n", sep = "")
  cat("Correction: ", x$correction, "\n", sep = "")

  if (!is.null(x$GET_summary)) {
    cat("\nGET summary:\n")
    print(x$GET_summary)
  }

  invisible(x)
}


#' @rdname pp_envelope-methods
#' @method plot pp_envelope
#' @export
plot.pp_envelope <- function(
    x,
    what = c("envelope", "residual"),
    ...
) {
  what <- match.arg(what)

  if (identical(what, "residual")) {
    plot(x$observed, what = "residual", ...)
    return(invisible(x))
  }

  if (identical(x$observed$type, "spatial")) {
    plotted <- plot(x$envelope, ...)

    # GET returns a ggplot object, which must be printed explicitly inside an
    # S3 method; otherwise nothing is drawn.
    if (inherits(plotted, "ggplot")) {
      plotted <- plotted + ggplot2::labs(x = "Spatial lag r", y = "K(r)")
      print(plotted)
    }

    return(invisible(x))
  }

  warning(
    paste0(
      "GET represents the spatio-temporal surface as a flattened vector. ",
      "Showing the observed residual surface instead; the full GET object is available in `$envelope`."
    ),
    call. = FALSE
  )
  plot(x$observed, what = "residual", ...)
  invisible(x)
}


#' @rdname pp_envelope-methods
#' @method print pp_diagnostic
#' @export
print.pp_diagnostic <- function(x, ...) {
  print.pp_envelope(x, ...)
}


#' @rdname pp_envelope-methods
#' @method summary pp_diagnostic
#' @export
summary.pp_diagnostic <- function(object, ...) {
  summary.pp_envelope(object, ...)
}


#' @rdname pp_envelope-methods
#' @method plot pp_diagnostic
#' @export
plot.pp_diagnostic <- function(x, ...) {
  plot.pp_envelope(x, ...)
}
