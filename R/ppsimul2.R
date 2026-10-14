# Poisson point-process simulation -----------------------------------------
#
# This module depends on the constructors and covariate tools defined in
# ppdata.R / ppcovariates.R, in particular:
# pp_window(), pp_pattern(), pp_grid(), evaluate_covariate(), and the
# pp_covariate classes.


# Internal validation and data helpers -------------------------------------

.pp_sim_validate_window <- function(window) {
  if (!inherits(window, "pp_window")) {
    stop(
      "`window` must be an object of class 'pp_window'.",
      call. = FALSE
    )
  }

  required <- c("bounds", "axis_info")
  missing_components <- setdiff(required, names(window))

  if (length(missing_components) > 0L) {
    stop(
      sprintf(
        "The supplied window is missing: %s.",
        paste(missing_components, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  axis_names <- window$axis_info$name

  if (!all(axis_names %in% colnames(window$bounds))) {
    stop(
      "The window bounds are inconsistent with its axis information.",
      call. = FALSE
    )
  }

  invisible(window)
}


.pp_sim_window_volume <- function(window) {
  widths <- window$bounds["upper", ] - window$bounds["lower", ]
  volume <- prod(as.numeric(widths))

  if (!is.finite(volume) || volume <= 0) {
    stop(
      "The observation-window hypervolume must be finite and strictly positive.",
      call. = FALSE
    )
  }

  volume
}


.pp_sim_coords_data_frame <- function(coords) {
  coords <- as.data.frame(coords, check.names = FALSE)
  rownames(coords) <- NULL
  coords
}


.pp_sim_validate_choice <- function(x, choices, arg_name) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    stop(
      sprintf("`%s` must be one character string.", arg_name),
      call. = FALSE
    )
  }

  x <- match.arg(x, choices)
  x
}


.pp_sim_validate_positive_integer <- function(x, arg_name, min_value = 1L) {
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


.pp_sim_validate_nonnegative_number <- function(x, arg_name) {
  if (!is.numeric(x) ||
      length(x) != 1L ||
      !is.finite(x) ||
      x < 0) {
    stop(
      sprintf(
        "`%s` must be one non-negative finite number.",
        arg_name
      ),
      call. = FALSE
    )
  }

  as.numeric(x)
}


.pp_sim_validate_positive_number <- function(x, arg_name, lower = 0) {
  if (!is.numeric(x) ||
      length(x) != 1L ||
      !is.finite(x) ||
      x <= lower) {
    comparator <- if (lower == 0) "strictly positive" else paste0("greater than ", lower)

    stop(
      sprintf(
        "`%s` must be one %s finite number.",
        arg_name,
        comparator
      ),
      call. = FALSE
    )
  }

  as.numeric(x)
}


.pp_sim_validate_seed <- function(seed) {
  if (is.null(seed)) {
    return(NULL)
  }

  if (!is.numeric(seed) ||
      length(seed) != 1L ||
      !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop(
      "`seed` must be NULL or one finite integer.",
      call. = FALSE
    )
  }

  as.integer(seed)
}


.pp_sim_uniform_candidates <- function(window, n) {
  axis_names <- window$axis_info$name
  d <- length(axis_names)

  coords <- matrix(
    NA_real_,
    nrow = n,
    ncol = d,
    dimnames = list(NULL, axis_names)
  )

  if (n == 0L) {
    return(coords)
  }

  for (j in seq_len(d)) {
    coords[, j] <- stats::runif(
      n = n,
      min = window$bounds["lower", axis_names[j]],
      max = window$bounds["upper", axis_names[j]]
    )
  }

  coords
}


.pp_sim_call_user_function <- function(fun, supplied_args, fun_name) {
  if (!is.function(fun)) {
    stop(
      sprintf("`%s` must be a function.", fun_name),
      call. = FALSE
    )
  }

  fun_formals <- formals(fun)

  if (is.null(fun_formals)) {
    return(do.call(fun, supplied_args))
  }

  formal_names <- names(fun_formals)
  accepts_dots <- "..." %in% formal_names

  args_to_use <- if (accepts_dots) {
    supplied_args
  } else {
    supplied_args[intersect(names(supplied_args), formal_names)]
  }

  tryCatch(
    do.call(fun, args_to_use),
    error = function(error) {
      stop(
        sprintf(
          "Evaluation of `%s` failed: %s",
          fun_name,
          conditionMessage(error)
        ),
        call. = FALSE
      )
    }
  )
}


# Covariate preparation ----------------------------------------------------

.pp_sim_normalize_covariates <- function(covariates) {
  if (is.null(covariates)) {
    return(list())
  }

  if (inherits(covariates, "pp_covariate")) {
    covariates <- list(covariates)
  }

  if (!is.list(covariates) || length(covariates) == 0L) {
    stop(
      "`covariates` must be NULL, one pp_covariate object, or a non-empty list of pp_covariate objects.",
      call. = FALSE
    )
  }

  for (covariate in covariates) {
    .pp_cov_validate_object(covariate)
  }

  covariate_names <- names(covariates)

  if (is.null(covariate_names) || any(!nzchar(covariate_names))) {
    covariate_names <- vapply(
      covariates,
      function(z) z$name,
      character(1)
    )
  }

  if (anyDuplicated(covariate_names)) {
    stop(
      "Covariates must have unique names, either in the list or in their `name` fields.",
      call. = FALSE
    )
  }

  names(covariates) <- covariate_names
  covariates
}


.pp_sim_normalize_covariate_methods <- function(covariates, methods) {
  if (length(covariates) == 0L) {
    if (!is.null(methods)) {
      stop(
        "`covariate_methods` can be supplied only when `covariates` are supplied.",
        call. = FALSE
      )
    }

    return(character(0))
  }

  covariate_names <- names(covariates)

  if (is.null(methods)) {
    methods <- vapply(
      covariates,
      .pp_cov_default_method,
      character(1)
    )
  } else {
    if (!is.character(methods)) {
      stop(
        "`covariate_methods` must be NULL or a character vector.",
        call. = FALSE
      )
    }

    if (!is.null(names(methods))) {
      missing_names <- setdiff(covariate_names, names(methods))
      extra_names <- setdiff(names(methods), covariate_names)

      if (length(missing_names) > 0L || length(extra_names) > 0L) {
        stop(
          "Named `covariate_methods` must match the covariate names exactly.",
          call. = FALSE
        )
      }

      methods <- methods[covariate_names]
    } else if (length(methods) != length(covariates)) {
      stop(
        "`covariate_methods` must have one entry per covariate.",
        call. = FALSE
      )
    }
  }

  methods <- vapply(
    seq_along(covariates),
    function(i) {
      .pp_cov_validate_method(
        covariate = covariates[[i]],
        method = methods[i]
      )
    },
    character(1)
  )

  stats::setNames(methods, covariate_names)
}


.pp_sim_normalize_covariate_args <- function(covariates, covariate_args) {
  if (length(covariates) == 0L) {
    if (!is.null(covariate_args)) {
      stop(
        "`covariate_args` can be supplied only when `covariates` are supplied.",
        call. = FALSE
      )
    }

    return(list())
  }

  covariate_names <- names(covariates)

  if (is.null(covariate_args)) {
    covariate_args <- list()
  }

  if (!is.list(covariate_args)) {
    stop(
      "`covariate_args` must be NULL or a named list.",
      call. = FALSE
    )
  }

  if (length(covariate_args) > 0L) {
    if (is.null(names(covariate_args)) || any(!nzchar(names(covariate_args)))) {
      stop(
        "`covariate_args` must be a named list.",
        call. = FALSE
      )
    }

    unknown_names <- setdiff(names(covariate_args), covariate_names)

    if (length(unknown_names) > 0L) {
      stop(
        sprintf(
          "`covariate_args` contains unknown covariate name(s): %s.",
          paste(unknown_names, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  allowed_names <- c("p", "k", "axis_scale", "chunk_size")
  out <- vector("list", length(covariates))
  names(out) <- covariate_names

  for (i in seq_along(covariates)) {
    name_i <- covariate_names[i]
    args_i <- covariate_args[[name_i]]

    if (is.null(args_i)) {
      args_i <- list()
    }

    if (!is.list(args_i)) {
      stop(
        sprintf("`covariate_args[[\"%s\"]]` must be a list.", name_i),
        call. = FALSE
      )
    }

    bad_names <- setdiff(names(args_i), allowed_names)

    if (length(bad_names) > 0L) {
      stop(
        sprintf(
          "`covariate_args[[\"%s\"]]` contains unsupported option(s): %s.",
          name_i,
          paste(bad_names, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    # Resolve scale only once. This avoids repeating the space-time scaling
    # warning for the bound grid and for every thinning restart.
    if ("axis_scale" %in% names(args_i)) {
      args_i$axis_scale <- .pp_cov_validate_axis_scale(
        covariate = covariates[[i]],
        axis_scale = args_i$axis_scale
      )
    } else {
      args_i$axis_scale <- .pp_cov_validate_axis_scale(
        covariate = covariates[[i]],
        axis_scale = NULL
      )
    }

    out[[i]] <- args_i
  }

  out
}


.pp_sim_evaluate_covariates <- function(
    covariates,
    coords,
    methods,
    covariate_args
) {
  if (length(covariates) == 0L) {
    return(NULL)
  }

  target <- .pp_sim_coords_data_frame(coords)
  out <- data.frame(row.names = seq_len(nrow(target)))

  for (i in seq_along(covariates)) {
    name_i <- names(covariates)[i]

    eval_args <- c(
      list(
        covariate = covariates[[i]],
        at = target,
        method = methods[i]
      ),
      covariate_args[[i]]
    )

    out[[name_i]] <- do.call(evaluate_covariate, eval_args)
  }

  rownames(out) <- NULL
  out
}


# Mark preparation ---------------------------------------------------------

.pp_sim_normalize_mark_generators <- function(mark_generators, axis_names, covariate_names) {
  if (is.null(mark_generators)) {
    return(list())
  }

  if (!is.list(mark_generators) || length(mark_generators) == 0L) {
    stop(
      "`mark_generators` must be NULL or a non-empty named list of functions.",
      call. = FALSE
    )
  }

  mark_names <- names(mark_generators)

  if (is.null(mark_names) ||
      any(!nzchar(mark_names)) ||
      anyDuplicated(mark_names)) {
    stop(
      "`mark_generators` must have unique non-empty names.",
      call. = FALSE
    )
  }

  if (!all(vapply(mark_generators, is.function, logical(1)))) {
    stop(
      "Every element of `mark_generators` must be a function.",
      call. = FALSE
    )
  }

  coordinate_conflicts <- intersect(mark_names, axis_names)
  covariate_conflicts <- intersect(mark_names, covariate_names)

  if (length(coordinate_conflicts) > 0L) {
    stop(
      sprintf(
        "Mark name(s) conflict with coordinate axes: %s.",
        paste(coordinate_conflicts, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (length(covariate_conflicts) > 0L) {
    stop(
      sprintf(
        "Mark name(s) conflict with covariate names: %s.",
        paste(covariate_conflicts, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  mark_generators
}


.pp_sim_validate_mark_value <- function(value, n, mark_name) {
  if (is.matrix(value) || is.data.frame(value) || is.list(value)) {
    stop(
      sprintf(
        "Mark generator `%s` must return an atomic vector, not a matrix, data frame, or list.",
        mark_name
      ),
      call. = FALSE
    )
  }

  if (length(value) != n) {
    stop(
      sprintf(
        "Mark generator `%s` returned %d values; %d were required.",
        mark_name,
        length(value),
        n
      ),
      call. = FALSE
    )
  }

  if (is.numeric(value)) {
    value <- as.numeric(value)

    if (any(!is.finite(value))) {
      stop(
        sprintf(
          "Continuous mark `%s` contains missing or non-finite values.",
          mark_name
        ),
        call. = FALSE
      )
    }

    return(value)
  }

  if (is.character(value)) {
    if (anyNA(value)) {
      stop(
        sprintf(
          "Categorical mark `%s` contains missing values.",
          mark_name
        ),
        call. = FALSE
      )
    }

    return(factor(value))
  }

  if (is.factor(value)) {
    if (anyNA(value)) {
      stop(
        sprintf(
          "Categorical mark `%s` contains missing values.",
          mark_name
        ),
        call. = FALSE
      )
    }

    return(value)
  }

  stop(
    sprintf(
      paste0(
        "Mark generator `%s` must return either a continuous numeric vector ",
        "or a factor/character vector."
      ),
      mark_name
    ),
    call. = FALSE
  )
}


.pp_sim_generate_marks <- function(
    mark_generators,
    n,
    coords,
    covariates,
    params
) {
  if (length(mark_generators) == 0L) {
    return(NULL)
  }

  coords <- .pp_sim_coords_data_frame(coords)
  out <- vector("list", length(mark_generators))
  names(out) <- names(mark_generators)

  for (i in seq_along(mark_generators)) {
    mark_name <- names(mark_generators)[i]

    value <- .pp_sim_call_user_function(
      fun = mark_generators[[i]],
      supplied_args = list(
        n = n,
        coords = coords,
        covariates = covariates,
        params = params
      ),
      fun_name = paste0("mark_generators$", mark_name)
    )

    out[[i]] <- .pp_sim_validate_mark_value(
      value = value,
      n = n,
      mark_name = mark_name
    )
  }

  out <- as.data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}


# Intensity and bound helpers ----------------------------------------------

.pp_sim_evaluate_log_intensity <- function(
    log_intensity,
    coords,
    covariates,
    marks,
    params
) {
  coords <- .pp_sim_coords_data_frame(coords)
  n <- nrow(coords)

  eta <- .pp_sim_call_user_function(
    fun = log_intensity,
    supplied_args = list(
      coords = coords,
      covariates = covariates,
      marks = marks,
      params = params
    ),
    fun_name = "log_intensity"
  )

  if (!is.numeric(eta) || length(eta) != n) {
    stop(
      paste0(
        "`log_intensity` must return a numeric vector with one value per ",
        "candidate location."
      ),
      call. = FALSE
    )
  }

  eta <- as.numeric(eta)

  if (anyNA(eta) || any(is.nan(eta)) || any(is.infinite(eta) & eta > 0)) {
    stop(
      paste0(
        "`log_intensity` returned NA, NaN, or +Inf. Only finite values and ",
        "-Inf (zero intensity) are allowed."
      ),
      call. = FALSE
    )
  }

  max_log <- log(.Machine$double.xmax)

  if (any(eta > max_log)) {
    stop(
      "`log_intensity` is too large to exponentiate safely.",
      call. = FALSE
    )
  }

  eta
}


.pp_sim_lambda_from_eta <- function(eta) {
  lambda <- numeric(length(eta))
  finite <- is.finite(eta)

  lambda[finite] <- exp(eta[finite])
  lambda[is.infinite(eta) & eta < 0] <- 0

  lambda
}


.pp_sim_resolve_bound_n <- function(window, bound_n, bound_max_points) {
  axis_names <- window$axis_info$name
  d <- length(axis_names)

  bound_max_points <- .pp_sim_validate_positive_integer(
    bound_max_points,
    arg_name = "bound_max_points",
    min_value = 2L
  )

  if (is.null(bound_n)) {
    per_axis <- max(
      2L,
      floor(bound_max_points^(1 / d))
    )

    n <- rep(as.integer(per_axis), d)
    names(n) <- axis_names

    return(n)
  }

  if (!is.numeric(bound_n) ||
      any(!is.finite(bound_n)) ||
      any(bound_n != as.integer(bound_n))) {
    stop(
      "`bound_n` must contain finite integer grid resolutions.",
      call. = FALSE
    )
  }

  if (!is.null(names(bound_n))) {
    if (anyDuplicated(names(bound_n))) {
      stop("`bound_n` cannot contain duplicated axis names.", call. = FALSE)
    }

    missing_axes <- setdiff(axis_names, names(bound_n))
    extra_axes <- setdiff(names(bound_n), axis_names)

    if (length(missing_axes) > 0L || length(extra_axes) > 0L) {
      stop(
        sprintf(
          "`bound_n` must be named exactly with these axes: %s.",
          paste(axis_names, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    bound_n <- bound_n[axis_names]
  } else if (length(bound_n) == 1L) {
    bound_n <- rep(bound_n, d)
  } else if (length(bound_n) != d) {
    stop(
      sprintf("`bound_n` must have length 1 or length %d.", d),
      call. = FALSE
    )
  }

  bound_n <- as.integer(bound_n)

  if (any(bound_n < 2L)) {
    stop("Every bound-grid axis must contain at least two points.", call. = FALSE)
  }

  names(bound_n) <- axis_names

  if (prod(as.double(bound_n)) > bound_max_points) {
    stop(
      sprintf(
        paste0(
          "The requested bound grid contains %.0f points, exceeding ",
          "`bound_max_points` (%.0f)."
        ),
        prod(as.double(bound_n)),
        bound_max_points
      ),
      call. = FALSE
    )
  }

  bound_n
}


.pp_sim_grid_bound <- function(
    window,
    log_intensity,
    covariates,
    covariate_methods,
    covariate_args,
    mark_generators,
    mark_mode,
    params,
    bound_n,
    bound_max_points,
    bound_safety,
    bound_mark_reps
) {
  resolution <- .pp_sim_resolve_bound_n(
    window = window,
    bound_n = bound_n,
    bound_max_points = bound_max_points
  )

  grid <- pp_grid(
    window = window,
    n = resolution,
    max_points = bound_max_points
  )

  grid_coords <- .pp_sim_coords_data_frame(grid$coords)

  grid_covariates <- .pp_sim_evaluate_covariates(
    covariates = covariates,
    coords = grid_coords,
    methods = covariate_methods,
    covariate_args = covariate_args
  )

  n_grid <- nrow(grid_coords)
  n_reps <- if (identical(mark_mode, "joint") && length(mark_generators) > 0L) {
    bound_mark_reps
  } else {
    1L
  }

  lambda_sum <- numeric(n_grid)
  lambda_max <- 0

  for (replicate_id in seq_len(n_reps)) {
    grid_marks <- if (identical(mark_mode, "joint")) {
      .pp_sim_generate_marks(
        mark_generators = mark_generators,
        n = n_grid,
        coords = grid_coords,
        covariates = grid_covariates,
        params = params
      )
    } else {
      NULL
    }

    eta <- .pp_sim_evaluate_log_intensity(
      log_intensity = log_intensity,
      coords = grid_coords,
      covariates = grid_covariates,
      marks = grid_marks,
      params = params
    )

    lambda <- .pp_sim_lambda_from_eta(eta)

    lambda_sum <- lambda_sum + lambda
    lambda_max <- max(lambda_max, max(lambda, 0))
  }

  mean_lambda <- lambda_sum / n_reps
  volume <- .pp_sim_window_volume(window)

  list(
    lambda_upper = bound_safety * lambda_max,
    estimated_expected_n = mean(mean_lambda) * volume,
    resolution = resolution,
    n_grid = n_grid,
    n_mark_reps = n_reps,
    max_lambda_on_grid = lambda_max
  )
}


.pp_sim_validate_lambda_upper <- function(lambda_upper, arg_name = "lambda_upper") {
  if (!is.numeric(lambda_upper) ||
      length(lambda_upper) != 1L ||
      !is.finite(lambda_upper) ||
      lambda_upper < 0) {
    stop(
      sprintf("`%s` must be one non-negative finite number.", arg_name),
      call. = FALSE
    )
  }

  as.numeric(lambda_upper)
}


# Pattern construction and output -----------------------------------------

.pp_sim_build_pattern <- function(
    coords,
    window,
    marks = NULL,
    covariates = NULL,
    include_marks = TRUE,
    include_covariates = TRUE,
    simulation_metadata = NULL
) {
  axis_names <- window$axis_info$name
  coords <- .pp_sim_coords_data_frame(coords)
  coords <- coords[, axis_names, drop = FALSE]

  data <- coords
  mark_names <- character(0)

  if (isTRUE(include_marks) && !is.null(marks) && ncol(marks) > 0L) {
    marks <- as.data.frame(marks, check.names = FALSE, stringsAsFactors = FALSE)
    rownames(marks) <- NULL

    if (nrow(marks) != nrow(coords)) {
      stop(
        "The number of mark rows does not match the number of simulated events.",
        call. = FALSE
      )
    }

    mark_names <- names(marks)
    data <- cbind(data, marks)
  }

  pattern <- pp_pattern(
    data = data,
    coords = axis_names,
    marks = if (length(mark_names) == 0L) NULL else mark_names,
    window = window
  )

  if (isTRUE(include_covariates) && !is.null(covariates)) {
    covariates <- as.data.frame(
      covariates,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    rownames(covariates) <- NULL

    if (nrow(covariates) != nrow(coords)) {
      stop(
        "The number of covariate rows does not match the number of simulated events.",
        call. = FALSE
      )
    }

    pattern$covariates <- covariates
  }

  if (!is.null(simulation_metadata)) {
    pattern$metadata$simulation <- simulation_metadata
  }

  pattern
}


.pp_sim_subset_or_null <- function(x, index) {
  if (is.null(x)) {
    return(NULL)
  }

  out <- x[index, , drop = FALSE]
  rownames(out) <- NULL
  out
}


.pp_sim_single <- function(
    window,
    mode,
    expected_n,
    lambda0,
    log_intensity,
    params,
    covariates,
    covariate_methods,
    covariate_args,
    mark_generators,
    mark_mode,
    lambda_upper,
    estimated_expected_n,
    bound_details,
    on_bound_violation,
    restart_safety,
    max_restarts,
    max_candidates,
    output
) {
  volume <- .pp_sim_window_volume(window)
  axis_names <- window$axis_info$name
  restarts <- 0L
  current_lambda_upper <- lambda_upper

  repeat {
    expected_candidates <- current_lambda_upper * volume

    if (!is.finite(expected_candidates) || expected_candidates > max_candidates) {
      stop(
        sprintf(
          paste0(
            "The thinning process would require an expected %.0f candidates, ",
            "exceeding `max_candidates` (%.0f)."
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
          paste0(
            "The thinning process generated %d candidates, exceeding ",
            "`max_candidates` (%d)."
          ),
          n_candidates,
          max_candidates
        ),
        call. = FALSE
      )
    }

    candidate_coords <- .pp_sim_uniform_candidates(
      window = window,
      n = n_candidates
    )

    candidate_coords_df <- .pp_sim_coords_data_frame(candidate_coords)

    candidate_covariates <- .pp_sim_evaluate_covariates(
      covariates = covariates,
      coords = candidate_coords_df,
      methods = covariate_methods,
      covariate_args = covariate_args
    )

    candidate_marks <- if (identical(mark_mode, "joint")) {
      .pp_sim_generate_marks(
        mark_generators = mark_generators,
        n = n_candidates,
        coords = candidate_coords_df,
        covariates = candidate_covariates,
        params = params
      )
    } else {
      NULL
    }

    if (identical(mode, "homogeneous")) {
      candidate_eta <- if (lambda0 == 0) {
        numeric(n_candidates)
      } else {
        rep(log(lambda0), n_candidates)
      }
      candidate_lambda <- rep(lambda0, n_candidates)
    } else {
      candidate_eta <- .pp_sim_evaluate_log_intensity(
        log_intensity = log_intensity,
        coords = candidate_coords_df,
        covariates = candidate_covariates,
        marks = candidate_marks,
        params = params
      )

      candidate_lambda <- .pp_sim_lambda_from_eta(candidate_eta)
    }

    tolerance <- .Machine$double.eps^0.5 * max(1, current_lambda_upper)
    violation <- candidate_lambda > current_lambda_upper + tolerance

    if (any(violation)) {
      observed_max <- max(candidate_lambda[violation])

      if (identical(on_bound_violation, "error")) {
        stop(
          sprintf(
            paste0(
              "The thinning bound was violated: the largest candidate intensity ",
              "was %.8g, above lambda_upper = %.8g. Supply a larger user bound, ",
              "increase `bound_n`/`bound_safety`, or use `on_bound_violation = 'restart'`."
            ),
            observed_max,
            current_lambda_upper
          ),
          call. = FALSE
        )
      }

      restarts <- restarts + 1L

      if (restarts > max_restarts) {
        stop(
          sprintf(
            "The thinning bound was violated more than `max_restarts` (%d) times.",
            max_restarts
          ),
          call. = FALSE
        )
      }

      current_lambda_upper <- max(
        current_lambda_upper * restart_safety,
        observed_max * restart_safety
      )

      next
    }

    if (current_lambda_upper == 0) {
      acceptance_probability <- numeric(n_candidates)
    } else {
      acceptance_probability <- pmin(
        1,
        candidate_lambda / current_lambda_upper
      )
    }

    uniform_draw <- stats::runif(n_candidates)
    accepted <- uniform_draw <= acceptance_probability

    accepted_coords <- candidate_coords[accepted, , drop = FALSE]
    accepted_covariates <- .pp_sim_subset_or_null(
      candidate_covariates,
      accepted
    )

    accepted_marks <- if (identical(mark_mode, "joint")) {
      .pp_sim_subset_or_null(candidate_marks, accepted)
    } else {
      .pp_sim_generate_marks(
        mark_generators = mark_generators,
        n = sum(accepted),
        coords = .pp_sim_coords_data_frame(accepted_coords),
        covariates = accepted_covariates,
        params = params
      )
    }

    simulation_metadata <- list(
      mode = mode,
      mark_mode = mark_mode,
      lambda_upper = current_lambda_upper,
      estimated_expected_n = estimated_expected_n,
      observed_n = sum(accepted),
      n_candidates = n_candidates,
      acceptance_rate = if (n_candidates == 0L) 0 else mean(accepted),
      thinning_restarts = restarts
    )

    enriched_pattern <- .pp_sim_build_pattern(
      coords = .pp_sim_coords_data_frame(accepted_coords),
      window = window,
      marks = accepted_marks,
      covariates = accepted_covariates,
      include_marks = TRUE,
      include_covariates = TRUE,
      simulation_metadata = simulation_metadata
    )

    coordinate_pattern <- .pp_sim_build_pattern(
      coords = .pp_sim_coords_data_frame(accepted_coords),
      window = window,
      include_marks = FALSE,
      include_covariates = FALSE,
      simulation_metadata = simulation_metadata
    )

    if (identical(output, "pattern")) {
      return(coordinate_pattern)
    }

    if (identical(output, "enriched")) {
      return(enriched_pattern)
    }

    diagnostic <- list(
      pattern = enriched_pattern,
      candidates = candidate_coords_df[, axis_names, drop = FALSE],
      candidate_covariates = candidate_covariates,
      candidate_marks = candidate_marks,
      log_intensity = candidate_eta,
      intensity = candidate_lambda,
      acceptance_probability = acceptance_probability,
      uniform_draw = uniform_draw,
      accepted = accepted,
      lambda_upper = current_lambda_upper,
      estimated_expected_n = estimated_expected_n,
      expected_candidates = expected_candidates,
      n_candidates = n_candidates,
      n_accepted = sum(accepted),
      acceptance_rate = if (n_candidates == 0L) 0 else mean(accepted),
      thinning_restarts = restarts,
      bound = bound_details,
      settings = list(
        mode = mode,
        expected_n = expected_n,
        mark_mode = mark_mode,
        output = output
      )
    )

    class(diagnostic) <- "pp_simulation"
    return(diagnostic)
  }
}


#' Simulate multidimensional Poisson point processes
#'
#' Simulates homogeneous or inhomogeneous Poisson point processes inside an
#' arbitrary-dimensional hyperrectangular observation window. Inhomogeneous
#' intensities may depend on event coordinates, externally supplied covariates,
#' generated marks, and user-defined model parameters. `simulate_pp()` is the
#' general Poisson-process simulator used by `ppnexus`.
#' A homogeneous process is specified through its expected number of events,
#' whereas an inhomogeneous process is specified through a function that
#' evaluates the log-intensity. Both cases are simulated by independent
#' thinning of a homogeneous Poisson candidate process.
#'
#' @param window A finite hyperrectangular observation window created by
#'   [pp_window()]. Its coordinate names, axis roles, units, temporal origin,
#'   and coordinate reference system are inherited by the simulated pattern.
#' @param expected_n A single non-negative finite number giving the expected
#'   number of events in `window` for a homogeneous Poisson process. Supply
#'   exactly one of `expected_n` and `log_intensity`.
#' @param log_intensity A function that returns one log-intensity value for each
#'   candidate location. The function may declare any subset of the named
#'   arguments `coords`, `covariates`, `marks`, and `params`. It must return a
#'   numeric vector of the same length as the number of candidate points.
#'   Finite values and `-Inf`, representing zero intensity, are allowed; `NA`,
#'   `NaN`, and `+Inf` are not allowed. Supply exactly one of `expected_n` and
#'   `log_intensity`.
#' @param params A list of user-defined parameters passed unchanged to
#'   `log_intensity` and to every function in `mark_generators`. `NULL` is
#'   converted to an empty list.
#' @param covariates `NULL`, one object created by [pp_covariate()], or a named
#'   list of such objects. When a list is supplied, its names identify the
#'   columns made available to `log_intensity` and to mark generators. Names
#'   must be unique; names used by `mark_generators` cannot duplicate them.
#' @param covariate_methods Optional character vector selecting the evaluation
#'   method for each supplied covariate. It may be named by covariate or have one
#'   entry per covariate. Continuous fields support `"idw"` and `"nearest"`;
#'   geometric supports use `"distance_min"`. When `NULL`, each covariate's
#'   default method is used. See [evaluate_covariate()].
#' @param covariate_args Optional named list with one sub-list per covariate.
#'   Each sub-list may contain `p`, `k`, `axis_scale`, and `chunk_size`, which
#'   are passed to [evaluate_covariate()]. Unspecified covariates use their
#'   default evaluation settings.
#' @param mark_generators `NULL` or a named list of functions used to generate
#'   event marks. Each generator may declare any subset of the arguments `n`,
#'   `coords`, `covariates`, and `params`, and must return an atomic vector of
#'   length `n`. Numeric vectors create continuous marks; character or factor
#'   vectors create categorical marks. Mark names must be unique and cannot
#'   duplicate coordinate or covariate names.
#' @param mark_mode Character string controlling when marks are generated.
#'   With `"posthoc"`, marks are generated only for retained events and cannot
#'   affect the intensity. With `"joint"`, marks are generated for every
#'   candidate before `log_intensity` is evaluated and may therefore affect the
#'   thinning probability.
#' @param bound Character string selecting the intensity-bound strategy. With
#'   `"user"`, `lambda_upper` must be supplied. With `"grid"`, the function
#'   estimates an upper bound numerically on a regular grid and multiplies the
#'   largest evaluated intensity by `bound_safety`.
#' @param lambda_upper A single non-negative finite number giving a verified
#'   upper bound for the intensity when `bound = "user"`. It must bound the
#'   intensity over the complete coordinate domain and, for joint marked
#'   models, over all marks that can be generated.
#' @param bound_n Optional positive integer resolution of the numerical bound
#'   grid. A scalar is recycled over all axes; alternatively, supply one value
#'   per axis or a vector named exactly with the window axes. Every axis must
#'   contain at least two grid points. When `NULL`, a balanced resolution is
#'   selected subject to `bound_max_points`.
#' @param bound_max_points Positive integer giving the maximum total number of
#'   locations allowed in the numerical bound grid. This controls memory and
#'   computation when `bound = "grid"`.
#' @param bound_safety Multiplicative safety factor applied to the largest
#'   intensity found on the numerical grid. It must be at least one. A larger
#'   value reduces the risk of a bound violation but lowers the thinning
#'   acceptance rate.
#' @param bound_mark_reps Positive integer giving the number of independent mark
#'   realizations generated at each bound-grid location when
#'   `mark_mode = "joint"`. It is ignored for post-hoc marks and models without
#'   generated marks.
#' @param on_bound_violation Character string specifying the action when a
#'   candidate intensity exceeds the current bound. `"error"` stops and reports
#'   the violation. `"restart"` discards the complete candidate realization,
#'   increases the bound, and starts the thinning run again.
#' @param restart_safety Multiplicative factor, strictly greater than one, used
#'   to enlarge the intensity bound after an observed violation when
#'   `on_bound_violation = "restart"`.
#' @param max_restarts Positive integer giving the maximum number of automatic
#'   bound-violation restarts permitted for each simulation.
#' @param output Character string controlling the returned object. `"pattern"`
#'   returns a coordinate-only `pp` object; `"enriched"` returns a `pp` object
#'   containing generated marks and evaluated covariates; `"diagnostic"`
#'   returns a `pp_simulation` object retaining the enriched pattern and the
#'   complete thinning diagnostics.
#' @param nsim Positive integer giving the number of independent simulations.
#'   A single object is returned when `nsim = 1`; otherwise a named list is
#'   returned.
#' @param seed `NULL` or one finite integer used to initialize R's random-number
#'   generator before all simulations. For `nsim > 1`, one reproducible random
#'   stream is used across the successive simulations.
#' @param max_candidates Positive integer limiting both the expected and the
#'   realized number of homogeneous candidate points in each thinning run. The
#'   simulation stops before allocating an excessively large candidate set.
#'
#' @details
#' Let \eqn{W} denote the observation window, \eqn{|W|} its hypervolume,
#' and \eqn{\mu_W = \mathrm{E}[N(W)]}. For a homogeneous model,
#' \eqn{\mu_W} is specified by `expected_n` and the constant intensity is
#'
#' \deqn{\lambda_0 = \frac{\mu_W}{|W|}.}
#'
#' Consequently, \eqn{N(W)} has a Poisson distribution with mean
#' `expected_n`, and conditional on \eqn{N(W)=n}, the event locations are
#' independent and uniformly distributed over \eqn{W}.
#'
#' For an inhomogeneous model with intensity \eqn{\lambda(u)}, the thinning
#' algorithm first simulates a homogeneous candidate process with constant
#' intensity `lambda_upper`. A candidate at \eqn{u} is retained independently
#' with probability
#'
#' \deqn{p(u)=\frac{\lambda(u)}{\lambda_{\max}}.}
#'
#' A user-supplied bound is the preferred option whenever a valid analytical
#' upper bound is available. The bound obtained with `bound = "grid"` is only a
#' numerical estimate: neither a fine grid nor a safety factor provides a
#' mathematical guarantee that unsampled intensity peaks are bounded. The
#' restart mechanism corrects observed violations by regenerating the complete
#' candidate process, but it cannot detect violations at locations that were
#' never proposed.
#'
#' In joint marked mode, the functions in `mark_generators` define the proposal
#' or base distribution of the marks. For every candidate location \eqn{u}, a
#' mark vector \eqn{m} is generated and retained with probability proportional
#' to \eqn{\lambda(u,m)}. Thus `log_intensity` reweights the generated mark
#' distribution, and `lambda_upper` must be valid over every mark value that the
#' generators can produce. In post-hoc mode, marks are conditionally generated
#' only after locations have been accepted and do not alter the point-process
#' intensity.
#'
#' External covariates are evaluated both on candidate locations and, when
#' needed, on the numerical bound grid. In a spatio-temporal setting, carefully
#' choose `axis_scale` because Euclidean interpolation combines spatial and
#' temporal coordinate differences.
#'
#' @section Diagnostic output:
#' With `output = "diagnostic"`, the returned `pp_simulation` object contains:
#'
#' * `pattern`: the accepted enriched point pattern;
#' * `candidates`: candidate coordinates from the final thinning run;
#' * `candidate_covariates` and `candidate_marks`: values evaluated or generated
#'   for those candidates;
#' * `log_intensity` and `intensity`: candidate log-intensities and intensities;
#' * `acceptance_probability`, `uniform_draw`, and `accepted`: the quantities
#'   defining each thinning decision;
#' * `lambda_upper`, `n_candidates`, `n_accepted`, `acceptance_rate`, and
#'   `thinning_restarts`: simulation diagnostics;
#' * `estimated_expected_n`: the requested homogeneous mean or, for a grid
#'   bound, a numerical integral estimate of the expected count;
#' * `bound` and `settings`: the bound construction and simulation settings.
#'
#' @return
#' If `nsim = 1`, the result depends on `output`:
#'
#' * `"pattern"`: an unmarked coordinate-only object inheriting from `"pp"`;
#' * `"enriched"`: an object inheriting from `"pp"` with generated marks and,
#'   when supplied, evaluated covariates;
#' * `"diagnostic"`: an object of class `"pp_simulation"` described above.
#'
#' If `nsim > 1`, a named list `sim_1`, ..., `sim_nsim` of the corresponding
#' objects is returned.
#'
#' @references
#' Lewis, P. A. W., and Shedler, G. S. (1979). Simulation of nonhomogeneous
#' Poisson processes by thinning. *Naval Research Logistics Quarterly*, 26(3),
#' 403--413. \doi{10.1002/nav.3800260304}.
#'
#' Møller, J., and Waagepetersen, R. P. (2004). *Statistical Inference and
#' Simulation for Spatial Point Processes*. Chapman & Hall/CRC.
#'
#' Daley, D. J., and Vere-Jones, D. (2003). *An Introduction to the Theory of
#' Point Processes. Volume I: Elementary Theory and Methods*, 2nd ed. Springer.
#'
#' @examples
#' # A common two-dimensional observation window
#' w <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1),
#'   axis_units = c(x = "km", y = "km")
#' )
#'
#' # Homogeneous Poisson process
#' x_hom <- simulate_pp(
#'   window = w,
#'   expected_n = 40,
#'   output = "pattern",
#'   seed = 123
#' )
#' summary(x_hom)
#'
#' # Inhomogeneous process depending on coordinates
#' # eta(x, y) = beta0 + beta_x * x + beta_y * y.
#' x_inhom <- simulate_pp(
#'   window = w,
#'   log_intensity = function(coords, params) {
#'     params$beta0 +
#'       params$beta_x * coords$x +
#'       params$beta_y * coords$y
#'   },
#'   params = list(
#'     beta0 = log(30),
#'     beta_x = 0.8,
#'     beta_y = -0.4
#'   ),
#'   bound = "user",
#'   lambda_upper = exp(log(30) + 0.8),
#'   output = "diagnostic",
#'   seed = 123
#' )
#' summary(x_inhom)
#'
#' # Post-hoc continuous and categorical marks
#' x_posthoc <- simulate_pp(
#'   window = w,
#'   expected_n = 35,
#'   mark_generators = list(
#'     magnitude = function(n) stats::runif(n, 1, 5),
#'     type = function(n) sample(c("A", "B"), n, replace = TRUE)
#'   ),
#'   mark_mode = "posthoc",
#'   output = "enriched",
#'   seed = 321
#' )
#' head(x_posthoc$marks)
#'
#' # Joint marked process: generated marks affect acceptance
#' x_joint <- simulate_pp(
#'   window = w,
#'   log_intensity = function(coords, marks, params) {
#'     params$beta0 + params$beta_x * coords$x +
#'       params$beta_m * marks$magnitude +
#'       ifelse(marks$type == "B", params$beta_B, 0)
#'   },
#'   params = list(
#'     beta0 = log(20), beta_x = 0.5,
#'     beta_m = 0.3, beta_B = 0.4
#'   ),
#'   mark_generators = list(
#'     magnitude = function(n) stats::runif(n),
#'     type = function(n) sample(c("A", "B"), n, replace = TRUE)
#'   ),
#'   mark_mode = "joint",
#'   bound = "user",
#'   lambda_upper = exp(log(20) + 0.5 + 0.3 + 0.4),
#'   output = "enriched",
#'   seed = 456
#' )
#' table(x_joint$marks$type)
#'
#' # Process driven by an external spatial covariate
#' cov_data <- expand.grid(
#'   x = seq(0, 1, length.out = 5),
#'   y = seq(0, 1, length.out = 5)
#' )
#' cov_data$temperature <- cov_data$x + cov_data$y
#'
#' temperature <- pp_covariate(
#'   data = cov_data,
#'   coords = c("x", "y"),
#'   value = "temperature",
#'   window = w,
#'   name = "temperature"
#' )
#'
#' x_covariate <- simulate_pp(
#'   window = w,
#'   log_intensity = function(covariates, params) {
#'     params$beta0 + params$beta_temperature * covariates$temperature
#'   },
#'   params = list(beta0 = log(15), beta_temperature = 0.5),
#'   covariates = list(temperature = temperature),
#'   covariate_methods = c(temperature = "idw"),
#'   covariate_args = list(temperature = list(k = 6)),
#'   bound = "user",
#'   lambda_upper = exp(log(15) + 0.5 * 2),
#'   output = "enriched",
#'   seed = 789
#' )
#' head(x_covariate$covariates)
#'
#' # Several independent simulations
#' sims <- simulate_pp(
#'   window = w,
#'   expected_n = 20,
#'   nsim = 3,
#'   seed = 42
#' )
#' vapply(sims, function(z) nrow(z$coords), numeric(1))
#'
#' @seealso
#' [pp_window()], [pp_pattern()], [pp_covariate()], [evaluate_covariate()],
#' [pp_grid()], [fit_pp()], and [simulate_lgcp_pp()].
#' @family Poisson-process simulation
#' @export
simulate_pp <- function(
    window,
    expected_n = NULL,
    log_intensity = NULL,
    params = list(),
    covariates = NULL,
    covariate_methods = NULL,
    covariate_args = NULL,
    mark_generators = NULL,
    mark_mode = c("posthoc", "joint"),
    bound = c("grid", "user"),
    lambda_upper = NULL,
    bound_n = NULL,
    bound_max_points = 50000L,
    bound_safety = 1.20,
    bound_mark_reps = 20L,
    on_bound_violation = c("error", "restart"),
    restart_safety = 1.20,
    max_restarts = 10L,
    output = c("pattern", "enriched", "diagnostic"),
    nsim = 1L,
    seed = NULL,
    max_candidates = 1e6
) {
  .pp_sim_validate_window(window)

  mark_mode <- match.arg(mark_mode)
  bound <- match.arg(bound)
  on_bound_violation <- match.arg(on_bound_violation)
  output <- match.arg(output)

  nsim <- .pp_sim_validate_positive_integer(
    nsim,
    arg_name = "nsim",
    min_value = 1L
  )

  max_candidates <- .pp_sim_validate_positive_integer(
    max_candidates,
    arg_name = "max_candidates",
    min_value = 1L
  )

  max_restarts <- .pp_sim_validate_positive_integer(
    max_restarts,
    arg_name = "max_restarts",
    min_value = 1L
  )

  bound_mark_reps <- .pp_sim_validate_positive_integer(
    bound_mark_reps,
    arg_name = "bound_mark_reps",
    min_value = 1L
  )

  bound_safety <- .pp_sim_validate_positive_number(
    bound_safety,
    arg_name = "bound_safety",
    lower = 0
  )

  if (bound_safety < 1) {
    stop(
      "`bound_safety` must be greater than or equal to 1.",
      call. = FALSE
    )
  }

  restart_safety <- .pp_sim_validate_positive_number(
    restart_safety,
    arg_name = "restart_safety",
    lower = 0
  )

  if (restart_safety <= 1) {
    stop(
      "`restart_safety` must be strictly greater than 1.",
      call. = FALSE
    )
  }

  if (is.null(params)) {
    params <- list()
  }

  if (!is.list(params)) {
    stop("`params` must be a list or NULL.", call. = FALSE)
  }

  has_expected_n <- !is.null(expected_n)
  has_log_intensity <- !is.null(log_intensity)

  if (identical(has_expected_n, has_log_intensity)) {
    stop(
      "Supply exactly one of `expected_n` (homogeneous) or `log_intensity` (inhomogeneous).",
      call. = FALSE
    )
  }

  mode <- if (has_expected_n) "homogeneous" else "inhomogeneous"

  if (has_expected_n) {
    expected_n <- .pp_sim_validate_nonnegative_number(
      expected_n,
      arg_name = "expected_n"
    )
  }

  if (has_log_intensity && !is.function(log_intensity)) {
    stop("`log_intensity` must be a function.", call. = FALSE)
  }

  covariates <- .pp_sim_normalize_covariates(covariates)
  covariate_methods <- .pp_sim_normalize_covariate_methods(
    covariates = covariates,
    methods = covariate_methods
  )
  covariate_args <- .pp_sim_normalize_covariate_args(
    covariates = covariates,
    covariate_args = covariate_args
  )

  mark_generators <- .pp_sim_normalize_mark_generators(
    mark_generators = mark_generators,
    axis_names = window$axis_info$name,
    covariate_names = names(covariates)
  )

  if (identical(mark_mode, "joint") &&
      !has_log_intensity &&
      length(mark_generators) > 0L) {
    warning(
      paste0(
        "`mark_mode = 'joint'` has no effect on a homogeneous process defined ",
        "only through `expected_n`; marks are generated before thinning but do ",
        "not alter the constant intensity."
      ),
      call. = FALSE
    )
  }

  seed <- .pp_sim_validate_seed(seed)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  volume <- .pp_sim_window_volume(window)

  lambda0 <- NULL
  bound_details <- NULL
  estimated_expected_n <- NULL

  if (identical(mode, "homogeneous")) {
    lambda0 <- expected_n / volume

    if (identical(bound, "user")) {
      if (is.null(lambda_upper)) {
        stop(
          "`lambda_upper` must be supplied when `bound = 'user'`.",
          call. = FALSE
        )
      }

      lambda_upper <- .pp_sim_validate_lambda_upper(lambda_upper)

      tolerance <- .Machine$double.eps^0.5 * max(1, lambda0)

      if (lambda_upper + tolerance < lambda0) {
        stop(
          paste0(
            "For a homogeneous process, `lambda_upper` must be at least ",
            "expected_n / volume(window)."
          ),
          call. = FALSE
        )
      }
    } else {
      lambda_upper <- lambda0
    }

    estimated_expected_n <- expected_n
    bound_details <- list(
      method = if (identical(bound, "user")) "user" else "exact_homogeneous",
      lambda_upper_initial = lambda_upper,
      estimated_expected_n = expected_n,
      resolution = NULL,
      max_lambda_on_grid = lambda0,
      bound_safety = NA_real_,
      bound_mark_reps = NA_integer_
    )
  } else if (identical(bound, "user")) {
    if (is.null(lambda_upper)) {
      stop(
        "`lambda_upper` must be supplied when `bound = 'user'`.",
        call. = FALSE
      )
    }

    lambda_upper <- .pp_sim_validate_lambda_upper(lambda_upper)
    estimated_expected_n <- NA_real_

    if (lambda_upper == 0) {
      warning(
        paste0(
          "`lambda_upper` is zero. The simulation will contain no candidates. ",
          "Ensure that the log-intensity is identically -Inf."
        ),
        call. = FALSE
      )
    }

    bound_details <- list(
      method = "user",
      lambda_upper_initial = lambda_upper,
      estimated_expected_n = NA_real_,
      resolution = NULL,
      max_lambda_on_grid = NA_real_,
      bound_safety = NA_real_,
      bound_mark_reps = NA_integer_
    )
  } else {
    warning(
      paste0(
        "`bound = 'grid'` uses a numerical upper-bound estimate and is not a ",
        "mathematical guarantee. For exact thinning, provide a verified ",
        "`lambda_upper` with `bound = 'user'`."
      ),
      call. = FALSE
    )

    grid_bound <- .pp_sim_grid_bound(
      window = window,
      log_intensity = log_intensity,
      covariates = covariates,
      covariate_methods = covariate_methods,
      covariate_args = covariate_args,
      mark_generators = mark_generators,
      mark_mode = mark_mode,
      params = params,
      bound_n = bound_n,
      bound_max_points = bound_max_points,
      bound_safety = bound_safety,
      bound_mark_reps = bound_mark_reps
    )

    lambda_upper <- grid_bound$lambda_upper
    estimated_expected_n <- grid_bound$estimated_expected_n

    bound_details <- list(
      method = "grid",
      lambda_upper_initial = lambda_upper,
      estimated_expected_n = estimated_expected_n,
      resolution = grid_bound$resolution,
      n_grid = grid_bound$n_grid,
      max_lambda_on_grid = grid_bound$max_lambda_on_grid,
      bound_safety = bound_safety,
      bound_mark_reps = grid_bound$n_mark_reps
    )
  }

  results <- vector("list", nsim)

  for (simulation_id in seq_len(nsim)) {
    results[[simulation_id]] <- .pp_sim_single(
      window = window,
      mode = mode,
      expected_n = expected_n,
      lambda0 = lambda0,
      log_intensity = log_intensity,
      params = params,
      covariates = covariates,
      covariate_methods = covariate_methods,
      covariate_args = covariate_args,
      mark_generators = mark_generators,
      mark_mode = mark_mode,
      lambda_upper = lambda_upper,
      estimated_expected_n = estimated_expected_n,
      bound_details = bound_details,
      on_bound_violation = on_bound_violation,
      restart_safety = restart_safety,
      max_restarts = max_restarts,
      max_candidates = max_candidates,
      output = output
    )
  }

  if (nsim == 1L) {
    return(results[[1]])
  }

  names(results) <- paste0("sim_", seq_len(nsim))
  results
}


#' Methods for Poisson thinning-simulation diagnostics
#'
#' Print and summarize the diagnostic object returned by
#' `simulate_pp(output = "diagnostic")`. `print.pp_simulation()` displays
#' the principal thinning diagnostics.
#' `summary.pp_simulation()` constructs a structured summary containing the
#' simulation settings, acceptance diagnostics, bound information, and a
#' summary of the accepted point pattern. The corresponding print method
#' formats that summary.
#'
#' @param x An object of class `"pp_simulation"` for
#'   `print.pp_simulation()`, or an object of class
#'   `"summary_pp_simulation"` for `print.summary_pp_simulation()`.
#' @param object An object of class `"pp_simulation"`, as returned by
#'   [simulate_pp()] with `output = "diagnostic"`.
#' @param digits Positive integer controlling the number of significant digits
#'   shown for rates, intensities, and expected counts.
#' @param ... Additional arguments. They are currently ignored.
#'
#' @return
#' `print.pp_simulation()` and `print.summary_pp_simulation()` invisibly return
#' `x`. `summary.pp_simulation()` returns a list of class
#' `"summary_pp_simulation"` containing simulation settings, accepted and
#' candidate counts, acceptance rate, final intensity bound, estimated expected
#' count, restart count, bound details, and the summary of the accepted pattern.
#'
#' @examples
#' w <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#'
#' sim <- simulate_pp(
#'   window = w,
#'   expected_n = 30,
#'   output = "diagnostic",
#'   seed = 123
#' )
#'
#' print(sim)
#' sim_summary <- summary(sim)
#' print(sim_summary)
#'
#' @seealso [simulate_pp()], [summary.pp()], and [print.pp()].
#' @family Poisson-process simulation
#' @rdname pp_simulation-methods
#' @method print pp_simulation
#' @export
print.pp_simulation <- function(x, digits = 4, ...) {
  cat("<pp_simulation>\n")
  cat("Poisson point-process simulation by thinning\n\n")

  cat("Mode: ", x$settings$mode, "\n", sep = "")
  cat("Accepted events: ", x$n_accepted, "\n", sep = "")
  cat("Candidates: ", x$n_candidates, "\n", sep = "")
  cat(
    "Acceptance rate: ",
    .pp_format_number(x$acceptance_rate, digits),
    "\n",
    sep = ""
  )
  cat(
    "Lambda upper: ",
    .pp_format_number(x$lambda_upper, digits),
    "\n",
    sep = ""
  )

  if (!is.null(x$estimated_expected_n) && is.finite(x$estimated_expected_n)) {
    cat(
      "Estimated expected events: ",
      .pp_format_number(x$estimated_expected_n, digits),
      "\n",
      sep = ""
    )
  }

  cat("Bound method: ", x$bound$method, "\n", sep = "")
  cat("Thinning restarts: ", x$thinning_restarts, "\n", sep = "")

  invisible(x)
}


#' @rdname pp_simulation-methods
#' @method summary pp_simulation
#' @export
summary.pp_simulation <- function(object, ...) {
  if (!inherits(object, "pp_simulation")) {
    stop("`object` must inherit from class 'pp_simulation'.", call. = FALSE)
  }

  out <- list(
    mode = object$settings$mode,
    mark_mode = object$settings$mark_mode,
    n_accepted = object$n_accepted,
    n_candidates = object$n_candidates,
    acceptance_rate = object$acceptance_rate,
    lambda_upper = object$lambda_upper,
    estimated_expected_n = if (is.null(object$estimated_expected_n)) {
      NA_real_
    } else {
      object$estimated_expected_n
    },
    thinning_restarts = object$thinning_restarts,
    bound = object$bound,
    pattern_summary = summary(object$pattern)
  )

  class(out) <- "summary_pp_simulation"
  out
}


#' @rdname pp_simulation-methods
#' @method print summary_pp_simulation
#' @export
print.summary_pp_simulation <- function(x, digits = 4, ...) {
  cat("Poisson point-process simulation by thinning\n\n")
  cat("Mode: ", x$mode, "\n", sep = "")
  cat("Mark mode: ", x$mark_mode, "\n", sep = "")
  cat("Accepted events: ", x$n_accepted, "\n", sep = "")
  cat("Candidates: ", x$n_candidates, "\n", sep = "")
  cat(
    "Acceptance rate: ",
    .pp_format_number(x$acceptance_rate, digits),
    "\n",
    sep = ""
  )
  cat(
    "Lambda upper: ",
    .pp_format_number(x$lambda_upper, digits),
    "\n",
    sep = ""
  )

  if (!is.null(x$estimated_expected_n) &&
      length(x$estimated_expected_n) == 1L &&
      is.finite(x$estimated_expected_n)) {
    cat(
      "Estimated expected events: ",
      .pp_format_number(x$estimated_expected_n, digits),
      "\n",
      sep = ""
    )
  }

  cat("Bound method: ", x$bound$method, "\n", sep = "")
  cat("Thinning restarts: ", x$thinning_restarts, "\n\n", sep = "")

  cat("Simulated-pattern summary:\n")
  print(x$pattern_summary, digits = digits)

  invisible(x)
}




