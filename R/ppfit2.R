# Poisson point-process fitting by multidimensional cubature ----------------
#
# This module depends on the data and covariate tools defined in ppdata.R,
# including pp_pattern(), pp_window(), pp_covariate(), and
# evaluate_covariate().  It requires mgcv in Imports.


# General utilities ----------------------------------------------------------

.pp_fit_validate_pattern <- function(pattern) {
  if (!inherits(pattern, "pp")) {
    stop("`pattern` must inherit from class 'pp'.", call. = FALSE)
  }

  if (is.null(pattern$coords) || !is.matrix(pattern$coords)) {
    stop("`pattern` must contain a coordinate matrix.", call. = FALSE)
  }

  if (is.null(pattern$window) || !inherits(pattern$window, "pp_window")) {
    stop(
      "`pattern` must have a supplied pp_window; fitting requires a finite observation window.",
      call. = FALSE
    )
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

  if (!all(axis_names %in% colnames(pattern$window$bounds))) {
    stop("The pattern window is inconsistent with its coordinate axes.", call. = FALSE)
  }

  invisible(pattern)
}


.pp_fit_validate_positive_integer <- function(x, arg_name, min_value = 1L) {
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


.pp_fit_validate_positive_number <- function(x, arg_name) {
  if (!is.numeric(x) ||
      length(x) != 1L ||
      !is.finite(x) ||
      x <= 0) {
    stop(
      sprintf("`%s` must be one strictly positive finite number.", arg_name),
      call. = FALSE
    )
  }

  as.numeric(x)
}


.pp_fit_validate_seed <- function(seed) {
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


.pp_fit_validate_formula <- function(formula) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula without a response, for example `~ x + y`.", call. = FALSE)
  }

  if (length(formula) != 2L) {
    stop(
      "`formula` must not contain a response; use a one-sided formula such as `~ x + y`.",
      call. = FALSE
    )
  }

  formula
}


.pp_fit_formula_variables <- function(formula) {
  unique(all.vars(formula))
}


.pp_fit_make_gam_formula <- function(formula) {
  formula <- .pp_fit_validate_formula(formula)

  stats::as.formula(
    call("~", as.name(".pp_response"), formula[[2L]]),
    env = environment(formula)
  )
}


.pp_fit_window_volume <- function(window) {
  axis_names <- window$axis_info$name
  widths <- window$bounds["upper", axis_names] - window$bounds["lower", axis_names]
  volume <- prod(as.numeric(widths))

  if (!is.finite(volume) || volume <= 0) {
    stop("The observation-window hypervolume must be finite and positive.", call. = FALSE)
  }

  volume
}


.pp_fit_same_units <- function(x, y) {
  x <- as.character(x)
  y <- as.character(y)

  same_missing <- is.na(x) & is.na(y)
  same_value <- !is.na(x) & !is.na(y) & x == y

  all(same_missing | same_value)
}


.pp_fit_call_user_function <- function(fun, supplied_args, fun_name) {
  if (!is.function(fun)) {
    stop(sprintf("`%s` must be a function.", fun_name), call. = FALSE)
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
        sprintf("Evaluation of `%s` failed: %s", fun_name, conditionMessage(error)),
        call. = FALSE
      )
    }
  )
}


# Joint continuous domain ----------------------------------------------------

.pp_fit_normalize_continuous_marks <- function(pattern, continuous_marks) {
  if (is.null(continuous_marks) ||
      (is.list(continuous_marks) && length(continuous_marks) == 0L)) {
    return(list())
  }

  if (!is.list(continuous_marks)) {
    stop(
      "`continuous_marks` must be NULL or a non-empty named list.",
      call. = FALSE
    )
  }

  mark_names <- names(continuous_marks)

  if (is.null(mark_names) ||
      any(!nzchar(mark_names)) ||
      anyDuplicated(mark_names)) {
    stop("`continuous_marks` must have unique non-empty names.", call. = FALSE)
  }

  if (is.null(pattern$marks) || !is.data.frame(pattern$marks)) {
    stop("`continuous_marks` were supplied but `pattern` has no marks.", call. = FALSE)
  }

  missing_marks <- setdiff(mark_names, names(pattern$marks))

  if (length(missing_marks) > 0L) {
    stop(
      sprintf(
        "Unknown continuous mark(s): %s.",
        paste(missing_marks, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  out <- vector("list", length(continuous_marks))
  names(out) <- mark_names

  for (i in seq_along(continuous_marks)) {
    mark_name <- mark_names[i]
    spec <- continuous_marks[[i]]

    if (!is.list(spec)) {
      stop(
        sprintf("`continuous_marks[[\"%s\"]]` must be a list.", mark_name),
        call. = FALSE
      )
    }

    allowed <- c("mode", "bounds", "method", "args")
    unknown <- setdiff(names(spec), allowed)

    if (length(unknown) > 0L) {
      stop(
        sprintf(
          "`continuous_marks[[\"%s\"]]` contains unsupported option(s): %s.",
          mark_name,
          paste(unknown, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    if (is.null(spec$mode) ||
        !is.character(spec$mode) ||
        length(spec$mode) != 1L ||
        !spec$mode %in% c("joint", "field")) {
      stop(
        sprintf(
          "`continuous_marks[[\"%s\"]]$mode` must be either 'joint' or 'field'.",
          mark_name
        ),
        call. = FALSE
      )
    }

    values <- pattern$marks[[mark_name]]

    if (!is.numeric(values) || any(!is.finite(values))) {
      stop(
        sprintf(
          "Continuous mark `%s` must be a finite numeric vector.",
          mark_name
        ),
        call. = FALSE
      )
    }

    values <- as.numeric(values)

    if (identical(spec$mode, "joint")) {
      bounds <- spec$bounds

      if (is.null(bounds)) {
        bounds <- range(values)

        warning(
          sprintf(
            paste0(
              "No bounds were supplied for joint continuous mark `%s`; ",
              "using its observed range [%.6g, %.6g]. This makes the joint ",
              "domain sample-dependent."
            ),
            mark_name,
            bounds[1],
            bounds[2]
          ),
          call. = FALSE
        )
      }

      if (!is.numeric(bounds) || length(bounds) != 2L || any(!is.finite(bounds))) {
        stop(
          sprintf(
            "`continuous_marks[[\"%s\"]]$bounds` must contain two finite numeric values.",
            mark_name
          ),
          call. = FALSE
        )
      }

      bounds <- sort(as.numeric(bounds))

      if (bounds[1] >= bounds[2]) {
        stop(
          sprintf(
            "Joint continuous mark `%s` requires bounds with lower < upper.",
            mark_name
          ),
          call. = FALSE
        )
      }

      tolerance <- .Machine$double.eps^0.5 * max(1, abs(bounds))

      if (any(values < bounds[1] - tolerance | values > bounds[2] + tolerance)) {
        stop(
          sprintf(
            "Observed values of joint continuous mark `%s` fall outside its declared bounds.",
            mark_name
          ),
          call. = FALSE
        )
      }

      if (!is.null(spec$method) || length(spec$args) > 0L) {
        warning(
          sprintf(
            "`method` and `args` are ignored for joint continuous mark `%s`.",
            mark_name
          ),
          call. = FALSE
        )
      }

      out[[i]] <- list(
        name = mark_name,
        mode = "joint",
        bounds = bounds,
        method = NULL,
        args = list(),
        values = values
      )
    } else {
      method <- spec$method

      if (is.null(method)) {
        method <- "idw"
      }

      if (!is.character(method) || length(method) != 1L ||
          !method %in% c("idw", "nearest")) {
        stop(
          sprintf(
            "Field continuous mark `%s` must use method 'idw' or 'nearest'.",
            mark_name
          ),
          call. = FALSE
        )
      }

      args <- spec$args

      if (is.null(args)) {
        args <- list()
      }

      if (!is.list(args)) {
        stop(
          sprintf(
            "`continuous_marks[[\"%s\"]]$args` must be a list.",
            mark_name
          ),
          call. = FALSE
        )
      }

      bad_args <- setdiff(names(args), c("p", "k", "axis_scale", "chunk_size"))

      if (length(bad_args) > 0L) {
        stop(
          sprintf(
            "Field continuous mark `%s` has unsupported argument(s): %s.",
            mark_name,
            paste(bad_args, collapse = ", ")
          ),
          call. = FALSE
        )
      }

      if (!is.null(spec$bounds)) {
        warning(
          sprintf("`bounds` are ignored for field continuous mark `%s`.", mark_name),
          call. = FALSE
        )
      }

      out[[i]] <- list(
        name = mark_name,
        mode = "field",
        bounds = NULL,
        method = method,
        args = args,
        values = values
      )
    }
  }

  out
}


.pp_fit_normalize_categorical_marks <- function(
    pattern,
    categorical_marks,
    categorical_levels = NULL
) {
  if (is.null(categorical_marks) ||
      (is.character(categorical_marks) && length(categorical_marks) == 0L)) {

    if (is.null(categorical_levels) ||
        (is.list(categorical_levels) && length(categorical_levels) == 0L)) {
      return(list())
    }

    stop(
      "`categorical_levels` can be supplied only when `categorical_marks` are supplied.",
      call. = FALSE
    )
  }

  if (!is.character(categorical_marks) ||
      length(categorical_marks) == 0L ||
      any(!nzchar(categorical_marks)) ||
      anyDuplicated(categorical_marks)) {
    stop(
      "`categorical_marks` must be NULL or a character vector of unique mark names.",
      call. = FALSE
    )
  }

  if (is.null(pattern$marks) || !is.data.frame(pattern$marks)) {
    stop(
      "`categorical_marks` were supplied but `pattern` has no marks.",
      call. = FALSE
    )
  }

  missing_marks <- setdiff(categorical_marks, names(pattern$marks))

  if (length(missing_marks) > 0L) {
    stop(
      sprintf(
        "Unknown categorical mark(s): %s.",
        paste(missing_marks, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!is.null(categorical_levels)) {
    if (!is.list(categorical_levels) ||
        is.null(names(categorical_levels)) ||
        any(!nzchar(names(categorical_levels))) ||
        anyDuplicated(names(categorical_levels))) {
      stop(
        "`categorical_levels` must be NULL or a named list.",
        call. = FALSE
      )
    }

    missing_levels <- setdiff(categorical_marks, names(categorical_levels))
    extra_levels <- setdiff(names(categorical_levels), categorical_marks)

    if (length(missing_levels) > 0L || length(extra_levels) > 0L) {
      stop(
        paste0(
          "`categorical_levels` must be named exactly with the categorical ",
          "mark names."
        ),
        call. = FALSE
      )
    }

    categorical_levels <- categorical_levels[categorical_marks]
  }

  out <- vector("list", length(categorical_marks))
  names(out) <- categorical_marks

  for (i in seq_along(categorical_marks)) {
    mark_name <- categorical_marks[i]
    values <- pattern$marks[[mark_name]]

    if (is.numeric(values)) {
      stop(
        sprintf(
          "Categorical mark `%s` is numeric. Declare it as a continuous mark instead.",
          mark_name
        ),
        call. = FALSE
      )
    }

    if (!(is.factor(values) || is.character(values))) {
      stop(
        sprintf(
          "Categorical mark `%s` must be a factor or character vector.",
          mark_name
        ),
        call. = FALSE
      )
    }

    if (anyNA(values)) {
      stop(
        sprintf("Categorical mark `%s` contains missing values.", mark_name),
        call. = FALSE
      )
    }

    values_character <- as.character(values)

    default_levels <- if (is.factor(values)) {
      as.character(levels(values))
    } else {
      sort(unique(values_character))
    }

    levels_i <- if (is.null(categorical_levels)) {
      default_levels
    } else {
      categorical_levels[[mark_name]]
    }

    if (is.factor(levels_i)) {
      levels_i <- as.character(levels_i)
    }

    if (!is.character(levels_i) ||
        length(levels_i) == 0L ||
        anyNA(levels_i) ||
        any(!nzchar(levels_i)) ||
        anyDuplicated(levels_i)) {
      stop(
        sprintf(
          "`categorical_levels[[\"%s\"]]` must contain unique non-empty character levels.",
          mark_name
        ),
        call. = FALSE
      )
    }

    missing_observed <- setdiff(unique(values_character), levels_i)

    if (length(missing_observed) > 0L) {
      stop(
        sprintf(
          paste0(
            "Observed values of categorical mark `%s` are absent from ",
            "`categorical_levels`: %s."
          ),
          mark_name,
          paste(missing_observed, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    out[[i]] <- list(
      name = mark_name,
      values = factor(values_character, levels = levels_i),
      levels = levels_i
    )
  }

  out
}


.pp_fit_prepare_joint_domain <- function(pattern, continuous_marks) {
  window_axes <- pattern$axis_info$name
  bounds <- pattern$window$bounds[, window_axes, drop = FALSE]

  axis_info <- pattern$axis_info
  rownames(axis_info) <- NULL

  joint_specs <- Filter(
    function(z) identical(z$mode, "joint"),
    continuous_marks
  )

  if (length(joint_specs) > 0L) {
    joint_names <- names(joint_specs)

    if (any(joint_names %in% window_axes)) {
      stop("Joint continuous-mark names cannot duplicate window-axis names.", call. = FALSE)
    }

    mark_bounds <- vapply(joint_specs, function(z) z$bounds, numeric(2L))
    mark_bounds <- as.matrix(mark_bounds)
    colnames(mark_bounds) <- joint_names
    rownames(mark_bounds) <- c("lower", "upper")

    bounds <- cbind(bounds, mark_bounds)

    axis_info <- rbind(
      axis_info,
      data.frame(
        name = joint_names,
        role = rep("mark", length(joint_names)),
        unit = rep(NA_character_, length(joint_names)),
        stringsAsFactors = FALSE
      )
    )
  }

  widths <- bounds["upper", ] - bounds["lower", ]
  hypervolume <- prod(as.numeric(widths))

  if (!is.finite(hypervolume) || hypervolume <= 0) {
    stop("The continuous joint domain must have positive finite hypervolume.", call. = FALSE)
  }

  list(
    bounds = bounds,
    axis_info = axis_info,
    window_axes = window_axes,
    joint_mark_names = names(joint_specs),
    dimension = ncol(bounds),
    hypervolume = hypervolume
  )
}


# Covariate sources ----------------------------------------------------------

.pp_fit_normalize_external_covariates <- function(covariates) {
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

  for (z in covariates) {
    if (!inherits(z, "pp_covariate")) {
      stop("Every external covariate must inherit from class 'pp_covariate'.", call. = FALSE)
    }

    .pp_cov_validate_object(z)
  }

  covariate_names <- names(covariates)

  if (is.null(covariate_names) || any(!nzchar(covariate_names))) {
    covariate_names <- vapply(covariates, function(z) z$name, character(1))
  }

  if (anyDuplicated(covariate_names)) {
    stop("External covariates must have unique names.", call. = FALSE)
  }

  names(covariates) <- covariate_names
  covariates
}


.pp_fit_validate_covariate_compatibility <- function(covariate, pattern) {
  pattern_axes <- pattern$axis_info$name
  cov_axes <- covariate$axis_info$name

  missing_axes <- setdiff(cov_axes, pattern_axes)

  if (length(missing_axes) > 0L) {
    stop(
      sprintf(
        "Covariate `%s` uses axis/axes absent from the pattern: %s.",
        covariate$name,
        paste(missing_axes, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  pattern_info <- pattern$axis_info[match(cov_axes, pattern_axes), , drop = FALSE]

  if (!identical(as.character(covariate$axis_info$role), as.character(pattern_info$role))) {
    stop(
      sprintf("Covariate `%s` has incompatible axis roles.", covariate$name),
      call. = FALSE
    )
  }

  if (!.pp_fit_same_units(covariate$axis_info$unit, pattern_info$unit)) {
    stop(
      sprintf("Covariate `%s` has incompatible axis units.", covariate$name),
      call. = FALSE
    )
  }

  if (!is.null(covariate$crs) &&
      !is.null(pattern$crs) &&
      !identical(covariate$crs, pattern$crs)) {
    stop(
      sprintf("Covariate `%s` has an incompatible CRS.", covariate$name),
      call. = FALSE
    )
  }

  cov_time <- covariate$time_info
  pattern_time <- pattern$time_info

  if (!is.null(cov_time) && !is.null(pattern_time)) {
    if (!identical(cov_time$axis, pattern_time$axis) ||
        !identical(cov_time$unit, pattern_time$unit) ||
        !identical(cov_time$origin, pattern_time$origin)) {
      stop(
        sprintf("Covariate `%s` has incompatible temporal metadata.", covariate$name),
        call. = FALSE
      )
    }
  }

  invisible(covariate)
}


.pp_fit_make_pattern_field <- function(pattern, values, name) {
  if (!is.numeric(values) || any(!is.finite(values))) {
    stop(
      sprintf("Pattern-based covariate `%s` must contain finite numeric values.", name),
      call. = FALSE
    )
  }

  axis_names <- pattern$axis_info$name
  field_data <- as.data.frame(pattern$coords[, axis_names, drop = FALSE], check.names = FALSE)
  field_data[[name]] <- as.numeric(values)

  pp_covariate(
    data = field_data,
    coords = axis_names,
    value = name,
    axis_roles = stats::setNames(pattern$axis_info$role, axis_names),
    axis_units = stats::setNames(pattern$axis_info$unit, axis_names),
    time_unit = if (is.null(pattern$time_info)) NULL else pattern$time_info$unit,
    time_origin = if (is.null(pattern$time_info)) NULL else pattern$time_info$origin,
    window = pattern$window,
    crs = pattern$crs,
    name = name
  )
}


.pp_fit_normalize_covariate_options <- function(
    active_names,
    methods,
    args,
    default_methods = NULL
) {
  if (length(active_names) == 0L) {
    if (!is.null(methods) || !is.null(args)) {
      stop(
        "Covariate options were supplied but no active covariates appear in the formula.",
        call. = FALSE
      )
    }

    return(list(methods = character(0), args = list()))
  }

  if (is.null(default_methods)) {
    default_methods <- stats::setNames(rep("idw", length(active_names)), active_names)
  }

  if (!is.character(default_methods) ||
      length(default_methods) != length(active_names) ||
      is.null(names(default_methods)) ||
      !identical(names(default_methods), active_names)) {
    stop("Internal error: invalid default covariate methods.", call. = FALSE)
  }

  if (is.null(methods)) {
    methods <- default_methods
  } else {
    if (!is.character(methods)) {
      stop("`covariate_methods` must be NULL or a character vector.", call. = FALSE)
    }

    if (!is.null(names(methods))) {
      missing_names <- setdiff(active_names, names(methods))
      extra_names <- setdiff(names(methods), active_names)

      if (length(missing_names) > 0L || length(extra_names) > 0L) {
        stop(
          "Named `covariate_methods` must match the active covariate names exactly.",
          call. = FALSE
        )
      }

      methods <- methods[active_names]
    } else if (length(methods) != length(active_names)) {
      stop("`covariate_methods` must have one entry per active covariate.", call. = FALSE)
    }
  }

  allowed_methods <- c("idw", "nearest", "distance_min")

  if (any(!methods %in% allowed_methods)) {
    stop(
      sprintf(
        "Covariate methods must be among: %s.",
        paste(allowed_methods, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  methods <- stats::setNames(as.character(methods), active_names)

  if (is.null(args)) {
    args <- list()
  }

  if (!is.list(args)) {
    stop("`covariate_args` must be NULL or a named list.", call. = FALSE)
  }

  if (length(args) > 0L) {
    if (is.null(names(args)) || any(!nzchar(names(args)))) {
      stop("`covariate_args` must be a named list.", call. = FALSE)
    }

    unknown <- setdiff(names(args), active_names)

    if (length(unknown) > 0L) {
      stop(
        sprintf(
          "`covariate_args` contains inactive or unknown covariate(s): %s.",
          paste(unknown, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  out_args <- vector("list", length(active_names))
  names(out_args) <- active_names

  allowed_args <- c("p", "k", "axis_scale", "chunk_size")

  for (name_i in active_names) {
    args_i <- args[[name_i]]

    if (is.null(args_i)) {
      args_i <- list()
    }

    if (!is.list(args_i)) {
      stop(sprintf("`covariate_args[[\"%s\"]]` must be a list.", name_i), call. = FALSE)
    }

    if (length(args_i) > 0L && (is.null(names(args_i)) || any(!nzchar(names(args_i))))) {
      stop(sprintf("`covariate_args[[\"%s\"]]` must be a named list.", name_i), call. = FALSE)
    }

    bad_args <- setdiff(names(args_i), allowed_args)

    if (length(bad_args) > 0L) {
      stop(
        sprintf(
          "`covariate_args[[\"%s\"]]` contains unsupported option(s): %s.",
          name_i,
          paste(bad_args, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    out_args[[name_i]] <- args_i
  }

  list(methods = methods, args = out_args)
}


.pp_fit_prepare_covariate_sources <- function(
    pattern,
    formula,
    covariates,
    continuous_marks,
    covariate_methods,
    covariate_args
) {
  formula_variables <- .pp_fit_formula_variables(formula)
  axis_names <- pattern$axis_info$name

  external <- .pp_fit_normalize_external_covariates(covariates)
  external_names <- names(external)

  if (length(external_names) > 0L) {
    unused_external <- setdiff(external_names, formula_variables)

    if (length(unused_external) > 0L) {
      stop(
        sprintf(
          paste0(
            "External covariate(s) must appear in `formula`: %s. ",
            "Remove them or include them explicitly."
          ),
          paste(unused_external, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    for (z in external) {
      .pp_fit_validate_covariate_compatibility(z, pattern)
    }
  }

  pattern_covariates <- pattern$covariates

  if (is.null(pattern_covariates)) {
    pattern_covariates <- data.frame(row.names = seq_len(nrow(pattern$coords)))
  }

  if (!is.data.frame(pattern_covariates)) {
    stop("`pattern$covariates` must be a data frame when present.", call. = FALSE)
  }

  if (nrow(pattern_covariates) != nrow(pattern$coords)) {
    stop("`pattern$covariates` must have one row per observed event.", call. = FALSE)
  }

  candidate_pattern_names <- intersect(names(pattern_covariates), formula_variables)
  candidate_field_marks <- names(Filter(function(z) identical(z$mode, "field"), continuous_marks))
  candidate_field_marks <- intersect(candidate_field_marks, formula_variables)

  overlap_external_pattern <- intersect(external_names, candidate_pattern_names)
  overlap_external_field <- intersect(external_names, candidate_field_marks)
  overlap_pattern_field <- intersect(candidate_pattern_names, candidate_field_marks)

  if (length(overlap_external_pattern) > 0L ||
      length(overlap_external_field) > 0L ||
      length(overlap_pattern_field) > 0L) {
    overlap <- unique(c(overlap_external_pattern, overlap_external_field, overlap_pattern_field))

    stop(
      sprintf(
        "A covariate name is supplied from more than one source: %s.",
        paste(overlap, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  option_names <- c(external_names, candidate_pattern_names)

  default_methods <- character(0)
  if (length(external_names) > 0L) {
    default_methods <- c(
      default_methods,
      vapply(
        external,
        function(z) if (.pp_cov_is_support(z)) "distance_min" else "idw",
        character(1)
      )
    )
  }
  if (length(candidate_pattern_names) > 0L) {
    default_methods <- c(
      default_methods,
      stats::setNames(rep("idw", length(candidate_pattern_names)), candidate_pattern_names)
    )
  }
  if (length(default_methods) > 0L) {
    default_methods <- default_methods[option_names]
  }

  options <- .pp_fit_normalize_covariate_options(
    active_names = option_names,
    methods = covariate_methods,
    args = covariate_args,
    default_methods = default_methods
  )

  active_names <- c(option_names, candidate_field_marks)
  sources <- vector("list", length(active_names))
  names(sources) <- active_names

  if (length(external_names) > 0L) {
    for (name_i in external_names) {
      method_i <- options$methods[[name_i]]

      if (.pp_cov_is_field(external[[name_i]]) && !method_i %in% c("idw", "nearest")) {
        stop(
          sprintf(
            "Field covariate `%s` must use method 'idw' or 'nearest'.",
            name_i
          ),
          call. = FALSE
        )
      }

      if (.pp_cov_is_support(external[[name_i]]) && !identical(method_i, "distance_min")) {
        stop(
          sprintf(
            "Geometric-support covariate `%s` must use method 'distance_min'.",
            name_i
          ),
          call. = FALSE
        )
      }

      sources[[name_i]] <- list(
        name = name_i,
        source = "external",
        object = external[[name_i]],
        method = method_i,
        args = options$args[[name_i]]
      )
    }
  }

  if (length(candidate_pattern_names) > 0L) {
    for (name_i in candidate_pattern_names) {
      object_i <- .pp_fit_make_pattern_field(
        pattern = pattern,
        values = pattern_covariates[[name_i]],
        name = name_i
      )

      method_i <- options$methods[[name_i]]

      if (!method_i %in% c("idw", "nearest")) {
        stop(
          sprintf(
            "Pattern-based covariate `%s` must use method 'idw' or 'nearest'.",
            name_i
          ),
          call. = FALSE
        )
      }

      sources[[name_i]] <- list(
        name = name_i,
        source = "pattern",
        object = object_i,
        method = method_i,
        args = options$args[[name_i]]
      )
    }
  }

  if (length(candidate_field_marks) > 0L) {
    for (name_i in candidate_field_marks) {
      spec_i <- continuous_marks[[name_i]]
      object_i <- .pp_fit_make_pattern_field(
        pattern = pattern,
        values = spec_i$values,
        name = name_i
      )

      sources[[name_i]] <- list(
        name = name_i,
        source = "continuous_mark_field",
        object = object_i,
        method = spec_i$method,
        args = spec_i$args
      )
    }
  }

  sources
}


.pp_fit_evaluate_sources <- function(sources, target_coords) {
  if (length(sources) == 0L) {
    return(data.frame(row.names = seq_len(nrow(target_coords))))
  }

  target <- as.data.frame(target_coords, check.names = FALSE)
  out <- data.frame(row.names = seq_len(nrow(target)))

  for (name_i in names(sources)) {
    source_i <- sources[[name_i]]

    eval_args <- c(
      list(
        covariate = source_i$object,
        at = target,
        method = source_i$method
      ),
      source_i$args
    )

    value_i <- do.call(evaluate_covariate, eval_args)

    if (!is.numeric(value_i) || length(value_i) != nrow(target) || any(!is.finite(value_i))) {
      stop(
        sprintf(
          "Evaluation of covariate `%s` did not produce finite numeric values.",
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


# Cubature construction ------------------------------------------------------

.pp_fit_resolve_nc <- function(domain, n_events, q, nc, max_tiles) {
  axis_names <- domain$axis_info$name
  d <- length(axis_names)

  max_tiles <- .pp_fit_validate_positive_integer(max_tiles, "max_tiles", min_value = 1L)

  if (is.null(nc)) {
    q <- .pp_fit_validate_positive_number(q, "q")
    target_tiles <- max(1L, ceiling(q * max(1L, n_events)))

    if (target_tiles > max_tiles) {
      stop(
        sprintf(
          paste0(
            "Automatic cubature requests approximately %d tiles, exceeding ",
            "`max_tiles = %d`. Reduce `q`, increase `max_tiles`, or supply `nc`."
          ),
          target_tiles,
          max_tiles
        ),
        call. = FALSE
      )
    }

    base <- max(1L, floor(target_tiles^(1 / d)))
    out <- rep(as.integer(base), d)
    names(out) <- axis_names

    while (prod(as.double(out)) < target_tiles) {
      candidate_products <- vapply(
        seq_len(d),
        function(j) prod(as.double(replace(out, j, out[j] + 1L))),
        numeric(1)
      )

      feasible <- which(candidate_products <= max_tiles)

      if (length(feasible) == 0L) {
        break
      }

      current_min <- min(out[feasible])
      eligible <- feasible[out[feasible] == current_min]
      chosen <- eligible[which.min(candidate_products[eligible])]
      out[chosen] <- out[chosen] + 1L
    }

    if (prod(as.double(out)) > max_tiles) {
      stop("Automatic cubature resolution exceeds `max_tiles`.", call. = FALSE)
    }

    return(out)
  }

  if (!is.numeric(nc) || any(!is.finite(nc)) || any(nc != as.integer(nc)) || any(nc < 1L)) {
    stop("`nc` must contain positive integer numbers of tiles.", call. = FALSE)
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
          "Named `nc` must match the joint-domain axes exactly: %s.",
          paste(axis_names, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    out <- as.integer(nc[axis_names])
  } else if (length(nc) == 1L) {
    out <- rep(as.integer(nc), d)
  } else if (length(nc) == d) {
    out <- as.integer(nc)
  } else {
    stop(
      sprintf("`nc` must be scalar, length %d, or named by joint-domain axes.", d),
      call. = FALSE
    )
  }

  names(out) <- axis_names

  n_tiles <- prod(as.double(out))

  if (n_tiles > max_tiles) {
    stop(
      sprintf("`nc` defines %s tiles, exceeding `max_tiles = %s`.", n_tiles, max_tiles),
      call. = FALSE
    )
  }

  out
}


.pp_fit_tile_index_grid <- function(nc) {
  index_list <- lapply(nc, seq_len)
  names(index_list) <- names(nc)

  out <- do.call(
    expand.grid,
    c(index_list, list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  )

  out <- as.matrix(out)
  storage.mode(out) <- "integer"
  colnames(out) <- names(nc)
  out
}


.pp_fit_generate_dummy_points <- function(domain, nc, dummy_layout) {
  dummy_layout <- match.arg(dummy_layout, c("grid", "random"))

  index_grid <- .pp_fit_tile_index_grid(nc)
  axis_names <- colnames(index_grid)
  n_tiles <- nrow(index_grid)

  lower <- domain$bounds["lower", axis_names]
  upper <- domain$bounds["upper", axis_names]
  widths <- upper - lower

  coords <- matrix(
    NA_real_,
    nrow = n_tiles,
    ncol = length(axis_names),
    dimnames = list(NULL, axis_names)
  )

  for (j in seq_along(axis_names)) {
    if (identical(dummy_layout, "grid")) {
      fraction <- index_grid[, j] - 0.5
    } else {
      fraction <- index_grid[, j] - 1 + stats::runif(n_tiles)
    }

    coords[, j] <- lower[j] + fraction * widths[j] / nc[j]
  }

  list(coords = coords, tile_index = index_grid)
}


.pp_fit_assign_tiles <- function(coords, domain, nc) {
  axis_names <- names(nc)
  lower <- domain$bounds["lower", axis_names]
  upper <- domain$bounds["upper", axis_names]
  widths <- upper - lower

  if (nrow(coords) == 0L) {
    return(integer(0))
  }

  indices <- matrix(
    NA_integer_,
    nrow = nrow(coords),
    ncol = length(axis_names)
  )

  for (j in seq_along(axis_names)) {
    value <- coords[, axis_names[j]]
    scaled <- floor((value - lower[j]) / widths[j] * nc[j]) + 1L
    scaled <- pmax(1L, pmin(nc[j], scaled))
    indices[, j] <- as.integer(scaled)
  }

  multipliers <- c(1L, cumprod(as.integer(nc[-length(nc)])))
  as.integer(1L + (indices - 1L) %*% multipliers)
}


.pp_fit_build_base_cubature <- function(pattern, domain, nc, dummy_layout) {
  axis_names <- domain$axis_info$name
  window_axes <- domain$window_axes
  joint_mark_names <- domain$joint_mark_names
  n_events <- nrow(pattern$coords)

  event_coords <- matrix(
    NA_real_,
    nrow = n_events,
    ncol = length(axis_names),
    dimnames = list(NULL, axis_names)
  )

  if (n_events > 0L) {
    event_coords[, window_axes] <- pattern$coords[, window_axes, drop = FALSE]

    if (length(joint_mark_names) > 0L) {
      for (name_i in joint_mark_names) {
        event_coords[, name_i] <- as.numeric(pattern$marks[[name_i]])
      }
    }
  }

  dummy <- .pp_fit_generate_dummy_points(
    domain = domain,
    nc = nc,
    dummy_layout = dummy_layout
  )

  coords <- rbind(event_coords, dummy$coords)
  event <- c(rep.int(1L, n_events), rep.int(0L, nrow(dummy$coords)))
  observed_id <- c(seq_len(n_events), rep.int(NA_integer_, nrow(dummy$coords)))

  tile_id <- .pp_fit_assign_tiles(
    coords = coords,
    domain = domain,
    nc = nc
  )

  n_tiles <- prod(as.double(nc))
  tile_counts <- tabulate(tile_id, nbins = n_tiles)
  tile_volume <- domain$hypervolume / n_tiles
  weight <- tile_volume / tile_counts[tile_id]

  if (any(!is.finite(weight)) || any(weight <= 0)) {
    stop("Failed to construct strictly positive cubature weights.", call. = FALSE)
  }

  out <- as.data.frame(coords, check.names = FALSE)
  out$.pp_event <- event
  out$.pp_weight <- as.numeric(weight)
  out$.pp_tile <- tile_id
  out$.pp_observed_id <- observed_id
  out$.pp_origin <- factor(ifelse(event == 1L, "event", "dummy"), levels = c("event", "dummy"))

  list(
    data = out,
    n_events = n_events,
    n_dummy = nrow(dummy$coords),
    n_tiles = n_tiles,
    tile_volume = tile_volume
  )
}


.pp_fit_replicate_categorical_marks <- function(
    base,
    categorical_marks,
    max_mark_combinations,
    max_quadrature_points
) {
  if (length(categorical_marks) == 0L) {
    return(base)
  }

  max_mark_combinations <- .pp_fit_validate_positive_integer(
    max_mark_combinations,
    "max_mark_combinations",
    min_value = 1L
  )

  level_list <- lapply(categorical_marks, function(z) z$levels)
  names(level_list) <- names(categorical_marks)

  combinations <- do.call(
    expand.grid,
    c(level_list, list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  )

  n_combinations <- nrow(combinations)

  if (n_combinations > max_mark_combinations) {
    stop(
      sprintf(
        paste0(
          "Categorical marks define %d level combinations, exceeding ",
          "`max_mark_combinations = %d`."
        ),
        n_combinations,
        max_mark_combinations
      ),
      call. = FALSE
    )
  }

  n_base <- nrow(base$data)

  if (n_base * n_combinations > max_quadrature_points) {
    stop(
      sprintf(
        paste0(
          "The cubature data would contain %d rows after categorical replication, ",
          "exceeding `max_quadrature_points = %d`."
        ),
        n_base * n_combinations,
        max_quadrature_points
      ),
      call. = FALSE
    )
  }

  index <- rep(seq_len(n_base), each = n_combinations)
  expanded <- base$data[index, , drop = FALSE]
  rownames(expanded) <- NULL

  for (name_i in names(categorical_marks)) {
    levels_i <- categorical_marks[[name_i]]$levels
    expanded[[name_i]] <- factor(
      rep(as.character(combinations[[name_i]]), times = n_base),
      levels = levels_i
    )
  }

  event_indicator <- expanded$.pp_event == 1L
  observed_id <- expanded$.pp_observed_id

  for (name_i in names(categorical_marks)) {
    observed_values <- as.character(categorical_marks[[name_i]]$values)
    target_values <- as.character(expanded[[name_i]])

    matches <- rep(FALSE, nrow(expanded))
    valid <- !is.na(observed_id)
    matches[valid] <- target_values[valid] == observed_values[observed_id[valid]]
    event_indicator <- event_indicator & matches
  }

  expanded$.pp_event <- as.integer(event_indicator)

  base$data <- expanded
  base$n_mark_combinations <- n_combinations
  base$categorical_levels <- lapply(categorical_marks, function(z) z$levels)
  base
}

# -------------------------------------------------------------------------
# Event weights for fractional-count Poisson cubature
# -------------------------------------------------------------------------

.pp_fit_validate_event_weights <- function(pattern, event_weights = NULL) {
  n_events <- nrow(pattern$coords)

  if (is.null(event_weights)) {
    return(rep.int(1, n_events))
  }

  if (!is.numeric(event_weights) ||
      is.matrix(event_weights) ||
      is.data.frame(event_weights) ||
      length(event_weights) != n_events) {
    stop(
      sprintf(
        "`event_weights` must be a numeric vector of length %d.",
        n_events
      ),
      call. = FALSE
    )
  }

  event_weights <- as.numeric(event_weights)

  if (any(!is.finite(event_weights))) {
    stop("`event_weights` must contain only finite values.", call. = FALSE)
  }

  if (any(event_weights < 0)) {
    stop("`event_weights` must be non-negative.", call. = FALSE)
  }

  event_weights
}


#' Construct a multidimensional Poisson cubature scheme
#'
#' Constructs the event--dummy cubature data required to approximate the
#' likelihood of a multidimensional Poisson point-process model. The continuous
#' integration domain is the observation window, optionally augmented by
#' continuous marks treated as joint-domain coordinates. Categorical marks are
#' handled by replicating the continuous cubature scheme over their admissible
#' level combinations.
#'
#' @param pattern A point pattern created by [pp_pattern()]. It must contain a
#'   finite observation window created by [pp_window()]. Its coordinates, marks,
#'   covariates, axis metadata, and window define the observed component of the
#'   cubature scheme.
#' @param formula A one-sided model formula, such as `~ x + y` or
#'   `~ s(x, y) + cov1`. Every declared external covariate and mark must
#'   occur in the formula. Formula variables may be coordinate names, declared
#'   marks, external covariates, or numeric columns stored in
#'   `pattern$covariates`.
#' @param covariates `NULL`, one object created by [pp_covariate()], or a named
#'   list of such objects. Each supplied covariate must occur in `formula` and
#'   must use axes, roles, units, temporal metadata, and a coordinate reference
#'   system compatible with `pattern`. Unnamed lists are named from the
#'   covariate objects.
#' @param covariate_methods Optional character vector selecting the evaluation
#'   method for each active external or pattern-stored covariate. It may be
#'   named exactly with the active covariate names or contain one entry per
#'   active covariate. Allowed values are `"idw"`, `"nearest"`, and
#'   `"distance_min"`. Continuous fields default to `"idw"`; geometric
#'   supports use `"distance_min"`. See [evaluate_covariate()].
#' @param covariate_args Optional named list with one sub-list for each active
#'   external or pattern-stored covariate requiring non-default evaluation
#'   options. Supported entries are `p`, `k`, `axis_scale`, and `chunk_size`, as
#'   described in [evaluate_covariate()]. Continuous marks declared as fields
#'   take their method and arguments from `continuous_marks` instead.
#' @param continuous_marks `NULL` or a named list describing numeric marks in
#'   `pattern$marks`. Each sub-list must contain `mode = "joint"` or
#'   `mode = "field"`:
#'
#'   * a joint mark becomes an additional continuous integration coordinate and
#'     should provide finite `bounds = c(lower, upper)`. If bounds are omitted,
#'     the observed range is used with a warning; `method` and `args` are
#'     ignored;
#'   * a field mark is interpolated from its observed event values and may
#'     specify `method = "idw"` or `"nearest"` and an `args` list containing
#'     `p`, `k`, `axis_scale`, or `chunk_size`. Its default method is `"idw"`,
#'     and `bounds` are ignored.
#'
#'   Every declared continuous mark must occur in `formula`.
#' @param categorical_marks `NULL` or a character vector naming unique factor
#'   or character marks in `pattern$marks`. The cubature data are replicated
#'   over the Cartesian product of their levels, and only the row matching an
#'   observed event's mark combination receives an event contribution. Every
#'   declared categorical mark must occur in `formula`.
#' @param event_weights `NULL` or a non-negative finite numeric vector with one
#'   value per observed event. The default assigns unit mass to every event.
#'   Supplied values replace the unit event contribution and permit
#'   fractional-count objectives, such as posterior background probabilities in
#'   EM or ETAS algorithms. Zero weights are allowed.
#' @param q A strictly positive number giving the target ratio of continuous
#'   cubature tiles to `max(1, n)`, where `n` is the number of observed events,
#'   when `nc` is `NULL`. The resulting balanced integer resolution is chosen
#'   automatically and is constrained by `max_tiles`.
#' @param dummy_layout Character string selecting the dummy-point position in
#'   each regular hyperrectangular tile. `"grid"` places the dummy at the tile
#'   centre; `"random"` draws one point uniformly inside each tile.
#' @param nc Optional cubature resolution for the continuous joint domain. It
#'   may be one positive integer recycled over all continuous axes, one value
#'   per continuous axis, or a named vector matching exactly the coordinate and
#'   joint-mark axes. When supplied, `nc` overrides the automatic resolution
#'   based on `q`.
#' @param categorical_levels Optional named list defining the full admissible
#'   levels of every mark in `categorical_marks`. Its names must match
#'   `categorical_marks` exactly, and all observed values must be included. If
#'   omitted, factor levels are retained and character values are converted to
#'   sorted unique levels. Explicit levels allow unobserved but admissible
#'   categories to remain in the discrete mark domain.
#' @param max_tiles Positive integer giving the maximum allowed number of
#'   continuous-domain tiles, equal to the product of the resolved entries of
#'   `nc`. It protects against excessive multidimensional grids.
#' @param max_quadrature_points Positive integer giving the maximum number of
#'   cubature rows after replication over categorical-mark combinations.
#' @param max_mark_combinations Positive integer giving the maximum allowed
#'   Cartesian product of categorical-mark levels.
#' @param seed `NULL` or one finite integer used to initialize R's random-number
#'   generator when `dummy_layout = "random"`. It has no effect for centred grid
#'   dummies.
#'
#' @details
#' Let \eqn{D} denote the continuous model domain. For an unmarked process,
#' \eqn{D=W}, where \eqn{W} is the observation window. If continuous marks
#' \eqn{m_1,\ldots,m_J} are declared with `mode = "joint"`, then
#'
#' \deqn{D = W \mathbin{\times} M_1 \mathbin{\times} \cdots
#'       \mathbin{\times} M_J,}
#'
#' where each \eqn{M_j} is the interval specified by its bounds. Categorical
#' marks are integrated with respect to counting measure by summing over their
#' level combinations.
#'
#' For log-intensity \eqn{\eta(u)=\log\lambda(u)} and optional event masses
#' \eqn{v_i}, the weighted Poisson objective is
#'
#' \deqn{\ell(\eta) = \sum_{i=1}^{n} v_i\eta(u_i)
#'       - \int_D \exp\{\eta(u)\}\,\mathrm{d}u.}
#'
#' The cubature scheme approximates the integral by
#'
#' \deqn{\int_D \lambda(u)\,\mathrm{d}u
#'       \simeq \sum_{j=1}^{J} a_j\lambda(u_j),}
#'
#' where the locations \eqn{u_j} include observed events and one dummy point per
#' tile. Within a tile, its hypervolume is divided equally among all event and
#' dummy rows assigned to that tile, so the corresponding weights \eqn{a_j}
#' sum to the tile hypervolume. The pseudo-response stored in `.pp_response` is
#' the event mass divided by the cubature weight. This representation allows
#' [fit_pp()] to fit the approximation as a weighted Poisson GAM with log link.
#'
#' A joint continuous mark is part of the modelled point-process domain and is
#' integrated in the likelihood. A field mark is instead treated as a numeric
#' field reconstructed from the observed events and evaluated at dummy
#' locations. These two modes therefore have different statistical meanings.
#'
#' The number of tiles grows as the product of the axis resolutions, and
#' categorical replication multiplies the number of cubature rows further.
#' `max_tiles`, `max_quadrature_points`, and `max_mark_combinations` are explicit
#' safeguards against this dimensional growth.
#'
#' @section Cubature object:
#' The returned `pp_cubature` object contains:
#'
#' * `data`: the complete event--dummy data frame. Internal columns include
#'   `.pp_event`, `.pp_weight`, `.pp_tile`, `.pp_observed_id`, `.pp_origin`,
#'   `.pp_event_weight`, and `.pp_response`;
#' * `formula` and `pattern`: the requested model formula and original pattern;
#' * `domain`: continuous bounds, axis metadata, dimension, and hypervolume;
#' * `continuous_marks` and `categorical_marks`: normalized mark specifications;
#' * `covariates`: the retained source, evaluation method, and arguments for
#'   each active covariate;
#' * `settings`: the resolved cubature controls;
#' * `summary`: event, dummy, tile, mark-combination, row, weight, and
#'   hypervolume summaries.
#'
#' @return An object of class `"pp_cubature"`.
#'
#' @examples
#' w <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#'
#' set.seed(1)
#' x <- pp_pattern(
#'   data = data.frame(x = runif(40), y = runif(40)),
#'   coords = c("x", "y"),
#'   window = w
#' )
#'
#' # Coordinate-only cubature
#' cub <- pp_cubature(
#'   pattern = x,
#'   formula = ~ x + y,
#'   q = 2
#' )
#' cub
#' summary(cub)
#'
#' # Include an externally observed covariate field
#' field_data <- expand.grid(
#'   x = seq(0, 1, length.out = 5),
#'   y = seq(0, 1, length.out = 5)
#' )
#' field_data$cov1 <- with(field_data, 15 + 3 * x - 2 * y)
#' cov1 <- pp_covariate(
#'   data = field_data,
#'   coords = c("x", "y"),
#'   value = "cov1",
#'   window = w,
#'   name = "cov1"
#' )
#'
#' cub_cov <- pp_cubature(
#'   pattern = x,
#'   formula = ~ x + cov1,
#'   covariates = list(cov1 = cov1),
#'   nc = c(x = 6, y = 6)
#' )
#'
#' # A joint continuous mark and a categorical mark
#' set.seed(2)
#' xm <- pp_pattern(
#'   data = data.frame(
#'     x = runif(30),
#'     y = runif(30),
#'     magnitude = runif(30, 0, 1),
#'     type = factor(sample(c("A", "B"), 30, replace = TRUE))
#'   ),
#'   coords = c("x", "y"),
#'   marks = c("magnitude", "type"),
#'   window = w
#' )
#'
#' cub_marked <- pp_cubature(
#'   pattern = xm,
#'   formula = ~ x + magnitude + type,
#'   continuous_marks = list(
#'     magnitude = list(mode = "joint", bounds = c(0, 1))
#'   ),
#'   categorical_marks = "type",
#'   nc = c(x = 4, y = 4, magnitude = 3)
#' )
#'
#' @references
#' Berman, M. and Turner, T. R. (1992). Approximating point process
#' likelihoods with GLIM. *Journal of the Royal Statistical Society: Series C
#' (Applied Statistics)*, **41**(1), 31--38.
#' \doi{10.2307/2347614}
#'
#' Baddeley, A. and Turner, R. (2000). Practical maximum pseudolikelihood for
#' spatial point patterns. *Australian & New Zealand Journal of Statistics*,
#' **42**(3), 283--322. \doi{10.1111/1467-842X.00128}
#'
#' @seealso
#' [fit_pp()] fits a Poisson point-process model from the cubature scheme;
#' [pp_pattern()] and [pp_window()] construct the observed pattern and domain;
#' [pp_covariate()] and [evaluate_covariate()] define and evaluate covariates.
#'
#' @family Poisson point-process models
#' @export
pp_cubature <- function(
    pattern,
    formula = ~ 1,
    covariates = NULL,
    covariate_methods = NULL,
    covariate_args = NULL,
    continuous_marks = NULL,
    categorical_marks = NULL,
    event_weights = NULL,
    q = 3,
    dummy_layout = c("grid", "random"),
    nc = NULL,
    categorical_levels = NULL,
    max_tiles = 1e6,
    max_quadrature_points = 1e6,
    max_mark_combinations = 100L,
    seed = NULL
) {
  .pp_fit_validate_pattern(pattern)
  formula <- .pp_fit_validate_formula(formula)
  dummy_layout <- match.arg(dummy_layout)

  event_weights_supplied <- !is.null(event_weights)

  event_weights <- .pp_fit_validate_event_weights(
    pattern = pattern,
    event_weights = event_weights
  )

  max_quadrature_points <- .pp_fit_validate_positive_integer(
    max_quadrature_points,
    "max_quadrature_points",
    min_value = 1L
  )

  seed <- .pp_fit_validate_seed(seed)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  continuous_marks <- .pp_fit_normalize_continuous_marks(
    pattern = pattern,
    continuous_marks = continuous_marks
  )

  categorical_marks <- .pp_fit_normalize_categorical_marks(
    pattern = pattern,
    categorical_marks = categorical_marks,
    categorical_levels = categorical_levels
  )

  overlap_marks <- intersect(names(continuous_marks), names(categorical_marks))

  if (length(overlap_marks) > 0L) {
    stop(
      sprintf(
        "A mark cannot be both continuous and categorical: %s.",
        paste(overlap_marks, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  formula_variables <- .pp_fit_formula_variables(formula)
  unused_continuous <- setdiff(names(continuous_marks), formula_variables)
  unused_categorical <- setdiff(names(categorical_marks), formula_variables)

  if (length(unused_continuous) > 0L || length(unused_categorical) > 0L) {
    unused_marks <- c(unused_continuous, unused_categorical)

    stop(
      sprintf(
        "Declared mark(s) must appear in `formula`: %s.",
        paste(unused_marks, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  domain <- .pp_fit_prepare_joint_domain(
    pattern = pattern,
    continuous_marks = continuous_marks
  )

  nc_resolved <- .pp_fit_resolve_nc(
    domain = domain,
    n_events = nrow(pattern$coords),
    q = q,
    nc = nc,
    max_tiles = max_tiles
  )

  base <- .pp_fit_build_base_cubature(
    pattern = pattern,
    domain = domain,
    nc = nc_resolved,
    dummy_layout = dummy_layout
  )

  sources <- .pp_fit_prepare_covariate_sources(
    pattern = pattern,
    formula = formula,
    covariates = covariates,
    continuous_marks = continuous_marks,
    covariate_methods = covariate_methods,
    covariate_args = covariate_args
  )

  covariate_values <- .pp_fit_evaluate_sources(
    sources = sources,
    target_coords = base$data[, domain$window_axes, drop = FALSE]
  )

  if (ncol(covariate_values) > 0L) {
    base$data <- cbind(base$data, covariate_values)
  }

  base <- .pp_fit_replicate_categorical_marks(
    base = base,
    categorical_marks = categorical_marks,
    max_mark_combinations = max_mark_combinations,
    max_quadrature_points = max_quadrature_points
  )

  available_variables <- names(base$data)

  unknown_formula_variables <- setdiff(
    .pp_fit_formula_variables(formula),
    available_variables
  )

  if (length(unknown_formula_variables) > 0L) {
    stop(
      sprintf(
        paste0(
          "`formula` refers to variable(s) unavailable in the cubature data: %s. ",
          "Provide them as coordinates, declared marks, external pp_covariates, ",
          "or columns of `pattern$covariates`."
        ),
        paste(unknown_formula_variables, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # -----------------------------------------------------------------------
  # Fractional event contribution.
  #
  # For ordinary fits, .pp_event_weight equals 1 for observed event rows
  # and 0 for dummies. With event_weights, it equals the corresponding
  # supplied weight, e.g. ETAS posterior background probability rho_i.
  # -----------------------------------------------------------------------

  base$data$.pp_event_weight <- 0

  event_rows <- base$data$.pp_event == 1L
  observed_id <- base$data$.pp_observed_id[event_rows]

  if (length(observed_id) > 0L) {
    if (anyNA(observed_id) ||
        any(observed_id < 1L) ||
        any(observed_id > length(event_weights))) {
      stop(
        "Internal error: invalid observed-event identifiers in cubature data.",
        call. = FALSE
      )
    }

    base$data$.pp_event_weight[event_rows] <-
      event_weights[observed_id]
  }

  base$data$.pp_response <-
    base$data$.pp_event_weight / base$data$.pp_weight

  if (any(!is.finite(base$data$.pp_response)) ||
      any(base$data$.pp_response < 0)) {
    stop("Failed to construct a valid cubature response.", call. = FALSE)
  }

  covariate_summary <- lapply(
    sources,
    function(z) {
      list(
        source = z$source,
        method = z$method,
        args = z$args,
        support_type = if (.pp_cov_is_field(z$object)) {
          "field"
        } else {
          "support"
        }
      )
    }
  )

  out <- list(
    data = base$data,
    formula = formula,
    pattern = pattern,
    domain = domain,

    continuous_marks = lapply(
      continuous_marks,
      function(z) z[c("name", "mode", "bounds", "method", "args")]
    ),

    categorical_marks = lapply(
      categorical_marks,
      function(z) z$levels
    ),

    covariates = covariate_summary,

    settings = list(
      q = q,
      dummy_layout = dummy_layout,
      nc = nc_resolved,
      max_tiles = max_tiles,
      max_quadrature_points = max_quadrature_points,
      max_mark_combinations = max_mark_combinations,
      seed = seed,
      event_weights_supplied = event_weights_supplied
    ),

    summary = list(
      n_events = nrow(pattern$coords),
      n_dummy = base$n_dummy,
      n_tiles = base$n_tiles,

      n_mark_combinations = if (is.null(base$n_mark_combinations)) {
        1L
      } else {
        base$n_mark_combinations
      },

      n_rows = nrow(base$data),

      continuous_dimension = domain$dimension,
      continuous_hypervolume = domain$hypervolume,

      sum_weights = sum(base$data$.pp_weight),
      weight_range = range(base$data$.pp_weight),

      n_event_rows = sum(base$data$.pp_event),
      n_positive_event_rows = sum(base$data$.pp_event_weight > 0),

      event_weights_supplied = event_weights_supplied,
      total_event_weight = sum(event_weights),
      event_weight_range = range(event_weights),
      n_positive_event_weights = sum(event_weights > 0)
    )
  )

  class(out) <- "pp_cubature"
  out
}


#' Methods for Poisson cubature objects
#'
#' Print and summarise the multidimensional cubature schemes returned by
#' [pp_cubature()].
#'
#' @details
#' `print.pp_cubature()` reports the number of observed events, dummy points,
#' continuous tiles, categorical combinations, final quadrature rows,
#' continuous-domain dimension and hypervolume, and the total cubature weight.
#'
#' `summary.pp_cubature()` returns the model formula, joint-domain axes,
#' cubature settings and numerical summary, covariate specifications, and mark
#' specifications. `print.summary_pp_cubature()` formats this information.
#'
#' @param x A `pp_cubature` object for `print.pp_cubature()`, or a
#'   `summary_pp_cubature` object for `print.summary_pp_cubature()`.
#' @param object A `pp_cubature` object to be summarised.
#' @param digits Number of significant digits used when printing numerical
#'   quantities.
#' @param ... Additional arguments. They are currently ignored.
#'
#' @return `print.pp_cubature()` and `print.summary_pp_cubature()` invisibly
#' return `x`. `summary.pp_cubature()` returns an object of class
#' `"summary_pp_cubature"`.
#'
#' @examples
#' w <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' set.seed(3)
#' x <- pp_pattern(
#'   data = data.frame(x = runif(20), y = runif(20)),
#'   coords = c("x", "y"),
#'   window = w
#' )
#' cub <- pp_cubature(x, formula = ~ x + y, nc = 4)
#'
#' print(cub)
#' sc <- summary(cub)
#' print(sc)
#'
#' @seealso [pp_cubature()] and [fit_pp()].
#' @family Poisson point-process models
#' @rdname pp_cubature-methods
#' @method print pp_cubature
#' @export
print.pp_cubature <- function(x, digits = 4, ...) {
  if (!inherits(x, "pp_cubature")) {
    stop("`x` must inherit from class 'pp_cubature'.", call. = FALSE)
  }

  cat("<pp_cubature>\n")
  cat("Multidimensional Poisson cubature scheme\n\n")
  cat("Observed events: ", x$summary$n_events, "\n", sep = "")
  cat("Dummy points: ", x$summary$n_dummy, "\n", sep = "")
  cat("Continuous tiles: ", x$summary$n_tiles, "\n", sep = "")
  cat("Categorical combinations: ", x$summary$n_mark_combinations, "\n", sep = "")
  cat("Quadrature rows: ", x$summary$n_rows, "\n", sep = "")
  cat("Continuous dimensions: ", x$summary$continuous_dimension, "\n", sep = "")
  cat(
    "Continuous hypervolume: ",
    format(signif(x$summary$continuous_hypervolume, digits), trim = TRUE),
    "\n",
    sep = ""
  )
  cat(
    "Sum of weights: ",
    format(signif(x$summary$sum_weights, digits), trim = TRUE),
    "\n",
    sep = ""
  )

  invisible(x)
}


#' @rdname pp_cubature-methods
#' @method summary pp_cubature
#' @export
summary.pp_cubature <- function(object, ...) {
  if (!inherits(object, "pp_cubature")) {
    stop("`object` must inherit from class 'pp_cubature'.", call. = FALSE)
  }

  out <- list(
    formula = object$formula,
    domain_axes = object$domain$axis_info,
    settings = object$settings,
    cubature = object$summary,
    covariates = object$covariates,
    continuous_marks = object$continuous_marks,
    categorical_marks = object$categorical_marks
  )

  class(out) <- "summary_pp_cubature"
  out
}


#' @rdname pp_cubature-methods
#' @method print summary_pp_cubature
#' @export
print.summary_pp_cubature <- function(x, digits = 4, ...) {
  cat("Cubature summary\n\n")
  cat("Formula: ")
  print(x$formula)
  cat("\n")

  cat("Joint continuous-domain axes:\n")
  domain_axes <- x$domain_axes
  if (!is.null(domain_axes$unit)) {
    domain_axes$unit[is.na(domain_axes$unit)] <- ""
  }
  print(domain_axes, row.names = FALSE)
  cat("\n")

  cat("Cubature:\n")
  print(as.data.frame(x$cubature), row.names = FALSE)

  if (length(x$covariates) > 0L) {
    cat("\nCovariates:\n")
    cov_names <- names(x$covariates)
    cov_source <- vapply(x$covariates, function(z) z$source, character(1))
    cov_method <- vapply(x$covariates, function(z) z$method, character(1))
    print(data.frame(name = cov_names, source = cov_source, method = cov_method), row.names = FALSE)
  }

  if (length(x$categorical_marks) > 0L) {
    cat("\nCategorical-mark levels:\n")
    for (name_i in names(x$categorical_marks)) {
      cat("  ", name_i, ": ", paste(x$categorical_marks[[name_i]], collapse = ", "), "\n", sep = "")
    }
  }

  invisible(x)
}


# Model fitting --------------------------------------------------------------

# Poisson family for Berman--Turner cubature responses ------------------------
#
# In the weighted Poisson representation of the point-process likelihood the
# working responses are y = z / w, which are generally non-integer. The
# estimating equations of the Poisson GLM/GAM are unaffected, but the `aic()`
# and `ls()` components of stats::poisson() evaluate stats::dpois(), which
# warns for every non-integer response. This family is identical to
# stats::poisson(link = "log") except that those two components use the
# continuous extension log p(y; mu) = y log(mu) - mu - lgamma(y + 1), which
# coincides with the Poisson log-probability for integer y.

.pp_log_poisson_density <- function(y, mu) {
  ifelse(y > 0, y * log(mu), 0) - mu - lgamma(y + 1)
}

.pp_poisson_family <- function() {
  family <- stats::poisson(link = "log")

  family$aic <- function(y, n, mu, wt, dev) {
    -2 * sum(.pp_log_poisson_density(y, mu) * wt)
  }

  family$ls <- function(y, w, n, scale) {
    c(sum(.pp_log_poisson_density(y, y) * w), 0, 0)
  }

  family
}


.pp_fit_canonical_continuous_marks <- function(continuous_marks) {
  if (length(continuous_marks) == 0L) {
    return(list())
  }

  out <- lapply(
    continuous_marks,
    function(z) {
      z[c("mode", "bounds", "method", "args")]
    }
  )

  names(out) <- names(continuous_marks)
  out
}


.pp_fit_make_model_spec <- function(
    pattern,
    formula,
    sources,
    continuous_marks,
    categorical_marks,
    q,
    dummy_layout,
    nc_requested,
    nc_resolved,
    max_tiles,
    max_quadrature_points,
    max_mark_combinations,
    seed,
    gam_args,
    keep_cubature
) {
  source_spec <- lapply(
    sources,
    function(z) {
      list(
        name = z$name,
        source = z$source,
        role = if (identical(z$source, "continuous_mark_field")) {
          "field_mark"
        } else {
          "covariate"
        },
        object = z$object,
        method = z$method,
        args = z$args
      )
    }
  )

  names(source_spec) <- names(sources)

  continuous_spec <- .pp_fit_canonical_continuous_marks(
    continuous_marks
  )

  categorical_spec <- lapply(
    categorical_marks,
    function(z) as.character(z$levels)
  )

  joint_mark_names <- names(
    Filter(
      function(z) identical(z$mode, "joint"),
      continuous_marks
    )
  )

  field_mark_names <- names(
    Filter(
      function(z) identical(z$mode, "field"),
      continuous_marks
    )
  )

  list(
    version = "1.0.0",
    formula = formula,
    formula_variables = .pp_fit_formula_variables(formula),

    window = pattern$window,
    axis_info = pattern$axis_info,
    time_info = pattern$time_info,
    crs = pattern$crs,

    covariate_sources = source_spec,

    continuous_marks = continuous_spec,
    categorical_levels = categorical_spec,

    mark_structure = list(
      joint_continuous = joint_mark_names,
      field_continuous = field_mark_names,
      categorical = names(categorical_marks)
    ),

    quadrature = list(
      q = q,
      dummy_layout = dummy_layout,
      nc_requested = nc_requested,
      nc_resolved = nc_resolved,
      max_tiles = max_tiles,
      max_quadrature_points = max_quadrature_points,
      max_mark_combinations = max_mark_combinations,
      seed = seed
    ),

    gam_args = gam_args,
    keep_cubature = keep_cubature
  )
}

.pp_fit_normalize_gam_args <- function(gam_args) {
  if (is.null(gam_args)) {
    gam_args <- list()
  }

  if (!is.list(gam_args)) {
    stop("`gam_args` must be NULL or a named list.", call. = FALSE)
  }

  if (length(gam_args) > 0L && (is.null(names(gam_args)) || any(!nzchar(names(gam_args))))) {
    stop("`gam_args` must be a named list.", call. = FALSE)
  }

  reserved <- c("formula", "family", "data", "weights")
  conflict <- intersect(names(gam_args), reserved)

  if (length(conflict) > 0L) {
    stop(
      sprintf(
        "`gam_args` cannot override internally defined argument(s): %s.",
        paste(conflict, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (is.null(gam_args$method)) {
    gam_args$method <- "REML"
  }

  gam_args
}


#' Fit a multidimensional Poisson point-process model
#'
#' Fits log-linear or additive Poisson point-process intensity models in an
#' arbitrary-dimensional spatial or spatio-temporal window. The point-process
#' likelihood is approximated by the event--dummy cubature scheme constructed
#' by [pp_cubature()] and fitted as a weighted Poisson model with
#' [mgcv::gam()]. Continuous and categorical marks can be included in either a
#' joint marked intensity or a field representation.
#'
#' @inheritParams pp_cubature
#' @param gam_args `NULL` or a named list of additional arguments passed to
#'   [mgcv::gam()]. The function internally fixes `formula`, `family`, `data`,
#'   and `weights`, so those names cannot be supplied. If `method` is omitted,
#'   `method = "REML"` is added. Terms such as `s()`, `te()`, and other
#'   `mgcv` smooth constructors may be used directly in `formula`.
#' @param keep_cubature Character string controlling how much cubature
#'   information is retained in the fitted object. `"full"` stores the complete
#'   `pp_cubature` object, including all event and dummy rows; `"summary"`
#'   stores its formula, domain, mark and covariate specifications, settings,
#'   and numerical summary; `"none"` omits the `cubature` component. The
#'   separate `cubature_summary` and the model specification are retained in all
#'   cases.
#'
#' @details
#' Let \eqn{D} be the continuous process domain and let
#' \eqn{\lambda(u;\theta)=\exp\{\eta(u;\theta)\}} be the model intensity. For
#' ordinary unit event contributions, the Poisson point-process log-likelihood
#' is
#'
#' \deqn{\ell(\theta)=\sum_{i=1}^{n}\eta(u_i;\theta)
#'       -\int_D\exp\{\eta(u;\theta)\}\,\mathrm{d}u.}
#'
#' With `event_weights = v`, the first term becomes
#' \eqn{\sum_i v_i\eta(u_i;\theta)}. The integral is replaced by the
#' multidimensional cubature approximation described in [pp_cubature()]. If
#' \eqn{a_j} is the cubature weight and \eqn{e_j} the event mass at row
#' \eqn{j}, the model is fitted through the equivalent weighted Poisson form
#'
#' \deqn{y_j = \frac{e_j}{a_j}, \qquad w_j=a_j, \qquad
#'       \log\{\mathrm{E}(y_j)\}=\eta(u_j;\theta).}
#'
#' This is a working weighted-Poisson representation rather than an additional
#' point-process assumption. The object stores
#' the directly evaluated approximation
#'
#' \deqn{\widetilde{\ell}(\theta)
#'       =\sum_j e_j\eta(u_j;\theta)
#'       -\sum_j a_j\exp\{\eta(u_j;\theta)\}.}
#'
#' The model formula is one-sided because the internal pseudo-response is
#' generated by the cubature routine. Coordinate variables, external
#' covariates, reconstructed pattern covariates, and declared marks can be used
#' as linear terms, interactions, transformations, or `mgcv` smooth terms.
#'
#' Continuous marks declared with `mode = "joint"` enlarge the integration
#' domain and define a joint location--mark intensity. Marks declared with
#' `mode = "field"` are interpolated from the observed points and enter as
#' explanatory fields. Categorical marks enlarge the domain with a discrete
#' counting measure by replicating the cubature rows over their level
#' combinations.
#'
#' Fractional `event_weights` produce a valid weighted cubature objective but
#' not an ordinary observed-data Poisson likelihood. In that case the returned
#' `aic` is the working-Poisson value calculated by `mgcv`; the
#' `aic_interpretation` component records that it should not be interpreted as a
#' conventional observed-data AIC or compared indiscriminately with unweighted
#' fits.
#'
#' @section Fitted object:
#' The returned `pp_fit` object contains:
#'
#' * `call`, `formula`, and the fitted `mgcv::gam` object in `gam`;
#' * `coefficients` and their covariance matrix `vcov`;
#' * `loglik_approx`, the evaluated cubature objective, and `aic`, the value
#'   returned by [stats::AIC()] for the fitted GAM;
#' * `pattern` and `domain`;
#' * the actual `event_weights`, flags describing whether they were supplied and
#'   whether they are fractional, and `aic_interpretation`;
#' * `cubature`, according to `keep_cubature`, and `cubature_summary`;
#' * `model_spec`, which retains the formula, window and axis metadata,
#'   covariate sources, mark structure, quadrature controls, and GAM arguments
#'   required by prediction, simulation, refitting, and diagnostics;
#' * `settings`, a compact record of fitting controls and total event weight.
#'
#' @return An object of class `"pp_fit"`.
#'
#' @examples
#' w <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#'
#' # Fit a coordinate-dependent intensity
#' x <- simulate_pp(
#'   window = w,
#'   log_intensity = function(coords) {
#'     log(80) + 0.8 * coords$x - 0.5 * coords$y
#'   },
#'   bound = "user",
#'   lambda_upper = exp(log(80) + 0.8),
#'   output = "pattern",
#'   seed = 10
#' )
#'
#' fit <- fit_pp(
#'   pattern = x,
#'   formula = ~ x + y,
#'   nc = c(x = 8, y = 8)
#' )
#' fit
#' summary(fit)
#'
#' # Fit an additive model with an external covariate field
#' field_data <- expand.grid(
#'   x = seq(0, 1, length.out = 6),
#'   y = seq(0, 1, length.out = 6)
#' )
#' field_data$cov1 <- with(field_data, 2 * x + sin(pi * y))
#' cov1 <- pp_covariate(
#'   data = field_data,
#'   coords = c("x", "y"),
#'   value = "cov1",
#'   window = w,
#'   name = "cov1"
#' )
#'
#' fit_cov <- fit_pp(
#'   pattern = x,
#'   formula = ~ x + y + cov1,
#'   covariates = list(cov1 = cov1),
#'   nc = 8
#' )
#'
#' # Joint continuous and categorical marks
#' set.seed(11)
#' xm <- pp_pattern(
#'   data = data.frame(
#'     x = runif(50),
#'     y = runif(50),
#'     magnitude = runif(50, 0, 1),
#'     type = factor(sample(c("A", "B"), 50, replace = TRUE))
#'   ),
#'   coords = c("x", "y"),
#'   marks = c("magnitude", "type"),
#'   window = w
#' )
#'
#' fit_marked <- fit_pp(
#'   pattern = xm,
#'   formula = ~ x + y + magnitude + type,
#'   continuous_marks = list(
#'     magnitude = list(mode = "joint", bounds = c(0, 1))
#'   ),
#'   categorical_marks = "type",
#'   nc = c(x = 5, y = 5, magnitude = 4),
#'   keep_cubature = "summary"
#' )
#'
#' @references
#' Berman, M. and Turner, T. R. (1992). Approximating point process
#' likelihoods with GLIM. *Journal of the Royal Statistical Society: Series C
#' (Applied Statistics)*, **41**(1), 31--38.
#' \doi{10.2307/2347614}
#'
#' Baddeley, A. and Turner, R. (2000). Practical maximum pseudolikelihood for
#' spatial point patterns. *Australian & New Zealand Journal of Statistics*,
#' **42**(3), 283--322. \doi{10.1111/1467-842X.00128}
#'
#' Wood, S. N. (2017). *Generalized Additive Models: An Introduction with R*,
#' 2nd ed. Chapman and Hall/CRC. \doi{10.1201/9781315370279}
#'
#' @seealso
#' [pp_cubature()] constructs the numerical integration scheme;
#' [simulate_pp()] simulates Poisson patterns; [predict_pp_fit()] evaluates
#' fitted joint or ground intensities; [K_fit_pp()] and [residuals_pp()] provide
#' second- and first-order diagnostics.
#'
#' @family Poisson point-process models
#' @export
fit_pp <- function(
    pattern,
    formula = ~ 1,
    covariates = NULL,
    covariate_methods = NULL,
    covariate_args = NULL,
    continuous_marks = NULL,
    categorical_marks = NULL,
    categorical_levels = NULL,
    event_weights = NULL,
    q = 3,
    dummy_layout = c("grid", "random"),
    nc = NULL,
    max_tiles = 1e6,
    max_quadrature_points = 1e6,
    max_mark_combinations = 100L,
    seed = NULL,
    gam_args = list(method = "REML"),
    keep_cubature = c("summary", "full", "none")
) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop(
      "Package 'mgcv' is required to fit Poisson point-process models.",
      call. = FALSE
    )
  }

  .pp_fit_validate_pattern(pattern)

  formula <- .pp_fit_validate_formula(formula)
  dummy_layout <- match.arg(dummy_layout)
  keep_cubature <- match.arg(keep_cubature)
  gam_args <- .pp_fit_normalize_gam_args(gam_args)

  if (nrow(pattern$coords) == 0L) {
    stop(
      "Cannot fit a Poisson point-process model to an empty pattern.",
      call. = FALSE
    )
  }

  event_weights_supplied <- !is.null(event_weights)

  event_weights <- .pp_fit_validate_event_weights(
    pattern = pattern,
    event_weights = event_weights
  )

  if (sum(event_weights) <= 0) {
    stop(
      paste0(
        "`event_weights` must have strictly positive total weight. ",
        "A Poisson point-process model cannot be fitted with no effective events."
      ),
      call. = FALSE
    )
  }

  is_fractionally_weighted <- any(event_weights != 1)

  continuous_marks_normalized <- .pp_fit_normalize_continuous_marks(
    pattern = pattern,
    continuous_marks = continuous_marks
  )

  categorical_marks_normalized <- .pp_fit_normalize_categorical_marks(
    pattern = pattern,
    categorical_marks = categorical_marks,
    categorical_levels = categorical_levels
  )

  overlap_marks <- intersect(
    names(continuous_marks_normalized),
    names(categorical_marks_normalized)
  )

  if (length(overlap_marks) > 0L) {
    stop(
      sprintf(
        "A mark cannot be both continuous and categorical: %s.",
        paste(overlap_marks, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  continuous_marks_canonical <- .pp_fit_canonical_continuous_marks(
    continuous_marks_normalized
  )

  categorical_names <- names(categorical_marks_normalized)

  categorical_levels_canonical <- lapply(
    categorical_marks_normalized,
    function(z) z$levels
  )

  sources <- .pp_fit_prepare_covariate_sources(
    pattern = pattern,
    formula = formula,
    covariates = covariates,
    continuous_marks = continuous_marks_normalized,
    covariate_methods = covariate_methods,
    covariate_args = covariate_args
  )

  cubature <- pp_cubature(
    pattern = pattern,
    formula = formula,
    covariates = covariates,
    covariate_methods = covariate_methods,
    covariate_args = covariate_args,
    continuous_marks = continuous_marks_canonical,
    categorical_marks = categorical_names,
    categorical_levels = categorical_levels_canonical,
    event_weights = event_weights,
    q = q,
    dummy_layout = dummy_layout,
    nc = nc,
    max_tiles = max_tiles,
    max_quadrature_points = max_quadrature_points,
    max_mark_combinations = max_mark_combinations,
    seed = seed
  )

  internal_formula <- .pp_fit_make_gam_formula(formula)

  gam_call <- c(
    list(
      formula = internal_formula,
      family = .pp_poisson_family(),
      data = cubature$data,
      weights = cubature$data$.pp_weight
    ),
    gam_args
  )

  gam_fit <- tryCatch(
    do.call(mgcv::gam, gam_call),
    error = function(error) {
      stop(
        sprintf(
          "`mgcv::gam()` failed: %s",
          conditionMessage(error)
        ),
        call. = FALSE
      )
    }
  )

  eta <- as.numeric(
    stats::predict(gam_fit, type = "link")
  )

  lambda <- exp(eta)

  event_mass <- cubature$data$.pp_event_weight
  weight <- cubature$data$.pp_weight

  loglik_approx <- sum(
    event_mass * eta -
      weight * lambda
  )

  aic_gam <- stats::AIC(gam_fit)

  cubature_value <- switch(
    keep_cubature,

    full = cubature,

    summary = list(
      formula = cubature$formula,
      domain = cubature$domain,
      continuous_marks = cubature$continuous_marks,
      categorical_marks = cubature$categorical_marks,
      covariates = cubature$covariates,
      settings = cubature$settings,
      summary = cubature$summary
    ),

    none = NULL
  )

  model_spec <- .pp_fit_make_model_spec(
    pattern = pattern,
    formula = formula,
    sources = sources,
    continuous_marks = continuous_marks_normalized,
    categorical_marks = categorical_marks_normalized,
    q = q,
    dummy_layout = dummy_layout,
    nc_requested = nc,
    nc_resolved = cubature$settings$nc,
    max_tiles = max_tiles,
    max_quadrature_points = max_quadrature_points,
    max_mark_combinations = max_mark_combinations,
    seed = seed,
    gam_args = gam_args,
    keep_cubature = keep_cubature
  )

  # The actual weights are data-specific and must not be reused automatically
  # by refit_pp() on a simulated or replacement pattern.
  model_spec$event_weights_supplied <- event_weights_supplied
  model_spec$event_weights_reusable <- FALSE

  out <- list(
    call = match.call(),

    formula = formula,
    gam = gam_fit,
    coefficients = stats::coef(gam_fit),
    vcov = stats::vcov(gam_fit),

    loglik_approx = loglik_approx,
    aic = aic_gam,

    pattern = pattern,
    domain = cubature$domain,

    event_weights = event_weights,
    event_weights_supplied = event_weights_supplied,
    is_fractionally_weighted = is_fractionally_weighted,

    aic_interpretation = if (is_fractionally_weighted) {
      paste0(
        "AIC is the mgcv working-Poisson value for fractional event weights; ",
        "it should not be interpreted as an ordinary observed-data AIC."
      )
    } else {
      "Ordinary Poisson cubature AIC."
    },

    cubature = cubature_value,
    cubature_summary = cubature$summary,

    model_spec = model_spec,

    settings = list(
      q = q,
      dummy_layout = dummy_layout,
      nc = cubature$settings$nc,
      gam_args = gam_args,
      keep_cubature = keep_cubature,

      event_weights_supplied = event_weights_supplied,
      is_fractionally_weighted = is_fractionally_weighted,
      total_event_weight = sum(event_weights)
    )
  )

  class(out) <- "pp_fit"

  out
}


#' Methods for fitted Poisson point-process models
#'
#' Print and summarise models returned by [fit_pp()].
#'
#' @details
#' `print.pp_fit()` displays the model formula, cubature size, approximate
#' log-likelihood, working GAM AIC, and estimated coefficients.
#'
#' `summary.pp_fit()` combines the point-process metadata with
#' `summary(object$gam)`, the coefficient covariance matrix, cubature summary,
#' and retained fitting settings. `print.summary_pp_fit()` formats that object.
#'
#' For fractionally weighted fits, the displayed AIC is the working-Poisson AIC
#' returned by `mgcv` and should be interpreted according to the
#' `aic_interpretation` component of the original `pp_fit` object.
#'
#' @param x A `pp_fit` object for `print.pp_fit()`, or a `summary_pp_fit` object
#'   for `print.summary_pp_fit()`.
#' @param object A `pp_fit` object to be summarised.
#' @param digits Number of significant digits used when printing numerical
#'   quantities and coefficients.
#' @param ... Additional arguments. They are currently ignored.
#'
#' @return `print.pp_fit()` and `print.summary_pp_fit()` invisibly return `x`.
#' `summary.pp_fit()` returns an object of class `"summary_pp_fit"`.
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
#'   seed = 12
#' )
#' fit <- fit_pp(x, formula = ~ x + y, nc = 6)
#'
#' print(fit)
#' sf <- summary(fit)
#' print(sf)
#'
#' @seealso [fit_pp()], [pp_cubature()], and [predict_pp_fit()].
#' @family Poisson point-process models
#' @rdname pp_fit-methods
#' @method print pp_fit
#' @export
print.pp_fit <- function(x, digits = 4, ...) {
  if (!inherits(x, "pp_fit")) {
    stop("`x` must inherit from class 'pp_fit'.", call. = FALSE)
  }

  cat("<pp_fit>\n")
  cat("Poisson point-process GAM fitted by cubature\n\n")
  cat("Formula: ")
  print(x$formula)
  cat("\n")
  cat("Observed events: ", x$cubature_summary$n_events, "\n", sep = "")
  cat("Quadrature rows: ", x$cubature_summary$n_rows, "\n", sep = "")
  cat(
    "Approximate log-likelihood: ",
    format(signif(x$loglik_approx, digits), trim = TRUE),
    "\n",
    sep = ""
  )
  cat("AIC (mgcv): ", format(signif(x$aic, digits), trim = TRUE), "\n\n", sep = "")

  cat("Coefficients:\n")
  print(signif(x$coefficients, digits))

  invisible(x)
}


#' @rdname pp_fit-methods
#' @method summary pp_fit
#' @export
summary.pp_fit <- function(object, ...) {
  if (!inherits(object, "pp_fit")) {
    stop("`object` must inherit from class 'pp_fit'.", call. = FALSE)
  }

  out <- list(
    call = object$call,
    formula = object$formula,
    gam_summary = summary(object$gam),
    coefficients = object$coefficients,
    vcov = object$vcov,
    loglik_approx = object$loglik_approx,
    aic = object$aic,
    cubature_summary = object$cubature_summary,
    settings = object$settings
  )

  class(out) <- "summary_pp_fit"
  out
}


#' @rdname pp_fit-methods
#' @method print summary_pp_fit
#' @export
print.summary_pp_fit <- function(x, digits = 4, ...) {
  cat("Poisson point-process GAM summary\n\n")
  cat("Formula: ")
  print(x$formula)
  cat("\n")

  cat("Cubature:\n")
  print(as.data.frame(x$cubature_summary), row.names = FALSE)
  cat("\n")

  cat(
    "Approximate log-likelihood: ",
    format(signif(x$loglik_approx, digits), trim = TRUE),
    "\n",
    sep = ""
  )
  cat("AIC (mgcv): ", format(signif(x$aic, digits), trim = TRUE), "\n\n", sep = "")

  print(x$gam_summary)

  invisible(x)
}
