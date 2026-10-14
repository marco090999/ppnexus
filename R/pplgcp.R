# Log-Gaussian Cox process simulation and fitting -----------------------------
#
# This module implements arbitrary-dimensional spatial and spatio-temporal
# log-Gaussian Cox processes (LGCPs) on the hyperrectangular pp_window objects
# used throughout ppnexus.
#
# The implementation is deliberately modular:
#   * first-order marginal intensity is fitted by fit_pp();
#   * second-order dependence is estimated from K_fit_pp() by minimum contrast;
#   * direct LGCP simulation reuses simulate_pp();
#   * simulation from a fitted LGCP reuses the fitted joint-domain machinery in
#     ppdiagnostic.R, so joint continuous marks, categorical marks, and field
#     marks remain coherent with the fitted pp_fit object.
#
# Required package modules:
#   ppdata.R       -> pp_window(), pp_pattern(), pp_covariate(), pp_grid(),
#                     evaluate_covariate()
#   ppsimul.R      -> simulate_pp() and simulation helpers
#   ppfit.R        -> fit_pp()
#   ppdiagnostic.R -> K_fit_pp(), prediction and fitted-model simulation helpers

.pp_lgcp_version <- "0.1.0"


# Dependencies and generic validation ----------------------------------------

.pp_lgcp_require <- function(functions, context = "LGCP modelling") {
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
          "Source ppdata.R, ppsimul.R, ppfit.R, and ppdiagnostic.R before pplgcp.R."
        ),
        context,
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.pp_lgcp_validate_logical <- function(x, arg_name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", arg_name), call. = FALSE)
  }

  isTRUE(x)
}


.pp_lgcp_validate_positive_integer <- function(x, arg_name, min_value = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != as.integer(x) || x < min_value) {
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


.pp_lgcp_validate_positive_number <- function(
    x,
    arg_name,
    lower = 0,
    inclusive = FALSE
) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop(sprintf("`%s` must be one finite numeric value.", arg_name), call. = FALSE)
  }

  valid <- if (inclusive) x >= lower else x > lower

  if (!valid) {
    comparator <- if (inclusive) "greater than or equal to" else "greater than"
    stop(
      sprintf("`%s` must be %s %s.", arg_name, comparator, format(lower)),
      call. = FALSE
    )
  }

  as.numeric(x)
}


.pp_lgcp_validate_seed <- function(seed) {
  if (is.null(seed)) {
    return(NULL)
  }

  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("`seed` must be NULL or one finite integer.", call. = FALSE)
  }

  as.integer(seed)
}



.pp_lgcp_validate_pattern <- function(pattern) {
  .pp_lgcp_require(c(".pp_diag_validate_pattern"), "LGCP fitting")
  .pp_diag_validate_pattern(pattern, require_window = TRUE)

  axes <- .pp_lgcp_axis_sets(pattern$axis_info)

  if (length(axes$space) < 1L) {
    stop("An LGCP requires at least one spatial coordinate axis.", call. = FALSE)
  }

  invisible(pattern)
}


.pp_lgcp_validate_window <- function(window) {
  .pp_lgcp_require(c(".pp_sim_validate_window"), "LGCP simulation")
  .pp_sim_validate_window(window)

  axes <- .pp_lgcp_axis_sets(window$axis_info)

  if (length(axes$space) < 1L) {
    stop("An LGCP requires at least one spatial coordinate axis.", call. = FALSE)
  }

  invisible(window)
}


.pp_lgcp_axis_sets <- function(axis_info) {
  if (is.null(axis_info) || !is.data.frame(axis_info) ||
      !all(c("name", "role") %in% names(axis_info))) {
    stop("Valid coordinate-axis information is required.", call. = FALSE)
  }

  space <- as.character(axis_info$name[axis_info$role == "space"])
  time <- as.character(axis_info$name[axis_info$role == "time"])

  if (length(time) > 1L) {
    stop("LGCP models support at most one temporal coordinate axis.", call. = FALSE)
  }

  list(
    space = space,
    time = time,
    is_st = length(time) == 1L,
    d_space = length(space)
  )
}


.pp_lgcp_window_widths <- function(window) {
  axis_names <- window$axis_info$name
  widths <- as.numeric(
    window$bounds["upper", axis_names] - window$bounds["lower", axis_names]
  )
  names(widths) <- axis_names

  if (any(!is.finite(widths)) || any(widths <= 0)) {
    stop("Every LGCP window axis must have positive finite width.", call. = FALSE)
  }

  widths
}


.pp_lgcp_window_volume <- function(window) {
  prod(.pp_lgcp_window_widths(window))
}


.pp_lgcp_call_user_function <- function(fun, supplied_args, fun_name) {
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


# Covariance models -----------------------------------------------------------

.pp_lgcp_allowed_models <- function(is_st) {
  if (isTRUE(is_st)) {
    c("separable_exponential", "gneiting", "custom")
  } else {
    c("exponential", "gaussian", "matern", "custom")
  }
}


.pp_lgcp_default_model <- function(is_st) {
  if (isTRUE(is_st)) "separable_exponential" else "exponential"
}


.pp_lgcp_normalize_cov_model <- function(
    cov_model,
    is_st,
    covariance_function = NULL
) {
  if (!is.null(covariance_function) && !is.function(covariance_function)) {
    stop("`covariance_function` must be NULL or a function.", call. = FALSE)
  }

  if (is.null(cov_model)) {
    cov_model <- if (is.null(covariance_function)) {
      .pp_lgcp_default_model(is_st)
    } else {
      "custom"
    }
  }

  if (!is.character(cov_model) || length(cov_model) != 1L || is.na(cov_model)) {
    stop("`cov_model` must be one character string.", call. = FALSE)
  }

  cov_model <- tolower(cov_model)
  allowed <- .pp_lgcp_allowed_models(is_st)

  if (!cov_model %in% allowed) {
    stop(
      sprintf(
        "For this process type, `cov_model` must be one of: %s.",
        paste(allowed, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (identical(cov_model, "custom") && is.null(covariance_function)) {
    stop(
      "`cov_model = 'custom'` requires `covariance_function`.",
      call. = FALSE
    )
  }

  if (!identical(cov_model, "custom") && !is.null(covariance_function)) {
    warning(
      "`covariance_function` is ignored because `cov_model` is not 'custom'.",
      call. = FALSE
    )
    covariance_function <- NULL
  }

  list(model = cov_model, covariance_function = covariance_function)
}


.pp_lgcp_parameter_names <- function(cov_model) {
  switch(
    cov_model,
    exponential = c("sigma2", "alpha"),
    gaussian = c("sigma2", "alpha"),
    matern = c("sigma2", "alpha", "nu"),
    separable_exponential = c("sigma2", "alpha_space", "alpha_time"),
    gneiting = c("sigma2", "alpha_space", "alpha_time", "beta", "gamma"),
    custom = character(0),
    stop("Internal error: unknown LGCP covariance model.", call. = FALSE)
  )
}


.pp_lgcp_as_named_numeric <- function(x, arg_name) {
  if (is.list(x)) {
    if (length(x) == 0L) {
      return(numeric(0))
    }

    x <- unlist(x, recursive = FALSE, use.names = TRUE)
  }

  if (!is.numeric(x) || is.null(names(x)) || any(!nzchar(names(x))) ||
      anyDuplicated(names(x)) || any(!is.finite(x))) {
    stop(
      sprintf("`%s` must be a finite named numeric vector or named list.", arg_name),
      call. = FALSE
    )
  }

  out <- as.numeric(x)
  names(out) <- names(x)
  out
}


.pp_lgcp_validate_model_params <- function(
    params,
    cov_model,
    is_st,
    arg_name = "covariance_params",
    exact = TRUE
) {
  params <- .pp_lgcp_as_named_numeric(params, arg_name)

  if (identical(cov_model, "custom")) {
    if (length(params) == 0L) {
      stop(
        sprintf("`%s` must contain at least one named parameter for a custom covariance.", arg_name),
        call. = FALSE
      )
    }
    return(params)
  }

  required <- .pp_lgcp_parameter_names(cov_model)
  missing <- setdiff(required, names(params))
  extra <- setdiff(names(params), required)

  if (length(missing) > 0L) {
    stop(
      sprintf(
        "`%s` is missing required parameter(s): %s.",
        arg_name,
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (isTRUE(exact) && length(extra) > 0L) {
    stop(
      sprintf(
        "`%s` contains unsupported parameter(s) for `%s`: %s.",
        arg_name,
        cov_model,
        paste(extra, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  params <- params[required]

  positive <- switch(
    cov_model,
    exponential = c("sigma2", "alpha"),
    gaussian = c("sigma2", "alpha"),
    matern = c("sigma2", "alpha", "nu"),
    separable_exponential = c("sigma2", "alpha_space", "alpha_time"),
    gneiting = c("sigma2", "alpha_space", "alpha_time"),
    character(0)
  )

  if (length(positive) > 0L && any(params[positive] <= 0)) {
    stop(
      sprintf(
        "The following `%s` parameter(s) must be strictly positive: %s.",
        arg_name,
        paste(positive[params[positive] <= 0], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (identical(cov_model, "gneiting")) {
    if (params[["beta"]] <= 0 || params[["beta"]] > 1) {
      stop("For the Gneiting model, `beta` must lie in (0, 1].", call. = FALSE)
    }

    if (params[["gamma"]] <= 0 || params[["gamma"]] > 1) {
      stop("For the Gneiting model, `gamma` must lie in (0, 1].", call. = FALSE)
    }
  }

  params
}


.pp_lgcp_call_covariance_function <- function(
    covariance_function,
    h,
    u,
    params,
    d_space
) {
  value <- .pp_lgcp_call_user_function(
    fun = covariance_function,
    supplied_args = list(
      h = h,
      u = u,
      params = params,
      d_space = d_space
    ),
    fun_name = "covariance_function"
  )

  if (!is.numeric(value)) {
    stop("`covariance_function` must return numeric values.", call. = FALSE)
  }

  if (length(value) == 1L && length(h) > 1L) {
    value <- rep(value, length(h))
    dim(value) <- dim(h)
  }

  if (length(value) != length(h)) {
    stop(
      "`covariance_function` must return one value per supplied lag.",
      call. = FALSE
    )
  }

  if (!is.null(dim(h))) {
    dim(value) <- dim(h)
    dimnames(value) <- dimnames(h)
  }

  if (any(!is.finite(value))) {
    stop("`covariance_function` returned non-finite values.", call. = FALSE)
  }

  value
}


.pp_lgcp_covariance <- function(
    h,
    u = NULL,
    params,
    cov_model,
    d_space,
    covariance_function = NULL
) {
  if (!is.numeric(h) || any(!is.finite(h)) || any(h < 0)) {
    stop("Spatial lags `h` must be finite and non-negative.", call. = FALSE)
  }

  is_st <- !is.null(u)

  if (is_st) {
    if (!is.numeric(u) || length(u) != length(h) ||
        any(!is.finite(u)) || any(u < 0)) {
      stop(
        "Temporal lags `u` must be finite, non-negative, and aligned with `h`.",
        call. = FALSE
      )
    }

    if (!is.null(dim(h)) && !identical(dim(h), dim(u))) {
      stop("`h` and `u` must have identical dimensions.", call. = FALSE)
    }
  }

  params <- .pp_lgcp_validate_model_params(
    params = params,
    cov_model = cov_model,
    is_st = is_st,
    arg_name = "params",
    exact = !identical(cov_model, "custom")
  )

  out <- switch(
    cov_model,

    exponential = {
      params[["sigma2"]] * exp(-h / params[["alpha"]])
    },

    gaussian = {
      params[["sigma2"]] * exp(-(h / params[["alpha"]])^2)
    },

    matern = {
      z <- h / params[["alpha"]]
      corr <- numeric(length(z))
      zero <- z == 0
      corr[zero] <- 1

      if (any(!zero)) {
        nu <- params[["nu"]]
        z_positive <- z[!zero]
        corr[!zero] <-
          2^(1 - nu) / gamma(nu) *
          z_positive^nu * besselK(z_positive, nu = nu)
      }

      if (!is.null(dim(h))) {
        dim(corr) <- dim(h)
      }

      params[["sigma2"]] * corr
    },

    separable_exponential = {
      params[["sigma2"]] *
        exp(-h / params[["alpha_space"]]) *
        exp(-u / params[["alpha_time"]])
    },

    gneiting = {
      # A valid Gneiting-type subclass:
      # psi(u^2) = (1 + (u / alpha_time)^2)^beta,
      # C(h,u) = sigma2 * psi^(-d/2) *
      #          exp(-((h / alpha_space)^2 / psi)^gamma).
      psi <- (1 + (u / params[["alpha_time"]])^2)^params[["beta"]]

      params[["sigma2"]] *
        psi^(-d_space / 2) *
        exp(
          -(
            (h / params[["alpha_space"]])^2 / psi
          )^params[["gamma"]]
        )
    },

    custom = {
      .pp_lgcp_call_covariance_function(
        covariance_function = covariance_function,
        h = h,
        u = u,
        params = params,
        d_space = d_space
      )
    },

    stop("Internal error: unsupported covariance model.", call. = FALSE)
  )

  if (any(!is.finite(out))) {
    stop("The covariance model produced non-finite values.", call. = FALSE)
  }

  out
}


.pp_lgcp_variance <- function(
    params,
    cov_model,
    is_st,
    d_space,
    covariance_function = NULL
) {
  variance <- .pp_lgcp_covariance(
    h = 0,
    u = if (isTRUE(is_st)) 0 else NULL,
    params = params,
    cov_model = cov_model,
    d_space = d_space,
    covariance_function = covariance_function
  )

  variance <- as.numeric(variance[1L])

  if (!is.finite(variance) || variance <= 0) {
    stop(
      paste0(
        "The covariance at zero lag must be strictly positive and finite. ",
        "This value is used to centre the latent Gaussian field."
      ),
      call. = FALSE
    )
  }

  variance
}


# Latent-field support and simulation ----------------------------------------

.pp_lgcp_default_field_n <- function(window, max_field_points) {
  max_field_points <- .pp_lgcp_validate_positive_integer(
    max_field_points,
    "max_field_points",
    min_value = 2L
  )

  d <- nrow(window$axis_info)
  target <- min(1200L, max_field_points)
  per_axis <- max(2L, floor(target^(1 / d)))

  out <- rep.int(as.integer(per_axis), d)
  names(out) <- window$axis_info$name

  while (prod(as.double(out)) > max_field_points) {
    j <- which.max(out)
    out[j] <- out[j] - 1L

    if (any(out < 2L)) {
      stop(
        "`max_field_points` is too small for a grid with at least two points per axis.",
        call. = FALSE
      )
    }
  }

  out
}


.pp_lgcp_validate_field_support <- function(
    field_support,
    window,
    field_n,
    max_field_points
) {
  .pp_lgcp_require(c("pp_grid"), "LGCP latent-field construction")

  axis_names <- window$axis_info$name

  if (is.null(field_support)) {
    if (is.null(field_n)) {
      field_n <- .pp_lgcp_default_field_n(window, max_field_points)
    }

    grid <- pp_grid(
      window = window,
      n = field_n,
      max_points = max_field_points
    )

    coords <- grid$coords
    resolution <- grid$resolution
    source <- "regular_grid"
  } else if (inherits(field_support, "pp_grid")) {
    if (!identical(field_support$axis_info$name, axis_names) ||
        !isTRUE(all.equal(field_support$window$bounds, window$bounds))) {
      stop(
        "`field_support` must use the same coordinate axes and window as `window`.",
        call. = FALSE
      )
    }

    coords <- field_support$coords
    resolution <- field_support$resolution
    source <- "pp_grid"
  } else if (is.data.frame(field_support) || is.matrix(field_support)) {
    support_df <- as.data.frame(field_support, check.names = FALSE)

    missing_axes <- setdiff(axis_names, names(support_df))

    if (length(missing_axes) > 0L) {
      stop(
        sprintf(
          "`field_support` is missing coordinate column(s): %s.",
          paste(missing_axes, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    coords <- as.matrix(support_df[, axis_names, drop = FALSE])
    storage.mode(coords) <- "double"
    colnames(coords) <- axis_names
    resolution <- NULL
    source <- "user"
  } else {
    stop(
      "`field_support` must be NULL, a pp_grid object, a data frame, or a matrix.",
      call. = FALSE
    )
  }

  if (!is.numeric(coords) || nrow(coords) < 2L || any(!is.finite(coords))) {
    stop(
      "The latent-field support must contain at least two finite coordinate locations.",
      call. = FALSE
    )
  }

  max_field_points <- .pp_lgcp_validate_positive_integer(
    max_field_points,
    "max_field_points",
    min_value = 2L
  )

  if (nrow(coords) > max_field_points) {
    stop(
      sprintf(
        "The latent-field support contains %d points, exceeding `max_field_points = %d`.",
        nrow(coords),
        max_field_points
      ),
      call. = FALSE
    )
  }

  lower <- window$bounds["lower", axis_names]
  upper <- window$bounds["upper", axis_names]
  tolerance <- .Machine$double.eps^0.5 * pmax(1, abs(upper - lower))

  below <- sweep(coords, 2L, lower - tolerance, FUN = "<")
  above <- sweep(coords, 2L, upper + tolerance, FUN = ">")

  if (any(below) || any(above)) {
    stop("Every latent-field support location must lie inside `window`.", call. = FALSE)
  }

  list(
    coords = coords,
    resolution = resolution,
    source = source
  )
}


.pp_lgcp_distance_matrices <- function(coords, axis_info) {
  axes <- .pp_lgcp_axis_sets(axis_info)
  n <- nrow(coords)

  h2 <- matrix(0, nrow = n, ncol = n)

  for (name_i in axes$space) {
    difference <- outer(coords[, name_i], coords[, name_i], FUN = "-")
    h2 <- h2 + difference^2
  }

  h <- sqrt(pmax(h2, 0))

  u <- if (axes$is_st) {
    abs(outer(coords[, axes$time], coords[, axes$time], FUN = "-"))
  } else {
    NULL
  }

  list(h = h, u = u, axes = axes)
}


.pp_lgcp_chol <- function(C, jitter, max_attempts) {
  if (!is.matrix(C) || nrow(C) != ncol(C) || any(!is.finite(C))) {
    stop("The latent-field covariance matrix must be finite and square.", call. = FALSE)
  }

  jitter <- .pp_lgcp_validate_positive_number(
    jitter,
    "jitter",
    lower = 0,
    inclusive = TRUE
  )

  max_attempts <- .pp_lgcp_validate_positive_integer(
    max_attempts,
    "max_jitter_attempts",
    min_value = 1L
  )

  base_scale <- max(1, max(diag(C)))
  current_jitter <- jitter * base_scale

  for (attempt in seq_len(max_attempts)) {
    C_attempt <- C

    if (current_jitter > 0) {
      diag(C_attempt) <- diag(C_attempt) + current_jitter
    }

    R <- try(chol(C_attempt), silent = TRUE)

    if (!inherits(R, "try-error")) {
      return(list(
        R = R,
        jitter_used = current_jitter,
        attempts = attempt
      ))
    }

    current_jitter <- if (current_jitter == 0) {
      .Machine$double.eps^0.5 * base_scale
    } else {
      current_jitter * 10
    }
  }

  stop(
    paste0(
      "Cholesky factorisation of the latent covariance matrix failed after ",
      max_attempts,
      " attempts. Check covariance parameters, custom covariance validity, ",
      "or increase `jitter` deliberately."
    ),
    call. = FALSE
  )
}


.pp_lgcp_simulate_field <- function(
    window,
    cov_model,
    covariance_params,
    covariance_function,
    field_support = NULL,
    field_n = NULL,
    max_field_points = 2000L,
    jitter = 1e-10,
    max_jitter_attempts = 8L
) {
  support <- .pp_lgcp_validate_field_support(
    field_support = field_support,
    window = window,
    field_n = field_n,
    max_field_points = max_field_points
  )

  distances <- .pp_lgcp_distance_matrices(
    coords = support$coords,
    axis_info = window$axis_info
  )

  params <- .pp_lgcp_validate_model_params(
    params = covariance_params,
    cov_model = cov_model,
    is_st = distances$axes$is_st,
    arg_name = "covariance_params",
    exact = !identical(cov_model, "custom")
  )

  C <- .pp_lgcp_covariance(
    h = distances$h,
    u = distances$u,
    params = params,
    cov_model = cov_model,
    d_space = distances$axes$d_space,
    covariance_function = covariance_function
  )

  # Numerical symmetry is enforced before Cholesky factorisation.
  C <- (C + t(C)) / 2

  chol_result <- .pp_lgcp_chol(
    C = C,
    jitter = jitter,
    max_attempts = max_jitter_attempts
  )

  z <- as.numeric(t(chol_result$R) %*% stats::rnorm(nrow(support$coords)))

  variance <- .pp_lgcp_variance(
    params = params,
    cov_model = cov_model,
    is_st = distances$axes$is_st,
    d_space = distances$axes$d_space,
    covariance_function = covariance_function
  )

  out <- list(
    coords = support$coords,
    value = z,
    centred_value = z - variance / 2,
    variance = variance,
    covariance_model = cov_model,
    covariance_params = params,
    covariance_function = covariance_function,
    support_source = support$source,
    support_resolution = support$resolution,
    jitter_used = chol_result$jitter_used,
    chol_attempts = chol_result$attempts,
    window = window
  )

  class(out) <- "pp_lgcp_field"
  out
}


.pp_lgcp_field_covariate <- function(field, name = ".pp_lgcp_latent") {
  .pp_lgcp_require(c("pp_covariate"), "LGCP latent-field interpolation")

  axis_names <- field$window$axis_info$name
  data <- as.data.frame(field$coords, check.names = FALSE)
  data[[name]] <- field$value

  pp_covariate(
    data = data,
    coords = axis_names,
    value = name,
    window = field$window,
    name = name
  )
}


.pp_lgcp_validate_field_method <- function(field_method) {
  if (!is.character(field_method) || length(field_method) != 1L || is.na(field_method)) {
    stop("`field_method` must be one character string.", call. = FALSE)
  }

  match.arg(field_method, c("idw", "nearest"))
}


.pp_lgcp_strip_latent_covariate <- function(covariates, latent_name) {
  if (is.null(covariates)) {
    return(NULL)
  }

  covariates <- as.data.frame(covariates, check.names = FALSE)
  keep <- setdiff(names(covariates), latent_name)

  if (length(keep) == 0L) {
    return(NULL)
  }

  covariates[, keep, drop = FALSE]
}


.pp_lgcp_wrap_mark_generators <- function(mark_generators, latent_name) {
  if (is.null(mark_generators)) {
    return(NULL)
  }

  if (!is.list(mark_generators)) {
    stop("`mark_generators` must be NULL or a named list of functions.", call. = FALSE)
  }

  if (length(mark_generators) > 0L &&
      (is.null(names(mark_generators)) || any(!nzchar(names(mark_generators))) ||
       anyDuplicated(names(mark_generators)))) {
    stop("`mark_generators` must have unique non-empty names.", call. = FALSE)
  }

  out <- vector("list", length(mark_generators))
  names(out) <- names(mark_generators)

  for (i in seq_along(mark_generators)) {
    generator <- mark_generators[[i]]
    generator_name <- names(mark_generators)[i]

    if (!is.function(generator)) {
      stop(
        sprintf("`mark_generators[[\"%s\"]]` must be a function.", generator_name),
        call. = FALSE
      )
    }

    out[[i]] <- local({
      fun_i <- generator
      name_i <- generator_name
      latent_i <- latent_name

      function(n, coords, covariates = NULL, params = list()) {
        .pp_lgcp_call_user_function(
          fun = fun_i,
          supplied_args = list(
            n = n,
            coords = coords,
            covariates = .pp_lgcp_strip_latent_covariate(covariates, latent_i),
            params = params
          ),
          fun_name = paste0("mark_generators[[\"", name_i, "\"]]")
        )
      }
    })
  }

  out
}


.pp_lgcp_attach_simulation_metadata <- function(
    x,
    field,
    field_method,
    field_args,
    keep_field,
    source
) {
  metadata <- list(
    source = source,
    covariance_model = field$covariance_model,
    covariance_params = field$covariance_params,
    marginal_variance = field$variance,
    centred = TRUE,
    field_method = field_method,
    field_args = field_args,
    support_source = field$support_source,
    support_resolution = field$support_resolution,
    n_field_points = nrow(field$coords),
    jitter_used = field$jitter_used,
    field = if (isTRUE(keep_field)) field else NULL
  )

  if (inherits(x, "pp")) {
    x$metadata$lgcp <- metadata
    return(x)
  }

  if (inherits(x, "pp_simulation") || inherits(x, "pp_fit_simulation")) {
    if (!is.null(x$pattern) && inherits(x$pattern, "pp")) {
      x$pattern$metadata$lgcp <- metadata
    }

    x$lgcp <- metadata
    class(x) <- unique(c("pp_lgcp_simulation", class(x)))
    return(x)
  }

  x
}


#' Simulate a multidimensional log-Gaussian Cox process
#'
#' Simulates a spatial or spatio-temporal log-Gaussian Cox process (LGCP) on an
#' arbitrary-dimensional hyperrectangular `pp_window`. The marginal first-order
#' intensity is specified exactly as in `simulate_pp()`, while second-order
#' clustering is induced by a latent Gaussian random field.
#'
#' @details
#' Conditional on a zero-mean Gaussian field `Z(u)` with covariance
#' `C_theta(u,v)`, the process intensity is parameterised as
#' \deqn{\Lambda(u) = \rho(u)\exp\{Z(u)-C_\theta(u,u)/2\},}
#' where `rho(u)` is the user-specified marginal intensity. Consequently,
#' \deqn{E\{\Lambda(u)\}=\rho(u)}
#' and the pair correlation function is
#' \deqn{g(u,v)=\exp\{C_\theta(u,v)\}.}
#'
#' The Gaussian field is simulated on a finite support using dense covariance
#' matrix factorisation and then interpolated to thinning candidates using
#' either inverse-distance weighting or nearest-neighbour evaluation. This is a
#' discretised approximation to a continuous LGCP. Each replicate receives an
#' independently simulated latent field, including when `nsim > 1`.
#'
#' Spatial models support `"exponential"`, `"gaussian"`, and `"matern"`
#' covariance functions. Spatio-temporal models support
#' `"separable_exponential"` and a nonseparable Gneiting-type model. A custom
#' stationary isotropic covariance may be supplied through
#' `cov_model = "custom"` and `covariance_function`.
#'
#' @param window A `pp_window` object with at least one spatial axis and at most
#'   one temporal axis.
#' @param expected_n Marginal expected number of events for a homogeneous
#'   first-order intensity. Supply exactly one of `expected_n` or
#'   `log_intensity`.
#' @param log_intensity Function returning the marginal log-intensity at
#'   candidate locations. It may use any subset of `coords`, `covariates`,
#'   `marks`, and `params`, exactly as in `simulate_pp()`.
#' @param params Optional list passed to `log_intensity` and mark generators.
#' @param covariates Optional external `pp_covariate` object or named list of
#'   such objects used by the marginal first-order intensity.
#' @param covariate_methods,covariate_args Covariate interpolation settings
#'   passed to `simulate_pp()`.
#' @param mark_generators Optional named list of mark-generating functions.
#' @param mark_mode Either `"posthoc"` or `"joint"`, as in `simulate_pp()`.
#' @param cov_model Latent covariance model. For purely spatial processes use
#'   `"exponential"`, `"gaussian"`, `"matern"`, or `"custom"`. For
#'   spatio-temporal processes use `"separable_exponential"`, `"gneiting"`, or
#'   `"custom"`.
#' @param covariance_params Named numeric vector or named list of covariance
#'   parameters. Spatial exponential/Gaussian models require `sigma2, alpha`;
#'   Matérn requires `sigma2, alpha, nu`; separable space-time exponential
#'   requires `sigma2, alpha_space, alpha_time`; the Gneiting model additionally
#'   requires `beta, gamma`, both in `(0, 1]`.
#' @param covariance_function Optional custom covariance function. It may use
#'   `h`, `u`, `params`, and `d_space`, and must return one covariance value per
#'   lag. For a purely spatial model `u` is `NULL`. The user is responsible for
#'   supplying a valid stationary isotropic positive-semidefinite covariance.
#' @param field_support Optional latent-field support. It may be a `pp_grid`, a
#'   data frame, or a matrix containing all window coordinate axes. If `NULL`,
#'   a regular grid is generated automatically.
#' @param field_n Optional number of support locations per axis when
#'   `field_support = NULL`. It follows the conventions of `pp_grid()`.
#' @param field_method Interpolation method used to evaluate the latent field at
#'   thinning candidates: `"idw"` or `"nearest"`.
#' @param field_args Optional list of interpolation arguments such as `p`, `k`,
#'   `axis_scale`, and `chunk_size`.
#' @param max_field_points Maximum permitted number of latent-field support
#'   points. Dense Gaussian simulation has cubic time complexity in this number.
#' @param jitter Initial non-negative diagonal jitter used for Cholesky
#'   factorisation of the covariance matrix.
#' @param max_jitter_attempts Maximum number of progressively larger jitter
#'   attempts when Cholesky factorisation fails numerically.
#' @param keep_field Logical; retain the simulated latent field in output
#'   metadata.
#' @param bound,lambda_upper,bound_n,bound_max_points,bound_safety Thinning
#'   bound settings passed to `simulate_pp()`.
#' @param bound_mark_reps,on_bound_violation Additional mark-bound and
#'   bound-violation settings passed to `simulate_pp()`.
#' @param restart_safety,max_restarts,max_candidates Thinning restart and
#'   candidate-limit settings passed to `simulate_pp()`. A user-supplied
#'   `lambda_upper` must bound the full conditional LGCP intensity for the
#'   realised latent field.
#' @param output One of `"pattern"`, `"enriched"`, or `"diagnostic"`, as in
#'   `simulate_pp()`.
#' @param nsim Number of independent LGCP simulations.
#' @param seed Optional integer random seed.
#'
#' @return For `nsim = 1`, an object following the selected `simulate_pp()`
#'   output mode with additional LGCP metadata. Diagnostic outputs inherit from
#'   `pp_lgcp_simulation`. For `nsim > 1`, a named list of independent outputs.
#'
#' @references
#' Møller, J., Syversveen, A. R., and Waagepetersen, R. P. (1998).
#' Log Gaussian Cox processes. *Scandinavian Journal of Statistics*, 25,
#' 451--482. doi:10.1111/1467-9469.00115.
#'
#' Gneiting, T. (2002). Nonseparable, stationary covariance functions for
#' space-time data. *Journal of the American Statistical Association*, 97,
#' 590--600. doi:10.1198/016214502760047113.
#'
#' Lewis, P. A. W. and Shedler, G. S. (1979). Simulation of nonhomogeneous
#' Poisson processes by thinning. *Naval Research Logistics Quarterly*, 26,
#' 403--413.
#'
#' @seealso `simulate_pp()`, `fit_lgcp_pp()`, `pp_grid()`, `pp_covariate()`
#'
#' @examples
#' w <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#'
#' sim <- simulate_lgcp_pp(
#'   window = w,
#'   expected_n = 80,
#'   cov_model = "exponential",
#'   covariance_params = c(sigma2 = 0.5, alpha = 0.15),
#'   field_n = 8,
#'   output = "pattern",
#'   seed = 1
#' )
#'
#' summary(sim)
#'
#' @export
simulate_lgcp_pp <- function(
    window,
    expected_n = NULL,
    log_intensity = NULL,
    params = list(),
    covariates = NULL,
    covariate_methods = NULL,
    covariate_args = NULL,
    mark_generators = NULL,
    mark_mode = c("posthoc", "joint"),
    cov_model = NULL,
    covariance_params,
    covariance_function = NULL,
    field_support = NULL,
    field_n = NULL,
    field_method = c("nearest", "idw"),
    field_args = list(),
    max_field_points = 2000L,
    jitter = 1e-10,
    max_jitter_attempts = 8L,
    keep_field = TRUE,
    bound = c("grid", "user"),
    lambda_upper = NULL,
    bound_n = NULL,
    bound_max_points = 50000L,
    bound_safety = 1.20,
    bound_mark_reps = 20L,
    on_bound_violation = c("restart", "error"),
    restart_safety = 1.20,
    max_restarts = 10L,
    output = c("pattern", "enriched", "diagnostic"),
    nsim = 1L,
    seed = NULL,
    max_candidates = 1e6
) {
  .pp_lgcp_require(
    c(
      "simulate_pp", ".pp_sim_normalize_covariates",
      ".pp_sim_normalize_covariate_methods", ".pp_sim_normalize_covariate_args"
    ),
    "simulate_lgcp_pp()"
  )

  .pp_lgcp_validate_window(window)

  has_expected_n <- !is.null(expected_n)
  has_log_intensity <- !is.null(log_intensity)

  if (identical(has_expected_n, has_log_intensity)) {
    stop(
      "Supply exactly one of `expected_n` or `log_intensity`.",
      call. = FALSE
    )
  }

  if (identical(field_method, "idw")) {
    warning(
      paste0(
        "IDW interpolation changes the marginal variance of the latent ",
        "Gaussian field. Consequently, centering by C(0) / 2 may not ",
        "preserve the requested marginal first-order intensity exactly. ",
        "Use `field_method = 'nearest'` for the default variance-preserving ",
        "finite-support approximation."
      ),
      call. = FALSE
    )
  }

  if (has_expected_n) {
    if (!is.numeric(expected_n) || length(expected_n) != 1L ||
        !is.finite(expected_n) || expected_n < 0) {
      stop("`expected_n` must be one finite non-negative number.", call. = FALSE)
    }
  } else if (!is.function(log_intensity)) {
    stop("`log_intensity` must be a function.", call. = FALSE)
  }

  if (is.null(params)) {
    params <- list()
  }

  if (!is.list(params)) {
    stop("`params` must be a list or NULL.", call. = FALSE)
  }

  nsim <- .pp_lgcp_validate_positive_integer(nsim, "nsim", min_value = 1L)
  seed <- .pp_lgcp_validate_seed(seed)
  keep_field <- .pp_lgcp_validate_logical(keep_field, "keep_field")
  field_method <- .pp_lgcp_validate_field_method(match.arg(field_method))
  field_args <- if (is.null(field_args)) list() else field_args

  if (!is.list(field_args)) {
    stop("`field_args` must be a list or NULL.", call. = FALSE)
  }

  mark_mode <- match.arg(mark_mode)
  bound <- match.arg(bound)
  on_bound_violation <- match.arg(on_bound_violation)
  output <- match.arg(output)

  axes <- .pp_lgcp_axis_sets(window$axis_info)
  covariance_spec <- .pp_lgcp_normalize_cov_model(
    cov_model = cov_model,
    is_st = axes$is_st,
    covariance_function = covariance_function
  )

  covariance_params <- .pp_lgcp_validate_model_params(
    params = covariance_params,
    cov_model = covariance_spec$model,
    is_st = axes$is_st,
    arg_name = "covariance_params",
    exact = !identical(covariance_spec$model, "custom")
  )

  base_covariates <- .pp_sim_normalize_covariates(covariates)
  base_methods <- .pp_sim_normalize_covariate_methods(
    covariates = base_covariates,
    methods = covariate_methods
  )
  base_args <- .pp_sim_normalize_covariate_args(
    covariates = base_covariates,
    covariate_args = covariate_args
  )

  latent_name <- ".pp_lgcp_latent"

  if (latent_name %in% names(base_covariates)) {
    stop(
      sprintf("The internal covariate name `%s` is reserved by LGCP simulation.", latent_name),
      call. = FALSE
    )
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  results <- vector("list", nsim)

  for (simulation_id in seq_len(nsim)) {
    field <- .pp_lgcp_simulate_field(
      window = window,
      cov_model = covariance_spec$model,
      covariance_params = covariance_params,
      covariance_function = covariance_spec$covariance_function,
      field_support = field_support,
      field_n = field_n,
      max_field_points = max_field_points,
      jitter = jitter,
      max_jitter_attempts = max_jitter_attempts
    )

    latent_covariate <- .pp_lgcp_field_covariate(field, name = latent_name)

    all_covariates <- c(base_covariates, stats::setNames(list(latent_covariate), latent_name))
    all_methods <- c(base_methods, stats::setNames(field_method, latent_name))
    all_args <- c(base_args, stats::setNames(list(field_args), latent_name))

    wrapped_mark_generators <- .pp_lgcp_wrap_mark_generators(
      mark_generators = mark_generators,
      latent_name = latent_name
    )

    marginal_log_intensity <- log_intensity
    expected_n_i <- expected_n
    variance_i <- field$variance
    volume_i <- .pp_lgcp_window_volume(window)

    lgcp_log_intensity <- local({
      user_fun <- marginal_log_intensity
      expected_count <- expected_n_i
      latent_i <- latent_name
      variance_value <- variance_i
      volume_value <- volume_i

      function(coords, covariates = NULL, marks = NULL, params = list()) {
        latent <- covariates[[latent_i]]

        base_covs <- .pp_lgcp_strip_latent_covariate(
          covariates = covariates,
          latent_name = latent_i
        )

        eta_base <- if (!is.null(expected_count)) {
          if (expected_count == 0) {
            rep(-Inf, nrow(coords))
          } else {
            rep(log(expected_count / volume_value), nrow(coords))
          }
        } else {
          .pp_lgcp_call_user_function(
            fun = user_fun,
            supplied_args = list(
              coords = coords,
              covariates = base_covs,
              marks = marks,
              params = params
            ),
            fun_name = "log_intensity"
          )
        }

        if (!is.numeric(eta_base) || length(eta_base) != nrow(coords)) {
          stop(
            "`log_intensity` must return one numeric value per candidate location.",
            call. = FALSE
          )
        }

        as.numeric(eta_base) + latent - variance_value / 2
      }
    })

    # For an LGCP the intensity contains a random latent field, so an exact
    # user bound cannot be supplied in general: the grid bound, combined with
    # automatic restarts on detected violations, is the intended mechanism.
    # The generic simulate_pp() warning about grid bounds is therefore not
    # actionable here and is muffled; any other warning is still signalled.
    result_i <- .pp_lgcp_muffle_grid_bound_warning(simulate_pp(
      window = window,
      expected_n = NULL,
      log_intensity = lgcp_log_intensity,
      params = params,
      covariates = all_covariates,
      covariate_methods = all_methods,
      covariate_args = all_args,
      mark_generators = wrapped_mark_generators,
      mark_mode = mark_mode,
      bound = bound,
      lambda_upper = lambda_upper,
      bound_n = bound_n,
      bound_max_points = bound_max_points,
      bound_safety = bound_safety,
      bound_mark_reps = bound_mark_reps,
      on_bound_violation = on_bound_violation,
      restart_safety = restart_safety,
      max_restarts = max_restarts,
      output = output,
      nsim = 1L,
      seed = NULL,
      max_candidates = max_candidates
    ))

    results[[simulation_id]] <- .pp_lgcp_attach_simulation_metadata(
      x = result_i,
      field = field,
      field_method = field_method,
      field_args = field_args,
      keep_field = keep_field,
      source = "simulate_lgcp_pp"
    )
  }

  if (nsim == 1L) {
    return(results[[1L]])
  }

  names(results) <- paste0("sim_", seq_len(nsim))
  results
}


# Theoretical K-functions -----------------------------------------------------

.pp_lgcp_unit_ball_volume <- function(d) {
  d <- .pp_lgcp_validate_positive_integer(d, "d", min_value = 1L)
  pi^(d / 2) / gamma(d / 2 + 1)
}


.pp_lgcp_unit_sphere_area <- function(d) {
  d * .pp_lgcp_unit_ball_volume(d)
}


.pp_lgcp_validate_lags <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) == 0L || any(!is.finite(x)) || any(x < 0)) {
    stop(
      sprintf("`%s` must contain finite non-negative lag values.", arg_name),
      call. = FALSE
    )
  }

  sort(unique(as.numeric(x)))
}


.pp_lgcp_cumtrapz <- function(x, y) {
  if (!is.numeric(x) || !is.numeric(y) || length(x) != length(y) || length(x) < 1L) {
    stop("Internal LGCP integration error.", call. = FALSE)
  }

  if (length(x) == 1L) {
    return(0)
  }

  increments <- diff(x) * (y[-length(y)] + y[-1L]) / 2
  c(0, cumsum(increments))
}


.pp_lgcp_theoretical_k_spatial <- function(
    r,
    d_space,
    params,
    cov_model,
    covariance_function,
    integration_n = 400L
) {
  r <- .pp_lgcp_validate_lags(r, "r")
  integration_n <- .pp_lgcp_validate_positive_integer(
    integration_n,
    "integration_n",
    min_value = 20L
  )

  max_r <- max(r)

  if (max_r == 0) {
    return(rep(0, length(r)))
  }

  h_grid <- sort(unique(c(
    0,
    r,
    seq(0, max_r, length.out = integration_n)
  )))

  C <- .pp_lgcp_covariance(
    h = h_grid,
    u = NULL,
    params = params,
    cov_model = cov_model,
    d_space = d_space,
    covariance_function = covariance_function
  )

  if (any(C > log(.Machine$double.xmax) - 10)) {
    return(rep(Inf, length(r)))
  }

  integrand <- h_grid^(d_space - 1) * exp(C)
  cumulative <- .pp_lgcp_cumtrapz(h_grid, integrand)

  area <- .pp_lgcp_unit_sphere_area(d_space)
  area * cumulative[match(r, h_grid)]
}


.pp_lgcp_theoretical_k_st <- function(
    r,
    u,
    d_space,
    params,
    cov_model,
    covariance_function,
    integration_n = c(100L, 100L)
) {
  r <- .pp_lgcp_validate_lags(r, "r")
  u <- .pp_lgcp_validate_lags(u, "u")

  if (!is.numeric(integration_n) || any(!is.finite(integration_n)) ||
      any(integration_n != as.integer(integration_n))) {
    stop("`integration_n` must contain positive integer resolutions.", call. = FALSE)
  }

  if (length(integration_n) == 1L) {
    integration_n <- rep.int(as.integer(integration_n), 2L)
  } else if (length(integration_n) == 2L) {
    integration_n <- as.integer(integration_n)
  } else {
    stop("For a spatio-temporal LGCP, `integration_n` must have length one or two.", call. = FALSE)
  }

  if (any(integration_n < 20L)) {
    stop("Every `integration_n` value must be at least 20.", call. = FALSE)
  }

  max_r <- max(r)
  max_u <- max(u)

  h_grid <- sort(unique(c(
    0,
    r,
    if (max_r > 0) seq(0, max_r, length.out = integration_n[1L]) else 0
  )))

  u_grid <- sort(unique(c(
    0,
    u,
    if (max_u > 0) seq(0, max_u, length.out = integration_n[2L]) else 0
  )))

  H <- outer(h_grid, u_grid, function(a, b) a)
  U <- outer(h_grid, u_grid, function(a, b) b)

  C <- .pp_lgcp_covariance(
    h = H,
    u = U,
    params = params,
    cov_model = cov_model,
    d_space = d_space,
    covariance_function = covariance_function
  )

  if (any(C > log(.Machine$double.xmax) - 10)) {
    return(matrix(Inf, nrow = length(r), ncol = length(u)))
  }

  F <- H^(d_space - 1) * exp(C)

  # First integrate cumulatively over non-negative temporal lags.
  cumulative_u <- matrix(0, nrow = length(h_grid), ncol = length(u_grid))

  if (length(u_grid) > 1L) {
    for (j in 2:length(u_grid)) {
      du <- u_grid[j] - u_grid[j - 1L]
      cumulative_u[, j] <- cumulative_u[, j - 1L] +
        du * (F[, j - 1L] + F[, j]) / 2
    }
  }

  # Then integrate cumulatively over spatial radius.
  cumulative_hu <- matrix(0, nrow = length(h_grid), ncol = length(u_grid))

  if (length(h_grid) > 1L) {
    for (i in 2:length(h_grid)) {
      dh <- h_grid[i] - h_grid[i - 1L]
      cumulative_hu[i, ] <- cumulative_hu[i - 1L, ] +
        dh * (cumulative_u[i - 1L, ] + cumulative_u[i, ]) / 2
    }
  }

  area <- .pp_lgcp_unit_sphere_area(d_space)
  out <- 2 * area * cumulative_hu[match(r, h_grid), match(u, u_grid), drop = FALSE]

  dimnames(out) <- list(as.character(r), as.character(u))
  out
}


.pp_lgcp_theoretical_k <- function(
    empirical_K,
    params,
    cov_model,
    covariance_function,
    integration_n
) {
  if (identical(empirical_K$type, "spatial")) {
    return(.pp_lgcp_theoretical_k_spatial(
      r = empirical_K$r,
      d_space = empirical_K$spatial_dimension,
      params = params,
      cov_model = cov_model,
      covariance_function = covariance_function,
      integration_n = integration_n
    ))
  }

  .pp_lgcp_theoretical_k_st(
    r = empirical_K$r,
    u = empirical_K$u,
    d_space = empirical_K$spatial_dimension,
    params = params,
    cov_model = cov_model,
    covariance_function = covariance_function,
    integration_n = integration_n
  )
}


# Minimum-contrast estimation -------------------------------------------------

.pp_lgcp_trapezoid_weights <- function(x) {
  x <- .pp_lgcp_validate_lags(x, "lag grid")
  n <- length(x)

  if (n == 1L) {
    return(1)
  }

  w <- numeric(n)
  w[1L] <- (x[2L] - x[1L]) / 2
  w[n] <- (x[n] - x[n - 1L]) / 2

  if (n > 2L) {
    w[2:(n - 1L)] <- (x[3:n] - x[1:(n - 2L)]) / 2
  }

  w
}


.pp_lgcp_resolve_contrast_weight <- function(
    contrast_weight,
    empirical_K,
    theoretical_K
) {
  Khat <- empirical_K$Khat

  user_weight <- if (is.null(contrast_weight)) {
    if (is.matrix(Khat)) matrix(1, nrow(Khat), ncol(Khat)) else rep(1, length(Khat))
  } else if (is.function(contrast_weight)) {
    .pp_lgcp_call_user_function(
      fun = contrast_weight,
      supplied_args = list(
        r = empirical_K$r,
        u = empirical_K$u,
        Khat = Khat,
        Kmodel = theoretical_K
      ),
      fun_name = "contrast_weight"
    )
  } else {
    contrast_weight
  }

  if (!is.numeric(user_weight) || length(user_weight) != length(Khat)) {
    stop(
      "`contrast_weight` must be NULL, a function, or numeric values aligned with the empirical K-function.",
      call. = FALSE
    )
  }

  if (is.matrix(Khat)) {
    dim(user_weight) <- dim(Khat)
  }

  if (any(!is.finite(user_weight)) || any(user_weight < 0)) {
    stop("Contrast weights must be finite and non-negative.", call. = FALSE)
  }

  geometry_weight <- if (is.matrix(Khat)) {
    outer(
      .pp_lgcp_trapezoid_weights(empirical_K$r),
      .pp_lgcp_trapezoid_weights(empirical_K$u),
      FUN = "*"
    )
  } else {
    .pp_lgcp_trapezoid_weights(empirical_K$r)
  }

  weight <- user_weight * geometry_weight

  if (!any(weight > 0)) {
    stop("The combined contrast weights must contain at least one positive value.", call. = FALSE)
  }

  weight
}


.pp_lgcp_default_start <- function(window, cov_model) {
  widths <- .pp_lgcp_window_widths(window)
  axes <- .pp_lgcp_axis_sets(window$axis_info)
  spatial_scale <- min(widths[axes$space]) / 5

  switch(
    cov_model,
    exponential = c(sigma2 = 1, alpha = spatial_scale),
    gaussian = c(sigma2 = 1, alpha = spatial_scale),
    matern = c(sigma2 = 1, alpha = spatial_scale, nu = 0.5),
    separable_exponential = c(
      sigma2 = 1,
      alpha_space = spatial_scale,
      alpha_time = widths[axes$time] / 5
    ),
    gneiting = c(
      sigma2 = 1,
      alpha_space = spatial_scale,
      alpha_time = widths[axes$time] / 5,
      beta = 0.5,
      gamma = 0.5
    ),
    custom = stop(
      "A custom covariance model requires explicit named `start` values.",
      call. = FALSE
    ),
    stop("Internal error: unsupported covariance model.", call. = FALSE)
  )
}


.pp_lgcp_default_bounds <- function(window, cov_model) {
  widths <- .pp_lgcp_window_widths(window)
  axes <- .pp_lgcp_axis_sets(window$axis_info)
  spatial_min <- min(widths[axes$space])
  spatial_max <- max(widths[axes$space])

  tiny_space <- max(.Machine$double.eps^0.25, spatial_min * 1e-5)
  max_space <- spatial_max * 10

  if (axes$is_st) {
    time_width <- widths[axes$time]
    tiny_time <- max(.Machine$double.eps^0.25, time_width * 1e-5)
    max_time <- time_width * 10
  }

  lower <- switch(
    cov_model,
    exponential = c(sigma2 = 1e-6, alpha = tiny_space),
    gaussian = c(sigma2 = 1e-6, alpha = tiny_space),
    matern = c(sigma2 = 1e-6, alpha = tiny_space, nu = 0.05),
    separable_exponential = c(
      sigma2 = 1e-6,
      alpha_space = tiny_space,
      alpha_time = tiny_time
    ),
    gneiting = c(
      sigma2 = 1e-6,
      alpha_space = tiny_space,
      alpha_time = tiny_time,
      beta = 0.05,
      gamma = 0.05
    ),
    custom = NULL
  )

  upper <- switch(
    cov_model,
    exponential = c(sigma2 = 20, alpha = max_space),
    gaussian = c(sigma2 = 20, alpha = max_space),
    matern = c(sigma2 = 20, alpha = max_space, nu = 5),
    separable_exponential = c(
      sigma2 = 20,
      alpha_space = max_space,
      alpha_time = max_time
    ),
    gneiting = c(
      sigma2 = 20,
      alpha_space = max_space,
      alpha_time = max_time,
      beta = 1,
      gamma = 1
    ),
    custom = NULL
  )

  list(lower = lower, upper = upper)
}


.pp_lgcp_merge_named_parameters <- function(default, user, arg_name, names_required) {
  if (is.null(user)) {
    return(default[names_required])
  }

  user <- .pp_lgcp_as_named_numeric(user, arg_name)
  unknown <- setdiff(names(user), names_required)

  if (length(unknown) > 0L) {
    stop(
      sprintf(
        "`%s` contains unknown parameter(s): %s.",
        arg_name,
        paste(unknown, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  out <- default
  out[names(user)] <- user
  out[names_required]
}


.pp_lgcp_prepare_optim_parameters <- function(
    window,
    cov_model,
    start,
    lower,
    upper
) {
  if (identical(cov_model, "custom")) {
    if (is.null(start)) {
      stop("A custom covariance model requires named `start` values.", call. = FALSE)
    }

    start <- .pp_lgcp_as_named_numeric(start, "start")
    parameter_names <- names(start)

    lower_default <- stats::setNames(rep(-Inf, length(start)), parameter_names)
    upper_default <- stats::setNames(rep(Inf, length(start)), parameter_names)

    lower <- .pp_lgcp_merge_named_parameters(
      default = lower_default,
      user = lower,
      arg_name = "lower",
      names_required = parameter_names
    )

    upper <- .pp_lgcp_merge_named_parameters(
      default = upper_default,
      user = upper,
      arg_name = "upper",
      names_required = parameter_names
    )
  } else {
    parameter_names <- .pp_lgcp_parameter_names(cov_model)
    start_default <- .pp_lgcp_default_start(window, cov_model)
    bounds_default <- .pp_lgcp_default_bounds(window, cov_model)

    start <- .pp_lgcp_merge_named_parameters(
      default = start_default,
      user = start,
      arg_name = "start",
      names_required = parameter_names
    )

    lower <- .pp_lgcp_merge_named_parameters(
      default = bounds_default$lower,
      user = lower,
      arg_name = "lower",
      names_required = parameter_names
    )

    upper <- .pp_lgcp_merge_named_parameters(
      default = bounds_default$upper,
      user = upper,
      arg_name = "upper",
      names_required = parameter_names
    )
  }

  if (any(lower >= upper)) {
    stop("Every `lower` bound must be strictly smaller than its `upper` bound.", call. = FALSE)
  }

  if (any(start < lower) || any(start > upper)) {
    stop("Every starting value must lie inside its optimisation bounds.", call. = FALSE)
  }

  list(start = start, lower = lower, upper = upper)
}


.pp_lgcp_mincontrast <- function(
    empirical_K,
    window,
    cov_model,
    covariance_function,
    start,
    lower,
    upper,
    contrast_q,
    contrast_p,
    contrast_weight,
    integration_n,
    optim_control,
    verbose
) {
  contrast_q <- .pp_lgcp_validate_positive_number(
    contrast_q,
    "contrast_q",
    lower = 0
  )
  contrast_p <- .pp_lgcp_validate_positive_number(
    contrast_p,
    "contrast_p",
    lower = 0
  )

  if (!is.list(optim_control)) {
    stop("`optim_control` must be a list.", call. = FALSE)
  }

  parameters <- .pp_lgcp_prepare_optim_parameters(
    window = window,
    cov_model = cov_model,
    start = start,
    lower = lower,
    upper = upper
  )

  parameter_names <- names(parameters$start)
  empirical_values <- empirical_K$Khat

  objective <- function(par) {
    names(par) <- parameter_names

    theoretical <- tryCatch(
      .pp_lgcp_theoretical_k(
        empirical_K = empirical_K,
        params = par,
        cov_model = cov_model,
        covariance_function = covariance_function,
        integration_n = integration_n
      ),
      error = function(error) NULL
    )

    if (is.null(theoretical) || any(!is.finite(theoretical)) || any(theoretical < 0)) {
      return(.Machine$double.xmax^0.25)
    }

    weight <- tryCatch(
      .pp_lgcp_resolve_contrast_weight(
        contrast_weight = contrast_weight,
        empirical_K = empirical_K,
        theoretical_K = theoretical
      ),
      error = function(error) NULL
    )

    if (is.null(weight)) {
      return(.Machine$double.xmax^0.25)
    }

    discrepancy <- abs(
      pmax(empirical_values, 0)^contrast_q -
        pmax(theoretical, 0)^contrast_q
    )^contrast_p

    valid <- is.finite(discrepancy) & is.finite(weight) & weight > 0

    if (!any(valid)) {
      return(.Machine$double.xmax^0.25)
    }

    sum(weight[valid] * discrepancy[valid]) / sum(weight[valid])
  }

  if (isTRUE(verbose)) {
    message(
      sprintf(
        "Estimating LGCP second-order parameters by minimum contrast (%s covariance).",
        cov_model
      )
    )
  }

  # L-BFGS-B needs finite box constraints: with unbounded parameters (the
  # default for a custom covariance) its first quasi-Newton step can leave the
  # valid region and the search stops at the starting values with a spurious
  # convergence code. Nelder-Mead is used when no bound is finite.
  if (any(is.finite(parameters$lower)) || any(is.finite(parameters$upper))) {
    optimisation <- stats::optim(
      par = parameters$start,
      fn = objective,
      method = "L-BFGS-B",
      lower = parameters$lower,
      upper = parameters$upper,
      control = optim_control
    )
  } else {
    optimisation <- stats::optim(
      par = parameters$start,
      fn = objective,
      method = "Nelder-Mead",
      control = optim_control
    )
  }

  theta <- as.numeric(optimisation$par)
  names(theta) <- parameter_names

  theoretical <- .pp_lgcp_theoretical_k(
    empirical_K = empirical_K,
    params = theta,
    cov_model = cov_model,
    covariance_function = covariance_function,
    integration_n = integration_n
  )

  weight <- .pp_lgcp_resolve_contrast_weight(
    contrast_weight = contrast_weight,
    empirical_K = empirical_K,
    theoretical_K = theoretical
  )

  list(
    theta = theta,
    theoretical_K = theoretical,
    contrast_value = objective(theta),
    weight = weight,
    optimisation = optimisation,
    start = parameters$start,
    lower = parameters$lower,
    upper = parameters$upper,
    contrast_q = contrast_q,
    contrast_p = contrast_p
  )
}


# Internal LGCP fitting engine. The documented user interface is
# fit_lgcp_pp() in ppinterface.R, which forwards first-order settings to
# fit_pp() and second-order settings from lgcp_control().
.pp_fit_lgcp_engine <- function(
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
    keep_cubature = c("summary", "full", "none"),
    cov_model = NULL,
    covariance_function = NULL,
    start = NULL,
    lower = NULL,
    upper = NULL,
    r = NULL,
    u = NULL,
    correction = c("translate", "none"),
    contrast_q = 1 / 4,
    contrast_p = 2,
    contrast_weight = NULL,
    integration_n = NULL,
    optim_control = list(maxit = 200L),
    k_chunk_size = 2000L,
    mark_n = NULL,
    max_ground_mark_combinations = 100000L,
    verbose = FALSE
) {
  .pp_lgcp_require(c("fit_pp", "K_fit_pp"), "fit_lgcp_pp()")
  .pp_lgcp_validate_pattern(pattern)

  verbose <- .pp_lgcp_validate_logical(verbose, "verbose")
  seed <- .pp_lgcp_validate_seed(seed)
  correction <- match.arg(correction)
  dummy_layout <- match.arg(dummy_layout)
  keep_cubature <- match.arg(keep_cubature)
  k_chunk_size <- .pp_lgcp_validate_positive_integer(
    k_chunk_size,
    "k_chunk_size",
    min_value = 1L
  )

  axes <- .pp_lgcp_axis_sets(pattern$axis_info)

  if (!axes$is_st && !is.null(u)) {
    stop("`u` can be supplied only for a spatio-temporal LGCP.", call. = FALSE)
  }

  covariance_spec <- .pp_lgcp_normalize_cov_model(
    cov_model = cov_model,
    is_st = axes$is_st,
    covariance_function = covariance_function
  )

  if (is.null(integration_n)) {
    integration_n <- if (axes$is_st) c(100L, 100L) else 400L
  }

  if (!is.null(event_weights) && is.numeric(event_weights) &&
      any(event_weights != 1, na.rm = TRUE)) {
    warning(
      paste0(
        "`event_weights` affect only the first-order `fit_pp()` stage. ",
        "The empirical second-order K-function is calculated from the full ",
        "observed point pattern without event weighting."
      ),
      call. = FALSE
    )
  }

  if (isTRUE(verbose)) {
    message("Fitting marginal first-order intensity with fit_pp().")
  }

  first_order <- fit_pp(
    pattern = pattern,
    formula = formula,
    covariates = covariates,
    covariate_methods = covariate_methods,
    covariate_args = covariate_args,
    continuous_marks = continuous_marks,
    categorical_marks = categorical_marks,
    categorical_levels = categorical_levels,
    event_weights = event_weights,
    q = q,
    dummy_layout = dummy_layout,
    nc = nc,
    max_tiles = max_tiles,
    max_quadrature_points = max_quadrature_points,
    max_mark_combinations = max_mark_combinations,
    seed = seed,
    gam_args = gam_args,
    keep_cubature = keep_cubature
  )

  if (isTRUE(verbose)) {
    message("Computing empirical inhomogeneous K-function from the fitted ground intensity.")
  }

  empirical_K <- K_fit_pp(
    fit = first_order,
    pattern = pattern,
    r = r,
    u = u,
    correction = correction,
    local = FALSE,
    chunk_size = k_chunk_size,
    mark_n = mark_n,
    max_mark_combinations = max_ground_mark_combinations
  )

  second_order <- .pp_lgcp_mincontrast(
    empirical_K = empirical_K,
    window = pattern$window,
    cov_model = covariance_spec$model,
    covariance_function = covariance_spec$covariance_function,
    start = start,
    lower = lower,
    upper = upper,
    contrast_q = contrast_q,
    contrast_p = contrast_p,
    contrast_weight = contrast_weight,
    integration_n = integration_n,
    optim_control = optim_control,
    verbose = verbose
  )

  theoretical_K <- empirical_K
  theoretical_K$Khat <- second_order$theoretical_K
  theoretical_K$residual <- empirical_K$Khat - second_order$theoretical_K
  theoretical_K$Ktheo <- second_order$theoretical_K
  theoretical_K$source <- "fitted_lgcp"

  model_spec <- list(
    covariance_model = covariance_spec$model,
    covariance_function = covariance_spec$covariance_function,
    covariance_params = second_order$theta,
    centred = TRUE,
    spatial_axes = axes$space,
    time_axis = if (axes$is_st) axes$time else NULL,
    spatial_dimension = axes$d_space,
    window = pattern$window,
    integration_n = integration_n,
    contrast_q = second_order$contrast_q,
    contrast_p = second_order$contrast_p,
    correction = correction
  )

  out <- list(
    call = match.call(),
    first_order = first_order,
    theta = second_order$theta,
    covariance_model = covariance_spec$model,
    covariance_function = covariance_spec$covariance_function,
    empirical_K = empirical_K,
    theoretical_K = theoretical_K,
    second_order = second_order,
    pattern = pattern,
    model_spec = model_spec,
    settings = list(
      second_method = "minimum_contrast_K",
      r = empirical_K$r,
      u = empirical_K$u,
      correction = correction,
      contrast_q = second_order$contrast_q,
      contrast_p = second_order$contrast_p,
      integration_n = integration_n,
      k_chunk_size = k_chunk_size,
      mark_n = mark_n,
      max_ground_mark_combinations = max_ground_mark_combinations
    )
  )

  class(out) <- "pp_lgcp_fit"
  out
}


# Simulation from a fitted LGCP ----------------------------------------------

.pp_lgcp_validate_fit <- function(object) {
  if (!inherits(object, "pp_lgcp_fit")) {
    stop("`object` must inherit from class 'pp_lgcp_fit'.", call. = FALSE)
  }

  if (is.null(object$first_order) || !inherits(object$first_order, "pp_fit")) {
    stop("The LGCP object does not contain a valid first-order `pp_fit`.", call. = FALSE)
  }

  if (is.null(object$model_spec) || !is.list(object$model_spec)) {
    stop("The LGCP object is missing its retained model specification.", call. = FALSE)
  }

  invisible(object)
}


.pp_lgcp_evaluate_field <- function(field, coords, field_method, field_args) {
  covariate <- .pp_lgcp_field_covariate(field)

  args <- c(
    list(
      covariate = covariate,
      at = coords,
      method = field_method
    ),
    field_args
  )

  value <- tryCatch(
    do.call(evaluate_covariate, args),
    error = function(error) {
      stop(
        sprintf("Latent-field interpolation failed: %s", conditionMessage(error)),
        call. = FALSE
      )
    }
  )

  as.numeric(value)
}


.pp_lgcp_simulate_fit_one <- function(
    object,
    field,
    field_method,
    field_args,
    bound,
    lambda_upper,
    bound_n,
    bound_max_points,
    bound_safety,
    on_bound_violation,
    restart_safety,
    max_restarts,
    max_candidates,
    output,
    keep_field
) {
  fit <- object$first_order
  domain <- .pp_diag_joint_domain(fit)
  variance <- field$variance

  if (identical(bound, "user")) {
    if (!is.numeric(lambda_upper) || length(lambda_upper) != 1L ||
        !is.finite(lambda_upper) || lambda_upper < 0) {
      stop(
        "`lambda_upper` must be one finite non-negative number when `bound = 'user'`.",
        call. = FALSE
      )
    }

    initial_bound <- as.numeric(lambda_upper)
    bound_details <- list(
      method = "user",
      lambda_upper_initial = initial_bound,
      base_lambda_upper = NA_real_,
      latent_factor_bound = NA_real_
    )
  } else {
    grid_bound <- .pp_diag_joint_bound_grid(
      fit = fit,
      bound_n = bound_n,
      bound_max_points = bound_max_points,
      bound_safety = bound_safety
    )

    # For nearest-neighbour and IDW interpolation, interpolated values cannot
    # exceed the maximum support value. This bounds the latent multiplier; any
    # remaining first-order grid-bound error is still caught by restart logic.
    latent_factor_bound <- exp(max(field$value) - variance / 2)
    initial_bound <- grid_bound$lambda_upper * latent_factor_bound

    bound_details <- c(
      list(
        method = "grid",
        lambda_upper_initial = initial_bound,
        base_lambda_upper = grid_bound$lambda_upper,
        latent_factor_bound = latent_factor_bound
      ),
      grid_bound[-1L]
    )
  }

  current_bound <- initial_bound
  restarts <- 0L

  repeat {
    expected_candidates <- current_bound * domain$augmented_volume

    if (!is.finite(expected_candidates) || expected_candidates > max_candidates) {
      stop(
        sprintf(
          paste0(
            "The fitted LGCP thinning process would require an expected %.0f candidates, ",
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
          "The fitted LGCP generated %d candidates, exceeding `max_candidates = %d`.",
          n_candidates,
          max_candidates
        ),
        call. = FALSE
      )
    }

    candidates <- .pp_diag_sample_joint_candidates(domain, n_candidates)

    base_lambda <- .pp_diag_predict_joint(
      fit = fit,
      coords = as.matrix(candidates$coords),
      marks = candidates$marks
    )

    latent <- .pp_lgcp_evaluate_field(
      field = field,
      coords = candidates$coords,
      field_method = field_method,
      field_args = field_args
    )

    lambda_candidate <- base_lambda * exp(latent - variance / 2)

    tolerance <- .Machine$double.eps^0.5 * max(1, current_bound)
    violation <- lambda_candidate > current_bound + tolerance

    if (any(violation)) {
      observed_max <- max(lambda_candidate[violation])

      if (identical(on_bound_violation, "error")) {
        stop(
          sprintf(
            paste0(
              "The fitted LGCP thinning bound was violated: %.8g exceeds ",
              "lambda_upper = %.8g. Supply a larger verified bound or use ",
              "`on_bound_violation = 'restart'`."
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
          "The fitted LGCP thinning bound was violated more often than `max_restarts`.",
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
      mode = "joint_thinning_from_pp_lgcp_fit",
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

    pattern <- .pp_lgcp_attach_simulation_metadata(
      x = pattern,
      field = field,
      field_method = field_method,
      field_args = field_args,
      keep_field = keep_field,
      source = "simulate.pp_lgcp_fit"
    )

    if (identical(output, "pattern")) {
      return(pattern)
    }

    out <- c(
      metadata,
      list(
        pattern = pattern,
        lgcp = pattern$metadata$lgcp,
        latent_candidate = latent,
        base_intensity = base_lambda,
        conditional_intensity = lambda_candidate
      )
    )

    class(out) <- "pp_lgcp_simulation"
    return(out)
  }
}


#' Simulate from a fitted log-Gaussian Cox process
#'
#' Generates independent realisations from a `pp_lgcp_fit`. Each simulation
#' draws a new latent Gaussian field using the estimated covariance parameters,
#' multiplies the fitted first-order joint intensity by the centred latent
#' factor, and performs joint-domain thinning. Fitted continuous joint marks,
#' categorical marks, and field marks are therefore simulated consistently with
#' the retained first-order `pp_fit` specification.
#'
#' @param object A `pp_lgcp_fit` object created by `fit_lgcp_pp()`.
#' @param nsim Number of independent simulations.
#' @param seed Optional integer random seed.
#' @param field_support,field_n,field_method,field_args,max_field_points,jitter
#'   Latent-field support, interpolation, size, and jitter settings as in
#'   `simulate_lgcp_pp()`.
#' @param max_jitter_attempts,keep_field Additional latent-field factorisation
#'   and retention settings as in `simulate_lgcp_pp()`.
#' @param bound Either `"grid"` or `"user"`.
#' @param lambda_upper User-supplied upper bound for the full conditional joint
#'   LGCP intensity when `bound = "user"`.
#' @param bound_n,bound_max_points,bound_safety Grid-bound settings applied to
#'   the fitted first-order joint intensity.
#' @param on_bound_violation Either `"restart"` or `"error"`.
#' @param restart_safety,max_restarts,max_candidates Thinning restart and safety
#'   settings.
#' @param output Either `"pattern"` or `"diagnostic"`.
#' @param ... Ignored.
#'
#' @return A simulated `pp` pattern or `pp_lgcp_simulation` diagnostic object;
#'   for `nsim > 1`, a named list of independent outputs.
#'
#' @references
#' Møller, J., Syversveen, A. R., and Waagepetersen, R. P. (1998).
#' Log Gaussian Cox processes. *Scandinavian Journal of Statistics*, 25,
#' 451--482. doi:10.1111/1467-9469.00115.
#'
#' Lewis, P. A. W. and Shedler, G. S. (1979). Simulation of nonhomogeneous
#' Poisson processes by thinning. *Naval Research Logistics Quarterly*, 26,
#' 403--413.
#'
#' @seealso `fit_lgcp_pp()`, `simulate_lgcp_pp()`, `simulate_pp_fit()`
#'
#' @examples
#' \donttest{
#' w <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#'
#' x <- simulate_lgcp_pp(
#'   window = w,
#'   expected_n = 60,
#'   cov_model = "exponential",
#'   covariance_params = c(sigma2 = 0.3, alpha = 0.15),
#'   field_n = 6,
#'   output = "pattern",
#'   seed = 10
#' )
#'
#' fit <- fit_lgcp_pp(
#'   pattern = x,
#'   formula = ~ 1,
#'   cov_model = "exponential",
#'   r = seq(0, 0.2, length.out = 8),
#'   nc = 6,
#'   control = lgcp_control(integration_n = 80)
#' )
#'
#' x_new <- simulate(fit, nsim = 1, field_n = 6, seed = 12)
#' summary(x_new)
#' }
#'
#' @method simulate pp_lgcp_fit
#' @export
simulate.pp_lgcp_fit <- function(
    object,
    nsim = 1L,
    seed = NULL,
    field_support = NULL,
    field_n = NULL,
    field_method = c("nearest", "idw"),
    field_args = list(),
    max_field_points = 2000L,
    jitter = 1e-10,
    max_jitter_attempts = 8L,
    keep_field = TRUE,
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
    ...
) {
  .pp_lgcp_require(
    c(
      ".pp_diag_joint_domain", ".pp_diag_joint_bound_grid",
      ".pp_diag_sample_joint_candidates", ".pp_diag_predict_joint",
      ".pp_diag_build_simulated_pattern", "evaluate_covariate"
    ),
    "simulate.pp_lgcp_fit()"
  )

  .pp_lgcp_validate_fit(object)

  nsim <- .pp_lgcp_validate_positive_integer(nsim, "nsim", min_value = 1L)
  seed <- .pp_lgcp_validate_seed(seed)
  field_method <- .pp_lgcp_validate_field_method(match.arg(field_method))
  keep_field <- .pp_lgcp_validate_logical(keep_field, "keep_field")
  bound <- match.arg(bound)
  on_bound_violation <- match.arg(on_bound_violation)
  output <- match.arg(output)

  if (is.null(field_args)) {
    field_args <- list()
  }

  if (identical(field_method, "idw")) {
    warning(
      paste0(
        "IDW interpolation changes the marginal variance of the latent ",
        "Gaussian field. Consequently, centering by C(0) / 2 may not ",
        "preserve the requested marginal first-order intensity exactly. ",
        "Use `field_method = 'nearest'` for the default variance-preserving ",
        "finite-support approximation."
      ),
      call. = FALSE
    )
  }

  if (!is.list(field_args)) {
    stop("`field_args` must be a list or NULL.", call. = FALSE)
  }

  bound_max_points <- .pp_lgcp_validate_positive_integer(
    bound_max_points,
    "bound_max_points",
    min_value = 2L
  )
  max_restarts <- .pp_lgcp_validate_positive_integer(
    max_restarts,
    "max_restarts",
    min_value = 1L
  )
  max_candidates <- .pp_lgcp_validate_positive_integer(
    max_candidates,
    "max_candidates",
    min_value = 1L
  )

  bound_safety <- .pp_lgcp_validate_positive_number(
    bound_safety,
    "bound_safety",
    lower = 1,
    inclusive = TRUE
  )
  restart_safety <- .pp_lgcp_validate_positive_number(
    restart_safety,
    "restart_safety",
    lower = 1,
    inclusive = FALSE
  )

  if (!is.null(seed)) {
    set.seed(seed)
  }

  results <- vector("list", nsim)
  spec <- object$model_spec

  for (i in seq_len(nsim)) {
    field <- .pp_lgcp_simulate_field(
      window = spec$window,
      cov_model = spec$covariance_model,
      covariance_params = spec$covariance_params,
      covariance_function = spec$covariance_function,
      field_support = field_support,
      field_n = field_n,
      max_field_points = max_field_points,
      jitter = jitter,
      max_jitter_attempts = max_jitter_attempts
    )

    results[[i]] <- .pp_lgcp_simulate_fit_one(
      object = object,
      field = field,
      field_method = field_method,
      field_args = field_args,
      bound = bound,
      lambda_upper = lambda_upper,
      bound_n = bound_n,
      bound_max_points = bound_max_points,
      bound_safety = bound_safety,
      on_bound_violation = on_bound_violation,
      restart_safety = restart_safety,
      max_restarts = max_restarts,
      max_candidates = max_candidates,
      output = output,
      keep_field = keep_field
    )
  }

  if (nsim == 1L) {
    return(results[[1L]])
  }

  names(results) <- paste0("sim_", seq_len(nsim))
  results
}


# Printing, summaries, and plotting ------------------------------------------

#' @rdname fit_lgcp_pp
#' @param x A `pp_lgcp_fit` object.
#' @param digits Number of significant digits used for printing.
#' @param ... Additional arguments passed to methods.
#' @method print pp_lgcp_fit
#' @export
print.pp_lgcp_fit <- function(x, digits = 4, ...) {
  .pp_lgcp_validate_fit(x)

  axes <- .pp_lgcp_axis_sets(x$pattern$axis_info)

  cat("<pp_lgcp_fit>\n")
  cat("Log-Gaussian Cox process fit\n\n")

  cat(
    "Process: ",
    if (axes$is_st) "spatio-temporal" else "spatial",
    "\n",
    sep = ""
  )
  cat("Spatial dimension: ", axes$d_space, "\n", sep = "")
  if (axes$is_st) cat("Time axis: ", axes$time, "\n", sep = "")
  cat("Observed points: ", nrow(x$pattern$coords), "\n", sep = "")
  cat("Covariance model: ", x$covariance_model, "\n", sep = "")
  cat("Second-order method: minimum contrast on K\n\n")

  cat("Estimated covariance parameters:\n")
  print(round(x$theta, digits = digits))

  cat(
    "\nMinimum-contrast value: ",
    format(x$second_order$contrast_value, digits = digits),
    "\n",
    sep = ""
  )
  cat("Optimizer convergence code: ", x$second_order$optimisation$convergence, "\n", sep = "")

  invisible(x)
}


#' @rdname fit_lgcp_pp
#' @param object A `pp_lgcp_fit` object.
#' @method summary pp_lgcp_fit
#' @export
summary.pp_lgcp_fit <- function(object, ...) {
  .pp_lgcp_validate_fit(object)

  axes <- .pp_lgcp_axis_sets(object$pattern$axis_info)

  out <- list(
    call = object$call,
    process_type = if (axes$is_st) "spatio-temporal" else "spatial",
    spatial_dimension = axes$d_space,
    spatial_axes = axes$space,
    time_axis = if (axes$is_st) axes$time else NULL,
    n_points = nrow(object$pattern$coords),
    covariance_model = object$covariance_model,
    theta = object$theta,
    contrast_value = object$second_order$contrast_value,
    contrast_q = object$second_order$contrast_q,
    contrast_p = object$second_order$contrast_p,
    convergence = object$second_order$optimisation$convergence,
    convergence_message = object$second_order$optimisation$message,
    first_order_coefficients = object$first_order$coefficients,
    first_order_loglik_approx = object$first_order$loglik_approx,
    correction = object$settings$correction,
    r = object$empirical_K$r,
    u = object$empirical_K$u
  )

  class(out) <- "summary_pp_lgcp_fit"
  out
}


#' @rdname fit_lgcp_pp
#' @method print summary_pp_lgcp_fit
#' @export
print.summary_pp_lgcp_fit <- function(x, digits = 4, ...) {
  cat("Summary of <pp_lgcp_fit>\n\n")
  cat("Process: ", x$process_type, "\n", sep = "")
  cat("Spatial dimension: ", x$spatial_dimension, "\n", sep = "")
  cat("Observed points: ", x$n_points, "\n", sep = "")
  cat("Covariance model: ", x$covariance_model, "\n\n", sep = "")

  cat("First-order coefficients:\n")
  print(round(x$first_order_coefficients, digits = digits))

  cat("\nSecond-order covariance parameters:\n")
  print(round(x$theta, digits = digits))

  cat(
    "\nMinimum-contrast value: ",
    format(x$contrast_value, digits = digits),
    "\n",
    sep = ""
  )
  cat("Optimizer convergence code: ", x$convergence, "\n", sep = "")

  if (!is.null(x$convergence_message) && nzchar(x$convergence_message)) {
    cat("Optimizer message: ", x$convergence_message, "\n", sep = "")
  }

  invisible(x)
}


#' @rdname fit_lgcp_pp
#' @param what Plot type. For a spatial fit, `"comparison"` overlays empirical
#'   and fitted `K`; `"residual"`, `"empirical"`, and `"fitted"` are also
#'   available. For a spatio-temporal fit, `"comparison"` draws empirical and
#'   fitted `K(r,u)` surfaces side by side.
#' @method plot pp_lgcp_fit
#' @export
plot.pp_lgcp_fit <- function(
    x,
    what = c("comparison", "residual", "empirical", "fitted"),
    ...
) {
  .pp_lgcp_validate_fit(x)
  what <- match.arg(what)

  empirical <- x$empirical_K$Khat
  fitted <- x$second_order$theoretical_K

  if (!is.matrix(empirical)) {
    r <- x$empirical_K$r

    if (identical(what, "comparison")) {
      graphics::plot(
        r,
        empirical,
        type = "l",
        xlab = "r",
        ylab = "K(r)",
        main = "Empirical and fitted LGCP K-function",
        ...
      )
      graphics::lines(r, fitted, lty = 2)
      graphics::legend(
        "topleft",
        legend = c("Empirical", "Fitted LGCP"),
        lty = c(1, 2),
        bty = "n"
      )
      return(invisible(x))
    }

    y <- switch(
      what,
      residual = empirical - fitted,
      empirical = empirical,
      fitted = fitted
    )

    ylab <- switch(
      what,
      residual = "Empirical - fitted K(r)",
      empirical = "Empirical K(r)",
      fitted = "Fitted LGCP K(r)"
    )

    graphics::plot(r, y, type = "l", xlab = "r", ylab = ylab, ...)

    if (identical(what, "residual")) {
      graphics::abline(h = 0, lty = 2)
    }

    return(invisible(x))
  }

  r <- x$empirical_K$r
  u <- x$empirical_K$u

  extra_args <- list(...)

  draw_surface <- function(z, main) {
    do.call(
      graphics::image,
      c(
        list(
          x = r,
          y = u,
          z = z,
          xlab = "Spatial lag r",
          ylab = "Temporal lag u",
          main = main
        ),
        extra_args
      )
    )
  }

  if (identical(what, "comparison")) {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)
    graphics::par(mfrow = c(1, 2))
    draw_surface(empirical, "Empirical K(r,u)")
    draw_surface(fitted, "Fitted LGCP K(r,u)")
    return(invisible(x))
  }

  z <- switch(
    what,
    residual = empirical - fitted,
    empirical = empirical,
    fitted = fitted
  )

  title <- switch(
    what,
    residual = "Empirical - fitted K(r,u)",
    empirical = "Empirical K(r,u)",
    fitted = "Fitted LGCP K(r,u)"
  )

  draw_surface(z, title)
  invisible(x)
}


#' @rdname simulate_lgcp_pp
#' @param x A `pp_lgcp_simulation` object.
#' @param digits Number of significant digits used for printing.
#' @param ... Additional arguments passed to methods.
#' @method print pp_lgcp_simulation
#' @export
print.pp_lgcp_simulation <- function(x, digits = 4, ...) {
  cat("<pp_lgcp_simulation>\n")
  cat("Log-Gaussian Cox process simulation\n\n")

  if (!is.null(x$n_candidates)) {
    cat("Candidates: ", x$n_candidates, "\n", sep = "")
  }
  if (!is.null(x$n_accepted)) {
    cat("Accepted points: ", x$n_accepted, "\n", sep = "")
  }
  if (!is.null(x$acceptance_rate)) {
    cat(
      "Acceptance rate: ",
      format(x$acceptance_rate, digits = digits),
      "\n",
      sep = ""
    )
  }

  lgcp <- x$lgcp

  if (!is.null(lgcp)) {
    cat("Covariance model: ", lgcp$covariance_model, "\n", sep = "")
    cat("Latent-field support points: ", lgcp$n_field_points, "\n", sep = "")
  }

  invisible(x)
}


.pp_lgcp_muffle_grid_bound_warning <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("uses a numerical upper-bound estimate", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}
