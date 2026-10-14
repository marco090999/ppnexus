# ETAS catalog simulation ----------------------------------------------------
#
# This module provides a package-oriented implementation of the validated
# `etas.par.sim_v4()` branching simulator used in the parametric-background
# ETAS simulation study. The simulation mechanism is intentionally specialised
# to two spatial coordinates plus time, matching `fit_etas_pp()`.
#
# Required package modules:
#   ppdata.R   -> pp_window(), pp_pattern()
#   ppetas.R   -> etas_catalog()
#
# No compiled etasFLP routine is required for simulation.

.pp_etas_sim_version <- "0.1.0"


# Dependencies and validation ------------------------------------------------

.pp_etas_sim_require <- function(functions, context = "ETAS simulation") {
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
          "Source ppdata.R and ppetas.R before ppetasim.R."
        ),
        context,
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.pp_etas_sim_validate_logical <- function(x, arg_name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", arg_name), call. = FALSE)
  }

  isTRUE(x)
}


.pp_etas_sim_validate_numeric_scalar <- function(
    x,
    arg_name,
    lower = -Inf,
    inclusive = TRUE
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


.pp_etas_sim_validate_positive_integer <- function(x, arg_name, min_value = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != as.integer(x) || x < min_value) {
    stop(
      sprintf(
        "`%s` must be an integer greater than or equal to %d.",
        arg_name,
        min_value
      ),
      call. = FALSE
    )
  }

  as.integer(x)
}


.pp_etas_sim_validate_seed <- function(seed) {
  if (is.null(seed)) {
    return(NULL)
  }

  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("`seed` must be NULL or one finite integer.", call. = FALSE)
  }

  as.integer(seed)
}


.pp_etas_sim_validate_range <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) != 2L || any(!is.finite(x))) {
    stop(
      sprintf("`%s` must be a finite numeric vector of length two.", arg_name),
      call. = FALSE
    )
  }

  x <- as.numeric(x)

  if (x[1L] >= x[2L]) {
    stop(
      sprintf("`%s[1]` must be strictly smaller than `%s[2]`.", arg_name, arg_name),
      call. = FALSE
    )
  }

  x
}


.pp_etas_sim_validate_params <- function(params) {
  required <- c("mu", "k0", "c", "p", "gamma", "d", "q")

  if (!is.numeric(params) || is.null(names(params)) ||
      !all(required %in% names(params))) {
    stop(
      paste0(
        "`params` must be a named numeric vector containing: ",
        paste(required, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  params <- as.numeric(params[required])
  names(params) <- required

  if (any(!is.finite(params))) {
    stop("All ETAS parameters must be finite.", call. = FALSE)
  }

  if (params[["mu"]] <= 0) stop("`mu` must be > 0.", call. = FALSE)
  if (params[["k0"]] <= 0) stop("`k0` must be > 0.", call. = FALSE)
  if (params[["c"]] <= 0) stop("`c` must be > 0.", call. = FALSE)
  if (params[["p"]] <= 1) stop("`p` must be > 1.", call. = FALSE)
  if (params[["d"]] <= 0) stop("`d` must be > 0.", call. = FALSE)
  if (params[["q"]] <= 1) stop("`q` must be > 1.", call. = FALSE)

  params
}


.pp_etas_sim_validate_linear_coefs <- function(coefs, what) {
  if (is.null(coefs) || !is.numeric(coefs) || is.null(names(coefs)) ||
      any(names(coefs) == "") || anyDuplicated(names(coefs)) ||
      any(!is.finite(coefs))) {
    stop(
      sprintf("%s: coefficients must be a finite named numeric vector.", what),
      call. = FALSE
    )
  }

  as.numeric_coefs <- as.numeric(coefs)
  names(as.numeric_coefs) <- names(coefs)
  as.numeric_coefs
}


# Feature construction and linear predictors --------------------------------

.pp_etas_sim_clamp <- function(x, limit = 20) {
  pmax(pmin(x, limit), -limit)
}


.pp_etas_sim_scale_centered <- function(x, xmin, xmax) {
  (x - (xmin + xmax) / 2) / (xmax - xmin)
}


.pp_etas_sim_build_features <- function(x, y, window_working, m_rel = NULL) {
  xmin <- window_working[["xmin"]]
  xmax <- window_working[["xmax"]]
  ymin <- window_working[["ymin"]]
  ymax <- window_working[["ymax"]]

  xc <- (xmin + xmax) / 2
  yc <- (ymin + ymax) / 2

  out <- data.frame(
    x = as.numeric(x),
    y = as.numeric(y),
    x_std = .pp_etas_sim_scale_centered(x, xmin, xmax),
    y_std = .pp_etas_sim_scale_centered(y, ymin, ymax),
    r2_center = (x - xc)^2 + (y - yc)^2,
    check.names = FALSE
  )

  if (!is.null(m_rel)) {
    out$m_rel <- as.numeric(m_rel)
  }

  out
}


.pp_etas_sim_eval_linear_lp <- function(df, coefs, what = "linear LP") {
  coefs <- .pp_etas_sim_validate_linear_coefs(coefs, what)
  missing <- setdiff(names(coefs), names(df))

  if (length(missing) > 0L) {
    stop(
      sprintf(
        "%s: missing variables in feature data frame: %s.",
        what,
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  as.numeric(
    as.matrix(df[, names(coefs), drop = FALSE]) %*% as.numeric(coefs)
  )
}


.pp_etas_sim_eval_custom_lp <- function(fun, df, what = "custom LP") {
  if (!is.function(fun)) {
    stop(sprintf("%s: a function is required.", what), call. = FALSE)
  }

  value <- fun(df)

  if (!is.numeric(value)) {
    stop(sprintf("%s: the function must return numeric values.", what), call. = FALSE)
  }

  if (length(value) == 1L) {
    value <- rep(value, nrow(df))
  }

  if (length(value) != nrow(df) || any(!is.finite(value))) {
    stop(
      sprintf("%s: the function must return one finite value per row.", what),
      call. = FALSE
    )
  }

  as.numeric(value)
}


.pp_etas_sim_append_spatial_covariates <- function(
    df_features,
    bg_spatial_cov_funs = NULL,
    bg_spatial_cov_generator = NULL,
    window_working,
    seed = NULL,
    stage = c("support", "events")
) {
  stage <- match.arg(stage)
  out <- df_features

  if (!is.null(bg_spatial_cov_generator)) {
    if (!is.function(bg_spatial_cov_generator)) {
      stop(
        paste0(
          "`bg_spatial_cov_generator` must be a function with signature ",
          "function(df_features, window_km, seed = NULL)."
        ),
        call. = FALSE
      )
    }

    generated <- bg_spatial_cov_generator(
      df_features = df_features,
      window_km = window_working,
      seed = seed
    )

    if (!is.data.frame(generated)) {
      stop("`bg_spatial_cov_generator` must return a data frame.", call. = FALSE)
    }

    if (nrow(generated) != nrow(df_features)) {
      stop("`bg_spatial_cov_generator` returned the wrong number of rows.", call. = FALSE)
    }

    if (is.null(names(generated)) || any(names(generated) == "") ||
        anyDuplicated(names(generated))) {
      stop(
        "`bg_spatial_cov_generator` must return uniquely named columns.",
        call. = FALSE
      )
    }

    overlap <- intersect(names(generated), names(out))

    if (length(overlap) > 0L) {
      warning(
        paste0(
          "`bg_spatial_cov_generator` overwrites columns: ",
          paste(overlap, collapse = ", "),
          "."
        ),
        call. = FALSE
      )
      out[, overlap] <- NULL
    }

    out <- cbind(out, generated)
  }

  if (!is.null(bg_spatial_cov_funs)) {
    if (!is.list(bg_spatial_cov_funs) || is.null(names(bg_spatial_cov_funs)) ||
        any(names(bg_spatial_cov_funs) == "") || anyDuplicated(names(bg_spatial_cov_funs))) {
      stop(
        "`bg_spatial_cov_funs` must be a uniquely named list of functions.",
        call. = FALSE
      )
    }

    for (name_i in names(bg_spatial_cov_funs)) {
      fun_i <- bg_spatial_cov_funs[[name_i]]

      if (!is.function(fun_i)) {
        stop(
          sprintf("`bg_spatial_cov_funs[['%s']]` is not a function.", name_i),
          call. = FALSE
        )
      }

      value_i <- fun_i(out)

      if (!is.numeric(value_i)) {
        stop(
          sprintf("Spatial covariate function `%s` must return numeric values.", name_i),
          call. = FALSE
        )
      }

      if (length(value_i) == 1L) {
        value_i <- rep(value_i, nrow(out))
      }

      if (length(value_i) != nrow(out) || any(!is.finite(value_i))) {
        stop(
          sprintf(
            "Spatial covariate function `%s` must return one finite value per row.",
            name_i
          ),
          call. = FALSE
        )
      }

      out[[name_i]] <- as.numeric(value_i)
    }
  }

  out
}


.pp_etas_sim_eval_background_lp <- function(
    features,
    bg_lp_type,
    bg_lp_coefs,
    bg_lp_fun,
    clip_eta_bg
) {
  eta <- switch(
    bg_lp_type,
    constant = rep.int(0, nrow(features)),
    linear = .pp_etas_sim_eval_linear_lp(
      features,
      bg_lp_coefs,
      what = "background linear LP"
    ),
    custom = .pp_etas_sim_eval_custom_lp(
      bg_lp_fun,
      features,
      what = "background custom LP"
    )
  )

  .pp_etas_sim_clamp(eta, clip_eta_bg)
}


.pp_etas_sim_eval_trigger_lp <- function(
    features,
    trig_lp_type,
    trig_lp_coefs,
    trig_lp_fun,
    clip_eta_trig
) {
  eta <- switch(
    trig_lp_type,
    linear = .pp_etas_sim_eval_linear_lp(
      features,
      trig_lp_coefs,
      what = "triggering linear LP"
    ),
    custom = .pp_etas_sim_eval_custom_lp(
      trig_lp_fun,
      features,
      what = "triggering custom LP"
    )
  )

  .pp_etas_sim_clamp(eta, clip_eta_trig)
}


# Sampling helpers ------------------------------------------------------------

.pp_etas_sim_sample_xy <- function(n, x0, y0, d, q) {
  if (n == 0L) {
    return(matrix(numeric(0), nrow = 0L, ncol = 2L))
  }

  theta <- stats::runif(n, 0, 2 * pi)
  u <- stats::runif(n)
  radius <- sqrt(d * (u^(1 / (1 - q)) - 1))

  cbind(
    x = x0 + radius * cos(theta),
    y = y0 + radius * sin(theta)
  )
}


.pp_etas_sim_rbind_fill <- function(...) {
  dots <- list(...)
  dots <- dots[!vapply(dots, is.null, logical(1))]

  if (length(dots) == 0L) {
    return(data.frame())
  }

  dots <- lapply(dots, as.data.frame, stringsAsFactors = FALSE)
  dots <- dots[vapply(dots, nrow, integer(1)) > 0L]

  if (length(dots) == 0L) {
    return(data.frame())
  }

  all_names <- unique(unlist(lapply(dots, names), use.names = FALSE))

  dots <- lapply(
    dots,
    function(data_i) {
      missing <- setdiff(all_names, names(data_i))

      if (length(missing) > 0L) {
        for (name_i in missing) {
          data_i[[name_i]] <- NA
        }
      }

      data_i[, all_names, drop = FALSE]
    }
  )

  out <- do.call(base::rbind, dots)
  rownames(out) <- NULL
  out
}


.pp_etas_sim_empty_catalog_data <- function(extra_covariate_names = character(0)) {
  out <- data.frame(
    event_id = integer(0),
    father_id = integer(0),
    lgen = integer(0),
    time = numeric(0),
    lat = numeric(0),
    long = numeric(0),
    z = numeric(0),
    magn1 = numeric(0),
    x_km = numeric(0),
    y_km = numeric(0),
    m_rel = numeric(0),
    x_true = numeric(0),
    y_true = numeric(0),
    x_std_true = numeric(0),
    y_std_true = numeric(0),
    r2_center_true = numeric(0),
    eta_bg_true = numeric(0),
    eta_trig_true = numeric(0),
    check.names = FALSE
  )

  if (length(extra_covariate_names) > 0L) {
    for (name_i in extra_covariate_names) {
      out[[name_i]] <- numeric(0)
    }
  }

  out
}


.pp_etas_sim_build_support <- function(
    n_support_obs,
    mult_support,
    window_working,
    bg_spatial_cov_funs,
    bg_spatial_cov_generator,
    bg_lp_type,
    bg_lp_coefs,
    bg_lp_fun,
    clip_eta_bg,
    seed
) {
  n_support_dummy <- as.integer(mult_support * n_support_obs)
  n_support_total <- as.integer(n_support_obs + n_support_dummy)

  if (n_support_obs < 1L || n_support_dummy < 1L) {
    stop(
      "Need `n_support_obs >= 1` and `mult_support * n_support_obs >= 1`.",
      call. = FALSE
    )
  }

  xmin <- window_working[["xmin"]]
  xmax <- window_working[["xmax"]]
  ymin <- window_working[["ymin"]]
  ymax <- window_working[["ymax"]]

  x_obs <- stats::runif(n_support_obs, xmin, xmax)
  y_obs <- stats::runif(n_support_obs, ymin, ymax)
  x_dummy <- stats::runif(n_support_dummy, xmin, xmax)
  y_dummy <- stats::runif(n_support_dummy, ymin, ymax)

  x <- c(x_obs, x_dummy)
  y <- c(y_obs, y_dummy)
  is_dummy <- c(rep(FALSE, n_support_obs), rep(TRUE, n_support_dummy))

  area <- (xmax - xmin) * (ymax - ymin)
  weight <- rep(area / n_support_total, n_support_total)

  features <- .pp_etas_sim_build_features(
    x = x,
    y = y,
    window_working = window_working
  )

  features <- .pp_etas_sim_append_spatial_covariates(
    df_features = features,
    bg_spatial_cov_funs = bg_spatial_cov_funs,
    bg_spatial_cov_generator = bg_spatial_cov_generator,
    window_working = window_working,
    seed = seed,
    stage = "support"
  )

  eta <- .pp_etas_sim_eval_background_lp(
    features = features,
    bg_lp_type = bg_lp_type,
    bg_lp_coefs = bg_lp_coefs,
    bg_lp_fun = bg_lp_fun,
    clip_eta_bg = clip_eta_bg
  )

  lambda_shape <- exp(eta)
  mass <- weight * lambda_shape
  mass_total <- sum(mass)

  if (!is.finite(mass_total) || mass_total <= 0) {
    stop("Invalid background support mass.", call. = FALSE)
  }

  density <- lambda_shape / sum(weight * lambda_shape)

  list(
    features = features,
    eta = eta,
    lambda_shape = lambda_shape,
    density = density,
    weight = weight,
    probability = mass / mass_total,
    is_dummy = is_dummy,
    n_support_obs = n_support_obs,
    n_support_dummy = n_support_dummy,
    n_support_total = n_support_total,
    area = area
  )
}


.pp_etas_sim_sample_background_xy <- function(
    n_need,
    window_working,
    eta_max,
    bg_spatial_cov_funs,
    bg_spatial_cov_generator,
    bg_lp_type,
    bg_lp_coefs,
    bg_lp_fun,
    clip_eta_bg,
    seed
) {
  if (n_need == 0L) {
    return(data.frame(x = numeric(0), y = numeric(0)))
  }

  xmin <- window_working[["xmin"]]
  xmax <- window_working[["xmax"]]
  ymin <- window_working[["ymin"]]
  ymax <- window_working[["ymax"]]

  out_x <- numeric(0)
  out_y <- numeric(0)
  batch <- max(1000L, ceiling(1.5 * n_need))
  guard <- 0L

  if (!is.finite(eta_max)) {
    eta_max <- 0
  }

  while (length(out_x) < n_need) {
    guard <- guard + 1L

    if (guard > 10000L) {
      stop("Background rejection sampler failed to terminate.", call. = FALSE)
    }

    x_candidate <- stats::runif(batch, xmin, xmax)
    y_candidate <- stats::runif(batch, ymin, ymax)

    features <- .pp_etas_sim_build_features(
      x = x_candidate,
      y = y_candidate,
      window_working = window_working
    )

    features <- .pp_etas_sim_append_spatial_covariates(
      df_features = features,
      bg_spatial_cov_funs = bg_spatial_cov_funs,
      bg_spatial_cov_generator = bg_spatial_cov_generator,
      window_working = window_working,
      seed = seed,
      stage = "events"
    )

    eta_candidate <- .pp_etas_sim_eval_background_lp(
      features = features,
      bg_lp_type = bg_lp_type,
      bg_lp_coefs = bg_lp_coefs,
      bg_lp_fun = bg_lp_fun,
      clip_eta_bg = clip_eta_bg
    )

    accepted <- stats::runif(batch) <= exp(eta_candidate - eta_max)

    if (any(accepted)) {
      out_x <- c(out_x, x_candidate[accepted])
      out_y <- c(out_y, y_candidate[accepted])
    }

    if (length(out_x) < n_need && guard %% 20L == 0L) {
      batch <- min(batch * 2L, 200000L)
    }
  }

  data.frame(
    x = out_x[seq_len(n_need)],
    y = out_y[seq_len(n_need)]
  )
}


# Output construction ---------------------------------------------------------

.pp_etas_sim_make_truth <- function(
    params,
    m0,
    b,
    bg_lp_type,
    bg_lp_coefs,
    trig_lp_type,
    trig_lp_coefs,
    window_working,
    tmin,
    tmax,
    longlat.to.km,
    sectoday
) {
  list(
    params = params,
    m0 = m0,
    b = b,
    bg = list(
      type = bg_lp_type,
      coefs = if (identical(bg_lp_type, "linear")) bg_lp_coefs else NULL,
      note = paste0(
        "Background intercept is not identifiable separately from mu when ",
        "the spatial shape is normalised."
      )
    ),
    trig = list(
      type = trig_lp_type,
      coefs = if (identical(trig_lp_type, "linear")) trig_lp_coefs else NULL
    ),
    window_km = window_working,
    window_working = window_working,
    t_window = c(tmin = tmin, tmax = tmax),
    longlat.to.km = longlat.to.km,
    sectoday = sectoday
  )
}


.pp_etas_sim_make_bg_info <- function(support, window_working, return_support) {
  support_df <- cbind(
    support$features,
    data.frame(
      w = support$weight,
      eta_bg_true = support$eta,
      lambda_shape = support$lambda_shape,
      dens_bg = support$density,
      is_dummy = support$is_dummy,
      check.names = FALSE
    )
  )

  list(
    window_km = window_working,
    window_working = window_working,
    support = if (isTRUE(return_support)) support_df else NULL,
    stations = NULL
  )
}


.pp_etas_sim_build_pattern_catalog <- function(
    cat_sim,
    long.range,
    lat.range,
    tmin,
    tmax,
    longlat.to.km,
    sectoday
) {
  .pp_etas_sim_require(
    c("pp_window", "pp_pattern", "etas_catalog"),
    context = "ETAS simulation output construction"
  )

  time_unit <- if (isTRUE(sectoday)) "secs" else "days"
  spatial_unit <- if (isTRUE(longlat.to.km)) "degrees" else "km"

  window <- pp_window(
    lower = c(long = long.range[1L], lat = lat.range[1L], time = tmin),
    upper = c(long = long.range[2L], lat = lat.range[2L], time = tmax),
    axis_roles = c(long = "space", lat = "space", time = "time"),
    axis_units = c(long = spatial_unit, lat = spatial_unit, time = time_unit),
    time_unit = time_unit
  )

  mark_names <- setdiff(names(cat_sim), c("long", "lat", "time"))

  pattern <- pp_pattern(
    data = cat_sim,
    coords = c("long", "lat", "time"),
    marks = mark_names,
    window = window,
    time_unit = time_unit
  )

  catalog <- etas_catalog(
    pattern = pattern,
    long_axis = "long",
    lat_axis = "lat",
    time_axis = "time",
    depth_mark = "z",
    magnitude_mark = "magn1",
    longlat.to.km = longlat.to.km,
    sectoday = sectoday
  )

  list(
    window = window,
    pattern = pattern,
    catalog = catalog
  )
}


.pp_etas_sim_validate_output <- function(cat_sim, tmin, tmax, long.range, lat.range, m0) {
  if (nrow(cat_sim) == 0L) {
    return(invisible(TRUE))
  }

  if (any(!is.finite(cat_sim$time)) ||
      any(cat_sim$time < tmin) || any(cat_sim$time > tmax)) {
    stop("Internal ETAS simulation error: event times lie outside the observation window.", call. = FALSE)
  }

  if (any(!is.finite(cat_sim$long)) ||
      any(cat_sim$long < long.range[1L]) || any(cat_sim$long > long.range[2L])) {
    stop("Internal ETAS simulation error: event long/x coordinates lie outside the observation window.", call. = FALSE)
  }

  if (any(!is.finite(cat_sim$lat)) ||
      any(cat_sim$lat < lat.range[1L]) || any(cat_sim$lat > lat.range[2L])) {
    stop("Internal ETAS simulation error: event lat/y coordinates lie outside the observation window.", call. = FALSE)
  }

  if (any(!is.finite(cat_sim$magn1)) || any(cat_sim$magn1 < m0)) {
    stop("Internal ETAS simulation error: invalid simulated magnitudes.", call. = FALSE)
  }

  if (any(cat_sim$event_id != seq_len(nrow(cat_sim)))) {
    stop("Internal ETAS simulation error: invalid event identifiers.", call. = FALSE)
  }

  if (any(cat_sim$father_id < 0L) || any(cat_sim$father_id > nrow(cat_sim))) {
    stop("Internal ETAS simulation error: invalid father identifiers.", call. = FALSE)
  }

  if (any(cat_sim$father_id > 0L & cat_sim$father_id >= cat_sim$event_id)) {
    stop(
      "Internal ETAS simulation error: a triggered event does not point to an earlier parent.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.pp_etas_sim_make_result <- function(
    call,
    cat_sim,
    cat_pois,
    n0,
    nson,
    exploded,
    truth,
    bg_info,
    output_objects,
    settings
) {
  genealogy <- cat_sim[, intersect(c("event_id", "father_id", "lgen"), names(cat_sim)), drop = FALSE]

  n_events <- nrow(cat_sim)
  max_generation <- if (n_events > 0L) max(cat_sim$lgen) else 0L

  out <- list(
    call = call,
    pattern = output_objects$pattern,
    catalog = output_objects$catalog,
    window = output_objects$window,

    # Legacy-compatible components
    cat = NULL,
    cat.pois = cat_pois,
    cat.sim = cat_sim,
    n0 = as.integer(n0),
    nson = as.integer(nson),
    exploded = isTRUE(exploded),
    truth = truth,
    bg_info = bg_info,

    # Package-native summaries
    genealogy = genealogy,
    n_events = as.integer(n_events),
    n_background = as.integer(n0),
    n_triggered = as.integer(nson),
    background_fraction = if (n_events > 0L) n0 / n_events else NA_real_,
    triggered_fraction = if (n_events > 0L) nson / n_events else NA_real_,
    max_generation = as.integer(max_generation),
    settings = settings
  )

  class(out) <- "pp_etas_simulation"
  out
}


# Main simulator --------------------------------------------------------------

#' Simulate a parametric ETAS earthquake catalog
#'
#' Simulates a branching epidemic-type aftershock sequence (ETAS) catalog on a
#' rectangular two-dimensional spatial domain and a finite temporal window.
#'
#' @details
#' The simulator uses the branching representation of an ETAS process. First,
#' background or immigrant events are generated from a Poisson process. The
#' expected number of immigrants is
#'
#' \deqn{\mu (t_{\max} - t_{\min}),}
#'
#' while their spatial density is proportional to
#' \eqn{\exp\{\eta_{\mathrm{bg}}(x,y)\}} and is numerically normalised over
#' the rectangular observation window. The background predictor may be
#' constant, linear in built-in or user-generated features, or supplied as a
#' custom function.
#'
#' Conditional on its magnitude and covariates, every event may produce direct
#' offspring. Magnitudes are generated above `m0` according to a
#' Gutenberg--Richter distribution with parameter `b`. The expected number of
#' direct offspring is controlled by `k0` and by the triggering linear
#' predictor, which by default depends on relative magnitude. Offspring times
#' follow a modified Omori law governed by `c` and `p`, whereas offspring
#' locations follow the isotropic spatial triggering kernel governed by `d`
#' and `q`. The branching construction is iterated generation by generation
#' until no new offspring are produced or the optional explosion guard is
#' reached.
#'
#' Complete genealogy is retained. `event_id` identifies each simulated event,
#' `father_id` is zero for an immigrant and otherwise identifies its direct
#' parent, and `lgen` records the branching generation. Thus `lgen = 0`
#' identifies background events, `lgen = 1` their direct offspring, and so on.
#'
#' When `longlat.to.km = TRUE`, the original longitude and latitude coordinates
#' are converted to the planar working coordinates used by the ETAS triggering
#' calculation. Both representations are retained in the returned pattern and
#' catalog. Similarly, `sectoday` controls the interpretation of raw simulated
#' times when the package-native ETAS catalog is constructed.
#'
#' The explosion guard is a numerical safety mechanism rather than part of the
#' stochastic model. If the cumulative number of events exceeds `max_events`,
#' branching is stopped and the returned object's `exploded` component is set
#' to `TRUE`.
#'
#' @param params Named numeric vector containing `mu`, `k0`, `c`, `p`,
#'   `gamma`, `d`, and `q`.
#' @param m0 Magnitude threshold.
#' @param b Gutenberg--Richter `b` value.
#' @param tmin Lower temporal bound.
#' @param t.lag Duration of the temporal observation window.
#' @param long.range,lat.range Spatial observation ranges. They are interpreted
#'   as longitude/latitude when `longlat.to.km = TRUE`, and otherwise as planar
#'   working coordinates.
#' @param longlat.to.km Logical. If TRUE, use the legacy conversion
#'   `6371.3 * pi / 180` between original coordinates and ETAS working
#'   coordinates.
#' @param sectoday Logical. If TRUE, raw simulated time is interpreted in
#'   seconds by the resulting `pp_etas_catalog` and converted to days for ETAS
#'   likelihood evaluation. The simulation itself uses the raw units supplied
#'   through `tmin`, `t.lag`, and `c`.
#' @param n_support_obs Number of observed-type support points used to represent
#'   the true background spatial shape.
#' @param mult_support Multiplier determining the number of dummy support
#'   points.
#' @param bg_lp_type Background linear-predictor type: `"constant"`,
#'   `"linear"`, or `"custom"`.
#' @param bg_lp_coefs Named coefficient vector for a linear background
#'   predictor. Available built-in features are `x`, `y`, `x_std`, `y_std`, and
#'   `r2_center`, together with any generated spatial covariates.
#' @param bg_lp_fun Custom function receiving a feature data frame and returning
#'   one numeric predictor value per row.
#' @param bg_spatial_cov_funs Optional named list of deterministic functions
#'   used to create additional spatial covariates.
#' @param bg_spatial_cov_generator Optional function with signature
#'   `function(df_features, window_km, seed = NULL)` returning a data frame of
#'   spatial covariates.
#' @param trig_lp_type Triggering predictor type: `"linear"` or `"custom"`.
#' @param trig_lp_coefs Named coefficient vector for a linear triggering
#'   predictor. The usual choice is `c(m_rel = alpha)`.
#' @param trig_lp_fun Custom function receiving event features and returning one
#'   numeric triggering predictor per row.
#' @param seed Optional integer random seed.
#' @param clip_eta_bg,clip_eta_trig Absolute clipping limits for background and
#'   triggering predictors.
#' @param max_events Maximum allowed number of simulated events when
#'   `explode_guard = TRUE`.
#' @param explode_guard Logical. Stop the branching loop and flag the simulation
#'   when `max_events` is exceeded.
#' @param return_bg_info Logical. Retain information about the true background
#'   spatial shape.
#' @param return_support Logical. Retain the full background support data frame
#'   inside `bg_info`.
#' @param plot_bg Logical. Plot the true simulated background spatial density
#'   after simulation, when background support information is retained.
#' @param plot_catalog Logical. Plot the simulated ETAS catalog after simulation.
#'   Background events are shown as blue open circles and triggered events as
#'   red open circles, with point sizes proportional to earthquake magnitude.
#'
#' @return An object of class `pp_etas_simulation` containing a `pp` pattern, a
#'   `pp_etas_catalog`, the complete simulated event table, immigrant events,
#'   genealogy, true parameters, background information, and simulation
#'   diagnostics.
#'
#' @references
#' Ogata, Y. (1988). Statistical models for earthquake occurrences and residual
#' analysis for point processes. *Journal of the American Statistical
#' Association*, 83, 9--27.
#'
#' Chiodi, M. and Adelfio, G. (2017). Mixed non-parametric and parametric
#' estimation techniques in R package etasFLP for earthquakes' description.
#' *Journal of Statistical Software*, 76(3), 1--29.
#'
#' @seealso [etas_catalog()] for constructing ETAS catalogs from observed point
#'   patterns and [fit_etas_pp()] for fitting a parametric-background ETAS model.
#'
#' @family ETAS modelling
#'
#' @examples
#' \donttest{
#' simulation <- simulate_etas_pp(
#'   params = c(
#'     mu = 0.5, k0 = 0.003, c = 0.05, p = 1.2,
#'     gamma = 0, d = 0.05, q = 2
#'   ),
#'   m0 = 2.5,
#'   b = 1,
#'   tmin = 0,
#'   t.lag = 10,
#'   long.range = c(0, 1),
#'   lat.range = c(0, 1),
#'   longlat.to.km = FALSE,
#'   n_support_obs = 50,
#'   mult_support = 1,
#'   bg_lp_type = "linear",
#'   bg_lp_coefs = c(x_std = 0.4, y_std = -0.3),
#'   trig_lp_coefs = c(m_rel = 0.3),
#'   seed = 19,
#'   max_events = 500,
#'   return_support = TRUE
#' )
#'
#' simulation
#' summary(simulation)
#' plot(simulation)
#' plot(simulation, what = "background")
#' }
#'
#' @export
simulate_etas_pp <- function(
    params = c(
      mu = 0.5,
      k0 = 0.012,
      c = 0.01,
      p = 1.12,
      gamma = 0,
      d = 1.6,
      q = 1.85
    ),
    m0 = 2.5,
    b = 1.05,
    tmin = 0,
    t.lag = 700,
    long.range = c(6, 20),
    lat.range = c(36, 45),
    longlat.to.km = TRUE,
    sectoday = FALSE,
    n_support_obs = 1500,
    mult_support = 4,
    bg_lp_type = c("constant", "linear", "custom"),
    bg_lp_coefs = NULL,
    bg_lp_fun = NULL,
    bg_spatial_cov_funs = NULL,
    bg_spatial_cov_generator = NULL,
    trig_lp_type = c("linear", "custom"),
    trig_lp_coefs = c(m_rel = 0.7),
    trig_lp_fun = NULL,
    seed = NULL,
    clip_eta_bg = 20,
    clip_eta_trig = 20,
    max_events = 50000,
    explode_guard = TRUE,
    return_bg_info = TRUE,
    return_support = TRUE,
    plot_bg = FALSE,
    plot_catalog = FALSE
) {

  .pp_etas_sim_require(
    c(
      "pp_window",
      "pp_pattern",
      "etas_catalog"
    ),
    context = "simulate_etas_pp()"
  )

  simulation_call <- match.call()

  bg_lp_type <- match.arg(bg_lp_type)
  trig_lp_type <- match.arg(trig_lp_type)

  params <- .pp_etas_sim_validate_params(params)

  m0 <- .pp_etas_sim_validate_numeric_scalar(
    m0,
    "m0"
  )

  b <- .pp_etas_sim_validate_numeric_scalar(
    b,
    "b",
    lower = 0,
    inclusive = FALSE
  )

  tmin <- .pp_etas_sim_validate_numeric_scalar(
    tmin,
    "tmin"
  )

  t.lag <- .pp_etas_sim_validate_numeric_scalar(
    t.lag,
    "t.lag",
    lower = 0,
    inclusive = FALSE
  )

  long.range <- .pp_etas_sim_validate_range(
    long.range,
    "long.range"
  )

  lat.range <- .pp_etas_sim_validate_range(
    lat.range,
    "lat.range"
  )

  longlat.to.km <- .pp_etas_sim_validate_logical(
    longlat.to.km,
    "longlat.to.km"
  )

  sectoday <- .pp_etas_sim_validate_logical(
    sectoday,
    "sectoday"
  )

  n_support_obs <- .pp_etas_sim_validate_positive_integer(
    n_support_obs,
    "n_support_obs",
    min_value = 1L
  )

  mult_support <- .pp_etas_sim_validate_numeric_scalar(
    mult_support,
    "mult_support",
    lower = 0,
    inclusive = FALSE
  )

  seed <- .pp_etas_sim_validate_seed(seed)

  clip_eta_bg <- .pp_etas_sim_validate_numeric_scalar(
    clip_eta_bg,
    "clip_eta_bg",
    lower = 0,
    inclusive = FALSE
  )

  clip_eta_trig <- .pp_etas_sim_validate_numeric_scalar(
    clip_eta_trig,
    "clip_eta_trig",
    lower = 0,
    inclusive = FALSE
  )

  max_events <- .pp_etas_sim_validate_positive_integer(
    max_events,
    "max_events",
    min_value = 1L
  )

  explode_guard <- .pp_etas_sim_validate_logical(
    explode_guard,
    "explode_guard"
  )

  return_bg_info <- .pp_etas_sim_validate_logical(
    return_bg_info,
    "return_bg_info"
  )

  return_support <- .pp_etas_sim_validate_logical(
    return_support,
    "return_support"
  )

  plot_bg <- .pp_etas_sim_validate_logical(
    plot_bg,
    "plot_bg"
  )

  plot_catalog <- .pp_etas_sim_validate_logical(
    plot_catalog,
    "plot_catalog"
  )

  if (identical(bg_lp_type, "linear")) {

    bg_lp_coefs <- .pp_etas_sim_validate_linear_coefs(
      bg_lp_coefs,
      "background linear LP"
    )
  }

  if (
    identical(bg_lp_type, "custom") &&
    !is.function(bg_lp_fun)
  ) {
    stop(
      "`bg_lp_fun` must be a function when `bg_lp_type = 'custom'`.",
      call. = FALSE
    )
  }

  if (identical(trig_lp_type, "linear")) {

    trig_lp_coefs <- .pp_etas_sim_validate_linear_coefs(
      trig_lp_coefs,
      "triggering linear LP"
    )
  }

  if (
    identical(trig_lp_type, "custom") &&
    !is.function(trig_lp_fun)
  ) {
    stop(
      "`trig_lp_fun` must be a function when `trig_lp_type = 'custom'`.",
      call. = FALSE
    )
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  tmax <- tmin + t.lag

  if (
    !is.finite(tmax) ||
    tmax <= tmin
  ) {
    stop(
      "Invalid time window: `tmax <= tmin`.",
      call. = FALSE
    )
  }

  unit <- if (isTRUE(longlat.to.km)) {
    6371.3 * pi / 180
  } else {
    1
  }

  xmin <- long.range[1L] * unit
  xmax <- long.range[2L] * unit
  ymin <- lat.range[1L] * unit
  ymax <- lat.range[2L] * unit

  if (
    xmin >= xmax ||
    ymin >= ymax
  ) {
    stop(
      "Invalid `long.range` or `lat.range` after coordinate conversion.",
      call. = FALSE
    )
  }

  window_working <- c(
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax
  )

  support <- .pp_etas_sim_build_support(
    n_support_obs = n_support_obs,
    mult_support = mult_support,
    window_working = window_working,
    bg_spatial_cov_funs = bg_spatial_cov_funs,
    bg_spatial_cov_generator = bg_spatial_cov_generator,
    bg_lp_type = bg_lp_type,
    bg_lp_coefs = bg_lp_coefs,
    bg_lp_fun = bg_lp_fun,
    clip_eta_bg = clip_eta_bg,
    seed = seed
  )

  truth <- .pp_etas_sim_make_truth(
    params = params,
    m0 = m0,
    b = b,
    bg_lp_type = bg_lp_type,
    bg_lp_coefs = bg_lp_coefs,
    trig_lp_type = trig_lp_type,
    trig_lp_coefs = trig_lp_coefs,
    window_working = window_working,
    tmin = tmin,
    tmax = tmax,
    longlat.to.km = longlat.to.km,
    sectoday = sectoday
  )

  bg_info <- if (isTRUE(return_bg_info)) {

    .pp_etas_sim_make_bg_info(
      support = support,
      window_working = window_working,
      return_support = return_support
    )

  } else {

    NULL
  }

  mu <- params[["mu"]]
  k0 <- params[["k0"]]
  c_ <- params[["c"]]
  p <- params[["p"]]
  gamma <- params[["gamma"]]
  d <- params[["d"]]
  q <- params[["q"]]

  expected_background <- mu * (tmax - tmin)

  n0 <- stats::rpois(
    1L,
    expected_background
  )

  base_feature_names <- c(
    "x",
    "y",
    "x_std",
    "y_std",
    "r2_center",
    "m_rel"
  )

  spatial_covariate_names <- setdiff(
    names(support$features),
    base_feature_names
  )

  # -----------------------------------------------------------------------
  # Empty-catalog case
  # -----------------------------------------------------------------------

  if (n0 == 0L) {

    cat_sim <- .pp_etas_sim_empty_catalog_data(
      spatial_covariate_names
    )

    output_objects <- .pp_etas_sim_build_pattern_catalog(
      cat_sim = cat_sim,
      long.range = long.range,
      lat.range = lat.range,
      tmin = tmin,
      tmax = tmax,
      longlat.to.km = longlat.to.km,
      sectoday = sectoday
    )

    out <- .pp_etas_sim_make_result(
      call = simulation_call,
      cat_sim = cat_sim,
      cat_pois = cat_sim,
      n0 = 0L,
      nson = 0L,
      exploded = FALSE,
      truth = truth,
      bg_info = bg_info,
      output_objects = output_objects,
      settings = list(
        seed = seed,
        expected_background = expected_background,
        bg_lp_type = bg_lp_type,
        trig_lp_type = trig_lp_type,
        n_support_obs = n_support_obs,
        mult_support = mult_support,
        clip_eta_bg = clip_eta_bg,
        clip_eta_trig = clip_eta_trig,
        max_events = max_events,
        explode_guard = explode_guard,
        longlat.to.km = longlat.to.km,
        sectoday = sectoday
      )
    )

    if (
      isTRUE(plot_bg) &&
      !is.null(out$bg_info) &&
      !is.null(out$bg_info$support)
    ) {
      graphics::plot(
        out,
        what = "background"
      )
    }

    if (isTRUE(plot_catalog)) {
      graphics::plot(
        out,
        what = "catalog"
      )
    }

    return(out)
  }

  # -----------------------------------------------------------------------
  # Simulate immigrant events
  # -----------------------------------------------------------------------

  xy0 <- .pp_etas_sim_sample_background_xy(
    n_need = n0,
    window_working = window_working,
    eta_max = max(
      support$eta,
      na.rm = TRUE
    ),
    bg_spatial_cov_funs = bg_spatial_cov_funs,
    bg_spatial_cov_generator = bg_spatial_cov_generator,
    bg_lp_type = bg_lp_type,
    bg_lp_coefs = bg_lp_coefs,
    bg_lp_fun = bg_lp_fun,
    clip_eta_bg = clip_eta_bg,
    seed = seed
  )

  x0 <- xy0$x
  y0 <- xy0$y

  t0 <- stats::runif(
    n0,
    tmin,
    tmax
  )

  beta_gr <- log(10) * b

  magnitude0 <- m0 + stats::rexp(
    n0,
    rate = beta_gr
  )

  mrel0 <- magnitude0 - m0

  features0 <- .pp_etas_sim_build_features(
    x = x0,
    y = y0,
    window_working = window_working,
    m_rel = mrel0
  )

  features0 <- .pp_etas_sim_append_spatial_covariates(
    df_features = features0,
    bg_spatial_cov_funs = bg_spatial_cov_funs,
    bg_spatial_cov_generator = bg_spatial_cov_generator,
    window_working = window_working,
    seed = seed,
    stage = "events"
  )

  eta_bg0 <- .pp_etas_sim_eval_background_lp(
    features = features0,
    bg_lp_type = bg_lp_type,
    bg_lp_coefs = bg_lp_coefs,
    bg_lp_fun = bg_lp_fun,
    clip_eta_bg = clip_eta_bg
  )

  eta_trig0 <- .pp_etas_sim_eval_trigger_lp(
    features = features0,
    trig_lp_type = trig_lp_type,
    trig_lp_coefs = trig_lp_coefs,
    trig_lp_fun = trig_lp_fun,
    clip_eta_trig = clip_eta_trig
  )

  if (isTRUE(longlat.to.km)) {

    long0 <- x0 / unit
    lat0 <- y0 / unit

  } else {

    long0 <- x0
    lat0 <- y0
  }

  cat_pois <- data.frame(
    event_id = seq_len(n0),
    father_id = 0L,
    lgen = 0L,
    time = t0,
    lat = lat0,
    long = long0,
    z = 0,
    magn1 = magnitude0,
    x_km = x0,
    y_km = y0,
    m_rel = mrel0,
    x_true = x0,
    y_true = y0,
    x_std_true = features0$x_std,
    y_std_true = features0$y_std,
    r2_center_true = features0$r2_center,
    eta_bg_true = eta_bg0,
    eta_trig_true = eta_trig0,
    check.names = FALSE
  )

  spatial_covariate_names <- setdiff(
    names(features0),
    c(
      "x",
      "y",
      "x_std",
      "y_std",
      "r2_center",
      "m_rel"
    )
  )

  if (length(spatial_covariate_names) > 0L) {

    for (name_i in spatial_covariate_names) {
      cat_pois[[name_i]] <- features0[[name_i]]
    }
  }

  cat_pois <- cat_pois[
    order(cat_pois$time),
    ,
    drop = FALSE
  ]

  rownames(cat_pois) <- NULL

  cat_pois$event_id <- seq_len(
    nrow(cat_pois)
  )

  cat_new <- cat_pois

  # Infinite-support temporal and spatial kernel integrals, exactly as in the
  # validated etas.par.sim_v4() implementation. Offspring outside the finite
  # observation window are subsequently rejected.

  ak <- k0 * c_^(1 - p) / (p - 1)

  sk <- pi * d^(1 - q) / (q - 1)

  i <- 0L
  exploded <- FALSE

  # -----------------------------------------------------------------------
  # Simulate ETAS branching generations
  # -----------------------------------------------------------------------

  while (i < nrow(cat_new)) {

    i <- i + 1L

    parent_features <- data.frame(
      x = cat_new$x_km[i],
      y = cat_new$y_km[i],
      x_std = cat_new$x_std_true[i],
      y_std = cat_new$y_std_true[i],
      r2_center = cat_new$r2_center_true[i],
      m_rel = cat_new$m_rel[i],
      check.names = FALSE
    )

    if (length(spatial_covariate_names) > 0L) {

      for (name_i in spatial_covariate_names) {
        parent_features[[name_i]] <- cat_new[[name_i]][i]
      }
    }

    eta_trig_i <- .pp_etas_sim_eval_trigger_lp(
      features = parent_features,
      trig_lp_type = trig_lp_type,
      trig_lp_coefs = trig_lp_coefs,
      trig_lp_fun = trig_lp_fun,
      clip_eta_trig = clip_eta_trig
    )

    cat_new$eta_trig_true[i] <- eta_trig_i

    n_expected_i <-
      ak *
      sk *
      exp(
        gamma * cat_new$m_rel[i] +
          eta_trig_i
      )

    if (
      !is.finite(n_expected_i) ||
      n_expected_i < 0
    ) {
      n_expected_i <- 0
    }

    n_child <- stats::rpois(
      1L,
      lambda = n_expected_i
    )

    if (n_child > 0L) {

      t_child <-
        c_ *
        stats::runif(n_child)^(
          -1 / (p - 1)
        ) -
        c_ +
        cat_new$time[i]

      xy_child <- .pp_etas_sim_sample_xy(
        n = n_child,
        x0 = cat_new$x_km[i],
        y0 = cat_new$y_km[i],
        d = d,
        q = q
      )

      inside <-
        t_child > tmin &
        t_child < tmax &
        xy_child[, 1L] > xmin &
        xy_child[, 1L] < xmax &
        xy_child[, 2L] > ymin &
        xy_child[, 2L] < ymax

      if (any(inside)) {

        n_inside <- sum(inside)

        x1 <- xy_child[inside, 1L]
        y1 <- xy_child[inside, 2L]
        t1 <- t_child[inside]

        magnitude1 <-
          m0 +
          stats::rexp(
            n_inside,
            rate = beta_gr
          )

        mrel1 <- magnitude1 - m0

        features1 <- .pp_etas_sim_build_features(
          x = x1,
          y = y1,
          window_working = window_working,
          m_rel = mrel1
        )

        features1 <- .pp_etas_sim_append_spatial_covariates(
          df_features = features1,
          bg_spatial_cov_funs = bg_spatial_cov_funs,
          bg_spatial_cov_generator = bg_spatial_cov_generator,
          window_working = window_working,
          seed = seed,
          stage = "events"
        )

        eta_bg1 <- .pp_etas_sim_eval_background_lp(
          features = features1,
          bg_lp_type = bg_lp_type,
          bg_lp_coefs = bg_lp_coefs,
          bg_lp_fun = bg_lp_fun,
          clip_eta_bg = clip_eta_bg
        )

        eta_trig1 <- .pp_etas_sim_eval_trigger_lp(
          features = features1,
          trig_lp_type = trig_lp_type,
          trig_lp_coefs = trig_lp_coefs,
          trig_lp_fun = trig_lp_fun,
          clip_eta_trig = clip_eta_trig
        )

        if (isTRUE(longlat.to.km)) {

          long1 <- x1 / unit
          lat1 <- y1 / unit

        } else {

          long1 <- x1
          lat1 <- y1
        }

        child_df <- data.frame(
          event_id = NA_integer_,
          father_id = as.integer(
            cat_new$event_id[i]
          ),
          lgen = as.integer(
            cat_new$lgen[i] + 1L
          ),
          time = t1,
          lat = lat1,
          long = long1,
          z = 0,
          magn1 = magnitude1,
          x_km = x1,
          y_km = y1,
          m_rel = mrel1,
          x_true = x1,
          y_true = y1,
          x_std_true = features1$x_std,
          y_std_true = features1$y_std,
          r2_center_true = features1$r2_center,
          eta_bg_true = eta_bg1,
          eta_trig_true = eta_trig1,
          check.names = FALSE
        )

        if (length(spatial_covariate_names) > 0L) {

          for (name_i in spatial_covariate_names) {
            child_df[[name_i]] <- features1[[name_i]]
          }
        }

        cat_new <- .pp_etas_sim_rbind_fill(
          cat_new,
          child_df
        )

        cat_new <- cat_new[
          order(cat_new$time),
          ,
          drop = FALSE
        ]

        rownames(cat_new) <- NULL

        old_ids <- cat_new$event_id

        cat_new$event_id <- seq_len(
          nrow(cat_new)
        )

        valid_old <- !is.na(old_ids)

        max_old <- max(
          old_ids[valid_old],
          na.rm = TRUE
        )

        id_map <- integer(max_old)

        id_map[
          old_ids[valid_old]
        ] <- which(valid_old)

        father_rows <- which(
          cat_new$father_id > 0L
        )

        if (length(father_rows) > 0L) {

          cat_new$father_id[
            father_rows
          ] <- id_map[
            cat_new$father_id[father_rows]
          ]
        }

        # Check the explosion guard only after the new offspring have been
        # sorted and given valid event and father identifiers, so that a
        # stopped simulation still returns an internally consistent catalog.
        if (
          isTRUE(explode_guard) &&
          nrow(cat_new) > max_events
        ) {

          warning(
            paste0(
              "Maximum number of events exceeded (`max_events = ",
              max_events,
              "`). Stopping branching simulation."
            ),
            call. = FALSE
          )

          exploded <- TRUE
          break
        }
      }
    }

    if (exploded) {
      break
    }
  }

  nson <- nrow(cat_new) - n0

  cat_pois_final <- cat_new[
    cat_new$father_id == 0L,
    ,
    drop = FALSE
  ]

  cat_pois_final <- cat_pois_final[
    order(cat_pois_final$time),
    ,
    drop = FALSE
  ]

  rownames(cat_pois_final) <- NULL

  # -----------------------------------------------------------------------
  # Validate and construct package-native output objects
  # -----------------------------------------------------------------------

  .pp_etas_sim_validate_output(
    cat_sim = cat_new,
    tmin = tmin,
    tmax = tmax,
    long.range = long.range,
    lat.range = lat.range,
    m0 = m0
  )

  output_objects <- .pp_etas_sim_build_pattern_catalog(
    cat_sim = cat_new,
    long.range = long.range,
    lat.range = lat.range,
    tmin = tmin,
    tmax = tmax,
    longlat.to.km = longlat.to.km,
    sectoday = sectoday
  )

  out <- .pp_etas_sim_make_result(
    call = simulation_call,
    cat_sim = cat_new,
    cat_pois = cat_pois_final,
    n0 = n0,
    nson = nson,
    exploded = exploded,
    truth = truth,
    bg_info = bg_info,
    output_objects = output_objects,
    settings = list(
      seed = seed,
      expected_background = expected_background,
      bg_lp_type = bg_lp_type,
      bg_lp_coefs = bg_lp_coefs,
      bg_lp_fun = bg_lp_fun,
      bg_spatial_cov_funs = bg_spatial_cov_funs,
      bg_spatial_cov_generator = bg_spatial_cov_generator,
      trig_lp_type = trig_lp_type,
      trig_lp_coefs = trig_lp_coefs,
      trig_lp_fun = trig_lp_fun,
      n_support_obs = n_support_obs,
      mult_support = mult_support,
      clip_eta_bg = clip_eta_bg,
      clip_eta_trig = clip_eta_trig,
      max_events = max_events,
      explode_guard = explode_guard,
      longlat.to.km = longlat.to.km,
      sectoday = sectoday
    )
  )

  # -----------------------------------------------------------------------
  # Optional immediate plots
  # -----------------------------------------------------------------------

  if (
    isTRUE(plot_bg) &&
    !is.null(out$bg_info) &&
    !is.null(out$bg_info$support)
  ) {
    graphics::plot(
      out,
      what = "background"
    )
  }

  if (isTRUE(plot_catalog)) {
    graphics::plot(
      out,
      what = "catalog"
    )
  }

  out
}


# Printed summaries and plots -------------------------------------------------

#' Methods for simulated ETAS catalogs
#'
#' Print, summarise, and plot objects returned by [simulate_etas_pp()].
#'
#' @details
#' The print method reports the total number of events, the background and
#' triggered counts, the maximum branching generation, the explosion-guard
#' status, and the true ETAS parameters. The summary method additionally
#' returns summaries of magnitudes and direct-offspring counts.
#'
#' The plot method can display either the simulated catalog or the true
#' background support. For `what = "catalog"`, background and triggered events
#' are distinguished by colour and point size represents magnitude. For
#' `what = "background"`, the retained support points are coloured according
#' to the normalised true background spatial density. The latter option requires
#' simulation with `return_bg_info = TRUE` and `return_support = TRUE`.
#'
#' @param x A `pp_etas_simulation` object for `print()` and `plot()`, or a
#'   `summary_pp_etas_simulation` object for its print method.
#' @param object A `pp_etas_simulation` object to summarise.
#' @param digits Number of significant digits used for printed numerical values.
#' @param what One of `"catalog"` or `"background"`. The first displays the
#'   simulated events; the second displays the retained true background support.
#' @param point_cex_range Numeric vector of length two giving the minimum and
#'   maximum point sizes used to represent magnitudes.
#' @param point_pch Plotting symbol. The default `1` gives open circles.
#' @param class_colours Named or unnamed vector of two colours for background
#'   and triggered events.
#' @param show_window Logical; draw the observation-window boundary.
#' @param show_class_legend Logical; display the background/triggered legend.
#' @param show_magnitude_legend Logical; display representative point sizes for
#'   earthquake magnitude.
#' @param cex_axis Size of axis tick labels.
#' @param cex_lab Size of axis labels.
#' @param cex_main Size of the main title.
#' @param legend_cex Size of legend text.
#' @param main Optional plot title.
#' @param background_palette HCL palette used by `what = "background"`.
#' @param ... Additional arguments. For `plot()`, these are passed to the
#'   underlying base-R plotting function; the print and summary methods ignore
#'   additional arguments.
#'
#' @return The print and plot methods invisibly return `x`.
#' `summary.pp_etas_simulation()` returns an object of class
#' `"summary_pp_etas_simulation"` containing event counts, fractions,
#' generation information, true parameters, magnitude and offspring summaries,
#' and retained simulation settings.
#'
#' @examples
#' \donttest{
#' simulation <- simulate_etas_pp(
#'   params = c(
#'     mu = 0.4, k0 = 0.002, c = 0.05, p = 1.2,
#'     gamma = 0, d = 0.05, q = 2
#'   ),
#'   t.lag = 8,
#'   long.range = c(0, 1),
#'   lat.range = c(0, 1),
#'   longlat.to.km = FALSE,
#'   n_support_obs = 40,
#'   mult_support = 1,
#'   seed = 27,
#'   max_events = 300,
#'   return_support = TRUE
#' )
#' print(simulation)
#' simulation_summary <- summary(simulation)
#' print(simulation_summary)
#' plot(simulation, what = "catalog")
#' plot(simulation, what = "background")
#' }
#'
#' @seealso [simulate_etas_pp()], [etas_catalog()], and [fit_etas_pp()].
#'
#' @family ETAS modelling
#' @rdname pp_etas_simulation-methods
#' @method print pp_etas_simulation
#' @export
print.pp_etas_simulation <- function(x, digits = 5, ...) {
  if (!inherits(x, "pp_etas_simulation")) {
    stop("`x` must inherit from class 'pp_etas_simulation'.", call. = FALSE)
  }

  cat("<pp_etas_simulation>\n")
  cat("Simulated parametric ETAS catalog\n\n")
  cat("Events: ", x$n_events, "\n", sep = "")
  cat("Background events: ", x$n_background, "\n", sep = "")
  cat("Triggered events: ", x$n_triggered, "\n", sep = "")

  if (is.finite(x$background_fraction)) {
    cat(
      "Background fraction: ",
      formatC(x$background_fraction, digits = digits, format = "fg"),
      "\n",
      sep = ""
    )
  }

  cat("Maximum generation: ", x$max_generation, "\n", sep = "")
  cat("Explosion guard triggered: ", x$exploded, "\n\n", sep = "")
  cat("True ETAS parameters:\n")
  print(round(x$truth$params, digits = digits))

  invisible(x)
}


#' @rdname pp_etas_simulation-methods
#' @method summary pp_etas_simulation
#' @export
summary.pp_etas_simulation <- function(object, ...) {
  if (!inherits(object, "pp_etas_simulation")) {
    stop("`object` must inherit from class 'pp_etas_simulation'.", call. = FALSE)
  }

  offspring_counts <- if (object$n_events > 0L) {
    tabulate(
      object$cat.sim$father_id[object$cat.sim$father_id > 0L],
      nbins = object$n_events
    )
  } else {
    integer(0)
  }

  out <- list(
    n_events = object$n_events,
    n_background = object$n_background,
    n_triggered = object$n_triggered,
    background_fraction = object$background_fraction,
    triggered_fraction = object$triggered_fraction,
    max_generation = object$max_generation,
    exploded = object$exploded,
    true_params = object$truth$params,
    background_type = object$truth$bg$type,
    triggering_type = object$truth$trig$type,
    magnitude_summary = if (object$n_events > 0L) summary(object$cat.sim$magn1) else NULL,
    offspring_summary = if (length(offspring_counts) > 0L) summary(offspring_counts) else NULL,
    settings = object$settings
  )

  class(out) <- "summary_pp_etas_simulation"
  out
}


#' @rdname pp_etas_simulation-methods
#' @method print summary_pp_etas_simulation
#' @export
print.summary_pp_etas_simulation <- function(x, digits = 5, ...) {
  cat("Parametric ETAS simulation summary\n\n")
  cat("Events: ", x$n_events, "\n", sep = "")
  cat("Background events: ", x$n_background, "\n", sep = "")
  cat("Triggered events: ", x$n_triggered, "\n", sep = "")

  if (is.finite(x$background_fraction)) {
    cat(
      "Background fraction: ",
      formatC(x$background_fraction, digits = digits, format = "fg"),
      "\n",
      sep = ""
    )
  }

  cat("Maximum generation: ", x$max_generation, "\n", sep = "")
  cat("Explosion guard triggered: ", x$exploded, "\n", sep = "")
  cat("Background type: ", x$background_type, "\n", sep = "")
  cat("Triggering LP type: ", x$triggering_type, "\n\n", sep = "")

  cat("True ETAS parameters:\n")
  print(round(x$true_params, digits = digits))

  if (!is.null(x$magnitude_summary)) {
    cat("\nMagnitude summary:\n")
    print(x$magnitude_summary, digits = digits)
  }

  invisible(x)
}


# Plot helpers ----------------------------------------------------------------

.pp_etas_sim_colour_values <- function(values, palette = "YlOrRd", n = 100L) {
  values <- as.numeric(values)
  colours <- grDevices::hcl.colors(n, palette = palette)

  if (length(values) == 0L) {
    return(list(colours = character(0), palette = colours, limits = c(0, 1)))
  }

  limits <- range(values, finite = TRUE)

  if (!all(is.finite(limits))) {
    limits <- c(0, 1)
  }

  if (diff(limits) <= sqrt(.Machine$double.eps)) {
    index <- rep.int(ceiling(n / 2), length(values))
    half_width <- max(abs(limits[1L]) * 0.01, 0.5)
    limits <- limits[1L] + c(-half_width, half_width)
  } else {
    scaled <- (values - limits[1L]) / diff(limits)
    index <- pmin(n, pmax(1L, floor(scaled * (n - 1L)) + 1L))
  }

  list(
    colours = colours[index],
    palette = colours,
    limits = limits
  )
}


.pp_etas_sim_draw_colorbar <- function(palette, limits, label, digits = 4) {
  n <- length(palette)
  old_mar <- graphics::par("mar")
  on.exit(graphics::par(mar = old_mar), add = TRUE)

  graphics::par(mar = c(5, 0.5, 4, 3.5))

  graphics::image(
    x = c(0, 1),
    y = seq(limits[1L], limits[2L], length.out = n),
    z = rbind(seq_len(n), seq_len(n)),
    col = palette,
    axes = FALSE,
    xlab = "",
    ylab = ""
  )

  ticks <- pretty(limits, n = 5)
  ticks <- ticks[ticks >= limits[1L] & ticks <= limits[2L]]

  graphics::axis(
    side = 4L,
    at = ticks,
    labels = format(signif(ticks, digits), trim = TRUE),
    las = 1
  )

  graphics::mtext(label, side = 4L, line = 2.4)
  invisible(NULL)
}


.pp_etas_sim_plot_spatial_base <- function(
    x,
    colours,
    main,
    point_size,
    point_pch,
    show_window,
    xlab,
    ylab,
    ...
) {
  cat_sim <- x$cat.sim

  xlim <- range(x$window$bounds[, "long"])
  ylim <- range(x$window$bounds[, "lat"])

  graphics::plot(
    NA_real_,
    NA_real_,
    type = "n",
    xlim = xlim,
    ylim = ylim,
    xlab = xlab,
    ylab = ylab,
    main = main,
    xaxs = "i",
    yaxs = "i",
    ...
  )

  if (isTRUE(show_window)) {
    graphics::rect(
      xleft = xlim[1L],
      ybottom = ylim[1L],
      xright = xlim[2L],
      ytop = ylim[2L],
      border = "grey35",
      lwd = 1.2
    )
  }

  if (nrow(cat_sim) > 0L) {
    graphics::points(
      cat_sim$long,
      cat_sim$lat,
      col = colours,
      pch = point_pch,
      cex = point_size
    )
  }

  invisible(NULL)
}



#' @rdname pp_etas_simulation-methods
#' @method plot pp_etas_simulation
#' @export
plot.pp_etas_simulation <- function(
    x,
    what = c("catalog", "background"),
    point_cex_range = c(0.9, 2.8),
    point_pch = 1,
    class_colours = c(
      Background = "blue",
      Triggered = "red"
    ),
    show_window = TRUE,
    show_class_legend = TRUE,
    show_magnitude_legend = TRUE,
    cex_axis = 1.3,
    cex_lab = 1.5,
    cex_main = 1.5,
    legend_cex = 1.15,
    main = NULL,
    background_palette = "YlOrRd",
    ...
) {

  # -----------------------------------------------------------------------
  # Validation
  # -----------------------------------------------------------------------

  if (!inherits(x, "pp_etas_simulation")) {
    stop(
      "`x` must inherit from class 'pp_etas_simulation'.",
      call. = FALSE
    )
  }

  what <- match.arg(what)

  show_window <- .pp_etas_sim_validate_logical(
    show_window,
    "show_window"
  )

  show_class_legend <- .pp_etas_sim_validate_logical(
    show_class_legend,
    "show_class_legend"
  )

  show_magnitude_legend <- .pp_etas_sim_validate_logical(
    show_magnitude_legend,
    "show_magnitude_legend"
  )

  if (
    !is.numeric(point_cex_range) ||
    length(point_cex_range) != 2L ||
    any(!is.finite(point_cex_range)) ||
    any(point_cex_range <= 0) ||
    point_cex_range[1L] >= point_cex_range[2L]
  ) {
    stop(
      paste0(
        "`point_cex_range` must contain two positive finite values, ",
        "with the first smaller than the second."
      ),
      call. = FALSE
    )
  }

  if (
    !is.character(class_colours) ||
    length(class_colours) != 2L
  ) {
    stop(
      "`class_colours` must contain exactly two colours.",
      call. = FALSE
    )
  }

  if (is.null(names(class_colours))) {
    names(class_colours) <- c(
      "Background",
      "Triggered"
    )
  }

  # -----------------------------------------------------------------------
  # Optional plot of the true background density
  # -----------------------------------------------------------------------

  if (identical(what, "background")) {

    if (
      is.null(x$bg_info) ||
      is.null(x$bg_info$support)
    ) {
      stop(
        paste0(
          "The true background support is unavailable. ",
          "Simulate with `return_bg_info = TRUE` and ",
          "`return_support = TRUE`."
        ),
        call. = FALSE
      )
    }

    support <- x$bg_info$support

    background_colours <- .pp_etas_sim_colour_values(
      support$dens_bg,
      palette = background_palette,
      n = 100L
    )

    old_par <- graphics::par(no.readonly = TRUE)

    on.exit(
      graphics::par(old_par),
      add = TRUE
    )

    graphics::par(
      mar = c(5.5, 6, 4.5, 2)
    )

    graphics::plot(
      support$x,
      support$y,
      pch = 16,
      cex = 0.7,
      col = background_colours$colours,
      xlab = if (isTRUE(x$truth$longlat.to.km)) {
        "x [km]"
      } else {
        "x"
      },
      ylab = if (isTRUE(x$truth$longlat.to.km)) {
        "y [km]"
      } else {
        "y"
      },
      main = if (is.null(main)) {
        "True background spatial density"
      } else {
        main
      },
      cex.axis = cex_axis,
      cex.lab = cex_lab,
      cex.main = cex_main,
      ...
    )

    return(invisible(x))
  }

  # -----------------------------------------------------------------------
  # Main ETAS catalog plot
  # -----------------------------------------------------------------------

  cat_sim <- x$cat.sim

  if (nrow(cat_sim) == 0L) {
    stop(
      "The simulated ETAS catalog contains no events.",
      call. = FALSE
    )
  }

  required_columns <- c(
    "long",
    "lat",
    "magn1",
    "father_id"
  )

  missing_columns <- setdiff(
    required_columns,
    names(cat_sim)
  )

  if (length(missing_columns) > 0L) {
    stop(
      sprintf(
        "The simulated catalog is missing required column(s): %s.",
        paste(
          missing_columns,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }

  # -----------------------------------------------------------------------
  # Event classification
  # -----------------------------------------------------------------------

  is_background <- cat_sim$father_id == 0L

  event_colours <- ifelse(
    is_background,
    class_colours[["Background"]],
    class_colours[["Triggered"]]
  )

  # -----------------------------------------------------------------------
  # Point size proportional to magnitude
  # -----------------------------------------------------------------------

  magnitudes <- cat_sim$magn1

  magnitude_range <- range(
    magnitudes,
    finite = TRUE
  )

  if (
    !all(is.finite(magnitude_range))
  ) {
    stop(
      "The simulated catalog contains invalid magnitudes.",
      call. = FALSE
    )
  }

  if (
    diff(magnitude_range) <=
    sqrt(.Machine$double.eps)
  ) {

    point_cex <- rep(
      mean(point_cex_range),
      length(magnitudes)
    )

  } else {

    point_cex <-
      point_cex_range[1L] +
      (
        magnitudes -
          magnitude_range[1L]
      ) /
      diff(magnitude_range) *
      diff(point_cex_range)
  }

  # -----------------------------------------------------------------------
  # Spatial limits
  # -----------------------------------------------------------------------

  if (
    !is.null(x$window) &&
    !is.null(x$window$bounds) &&
    all(c("long", "lat") %in%
        colnames(x$window$bounds))
  ) {

    xlim <- range(
      x$window$bounds[, "long"]
    )

    ylim <- range(
      x$window$bounds[, "lat"]
    )

  } else {

    xlim <- range(
      cat_sim$long,
      finite = TRUE
    )

    ylim <- range(
      cat_sim$lat,
      finite = TRUE
    )
  }

  # -----------------------------------------------------------------------
  # Axis labels
  # -----------------------------------------------------------------------

  xlab <- if (
    isTRUE(x$truth$longlat.to.km)
  ) {
    "Longitude"
  } else {
    "x"
  }

  ylab <- if (
    isTRUE(x$truth$longlat.to.km)
  ) {
    "Latitude"
  } else {
    "y"
  }

  plot_title <- if (is.null(main)) {
    "Simulated ETAS catalog"
  } else {
    main
  }

  # -----------------------------------------------------------------------
  # Plot
  # -----------------------------------------------------------------------

  old_par <- graphics::par(
    no.readonly = TRUE
  )

  on.exit(
    graphics::par(old_par),
    add = TRUE
  )

  graphics::par(
    mar = c(5.8, 6.2, 4.5, 2)
  )

  graphics::plot(
    NA_real_,
    NA_real_,
    type = "n",
    xlim = xlim,
    ylim = ylim,
    xlab = xlab,
    ylab = ylab,
    main = plot_title,
    xaxs = "i",
    yaxs = "i",
    cex.axis = cex_axis,
    cex.lab = cex_lab,
    cex.main = cex_main,
    ...
  )

  # Observation window
  if (isTRUE(show_window)) {

    graphics::rect(
      xleft = xlim[1L],
      ybottom = ylim[1L],
      xright = xlim[2L],
      ytop = ylim[2L],
      border = "grey35",
      lwd = 1.8
    )
  }

  # Plot triggered events first, then background events.
  # This makes the two classes easier to distinguish when circles overlap.

  triggered_index <- which(
    !is_background
  )

  background_index <- which(
    is_background
  )

  if (length(triggered_index) > 0L) {

    graphics::points(
      cat_sim$long[triggered_index],
      cat_sim$lat[triggered_index],
      pch = point_pch,
      cex = point_cex[triggered_index],
      col = event_colours[triggered_index],
      lwd = 1.5
    )
  }

  if (length(background_index) > 0L) {

    graphics::points(
      cat_sim$long[background_index],
      cat_sim$lat[background_index],
      pch = point_pch,
      cex = point_cex[background_index],
      col = event_colours[background_index],
      lwd = 1.5
    )
  }

  # -----------------------------------------------------------------------
  # Classification legend
  # -----------------------------------------------------------------------

  if (isTRUE(show_class_legend)) {

    graphics::legend(
      "topright",
      legend = c(
        "Background",
        "Triggered"
      ),
      col = c(
        class_colours[["Background"]],
        class_colours[["Triggered"]]
      ),
      pch = point_pch,
      pt.cex = 1.5,
      pt.lwd = 1.7,
      cex = legend_cex,
      bty = "n",
      title = "Event type"
    )
  }

  # -----------------------------------------------------------------------
  # Magnitude-size legend
  # -----------------------------------------------------------------------

  if (isTRUE(show_magnitude_legend)) {

    representative_magnitudes <- unique(
      as.numeric(
        stats::quantile(
          magnitudes,
          probs = c(
            0.1,
            0.5,
            0.9
          ),
          names = FALSE,
          na.rm = TRUE
        )
      )
    )

    if (
      diff(magnitude_range) <=
      sqrt(.Machine$double.eps)
    ) {

      representative_cex <- rep(
        mean(point_cex_range),
        length(representative_magnitudes)
      )

    } else {

      representative_cex <-
        point_cex_range[1L] +
        (
          representative_magnitudes -
            magnitude_range[1L]
        ) /
        diff(magnitude_range) *
        diff(point_cex_range)
    }

    graphics::legend(
      "bottomleft",
      legend = format(
        round(
          representative_magnitudes,
          2
        ),
        nsmall = 2
      ),
      col = "grey20",
      pch = point_pch,
      pt.cex = representative_cex,
      pt.lwd = 1.5,
      cex = legend_cex,
      bty = "n",
      title = "Magnitude"
    )
  }

  invisible(x)
}
