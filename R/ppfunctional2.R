# Functional covariates for Poisson point-process models ---------------------
#
# This module converts densely observed functional covariates into scalar
# functional principal-component scores. Each retained score is represented as
# a standard pp_covariate field and can therefore be supplied directly to
# fit_pp(), predict_pp_fit(), simulation, refitting, and diagnostic functions.
#
# Required package module:
#   ppdata.R -> pp_covariate()
#
# The module deliberately does not define a separate fitting wrapper. Poisson
# point-process models with functional covariates are fitted by fit_pp().

.pp_functional_version <- "0.1.0"


# Validation and numerical helpers ------------------------------------------

.pp_functional_require <- function(functions, context = "Functional covariates") {
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
          "Source ppdata.R before ppfunctional.R."
        ),
        context,
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.pp_functional_validate_name <- function(x, arg_name = "name") {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(trimws(x))) {
    stop(sprintf("`%s` must be one non-empty character string.", arg_name), call. = FALSE)
  }

  trimws(x)
}


.pp_functional_validate_probability <- function(x, arg_name = "pve") {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x > 1) {
    stop(sprintf("`%s` must be one finite number in (0, 1].", arg_name), call. = FALSE)
  }

  as.numeric(x)
}


.pp_functional_validate_positive_integer <- function(
    x,
    arg_name,
    min_value = 1L,
    allow_null = FALSE
) {
  if (is.null(x) && isTRUE(allow_null)) {
    return(NULL)
  }

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


.pp_functional_validate_grid <- function(grid) {
  if (!is.numeric(grid) || length(grid) < 2L || any(!is.finite(grid))) {
    stop(
      "`grid` must be a finite numeric vector containing at least two values.",
      call. = FALSE
    )
  }

  grid <- as.numeric(grid)

  if (is.unsorted(grid, strictly = TRUE)) {
    stop("`grid` must be strictly increasing.", call. = FALSE)
  }

  grid
}


.pp_functional_trapezoid_weights <- function(grid) {
  grid <- .pp_functional_validate_grid(grid)
  increments <- diff(grid)
  m <- length(grid)

  weights <- numeric(m)
  weights[1L] <- increments[1L] / 2
  weights[m] <- increments[length(increments)] / 2

  if (m > 2L) {
    weights[2L:(m - 1L)] <-
      (increments[seq_len(m - 2L)] + increments[2L:(m - 1L)]) / 2
  }

  weights
}


.pp_functional_as_curve_matrix <- function(curves, grid, na_action, training_mean = NULL) {
  if (!is.matrix(curves) && !is.data.frame(curves)) {
    stop("`curves` must be a numeric matrix or data frame.", call. = FALSE)
  }

  if (is.data.frame(curves)) {
    numeric_columns <- vapply(curves, is.numeric, logical(1))

    if (!all(numeric_columns)) {
      bad <- names(curves)[!numeric_columns]
      stop(
        sprintf(
          "Every column of `curves` must be numeric. Non-numeric column(s): %s.",
          paste(bad, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  curves <- as.matrix(curves)
  storage.mode(curves) <- "double"

  if (nrow(curves) == 0L) {
    stop("`curves` must contain at least one curve.", call. = FALSE)
  }

  if (ncol(curves) != length(grid)) {
    stop(
      sprintf(
        "`curves` must have one column for each of the %d values in `grid`.",
        length(grid)
      ),
      call. = FALSE
    )
  }

  if (any(is.infinite(curves))) {
    stop("`curves` cannot contain infinite values.", call. = FALSE)
  }

  na_action <- match.arg(na_action, c("fail", "mean"))

  if (anyNA(curves)) {
    if (identical(na_action, "fail")) {
      stop(
        "`curves` contains missing values; use `na_action = \"mean\"` to impute them.",
        call. = FALSE
      )
    }

    imputation_mean <- training_mean

    if (is.null(imputation_mean)) {
      all_missing <- colSums(!is.na(curves)) == 0L

      if (any(all_missing)) {
        stop(
          sprintf(
            "Mean imputation is impossible because curve-grid column(s) %s are entirely missing.",
            paste(which(all_missing), collapse = ", ")
          ),
          call. = FALSE
        )
      }

      imputation_mean <- colMeans(curves, na.rm = TRUE)
    }

    if (!is.numeric(imputation_mean) ||
        length(imputation_mean) != ncol(curves) ||
        any(!is.finite(imputation_mean))) {
      stop("The mean function used for imputation is invalid.", call. = FALSE)
    }

    missing_index <- which(is.na(curves), arr.ind = TRUE)
    curves[missing_index] <- imputation_mean[missing_index[, "col"]]
  }

  curves
}


.pp_functional_resolve_columns <- function(data, columns, arg_name) {
  if (is.character(columns)) {
    if (length(columns) == 0L || anyNA(columns) || any(!nzchar(columns))) {
      stop(sprintf("`%s` must contain non-empty column names.", arg_name), call. = FALSE)
    }

    missing <- setdiff(columns, names(data))

    if (length(missing) > 0L) {
      stop(
        sprintf(
          "`%s` contains unknown column(s): %s.",
          arg_name,
          paste(missing, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    return(match(columns, names(data)))
  }

  if (is.numeric(columns)) {
    if (length(columns) == 0L || any(!is.finite(columns)) ||
        any(columns != as.integer(columns))) {
      stop(sprintf("`%s` must contain valid integer column indices.", arg_name), call. = FALSE)
    }

    columns <- as.integer(columns)

    if (any(columns < 1L | columns > ncol(data))) {
      stop(sprintf("`%s` contains indices outside `support`.", arg_name), call. = FALSE)
    }

    return(columns)
  }

  stop(
    sprintf("`%s` must contain column names or integer column positions.", arg_name),
    call. = FALSE
  )
}


.pp_functional_support_spec <- function(
    support,
    coords,
    d,
    axis_roles,
    axis_units,
    time_unit,
    time_origin,
    window,
    crs
) {
  if (inherits(support, "pp")) {
    if (!is.null(coords) || !is.null(d)) {
      stop("Do not supply `coords` or `d` when `support` is a `pp` object.", call. = FALSE)
    }

    supplied_metadata <- c(
      !is.null(axis_roles), !is.null(axis_units), !is.null(time_unit),
      !is.null(time_origin), !is.null(window), !is.null(crs)
    )

    if (any(supplied_metadata)) {
      stop(
        paste0(
          "Coordinate roles, units, time metadata, window, and CRS are inherited ",
          "from a `pp` support and must not be supplied again."
        ),
        call. = FALSE
      )
    }

    if (is.null(support$coords) || !is.matrix(support$coords) ||
        is.null(support$axis_info) || !is.data.frame(support$axis_info)) {
      stop("`support` is not a valid `pp` object.", call. = FALSE)
    }

    axis_names <- as.character(support$axis_info$name)
    support_data <- as.data.frame(
      support$coords[, axis_names, drop = FALSE],
      check.names = FALSE
    )

    return(list(
      data = support_data,
      coords = axis_names,
      d = NULL,
      axis_names = axis_names,
      axis_roles = stats::setNames(support$axis_info$role, axis_names),
      axis_units = stats::setNames(support$axis_info$unit, axis_names),
      time_unit = if (is.null(support$time_info)) NULL else support$time_info$unit,
      time_origin = if (is.null(support$time_info)) NULL else support$time_info$origin,
      window = support$window,
      crs = support$crs,
      source = "pp"
    ))
  }

  if (!is.data.frame(support) && !is.matrix(support)) {
    stop("`support` must be a `pp` object, data frame, or matrix.", call. = FALSE)
  }

  support_data <- as.data.frame(support, stringsAsFactors = FALSE)

  if (ncol(support_data) == 0L) {
    stop("`support` must contain coordinate columns.", call. = FALSE)
  }

  if (is.null(names(support_data)) || any(!nzchar(names(support_data)))) {
    names(support_data) <- paste0("V", seq_len(ncol(support_data)))
  }

  if (anyDuplicated(names(support_data))) {
    stop("`support` must have unique column names.", call. = FALSE)
  }

  if (!is.null(coords) && !is.null(d)) {
    stop("Supply either `coords` or `d`, not both.", call. = FALSE)
  }

  if (is.null(coords) && is.null(d)) {
    stop("Specify support coordinate columns through `coords` or `d`.", call. = FALSE)
  }

  if (!is.null(d)) {
    d <- .pp_functional_validate_positive_integer(d, "d")

    if (d > ncol(support_data)) {
      stop("`d` cannot exceed the number of columns in `support`.", call. = FALSE)
    }

    coordinate_index <- seq_len(d)
  } else {
    coordinate_index <- .pp_functional_resolve_columns(support_data, coords, "coords")
  }

  if (anyDuplicated(coordinate_index)) {
    stop("Support coordinate columns cannot be duplicated.", call. = FALSE)
  }

  list(
    data = support_data,
    coords = if (is.null(d)) names(support_data)[coordinate_index] else NULL,
    d = if (is.null(d)) NULL else d,
    axis_names = names(support_data)[coordinate_index],
    axis_roles = axis_roles,
    axis_units = axis_units,
    time_unit = time_unit,
    time_origin = time_origin,
    window = window,
    crs = crs,
    source = "data"
  )
}


.pp_functional_component_indices <- function(object, components = NULL) {
  available <- colnames(object$scores)

  if (is.null(components)) {
    return(seq_along(available))
  }

  if (is.character(components)) {
    missing <- setdiff(components, available)

    if (length(missing) > 0L) {
      stop(
        sprintf(
          "Unknown functional component(s): %s.",
          paste(missing, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    index <- match(components, available)
  } else if (is.numeric(components)) {
    if (length(components) == 0L || any(!is.finite(components)) ||
        any(components != as.integer(components))) {
      stop("`components` must contain valid component names or integer indices.", call. = FALSE)
    }

    index <- as.integer(components)

    if (any(index < 1L | index > length(available))) {
      stop("`components` contains indices outside the retained FPCA basis.", call. = FALSE)
    }
  } else {
    stop("`components` must be NULL, character, or numeric.", call. = FALSE)
  }

  if (anyDuplicated(index)) {
    stop("`components` cannot contain duplicates.", call. = FALSE)
  }

  index
}


.pp_functional_check_object <- function(object) {
  if (!inherits(object, "pp_functional_covariate")) {
    stop("`object` must be created by pp_functional_covariate().", call. = FALSE)
  }

  required <- c(
    "name", "grid", "weights", "mean_function", "eigenfunctions",
    "scores", "score_covariates", "pve", "cumulative_pve", "settings"
  )

  missing <- setdiff(required, names(object))

  if (length(missing) > 0L) {
    stop(
      sprintf(
        "The functional-covariate object is missing component(s): %s.",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(object)
}


# Construction and accessors -------------------------------------------------

#' Construct a functional covariate for point-process modelling
#'
#' Converts curves observed at spatial or spatio-temporal support locations into
#' functional principal-component scores and represents each retained score as
#' a standard [pp_covariate()] field. The resulting score fields can be passed
#' directly to [fit_pp()]; no separate functional-model fitting function is
#' required.
#'
#' @description
#' `pp_functional_covariate()` is the data-preparation entry point for Poisson
#' point-process models with a scalar-on-function intensity component. It
#' performs a weighted functional principal component analysis (FPCA) on curves
#' observed over a common one-dimensional argument grid, retains a selected
#' number of components, and creates one `pp_covariate` object for every retained
#' score.
#'
#' @details
#' Let \eqn{X_i(u)} denote the curve observed at support location \eqn{s_i}. The
#' curves are centred by the empirical mean \eqn{\bar X(u)} and approximated as
#'
#' \deqn{X_i(u)-\bar X(u) \approx \sum_{k=1}^{K} \xi_{ik}\phi_k(u),}
#'
#' where the eigenfunctions are orthonormal under the numerical inner product
#' induced by trapezoidal integration weights. The implementation obtains the
#' FPCA through a singular-value decomposition of the weighted centred curve
#' matrix. Component signs are oriented deterministically so that the largest
#' absolute loading of each eigenfunction is positive.
#'
#' In a log-linear point-process intensity,
#'
#' \deqn{\log\lambda(s_i) = \eta_0(s_i) +
#'       \int \beta(u)\{X_i(u)-\bar X(u)\}\,du,}
#'
#' truncation at \eqn{K} components gives
#'
#' \deqn{\log\lambda(s_i) \approx \eta_0(s_i) +
#'       \sum_{k=1}^{K}\gamma_k\xi_{ik}.}
#'
#' Consequently, the retained scores are ordinary scalar covariates. Use
#' [as_pp_covariates()] to obtain their `pp_covariate` representations, include
#' all selected score names as linear terms in the formula, and estimate the
#' model with [fit_pp()]. The score fields may use any interpolation method
#' supported by [evaluate_covariate()], typically inverse-distance weighting or
#' nearest-neighbour assignment.
#'
#' The curves must share a common grid. Missing values are rejected by default.
#' With `na_action = "mean"`, missing entries are replaced by the empirical mean
#' at the corresponding grid location before FPCA. This simple imputation is
#' intended for sparse isolated missing values, not for heavily fragmented or
#' irregular functional observations.
#'
#' If `ncomp` is supplied, exactly that many components are retained. Otherwise,
#' the smallest number reaching the cumulative proportion of variance `pve` is
#' selected, subject to `max_components`. The numerical rank cannot exceed
#' `min(n - 1, m)`, where `n` is the number of curves and `m` the grid size.
#'
#' @param support Locations at which the curves are observed. Supply either a
#'   `pp` object, whose coordinates and metadata are inherited, or a data frame
#'   or matrix containing coordinate columns.
#' @param curves Numeric matrix or data frame with one row per support location
#'   and one column per value of `grid`. Row order must agree with `support`.
#' @param grid Strictly increasing finite numeric vector giving the common
#'   functional argument values associated with the columns of `curves`.
#' @param coords Coordinate-column names or positions in `support`. Required
#'   when `support` is a data frame or matrix unless `d` is supplied. Do not use
#'   with a `pp` support.
#' @param d Optional positive integer indicating that the first `d` columns of a
#'   data-frame or matrix `support` are coordinates. Cannot be combined with
#'   `coords` and must not be used with a `pp` support.
#' @param ids Optional vector of unique identifiers, one per curve. If omitted,
#'   row names of `curves` are used when available; otherwise sequential
#'   character identifiers are generated.
#' @param name Non-empty name identifying the functional covariate. It is also
#'   used to construct default score names.
#' @param score_prefix Optional syntactically valid prefix for score covariates.
#'   `NULL` uses `paste0(make.names(name), "_pc")`; retained components are then
#'   named by appending `1`, `2`, and so forth.
#' @param ncomp Optional positive integer fixing the number of retained
#'   components. When supplied, it takes precedence over `pve`.
#' @param pve Cumulative proportion of weighted curve variation targeted when
#'   `ncomp = NULL`. Must lie in `(0, 1]`.
#' @param max_components Optional positive integer limiting the number of
#'   components considered. `NULL` uses the full numerical rank.
#' @param na_action Missing-value treatment, either `"fail"` or `"mean"`.
#'   The latter imputes a missing curve value by the empirical mean at the same
#'   functional-grid location.
#' @param axis_roles Optional coordinate roles passed to [pp_covariate()] when
#'   `support` is a data frame or matrix. Entries must be `"space"` or, for at
#'   most one axis, `"time"`.
#' @param axis_units Optional coordinate units passed to [pp_covariate()] for a
#'   data-frame or matrix support.
#' @param time_unit Optional temporal unit passed to [pp_covariate()] when one
#'   support coordinate has role `"time"`.
#' @param time_origin Optional temporal origin passed to [pp_covariate()].
#' @param window Optional [pp_window()] associated with a data-frame or matrix
#'   support. It is inherited automatically from a `pp` support.
#' @param crs Optional coordinate-reference-system string for a data-frame or
#'   matrix support. It is inherited automatically from a `pp` support.
#' @param keep_curves Logical value indicating whether the completed input curve
#'   matrix should be retained in the returned object. Keeping it is useful for
#'   inspection and reconstruction but may increase object size.
#'
#' @return An object of class `"pp_functional_covariate"`. Important components
#'   include:
#'
#' * `mean_function`: estimated functional mean;
#' * `eigenfunctions`: retained weighted-orthonormal FPCA basis;
#' * `scores`: retained score matrix at the observed support;
#' * `score_covariates`: named list of standard `pp_covariate` fields;
#' * `pve` and `cumulative_pve`: proportions of variation explained;
#' * `grid` and `weights`: functional grid and trapezoidal integration weights;
#' * `curves`: completed input curves when `keep_curves = TRUE`.
#'
#' @references
#' Ramsay, J. O. and Silverman, B. W. (2005). *Functional Data Analysis*,
#' 2nd ed. Springer.
#'
#' Hall, P., Muller, H.-G. and Wang, J.-L. (2006). Properties of principal
#' component methods for functional and longitudinal data analysis. *The Annals
#' of Statistics*, **34**(3), 1493--1517.
#'
#' Reiss, P. T., Goldsmith, J., Shang, H. L. and Ogden, R. T. (2017). Methods
#' for scalar-on-function regression. *International Statistical Review*,
#' **85**(2), 228--249.
#'
#' @examples
#' set.seed(10)
#' w <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' support <- data.frame(
#'   x = runif(45),
#'   y = runif(45)
#' )
#' u <- seq(0, 1, length.out = 31)
#' score1 <- scale(support$x + support$y)[, 1]
#' score2 <- scale(support$x - support$y)[, 1]
#' curves <- outer(score1, sin(pi * u)) +
#'   0.6 * outer(score2, cos(2 * pi * u)) +
#'   matrix(rnorm(45 * length(u), sd = 0.08), 45, length(u))
#'
#' functional <- pp_functional_covariate(
#'   support = support,
#'   curves = curves,
#'   grid = u,
#'   coords = c("x", "y"),
#'   window = w,
#'   name = "profile",
#'   ncomp = 2
#' )
#' functional
#' summary(functional)
#' names(as_pp_covariates(functional))
#' head(functional_scores(functional))
#'
#' @seealso
#' [as_pp_covariates()] extracts score fields for [fit_pp()],
#' [project_functional_scores()] projects additional curves, and
#' [functional_coefficient()] reconstructs the fitted coefficient function.
#'
#' @family functional covariates
#' @export
pp_functional_covariate <- function(
    support,
    curves,
    grid,
    coords = NULL,
    d = NULL,
    ids = NULL,
    name = "functional",
    score_prefix = NULL,
    ncomp = NULL,
    pve = 0.95,
    max_components = NULL,
    na_action = c("fail", "mean"),
    axis_roles = NULL,
    axis_units = NULL,
    time_unit = NULL,
    time_origin = NULL,
    window = NULL,
    crs = NULL,
    keep_curves = TRUE
) {
  .pp_functional_require("pp_covariate", "pp_functional_covariate()")

  grid <- .pp_functional_validate_grid(grid)
  na_action <- match.arg(na_action)
  name <- .pp_functional_validate_name(name)

  if (!is.logical(keep_curves) || length(keep_curves) != 1L || is.na(keep_curves)) {
    stop("`keep_curves` must be TRUE or FALSE.", call. = FALSE)
  }

  if (is.null(score_prefix)) {
    score_prefix <- paste0(make.names(name), "_pc")
  }

  score_prefix <- .pp_functional_validate_name(score_prefix, "score_prefix")

  if (!identical(make.names(score_prefix), score_prefix)) {
    stop("`score_prefix` must be a syntactically valid R name.", call. = FALSE)
  }

  support_spec <- .pp_functional_support_spec(
    support = support,
    coords = coords,
    d = d,
    axis_roles = axis_roles,
    axis_units = axis_units,
    time_unit = time_unit,
    time_origin = time_origin,
    window = window,
    crs = crs
  )

  curves <- .pp_functional_as_curve_matrix(
    curves = curves,
    grid = grid,
    na_action = na_action
  )

  n <- nrow(curves)
  m <- ncol(curves)

  if (nrow(support_spec$data) != n) {
    stop(
      sprintf(
        "`support` has %d rows but `curves` contains %d curves.",
        nrow(support_spec$data),
        n
      ),
      call. = FALSE
    )
  }

  if (n < 2L) {
    stop("At least two curves are required for FPCA.", call. = FALSE)
  }

  if (is.null(ids)) {
    ids <- rownames(curves)

    if (is.null(ids) || any(!nzchar(ids))) {
      ids <- as.character(seq_len(n))
    }
  }

  if (length(ids) != n || anyNA(ids)) {
    stop("`ids` must contain one non-missing identifier per curve.", call. = FALSE)
  }

  ids <- as.character(ids)

  if (any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("`ids` must be unique non-empty identifiers.", call. = FALSE)
  }

  ncomp <- .pp_functional_validate_positive_integer(
    ncomp,
    "ncomp",
    allow_null = TRUE
  )

  max_components <- .pp_functional_validate_positive_integer(
    max_components,
    "max_components",
    allow_null = TRUE
  )

  pve <- .pp_functional_validate_probability(pve)

  weights <- .pp_functional_trapezoid_weights(grid)
  mean_function <- colMeans(curves)
  centred <- sweep(curves, 2L, mean_function, FUN = "-")
  weighted <- sweep(centred, 2L, sqrt(weights), FUN = "*")

  decomposition <- svd(
    weighted,
    nu = min(n, m),
    nv = min(n, m)
  )

  if (length(decomposition$d) == 0L || max(decomposition$d) <= 0) {
    stop("The supplied curves have no positive functional variation.", call. = FALSE)
  }

  tolerance <- max(dim(weighted)) * max(decomposition$d) * .Machine$double.eps
  numerical_rank <- sum(decomposition$d > tolerance)
  theoretical_rank <- min(n - 1L, m)
  available_rank <- min(numerical_rank, theoretical_rank)

  if (available_rank < 1L) {
    stop("The centred curve matrix has numerical rank zero.", call. = FALSE)
  }

  if (is.null(max_components)) {
    max_components <- available_rank
  } else if (max_components > available_rank) {
    warning(
      sprintf(
        "`max_components` was reduced from %d to the available numerical rank %d.",
        max_components,
        available_rank
      ),
      call. = FALSE
    )
    max_components <- available_rank
  }

  variation <- decomposition$d[seq_len(available_rank)]^2
  pve_all <- variation / sum(variation)
  cumulative_all <- cumsum(pve_all)

  if (!is.null(ncomp)) {
    if (ncomp > max_components) {
      stop(
        sprintf(
          "`ncomp = %d` exceeds the maximum available component count %d.",
          ncomp,
          max_components
        ),
        call. = FALSE
      )
    }

    retained <- ncomp
    selection <- "fixed"
  } else {
    reached <- which(cumulative_all[seq_len(max_components)] >= pve)

    if (length(reached) > 0L) {
      retained <- reached[1L]
    } else {
      retained <- max_components
      warning(
        sprintf(
          paste0(
            "The requested cumulative PVE %.4f was not reached within ",
            "`max_components = %d`; retaining %d components with PVE %.4f."
          ),
          pve,
          max_components,
          retained,
          cumulative_all[retained]
        ),
        call. = FALSE
      )
    }

    selection <- "pve"
  }

  component_index <- seq_len(retained)
  score_names <- paste0(score_prefix, component_index)

  if (any(score_names %in% support_spec$axis_names)) {
    conflict <- intersect(score_names, support_spec$axis_names)
    stop(
      sprintf(
        "Generated score name(s) conflict with support coordinates: %s.",
        paste(conflict, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  existing_conflict <- intersect(score_names, names(support_spec$data))

  if (length(existing_conflict) > 0L) {
    stop(
      sprintf(
        "Generated score name(s) already occur in `support`: %s. Choose another `score_prefix`.",
        paste(existing_conflict, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  scores <- decomposition$u[, component_index, drop = FALSE] %*%
    diag(decomposition$d[component_index], nrow = retained)

  eigenfunctions <- sweep(
    decomposition$v[, component_index, drop = FALSE],
    1L,
    sqrt(weights),
    FUN = "/"
  )

  for (k in component_index) {
    orientation_index <- which.max(abs(eigenfunctions[, k]))

    if (eigenfunctions[orientation_index, k] < 0) {
      eigenfunctions[, k] <- -eigenfunctions[, k]
      scores[, k] <- -scores[, k]
    }
  }

  rownames(scores) <- ids
  colnames(scores) <- score_names
  rownames(eigenfunctions) <- NULL
  colnames(eigenfunctions) <- score_names

  score_covariates <- vector("list", retained)
  names(score_covariates) <- score_names

  for (k in component_index) {
    field_data <- support_spec$data
    field_data[[score_names[k]]] <- as.numeric(scores[, k])

    score_covariates[[k]] <- pp_covariate(
      data = field_data,
      coords = support_spec$coords,
      d = support_spec$d,
      value = score_names[k],
      axis_roles = support_spec$axis_roles,
      axis_units = support_spec$axis_units,
      time_unit = support_spec$time_unit,
      time_origin = support_spec$time_origin,
      window = support_spec$window,
      crs = support_spec$crs,
      name = score_names[k]
    )
  }

  eigenvalues <- decomposition$d[component_index]^2 / (n - 1L)
  names(eigenvalues) <- score_names

  out <- list(
    name = name,
    ids = ids,
    grid = grid,
    weights = weights,
    mean_function = mean_function,
    eigenfunctions = eigenfunctions,
    scores = scores,
    score_covariates = score_covariates,
    singular_values = stats::setNames(
      decomposition$d[component_index],
      score_names
    ),
    eigenvalues = eigenvalues,
    pve = stats::setNames(pve_all[component_index], score_names),
    cumulative_pve = stats::setNames(cumulative_all[component_index], score_names),
    numerical_rank = available_rank,
    support = list(
      coords = score_covariates[[1L]]$coords,
      axis_info = score_covariates[[1L]]$axis_info,
      time_info = score_covariates[[1L]]$time_info,
      window = score_covariates[[1L]]$window,
      crs = score_covariates[[1L]]$crs
    ),
    curves = if (isTRUE(keep_curves)) curves else NULL,
    settings = list(
      selection = selection,
      requested_ncomp = ncomp,
      requested_pve = pve,
      max_components = max_components,
      retained_components = retained,
      na_action = na_action,
      score_prefix = score_prefix,
      keep_curves = isTRUE(keep_curves)
    ),
    call = match.call()
  )

  class(out) <- "pp_functional_covariate"
  out
}


#' Extract scalar score fields from a functional covariate
#'
#' Returns the retained FPCA score fields as ordinary `pp_covariate` objects.
#' These fields are designed to be supplied directly through the `covariates`
#' argument of [fit_pp()].
#'
#' @param x An object to convert.
#' @param ... Additional arguments passed to methods.
#'
#' @return A named list of `pp_covariate` objects.
#'
#' @examples
#' set.seed(4)
#' support <- data.frame(x = runif(25), y = runif(25))
#' u <- seq(0, 1, length.out = 21)
#' curves <- outer(scale(support$x)[, 1], sin(pi * u)) +
#'   matrix(rnorm(25 * length(u), sd = 0.1), 25)
#' fc <- pp_functional_covariate(
#'   support, curves, u,
#'   coords = c("x", "y"),
#'   name = "signal",
#'   ncomp = 2
#' )
#' fields <- as_pp_covariates(fc)
#' names(fields)
#'
#' @seealso [pp_functional_covariate()] and [fit_pp()].
#'
#' @family functional covariates
#' @export
as_pp_covariates <- function(x, ...) {
  UseMethod("as_pp_covariates")
}


#' @rdname as_pp_covariates
#' @param components Optional component names or one-based indices. `NULL`
#'   returns all retained score fields.
#' @method as_pp_covariates pp_functional_covariate
#' @export
as_pp_covariates.pp_functional_covariate <- function(
    x,
    components = NULL,
    ...
) {
  .pp_functional_check_object(x)
  index <- .pp_functional_component_indices(x, components)
  x$score_covariates[index]
}


#' Extract functional principal-component scores
#'
#' Retrieves the scalar FPCA scores stored at the observed functional-covariate
#' support. The returned columns have exactly the names used by
#' [as_pp_covariates()] and therefore by formulas passed to [fit_pp()].
#'
#' @param object A `pp_functional_covariate` object.
#' @param components Optional retained component names or one-based indices.
#'   `NULL` selects all components.
#' @param include_id Logical value indicating whether an initial `id` column
#'   should be included.
#'
#' @return A data frame with one row per observed curve and one column per
#'   selected score, optionally preceded by `id`.
#'
#' @examples
#' set.seed(8)
#' support <- data.frame(x = runif(20), y = runif(20))
#' u <- seq(0, 1, length.out = 15)
#' curves <- outer(rnorm(20), sin(pi * u)) +
#'   outer(rnorm(20), cos(2 * pi * u))
#' fc <- pp_functional_covariate(
#'   support, curves, u,
#'   coords = c("x", "y"),
#'   name = "curve",
#'   ncomp = 2
#' )
#' head(functional_scores(fc, include_id = TRUE))
#'
#' @family functional covariates
#' @export
functional_scores <- function(
    object,
    components = NULL,
    include_id = FALSE
) {
  .pp_functional_check_object(object)

  if (!is.logical(include_id) || length(include_id) != 1L || is.na(include_id)) {
    stop("`include_id` must be TRUE or FALSE.", call. = FALSE)
  }

  index <- .pp_functional_component_indices(object, components)
  out <- as.data.frame(object$scores[, index, drop = FALSE], check.names = FALSE)
  rownames(out) <- NULL

  if (isTRUE(include_id)) {
    out <- data.frame(id = object$ids, out, check.names = FALSE)
  }

  out
}


#' Project additional curves onto a fitted functional basis
#'
#' Centres new curves using the mean function estimated by
#' [pp_functional_covariate()] and computes their scores under the retained
#' weighted-orthonormal FPCA basis.
#'
#' @details
#' For a new curve \eqn{X_{new}(u)}, score \eqn{k} is approximated by
#'
#' \deqn{\xi_{new,k} = \sum_j w_j
#'       \{X_{new}(u_j)-\bar X(u_j)\}\phi_k(u_j),}
#'
#' where \eqn{w_j} are the trapezoidal weights stored in `object`. New curves
#' must be observed on exactly the same grid used to construct the object.
#' Missing entries treated with `na_action = "mean"` are replaced by the
#' training mean function, so their centred contribution is zero.
#'
#' @param object A `pp_functional_covariate` object defining the training mean,
#'   functional grid, integration weights, and retained eigenfunctions.
#' @param curves Numeric matrix or data frame containing additional curves, one
#'   per row and on the same grid as `object`.
#' @param components Optional retained component names or one-based indices.
#' @param ids Optional unique identifiers for the new curves. Row names are used
#'   when available; otherwise sequential identifiers are generated.
#' @param na_action Optional missing-value rule, `"fail"` or `"mean"`. `NULL`
#'   reuses the rule stored in `object`.
#' @param include_id Logical value indicating whether the returned data frame
#'   should contain an initial `id` column.
#'
#' @return A data frame containing the projected scores.
#'
#' @examples
#' set.seed(11)
#' support <- data.frame(x = runif(30), y = runif(30))
#' u <- seq(0, 1, length.out = 25)
#' curves <- outer(rnorm(30), sin(pi * u)) +
#'   outer(rnorm(30), cos(2 * pi * u))
#' fc <- pp_functional_covariate(
#'   support, curves, u,
#'   coords = c("x", "y"),
#'   ncomp = 2
#' )
#' new_curves <- curves[1:3, ] + 0.05
#' project_functional_scores(fc, new_curves, include_id = TRUE)
#'
#' @references
#' Ramsay, J. O. and Silverman, B. W. (2005). *Functional Data Analysis*,
#' 2nd ed. Springer.
#'
#' @seealso [pp_functional_covariate()] and [functional_scores()].
#'
#' @family functional covariates
#' @export
project_functional_scores <- function(
    object,
    curves,
    components = NULL,
    ids = NULL,
    na_action = NULL,
    include_id = FALSE
) {
  .pp_functional_check_object(object)

  if (is.null(na_action)) {
    na_action <- object$settings$na_action
  }

  na_action <- match.arg(na_action, c("fail", "mean"))

  if (!is.logical(include_id) || length(include_id) != 1L || is.na(include_id)) {
    stop("`include_id` must be TRUE or FALSE.", call. = FALSE)
  }

  curves <- .pp_functional_as_curve_matrix(
    curves = curves,
    grid = object$grid,
    na_action = na_action,
    training_mean = object$mean_function
  )

  n <- nrow(curves)

  if (is.null(ids)) {
    ids <- rownames(curves)

    if (is.null(ids) || any(!nzchar(ids))) {
      ids <- as.character(seq_len(n))
    }
  }

  if (length(ids) != n || anyNA(ids)) {
    stop("`ids` must contain one non-missing identifier per new curve.", call. = FALSE)
  }

  ids <- as.character(ids)

  if (any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("`ids` must be unique non-empty identifiers.", call. = FALSE)
  }

  index <- .pp_functional_component_indices(object, components)
  centred <- sweep(curves, 2L, object$mean_function, FUN = "-")
  weighted_basis <- sweep(
    object$eigenfunctions[, index, drop = FALSE],
    1L,
    object$weights,
    FUN = "*"
  )
  scores <- centred %*% weighted_basis
  colnames(scores) <- colnames(object$eigenfunctions)[index]

  out <- as.data.frame(scores, check.names = FALSE)
  rownames(out) <- NULL

  if (isTRUE(include_id)) {
    out <- data.frame(id = ids, out, check.names = FALSE)
  }

  out
}


# Coefficient-function reconstruction ----------------------------------------

#' Reconstruct a fitted functional coefficient
#'
#' Reconstructs the coefficient function associated with a functional
#' covariate from the linear coefficients of a model fitted by [fit_pp()].
#'
#' @details
#' Suppose the score terms enter the fitted log intensity linearly with
#' coefficients \eqn{\widehat\gamma_1,\ldots,\widehat\gamma_K}. The truncated
#' coefficient function is reconstructed as
#'
#' \deqn{\widehat\beta_K(u) =
#'       \sum_{k=1}^{K}\widehat\gamma_k\widehat\phi_k(u).}
#'
#' Pointwise standard errors condition on the estimated FPCA basis. If
#' \eqn{\widehat V_\gamma} denotes the covariance matrix of the score
#' coefficients, then
#'
#' \deqn{\widehat{se}\{\widehat\beta_K(u)\} =
#'       \sqrt{\boldsymbol\phi(u)^T
#'       \widehat V_\gamma\boldsymbol\phi(u)}.}
#'
#' These intervals reflect uncertainty in the point-process regression
#' coefficients but not uncertainty from estimating the functional mean and
#' eigenfunctions. Bootstrap methods are required when basis-estimation
#' uncertainty must also be represented.
#'
#' Every selected score must occur as an untransformed linear parametric term in
#' the formula supplied to [fit_pp()]. Smooth score effects, interactions, or
#' transformed terms do not define the coefficient function above and are
#' therefore rejected.
#'
#' @param fit A fitted Poisson point-process model of class `"pp_fit"` returned
#'   by [fit_pp()].
#' @param covariate A `pp_functional_covariate` object whose score fields were
#'   used in `fit`.
#' @param components Optional retained component names or one-based indices.
#'   `NULL` uses every retained component and requires all of them in `fit`.
#' @param level Confidence level for pointwise normal intervals.
#' @param interval Interval type, either `"pointwise"` or `"none"`.
#' @param unconditional Logical value passed to `vcov()` for an `mgcv::gam`
#'   object. If `TRUE`, smoothing-parameter uncertainty is included when the
#'   fitted GAM supports the unconditional covariance correction.
#'
#' @return A data frame of class `"pp_functional_coefficient"` with columns
#'   `grid`, `estimate`, `standard_error`, `lower`, and `upper`. Attributes retain
#'   the functional-covariate name, component names, confidence level, and
#'   interval type.
#'
#' @examples
#' set.seed(15)
#' w <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' locations <- data.frame(x = runif(50), y = runif(50))
#' pattern <- pp_pattern(locations, coords = c("x", "y"), window = w)
#' u <- seq(0, 1, length.out = 21)
#' curves <- outer(scale(locations$x)[, 1], sin(pi * u)) +
#'   outer(scale(locations$y)[, 1], cos(2 * pi * u)) +
#'   matrix(rnorm(50 * length(u), sd = 0.1), 50)
#' fc <- pp_functional_covariate(
#'   support = pattern,
#'   curves = curves,
#'   grid = u,
#'   name = "profile",
#'   ncomp = 2
#' )
#' fields <- as_pp_covariates(fc)
#' model_formula <- stats::reformulate(names(fields))
#' fit <- fit_pp(
#'   pattern = pattern,
#'   formula = model_formula,
#'   covariates = fields,
#'   nc = c(x = 5, y = 5),
#'   q = 2
#' )
#' beta_hat <- functional_coefficient(fit, fc)
#' head(beta_hat)
#'
#' @references
#' Reiss, P. T., Goldsmith, J., Shang, H. L. and Ogden, R. T. (2017). Methods
#' for scalar-on-function regression. *International Statistical Review*,
#' **85**(2), 228--249.
#'
#' @seealso [fit_pp()], [pp_functional_covariate()], and
#' [as_pp_covariates()].
#'
#' @family functional covariates
#' @export
functional_coefficient <- function(
    fit,
    covariate,
    components = NULL,
    level = 0.95,
    interval = c("pointwise", "none"),
    unconditional = FALSE
) {
  .pp_functional_check_object(covariate)

  if (!inherits(fit, "pp_fit")) {
    stop("`fit` must be an object returned by fit_pp().", call. = FALSE)
  }

  if (is.null(fit$gam) || !inherits(fit$gam, "gam")) {
    stop("`fit` does not contain a valid `mgcv::gam` object.", call. = FALSE)
  }

  level <- .pp_functional_validate_probability(level, "level")
  interval <- match.arg(interval)

  if (!is.logical(unconditional) || length(unconditional) != 1L || is.na(unconditional)) {
    stop("`unconditional` must be TRUE or FALSE.", call. = FALSE)
  }

  index <- .pp_functional_component_indices(covariate, components)
  component_names <- colnames(covariate$scores)[index]
  coefficients <- stats::coef(fit$gam)
  missing <- setdiff(component_names, names(coefficients))

  if (length(missing) > 0L) {
    stop(
      sprintf(
        paste0(
          "The following score coefficients were not found as linear terms in `fit`: %s. ",
          "Include the selected score names directly in the fit_pp() formula."
        ),
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  gamma <- as.numeric(coefficients[component_names])

  if (any(!is.finite(gamma))) {
    stop("Some selected score coefficients are non-finite.", call. = FALSE)
  }

  covariance <- tryCatch(
    stats::vcov(fit$gam, unconditional = unconditional),
    error = function(error) {
      stop(
        sprintf("Could not extract the fitted coefficient covariance: %s", conditionMessage(error)),
        call. = FALSE
      )
    }
  )

  if (!all(component_names %in% rownames(covariance)) ||
      !all(component_names %in% colnames(covariance))) {
    stop("The fitted covariance matrix does not contain all score coefficients.", call. = FALSE)
  }

  covariance <- covariance[component_names, component_names, drop = FALSE]
  basis <- covariate$eigenfunctions[, index, drop = FALSE]
  estimate <- as.numeric(basis %*% gamma)
  variance <- rowSums((basis %*% covariance) * basis)
  tolerance <- .Machine$double.eps^0.5

  if (any(variance < -tolerance)) {
    stop("The reconstructed coefficient variance contains negative values.", call. = FALSE)
  }

  variance <- pmax(variance, 0)
  standard_error <- sqrt(variance)

  if (identical(interval, "pointwise")) {
    critical <- stats::qnorm(1 - (1 - level) / 2)
    lower <- estimate - critical * standard_error
    upper <- estimate + critical * standard_error
  } else {
    lower <- rep(NA_real_, length(estimate))
    upper <- rep(NA_real_, length(estimate))
  }

  out <- data.frame(
    grid = covariate$grid,
    estimate = estimate,
    standard_error = standard_error,
    lower = lower,
    upper = upper
  )

  attr(out, "functional_covariate") <- covariate$name
  attr(out, "components") <- component_names
  attr(out, "level") <- level
  attr(out, "interval") <- interval
  attr(out, "unconditional") <- isTRUE(unconditional)
  class(out) <- c("pp_functional_coefficient", "data.frame")
  out
}


# Print, summary, and plot methods -------------------------------------------

#' Methods for functional covariates
#'
#' Print, summarise, and plot objects created by
#' [pp_functional_covariate()].
#'
#' @details
#' The print method gives a concise description of the number of curves, grid
#' size, retained components, and cumulative proportion of variation explained.
#' The summary method returns a structured component table. The plot method can
#' display the estimated mean function, selected eigenfunctions, or cumulative
#' proportion of variation explained.
#'
#' @param x A `pp_functional_covariate` object for `print()` and `plot()`, or a
#'   `summary_pp_functional_covariate` object for its print method.
#' @param object A `pp_functional_covariate` object to summarise.
#' @param digits Number of significant digits used for printing.
#' @param type Plot type: `"mean"`, `"components"`, or `"pve"`.
#' @param components Optional component names or indices used when
#'   `type = "components"`.
#' @param xlab,ylab,main Optional plot annotations.
#' @param ... Additional graphical arguments passed to the underlying base-R
#'   plotting function.
#'
#' @return Print and plot methods invisibly return `x`.
#' `summary.pp_functional_covariate()` returns an object of class
#' `"summary_pp_functional_covariate"`.
#'
#' @examples
#' set.seed(20)
#' support <- data.frame(x = runif(25), y = runif(25))
#' u <- seq(0, 1, length.out = 21)
#' curves <- outer(rnorm(25), sin(pi * u)) +
#'   outer(rnorm(25), cos(2 * pi * u))
#' fc <- pp_functional_covariate(
#'   support, curves, u,
#'   coords = c("x", "y"),
#'   name = "curve",
#'   ncomp = 2
#' )
#' fc
#' summary(fc)
#' plot(fc, type = "components")
#'
#' @family functional covariates
#' @rdname pp_functional_covariate-methods
#' @method print pp_functional_covariate
#' @export
print.pp_functional_covariate <- function(x, digits = 4, ...) {
  .pp_functional_check_object(x)
  cat("Functional covariate: ", x$name, "\n", sep = "")
  cat("Observed curves: ", nrow(x$scores), "\n", sep = "")
  cat("Functional grid points: ", length(x$grid), "\n", sep = "")
  cat("Retained components: ", ncol(x$scores), "\n", sep = "")
  cat(
    "Cumulative PVE: ",
    formatC(tail(x$cumulative_pve, 1L), digits = digits, format = "fg"),
    "\n",
    sep = ""
  )
  cat("Score fields: ", paste(colnames(x$scores), collapse = ", "), "\n", sep = "")
  invisible(x)
}


#' @rdname pp_functional_covariate-methods
#' @method summary pp_functional_covariate
#' @export
summary.pp_functional_covariate <- function(object, ...) {
  .pp_functional_check_object(object)

  component_table <- data.frame(
    component = colnames(object$scores),
    singular_value = as.numeric(object$singular_values),
    eigenvalue = as.numeric(object$eigenvalues),
    pve = as.numeric(object$pve),
    cumulative_pve = as.numeric(object$cumulative_pve),
    row.names = NULL
  )

  out <- list(
    name = object$name,
    n_curves = nrow(object$scores),
    grid_size = length(object$grid),
    grid_range = range(object$grid),
    retained_components = ncol(object$scores),
    numerical_rank = object$numerical_rank,
    selection = object$settings$selection,
    requested_pve = object$settings$requested_pve,
    na_action = object$settings$na_action,
    axis_info = object$support$axis_info,
    component_table = component_table
  )

  class(out) <- "summary_pp_functional_covariate"
  out
}


#' @rdname pp_functional_covariate-methods
#' @method print summary_pp_functional_covariate
#' @export
print.summary_pp_functional_covariate <- function(x, digits = 4, ...) {
  cat("Functional covariate: ", x$name, "\n", sep = "")
  cat("Observed curves: ", x$n_curves, "\n", sep = "")
  cat("Functional grid size: ", x$grid_size, "\n", sep = "")
  cat(
    "Functional grid range: [",
    formatC(x$grid_range[1L], digits = digits, format = "fg"),
    ", ",
    formatC(x$grid_range[2L], digits = digits, format = "fg"),
    "]\n",
    sep = ""
  )
  cat("Numerical rank: ", x$numerical_rank, "\n", sep = "")
  cat("Retained components: ", x$retained_components, "\n", sep = "")
  cat("Selection rule: ", x$selection, "\n", sep = "")
  cat("Missing-value rule: ", x$na_action, "\n", sep = "")
  cat("\nComponents:\n")

  table_to_print <- x$component_table
  numeric_columns <- vapply(table_to_print, is.numeric, logical(1))
  table_to_print[numeric_columns] <- lapply(
    table_to_print[numeric_columns],
    function(z) signif(z, digits)
  )
  print(table_to_print, row.names = FALSE)
  invisible(x)
}


#' @rdname pp_functional_covariate-methods
#' @method plot pp_functional_covariate
#' @export
plot.pp_functional_covariate <- function(
    x,
    type = c("mean", "components", "pve"),
    components = NULL,
    xlab = NULL,
    ylab = NULL,
    main = NULL,
    ...
) {
  .pp_functional_check_object(x)
  type <- match.arg(type)

  if (identical(type, "mean")) {
    if (is.null(xlab)) xlab <- "Functional argument"
    if (is.null(ylab)) ylab <- "Mean function"
    if (is.null(main)) main <- paste("Mean of", x$name)

    graphics::plot(
      x$grid,
      x$mean_function,
      type = "l",
      xlab = xlab,
      ylab = ylab,
      main = main,
      ...
    )
  } else if (identical(type, "components")) {
    index <- .pp_functional_component_indices(x, components)
    if (is.null(xlab)) xlab <- "Functional argument"
    if (is.null(ylab)) ylab <- "Eigenfunction"
    if (is.null(main)) main <- paste("FPCA components of", x$name)

    graphics::matplot(
      x$grid,
      x$eigenfunctions[, index, drop = FALSE],
      type = "l",
      lty = seq_along(index),
      xlab = xlab,
      ylab = ylab,
      main = main,
      ...
    )
    graphics::legend(
      "topright",
      legend = colnames(x$eigenfunctions)[index],
      lty = seq_along(index),
      bty = "n"
    )
  } else {
    if (is.null(xlab)) xlab <- "Number of components"
    if (is.null(ylab)) ylab <- "Cumulative PVE"
    if (is.null(main)) main <- paste("Variation explained by", x$name)

    graphics::plot(
      seq_along(x$cumulative_pve),
      as.numeric(x$cumulative_pve),
      type = "b",
      ylim = c(0, 1),
      xlab = xlab,
      ylab = ylab,
      main = main,
      ...
    )
  }

  invisible(x)
}


#' Methods for reconstructed functional coefficients
#'
#' Print and plot objects returned by [functional_coefficient()].
#'
#' @param x A `pp_functional_coefficient` object.
#' @param digits Number of significant digits used when printing ranges.
#' @param xlab,ylab,main Optional plot annotations.
#' @param show_interval Logical value indicating whether available pointwise
#'   confidence limits should be displayed.
#' @param interval_col Fill colour for the confidence band.
#' @param interval_border Border colour for the confidence band; `NA` omits it.
#' @param ... Additional graphical arguments passed to `graphics::plot()`.
#'
#' @return Both methods invisibly return `x`.
#'
#' @examples
#' \donttest{
#' set.seed(31)
#' window <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' locations <- data.frame(x = runif(40), y = runif(40))
#' pattern <- pp_pattern(
#'   locations,
#'   coords = c("x", "y"),
#'   window = window
#' )
#' u <- seq(0, 1, length.out = 17)
#' curves <- outer(scale(locations$x)[, 1], sin(pi * u)) +
#'   outer(scale(locations$y)[, 1], cos(2 * pi * u)) +
#'   matrix(rnorm(40 * length(u), sd = 0.08), 40)
#' functional_covariate <- pp_functional_covariate(
#'   support = pattern,
#'   curves = curves,
#'   grid = u,
#'   name = "profile",
#'   ncomp = 2
#' )
#' score_fields <- as_pp_covariates(functional_covariate)
#' fitted_model <- fit_pp(
#'   pattern = pattern,
#'   formula = stats::reformulate(names(score_fields)),
#'   covariates = score_fields,
#'   nc = c(x = 4, y = 4),
#'   q = 2
#' )
#' coefficient <- functional_coefficient(
#'   fitted_model,
#'   functional_covariate
#' )
#' print(coefficient)
#' plot(coefficient)
#' }
#'
#' @family functional covariates
#' @rdname pp_functional_coefficient-methods
#' @method print pp_functional_coefficient
#' @export
print.pp_functional_coefficient <- function(x, digits = 4, ...) {
  name <- attr(x, "functional_covariate")
  components <- attr(x, "components")
  interval <- attr(x, "interval")
  level <- attr(x, "level")

  cat("Functional coefficient: ", name, "\n", sep = "")
  cat("Components: ", paste(components, collapse = ", "), "\n", sep = "")
  cat(
    "Estimate range: [",
    formatC(min(x$estimate), digits = digits, format = "fg"),
    ", ",
    formatC(max(x$estimate), digits = digits, format = "fg"),
    "]\n",
    sep = ""
  )

  if (identical(interval, "pointwise")) {
    cat(
      "Pointwise confidence level: ",
      formatC(level, digits = digits, format = "fg"),
      "\n",
      sep = ""
    )
  } else {
    cat("Confidence interval: none\n")
  }

  invisible(x)
}


#' @rdname pp_functional_coefficient-methods
#' @method plot pp_functional_coefficient
#' @export
plot.pp_functional_coefficient <- function(
    x,
    xlab = "Functional argument",
    ylab = "Coefficient",
    main = NULL,
    show_interval = TRUE,
    interval_col = "grey85",
    interval_border = NA,
    ...
) {
  if (!inherits(x, "pp_functional_coefficient")) {
    stop("`x` must be created by functional_coefficient().", call. = FALSE)
  }

  if (!is.logical(show_interval) || length(show_interval) != 1L || is.na(show_interval)) {
    stop("`show_interval` must be TRUE or FALSE.", call. = FALSE)
  }

  if (is.null(main)) {
    main <- paste("Functional coefficient:", attr(x, "functional_covariate"))
  }

  has_interval <- isTRUE(show_interval) &&
    all(is.finite(x$lower)) && all(is.finite(x$upper))

  y_range <- range(
    if (has_interval) c(x$lower, x$upper) else x$estimate,
    finite = TRUE
  )

  graphics::plot(
    x$grid,
    x$estimate,
    type = "n",
    ylim = y_range,
    xlab = xlab,
    ylab = ylab,
    main = main,
    ...
  )

  if (has_interval) {
    graphics::polygon(
      c(x$grid, rev(x$grid)),
      c(x$lower, rev(x$upper)),
      col = interval_col,
      border = interval_border
    )
  }

  graphics::lines(x$grid, x$estimate, lwd = 2)
  graphics::abline(h = 0, lty = 3)
  invisible(x)
}
