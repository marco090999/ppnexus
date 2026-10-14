# Parametric-background ETAS models ------------------------------------------
#
# This module implements a two-dimensional spatial ETAS model with a
# temporally homogeneous background component estimated through the package's
# weighted Poisson point-process machinery. The triggering likelihood follows
# the parametrisation used by etasFLP::etasclass and evaluates the eventwise
# triggering contribution through a Fortran routine adapted from etasFLP
# (src/etas_trigger.f90).
#
# Required modules:
#   ppdata.R
#   ppfit.R   (with event_weights support in pp_cubature() and fit_pp())
#
# Compiled code:
#   src/etas_trigger.f90 (adapted from etasFLP, GPL >= 2)
#
# The ETAS model is intentionally specialised: it requires two spatial axes,
# one temporal axis, a depth mark, and a magnitude mark. The Poisson background
# is fitted only on the working spatial plane x,y. Other package models remain
# available in arbitrary spatial and spatio-temporal dimension.

.pp_etas_version <- "0.1.0"


# Shared validation -----------------------------------------------------------

.pp_etas_require <- function(functions, context) {
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
        "%s requires function(s) not currently available: %s. Source ppdata.R and ppfit.R before ppetas.R.",
        context,
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.pp_etas_validate_logical <- function(x, arg_name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", arg_name), call. = FALSE)
  }

  isTRUE(x)
}


.pp_etas_validate_positive_integer <- function(x, arg_name, min_value = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != as.integer(x) || x < min_value) {
    stop(
      sprintf("`%s` must be an integer greater than or equal to %d.", arg_name, min_value),
      call. = FALSE
    )
  }

  as.integer(x)
}


.pp_etas_validate_probability <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 || x > 1) {
    stop(sprintf("`%s` must be a finite number in [0, 1].", arg_name), call. = FALSE)
  }

  as.numeric(x)
}


.pp_etas_validate_numeric_scalar <- function(
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


.pp_etas_validate_range <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) != 2L || any(!is.finite(x))) {
    stop(sprintf("`%s` must be a finite numeric vector of length two.", arg_name), call. = FALSE)
  }

  x <- as.numeric(x)

  if (x[1L] >= x[2L]) {
    stop(sprintf("`%s[1]` must be strictly smaller than `%s[2]`.", arg_name, arg_name), call. = FALSE)
  }

  x
}


.pp_etas_validate_mark_name <- function(pattern, mark_name, arg_name) {
  if (!is.character(mark_name) || length(mark_name) != 1L ||
      is.na(mark_name) || !nzchar(mark_name)) {
    stop(sprintf("`%s` must name one mark stored in the input pattern.", arg_name), call. = FALSE)
  }

  if (is.null(pattern$marks) || !mark_name %in% names(pattern$marks)) {
    stop(
      sprintf("`%s = \"%s\"` is not a mark stored in `pattern`.", arg_name, mark_name),
      call. = FALSE
    )
  }

  values <- pattern$marks[[mark_name]]

  if (!is.numeric(values) || any(!is.finite(values))) {
    stop(
      sprintf("The `%s` mark `%s` must be numeric and finite.", arg_name, mark_name),
      call. = FALSE
    )
  }

  as.numeric(values)
}


.pp_etas_validate_catalog <- function(catalog) {
  if (!inherits(catalog, "pp_etas_catalog")) {
    stop("`catalog` must be created by etas_catalog().", call. = FALSE)
  }

  required <- c(
    "data", "marks", "window_original", "window_working",
    "longlat.to.km", "sectoday", "axis_map", "source_pattern"
  )

  missing <- setdiff(required, names(catalog))

  if (length(missing) > 0L) {
    stop(
      sprintf(
        "The supplied ETAS catalog is missing required component(s): %s.",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  data_required <- c(
    ".pp_source_id", "time", "long", "lat", "z", "magn1",
    "time.work", "xcat.work", "ycat.work"
  )

  missing_data <- setdiff(data_required, names(catalog$data))

  if (length(missing_data) > 0L) {
    stop(
      sprintf(
        "The ETAS catalog data are missing required column(s): %s.",
        paste(missing_data, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(catalog)
}


.pp_etas_require_native <- function() {
  # The ETAS triggering routine is compiled into ppnexus (src/etas_trigger.f90);
  # no external package is required.
  invisible(TRUE)
}


.pp_etas_normalize_event_weights <- function(weights, n, arg_name = "wp") {
  if (is.null(weights)) {
    return(rep.int(1, n))
  }

  if (!is.numeric(weights) || is.matrix(weights) || is.data.frame(weights) ||
      length(weights) != n) {
    stop(
      sprintf("`%s` must be a numeric vector of length %d.", arg_name, n),
      call. = FALSE
    )
  }

  weights <- as.numeric(weights)

  if (any(!is.finite(weights)) || any(weights < 0) || any(weights > 1)) {
    stop(
      sprintf("`%s` must contain finite background probabilities in [0, 1].", arg_name),
      call. = FALSE
    )
  }

  if (sum(weights) <= 0) {
    stop(sprintf("`%s` must have positive total weight.", arg_name), call. = FALSE)
  }

  weights
}


# ETAS catalog construction ---------------------------------------------------

#' Construct an ETAS earthquake catalog
#'
#' Converts a marked spatio-temporal point pattern into the canonical catalog
#' representation used by [fit_etas_pp()]. The returned object retains both the
#' original coordinates and the working coordinates required by the ETAS
#' likelihood implemented in this package.
#'
#' @description
#' `etas_catalog()` is the data-preparation entry point for ETAS modelling in
#' `ppnexus`. The input must contain exactly two spatial coordinates, one
#' temporal coordinate, a numeric depth mark, and a numeric magnitude mark.
#' Additional marks are retained in the catalog and may subsequently be used as
#' triggering covariates or as covariates in the parametric background model.
#'
#' @details
#' ETAS models are conditional-intensity models for earthquake occurrence. At
#' time \eqn{t}, the intensity is conventionally decomposed as
#'
#' \deqn{\lambda(t, \mathbf{s} \mid \mathcal{H}_t)
#' = \mu b(\mathbf{s})
#' + \sum_{t_j < t} g(t-t_j, m_j) f(\mathbf{s}-\mathbf{s}_j, m_j),}
#'
#' where \eqn{b(\mathbf{s})} is the background spatial density and the sum is
#' the contribution of earthquakes occurring before \eqn{t}. `etas_catalog()`
#' does not fit this model; it validates and reorganizes the input data so that
#' [fit_etas_pp()] can evaluate the conditional likelihood.
#'
#' The canonical columns are:
#'
#' - `time`, `long`, and `lat`: original temporal and spatial coordinates;
#' - `z`: depth;
#' - `magn1`: magnitude;
#' - `time.work`, `xcat.work`, and `ycat.work`: coordinates on the numerical
#'   scales used by the ETAS likelihood;
#' - `.pp_source_id`: row identifier inherited from the input pattern.
#'
#' When `longlat.to.km = TRUE`, both spatial coordinates are multiplied by
#' `6371.3 * pi / 180`, following the convention used by `etasFLP`. This is a
#' local angular-to-distance scaling and is not a general map projection. It is
#' appropriate only when the supplied spatial coordinates are longitude and
#' latitude in degrees and the study region is sufficiently limited for this
#' approximation to be meaningful.
#'
#' When `sectoday = TRUE`, the temporal coordinate is divided by 86400 before
#' likelihood evaluation. The original values remain available in `time` and
#' the converted values are stored in `time.work`.
#'
#' Event ordering is not required at construction time. The selected catalog is
#' ordered internally by [fit_etas_pp()] before the conditional likelihood is
#' evaluated.
#'
#' @param pattern A marked spatio-temporal object inheriting from class `"pp"`
#'   and carrying a finite [pp_window()]. It must contain exactly two axes with
#'   role `"space"`, one axis with role `"time"`, and the numeric marks named
#'   by `depth_mark` and `magnitude_mark`.
#' @param long_axis,lat_axis Optional character strings identifying the two
#'   spatial axes. If omitted, the first and second spatial axes stored in
#'   `pattern$axis_info` are used. The two names must be distinct and must
#'   identify all spatial axes of `pattern`.
#' @param time_axis Optional character string identifying the temporal axis. If
#'   omitted, the unique axis with role `"time"` is used.
#' @param depth_mark Character string naming the numeric finite depth mark. The
#'   default is `"z"`.
#' @param magnitude_mark Character string naming the numeric finite magnitude
#'   mark. The default is `"magn1"`.
#' @param longlat.to.km Logical scalar. If `TRUE`, the two original spatial
#'   coordinates are multiplied by `6371.3 * pi / 180` to obtain the working
#'   coordinates used in the likelihood. If `FALSE`, the original numerical
#'   scale is retained.
#' @param sectoday Logical scalar. If `TRUE`, original numerical times are
#'   interpreted as seconds and divided by 86400 to obtain working times in
#'   days. If `FALSE`, the original numerical scale is retained.
#'
#' @return An object of class `"pp_etas_catalog"`. Its main components are:
#'
#' - `data`: canonical event data and any additional retained marks;
#' - `marks`: the complete mark data frame inherited from `pattern`;
#' - `source_pattern`: the original point pattern;
#' - `window_original`: temporal and spatial bounds on the original scales;
#' - `window_working`: corresponding bounds on the ETAS working scales;
#' - `axis_map`: names of the source axes and marks;
#' - `axis_units`: original and working units;
#' - `longlat.to.km`, `sectoday`, `spatial_scale`, and `time_scale`: conversion
#'   metadata.
#'
#' @examples
#' set.seed(42)
#' n <- 30
#'
#' w <- pp_window(
#'   lower = c(long = 0, lat = 0, time = 0),
#'   upper = c(long = 1, lat = 1, time = 100),
#'   axis_roles = c(long = "space", lat = "space", time = "time"),
#'   axis_units = c(long = "degrees", lat = "degrees", time = "days")
#' )
#'
#' earthquakes <- pp_pattern(
#'   data = data.frame(
#'     long = runif(n),
#'     lat = runif(n),
#'     time = sort(runif(n, 0, 100)),
#'     z = runif(n, 1, 20),
#'     magn1 = 2.5 + rexp(n, rate = 2)
#'   ),
#'   coords = c("long", "lat", "time"),
#'   marks = c("z", "magn1"),
#'   window = w
#' )
#'
#' catalog <- etas_catalog(
#'   earthquakes,
#'   longlat.to.km = FALSE,
#'   sectoday = FALSE
#' )
#' catalog
#' summary(catalog)
#'
#' @references
#' Ogata, Y. (1988). Statistical models for earthquake occurrences and residual
#' analysis for point processes. *Journal of the American Statistical
#' Association*, 83(401), 9--27.
#'
#' Chiodi, M. and Adelfio, G. (2017). Mixed non-parametric and parametric
#' estimation techniques in R package etasFLP for earthquakes' description.
#' *Journal of Statistical Software*, 76(3), 1--29.
#'
#' @seealso
#' [pp_pattern()] and [pp_window()] for constructing the input data,
#' [as_etas_catalog()] for coercion, [fit_etas_pp()] for estimation,
#' [decluster_etas_pp()] for posterior event classification, and
#' [simulate_etas_pp()] for ETAS simulation.
#'
#' @family ETAS modelling
#' @export
etas_catalog <- function(
    pattern,
    long_axis = NULL,
    lat_axis = NULL,
    time_axis = NULL,
    depth_mark = "z",
    magnitude_mark = "magn1",
    longlat.to.km = TRUE,
    sectoday = FALSE
) {
  .pp_etas_require(c("pp_pattern", "pp_window"), "etas_catalog()")

  if (!inherits(pattern, "pp")) {
    stop("`pattern` must inherit from class 'pp'.", call. = FALSE)
  }

  if (is.null(pattern$window) || !inherits(pattern$window, "pp_window")) {
    stop("`pattern` must carry a `pp_window`.", call. = FALSE)
  }

  longlat.to.km <- .pp_etas_validate_logical(longlat.to.km, "longlat.to.km")
  sectoday <- .pp_etas_validate_logical(sectoday, "sectoday")

  axis_info <- pattern$axis_info
  spatial_axes <- axis_info$name[axis_info$role == "space"]
  temporal_axes <- axis_info$name[axis_info$role == "time"]

  if (length(spatial_axes) != 2L || length(temporal_axes) != 1L) {
    stop(
      paste0(
        "ETAS requires exactly two spatial axes and one temporal axis. ",
        "The input pattern has ", length(spatial_axes), " spatial and ",
        length(temporal_axes), " temporal axes."
      ),
      call. = FALSE
    )
  }

  if (is.null(long_axis)) {
    long_axis <- spatial_axes[1L]
  }

  if (is.null(lat_axis)) {
    lat_axis <- spatial_axes[2L]
  }

  if (is.null(time_axis)) {
    time_axis <- temporal_axes[1L]
  }

  requested_axes <- c(long_axis, lat_axis, time_axis)

  if (!is.character(requested_axes) || length(requested_axes) != 3L ||
      anyNA(requested_axes) || any(!nzchar(requested_axes)) ||
      anyDuplicated(requested_axes)) {
    stop(
      "`long_axis`, `lat_axis`, and `time_axis` must identify three distinct pattern axes.",
      call. = FALSE
    )
  }

  if (!identical(c(long_axis, lat_axis), spatial_axes) &&
      !setequal(c(long_axis, lat_axis), spatial_axes)) {
    stop("`long_axis` and `lat_axis` must identify the two spatial axes.", call. = FALSE)
  }

  if (!identical(time_axis, temporal_axes[1L])) {
    stop("`time_axis` must identify the temporal axis.", call. = FALSE)
  }

  depth <- .pp_etas_validate_mark_name(pattern, depth_mark, "depth_mark")
  magnitude <- .pp_etas_validate_mark_name(pattern, magnitude_mark, "magnitude_mark")

  coords <- pattern$coords
  bounds <- pattern$window$bounds

  long <- as.numeric(coords[, long_axis])
  lat <- as.numeric(coords[, lat_axis])
  time <- as.numeric(coords[, time_axis])

  if (any(!is.finite(long)) || any(!is.finite(lat)) || any(!is.finite(time))) {
    stop("ETAS coordinate values must be finite.", call. = FALSE)
  }

  spatial_scale <- if (isTRUE(longlat.to.km)) 6371.3 * pi / 180 else 1
  time_scale <- if (isTRUE(sectoday)) 86400 else 1

  x_work <- spatial_scale * long
  y_work <- spatial_scale * lat
  time_work <- time / time_scale

  mark_data <- if (is.null(pattern$marks)) {
    data.frame(row.names = seq_len(nrow(coords)))
  } else {
    as.data.frame(pattern$marks, stringsAsFactors = FALSE, check.names = FALSE)
  }

  extra_names <- setdiff(names(mark_data), c(depth_mark, magnitude_mark))

  catalog_data <- data.frame(
    .pp_source_id = seq_len(nrow(coords)),
    time = time,
    long = long,
    lat = lat,
    z = depth,
    magn1 = magnitude,
    time.work = time_work,
    xcat.work = x_work,
    ycat.work = y_work,
    check.names = FALSE
  )

  if (length(extra_names) > 0L) {
    catalog_data <- cbind(
      catalog_data,
      mark_data[, extra_names, drop = FALSE]
    )
  }

  original_bounds <- rbind(
    lower = c(
      time = bounds["lower", time_axis],
      long = bounds["lower", long_axis],
      lat = bounds["lower", lat_axis]
    ),
    upper = c(
      time = bounds["upper", time_axis],
      long = bounds["upper", long_axis],
      lat = bounds["upper", lat_axis]
    )
  )

  working_bounds <- rbind(
    lower = c(
      time = original_bounds["lower", "time"] / time_scale,
      x = spatial_scale * original_bounds["lower", "long"],
      y = spatial_scale * original_bounds["lower", "lat"]
    ),
    upper = c(
      time = original_bounds["upper", "time"] / time_scale,
      x = spatial_scale * original_bounds["upper", "long"],
      y = spatial_scale * original_bounds["upper", "lat"]
    )
  )

  units <- axis_info$unit
  names(units) <- axis_info$name

  out <- list(
    call = match.call(),
    data = catalog_data,
    marks = mark_data,
    source_pattern = pattern,
    window_original = original_bounds,
    window_working = working_bounds,
    axis_map = list(
      long_axis = long_axis,
      lat_axis = lat_axis,
      time_axis = time_axis,
      depth_mark = depth_mark,
      magnitude_mark = magnitude_mark
    ),
    axis_units = list(
      long = units[[long_axis]],
      lat = units[[lat_axis]],
      time = units[[time_axis]],
      working_space = if (isTRUE(longlat.to.km)) "km" else units[[long_axis]],
      working_time = if (isTRUE(sectoday)) "days" else units[[time_axis]]
    ),
    longlat.to.km = longlat.to.km,
    sectoday = sectoday,
    spatial_scale = spatial_scale,
    time_scale = time_scale,
    crs = pattern$crs
  )

  class(out) <- "pp_etas_catalog"
  out
}


#' Coerce an object to an ETAS catalog
#'
#' Validates an existing `pp_etas_catalog` object or converts an eligible
#' marked spatio-temporal point pattern using [etas_catalog()].
#'
#' @param x An object of class `"pp_etas_catalog"`, which is returned unchanged
#'   after validation, or an eligible object inheriting from class `"pp"`.
#' @param ... Additional arguments passed to [etas_catalog()] when `x` is a
#'   point pattern. They are ignored when `x` is already an ETAS catalog.
#'
#' @return A validated object of class `"pp_etas_catalog"`.
#'
#' @examples
#' set.seed(1)
#' w <- pp_window(
#'   lower = c(x = 0, y = 0, time = 0),
#'   upper = c(x = 1, y = 1, time = 20),
#'   axis_roles = c(x = "space", y = "space", time = "time")
#' )
#'
#' pattern <- pp_pattern(
#'   data = data.frame(
#'     x = runif(10),
#'     y = runif(10),
#'     time = sort(runif(10, 0, 20)),
#'     depth = runif(10, 1, 10),
#'     magnitude = runif(10, 2.5, 4)
#'   ),
#'   coords = c("x", "y", "time"),
#'   marks = c("depth", "magnitude"),
#'   window = w
#' )
#'
#' catalog <- as_etas_catalog(
#'   pattern,
#'   depth_mark = "depth",
#'   magnitude_mark = "magnitude",
#'   longlat.to.km = FALSE
#' )
#'
#' identical(as_etas_catalog(catalog), catalog)
#'
#' @seealso [etas_catalog()] for all conversion rules and returned components.
#'
#' @family ETAS modelling
#' @export
as_etas_catalog <- function(x, ...) {
  if (inherits(x, "pp_etas_catalog")) {
    .pp_etas_validate_catalog(x)
    return(x)
  }

  if (inherits(x, "pp")) {
    return(etas_catalog(x, ...))
  }

  stop("`x` must be a `pp_etas_catalog` or an eligible `pp` object.", call. = FALSE)
}


#' Print and summarize ETAS catalogs
#'
#' Methods for compactly displaying an object returned by [etas_catalog()] and
#' summarizing its observation windows, conversion settings, magnitudes, and
#' depths.
#'
#' @param x An object of class `"pp_etas_catalog"` for `print()`, or an object
#'   of class `"summary_pp_etas_catalog"` for printing a summary.
#' @param object An object of class `"pp_etas_catalog"`.
#' @param ... Additional arguments. They are currently unused and are included
#'   for compatibility with the base `print()` and `summary()` generics.
#'
#' @return `print.pp_etas_catalog()` and
#'   `print.summary_pp_etas_catalog()` invisibly return `x`.
#'   `summary.pp_etas_catalog()` returns an object of class
#'   `"summary_pp_etas_catalog"` containing the number of events, original and
#'   working windows, magnitude and depth summaries, and conversion flags.
#'
#' @examples
#' set.seed(2)
#' w <- pp_window(
#'   lower = c(x = 0, y = 0, time = 0),
#'   upper = c(x = 1, y = 1, time = 10),
#'   axis_roles = c(x = "space", y = "space", time = "time")
#' )
#' pattern <- pp_pattern(
#'   data = data.frame(
#'     x = runif(8), y = runif(8), time = sort(runif(8, 0, 10)),
#'     z = runif(8, 1, 15), magn1 = runif(8, 2.5, 4)
#'   ),
#'   coords = c("x", "y", "time"),
#'   marks = c("z", "magn1"),
#'   window = w
#' )
#' catalog <- etas_catalog(pattern, longlat.to.km = FALSE)
#' print(catalog)
#' summary(catalog)
#'
#' @seealso [etas_catalog()] and [fit_etas_pp()].
#'
#' @family ETAS modelling
#' @rdname pp_etas_catalog-methods
#' @method print pp_etas_catalog
#' @export
print.pp_etas_catalog <- function(x, ...) {
  .pp_etas_validate_catalog(x)

  cat("<pp_etas_catalog>\n")
  cat("ETAS earthquake catalog\n\n")
  cat("Events: ", nrow(x$data), "\n", sep = "")
  cat("Spatial conversion: ", if (x$longlat.to.km) "longitude/latitude to km" else "none", "\n", sep = "")
  cat("Time conversion: ", if (x$sectoday) "seconds to days" else "none", "\n", sep = "")
  cat(
    "Original time window: [",
    formatC(x$window_original["lower", "time"], digits = 6, format = "fg"),
    ", ",
    formatC(x$window_original["upper", "time"], digits = 6, format = "fg"),
    "]\n",
    sep = ""
  )

  invisible(x)
}


#' @rdname pp_etas_catalog-methods
#' @method summary pp_etas_catalog
#' @export
summary.pp_etas_catalog <- function(object, ...) {
  .pp_etas_validate_catalog(object)

  out <- list(
    n_events = nrow(object$data),
    original_window = object$window_original,
    working_window = object$window_working,
    magnitude_summary = summary(object$data$magn1),
    depth_summary = summary(object$data$z),
    longlat.to.km = object$longlat.to.km,
    sectoday = object$sectoday
  )

  class(out) <- "summary_pp_etas_catalog"
  out
}


#' @rdname pp_etas_catalog-methods
#' @method print summary_pp_etas_catalog
#' @export
print.summary_pp_etas_catalog <- function(x, ...) {
  cat("ETAS catalog summary\n\n")
  cat("Events: ", x$n_events, "\n", sep = "")
  cat("Longitude/latitude conversion: ", x$longlat.to.km, "\n", sep = "")
  cat("Seconds-to-days conversion: ", x$sectoday, "\n\n", sep = "")
  cat("Magnitude summary:\n")
  print(x$magnitude_summary)
  invisible(x)
}


# Preparation and geometry ----------------------------------------------------

.pp_etas_prepare_fit_catalog <- function(
    catalog,
    magn.threshold,
    tmax = NULL,
    long.range = NULL,
    lat.range = NULL
) {
  .pp_etas_validate_catalog(catalog)

  magn.threshold <- .pp_etas_validate_numeric_scalar(
    magn.threshold,
    "magn.threshold",
    lower = -Inf
  )

  original_bounds <- catalog$window_original

  if (is.null(tmax)) {
    tmax <- original_bounds["upper", "time"]
  }

  if (is.null(long.range)) {
    long.range <- as.numeric(original_bounds[, "long"])
  }

  if (is.null(lat.range)) {
    lat.range <- as.numeric(original_bounds[, "lat"])
  }

  tmax <- .pp_etas_validate_numeric_scalar(
    tmax,
    "tmax",
    lower = original_bounds["lower", "time"],
    inclusive = FALSE
  )

  if (tmax > original_bounds["upper", "time"]) {
    stop("`tmax` must lie inside the catalog observation window.", call. = FALSE)
  }

  long.range <- .pp_etas_validate_range(long.range, "long.range")
  lat.range <- .pp_etas_validate_range(lat.range, "lat.range")

  if (long.range[1L] < original_bounds["lower", "long"] ||
      long.range[2L] > original_bounds["upper", "long"]) {
    stop("`long.range` must lie inside the catalog observation window.", call. = FALSE)
  }

  if (lat.range[1L] < original_bounds["lower", "lat"] ||
      lat.range[2L] > original_bounds["upper", "lat"]) {
    stop("`lat.range` must lie inside the catalog observation window.", call. = FALSE)
  }

  source <- catalog$data

  selected <- source$magn1 >= magn.threshold &
    source$time <= tmax &
    source$long >= long.range[1L] &
    source$long <= long.range[2L] &
    source$lat >= lat.range[1L] &
    source$lat <= lat.range[2L]

  if (sum(selected) < 2L) {
    stop(
      "ETAS fitting requires at least two events after magnitude, time, and spatial selection.",
      call. = FALSE
    )
  }

  source_id <- source$.pp_source_id[selected]
  cat <- source[selected, , drop = FALSE]
  cat <- cat[order(cat$time.work, cat$.pp_source_id), , drop = FALSE]
  rownames(cat) <- NULL

  mark_data <- catalog$marks[source_id, , drop = FALSE]
  mark_data <- mark_data[match(cat$.pp_source_id, source_id), , drop = FALSE]
  rownames(mark_data) <- NULL

  cat$magnitude <- cat$magn1 - magn.threshold
  cat$ind <- TRUE
  cat$ord <- seq_len(nrow(cat))

  spatial_scale <- catalog$spatial_scale
  time_scale <- catalog$time_scale

  working_bounds <- rbind(
    lower = c(
      x = spatial_scale * long.range[1L],
      y = spatial_scale * lat.range[1L],
      time = catalog$window_original["lower", "time"] / time_scale
    ),
    upper = c(
      x = spatial_scale * long.range[2L],
      y = spatial_scale * lat.range[2L],
      time = tmax / time_scale
    )
  )

  if (working_bounds["lower", "time"] >= working_bounds["upper", "time"]) {
    stop("The selected ETAS temporal window has non-positive duration.", call. = FALSE)
  }

  list(
    data = cat,
    marks = mark_data,
    source_id = cat$.pp_source_id,
    magn.threshold = magn.threshold,
    tmax_original = tmax,
    long.range = long.range,
    lat.range = lat.range,
    working_bounds = working_bounds,
    time_span = working_bounds["upper", "time"] - working_bounds["lower", "time"],
    spatial_area = prod(
      working_bounds["upper", c("x", "y")] -
        working_bounds["lower", c("x", "y")]
    ),
    catalog = catalog
  )
}


.pp_etas_make_background_pattern <- function(prepared) {
  .pp_etas_require(c("pp_pattern", "pp_window"), "ETAS Poisson background")

  bounds <- prepared$working_bounds[, c("x", "y"), drop = FALSE]

  # A geographic CRS is no longer valid after the legacy longitude/latitude
  # to kilometre conversion. Retain the original CRS only for already-planar
  # working coordinates.
  background_crs <- if (isTRUE(prepared$catalog$longlat.to.km)) {
    NULL
  } else {
    prepared$catalog$crs
  }

  window <- pp_window(
    lower = bounds["lower", ],
    upper = bounds["upper", ],
    axis_roles = c(x = "space", y = "space"),
    axis_units = c(
      x = prepared$catalog$axis_units$working_space,
      y = prepared$catalog$axis_units$working_space
    ),
    crs = background_crs
  )

  data <- data.frame(
    x = prepared$data$xcat.work,
    y = prepared$data$ycat.work,
    check.names = FALSE
  )

  if (ncol(prepared$marks) > 0L) {
    data <- cbind(data, prepared$marks)
  }

  pp_pattern(
    data = data,
    coords = c("x", "y"),
    marks = if (ncol(prepared$marks) > 0L) "auto" else NULL,
    window = window,
    crs = background_crs
  )
}


.pp_etas_validate_params_ind <- function(params.ind) {
  if (!is.logical(params.ind) || length(params.ind) != 7L || anyNA(params.ind)) {
    stop("`params.ind` must be a logical vector of length seven.", call. = FALSE)
  }

  names(params.ind) <- c("mu", "k0", "c", "p", "gamma", "d", "q")
  params.ind
}


.pp_etas_prepare_trigger_model <- function(cat, formula1, offset) {
  if (!is.character(formula1) || length(formula1) != 1L ||
      is.na(formula1) || !nzchar(formula1)) {
    stop("`formula1` must be one non-empty character formula.", call. = FALSE)
  }

  trigger_formula <- tryCatch(
    stats::as.formula(formula1),
    error = function(error) {
      stop(sprintf("Invalid `formula1`: %s", conditionMessage(error)), call. = FALSE)
    }
  )

  lhs <- all.vars(stats::update(trigger_formula, . ~ 0))

  if (length(lhs) != 1L || !identical(lhs, "time")) {
    stop(
      "`formula1` must have `time` as its left-hand side, as in \"time ~ magnitude - 1\".",
      call. = FALSE
    )
  }

  trigger_formula <- stats::update(trigger_formula, . ~ . - 1)

  frame <- tryCatch(
    stats::model.frame(trigger_formula, data = cat, na.action = stats::na.fail),
    error = function(error) {
      stop(
        sprintf("Could not build the triggering model frame from `formula1`: %s", conditionMessage(error)),
        call. = FALSE
      )
    }
  )

  cov.matrix <- tryCatch(
    stats::model.matrix(trigger_formula, data = frame),
    error = function(error) {
      stop(
        sprintf("Could not build the triggering design matrix: %s", conditionMessage(error)),
        call. = FALSE
      )
    }
  )

  if (!is.matrix(cov.matrix) || nrow(cov.matrix) != nrow(cat)) {
    stop("Internal error: invalid triggering design matrix.", call. = FALSE)
  }

  offset_formula <- stats::model.offset(frame)

  if (!is.null(offset_formula)) {
    offset_value <- as.numeric(offset_formula)
  } else if (is.null(offset)) {
    offset_value <- rep.int(0, nrow(cat))
  } else if (is.numeric(offset) && length(offset) == 1L && is.finite(offset)) {
    offset_value <- rep.int(as.numeric(offset), nrow(cat))
  } else if (is.numeric(offset) && length(offset) == nrow(cat) &&
             all(is.finite(offset))) {
    offset_value <- as.numeric(offset)
  } else {
    stop(
      "`offset` must be NULL, one finite number, or a finite vector with one entry per selected event.",
      call. = FALSE
    )
  }

  list(
    formula = trigger_formula,
    cov.matrix = cov.matrix,
    offset = offset_value
  )
}


.pp_etas_angle_positive <- function(theta) {
  theta %% (2 * pi)
}


.pp_etas_cartesian_to_polar <- function(vertices, origin) {
  centered <- sweep(vertices, 2L, origin, FUN = "-")

  cbind(
    rho = sqrt(rowSums(centered^2)),
    theta = .pp_etas_angle_positive(atan2(centered[, 2L], centered[, 1L]))
  )
}


.pp_etas_triangle_c <- function(a, b, gamma) {
  sqrt(a^2 + b^2 - 2 * a * b * cos(gamma))
}


.pp_etas_triangle_gamma <- function(a, b, c) {
  denominator <- 2 * a * b
  value <- (a^2 + b^2 - c^2) / denominator
  acos(pmin(1, pmax(-1, value)))
}


.pp_etas_region_radii <- function(region, origin, ntheta) {

  # Convert rectangular vertices to polar coordinates
  # relative to the earthquake location.
  vert0 <- .pp_etas_cartesian_to_polar(
    vertices = region,
    origin = origin
  )

  # Sort vertices by polar angle.
  ind_theta <- order(vert0[, "theta"])

  vert0 <- vert0[
    ind_theta,
    ,
    drop = FALSE
  ]

  n_vert0 <- nrow(vert0)

  # Repeat the first vertex at the end to close the polygon.
  vert <- rbind(
    vert0,
    vert0[1L, , drop = FALSE]
  )

  vert_theta <- vert[, "theta"]
  vert_rho <- vert[, "rho"]

  vert_seq <- seq_len(n_vert0)

  # Angular width between consecutive vertices.
  vert_gamma <- .pp_etas_angle_positive(
    diff(vert_theta)
  )

  # Length of the polygon side between consecutive vertices.
  side_c <- sqrt(
    vert_rho[vert_seq + 1L]^2 +
      vert_rho[vert_seq]^2 -
      2 *
      vert_rho[vert_seq + 1L] *
      vert_rho[vert_seq] *
      cos(vert_gamma)
  )

  numerator <-
    vert_rho[vert_seq]^2 +
    side_c^2 -
    vert_rho[vert_seq + 1L]^2

  denominator <-
    2 *
    vert_rho[vert_seq] *
    side_c

  cos_alpha <- numerator / denominator

  # Numerical safeguard.
  cos_alpha <- pmin(
    1,
    pmax(-1, cos_alpha)
  )

  vert_alpha <- acos(cos_alpha)

  # Angular directions used for polar integration.
  theta_grid <- .pp_etas_angle_positive(
    seq(
      0,
      2 * pi,
      length.out = ntheta + 1L
    ) +
      vert_theta[1L] -
      2 * pi
  )

  theta_grid <- theta_grid[
    2L:(ntheta + 1L)
  ]

  # Identify the polygon side intersected by each angular ray.
  ind1 <- findInterval(
    theta_grid,
    vert_theta[vert_seq]
  )

  ind1[ind1 == 0L] <- max(vert_seq)

  alpha_inf <- vert_alpha[ind1]

  theta_inf <- .pp_etas_angle_positive(
    theta_grid -
      vert_theta[ind1]
  )

  rho_inf <- vert_rho[ind1]

  denominator_rho <- sin(
    theta_inf +
      alpha_inf
  )

  # Very small numerical values should not generate division by zero.
  tiny <- sqrt(.Machine$double.eps)

  if (any(abs(denominator_rho) < tiny)) {
    denominator_rho[
      abs(denominator_rho) < tiny
    ] <- tiny
  }

  rho_grid <-
    rho_inf *
    sin(alpha_inf) /
    denominator_rho

  if (
    any(!is.finite(rho_grid)) ||
    any(rho_grid <= 0)
  ) {
    stop(
      paste0(
        "ETAS spatial boundary radii are invalid. ",
        "This usually indicates an event on or outside the working spatial window."
      ),
      call. = FALSE
    )
  }

  rho_grid
}


.pp_etas_build_rho_s2 <- function(cat, bounds, ntheta) {
  ntheta <- .pp_etas_validate_positive_integer(ntheta, "ntheta", min_value = 4L)

  region <- rbind(
    c(bounds["lower", "x"], bounds["lower", "y"]),
    c(bounds["upper", "x"], bounds["lower", "y"]),
    c(bounds["upper", "x"], bounds["upper", "y"]),
    c(bounds["lower", "x"], bounds["upper", "y"])
  )

  out <- matrix(
    NA_real_,
    nrow = nrow(cat),
    ncol = ntheta
  )

  for (i in seq_len(nrow(cat))) {
    out[i, ] <- .pp_etas_region_radii(
      region = region,
      origin = c(cat$xcat.work[i], cat$ycat.work[i]),
      ntheta = ntheta
    )
  }

  out
}


.pp_etas_starting <- function(
    cat,
    time_upper,
    p.start,
    gamma.start,
    q.start,
    betacov.start,
    onlytime
) {
  n <- nrow(cat)

  if (n < 2L) {
    stop("At least two events are required to construct ETAS starting values.", call. = FALSE)
  }

  time <- cat$time.work
  x <- cat$xcat.work
  y <- cat$ycat.work

  dt <- diff(time)
  ds <- diff(x)^2 + diff(y)^2

  time_span <- max(time_upper - min(time), .Machine$double.eps)
  spatial_scale <- max(diff(range(x)), diff(range(y)), 1)

  c.start <- max(
    as.numeric(stats::quantile(dt, probs = 0.25, names = FALSE, type = 8)),
    time_span * 1e-8,
    .Machine$double.eps
  )

  d.start <- max(
    as.numeric(stats::quantile(ds, probs = 0.05, names = FALSE, type = 8)),
    spatial_scale^2 * 1e-12,
    .Machine$double.eps
  )

  mu.start <- n * 0.5 / time_span

  if (isTRUE(onlytime)) {
    gamma.start <- 0
    q.start <- 0
    d.start <- 0
  }

  temporal_integral <- .pp_etas_temporal_integral(
    c = c.start,
    p = p.start,
    t = time,
    tmax = time_upper
  )

  magnitude_factor <- exp(as.numeric(betacov.start) * cat$magnitude)

  if (isTRUE(onlytime)) {
    denominator <- sum(magnitude_factor * temporal_integral)
  } else {
    q_safe <- max(q.start, 1 + 1e-6)
    spatial_integral <- d.start^(1 - q_safe) *
      exp(gamma.start * cat$magnitude) * pi / (q_safe - 1)

    denominator <- sum(magnitude_factor * temporal_integral * spatial_integral)
  }

  k0.start <- if (is.finite(denominator) && denominator > 0) {
    n * 0.5 / denominator
  } else {
    0.01
  }

  list(
    mu.start = mu.start,
    k0.start = k0.start,
    c.start = c.start,
    p.start = p.start,
    gamma.start = gamma.start,
    d.start = d.start,
    q.start = q.start,
    betacov.start = betacov.start
  )
}


# Triggering likelihood -------------------------------------------------------

.pp_etas_native_trigger <- function(
    cat,
    params,
    predictor,
    onlytime = FALSE
) {
  .pp_etas_require_native()

  n <- nrow(cat)

  if (length(predictor) != n || any(!is.finite(predictor))) {
    stop("The ETAS triggering predictor is invalid.", call. = FALSE)
  }

  required <- c("mu", "k0", "c", "p", "gamma", "d", "q")

  if (!identical(names(params), required)) {
    params <- params[required]
  }

if (
  length(params) != 7L ||
  any(!is.finite(params))
) {
  return(NULL)
}

  result <- tryCatch(
    .Fortran(
      C_pp_etas_trigger,
      NAOK = TRUE,
      tflag = as.integer(isTRUE(onlytime)),
      n = as.integer(n),
      mu = as.double(params[["mu"]]),
      k = as.double(params[["k0"]]),
      c = as.double(params[["c"]]),
      p = as.double(params[["p"]]),
      g = as.double(params[["gamma"]]),
      d = as.double(params[["d"]]),
      q = as.double(params[["q"]]),
      x = as.double(cat$xcat.work),
      y = as.double(cat$ycat.work),
      t = as.double(cat$time.work),
      m = as.double(cat$magnitude),
      predictor = as.double(predictor),
      l = double(n)
    ),
    error = function(error) {
      stop(
        paste0(
          "The compiled ETAS triggering routine could not be called. ",
          "Original error: ",
          conditionMessage(error)
        ),
        call. = FALSE
      )
    }
  )

  triggering <- as.numeric(result$l)

  if (length(triggering) != n || any(!is.finite(triggering)) ||
      any(triggering < 0)) {
    stop("The ETAS triggering routine returned an invalid triggering-intensity vector.", call. = FALSE)
  }

  triggering
}


.pp_etas_temporal_integral <- function(c, p, t, tmax) {
  delta <- pmax(as.numeric(tmax) - as.numeric(t), 0)

  if (abs(p - 1) <= 1e-8) {
    return(log1p(delta / c))
  }

  ((c + delta)^(1 - p) - c^(1 - p)) / (1 - p)
}


.pp_etas_trigger_integral <- function(
    cat,
    params,
    predictor,
    tmax,
    rho.s2,
    ntheta,
    onlytime = FALSE
) {
  temporal <- .pp_etas_temporal_integral(
    c = params[["c"]],
    p = params[["p"]],
    t = cat$time.work,
    tmax = tmax
  )

  if (isTRUE(onlytime)) {
    return(params[["k0"]] * sum(exp(predictor) * temporal))
  }

  q <- params[["q"]]
  d <- params[["d"]]
  gamma <- params[["gamma"]]

  if (q <= 1 || d <= 0) {
    stop("Spatial ETAS integration requires q > 1 and d > 0.", call. = FALSE)
  }

  magnitude_scale <- exp(gamma * cat$magnitude)
  productivity <- exp(predictor) * magnitude_scale

  boundary_term <- sweep(
    rho.s2^2,
    1L,
    magnitude_scale,
    FUN = "/"
  ) + d

  radial_integrand <- boundary_term^(1 - q) - d^(1 - q)

  spatial_integral <- pi / ((1 - q) * ntheta) *
    rowSums(sweep(radial_integrand, 1L, productivity, FUN = "*"))

  if (any(!is.finite(spatial_integral)) || any(spatial_integral < 0)) {
    stop("The integrated ETAS spatial triggering term is invalid.", call. = FALSE)
  }

  sum(params[["k0"]] * temporal * spatial_integral)
}


.pp_etas_unpack_parameters <- function(
    candidate,
    params.fix,
    params.ind,
    params.lim,
    n_covariates
) {
  n_etas <- sum(params.ind)

  if (length(candidate) != n_etas + n_covariates) {
    stop("Internal error: ETAS optimisation vector has invalid length.", call. = FALSE)
  }

  transformed <- candidate[seq_len(n_etas)]

  params <- params.fix
  params[params.ind] <- exp(transformed) + params.lim[params.ind]

  betacov <- if (n_covariates > 0L) {
    candidate[(n_etas + 1L):(n_etas + n_covariates)]
  } else {
    numeric(0)
  }

  list(
    params = params,
    betacov = as.numeric(betacov)
  )
}


.pp_etas_objective <- function(
    candidate,
    state,
    trace = FALSE
) {

  # Penalità finita per regioni numericamente non valide.
  # Non modifica la likelihood nelle regioni ammissibili.
  penalty <- 1e50

  # -----------------------------------------------------------------------
  # 1. Il vettore candidato deve essere finito
  # -----------------------------------------------------------------------

  if (
    !is.numeric(candidate) ||
    any(!is.finite(candidate))
  ) {
    return(penalty)
  }

  # -----------------------------------------------------------------------
  # 2. Trasformazione dei parametri
  # -----------------------------------------------------------------------

  unpacked <- tryCatch(
    .pp_etas_unpack_parameters(
      candidate = candidate,
      params.fix = state$params.fix,
      params.ind = state$params.ind,
      params.lim = state$params.lim,
      n_covariates = ncol(state$cov.matrix)
    ),
    error = function(e) NULL
  )

  if (is.null(unpacked)) {
    return(penalty)
  }

  params <- unpacked$params

  # Non introduciamo nuovi vincoli teorici:
  # controlliamo soltanto che i valori ottenuti siano finiti.
  if (
    length(params) != 7L ||
    any(!is.finite(params)) ||
    any(!is.finite(unpacked$betacov))
  ) {
    return(penalty)
  }

  # -----------------------------------------------------------------------
  # 3. Predittore del triggering
  # -----------------------------------------------------------------------

  predictor <- tryCatch(
    as.vector(
      state$cov.matrix %*%
        unpacked$betacov +
        state$offset
    ),
    error = function(e) NULL
  )

  if (
    is.null(predictor) ||
    length(predictor) != nrow(state$cat) ||
    any(!is.finite(predictor))
  ) {
    return(penalty)
  }

  # -----------------------------------------------------------------------
  # 4. Componente triggering tramite routine Fortran (src/etas_trigger.f90)
  # -----------------------------------------------------------------------

  triggering <- tryCatch(
    .pp_etas_native_trigger(
      cat = state$cat,
      params = params,
      predictor = predictor,
      onlytime = state$onlytime
    ),
    error = function(e) NULL
  )

  if (
    is.null(triggering) ||
    length(triggering) != nrow(state$cat) ||
    any(!is.finite(triggering)) ||
    any(triggering < 0)
  ) {
    return(penalty)
  }

  # -----------------------------------------------------------------------
  # 5. Intensità condizionata totale
  # -----------------------------------------------------------------------

  intensity <-
    params[["mu"]] *
    state$back.dens +
    triggering

  if (
    any(!is.finite(intensity)) ||
    any(intensity <= 0)
  ) {
    return(penalty)
  }

  # -----------------------------------------------------------------------
  # 6. Integrale della componente triggering
  # -----------------------------------------------------------------------

  trigger_integral <- tryCatch(
    .pp_etas_trigger_integral(
      cat = state$cat,
      params = params,
      predictor = predictor,
      tmax = state$tmax,
      rho.s2 = state$rho.s2,
      ntheta = state$ntheta,
      onlytime = state$onlytime
    ),
    error = function(e) NA_real_
  )

  if (
    length(trigger_integral) != 1L ||
    !is.finite(trigger_integral)
  ) {
    return(penalty)
  }

  # -----------------------------------------------------------------------
  # 7. Integrale totale ETAS
  # -----------------------------------------------------------------------

  integral <-
    trigger_integral +
    params[["mu"]] *
    state$time_span *
    state$back.integral

  if (!is.finite(integral)) {
    return(penalty)
  }

  # -----------------------------------------------------------------------
  # 8. Negative log-likelihood
  # -----------------------------------------------------------------------

  negative_loglik <-
    -sum(log(intensity)) +
    integral

  if (!is.finite(negative_loglik)) {
    return(penalty)
  }

  # -----------------------------------------------------------------------
  # 9. Attributi necessari dopo l'ottimizzazione
  # -----------------------------------------------------------------------

  attr(
    negative_loglik,
    "etas.vec"
  ) <- triggering

  attr(
    negative_loglik,
    "lambda.vec"
  ) <- intensity

  attr(
    negative_loglik,
    "integraltot"
  ) <- integral

  attr(
    negative_loglik,
    "trigger.integral"
  ) <- trigger_integral

  attr(
    negative_loglik,
    "predictor"
  ) <- predictor

  attr(
    negative_loglik,
    "params"
  ) <- params

  attr(
    negative_loglik,
    "betacov"
  ) <- unpacked$betacov

  if (isTRUE(trace)) {
    cat("*")
  }

  negative_loglik
}


.pp_etas_optimize <- function(
    state,
    hessian = TRUE,
    iterlim = 50L,
    trace = TRUE
) {

  n <- nrow(state$cat)

  n_etas <- sum(
    state$params.ind
  )

  start <- c(
    log(
      (
        state$params.fix -
          state$params.lim
      )[
        state$params.ind
      ]
    ),
    state$betacov
  )

  objective <- function(par) {

    .pp_etas_objective(
      candidate = par,
      state = state,
      trace = trace
    )
  }

  if (isTRUE(state$usenlm)) {

    result <- stats::nlm(
      f = objective,
      p = start,

      # Exact legacy behaviour
      typsize = abs(start),

      hessian = hessian,
      iterlim = iterlim
    )

    parameters <- result$estimate
    value <- result$minimum

  } else {

    result <- stats::optim(
      par = start,
      fn = objective,
      method = state$method,
      hessian = hessian,
      control = list(
        trace = if (isTRUE(trace)) 2 else 0,
        maxit = iterlim,

        # Same scaling principle as original code
        fnscale =
          n /
          state$time_span,

        parscale =
          sqrt(exp(start))
      )
    )

    parameters <- result$par
    value <- result$value
  }

  list(
    parameters = as.numeric(parameters),
    value = as.numeric(value),
    result = result,
    state = state,
    n_etas = n_etas
  )
}


.pp_etas_standard_errors <- function(
    optimisation,
    params.ind,
    n_covariates
) {
  n_parameters <- sum(params.ind) + n_covariates
  out <- rep(NA_real_, n_parameters)

  hessian <- optimisation$result$hessian

  if (is.null(hessian) || !is.matrix(hessian) ||
      !identical(dim(hessian), c(n_parameters, n_parameters)) ||
      any(!is.finite(hessian))) {
    return(out)
  }

  inverse <- try(solve(hessian), silent = TRUE)

  if (inherits(inverse, "try-error")) {
    return(out)
  }

  diagonal <- diag(inverse)

  if (any(!is.finite(diagonal)) || any(diagonal < 0)) {
    return(out)
  }

  scale <- c(
    exp(optimisation$parameters[seq_len(sum(params.ind))]),
    rep.int(1, n_covariates)
  )

  sqrt(diagonal) * scale
}


# Parametric Poisson background ----------------------------------------------

.pp_etas_fit_background <- function(
    prepared,
    background_pattern,
    event_weights,
    formula.bg,
    covariates.bg,
    covariate_methods.bg,
    covariate_args.bg,
    continuous_marks.bg,
    categorical_marks.bg,
    categorical_levels.bg,
    q.bg,
    dummy_layout.bg,
    nc.bg,
    max_tiles.bg,
    max_quadrature_points.bg,
    max_mark_combinations.bg,
    seed.bg,
    gam_args.bg,
    mark_n.bg,
    force_constant = FALSE,
    onlytime = FALSE
) {
  n <- nrow(prepared$data)

  if (isTRUE(onlytime)) {
    return(list(
      fit = NULL,
      density_events = rep.int(1, n),
      integral = 1,
      kappa = 1,
      edf = 0,
      mode = "onlytime"
    ))
  }

  if (isTRUE(force_constant)) {
    density <- 1 / prepared$spatial_area

    return(list(
      fit = NULL,
      density_events = rep.int(density, n),
      integral = 1,
      kappa = 1,
      edf = 0,
      mode = "constant"
    ))
  }

  if (!is.null(categorical_marks.bg) && length(categorical_marks.bg) > 0L) {
    warning(
      paste0(
        "Categorical marks are included in the Poisson background and then ",
        "marginalised before ETAS likelihood evaluation. They do not modify ",
        "the triggering mechanism or the ETAS mark distribution."
      ),
      call. = FALSE
    )
  }

  background_fit <- fit_pp(
    pattern = background_pattern,
    formula = formula.bg,
    covariates = covariates.bg,
    covariate_methods = covariate_methods.bg,
    covariate_args = covariate_args.bg,
    continuous_marks = continuous_marks.bg,
    categorical_marks = categorical_marks.bg,
    categorical_levels = categorical_levels.bg,
    event_weights = event_weights,
    q = q.bg,
    dummy_layout = dummy_layout.bg,
    nc = nc.bg,
    max_tiles = max_tiles.bg,
    max_quadrature_points = max_quadrature_points.bg,
    max_mark_combinations = max_mark_combinations.bg,
    seed = seed.bg,
    gam_args = gam_args.bg,
    keep_cubature = "full"
  )

  if (is.null(background_fit$cubature) ||
      is.null(background_fit$cubature$data)) {
    stop(
      "Internal ETAS error: the weighted background fit must retain full cubature.",
      call. = FALSE
    )
  }

  quad_data <- background_fit$cubature$data
  lambda_quad <- exp(stats::predict(background_fit$gam, newdata = quad_data, type = "link"))
  quad_weights <- quad_data$.pp_weight

  kappa <- sum(quad_weights * lambda_quad)

  if (!is.finite(kappa) || kappa <= 0) {
    stop("Background normalisation failed: `kappa` is non-positive or non-finite.", call. = FALSE)
  }

  event_coordinates <- data.frame(
    x = prepared$data$xcat.work,
    y = prepared$data$ycat.work,
    check.names = FALSE
  )

  lambda_events <- predict_pp_fit(
    object = background_fit,
    at = event_coordinates,
    type = "ground",
    mark_n = mark_n.bg
  )

  density_events <- lambda_events / kappa
  density_events <- pmax(density_events, .Machine$double.eps)

  normalised_integral <- sum(quad_weights * (lambda_quad / kappa))

  list(
    fit = background_fit,
    density_events = density_events,
    integral = normalised_integral,
    kappa = kappa,
    edf = if (is.null(background_fit$gam$edf)) 0 else sum(background_fit$gam$edf),
    mode = "poisson"
  )
}


# Main ETAS fit ---------------------------------------------------------------

# Internal ETAS fitting engine. The documented user interface is
# fit_etas_pp() in ppinterface.R, which assembles starting values,
# background settings and etas_control() options for this engine.
.pp_fit_etas_engine <- function(
    catalog,
    time.update = FALSE,
    magn.threshold = 2.5,
    magn.threshold.back = magn.threshold + 2,
    tmax = NULL,
    long.range = NULL,
    lat.range = NULL,
    mu = 1,
    k0 = 1,
    c = 0.5,
    p = 1.01,
    gamma = 0,
    d = 1,
    q = 1.5,
    betacov = 0.7,
    params.ind = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
    formula1 = "time ~ magnitude - 1",
    offset = 0,
    hdef = c(1, 1),
    wp = NULL,
    hvarx = NULL,
    hvary = NULL,
    declustering = TRUE,
    thinning = FALSE,
    flp = FALSE,
    m1 = NULL,
    ndeclust = 5,
    n.iterweight = 1,
    onlytime = FALSE,
    is.backconstant = FALSE,
    description = "",
    cat.back = NULL,
    back.smooth = 1,
    sectoday = NULL,
    longlat.to.km = NULL,
    usenlm = TRUE,
    method = "BFGS",
    compsqm = TRUE,
    epsmax = 1e-4,
    iterlim = 50,
    ntheta = 36,
    formula.bg = ~ 1,
    covariates.bg = NULL,
    covariate_methods.bg = NULL,
    covariate_args.bg = NULL,
    continuous_marks.bg = NULL,
    categorical_marks.bg = NULL,
    categorical_levels.bg = NULL,
    q.bg = 3,
    dummy_layout.bg = c("grid", "random"),
    nc.bg = NULL,
    max_tiles.bg = 1e6,
    max_quadrature_points.bg = 1e6,
    max_mark_combinations.bg = 100L,
    seed.bg = NULL,
    gam_args.bg = list(method = "REML"),
    mark_n.bg = NULL,
    min.declust.iter = 2,
    tol.rho = 1e-3,
    tol.par = 1e-3,
    stop.on.increasing.AIC = FALSE,
    verbose.convergence = TRUE,
    verbose = TRUE,
    trace = TRUE
) {
  fit_start_time <- Sys.time()

  .pp_etas_require(
    c("pp_pattern", "pp_window", "fit_pp", "predict_pp_fit"),
    "fit_etas_pp()"
  )
  .pp_etas_require_native()
  .pp_etas_validate_catalog(catalog)

  time.update <- .pp_etas_validate_logical(time.update, "time.update")
  declustering <- .pp_etas_validate_logical(declustering, "declustering")
  thinning <- .pp_etas_validate_logical(thinning, "thinning")
  flp <- .pp_etas_validate_logical(flp, "flp")
  onlytime <- .pp_etas_validate_logical(onlytime, "onlytime")
  is.backconstant <- .pp_etas_validate_logical(is.backconstant, "is.backconstant")
  usenlm <- .pp_etas_validate_logical(usenlm, "usenlm")
  compsqm <- .pp_etas_validate_logical(compsqm, "compsqm")
  stop.on.increasing.AIC <- .pp_etas_validate_logical(
    stop.on.increasing.AIC,
    "stop.on.increasing.AIC"
  )
  verbose.convergence <- .pp_etas_validate_logical(
    verbose.convergence,
    "verbose.convergence"
  )
  verbose <- .pp_etas_validate_logical(verbose, "verbose")
  trace <- .pp_etas_validate_logical(trace, "trace")

  if (isTRUE(time.update)) {
    stop(
      "`time.update = TRUE` is not yet implemented for pp_etas_fit objects.",
      call. = FALSE
    )
  }

  if (!is.null(cat.back)) {
    stop("`cat.back` is not supported; the parametric background is fitted internally.", call. = FALSE)
  }

  if (isTRUE(flp)) {
    stop(
      paste0(
        "`flp = TRUE` requests the classical FLP kernel-background path. ",
        "ppetas.R implements the parametric Poisson-background ETAS model, ",
        "so use `flp = FALSE`."
      ),
      call. = FALSE
    )
  }

  if (isTRUE(thinning)) {
    warning(
      "`thinning = TRUE` is ignored for the parametric Poisson background; fractional posterior weights are used.",
      call. = FALSE
    )
  }

  if (!is.null(sectoday) &&
      (!is.logical(sectoday) || length(sectoday) != 1L ||
       is.na(sectoday) || !identical(isTRUE(sectoday), catalog$sectoday))) {
    stop(
      "`sectoday` must match the value used when creating `catalog` with etas_catalog().",
      call. = FALSE
    )
  }

  if (!is.null(longlat.to.km) &&
      (!is.logical(longlat.to.km) || length(longlat.to.km) != 1L ||
       is.na(longlat.to.km) || !identical(isTRUE(longlat.to.km), catalog$longlat.to.km))) {
    stop(
      "`longlat.to.km` must match the value used when creating `catalog` with etas_catalog().",
      call. = FALSE
    )
  }

  if (!inherits(formula.bg, "formula")) {
    formula.bg <- tryCatch(
      stats::as.formula(formula.bg),
      error = function(error) {
        stop(sprintf("Invalid `formula.bg`: %s", conditionMessage(error)), call. = FALSE)
      }
    )
  }

  if (length(all.vars(formula.bg)) == 0L && length(formula.bg) == 0L) {
    stop("`formula.bg` is invalid.", call. = FALSE)
  }

  dummy_layout.bg <- match.arg(dummy_layout.bg)

  ndeclust <- .pp_etas_validate_positive_integer(ndeclust, "ndeclust", min_value = 1L)
  n.iterweight <- .pp_etas_validate_positive_integer(n.iterweight, "n.iterweight", min_value = 1L)
  min.declust.iter <- .pp_etas_validate_positive_integer(
    min.declust.iter,
    "min.declust.iter",
    min_value = 1L
  )
  iterlim <- .pp_etas_validate_positive_integer(iterlim, "iterlim", min_value = 1L)
  ntheta <- .pp_etas_validate_positive_integer(ntheta, "ntheta", min_value = 4L)

  epsmax <- .pp_etas_validate_numeric_scalar(epsmax, "epsmax", lower = 0, inclusive = FALSE)
  tol.rho <- .pp_etas_validate_numeric_scalar(tol.rho, "tol.rho", lower = 0, inclusive = FALSE)
  tol.par <- .pp_etas_validate_numeric_scalar(tol.par, "tol.par", lower = 0, inclusive = FALSE)

  if (!is.character(method) || length(method) != 1L || is.na(method)) {
    stop("`method` must be one optimisation-method string.", call. = FALSE)
  }

  params.ind <- .pp_etas_validate_params_ind(params.ind)

  mu_missing <- missing(mu)
  k0_missing <- missing(k0)
  c_missing <- missing(c)
  d_missing <- missing(d)

  prepared <- .pp_etas_prepare_fit_catalog(
    catalog = catalog,
    magn.threshold = magn.threshold,
    tmax = tmax,
    long.range = long.range,
    lat.range = lat.range
  )

  cat <- prepared$data
  n <- nrow(cat)

  if (isTRUE(onlytime)) {
    is.backconstant <- TRUE
    declustering <- FALSE
    params.ind[c("gamma", "d", "q")] <- FALSE
    gamma <- 0
    d <- 0
    q <- 0
  }

  if (isTRUE(is.backconstant)) {
    declustering <- FALSE
  }

  if (!isTRUE(declustering)) {
    ndeclust <- 1L
  }

  trigger_model <- .pp_etas_prepare_trigger_model(
    cat = cat,
    formula1 = formula1,
    offset = offset
  )

  cov.matrix <- trigger_model$cov.matrix
  n_covariates <- ncol(cov.matrix)

  if (length(betacov) != n_covariates) {
    if (isTRUE(verbose)) {
      message(
        "Wrong number of `betacov` starting values; using zero for every triggering coefficient."
      )
    }

    betacov <- rep.int(0, n_covariates)
  }

  betacov <- as.numeric(betacov)

  if (any(!is.finite(betacov))) {
    stop("`betacov` must contain finite values.", call. = FALSE)
  }

  starting <- .pp_etas_starting(
    cat = cat,
    time_upper = prepared$working_bounds["upper", "time"],
    p.start = p,
    gamma.start = gamma,
    q.start = q,
    betacov.start = betacov[1L],
    onlytime = onlytime
  )

  if (isTRUE(mu_missing)) mu <- starting$mu.start
  if (isTRUE(k0_missing)) k0 <- starting$k0.start
  if (isTRUE(c_missing)) c <- starting$c.start
  if (isTRUE(d_missing) && !isTRUE(onlytime)) d <- starting$d.start

  namespar <- c("mu", "k0", "c", "p", "gamma", "d", "q")
  params.fix <- c(mu = mu, k0 = k0, c = c, p = p, gamma = gamma, d = d, q = q)
  params.lim <- c(mu = 0, k0 = 0, c = 0, p = 1, gamma = 0, d = 0, q = 1)

  if (isTRUE(onlytime)) {
    valid_fixed <- params.fix[c("mu", "k0", "c", "p")] >=
      params.lim[c("mu", "k0", "c", "p")]
  } else {
    valid_fixed <- params.fix >= params.lim
  }

  if (!all(valid_fixed)) {
    stop(
      "Invalid ETAS starting/fixed parameters relative to their lower bounds.",
      call. = FALSE
    )
  }

  if (any(params.ind & (params.fix <= params.lim))) {
    offending <- names(params.fix)[params.ind & (params.fix <= params.lim)]

    stop(
      sprintf(
        "Estimated parameter(s) must be strictly above their lower bounds: %s.",
        paste(offending, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!isTRUE(onlytime) &&
      (params.fix[["q"]] <= 1 || params.fix[["d"]] <= 0)) {
    stop(
      "The spatial ETAS model requires starting/fixed values q > 1 and d > 0.",
      call. = FALSE
    )
  }

  background_pattern <- .pp_etas_make_background_pattern(prepared)

  if (is.null(wp)) {
    background_weights <- rep.int(1, n)
  } else {
    background_weights <- .pp_etas_normalize_event_weights(wp, n, "wp")
  }

  rho.s2 <- if (isTRUE(onlytime)) {
    matrix(0, nrow = n, ncol = 1L)
  } else {
    .pp_etas_build_rho_s2(
      cat = cat,
      bounds = prepared$working_bounds,
      ntheta = ntheta
    )
  }

  fit_background <- function(weights, force_constant = FALSE) {
    .pp_etas_fit_background(
      prepared = prepared,
      background_pattern = background_pattern,
      event_weights = weights,
      formula.bg = if (isTRUE(force_constant)) ~ 1 else formula.bg,
      covariates.bg = if (isTRUE(force_constant)) NULL else covariates.bg,
      covariate_methods.bg = if (isTRUE(force_constant)) NULL else covariate_methods.bg,
      covariate_args.bg = if (isTRUE(force_constant)) NULL else covariate_args.bg,
      continuous_marks.bg = if (isTRUE(force_constant)) NULL else continuous_marks.bg,
      categorical_marks.bg = if (isTRUE(force_constant)) NULL else categorical_marks.bg,
      categorical_levels.bg = if (isTRUE(force_constant)) NULL else categorical_levels.bg,
      q.bg = q.bg,
      dummy_layout.bg = dummy_layout.bg,
      nc.bg = nc.bg,
      max_tiles.bg = max_tiles.bg,
      max_quadrature_points.bg = max_quadrature_points.bg,
      max_mark_combinations.bg = max_mark_combinations.bg,
      seed.bg = seed.bg,
      gam_args.bg = gam_args.bg,
      mark_n.bg = mark_n.bg,
      force_constant = force_constant,
      onlytime = onlytime
    )
  }

  background_state <- fit_background(
    weights = background_weights,
    force_constant = isTRUE(is.backconstant)
  )

  if (isTRUE(verbose)) {
    initial <- params.fix
    message(
      "Initial ETAS parameters: ",
      paste(
        paste0(names(initial), " = ", formatC(initial, digits = 5, format = "fg")),
        collapse = ", "
      )
    )
  }

  state <- list(
    cat = cat,
    params.fix = params.fix,
    params.ind = params.ind,
    params.lim = params.lim,
    betacov = betacov,
    cov.matrix = cov.matrix,
    offset = trigger_model$offset,
    back.dens = background_state$density_events,
    back.integral = background_state$integral,
    time_span = prepared$time_span,
    tmax = prepared$working_bounds["upper", "time"],
    rho.s2 = rho.s2,
    ntheta = ntheta,
    onlytime = onlytime,
    usenlm = usenlm,
    method = method
  )

  history_parameters <- list()
  history_se <- list()
  history_rho <- list()
  history_background_density <- list()
  history_aic_ml <- numeric(0)
  history_aic_post <- numeric(0)
  history_eps_rho <- numeric(0)
  history_eps_par <- numeric(0)
  history_eps_bg <- numeric(0)

  current_params_total <- NULL
  current_se_total <- NULL
  current_optimisation <- NULL
  current_objective <- NULL
  current_predictor <- NULL
  current_trigger <- NULL
  current_intensity <- NULL
  current_rho <- background_weights
  converged <- FALSE
  stop_reason <- NULL

  for (iteration in seq_len(ndeclust)) {
    if (iteration > min.declust.iter &&
        length(history_eps_rho) > 0L &&
        tail(history_eps_rho, 1L) <= tol.rho &&
        tail(history_eps_par, 1L) <= tol.par) {
      converged <- TRUE
      stop_reason <- "Convergence thresholds reached."
      break
    }

    if (isTRUE(verbose)) {
      message("Start ETAS maximum-likelihood step ", iteration, " of ", ndeclust)
    }

    optimisation <- .pp_etas_optimize(
      state = state,
      hessian = TRUE,
      iterlim = iterlim,
      trace = trace
    )

    n_etas <- sum(params.ind)

    unpacked <- .pp_etas_unpack_parameters(
      candidate = optimisation$parameters,
      params.fix = state$params.fix,
      params.ind = params.ind,
      params.lim = params.lim,
      n_covariates = n_covariates
    )

    params_current <- unpacked$params
    betacov_current <- unpacked$betacov

    objective_ml <- .pp_etas_objective(
      candidate = optimisation$parameters,
      state = state,
      trace = FALSE
    )

    lambda_ml <- attr(objective_ml, "lambda.vec")

    if (is.null(lambda_ml) || any(!is.finite(lambda_ml)) || any(lambda_ml <= 0)) {
      stop("The ETAS likelihood returned invalid total intensities.", call. = FALSE)
    }

    rho_raw <- params_current[["mu"]] * state$back.dens / lambda_ml

    if (any(rho_raw < -1e-8 | rho_raw > 1 + 1e-8, na.rm = TRUE)) {
      warning(
        "Some posterior background probabilities were outside [0, 1] and were truncated.",
        call. = FALSE
      )
    }

    rho_current <- pmin(pmax(rho_raw, 0), 1)

    post_background <- state
    post_background_state <- background_state

    if (isTRUE(declustering) && !isTRUE(is.backconstant) && !isTRUE(onlytime)) {
      for (weight_step in seq_len(n.iterweight)) {
        if (isTRUE(verbose)) {
          message("Parametric background update ", weight_step, " of ", n.iterweight)
        }

        post_background_state <- fit_background(
          weights = rho_current,
          force_constant = FALSE
        )

        post_background$back.dens <- post_background_state$density_events
        post_background$back.integral <- post_background_state$integral
      }
    }

    objective_post <- .pp_etas_objective(
      candidate = optimisation$parameters,
      state = post_background,
      trace = FALSE
    )

    lambda_post <- attr(objective_post, "lambda.vec")
    trigger_post <- attr(objective_post, "etas.vec")
    predictor_post <- attr(objective_post, "predictor")

    se <- .pp_etas_standard_errors(
      optimisation = optimisation,
      params.ind = params.ind,
      n_covariates = n_covariates
    )

    se_etas <- rep(NA_real_, 7L)
    names(se_etas) <- namespar
    se_etas[params.ind] <- se[seq_len(n_etas)]
    se_etas[!params.ind] <- 0

    se_covariates <- if (n_covariates > 0L) {
      se[(n_etas + 1L):(n_etas + n_covariates)]
    } else {
      numeric(0)
    }

    params_total <- c(params_current, betacov_current)
    se_total <- c(se_etas, se_covariates)

    if (n_covariates > 0L) {
      names(betacov_current) <- colnames(cov.matrix)
      names(se_covariates) <- colnames(cov.matrix)
      params_total <- c(params_current, betacov_current)
      se_total <- c(se_etas, se_covariates)
    }

    edf_background <- post_background_state$edf
    aic_ml <- 2 * optimisation$value + 2 * (n_etas + n_covariates + edf_background)
    aic_post <- 2 * as.numeric(objective_post) + 2 * (n_etas + n_covariates + edf_background)

    eps_rho <- if (length(history_rho) == 0L) {
      Inf
    } else {
      mean(abs(rho_current - history_rho[[length(history_rho)]]), na.rm = TRUE)
    }

    eps_par <- if (length(history_parameters) == 0L) {
      Inf
    } else {
      previous <- history_parameters[[length(history_parameters)]]
      max(
        abs(params_total - previous) / (abs(previous) + 1e-8),
        na.rm = TRUE
      )
    }

    eps_bg <- if (length(history_background_density) == 0L) {
      Inf
    } else {
      previous <- history_background_density[[length(history_background_density)]]
      mean(
        abs(post_background$back.dens - previous) / (abs(previous) + 1e-8),
        na.rm = TRUE
      )
    }

    if (isTRUE(stop.on.increasing.AIC) &&
        length(history_aic_ml) > 0L &&
        aic_ml > tail(history_aic_ml, 1L)) {
      stop_reason <- "Working AIC increased."
      break
    }

    history_parameters[[length(history_parameters) + 1L]] <- params_total
    history_se[[length(history_se) + 1L]] <- se_total
    history_rho[[length(history_rho) + 1L]] <- rho_current
    history_background_density[[length(history_background_density) + 1L]] <- post_background$back.dens
    history_aic_ml <- c(history_aic_ml, aic_ml)
    history_aic_post <- c(history_aic_post, aic_post)
    history_eps_rho <- c(history_eps_rho, eps_rho)
    history_eps_par <- c(history_eps_par, eps_par)
    history_eps_bg <- c(history_eps_bg, eps_bg)

    current_params_total <- params_total
    current_se_total <- se_total
    current_optimisation <- optimisation$result
    current_objective <- objective_post
    current_predictor <- predictor_post
    current_trigger <- trigger_post
    current_intensity <- lambda_post
    current_rho <- rho_current

    state <- post_background
    state$params.fix <- params_current
    state$betacov <- betacov_current
    background_state <- post_background_state

    if (isTRUE(verbose.convergence)) {
      message(
        "Iteration ", iteration,
        ": eps.rho = ", signif(eps_rho, 4),
        "; eps.par = ", signif(eps_par, 4),
        "; eps.bg = ", signif(eps_bg, 4),
        "; mean rho = ", signif(mean(rho_current), 4),
        "; sum rho = ", signif(sum(rho_current), 4)
      )
    }

    if (iteration >= min.declust.iter &&
        eps_rho <= tol.rho &&
        eps_par <= tol.par) {
      converged <- TRUE
      stop_reason <- "Convergence thresholds reached."
      break
    }
  }

  if (is.null(current_params_total)) {
    stop("ETAS fitting stopped before completing one estimation iteration.", call. = FALSE)
  }

  parameter_history <- do.call(rbind, history_parameters)
  se_history <- do.call(rbind, history_se)
  rho_history <- do.call(rbind, history_rho)

  rownames(parameter_history) <- seq_len(nrow(parameter_history))
  rownames(se_history) <- seq_len(nrow(se_history))
  rownames(rho_history) <- seq_len(nrow(rho_history))

  if (is.null(stop_reason)) {
    stop_reason <- if (converged) "Convergence thresholds reached." else "Maximum declustering iterations reached."
  }

  if (any(is.na(current_se_total))) {
    warning(
      paste0(
        "Standard errors of the triggering parameters could not be computed: ",
        "the Hessian of the log-likelihood is singular or not positive definite ",
        "(estimated k0 = ", format(signif(current_params_total[["k0"]], 3)), "). ",
        "This usually means that the catalogue contains too little triggering ",
        "to identify the triggering parameters."
      ),
      call. = FALSE
    )
  }

  out <- list(
    call = match.call(),
    description = description,

    catalog = catalog,
    cat = cat,
    selected_catalog = prepared,

    parameters = current_params_total,
    params = current_params_total[namespar],
    betacov = current_params_total[setdiff(names(current_params_total), namespar)],
    standard_errors = current_se_total,
    sqm = current_se_total,

    params.ind = params.ind,
    params.lim = params.lim,
    formula1 = trigger_model$formula,
    cov.matrix = cov.matrix,
    offset = trigger_model$offset,
    predictor = current_predictor,

    background_fit = background_state$fit,
    model.bg = background_state$fit,
    background_density = state$back.dens,
    back.dens = state$back.dens,
    background_integral = state$back.integral,
    back.integral = state$back.integral,
    background_normalisation = background_state$kappa,

    triggering_intensity = current_trigger,
    total_intensity = current_intensity,
    lambda = current_intensity,
    background_probability = current_rho,
    rho.weights = current_rho,
    triggered_probability = 1 - current_rho,

    negative_loglik = as.numeric(current_objective),
    loglik = -as.numeric(current_objective),
    integral = attr(current_objective, "integraltot"),
    trigger_integral = attr(current_objective, "trigger.integral"),

    aic_iter = history_aic_ml,
    aic_post_background_iter = history_aic_post,
    AIC.iter = history_aic_ml,
    AIC.iter2 = history_aic_post,

    params.iter = parameter_history,
    sqm.iter = se_history,
    rho.weights.iter = rho_history,

    history = list(
      parameters = parameter_history,
      standard_errors = se_history,
      background_probability = rho_history,
      background_density = history_background_density,
      aic_ml = history_aic_ml,
      aic_post_background = history_aic_post,
      eps_rho = history_eps_rho,
      eps_par = history_eps_par,
      eps_bg = history_eps_bg
    ),

    optimisation = current_optimisation,
    convergence = list(
      converged = converged,
      reason = stop_reason,
      iterations = length(history_parameters),
      epsmax = epsmax,
      tol_rho = tol.rho,
      tol_par = tol.par,
      final_eps_rho = tail(history_eps_rho, 1L),
      final_eps_par = tail(history_eps_par, 1L),
      final_eps_bg = tail(history_eps_bg, 1L)
    ),

    settings = list(
      magn.threshold = magn.threshold,
      magn.threshold.back = magn.threshold.back,
      tmax = prepared$tmax_original,
      long.range = prepared$long.range,
      lat.range = prepared$lat.range,
      longlat.to.km = catalog$longlat.to.km,
      sectoday = catalog$sectoday,
      onlytime = onlytime,
      is.backconstant = is.backconstant,
      declustering = declustering,
      ndeclust = ndeclust,
      n.iterweight = n.iterweight,
      hdef = hdef,
      hvarx = hvarx,
      hvary = hvary,
      m1 = m1,
      back.smooth = back.smooth,
      formula.bg = formula.bg,
      q.bg = q.bg,
      dummy_layout.bg = dummy_layout.bg,
      nc.bg = nc.bg,
      seed.bg = seed.bg,
      mark_n.bg = mark_n.bg,
      usenlm = usenlm,
      method = method,
      compsqm = compsqm,
      iterlim = iterlim,
      ntheta = ntheta
    ),

    time.start = fit_start_time,
    time.end = Sys.time()
  )

  class(out) <- "pp_etas_fit"
  out
}


# Fitted-object methods -------------------------------------------------------

#' Extract eventwise ETAS declustering probabilities
#'
#' Returns the fitted background and triggering components for each selected
#' earthquake together with its posterior probabilities of being a background
#' event or a triggered event.
#'
#' @details
#' For event \eqn{i}, the posterior background probability is
#'
#' \deqn{\rho_i = \frac{\mu b(\mathbf{s}_i)}
#' {\mu b(\mathbf{s}_i) + \lambda_{\mathrm{trig},i}},}
#'
#' and the posterior triggered probability is \eqn{1-\rho_i}. These are soft
#' classifications: the function does not randomly thin the catalog and does
#' not assign a deterministic binary label.
#'
#' @param fit An object of class `"pp_etas_fit"` returned by [fit_etas_pp()].
#' @param sort_by Character string controlling row order. `"time"` orders by
#'   event time and source identifier; `"background_probability"` orders from
#'   highest to lowest posterior background probability.
#'
#' @return A data frame with source identifier, original catalog coordinates,
#'   depth, magnitude, fitted background density, triggering and total
#'   intensities, and posterior background and triggered probabilities.
#'
#' @examples
#' \donttest{
#' sim <- simulate_etas_pp(
#'   params = c(mu = 0.5, k0 = 1.5e-5, c = 0.05, p = 1.2,
#'              gamma = 0, d = 0.001, q = 2),
#'   m0 = 2.5, b = 1, tmin = 0, t.lag = 500,
#'   long.range = c(0, 1), lat.range = c(0, 1), longlat.to.km = FALSE,
#'   n_support_obs = 50, mult_support = 1,
#'   trig_lp_coefs = c(m_rel = 0.3),
#'   seed = 19
#' )
#' fit <- fit_etas_pp(sim$catalog)
#'
#' event_probabilities <- decluster_etas_pp(fit)
#' head(event_probabilities)
#'
#' most_likely_background <- decluster_etas_pp(
#'   fit,
#'   sort_by = "background_probability"
#' )
#' head(most_likely_background)
#' }
#'
#' @references
#' Zhuang, J., Ogata, Y. and Vere-Jones, D. (2002). Stochastic declustering of
#' space-time earthquake occurrences. *Journal of the American Statistical
#' Association*, 97(458), 369--380.
#'
#' @seealso [fit_etas_pp()] and [plot.pp_etas_fit()].
#'
#' @family ETAS modelling
#' @export
decluster_etas_pp <- function(
    fit,
    sort_by = c("time", "background_probability")
) {
  if (!inherits(fit, "pp_etas_fit")) {
    stop("`fit` must inherit from class 'pp_etas_fit'.", call. = FALSE)
  }

  sort_by <- match.arg(sort_by)

  out <- data.frame(
    source_id = fit$cat$.pp_source_id,
    time = fit$cat$time,
    long = fit$cat$long,
    lat = fit$cat$lat,
    z = fit$cat$z,
    magnitude = fit$cat$magn1,
    background_density = fit$background_density,
    triggering_intensity = fit$triggering_intensity,
    total_intensity = fit$total_intensity,
    background_probability = fit$background_probability,
    triggered_probability = fit$triggered_probability,
    check.names = FALSE
  )

  if (identical(sort_by, "background_probability")) {
    out <- out[order(out$background_probability, decreasing = TRUE), , drop = FALSE]
  } else {
    out <- out[order(out$time, out$source_id), , drop = FALSE]
  }

  rownames(out) <- NULL
  out
}


#' Methods for fitted ETAS models
#'
#' Print, summarize, and plot an object returned by [fit_etas_pp()].
#'
#' @description
#' The print method reports the number of fitted events, alternating-fit
#' convergence, log likelihood, and parameter estimates. The summary method
#' additionally returns standard errors, the integrated intensity, summaries of
#' the posterior probabilities, and retained settings. The plot method displays
#' posterior background probabilities over time or space.
#'
#' @param x An object of class `"pp_etas_fit"` for `print()` and `plot()`, or an
#'   object of class `"summary_pp_etas_fit"` when printing a summary.
#' @param object An object of class `"pp_etas_fit"`.
#' @param digits Positive integer controlling the number of significant digits
#'   printed for likelihoods and parameter estimates.
#' @param what Character string selecting the graphical display:
#'   `"declustering"` and `"space"` both show the spatial event pattern colored
#'   by posterior background probability, whereas `"time"` plots posterior
#'   background probability against event time.
#' @param main Optional character string used as the plot title. A
#'   display-specific title is generated when `NULL`.
#' @param ... Additional graphical parameters passed to [graphics::plot()] by
#'   the plot method. Additional arguments supplied to print or summary methods
#'   are currently unused.
#'
#' @return `print.pp_etas_fit()`, `print.summary_pp_etas_fit()`, and
#'   `plot.pp_etas_fit()` invisibly return `x`. `summary.pp_etas_fit()` returns
#'   an object of class `"summary_pp_etas_fit"` containing the event count,
#'   log likelihood, integrated intensity, parameter table, summaries of
#'   background and triggered probabilities, convergence information, and
#'   settings.
#'
#' @examples
#' \donttest{
#' sim <- simulate_etas_pp(
#'   params = c(mu = 0.5, k0 = 1.5e-5, c = 0.05, p = 1.2,
#'              gamma = 0, d = 0.001, q = 2),
#'   m0 = 2.5, b = 1, tmin = 0, t.lag = 500,
#'   long.range = c(0, 1), lat.range = c(0, 1), longlat.to.km = FALSE,
#'   n_support_obs = 50, mult_support = 1,
#'   trig_lp_coefs = c(m_rel = 0.3),
#'   seed = 19
#' )
#' fit <- fit_etas_pp(sim$catalog)
#'
#' print(fit)
#' summary(fit)
#' plot(fit, what = "declustering")
#' plot(fit, what = "time")
#' plot(fit, what = "space")
#' }
#'
#' @seealso [fit_etas_pp()] and [decluster_etas_pp()].
#'
#' @family ETAS modelling
#' @rdname pp_etas_fit-methods
#' @method print pp_etas_fit
#' @export
print.pp_etas_fit <- function(x, digits = 5, ...) {
  if (!inherits(x, "pp_etas_fit")) {
    stop("`x` must inherit from class 'pp_etas_fit'.", call. = FALSE)
  }

  cat("<pp_etas_fit>\n")
  cat("Parametric-background ETAS model\n\n")
  cat("Events: ", nrow(x$cat), "\n", sep = "")
  cat("Declustering iterations: ", x$convergence$iterations, "\n", sep = "")
  cat("Converged: ", x$convergence$converged, "\n", sep = "")
  cat("Stop reason: ", x$convergence$reason, "\n", sep = "")
  cat("Log-likelihood: ", formatC(x$loglik, digits = digits, format = "fg"), "\n\n", sep = "")
  cat("Parameters:\n")
  print(round(x$parameters, digits = digits))
  invisible(x)
}


#' @rdname pp_etas_fit-methods
#' @method summary pp_etas_fit
#' @export
summary.pp_etas_fit <- function(object, ...) {
  if (!inherits(object, "pp_etas_fit")) {
    stop("`object` must inherit from class 'pp_etas_fit'.", call. = FALSE)
  }

  parameter_table <- data.frame(
    estimate = object$parameters,
    std_error = object$standard_errors[names(object$parameters)],
    row.names = names(object$parameters),
    check.names = FALSE
  )

  out <- list(
    n_events = nrow(object$cat),
    loglik = object$loglik,
    integral = object$integral,
    parameters = parameter_table,
    background_probability = summary(object$background_probability),
    triggered_probability = summary(object$triggered_probability),
    convergence = object$convergence,
    settings = object$settings
  )

  class(out) <- "summary_pp_etas_fit"
  out
}


#' @rdname pp_etas_fit-methods
#' @method print summary_pp_etas_fit
#' @export
print.summary_pp_etas_fit <- function(x, digits = 5, ...) {
  cat("Parametric-background ETAS summary\n\n")
  cat("Events: ", x$n_events, "\n", sep = "")
  cat("Log-likelihood: ", formatC(x$loglik, digits = digits, format = "fg"), "\n", sep = "")
  cat("Integrated intensity: ", formatC(x$integral, digits = digits, format = "fg"), "\n", sep = "")
  cat("Converged: ", x$convergence$converged, "\n", sep = "")
  cat("Stop reason: ", x$convergence$reason, "\n\n", sep = "")

  print(x$parameters, digits = digits)
  invisible(x)
}


#' @rdname pp_etas_fit-methods
#' @method plot pp_etas_fit
#' @export
plot.pp_etas_fit <- function(
    x,
    what = c("declustering", "space", "time"),
    main = NULL,
    ...
) {
  if (!inherits(x, "pp_etas_fit")) {
    stop("`x` must inherit from class 'pp_etas_fit'.", call. = FALSE)
  }

  what <- match.arg(what)

  if (identical(what, "time")) {
    graphics::plot(
      x$cat$time,
      x$background_probability,
      xlab = "time",
      ylab = "posterior background probability",
      main = if (is.null(main)) "ETAS declustering over time" else main,
      ylim = c(0, 1),
      pch = 19,
      ...
    )
    graphics::abline(h = c(0, 0.5, 1), lty = c(1, 2, 1))
    return(invisible(x))
  }

  if (identical(what, "space")) {
    probability <- x$background_probability
    pal <- grDevices::hcl.colors(101L, palette = "YlOrRd", rev = TRUE)
    index <- pmin(101L, pmax(1L, floor(probability * 100) + 1L))

    graphics::plot(
      x$cat$long,
      x$cat$lat,
      xlab = "longitude / x",
      ylab = "latitude / y",
      main = if (is.null(main)) "Posterior background probability in space" else main,
      col = pal[index],
      pch = 19,
      ...
    )

    graphics::legend(
      "topright",
      legend = c("mostly triggered", "mixed", "mostly background"),
      col = pal[c(1L, 51L, 101L)],
      pch = 19,
      bty = "o",
      bg = grDevices::adjustcolor("white", alpha.f = 0.85),
      box.col = "grey80"
    )

    return(invisible(x))
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::par(mfrow = c(1, 2))

  graphics::plot(
    x$cat$time,
    x$background_probability,
    xlab = "time",
    ylab = "background probability",
    ylim = c(0, 1),
    main = "Declustering over time",
    pch = 19
  )
  graphics::abline(h = 0.5, lty = 2)

  pal <- grDevices::hcl.colors(101L, palette = "YlOrRd", rev = TRUE)
  index <- pmin(101L, pmax(1L, floor(x$background_probability * 100) + 1L))

  graphics::plot(
    x$cat$long,
    x$cat$lat,
    xlab = "longitude / x",
    ylab = "latitude / y",
    main = "Declustering in space",
    col = pal[index],
    pch = 19
  )

  graphics::legend(
    "topright",
    legend = c("triggered", "background"),
    col = pal[c(1L, 101L)],
    pch = 19,
    bty = "o",
    bg = grDevices::adjustcolor("white", alpha.f = 0.85),
    box.col = "grey80"
  )

  invisible(x)
}
