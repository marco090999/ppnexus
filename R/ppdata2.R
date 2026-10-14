

# Internal helpers ---------------------------------------------------------

.pp_resolve_columns <- function(data, columns, arg_name) {
  nms <- names(data)

  if (is.character(columns)) {
    if (length(columns) == 0L) {
      stop(sprintf("`%s` cannot be empty.", arg_name), call. = FALSE)
    }

    missing_cols <- setdiff(columns, nms)
    if (length(missing_cols) > 0L) {
      stop(
        sprintf(
          "`%s` contains unknown column(s): %s.",
          arg_name,
          paste(missing_cols, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    return(match(columns, nms))
  }

  if (is.numeric(columns)) {
    if (length(columns) == 0L ||
        any(!is.finite(columns)) ||
        any(columns != as.integer(columns))) {
      stop(
        sprintf("`%s` must contain valid integer column indices.", arg_name),
        call. = FALSE
      )
    }

    columns <- as.integer(columns)

    if (any(columns < 1L | columns > ncol(data))) {
      stop(
        sprintf("`%s` contains column indices outside the data.", arg_name),
        call. = FALSE
      )
    }

    return(columns)
  }

  stop(
    sprintf("`%s` must be a character vector of names or numeric column indices.",
            arg_name),
    call. = FALSE
  )
}


.pp_align_axis_argument <- function(x, axis_names, arg_name) {
  if (is.null(x)) {
    return(NULL)
  }

  if (!is.character(x)) {
    stop(sprintf("`%s` must be a character vector.", arg_name), call. = FALSE)
  }

  if (!is.null(names(x))) {
    if (anyDuplicated(names(x))) {
      stop(
        sprintf("`%s` cannot contain duplicated axis names.", arg_name),
        call. = FALSE
      )
    }

    missing_axes <- setdiff(axis_names, names(x))
    extra_axes <- setdiff(names(x), axis_names)

    if (length(missing_axes) > 0L || length(extra_axes) > 0L) {
      stop(
        sprintf(
          "`%s` must be named exactly with these axes: %s.",
          arg_name,
          paste(axis_names, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    x <- x[axis_names]
  } else if (length(x) != length(axis_names)) {
    stop(
      sprintf(
        "`%s` must have length %d or be named using the coordinate names.",
        arg_name,
        length(axis_names)
      ),
      call. = FALSE
    )
  }

  unname(x)
}


.pp_convert_time_coordinate <- function(x, time_unit, time_origin, axis_name) {
  original_class <- class(x)[1]

  if (inherits(x, "Date")) {
    if (is.null(time_origin)) {
      if (length(x) == 0L) {
        stop(
          sprintf(
            "An empty Date time coordinate `%s` requires `time_origin`.",
            axis_name
          ),
          call. = FALSE
        )
      }

      time_origin <- min(x, na.rm = TRUE)
    }

    if (!inherits(time_origin, "Date")) {
      stop(
        sprintf(
          "`time_origin` must be a Date object for the time coordinate `%s`.",
          axis_name
        ),
        call. = FALSE
      )
    }

    values <- as.numeric(difftime(x, time_origin, units = time_unit))
  } else if (inherits(x, c("POSIXct", "POSIXt"))) {
    if (is.null(time_origin)) {
      if (length(x) == 0L) {
        stop(
          sprintf(
            "An empty POSIXct time coordinate `%s` requires `time_origin`.",
            axis_name
          ),
          call. = FALSE
        )
      }

      time_origin <- min(x, na.rm = TRUE)
    }

    if (!inherits(time_origin, c("POSIXct", "POSIXt"))) {
      stop(
        sprintf(
          "`time_origin` must be POSIXct/POSIXt for the time coordinate `%s`.",
          axis_name
        ),
        call. = FALSE
      )
    }

    values <- as.numeric(difftime(x, time_origin, units = time_unit))
  } else if (is.numeric(x)) {
    if (is.null(time_origin)) {
      time_origin <- 0
    }

    if (is.numeric(time_origin)) {
      if (length(time_origin) != 1L || !is.finite(time_origin)) {
        stop(
          sprintf(
            "`time_origin` must be a single finite numeric value for `%s`.",
            axis_name
          ),
          call. = FALSE
        )
      }

      values <- as.numeric(x) - time_origin
    } else if (inherits(time_origin, "Date") ||
               inherits(time_origin, c("POSIXct", "POSIXt"))) {
      # Numeric input is interpreted as already expressed relative to the
      # supplied Date/POSIXct origin.
      values <- as.numeric(x)
    } else {
      stop(
        sprintf(
          "`time_origin` must be numeric, Date, or POSIXct/POSIXt for `%s`.",
          axis_name
        ),
        call. = FALSE
      )
    }
  } else {
    stop(
      sprintf(
        "The time coordinate `%s` must be numeric, Date, or POSIXct/POSIXt.",
        axis_name
      ),
      call. = FALSE
    )
  }

  if (any(!is.finite(values))) {
    stop(
      sprintf("The time coordinate `%s` contains missing or non-finite values.",
              axis_name),
      call. = FALSE
    )
  }

  list(
    values = values,
    origin = time_origin,
    original_class = original_class,
    unit = time_unit
  )
}


.pp_same_time_origin <- function(x, y) {
  if (is.null(x) || is.null(y)) {
    return(is.null(x) && is.null(y))
  }

  if (inherits(x, "Date") && inherits(y, "Date")) {
    return(identical(as.numeric(x), as.numeric(y)))
  }

  if (inherits(x, c("POSIXct", "POSIXt")) &&
      inherits(y, c("POSIXct", "POSIXt"))) {
    return(isTRUE(all.equal(
      as.numeric(x),
      as.numeric(y),
      tolerance = .Machine$double.eps^0.5
    )))
  }

  if (is.numeric(x) && is.numeric(y) &&
      length(x) == 1L && length(y) == 1L) {
    return(isTRUE(all.equal(
      as.numeric(x),
      as.numeric(y),
      tolerance = .Machine$double.eps^0.5
    )))
  }

  FALSE
}


.pp_validate_time_unit <- function(time_unit, arg_name = "time_unit") {
  allowed_units <- c("days", "hours", "mins", "secs", "weeks")

  if (!is.character(time_unit) ||
      length(time_unit) != 1L ||
      is.na(time_unit) ||
      !time_unit %in% allowed_units) {
    stop(
      sprintf(
        "`%s` must be one of: %s.",
        arg_name,
        paste(allowed_units, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  time_unit
}


#' Create a multidimensional point-pattern object
#'
#' Constructs a point-pattern object from a data frame or matrix. A pattern may
#' be purely spatial, with any positive number of spatial dimensions, or
#' spatio-temporal, with one temporal dimension and at least one spatial
#' dimension. One or more columns may be retained as event marks.
#'
#' @details
#' Coordinate columns are converted to a numeric matrix and stored in the order
#' supplied by `coords`, except when `window` is provided. In that case, the
#' order, roles, units, temporal origin, and coordinate reference system stored
#' in the window are authoritative. All observed events must lie inside the
#' supplied window.
#'
#' A temporal coordinate may be numeric, `Date`, or `POSIXct`/`POSIXt`. It is
#' represented internally as a numeric displacement from `time_origin`, in the
#' unit specified by `time_unit`. Numeric temporal values are shifted by a
#' numeric origin when one is supplied. If a `Date` or `POSIXct` origin is
#' supplied with numeric values, those values are treated as already expressed
#' relative to that origin.
#'
#' Marks are stored in a data frame without changing their original classes.
#' Columns that are neither coordinates nor selected marks are recorded in the
#' object's metadata but are not otherwise used.
#'
#' The returned object inherits from `"pp"` and from one of the following
#' subclasses:
#'
#' - `"pp_spatial"`: unmarked spatial pattern;
#' - `"pp_spatial_marked"`: marked spatial pattern;
#' - `"pp_st"`: unmarked spatio-temporal pattern;
#' - `"pp_st_marked"`: marked spatio-temporal pattern.
#'
#' @param data A data frame or matrix containing event coordinates and,
#'   optionally, event marks. Coordinate columns must be finite and numeric,
#'   except for a temporal coordinate, which may also be `Date` or
#'   `POSIXct`/`POSIXt`.
#' @param coords Optional coordinate-column specification. Supply a character
#'   vector of column names or a numeric vector of column positions. If omitted
#'   and `window` is supplied, the axes stored in `window` are used. Cannot be
#'   supplied together with `d`.
#' @param d Optional positive integer indicating that the first `d` columns of
#'   `data` are coordinates. Cannot be supplied together with `coords`.
#' @param axis_roles Optional character vector describing the role of each
#'   coordinate. Allowed values are `"space"` and `"time"`. It may be unnamed,
#'   in coordinate order, or named with the coordinate names. If omitted, roles
#'   are inherited from `window` when available; otherwise all coordinates are
#'   spatial. At most one temporal axis is allowed, and a spatio-temporal
#'   pattern must contain at least one spatial axis.
#' @param marks Optional mark-column specification. Use `NULL` for an unmarked
#'   pattern, `"auto"` to retain every non-coordinate column as a mark, or a
#'   character/numeric vector identifying selected mark columns. Coordinate and
#'   mark columns must be distinct.
#' @param window Optional observation window created by [pp_window()]. Its axes
#'   must be compatible with the selected coordinates, and all observed points
#'   must lie within its bounds. If omitted, no statistical observation window
#'   is inferred; the observed bounding box is used only for descriptive plots
#'   and summaries.
#' @param axis_units Optional character vector of coordinate units. It may be
#'   unnamed, in coordinate order, or named with the coordinate names. If
#'   omitted, units are inherited from `window` when available.
#' @param time_unit Optional unit used to represent the temporal coordinate
#'   internally. Allowed values are `"days"`, `"hours"`, `"mins"`, `"secs"`,
#'   and `"weeks"`. If omitted, it is inherited from `window`, then from the
#'   temporal entry of `axis_units`, and otherwise defaults to `"days"`.
#' @param time_origin Optional temporal origin. For numeric time it must be a
#'   single finite numeric value, `Date`, or `POSIXct`/`POSIXt`; a numeric origin
#'   is subtracted from the observed values. For `Date` and `POSIXct` input it
#'   must have the corresponding class. If omitted, it is inherited from
#'   `window` when available; otherwise the minimum observed date/time is used
#'   for date-time input and zero is used for numeric input.
#' @param crs Optional non-empty character string describing the coordinate
#'   reference system of the spatial coordinates. It is inherited from
#'   `window` when available and must agree with it when supplied explicitly.
#'
#' @return An object inheriting from `"pp"`. Its principal components are:
#'
#' - `coords`: numeric coordinate matrix;
#' - `marks`: data frame of marks, or `NULL`;
#' - `window`: the supplied [pp_window()] object, or `NULL`;
#' - `axis_info`: data frame containing axis names, roles, and units;
#' - `time_info`: temporal-axis metadata, or `NULL` for a purely spatial
#'   pattern;
#' - `crs`: coordinate reference system, or `NULL`;
#' - `metadata`: original column names and columns not retained as coordinates
#'   or marks.
#'
#' @examples
#' # A marked two-dimensional spatial pattern with an observation window
#' w <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1),
#'   axis_units = c(x = "km", y = "km")
#' )
#'
#' set.seed(1)
#' dat <- data.frame(
#'   x = runif(40),
#'   y = runif(40),
#'   magnitude = runif(40, 1, 4),
#'   type = sample(c("A", "B"), 40, replace = TRUE)
#' )
#'
#' x <- pp_pattern(
#'   data = dat,
#'   coords = c("x", "y"),
#'   marks = c("magnitude", "type"),
#'   window = w
#' )
#' x
#' summary(x)
#'
#' # An unmarked three-dimensional spatial pattern
#' x3 <- pp_pattern(
#'   data = data.frame(
#'     x = runif(20), y = runif(20), z = runif(20)
#'   ),
#'   coords = c("x", "y", "z")
#' )
#'
#' # A marked spatio-temporal pattern with calendar dates
#' origin <- as.Date("2020-01-01")
#' wst <- pp_window(
#'   lower = c(x = 0, y = 0, time = 0),
#'   upper = c(x = 1, y = 1, time = 30),
#'   axis_roles = c(x = "space", y = "space", time = "time"),
#'   axis_units = c(x = "km", y = "km", time = "days"),
#'   time_origin = origin
#' )
#'
#' xst <- pp_pattern(
#'   data = data.frame(
#'     x = runif(20),
#'     y = runif(20),
#'     time = origin + sample(0:30, 20, replace = TRUE),
#'     magnitude = runif(20, 1, 4)
#'   ),
#'   coords = c("x", "y", "time"),
#'   marks = "magnitude",
#'   window = wst
#' )
#'
#' @seealso
#' [pp_window()] for observation windows, [pp_covariate()] for covariate
#' supports, [pp_grid()] for regular evaluation grids, [simulate_pp()] for
#' simulation, and [fit_pp()] for Poisson point-process fitting.
#'
#' @family point-pattern data
#' @export
pp_pattern <- function(
    data,
    coords = NULL,
    d = NULL,
    axis_roles = NULL,
    marks = NULL,
    window = NULL,
    axis_units = NULL,
    time_unit = NULL,
    time_origin = NULL,
    crs = NULL
) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data frame or matrix.", call. = FALSE)
  }

  data <- as.data.frame(data, stringsAsFactors = FALSE)

  if (ncol(data) == 0L) {
    stop("`data` must contain at least one column.", call. = FALSE)
  }

  if (is.null(names(data)) || any(!nzchar(names(data)))) {
    names(data) <- paste0("V", seq_len(ncol(data)))
  }

  if (anyDuplicated(names(data))) {
    stop("`data` must have unique column names.", call. = FALSE)
  }

  if (!is.null(window) && !inherits(window, "pp_window")) {
    stop(
      "`window` must be NULL or an object of class 'pp_window'.",
      call. = FALSE
    )
  }

  if (!is.null(coords) && !is.null(d)) {
    stop("Supply either `coords` or `d`, not both.", call. = FALSE)
  }

  if (is.null(coords) && is.null(d) && is.null(window)) {
    stop(
      "Specify coordinate columns through `coords` or provide their number using `d`.",
      call. = FALSE
    )
  }

  if (is.null(coords) && is.null(d) && !is.null(window)) {
    coords <- window$axis_info$name
  }

  if (!is.null(d)) {
    if (!is.numeric(d) ||
        length(d) != 1L ||
        !is.finite(d) ||
        d != as.integer(d) ||
        d < 1L ||
        d > ncol(data)) {
      stop(
        "`d` must be a positive integer not larger than the number of columns.",
        call. = FALSE
      )
    }

    coord_idx <- seq_len(as.integer(d))
  } else {
    coord_idx <- .pp_resolve_columns(data, coords, "coords")
  }

  if (anyDuplicated(coord_idx)) {
    stop("`coords` cannot contain duplicated columns.", call. = FALSE)
  }

  coord_names <- names(data)[coord_idx]

  if (!is.null(window)) {
    window_axes <- window$axis_info$name

    if (!setequal(coord_names, window_axes)) {
      stop(
        paste0(
          "The selected coordinate axes must match the axes of `window`. ",
          "Pattern axes: ", paste(coord_names, collapse = ", "),
          ". Window axes: ", paste(window_axes, collapse = ", "), "."
        ),
        call. = FALSE
      )
    }

    # The stored point-pattern coordinates follow the order defined by window.
    coord_idx <- coord_idx[match(window_axes, coord_names)]
    coord_names <- names(data)[coord_idx]
  }

  d <- length(coord_idx)

  axis_roles <- .pp_align_axis_argument(
    x = axis_roles,
    axis_names = coord_names,
    arg_name = "axis_roles"
  )

  if (is.null(axis_roles) && !is.null(window)) {
    axis_roles <- window$axis_info$role[
      match(coord_names, window$axis_info$name)
    ]
  }

  if (is.null(axis_roles)) {
    axis_roles <- rep("space", d)
  }

  axis_roles <- tolower(trimws(axis_roles))

  if (!all(axis_roles %in% c("space", "time"))) {
    stop(
      "`axis_roles` may contain only 'space' and 'time'.",
      call. = FALSE
    )
  }

  if (!is.null(window)) {
    window_roles <- window$axis_info$role[
      match(coord_names, window$axis_info$name)
    ]

    if (!identical(axis_roles, window_roles)) {
      stop(
        "`axis_roles` is not consistent with the roles stored in `window`.",
        call. = FALSE
      )
    }
  }

  n_time <- sum(axis_roles == "time")

  if (n_time > 1L) {
    stop(
      "A spatio-temporal point pattern currently supports exactly one time axis.",
      call. = FALSE
    )
  }

  if (n_time == 1L && sum(axis_roles == "space") == 0L) {
    stop(
      "A spatio-temporal point pattern must contain at least one spatial axis.",
      call. = FALSE
    )
  }

  if (n_time == 0L &&
      (!is.null(time_unit) || !is.null(time_origin))) {
    stop(
      "`time_unit` and `time_origin` can be supplied only when one axis has role 'time'.",
      call. = FALSE
    )
  }

  axis_units <- .pp_align_axis_argument(
    x = axis_units,
    axis_names = coord_names,
    arg_name = "axis_units"
  )

  if (is.null(axis_units) && !is.null(window)) {
    axis_units <- window$axis_info$unit[
      match(coord_names, window$axis_info$name)
    ]
  }

  if (is.null(axis_units)) {
    axis_units <- rep(NA_character_, d)
  }

  if (n_time == 1L) {
    time_idx <- which(axis_roles == "time")
    window_time_info <- if (!is.null(window)) window$time_info else NULL

    if (is.null(time_unit) && !is.null(window_time_info)) {
      time_unit <- window_time_info$unit
    }

    if (is.null(time_unit) &&
        !is.na(axis_units[time_idx]) &&
        nzchar(axis_units[time_idx])) {
      time_unit <- axis_units[time_idx]
    }

    if (is.null(time_unit)) {
      time_unit <- "days"
    }

    time_unit <- .pp_validate_time_unit(time_unit)

    if (!is.na(axis_units[time_idx]) &&
        nzchar(axis_units[time_idx]) &&
        !identical(axis_units[time_idx], time_unit)) {
      stop(
        "`axis_units` for the time axis must agree with `time_unit`.",
        call. = FALSE
      )
    }

    axis_units[time_idx] <- time_unit

    if (!is.null(window_time_info)) {
      if (is.null(time_origin)) {
        time_origin <- window_time_info$origin
      } else if (!.pp_same_time_origin(
        time_origin,
        window_time_info$origin
      )) {
        stop(
          "`time_origin` is not consistent with the time origin stored in `window`.",
          call. = FALSE
        )
      }
    }
  }

  if (!is.null(window)) {
    window_units <- window$axis_info$unit[
      match(coord_names, window$axis_info$name)
    ]

    if (!identical(axis_units, window_units)) {
      stop(
        "`axis_units` is not consistent with the units stored in `window`.",
        call. = FALSE
      )
    }
  }

  if (!is.null(crs)) {
    if (!is.character(crs) || length(crs) != 1L || is.na(crs) || !nzchar(crs)) {
      stop(
        "`crs` must be NULL or a single non-empty character string.",
        call. = FALSE
      )
    }
  }

  if (!is.null(window)) {
    if (is.null(window$crs)) {
      if (!is.null(crs)) {
        stop(
          "The supplied `window` has no CRS. Define its CRS in `pp_window()` before attaching it to a pattern.",
          call. = FALSE
        )
      }
    } else if (is.null(crs)) {
      crs <- window$crs
    } else if (!identical(crs, window$crs)) {
      stop(
        "`crs` is not consistent with the CRS stored in `window`.",
        call. = FALSE
      )
    }
  }

  coords_matrix <- matrix(
    NA_real_,
    nrow = nrow(data),
    ncol = d,
    dimnames = list(NULL, coord_names)
  )

  time_info <- NULL

  for (j in seq_len(d)) {
    variable <- data[[coord_idx[j]]]

    if (axis_roles[j] == "time") {
      converted <- .pp_convert_time_coordinate(
        x = variable,
        time_unit = time_unit,
        time_origin = time_origin,
        axis_name = coord_names[j]
      )

      coords_matrix[, j] <- converted$values

      time_info <- list(
        axis = coord_names[j],
        origin = converted$origin,
        original_class = converted$original_class,
        unit = converted$unit
      )
    } else {
      if (!is.numeric(variable)) {
        stop(
          sprintf(
            "Spatial coordinate `%s` must be numeric.",
            coord_names[j]
          ),
          call. = FALSE
        )
      }

      variable <- as.numeric(variable)

      if (any(!is.finite(variable))) {
        stop(
          sprintf(
            "Spatial coordinate `%s` contains missing or non-finite values.",
            coord_names[j]
          ),
          call. = FALSE
        )
      }

      coords_matrix[, j] <- variable
    }
  }

  if (!is.null(window) && nrow(coords_matrix) > 0L) {
    window_bounds <- window$bounds[, coord_names, drop = FALSE]

    outside <- sweep(
      coords_matrix,
      2L,
      window_bounds["lower", ],
      FUN = "<"
    ) |
      sweep(
        coords_matrix,
        2L,
        window_bounds["upper", ],
        FUN = ">"
      )

    if (any(outside)) {
      bad_rows <- which(rowSums(outside) > 0L)
      bad_axes <- colnames(coords_matrix)[colSums(outside) > 0L]

      stop(
        paste0(
          "Some observed points lie outside the supplied observation window. ",
          "Affected axes: ", paste(bad_axes, collapse = ", "),
          ". First affected row(s): ",
          paste(utils::head(bad_rows, 5L), collapse = ", "),
          "."
        ),
        call. = FALSE
      )
    }
  }

  if (is.null(marks)) {
    mark_idx <- integer(0)
  } else if (is.character(marks) &&
             length(marks) == 1L &&
             identical(marks, "auto")) {
    mark_idx <- setdiff(seq_len(ncol(data)), coord_idx)
  } else {
    mark_idx <- .pp_resolve_columns(data, marks, "marks")
  }

  if (length(mark_idx) > 0L && any(mark_idx %in% coord_idx)) {
    overlap <- names(data)[intersect(mark_idx, coord_idx)]

    stop(
      sprintf(
        "Coordinate and mark columns must be distinct. Overlap: %s.",
        paste(overlap, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (length(mark_idx) > 0L && anyDuplicated(mark_idx)) {
    stop("`marks` cannot contain duplicated columns.", call. = FALSE)
  }

  marks_data <- NULL

  if (length(mark_idx) > 0L) {
    marks_data <- data[, mark_idx, drop = FALSE]
    rownames(marks_data) <- NULL
  }

  axis_info <- data.frame(
    name = coord_names,
    role = axis_roles,
    unit = axis_units,
    stringsAsFactors = FALSE
  )

  has_marks <- !is.null(marks_data) && ncol(marks_data) > 0L
  is_st <- n_time == 1L

  subclass <- if (is_st && has_marks) {
    "pp_st_marked"
  } else if (is_st) {
    "pp_st"
  } else if (has_marks) {
    "pp_spatial_marked"
  } else {
    "pp_spatial"
  }

  ignored_columns <- setdiff(
    names(data),
    c(coord_names, if (is.null(marks_data)) character(0) else names(marks_data))
  )

  out <- list(
    coords = coords_matrix,
    marks = marks_data,
    window = window,
    axis_info = axis_info,
    time_info = time_info,
    crs = crs,
    metadata = list(
      ignored_columns = ignored_columns,
      original_column_names = names(data)
    )
  )

  class(out) <- c(subclass, "pp")

  out
}


# Internal helpers ---------------------------------------------------------

.pp_pattern_info <- function(x) {
  is_st <- inherits(x, "pp_st") || inherits(x, "pp_st_marked")
  is_marked <- inherits(x, "pp_spatial_marked") ||
    inherits(x, "pp_st_marked")

  label <- if (is_st && is_marked) {
    "Marked spatio-temporal point pattern"
  } else if (is_st) {
    "Spatio-temporal point pattern"
  } else if (is_marked) {
    "Marked spatial point pattern"
  } else {
    "Spatial point pattern"
  }

  list(
    is_st = is_st,
    is_marked = is_marked,
    label = label,
    subclass = class(x)[1]
  )
}


.pp_format_number <- function(x, digits = 4) {
  if (length(x) == 0L || is.na(x)) {
    return("NA")
  }

  formatC(x, format = "fg", digits = digits, flag = "#")
}


.pp_axis_label <- function(name, unit = NA_character_) {
  if (is.na(unit) || !nzchar(unit)) {
    return(name)
  }

  paste0(name, " [", unit, "]")
}


.pp_observed_bounds <- function(x) {
  coords <- x$coords
  axis_names <- colnames(coords)

  if (nrow(coords) == 0L) {
    out <- matrix(
      NA_real_,
      nrow = 2L,
      ncol = ncol(coords),
      dimnames = list(c("lower", "upper"), axis_names)
    )

    return(out)
  }

  out <- rbind(
    lower = apply(coords, 2L, min),
    upper = apply(coords, 2L, max)
  )

  storage.mode(out) <- "double"
  out
}


.pp_window_bounds <- function(window, axis_names) {
  if (is.null(window)) {
    return(NULL)
  }

  bounds <- window$bounds

  if (is.null(bounds)) {
    stop(
      "The supplied `pp_window` object must contain a `bounds` component.",
      call. = FALSE
    )
  }

  bounds <- as.matrix(bounds)

  if (!is.numeric(bounds) || nrow(bounds) != 2L) {
    stop(
      "`window$bounds` must be a numeric matrix with two rows: lower and upper.",
      call. = FALSE
    )
  }

  if (is.null(colnames(bounds))) {
    if (ncol(bounds) != length(axis_names)) {
      stop(
        "`window$bounds` has an incompatible number of dimensions.",
        call. = FALSE
      )
    }

    colnames(bounds) <- axis_names
  }

  missing_axes <- setdiff(axis_names, colnames(bounds))

  if (length(missing_axes) > 0L) {
    stop(
      sprintf(
        "`window$bounds` is missing coordinate axes: %s.",
        paste(missing_axes, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  bounds <- bounds[, axis_names, drop = FALSE]
  rownames(bounds) <- c("lower", "upper")

  if (any(!is.finite(bounds))) {
    stop("`window$bounds` must contain finite values.", call. = FALSE)
  }

  if (any(bounds["lower", ] > bounds["upper", ])) {
    stop(
      "Every lower window bound must be less than or equal to its upper bound.",
      call. = FALSE
    )
  }

  bounds
}


.pp_effective_bounds <- function(x) {
  axis_names <- x$axis_info$name

  if (is.null(x$window)) {
    return(list(
      bounds = .pp_observed_bounds(x),
      supplied = FALSE,
      source = "Observed-point bounding box"
    ))
  }

  list(
    bounds = .pp_window_bounds(x$window, axis_names),
    supplied = TRUE,
    source = "Specified observation window"
  )
}


.pp_time_description <- function(x) {
  time_idx <- which(x$axis_info$role == "time")

  if (length(time_idx) != 1L) {
    return(NULL)
  }

  time_name <- x$axis_info$name[time_idx]
  time_unit <- x$axis_info$unit[time_idx]

  if (is.null(x$time_info)) {
    return(.pp_axis_label(time_name, time_unit))
  }

  origin <- x$time_info$origin

  if (inherits(origin, "Date")) {
    origin_label <- format(origin, "%Y-%m-%d")
    return(paste0(time_name, " [", time_unit, " since ", origin_label, "]"))
  }

  if (inherits(origin, c("POSIXct", "POSIXt"))) {
    origin_label <- format(origin, "%Y-%m-%d %H:%M:%S")
    return(paste0(time_name, " [", time_unit, " since ", origin_label, "]"))
  }

  if (is.numeric(origin) && length(origin) == 1L) {
    return(paste0(time_name, " [", time_unit, "; origin = ", origin, "]"))
  }

  .pp_axis_label(time_name, time_unit)
}


.pp_restore_time <- function(values, time_info) {
  if (is.null(time_info) || length(values) == 0L) {
    return(NULL)
  }

  origin <- time_info$origin
  unit <- time_info$unit

  seconds_per_unit <- switch(
    unit,
    secs = 1,
    mins = 60,
    hours = 3600,
    days = 86400,
    weeks = 604800,
    NA_real_
  )

  if (inherits(origin, "Date")) {
    days_per_unit <- seconds_per_unit / 86400
    return(origin + values * days_per_unit)
  }

  if (inherits(origin, c("POSIXct", "POSIXt"))) {
    return(origin + values * seconds_per_unit)
  }

  NULL
}


.pp_mark_type <- function(x) {
  if (is.numeric(x)) {
    return("numeric")
  }

  if (inherits(x, "Date")) {
    return("date")
  }

  if (inherits(x, c("POSIXct", "POSIXt"))) {
    return("datetime")
  }

  if (is.factor(x) || is.character(x) || is.logical(x)) {
    return("categorical")
  }

  class(x)[1]
}


.pp_mark_brief <- function(marks) {
  if (is.null(marks) || ncol(marks) == 0L) {
    return(NULL)
  }

  paste(
    vapply(
      names(marks),
      function(nm) {
        paste0(nm, ": ", .pp_mark_type(marks[[nm]]))
      },
      character(1)
    ),
    collapse = "; "
  )
}


.pp_coordinate_summary <- function(coords, axis_info) {
  safe_stats <- function(v) {
    if (length(v) == 0L) {
      return(c(
        min = NA_real_,
        q1 = NA_real_,
        median = NA_real_,
        mean = NA_real_,
        sd = NA_real_,
        q3 = NA_real_,
        max = NA_real_
      ))
    }

    q <- stats::quantile(
      v,
      probs = c(0.25, 0.50, 0.75),
      names = FALSE,
      na.rm = TRUE
    )

    c(
      min = min(v, na.rm = TRUE),
      q1 = q[1],
      median = q[2],
      mean = mean(v, na.rm = TRUE),
      sd = stats::sd(v, na.rm = TRUE),
      q3 = q[3],
      max = max(v, na.rm = TRUE)
    )
  }

  stats_matrix <- t(
    vapply(
      seq_len(ncol(coords)),
      function(j) safe_stats(coords[, j]),
      numeric(7)
    )
  )

  out <- data.frame(
    axis = axis_info$name,
    role = axis_info$role,
    unit = axis_info$unit,
    stats_matrix,
    row.names = NULL,
    check.names = FALSE
  )

  out
}


.pp_mark_summary <- function(x) {
  mark_type <- .pp_mark_type(x)

  if (identical(mark_type, "numeric")) {
    if (length(x) == 0L) {
      stats_out <- c(
        n = 0,
        missing = 0,
        min = NA_real_,
        q1 = NA_real_,
        median = NA_real_,
        mean = NA_real_,
        q3 = NA_real_,
        max = NA_real_
      )
    } else {
      q <- stats::quantile(
        x,
        probs = c(0.25, 0.50, 0.75),
        names = FALSE,
        na.rm = TRUE
      )

      stats_out <- c(
        n = length(x),
        missing = sum(is.na(x)),
        min = min(x, na.rm = TRUE),
        q1 = q[1],
        median = q[2],
        mean = mean(x, na.rm = TRUE),
        q3 = q[3],
        max = max(x, na.rm = TRUE)
      )
    }

    return(list(
      type = mark_type,
      statistics = stats_out
    ))
  }

  if (mark_type %in% c("date", "datetime")) {
    return(list(
      type = mark_type,
      n = length(x),
      missing = sum(is.na(x)),
      min = if (all(is.na(x))) NA else min(x, na.rm = TRUE),
      max = if (all(is.na(x))) NA else max(x, na.rm = TRUE)
    ))
  }

  if (mark_type == "categorical") {
    values <- as.character(x)
    values[is.na(values)] <- "<NA>"

    counts <- sort(table(values), decreasing = TRUE)

    frequency <- data.frame(
      level = names(counts),
      count = as.integer(counts),
      percentage = as.numeric(counts) / sum(counts) * 100,
      row.names = NULL,
      check.names = FALSE
    )

    return(list(
      type = mark_type,
      frequency = frequency
    ))
  }

  list(
    type = mark_type,
    message = "No automatic summary is available for this mark type."
  )
}


.pp_print_bounds <- function(bounds, axis_info, digits = 4) {
  for (j in seq_len(ncol(bounds))) {
    axis_label <- .pp_axis_label(
      axis_info$name[j],
      axis_info$unit[j]
    )

    lower <- .pp_format_number(bounds["lower", j], digits)
    upper <- .pp_format_number(bounds["upper", j], digits)

    cat(sprintf("  %s: [%s, %s]\n", axis_label, lower, upper))
  }
}


#' Methods for multidimensional point-pattern objects
#'
#' Print, summarise, and plot objects created by [pp_pattern()]. The methods
#' adapt their output to the spatial dimension, the presence of a temporal
#' coordinate, and the number and type of marks.
#'
#' @details
#' `print.pp()` reports the principal structural information. `summary.pp()`
#' computes coordinate summaries, temporal ranges, mark summaries, and window
#' information, returning an object of class `"summary_pp"`.
#' `print.summary_pp()` formats that object.
#'
#' `plot.pp()` produces exploratory displays. Two-dimensional panels use base
#' graphics, while three-dimensional panels use `plot3D::scatter3D()`. For a
#' spatio-temporal pattern, spatial and temporal projections are displayed. A
#' continuous mark is represented with a colour scale and a categorical mark
#' with a legend. When more than three spatial axes are available, `axes` must
#' select the two or three axes to display.
#'
#' If no observation window is stored, printing and plotting use the observed
#' bounding box only as a descriptive range; it is not treated as a statistical
#' observation window.
#'
#' @param x A `pp` object for `print()` and `plot()`, or a `summary_pp` object
#'   for `print.summary_pp()`.
#' @param object A `pp` object to be summarised.
#' @param digits Number of significant digits used when printing numeric values.
#' @param max_mark_levels Maximum number of levels displayed for each
#'   categorical mark. Additional levels are omitted from the printed output.
#' @param mark Optional mark name or one-based numeric mark position. If `NULL`,
#'   every mark is displayed, with one row of panels per mark where applicable.
#' @param axes Optional character or numeric vector selecting two or three
#'   spatial axes. Selection is required when more than three spatial axes are
#'   stored. Temporal axes are handled separately and must not be selected here.
#' @param palette_continuous Name of an HCL palette used for continuous marks.
#' @param palette_discrete Name of an HCL palette used for categorical marks.
#' @param missing_col Colour assigned to missing mark values.
#' @param point_col Colour used for unmarked events.
#' @param point_size Positive point-size multiplier.
#' @param point_pch Plotting symbol used for events.
#' @param theta,phi Viewing angles, in degrees, passed to three-dimensional
#'   displays.
#' @param show_window Logical value indicating whether the supplied observation
#'   window, or the descriptive observed bounding box when no window is stored,
#'   should be drawn.
#' @param cex_axis,cex_lab,cex_main,cex_legend Positive scale factors for axis
#'   annotations, axis labels, panel titles, and legends/colour bars.
#' @param colorbar_ticks Approximate positive number of ticks displayed on a
#'   continuous colour bar.
#' @param colorbar_digits Optional non-negative integer fixing the number of
#'   decimal places in continuous colour-bar labels. If `NULL`, a suitable
#'   value is selected automatically.
#' @param ... Additional arguments. They are currently ignored.
#'
#' @return `print.pp()`, `print.summary_pp()`, and `plot.pp()` invisibly return
#'   `x`. `summary.pp()` returns an object of class `"summary_pp"`.
#'
#' @examples
#' w <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' set.seed(2)
#' x <- pp_pattern(
#'   data = data.frame(
#'     x = runif(30),
#'     y = runif(30),
#'     mark = rnorm(30)
#'   ),
#'   coords = c("x", "y"),
#'   marks = "mark",
#'   window = w
#' )
#'
#' print(x)
#' sx <- summary(x)
#' print(sx)
#' plot(x, mark = "mark")
#'
#' @seealso [pp_pattern()] and [pp_window()].
#' @family point-pattern data
#' @rdname pp-methods
#' @method print pp
#' @export
print.pp <- function(x, digits = 4, ...) {
  info <- .pp_pattern_info(x)
  axis_info <- x$axis_info
  bounds_info <- .pp_effective_bounds(x)

  spatial_axes <- axis_info$name[axis_info$role == "space"]
  spatial_labels <- vapply(
    which(axis_info$role == "space"),
    function(j) .pp_axis_label(axis_info$name[j], axis_info$unit[j]),
    character(1)
  )

  cat("<", info$subclass, ">\n", sep = "")
  cat(info$label, "\n\n", sep = "")

  cat("Points: ", nrow(x$coords), "\n", sep = "")
  cat(
    "Spatial dimensions: ",
    length(spatial_axes),
    "\n",
    sep = ""
  )
  cat("Spatial axes: ", paste(spatial_labels, collapse = ", "), "\n", sep = "")

  if (info$is_st) {
    cat("Time axis: ", .pp_time_description(x), "\n", sep = "")
  }

  if (info$is_marked) {
    cat(
      "Marks: ",
      ncol(x$marks),
      " (",
      .pp_mark_brief(x$marks),
      ")\n",
      sep = ""
    )
  }

  if (!is.null(x$crs)) {
    cat("CRS: ", x$crs, "\n", sep = "")
  }

  cat("\n")

  if (bounds_info$supplied) {
    cat("Observation window: specified hyperrectangular window\n")
  } else {
    cat("Observation window: not supplied\n")
    cat("Enclosing bounding box of observed points:\n")
  }

  .pp_print_bounds(
    bounds = bounds_info$bounds,
    axis_info = axis_info,
    digits = digits
  )

  invisible(x)
}


#' @rdname pp-methods
#' @method summary pp
#' @export
summary.pp <- function(object, ...) {
  info <- .pp_pattern_info(object)
  axis_info <- object$axis_info
  bounds_info <- .pp_effective_bounds(object)

  marks_summary <- NULL

  if (info$is_marked) {
    marks_summary <- lapply(object$marks, .pp_mark_summary)
  }

  time_summary <- NULL

  if (info$is_st) {
    time_idx <- which(axis_info$role == "time")
    time_values <- object$coords[, time_idx]

    original_range <- NULL

    if (!is.null(object$time_info) && length(time_values) > 0L) {
      restored <- .pp_restore_time(
        values = range(time_values),
        time_info = object$time_info
      )

      if (!is.null(restored)) {
        original_range <- restored
      }
    }

    time_summary <- list(
      axis = axis_info$name[time_idx],
      unit = axis_info$unit[time_idx],
      numeric_range = if (length(time_values) == 0L) {
        c(NA_real_, NA_real_)
      } else {
        range(time_values)
      },
      original_range = original_range,
      origin = if (is.null(object$time_info)) NULL else object$time_info$origin
    )
  }

  window_volume <- NA_real_

  if (bounds_info$supplied) {
    window_volume <- prod(
      bounds_info$bounds["upper", ] - bounds_info$bounds["lower", ]
    )
  }

  out <- list(
    subclass = info$subclass,
    label = info$label,
    n_points = nrow(object$coords),
    dimension = ncol(object$coords),
    spatial_dimension = sum(axis_info$role == "space"),
    axis_info = axis_info,
    coordinate_summary = .pp_coordinate_summary(object$coords, axis_info),
    window = list(
      supplied = bounds_info$supplied,
      source = bounds_info$source,
      bounds = bounds_info$bounds,
      hypervolume = window_volume
    ),
    time = time_summary,
    marks = marks_summary,
    crs = object$crs
  )

  class(out) <- "summary_pp"
  out
}


#' @rdname pp-methods
#' @method print summary_pp
#' @export
print.summary_pp <- function(
    x,
    digits = 4,
    max_mark_levels = 10,
    ...
) {
  cat(x$label, "\n\n", sep = "")

  cat("Points: ", x$n_points, "\n", sep = "")
  cat("Total dimensions: ", x$dimension, "\n", sep = "")
  cat("Spatial dimensions: ", x$spatial_dimension, "\n", sep = "")

  if (!is.null(x$crs)) {
    cat("CRS: ", x$crs, "\n", sep = "")
  }

  cat("\nCoordinate summaries:\n")

  coord_summary <- x$coordinate_summary

  numeric_cols <- c("min", "q1", "median", "mean", "sd", "q3", "max")

  for (nm in numeric_cols) {
    coord_summary[[nm]] <- vapply(
      coord_summary[[nm]],
      .pp_format_number,
      character(1),
      digits = digits
    )
  }

  if (!is.null(coord_summary$unit)) {
    coord_summary$unit[is.na(coord_summary$unit)] <- ""
  }
  print(coord_summary, row.names = FALSE)

  cat("\n")

  if (x$window$supplied) {
    cat("Observation window: specified hyperrectangular window\n")
    .pp_print_bounds(
      bounds = x$window$bounds,
      axis_info = x$axis_info,
      digits = digits
    )

    cat(
      "Window hypervolume: ",
      .pp_format_number(x$window$hypervolume, digits),
      "\n",
      sep = ""
    )
  } else {
    cat("Observation window: not supplied\n")
    cat("Observed-point bounding box:\n")

    .pp_print_bounds(
      bounds = x$window$bounds,
      axis_info = x$axis_info,
      digits = digits
    )
  }

  if (!is.null(x$time)) {
    cat("\nTime information:\n")
    cat("  Axis: ", x$time$axis, "\n", sep = "")
    cat(
      "  Numeric range: [",
      .pp_format_number(x$time$numeric_range[1], digits),
      ", ",
      .pp_format_number(x$time$numeric_range[2], digits),
      "] ",
      x$time$unit,
      "\n",
      sep = ""
    )

    if (!is.null(x$time$original_range)) {
      cat(
        "  Original range: ",
        format(x$time$original_range[1]),
        " to ",
        format(x$time$original_range[2]),
        "\n",
        sep = ""
      )
    }
  }

  if (is.null(x$marks)) {
    cat("\nMarks: none\n")
    return(invisible(x))
  }

  cat("\nMarks:\n")

  for (nm in names(x$marks)) {
    mark_summary <- x$marks[[nm]]

    cat("\n", nm, " (", mark_summary$type, ")\n", sep = "")

    if (identical(mark_summary$type, "numeric")) {
      statistics <- mark_summary$statistics

      stats_out <- data.frame(
        n = statistics["n"],
        missing = statistics["missing"],
        min = .pp_format_number(statistics["min"], digits),
        q1 = .pp_format_number(statistics["q1"], digits),
        median = .pp_format_number(statistics["median"], digits),
        mean = .pp_format_number(statistics["mean"], digits),
        q3 = .pp_format_number(statistics["q3"], digits),
        max = .pp_format_number(statistics["max"], digits),
        row.names = NULL,
        check.names = FALSE
      )

      print(stats_out, row.names = FALSE)
    } else if (mark_summary$type %in% c("date", "datetime")) {
      cat("  Observations: ", mark_summary$n, "\n", sep = "")
      cat("  Missing: ", mark_summary$missing, "\n", sep = "")
      cat("  Range: ", format(mark_summary$min), " to ",
          format(mark_summary$max), "\n", sep = "")
    } else if (identical(mark_summary$type, "categorical")) {
      frequency <- mark_summary$frequency

      if (nrow(frequency) == 0L) {
        cat("  No observed mark values.\n")
      } else {
        shown <- utils::head(frequency, max_mark_levels)

        shown$percentage <- vapply(
          shown$percentage,
          .pp_format_number,
          character(1),
          digits = digits
        )

        print(shown, row.names = FALSE)

        if (nrow(frequency) > max_mark_levels) {
          cat(
            "  ... ",
            nrow(frequency) - max_mark_levels,
            " additional level(s) not displayed.\n",
            sep = ""
          )
        }
      }
    } else {
      cat("  ", mark_summary$message, "\n", sep = "")
    }
  }

  invisible(x)
}



#' Create a multidimensional observation window
#'
#' Creates a finite hyperrectangular observation window for a spatial or
#' spatio-temporal point process. Each coordinate axis is represented by one
#' lower and one upper bound.
#'
#' @details
#' The window is stored as a numeric matrix with rows `"lower"` and `"upper"`.
#' Its hypervolume is
#'
#' \deqn{|W| = \prod_{j=1}^{d} (b_j-a_j),}
#'
#' where \eqn{a_j} and \eqn{b_j} are the lower and upper bounds of axis \eqn{j}.
#' This quantity is a spatial volume for a purely spatial window and a
#' space-time hypervolume for a spatio-temporal window.
#'
#' Temporal bounds are always numeric and are interpreted on the internal scale
#' defined by `time_unit` and `time_origin`. For example, bounds `c(0, 30)` with
#' a `Date` origin and `time_unit = "days"` represent the first 30 days after
#' that origin. At most one temporal axis is supported, and a temporal window
#' must also contain at least one spatial axis.
#'
#' The window provides authoritative metadata when passed to [pp_pattern()] or
#' [pp_covariate()]. Those objects must use compatible axis names, roles, units,
#' temporal origin, and coordinate reference system.
#'
#' @param bounds Optional numeric matrix or data frame with exactly two rows and
#'   at least one column. The first row contains lower bounds and the second row
#'   contains upper bounds. Columns identify coordinate axes. Supply either
#'   `bounds` or both `lower` and `upper`.
#' @param lower,upper Optional numeric vectors of equal positive length. Every
#'   lower bound must be finite and strictly smaller than the corresponding
#'   upper bound. Supply both arguments together and do not combine them with
#'   `bounds`.
#' @param axes Optional character vector of unique coordinate-axis names. It is
#'   required when the supplied bounds are unnamed and may be used to reorder
#'   named bounds. If omitted, existing names are retained.
#' @param axis_roles Optional character vector describing each axis. Allowed
#'   values are `"space"` and `"time"`. It may be unnamed, in axis order, or
#'   named with `axes`. By default all axes are spatial. At most one temporal
#'   axis is allowed, and a temporal window must contain at least one spatial
#'   axis.
#' @param axis_units Optional character vector containing axis units. It may be
#'   unnamed, in axis order, or named with `axes`. Missing units may be supplied
#'   as `NA_character_`.
#' @param time_unit Optional unit of the temporal axis. Allowed values are
#'   `"days"`, `"hours"`, `"mins"`, `"secs"`, and `"weeks"`. If omitted for a
#'   temporal window, the temporal entry of `axis_units` is used when available;
#'   otherwise the default is `"days"`.
#' @param time_origin Optional single finite numeric value, `Date`, or
#'   `POSIXct`/`POSIXt` object defining the origin of the internal temporal
#'   scale. It defaults to zero for a temporal window. It cannot be supplied for
#'   a purely spatial window.
#' @param crs Optional non-empty character string describing the coordinate
#'   reference system of the spatial axes.
#'
#' @return An object of class `"pp_window"` with components:
#'
#' - `bounds`: the two-row numeric bounds matrix;
#' - `dimension`: number of coordinate axes;
#' - `axis_info`: axis names, roles, and units;
#' - `time_info`: temporal metadata, or `NULL` for a spatial window;
#' - `hypervolume`: product of the axis widths;
#' - `crs`: coordinate reference system, or `NULL`.
#'
#' @examples
#' # Construct a two-dimensional window from lower and upper vectors
#' w2 <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 10, y = 5),
#'   axis_units = c(x = "km", y = "km")
#' )
#' w2$hypervolume
#'
#' # The same interface also supports higher-dimensional spatial windows
#' w3 <- pp_window(
#'   bounds = rbind(
#'     lower = c(x = 0, y = 0, z = 0),
#'     upper = c(x = 1, y = 1, z = 2)
#'   )
#' )
#'
#' # A two-dimensional space-time window measured in days from a date origin
#' wst <- pp_window(
#'   lower = c(x = 0, y = 0, time = 0),
#'   upper = c(x = 100, y = 100, time = 365),
#'   axis_roles = c(x = "space", y = "space", time = "time"),
#'   axis_units = c(x = "km", y = "km", time = "days"),
#'   time_origin = as.Date("2020-01-01"),
#'   crs = "local Cartesian coordinates"
#' )
#'
#' @seealso
#' [pp_pattern()] for event patterns, [pp_covariate()] for covariate supports,
#' and [pp_grid()] for regular grids over a window.
#'
#' @family point-pattern data
#' @export
pp_window <- function(
    bounds = NULL,
    lower = NULL,
    upper = NULL,
    axes = NULL,
    axis_roles = NULL,
    axis_units = NULL,
    time_unit = NULL,
    time_origin = NULL,
    crs = NULL
) {
  if (!is.null(bounds) && (!is.null(lower) || !is.null(upper))) {
    stop(
      "Supply either `bounds` or both `lower` and `upper`, not both.",
      call. = FALSE
    )
  }

  if (is.null(bounds)) {
    if (is.null(lower) || is.null(upper)) {
      stop(
        "Supply either `bounds` or both `lower` and `upper`.",
        call. = FALSE
      )
    }

    if (!is.numeric(lower) || !is.numeric(upper)) {
      stop("`lower` and `upper` must be numeric vectors.", call. = FALSE)
    }

    if (length(lower) == 0L || length(upper) == 0L) {
      stop("`lower` and `upper` cannot be empty.", call. = FALSE)
    }

    if (length(lower) != length(upper)) {
      stop(
        "`lower` and `upper` must have the same length.",
        call. = FALSE
      )
    }

    lower_names <- names(lower)
    upper_names <- names(upper)

    if (!is.null(lower_names) && !is.null(upper_names)) {
      if (!setequal(lower_names, upper_names)) {
        stop(
          "`lower` and `upper` must have identical names.",
          call. = FALSE
        )
      }

      upper <- upper[lower_names]
      bounds <- rbind(lower = lower, upper = upper)
    } else {
      bounds <- rbind(lower = lower, upper = upper)
    }
  } else {
    if (!is.matrix(bounds) && !is.data.frame(bounds)) {
      stop(
        "`bounds` must be a numeric matrix or data frame.",
        call. = FALSE
      )
    }

    bounds <- as.matrix(bounds)

    if (!is.numeric(bounds)) {
      stop("`bounds` must contain numeric values only.", call. = FALSE)
    }

    if (nrow(bounds) != 2L) {
      stop(
        "`bounds` must have exactly two rows: lower and upper.",
        call. = FALSE
      )
    }

    if (ncol(bounds) == 0L) {
      stop("`bounds` must contain at least one dimension.", call. = FALSE)
    }
  }

  if (any(!is.finite(bounds))) {
    stop("Window bounds must be finite.", call. = FALSE)
  }

  d <- ncol(bounds)

  if (!is.null(axes)) {
    if (!is.character(axes) || length(axes) != d) {
      stop(
        sprintf("`axes` must be a character vector of length %d.", d),
        call. = FALSE
      )
    }

    if (any(!nzchar(axes)) || anyDuplicated(axes)) {
      stop(
        "`axes` must contain unique, non-empty names.",
        call. = FALSE
      )
    }
  }

  bound_names <- colnames(bounds)

  if (!is.null(bound_names) && any(!nzchar(bound_names))) {
    bound_names <- NULL
  }

  if (!is.null(axes) && !is.null(bound_names)) {
    if (!setequal(axes, bound_names)) {
      stop(
        "`axes` and the column names of `bounds` must identify the same axes.",
        call. = FALSE
      )
    }

    bounds <- bounds[, axes, drop = FALSE]
    axis_names <- axes
  } else if (!is.null(axes)) {
    axis_names <- axes
    colnames(bounds) <- axis_names
  } else if (!is.null(bound_names)) {
    if (anyDuplicated(bound_names)) {
      stop(
        "The columns of `bounds` must have unique names.",
        call. = FALSE
      )
    }

    axis_names <- bound_names
  } else {
    axis_names <- paste0("d", seq_len(d))
    colnames(bounds) <- axis_names
  }

  rownames(bounds) <- c("lower", "upper")

  if (any(bounds["lower", ] >= bounds["upper", ])) {
    invalid_axes <- axis_names[
      bounds["lower", ] >= bounds["upper", ]
    ]

    stop(
      sprintf(
        "Each lower bound must be strictly smaller than its upper bound. Problematic axis/axes: %s.",
        paste(invalid_axes, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  axis_roles <- .pp_align_axis_argument(
    x = axis_roles,
    axis_names = axis_names,
    arg_name = "axis_roles"
  )

  if (is.null(axis_roles)) {
    axis_roles <- rep("space", d)
  }

  axis_roles <- tolower(trimws(axis_roles))

  if (!all(axis_roles %in% c("space", "time"))) {
    stop(
      "`axis_roles` may contain only 'space' and 'time'.",
      call. = FALSE
    )
  }

  n_time <- sum(axis_roles == "time")

  if (n_time > 1L) {
    stop(
      "A `pp_window` currently supports at most one temporal axis.",
      call. = FALSE
    )
  }

  if (n_time == 1L && sum(axis_roles == "space") == 0L) {
    stop(
      "A spatio-temporal window must contain at least one spatial axis.",
      call. = FALSE
    )
  }

  axis_units <- .pp_align_axis_argument(
    x = axis_units,
    axis_names = axis_names,
    arg_name = "axis_units"
  )

  if (is.null(axis_units)) {
    axis_units <- rep(NA_character_, d)
  }

  time_info <- NULL

  if (n_time == 1L) {
    time_idx <- which(axis_roles == "time")

    if (is.null(time_unit) &&
        !is.na(axis_units[time_idx]) &&
        nzchar(axis_units[time_idx])) {
      time_unit <- axis_units[time_idx]
    }

    if (is.null(time_unit)) {
      time_unit <- "days"
    }

    time_unit <- .pp_validate_time_unit(time_unit)

    if (!is.na(axis_units[time_idx]) &&
        nzchar(axis_units[time_idx]) &&
        !identical(axis_units[time_idx], time_unit)) {
      stop(
        "`axis_units` for the time axis must agree with `time_unit`.",
        call. = FALSE
      )
    }

    axis_units[time_idx] <- time_unit

    if (is.null(time_origin)) {
      time_origin <- 0
    }

    if (is.numeric(time_origin)) {
      if (length(time_origin) != 1L || !is.finite(time_origin)) {
        stop(
          "`time_origin` must be a single finite numeric value, Date, or POSIXct/POSIXt.",
          call. = FALSE
        )
      }
    } else if (inherits(time_origin, "Date")) {
      if (length(time_origin) != 1L || is.na(time_origin)) {
        stop(
          "`time_origin` must be a non-missing single Date value.",
          call. = FALSE
        )
      }
    } else if (inherits(time_origin, c("POSIXct", "POSIXt"))) {
      if (length(time_origin) != 1L || is.na(time_origin)) {
        stop(
          "`time_origin` must be a non-missing single POSIXct/POSIXt value.",
          call. = FALSE
        )
      }
    } else {
      stop(
        "`time_origin` must be numeric, Date, or POSIXct/POSIXt.",
        call. = FALSE
      )
    }

    time_info <- list(
      axis = axis_names[time_idx],
      origin = time_origin,
      original_class = class(time_origin)[1],
      unit = time_unit
    )
  } else if (!is.null(time_unit) || !is.null(time_origin)) {
    stop(
      "`time_unit` and `time_origin` can be supplied only when one axis has role 'time'.",
      call. = FALSE
    )
  }

  if (!is.null(crs)) {
    if (!is.character(crs) || length(crs) != 1L || is.na(crs) || !nzchar(crs)) {
      stop(
        "`crs` must be NULL or a single non-empty character string.",
        call. = FALSE
      )
    }
  }

  axis_info <- data.frame(
    name = axis_names,
    role = axis_roles,
    unit = axis_units,
    stringsAsFactors = FALSE
  )

  out <- list(
    bounds = bounds,
    dimension = d,
    axis_info = axis_info,
    time_info = time_info,
    hypervolume = prod(bounds["upper", ] - bounds["lower", ]),
    crs = crs
  )

  class(out) <- "pp_window"

  out
}



# Internal plotting helpers -------------------------------------------------

.pp_plot_limits <- function(x) {
  x <- sort(as.numeric(x))

  if (length(x) != 2L || any(!is.finite(x))) {
    stop("Plot limits must contain two finite values.", call. = FALSE)
  }

  if (x[1] == x[2]) {
    pad <- max(abs(x[1]) * 0.05, 1)
    x <- x + c(-pad, pad)
  }

  x
}


.pp_plot_bounds <- function(x, axes) {
  out <- .pp_effective_bounds(x)$bounds[, axes, drop = FALSE]

  if (any(!is.finite(out))) {
    stop(
      "A non-empty observation window is required to plot an empty pattern.",
      call. = FALSE
    )
  }

  for (j in seq_len(ncol(out))) {
    out[, j] <- .pp_plot_limits(out[, j])
  }

  out
}


.pp_plot_axis_labels <- function(x, axes) {
  idx <- match(axes, x$axis_info$name)

  vapply(
    idx,
    function(j) {
      .pp_axis_label(
        name = x$axis_info$name[j],
        unit = x$axis_info$unit[j]
      )
    },
    character(1)
  )
}


.pp_plot_resolve_spatial_axes <- function(x, axes = NULL) {
  spatial_axes <- x$axis_info$name[
    x$axis_info$role == "space"
  ]

  n_spatial <- length(spatial_axes)

  if (n_spatial == 0L) {
    stop("The pattern contains no spatial axes.", call. = FALSE)
  }

  n_required <- min(n_spatial, 3L)

  if (is.null(axes)) {
    if (n_spatial > 3L) {
      stop(
        paste0(
          "The pattern has more than three spatial dimensions. ",
          "Specify the three axes to display through `axes`."
        ),
        call. = FALSE
      )
    }

    return(spatial_axes)
  }

  if (is.character(axes)) {
    bad_axes <- setdiff(axes, spatial_axes)

    if (length(bad_axes) > 0L) {
      stop(
        sprintf(
          "`axes` must contain only spatial axes. Invalid axes: %s.",
          paste(bad_axes, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    chosen <- axes
  } else if (is.numeric(axes)) {
    axes <- as.integer(axes)

    if (any(axes < 1L | axes > n_spatial)) {
      stop(
        "`axes` contains indices outside the available spatial axes.",
        call. = FALSE
      )
    }

    chosen <- spatial_axes[axes]
  } else {
    stop(
      "`axes` must be NULL, a character vector, or numeric positions.",
      call. = FALSE
    )
  }

  if (length(chosen) != n_required) {
    stop(
      sprintf(
        "Exactly %d spatial axis/axes must be selected.",
        n_required
      ),
      call. = FALSE
    )
  }

  if (anyDuplicated(chosen)) {
    stop("`axes` cannot contain duplicated axes.", call. = FALSE)
  }

  chosen
}


.pp_plot_resolve_marks <- function(x, mark = NULL) {
  if (is.null(x$marks) || ncol(x$marks) == 0L) {
    if (!is.null(mark)) {
      stop(
        "`mark` can be used only with marked point patterns.",
        call. = FALSE
      )
    }

    return(character(0))
  }

  available <- names(x$marks)

  if (is.null(mark)) {
    return(available)
  }

  if (is.character(mark)) {
    missing_marks <- setdiff(mark, available)

    if (length(missing_marks) > 0L) {
      stop(
        sprintf(
          "Unknown mark(s): %s.",
          paste(missing_marks, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    return(mark)
  }

  if (is.numeric(mark)) {
    mark <- as.integer(mark)

    if (any(mark < 1L | mark > length(available))) {
      stop(
        "`mark` contains positions outside the available marks.",
        call. = FALSE
      )
    }

    return(available[mark])
  }

  stop(
    "`mark` must be NULL, a mark name, or a numeric mark position.",
    call. = FALSE
  )
}


.pp_plot_pretty_ticks <- function(limits, n = 5L) {
  limits <- .pp_plot_limits(limits)

  ticks <- pretty(limits, n = n)

  tol <- diff(limits) * 1e-9

  ticks <- ticks[
    ticks >= limits[1] - tol &
      ticks <= limits[2] + tol
  ]

  if (length(ticks) < 2L) {
    ticks <- seq(limits[1], limits[2], length.out = n)
  }

  unique(ticks)
}


.pp_plot_tick_digits <- function(ticks) {
  ticks <- sort(unique(as.numeric(ticks)))

  if (length(ticks) < 2L) {
    return(2L)
  }

  step <- min(diff(ticks))

  if (!is.finite(step) || step <= 0) {
    return(2L)
  }

  if (step >= 10) {
    return(0L)
  }

  if (step >= 1) {
    return(1L)
  }

  max(1L, ceiling(-log10(step)) + 1L)
}


.pp_plot_tick_labels <- function(ticks, digits = NULL) {
  if (is.null(digits)) {
    digits <- .pp_plot_tick_digits(ticks)
  }

  lab <- formatC(
    ticks,
    format = "f",
    digits = digits
  )

  lab <- sub("0+$", "", lab)
  lab <- sub("\\.$", "", lab)

  lab
}


.pp_plot_mark_spec <- function(
    values,
    name,
    palette_continuous = "YlOrRd",
    palette_discrete = "Dark 3",
    missing_col = "grey75",
    colorbar_ticks = 5L,
    colorbar_digits = NULL
) {
  is_continuous <- is.numeric(values) ||
    inherits(values, "Date") ||
    inherits(values, c("POSIXct", "POSIXt"))

  if (is_continuous) {
    values_num <- as.numeric(values)
    observed <- is.finite(values_num)

    palette <- grDevices::hcl.colors(
      n = 80L,
      palette = palette_continuous
    )

    point_colours <- rep(missing_col, length(values_num))

    if (any(observed)) {
      limits <- range(values_num[observed])

      if (limits[1] == limits[2]) {
        pad <- max(abs(limits[1]) * 0.05, 1)
        limits <- limits + c(-pad, pad)
      }

      scaled <- (values_num[observed] - limits[1]) / diff(limits)

      index <- floor(
        scaled * (length(palette) - 1L)
      ) + 1L

      index <- pmax(1L, pmin(length(palette), index))

      point_colours[observed] <- palette[index]

      ticks <- .pp_plot_pretty_ticks(
        limits = limits,
        n = colorbar_ticks
      )

      tick_labels <- .pp_plot_tick_labels(
        ticks = ticks,
        digits = colorbar_digits
      )
    } else {
      limits <- NULL
      ticks <- NULL
      tick_labels <- NULL
    }

    return(list(
      type = "continuous",
      name = name,
      numeric = values_num,
      colours = point_colours,
      limits = limits,
      palette = palette,
      ticks = ticks,
      tick_labels = tick_labels,
      has_missing = any(!observed),
      missing_col = missing_col
    ))
  }

  missing_label <- "Missing"

  labels <- as.character(values)
  labels[is.na(labels)] <- missing_label

  if (is.factor(values)) {
    observed_levels <- intersect(
      levels(values),
      unique(labels[labels != missing_label])
    )
  } else {
    observed_levels <- unique(labels[labels != missing_label])
  }

  palette <- grDevices::hcl.colors(
    n = max(3L, length(observed_levels)),
    palette = palette_discrete
  )[seq_along(observed_levels)]

  names(palette) <- observed_levels

  all_levels <- observed_levels

  if (any(labels == missing_label)) {
    all_levels <- c(all_levels, missing_label)
  }

  level_colours <- c(
    palette,
    setNames(missing_col, missing_label)
  )

  level_colours <- level_colours[all_levels]

  list(
    type = "discrete",
    name = name,
    colours = unname(level_colours[labels]),
    levels = all_levels,
    level_colours = level_colours,
    has_missing = any(labels == missing_label),
    missing_col = missing_col
  )
}


.pp_plot_colourbar_panel <- function(
    spec,
    cex_axis = 1.25,
    cex_main = 1.35
) {
  if (!identical(spec$type, "continuous")) {
    stop("`spec` must describe a continuous mark.", call. = FALSE)
  }

  graphics::par(
    mar = c(3.0, 0.3, 3.8, 4.8)
  )

  graphics::plot.new()

  if (is.null(spec$limits)) {
    graphics::plot.window(
      xlim = c(0, 1),
      ylim = c(0, 1)
    )

    graphics::text(
      x = 0.5,
      y = 0.5,
      labels = "All values\nmissing",
      cex = cex_axis
    )

    graphics::title(
      main = spec$name,
      cex.main = 0.8 * cex_main
    )

    return(invisible(NULL))
  }

  limits <- spec$limits
  range_y <- diff(limits)

  extra_bottom <- if (spec$has_missing) 0.18 * range_y else 0

  graphics::plot.window(
    xlim = c(0, 1),
    ylim = c(limits[1] - extra_bottom, limits[2]),
    xaxs = "i",
    yaxs = "i"
  )

  y_seq <- seq(
    limits[1],
    limits[2],
    length.out = length(spec$palette) + 1L
  )

  for (i in seq_along(spec$palette)) {
    graphics::rect(
      xleft = 0.20,
      ybottom = y_seq[i],
      xright = 0.55,
      ytop = y_seq[i + 1L],
      col = spec$palette[i],
      border = NA
    )
  }

  graphics::rect(
    xleft = 0.20,
    ybottom = limits[1],
    xright = 0.55,
    ytop = limits[2],
    border = "grey25",
    lwd = 1
  )

  graphics::axis(
    side = 4,
    at = spec$ticks,
    labels = spec$tick_labels,
    las = 1,
    cex.axis = cex_axis,
    lwd = 1,
    lwd.ticks = 1
  )

  graphics::title(
    main = spec$name,
    cex.main = 0.75 * cex_main,
    font.main = 1
  )

  if (spec$has_missing) {
    y_missing <- limits[1] - 0.09 * range_y

    graphics::points(
      x = 0.38,
      y = y_missing,
      pch = 15,
      cex = 1.4,
      col = spec$missing_col
    )

    graphics::text(
      x = 0.62,
      y = y_missing,
      labels = "Missing",
      adj = c(0, 0.5),
      cex = cex_axis
    )
  }

  invisible(NULL)
}


.pp_plot_legend_panel <- function(
    spec,
    point_pch = 19,
    cex_legend = 1.2,
    cex_main = 1.35
) {
  if (!identical(spec$type, "discrete")) {
    stop("`spec` must describe a categorical mark.", call. = FALSE)
  }

  graphics::par(
    mar = c(2.0, 0.2, 3.8, 0.2)
  )

  graphics::plot.new()

  graphics::plot.window(
    xlim = c(0, 1),
    ylim = c(0, 1)
  )

  graphics::title(
    main = spec$name,
    cex.main = cex_main,
    font.main = 2
  )

  n_levels <- length(spec$levels)

  legend_cols <- if (n_levels > 8L) 2L else 1L

  graphics::legend(
    "center",
    legend = spec$levels,
    col = unname(spec$level_colours),
    pch = point_pch,
    pt.cex = 1.35,
    cex = cex_legend,
    ncol = legend_cols,
    bty = "n"
  )

  invisible(NULL)
}


.pp_plot_2d_panel <- function(
    coords,
    bounds,
    axis_labels,
    main,
    colours = "black",
    point_size = 1.35,
    point_pch = 19,
    equal_axes = TRUE,
    show_window = TRUE,
    cex_axis = 1.3,
    cex_lab = 1.55,
    cex_main = 1.45
) {
  xlim <- .pp_plot_limits(bounds[, 1])
  ylim <- .pp_plot_limits(bounds[, 2])

  graphics::par(
    mar = c(5.5, 5.6, 4.0, 1.0)
  )

  graphics::plot(
    NA_real_,
    NA_real_,
    type = "n",
    xlim = xlim,
    ylim = ylim,
    axes = FALSE,
    ann = FALSE,
    bty = "n",
    xaxs = "i",
    yaxs = "i",
    asp = if (equal_axes) 1 else NA
  )

  graphics::axis(
    side = 1,
    cex.axis = cex_axis,
    lwd = 1,
    lwd.ticks = 1
  )

  graphics::axis(
    side = 2,
    cex.axis = cex_axis,
    las = 1,
    lwd = 1,
    lwd.ticks = 1
  )

  graphics::mtext(
    text = axis_labels[1],
    side = 1,
    line = 3.4,
    cex = cex_lab
  )

  graphics::mtext(
    text = axis_labels[2],
    side = 2,
    line = 3.7,
    cex = cex_lab
  )

  graphics::title(
    main = main,
    cex.main = cex_main,
    font.main = 2
  )

  if (isTRUE(show_window)) {
    graphics::rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = ylim[2],
      border = "grey35",
      lwd = 1.3
    )
  }

  if (nrow(coords) > 0L) {
    graphics::points(
      x = coords[, 1],
      y = coords[, 2],
      pch = point_pch,
      cex = point_size,
      col = colours
    )
  }

  invisible(NULL)
}


.pp_plot_time_panel <- function(
    x,
    time_axis,
    main = "Cumulative number of events",
    show_window = TRUE,
    cex_axis = 1.3,
    cex_lab = 1.55,
    cex_main = 1.45
) {
  bounds <- .pp_plot_bounds(
    x = x,
    axes = time_axis
  )

  xlim <- .pp_plot_limits(bounds[, 1])
  times <- sort(x$coords[, time_axis])

  graphics::par(
    mar = c(5.5, 5.8, 4.0, 1.0)
  )

  graphics::plot(
    NA_real_,
    NA_real_,
    type = "n",
    xlim = xlim,
    ylim = c(0, max(1, length(times))),
    axes = FALSE,
    ann = FALSE,
    bty = "n",
    xaxs = "i",
    yaxs = "i"
  )

  graphics::axis(
    side = 1,
    cex.axis = cex_axis,
    lwd = 1,
    lwd.ticks = 1
  )

  graphics::axis(
    side = 2,
    cex.axis = cex_axis,
    las = 1,
    lwd = 1,
    lwd.ticks = 1
  )

  graphics::mtext(
    text = .pp_time_description(x),
    side = 1,
    line = 3.4,
    cex = cex_lab
  )

  graphics::mtext(
    text = "Cumulative number of events",
    side = 2,
    line = 3.8,
    cex = cex_lab
  )

  graphics::title(
    main = main,
    cex.main = cex_main,
    font.main = 2
  )

  if (isTRUE(show_window)) {
    graphics::segments(
      x0 = xlim[1],
      y0 = 0,
      x1 = xlim[2],
      y1 = 0,
      col = "grey35",
      lwd = 1.3
    )

    graphics::segments(
      x0 = xlim[1],
      y0 = 0,
      x1 = xlim[1],
      y1 = length(times),
      col = "grey65",
      lty = 2
    )

    graphics::segments(
      x0 = xlim[2],
      y0 = 0,
      x1 = xlim[2],
      y1 = length(times),
      col = "grey65",
      lty = 2
    )
  }

  if (length(times) > 0L) {
    graphics::lines(
      x = c(xlim[1], times),
      y = c(0, seq_along(times)),
      type = "s",
      lwd = 2.2
    )
  }

  invisible(NULL)
}


.pp_plot_3d_panel <- function(
    coords,
    bounds,
    axis_labels,
    main,
    colours = "black",
    point_size = 1.35,
    point_pch = 19,
    theta = -45,
    phi = 20,
    show_window = TRUE,
    cex_axis = 1.15,
    cex_lab = 1.3,
    cex_main = 1.35
) {
  xlim <- .pp_plot_limits(bounds[, 1])
  ylim <- .pp_plot_limits(bounds[, 2])
  zlim <- .pp_plot_limits(bounds[, 3])

  empty_pattern <- nrow(coords) == 0L

  if (empty_pattern) {
    coords <- matrix(
      c(mean(xlim), mean(ylim), mean(zlim)),
      nrow = 1
    )

    colours <- NA
  }

  graphics::par(
    mar = c(4.0, 3.8, 3.8, 2.0)
  )

  plot3D::scatter3D(
    x = coords[, 1],
    y = coords[, 2],
    z = coords[, 3],
    xlim = xlim,
    ylim = ylim,
    zlim = zlim,
    xlab = axis_labels[1],
    ylab = axis_labels[2],
    zlab = axis_labels[3],
    main = main,
    col = colours,
    pch = if (empty_pattern) NA else point_pch,
    cex = if (empty_pattern) 0 else point_size,
    theta = theta,
    phi = phi,
    ticktype = "detailed",
    bty = if (show_window) "b" else "n",
    colkey = FALSE,
    cex.axis = cex_axis,
    cex.lab = cex_lab,
    cex.main = cex_main
  )

  invisible(NULL)
}


.pp_plot_joint_mark_panel <- function(
    x,
    spatial_axes,
    mark_name,
    mark_spec,
    point_size,
    point_pch,
    theta,
    phi,
    show_window,
    cex_axis,
    cex_lab,
    cex_main
) {
  xy_bounds <- .pp_plot_bounds(
    x = x,
    axes = spatial_axes
  )

  zlim <- mark_spec$limits

  if (is.null(zlim)) {
    zlim <- c(0, 1)
  }

  z_values <- mark_spec$numeric

  z_values[!is.finite(z_values)] <- zlim[1]

  coords <- cbind(
    x$coords[, spatial_axes, drop = FALSE],
    z_values
  )

  bounds <- cbind(
    xy_bounds,
    mark = zlim
  )

  colnames(bounds) <- c(spatial_axes, mark_name)

  .pp_plot_3d_panel(
    coords = coords,
    bounds = bounds,
    axis_labels = c(
      .pp_plot_axis_labels(x, spatial_axes),
      mark_name
    ),
    main = paste0(
      "Location and ",
      mark_name
    ),
    colours = mark_spec$colours,
    point_size = point_size,
    point_pch = point_pch,
    theta = theta,
    phi = phi,
    show_window = show_window,
    cex_axis = cex_axis,
    cex_lab = cex_lab,
    cex_main = cex_main
  )
}


#' @rdname pp-methods
#' @method plot pp
#' @export
plot.pp <- function(
    x,
    mark = NULL,
    axes = NULL,
    palette_continuous = "YlOrRd",
    palette_discrete = "Dark 3",
    missing_col = "grey75",
    point_col = "black",
    point_size = 1.35,
    point_pch = 19,
    theta = -45,
    phi = 20,
    show_window = TRUE,
    cex_axis = 1.3,
    cex_lab = 1.55,
    cex_main = 1.45,
    cex_legend = 1.2,
    colorbar_ticks = 5L,
    colorbar_digits = NULL,
    ...
) {
  if (!inherits(x, "pp")) {
    stop("`x` must inherit from class 'pp'.", call. = FALSE)
  }

  info <- .pp_pattern_info(x)

  spatial_axes <- .pp_plot_resolve_spatial_axes(
    x = x,
    axes = axes
  )

  n_spatial <- length(spatial_axes)

  time_axis <- x$axis_info$name[
    x$axis_info$role == "time"
  ]

  mark_names <- .pp_plot_resolve_marks(
    x = x,
    mark = mark
  )

  old_par <- graphics::par(no.readonly = TRUE)

  on.exit({
    graphics::layout(matrix(1L, 1L, 1L))
    graphics::par(old_par)
  }, add = TRUE)

  make_mark_spec <- function(mark_name) {
    .pp_plot_mark_spec(
      values = x$marks[[mark_name]],
      name = mark_name,
      palette_continuous = palette_continuous,
      palette_discrete = palette_discrete,
      missing_col = missing_col,
      colorbar_ticks = colorbar_ticks,
      colorbar_digits = colorbar_digits
    )
  }

  # -------------------------------------------------------------------------
  # Purely spatial patterns
  # -------------------------------------------------------------------------

  if (!info$is_st) {
    if (n_spatial == 1L) {
      axis_name <- spatial_axes[1]
      bounds <- .pp_plot_bounds(x, axis_name)

      graphics::layout(matrix(1L, 1L, 1L))

      graphics::par(mar = c(5.5, 2.0, 4.0, 1.0))

      graphics::plot(
        NA_real_,
        NA_real_,
        type = "n",
        xlim = .pp_plot_limits(bounds[, 1]),
        ylim = c(0.5, 1.5),
        axes = FALSE,
        ann = FALSE,
        bty = "n",
        xaxs = "i",
        yaxs = "i"
      )

      graphics::axis(
        1,
        cex.axis = cex_axis
      )

      graphics::mtext(
        .pp_plot_axis_labels(x, axis_name),
        side = 1,
        line = 3.4,
        cex = cex_lab
      )

      graphics::title(
        main = "One-dimensional spatial point pattern",
        cex.main = cex_main,
        font.main = 2
      )

      if (show_window) {
        graphics::segments(
          x0 = bounds[1, 1],
          y0 = 1,
          x1 = bounds[2, 1],
          y1 = 1,
          lwd = 1.3,
          col = "grey35"
        )
      }

      graphics::points(
        x$coords[, axis_name],
        rep(1, nrow(x$coords)),
        pch = point_pch,
        cex = point_size,
        col = point_col
      )

      return(invisible(x))
    }

    if (n_spatial == 2L) {
      bounds <- .pp_plot_bounds(x, spatial_axes)
      labels <- .pp_plot_axis_labels(x, spatial_axes)
      coords <- x$coords[, spatial_axes, drop = FALSE]

      if (!info$is_marked) {
        graphics::layout(matrix(1L, 1L, 1L))

        .pp_plot_2d_panel(
          coords = coords,
          bounds = bounds,
          axis_labels = labels,
          main = "Spatial point pattern",
          colours = point_col,
          point_size = point_size,
          point_pch = point_pch,
          show_window = show_window,
          cex_axis = cex_axis,
          cex_lab = cex_lab,
          cex_main = cex_main
        )

        return(invisible(x))
      }

      specs <- lapply(mark_names, make_mark_spec)

      any_continuous <- any(
        vapply(
          specs,
          function(z) identical(z$type, "continuous"),
          logical(1)
        )
      )

      n_cols <- if (any_continuous) 4L else 3L

      layout_mat <- matrix(
        seq_len(length(mark_names) * n_cols),
        nrow = length(mark_names),
        byrow = TRUE
      )

      graphics::layout(
        mat = layout_mat,
        widths = if (any_continuous) {
          c(4.5, 4.5, 1.6, 5.2)
        } else {
          c(4.8, 4.8, 2.3)
        }
      )

      for (i in seq_along(mark_names)) {
        mark_name <- mark_names[i]
        spec <- specs[[i]]

        .pp_plot_2d_panel(
          coords = coords,
          bounds = bounds,
          axis_labels = labels,
          main = "Spatial point pattern",
          colours = point_col,
          point_size = point_size,
          point_pch = point_pch,
          show_window = show_window,
          cex_axis = cex_axis,
          cex_lab = cex_lab,
          cex_main = cex_main
        )

        .pp_plot_2d_panel(
          coords = coords,
          bounds = bounds,
          axis_labels = labels,
          main = paste0("Coloured by ", mark_name),
          colours = spec$colours,
          point_size = point_size,
          point_pch = point_pch,
          show_window = show_window,
          cex_axis = cex_axis,
          cex_lab = cex_lab,
          cex_main = cex_main
        )

        if (identical(spec$type, "continuous")) {
          .pp_plot_colourbar_panel(
            spec = spec,
            cex_axis = cex_legend,
            cex_main = cex_main
          )

          .pp_plot_joint_mark_panel(
            x = x,
            spatial_axes = spatial_axes,
            mark_name = mark_name,
            mark_spec = spec,
            point_size = point_size,
            point_pch = point_pch,
            theta = theta,
            phi = phi,
            show_window = show_window,
            cex_axis = cex_axis,
            cex_lab = cex_lab,
            cex_main = cex_main
          )
        } else {
          .pp_plot_legend_panel(
            spec = spec,
            point_pch = point_pch,
            cex_legend = cex_legend,
            cex_main = cex_main
          )

          if (any_continuous) {
            graphics::plot.new()
          }
        }
      }

      return(invisible(x))
    }

    bounds <- .pp_plot_bounds(x, spatial_axes)
    labels <- .pp_plot_axis_labels(x, spatial_axes)
    coords <- x$coords[, spatial_axes, drop = FALSE]

    if (!info$is_marked) {
      graphics::layout(matrix(1L, 1L, 1L))

      .pp_plot_3d_panel(
        coords = coords,
        bounds = bounds,
        axis_labels = labels,
        main = "Three-dimensional spatial point pattern",
        colours = point_col,
        point_size = point_size,
        point_pch = point_pch,
        theta = theta,
        phi = phi,
        show_window = show_window,
        cex_axis = cex_axis,
        cex_lab = cex_lab,
        cex_main = cex_main
      )

      return(invisible(x))
    }

    graphics::layout(
      matrix(
        seq_len(length(mark_names) * 2L),
        nrow = length(mark_names),
        byrow = TRUE
      ),
      widths = c(5.8, 2.1)
    )

    for (mark_name in mark_names) {
      spec <- make_mark_spec(mark_name)

      .pp_plot_3d_panel(
        coords = coords,
        bounds = bounds,
        axis_labels = labels,
        main = paste0(
          "Three-dimensional spatial pattern coloured by ",
          mark_name
        ),
        colours = spec$colours,
        point_size = point_size,
        point_pch = point_pch,
        theta = theta,
        phi = phi,
        show_window = show_window,
        cex_axis = cex_axis,
        cex_lab = cex_lab,
        cex_main = cex_main
      )

      if (identical(spec$type, "continuous")) {
        .pp_plot_colourbar_panel(
          spec = spec,
          cex_axis = cex_legend,
          cex_main = cex_main
        )
      } else {
        .pp_plot_legend_panel(
          spec = spec,
          point_pch = point_pch,
          cex_legend = cex_legend,
          cex_main = cex_main
        )
      }
    }

    return(invisible(x))
  }

  # -------------------------------------------------------------------------
  # Spatio-temporal patterns
  # -------------------------------------------------------------------------

  if (length(time_axis) != 1L) {
    stop(
      "A spatio-temporal pattern must contain exactly one time axis.",
      call. = FALSE
    )
  }

  if (n_spatial == 1L) {
    st_axes <- c(spatial_axes, time_axis)
    bounds <- .pp_plot_bounds(x, st_axes)
    labels <- .pp_plot_axis_labels(x, st_axes)
    coords <- x$coords[, st_axes, drop = FALSE]

    if (!info$is_marked) {
      graphics::layout(
        matrix(1:2, nrow = 1),
        widths = c(5.0, 5.0)
      )

      .pp_plot_2d_panel(
        coords = coords,
        bounds = bounds,
        axis_labels = labels,
        main = "Space-time point pattern",
        colours = point_col,
        point_size = point_size,
        point_pch = point_pch,
        equal_axes = FALSE,
        show_window = show_window,
        cex_axis = cex_axis,
        cex_lab = cex_lab,
        cex_main = cex_main
      )

      .pp_plot_time_panel(
        x = x,
        time_axis = time_axis,
        show_window = show_window,
        cex_axis = cex_axis,
        cex_lab = cex_lab,
        cex_main = cex_main
      )

      return(invisible(x))
    }

    graphics::layout(
      matrix(
        seq_len(length(mark_names) * 3L),
        nrow = length(mark_names),
        byrow = TRUE
      ),
      widths = c(5.0, 2.0, 5.0)
    )

    for (mark_name in mark_names) {
      spec <- make_mark_spec(mark_name)

      .pp_plot_2d_panel(
        coords = coords,
        bounds = bounds,
        axis_labels = labels,
        main = paste0(
          "Space-time pattern coloured by ",
          mark_name
        ),
        colours = spec$colours,
        point_size = point_size,
        point_pch = point_pch,
        equal_axes = FALSE,
        show_window = show_window,
        cex_axis = cex_axis,
        cex_lab = cex_lab,
        cex_main = cex_main
      )

      if (identical(spec$type, "continuous")) {
        .pp_plot_colourbar_panel(
          spec = spec,
          cex_axis = cex_legend,
          cex_main = cex_main
        )
      } else {
        .pp_plot_legend_panel(
          spec = spec,
          point_pch = point_pch,
          cex_legend = cex_legend,
          cex_main = cex_main
        )
      }

      .pp_plot_time_panel(
        x = x,
        time_axis = time_axis,
        show_window = show_window,
        cex_axis = cex_axis,
        cex_lab = cex_lab,
        cex_main = cex_main
      )
    }

    return(invisible(x))
  }

  if (n_spatial == 2L) {
    spatial_bounds <- .pp_plot_bounds(x, spatial_axes)
    spatial_labels <- .pp_plot_axis_labels(x, spatial_axes)
    spatial_coords <- x$coords[, spatial_axes, drop = FALSE]

    st_axes <- c(spatial_axes, time_axis)

    st_bounds <- .pp_plot_bounds(x, st_axes)
    st_labels <- .pp_plot_axis_labels(x, st_axes)
    st_coords <- x$coords[, st_axes, drop = FALSE]

    if (!info$is_marked) {
      graphics::layout(
        matrix(1:3, nrow = 1),
        widths = c(4.8, 4.8, 5.5)
      )

      .pp_plot_2d_panel(
        coords = spatial_coords,
        bounds = spatial_bounds,
        axis_labels = spatial_labels,
        main = "Spatial point pattern",
        colours = point_col,
        point_size = point_size,
        point_pch = point_pch,
        show_window = show_window,
        cex_axis = cex_axis,
        cex_lab = cex_lab,
        cex_main = cex_main
      )

      .pp_plot_time_panel(
        x = x,
        time_axis = time_axis,
        show_window = show_window,
        cex_axis = cex_axis,
        cex_lab = cex_lab,
        cex_main = cex_main
      )

      .pp_plot_3d_panel(
        coords = st_coords,
        bounds = st_bounds,
        axis_labels = st_labels,
        main = "Space-time point pattern",
        colours = point_col,
        point_size = point_size,
        point_pch = point_pch,
        theta = theta,
        phi = phi,
        show_window = show_window,
        cex_axis = cex_axis,
        cex_lab = cex_lab,
        cex_main = cex_main
      )

      return(invisible(x))
    }

    graphics::layout(
      matrix(
        seq_len(length(mark_names) * 4L),
        nrow = length(mark_names),
        byrow = TRUE
      ),
      widths = c(4.6, 1.8, 4.5, 5.4)
    )

    for (mark_name in mark_names) {
      spec <- make_mark_spec(mark_name)

      .pp_plot_2d_panel(
        coords = spatial_coords,
        bounds = spatial_bounds,
        axis_labels = spatial_labels,
        main = paste0(
          "Spatial pattern coloured by ",
          mark_name
        ),
        colours = spec$colours,
        point_size = point_size,
        point_pch = point_pch,
        show_window = show_window,
        cex_axis = cex_axis,
        cex_lab = cex_lab,
        cex_main = cex_main
      )

      if (identical(spec$type, "continuous")) {
        .pp_plot_colourbar_panel(
          spec = spec,
          cex_axis = cex_legend,
          cex_main = cex_main
        )
      } else {
        .pp_plot_legend_panel(
          spec = spec,
          point_pch = point_pch,
          cex_legend = cex_legend,
          cex_main = cex_main
        )
      }

      .pp_plot_time_panel(
        x = x,
        time_axis = time_axis,
        show_window = show_window,
        cex_axis = cex_axis,
        cex_lab = cex_lab,
        cex_main = cex_main
      )

      .pp_plot_3d_panel(
        coords = st_coords,
        bounds = st_bounds,
        axis_labels = st_labels,
        main = paste0(
          "Space-time pattern coloured by ",
          mark_name
        ),
        colours = spec$colours,
        point_size = point_size,
        point_pch = point_pch,
        theta = theta,
        phi = phi,
        show_window = show_window,
        cex_axis = cex_axis,
        cex_lab = cex_lab,
        cex_main = cex_main
      )
    }

    return(invisible(x))
  }

  # Three spatial axes plus time --------------------------------------------

  xyz_axes <- spatial_axes
  xyt_axes <- c(spatial_axes[1:2], time_axis)

  xyz_bounds <- .pp_plot_bounds(x, xyz_axes)
  xyt_bounds <- .pp_plot_bounds(x, xyt_axes)

  xyz_labels <- .pp_plot_axis_labels(x, xyz_axes)
  xyt_labels <- .pp_plot_axis_labels(x, xyt_axes)

  xyz_coords <- x$coords[, xyz_axes, drop = FALSE]
  xyt_coords <- x$coords[, xyt_axes, drop = FALSE]

  time_spec <- .pp_plot_mark_spec(
    values = x$coords[, time_axis],
    name = .pp_time_description(x),
    palette_continuous = palette_continuous,
    palette_discrete = palette_discrete,
    missing_col = missing_col,
    colorbar_ticks = colorbar_ticks,
    colorbar_digits = colorbar_digits
  )

  z_axis <- xyz_axes[3]

  z_spec <- .pp_plot_mark_spec(
    values = x$coords[, z_axis],
    name = .pp_plot_axis_labels(x, z_axis),
    palette_continuous = palette_continuous,
    palette_discrete = palette_discrete,
    missing_col = missing_col,
    colorbar_ticks = colorbar_ticks,
    colorbar_digits = colorbar_digits
  )

  if (!info$is_marked) {
    graphics::layout(
      matrix(1:4, nrow = 1),
      widths = c(5.0, 1.7, 5.0, 1.7)
    )

    .pp_plot_3d_panel(
      coords = xyz_coords,
      bounds = xyz_bounds,
      axis_labels = xyz_labels,
      main = "Spatial 3-D projection coloured by time",
      colours = time_spec$colours,
      point_size = point_size,
      point_pch = point_pch,
      theta = theta,
      phi = phi,
      show_window = show_window,
      cex_axis = cex_axis,
      cex_lab = cex_lab,
      cex_main = cex_main
    )

    .pp_plot_colourbar_panel(
      spec = time_spec,
      cex_axis = cex_legend,
      cex_main = cex_main
    )

    .pp_plot_3d_panel(
      coords = xyt_coords,
      bounds = xyt_bounds,
      axis_labels = xyt_labels,
      main = paste0("Space-time projection coloured by ", z_axis),
      colours = z_spec$colours,
      point_size = point_size,
      point_pch = point_pch,
      theta = theta,
      phi = phi,
      show_window = show_window,
      cex_axis = cex_axis,
      cex_lab = cex_lab,
      cex_main = cex_main
    )

    .pp_plot_colourbar_panel(
      spec = z_spec,
      cex_axis = cex_legend,
      cex_main = cex_main
    )

    return(invisible(x))
  }

  layout_mat <- matrix(
    0L,
    nrow = 2L * length(mark_names),
    ncol = 4L
  )

  panel_id <- 1L

  for (i in seq_along(mark_names)) {
    row_top <- 2L * i - 1L
    row_bottom <- 2L * i

    layout_mat[row_top, ] <- panel_id:(panel_id + 3L)
    panel_id <- panel_id + 4L

    layout_mat[row_bottom, 1:2] <- panel_id:(panel_id + 1L)
    panel_id <- panel_id + 2L
  }

  graphics::layout(
    mat = layout_mat,
    widths = c(5.0, 1.7, 5.0, 1.7)
  )

  for (mark_name in mark_names) {
    spec <- make_mark_spec(mark_name)

    .pp_plot_3d_panel(
      coords = xyz_coords,
      bounds = xyz_bounds,
      axis_labels = xyz_labels,
      main = "Spatial 3-D projection coloured by time",
      colours = time_spec$colours,
      point_size = point_size,
      point_pch = point_pch,
      theta = theta,
      phi = phi,
      show_window = show_window,
      cex_axis = cex_axis,
      cex_lab = cex_lab,
      cex_main = cex_main
    )

    .pp_plot_colourbar_panel(
      spec = time_spec,
      cex_axis = cex_legend,
      cex_main = cex_main
    )

    .pp_plot_3d_panel(
      coords = xyt_coords,
      bounds = xyt_bounds,
      axis_labels = xyt_labels,
      main = paste0("Space-time projection coloured by ", z_axis),
      colours = z_spec$colours,
      point_size = point_size,
      point_pch = point_pch,
      theta = theta,
      phi = phi,
      show_window = show_window,
      cex_axis = cex_axis,
      cex_lab = cex_lab,
      cex_main = cex_main
    )

    .pp_plot_colourbar_panel(
      spec = z_spec,
      cex_axis = cex_legend,
      cex_main = cex_main
    )

    .pp_plot_3d_panel(
      coords = xyz_coords,
      bounds = xyz_bounds,
      axis_labels = xyz_labels,
      main = paste0(
        "Spatial 3-D projection coloured by ",
        mark_name
      ),
      colours = spec$colours,
      point_size = point_size,
      point_pch = point_pch,
      theta = theta,
      phi = phi,
      show_window = show_window,
      cex_axis = cex_axis,
      cex_lab = cex_lab,
      cex_main = cex_main
    )

    if (identical(spec$type, "continuous")) {
      .pp_plot_colourbar_panel(
        spec = spec,
        cex_axis = cex_legend,
        cex_main = cex_main
      )
    } else {
      .pp_plot_legend_panel(
        spec = spec,
        point_pch = point_pch,
        cex_legend = cex_legend,
        cex_main = cex_main
      )
    }
  }

  invisible(x)
}



# Covariate helpers ---------------------------------------------------------

.pp_cov_validate_name <- function(name, default) {
  if (is.null(name)) {
    return(default)
  }

  if (!is.character(name) ||
      length(name) != 1L ||
      is.na(name) ||
      !nzchar(name)) {
    stop(
      "`name` must be NULL or one non-empty character string.",
      call. = FALSE
    )
  }

  name
}


.pp_cov_validate_crs <- function(crs) {
  if (is.null(crs)) {
    return(NULL)
  }

  if (!is.character(crs) ||
      length(crs) != 1L ||
      is.na(crs) ||
      !nzchar(crs)) {
    stop(
      "`crs` must be NULL or one non-empty character string.",
      call. = FALSE
    )
  }

  crs
}


.pp_cov_observed_bounds <- function(coords) {
  if (nrow(coords) == 0L) {
    return(
      matrix(
        NA_real_,
        nrow = 2L,
        ncol = ncol(coords),
        dimnames = list(
          c("lower", "upper"),
          colnames(coords)
        )
      )
    )
  }

  out <- rbind(
    lower = apply(coords, 2L, min),
    upper = apply(coords, 2L, max)
  )

  storage.mode(out) <- "double"
  out
}


.pp_cov_align_numeric_argument <- function(x, axis_names, arg_name) {
  if (is.null(x)) {
    return(NULL)
  }

  if (!is.numeric(x) || any(!is.finite(x))) {
    stop(
      sprintf("`%s` must contain finite numeric values.", arg_name),
      call. = FALSE
    )
  }

  if (!is.null(names(x))) {
    if (anyDuplicated(names(x))) {
      stop(
        sprintf("`%s` cannot contain duplicated axis names.", arg_name),
        call. = FALSE
      )
    }

    missing_axes <- setdiff(axis_names, names(x))
    extra_axes <- setdiff(names(x), axis_names)

    if (length(missing_axes) > 0L || length(extra_axes) > 0L) {
      stop(
        sprintf(
          "`%s` must be named exactly with these axes: %s.",
          arg_name,
          paste(axis_names, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    x <- x[axis_names]
  } else if (length(x) != length(axis_names)) {
    stop(
      sprintf(
        "`%s` must have length %d or be named using the axis names.",
        arg_name,
        length(axis_names)
      ),
      call. = FALSE
    )
  }

  unname(as.numeric(x))
}


.pp_cov_validate_axis_scale <- function(covariate, axis_scale) {
  axis_names <- covariate$axis_info$name

  if (is.null(axis_scale)) {
    has_space <- any(covariate$axis_info$role == "space")
    has_time <- any(covariate$axis_info$role == "time")

    if (has_space && has_time) {
      warning(
        paste0(
          "The covariate combines spatial and temporal axes, but ",
          "`axis_scale` was not supplied. Distances are therefore computed ",
          "with scale 1 on every axis. Consider specifying `axis_scale`, ",
          "especially for the temporal axis."
        ),
        call. = FALSE
      )
    }

    return(rep(1, length(axis_names)))
  }

  axis_scale <- .pp_cov_align_numeric_argument(
    x = axis_scale,
    axis_names = axis_names,
    arg_name = "axis_scale"
  )

  if (any(axis_scale <= 0)) {
    stop(
      "`axis_scale` must contain strictly positive values.",
      call. = FALSE
    )
  }

  axis_scale
}


.pp_cov_validate_object <- function(x) {
  if (!inherits(x, "pp_covariate")) {
    stop(
      "`covariate` must inherit from class 'pp_covariate'.",
      call. = FALSE
    )
  }

  required <- c("coords", "axis_info", "window", "crs", "name")

  missing_components <- setdiff(required, names(x))

  if (length(missing_components) > 0L) {
    stop(
      sprintf(
        "The covariate object is missing: %s.",
        paste(missing_components, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(x)
}


.pp_cov_is_field <- function(x) {
  inherits(x, "pp_covariate_field")
}


.pp_cov_is_support <- function(x) {
  inherits(x, "pp_covariate_support")
}


.pp_cov_compatible_units <- function(cov_units, target_units) {
  known <- !is.na(cov_units) &
    nzchar(cov_units) &
    !is.na(target_units) &
    nzchar(target_units)

  if (!any(known)) {
    return(TRUE)
  }

  all(cov_units[known] == target_units[known])
}


.pp_cov_check_target_metadata <- function(covariate, target) {
  cov_axes <- covariate$axis_info$name
  target_axes <- target$axis_info$name

  missing_axes <- setdiff(cov_axes, target_axes)

  if (length(missing_axes) > 0L) {
    stop(
      sprintf(
        "The target is missing covariate axis/axes: %s.",
        paste(missing_axes, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  target_axis_info <- target$axis_info[
    match(cov_axes, target_axes),
    ,
    drop = FALSE
  ]

  if (!identical(
    covariate$axis_info$role,
    target_axis_info$role
  )) {
    stop(
      "Covariate and target have incompatible axis roles.",
      call. = FALSE
    )
  }

  if (!.pp_cov_compatible_units(
    cov_units = covariate$axis_info$unit,
    target_units = target_axis_info$unit
  )) {
    stop(
      "Covariate and target have incompatible axis units.",
      call. = FALSE
    )
  }

  if (!is.null(covariate$crs) &&
      !is.null(target$crs) &&
      !identical(covariate$crs, target$crs)) {
    stop(
      "Covariate and target have incompatible CRS values.",
      call. = FALSE
    )
  }

  cov_time_idx <- which(covariate$axis_info$role == "time")

  if (length(cov_time_idx) == 1L &&
      !is.null(covariate$time_info) &&
      !is.null(target$time_info)) {
    if (!.pp_same_time_origin(
      covariate$time_info$origin,
      target$time_info$origin
    )) {
      stop(
        "Covariate and target have incompatible temporal origins.",
        call. = FALSE
      )
    }
  }

  invisible(NULL)
}


.pp_cov_target_from_object <- function(covariate, at) {
  .pp_cov_check_target_metadata(
    covariate = covariate,
    target = at
  )

  cov_axes <- covariate$axis_info$name

  list(
    coords = at$coords[, cov_axes, drop = FALSE],
    output_coords = as.data.frame(
      at$coords,
      check.names = FALSE
    )
  )
}


.pp_cov_target_from_data <- function(covariate, at) {
  if (!is.data.frame(at) && !is.matrix(at)) {
    stop(
      "`at` must be a pp object, pp_grid object, data frame, or matrix.",
      call. = FALSE
    )
  }

  at <- as.data.frame(at, stringsAsFactors = FALSE)

  if (is.null(names(at)) ||
      any(!nzchar(names(at))) ||
      anyDuplicated(names(at))) {
    stop(
      "A data-frame or matrix target must have unique non-empty column names.",
      call. = FALSE
    )
  }

  axis_names <- covariate$axis_info$name

  missing_axes <- setdiff(axis_names, names(at))

  if (length(missing_axes) > 0L) {
    stop(
      sprintf(
        "The target data are missing covariate axis/axes: %s.",
        paste(missing_axes, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  coords <- matrix(
    NA_real_,
    nrow = nrow(at),
    ncol = length(axis_names),
    dimnames = list(NULL, axis_names)
  )

  for (j in seq_along(axis_names)) {
    axis_name <- axis_names[j]
    role <- covariate$axis_info$role[j]
    values <- at[[axis_name]]

    if (role == "space") {
      if (!is.numeric(values)) {
        stop(
          sprintf(
            "Target spatial coordinate `%s` must be numeric.",
            axis_name
          ),
          call. = FALSE
        )
      }

      values <- as.numeric(values)

      if (any(!is.finite(values))) {
        stop(
          sprintf(
            "Target spatial coordinate `%s` contains missing or non-finite values.",
            axis_name
          ),
          call. = FALSE
        )
      }

      coords[, j] <- values
    } else {
      converted <- .pp_convert_time_coordinate(
        x = values,
        time_unit = covariate$time_info$unit,
        time_origin = covariate$time_info$origin,
        axis_name = axis_name
      )

      coords[, j] <- converted$values
    }
  }

  list(
    coords = coords,
    output_coords = at[, axis_names, drop = FALSE]
  )
}


.pp_cov_target_coordinates <- function(covariate, at) {
  if (inherits(at, "pp") || inherits(at, "pp_grid")) {
    return(
      .pp_cov_target_from_object(
        covariate = covariate,
        at = at
      )
    )
  }

  .pp_cov_target_from_data(
    covariate = covariate,
    at = at
  )
}


.pp_cov_distance_matrix <- function(target, source, axis_scale) {
  n_target <- nrow(target)
  n_source <- nrow(source)

  if (n_target == 0L || n_source == 0L) {
    return(
      matrix(
        numeric(0),
        nrow = n_target,
        ncol = n_source
      )
    )
  }

  distance_squared <- matrix(
    0,
    nrow = n_target,
    ncol = n_source
  )

  for (j in seq_len(ncol(target))) {
    difference <- outer(
      target[, j],
      source[, j],
      FUN = "-"
    ) / axis_scale[j]

    distance_squared <- distance_squared + difference^2
  }

  sqrt(distance_squared)
}


.pp_cov_nearest_values <- function(distance_matrix, values) {
  out <- numeric(nrow(distance_matrix))

  for (i in seq_len(nrow(distance_matrix))) {
    d <- distance_matrix[i, ]
    min_d <- min(d)

    tolerance <- .Machine$double.eps^0.5 * max(1, min_d)

    nearest <- which(abs(d - min_d) <= tolerance)

    out[i] <- mean(values[nearest])
  }

  out
}


.pp_cov_idw_values <- function(
    distance_matrix,
    values,
    p,
    k = NULL
) {
  out <- numeric(nrow(distance_matrix))

  for (i in seq_len(nrow(distance_matrix))) {
    d <- distance_matrix[i, ]

    zero_tolerance <- .Machine$double.eps^0.5

    zero_distance <- which(d <= zero_tolerance)

    # Se il punto target coincide con uno o più punti di supporto,
    # restituisce la media dei valori nei punti coincidenti.
    if (length(zero_distance) > 0L) {
      out[i] <- mean(values[zero_distance])
      next
    }

    # Se richiesto, usa soltanto i k supporti più vicini.
    if (!is.null(k) && k < length(d)) {
      selected <- order(d)[seq_len(k)]

      d_selected <- d[selected]
      values_selected <- values[selected]
    } else {
      d_selected <- d
      values_selected <- values
    }

    weights <- d_selected^(-p)

    out[i] <- sum(weights * values_selected) / sum(weights)
  }

  out
}


.pp_cov_chunk_indices <- function(n, chunk_size) {
  starts <- seq.int(1L, n, by = chunk_size)

  lapply(
    starts,
    function(start) {
      start:min(n, start + chunk_size - 1L)
    }
  )
}


.pp_cov_default_method <- function(covariate) {
  if (.pp_cov_is_field(covariate)) {
    return("idw")
  }

  "distance_min"
}


.pp_cov_validate_method <- function(covariate, method) {
  allowed <- if (.pp_cov_is_field(covariate)) {
    c("idw", "nearest")
  } else {
    "distance_min"
  }

  if (is.null(method)) {
    return(.pp_cov_default_method(covariate))
  }

  if (!is.character(method) ||
      length(method) != 1L ||
      is.na(method)) {
    stop(
      "`method` must be one character string.",
      call. = FALSE
    )
  }

  method <- tolower(trimws(method))

  if (!method %in% allowed) {
    stop(
      sprintf(
        "Method `%s` is not available for this covariate. Allowed method(s): %s.",
        method,
        paste(allowed, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  method
}


.pp_cov_output_coordinates <- function(at, required_axes) {
  if (inherits(at, "pp") || inherits(at, "pp_grid")) {
    return(
      as.data.frame(
        at$coords,
        check.names = FALSE
      )
    )
  }

  at <- as.data.frame(at, stringsAsFactors = FALSE)

  missing_axes <- setdiff(required_axes, names(at))

  if (length(missing_axes) > 0L) {
    stop(
      sprintf(
        "The target data are missing coordinate axis/axes: %s.",
        paste(missing_axes, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  at[, required_axes, drop = FALSE]
}


#' Create a multidimensional covariate object
#'
#' Creates either a continuous numerical covariate field or a geometric support
#' over spatial or spatio-temporal coordinates. Covariate objects can be
#' evaluated at event locations, regular grids, or arbitrary target points with
#' [evaluate_covariate()].
#'
#' @details
#' If `value` identifies a numeric column, the result is a
#' `"pp_covariate_field"`: each support location carries one observed continuous
#' value. Such a field can be evaluated by inverse-distance weighting or
#' nearest-neighbour assignment.
#'
#' If `value = NULL`, the result is a `"pp_covariate_support"`: only support
#' coordinates are retained. It represents a geometric feature such as fault
#' locations, facilities, stations, or other sources, and is evaluated as the
#' distance to the nearest support location.
#'
#' The covariate may use all or a subset of the axes of `window`. When a window
#' is supplied, coordinate order and metadata are aligned with the corresponding
#' window axes, and every support location must lie within the relevant bounds.
#' A spatio-temporal covariate may contain at most one temporal axis and must
#' contain at least one spatial axis.
#'
#' @param data A data frame or matrix containing support coordinates and,
#'   optionally, one continuous numeric covariate column. Coordinate values must
#'   be finite; a temporal coordinate may be numeric, `Date`, or
#'   `POSIXct`/`POSIXt`.
#' @param coords Coordinate-column specification as a character vector of names
#'   or numeric vector of positions. It is required unless `d` is supplied and
#'   cannot be used together with `d`.
#' @param d Optional positive integer indicating that the first `d` columns of
#'   `data` are coordinates. Cannot be supplied together with `coords`.
#' @param value Optional name or one-based position of exactly one continuous
#'   numeric covariate column. Use `NULL` to create a geometric support without
#'   observed values.
#' @param axis_roles Optional character vector containing `"space"` and at most
#'   one `"time"` entry. It may be unnamed, in coordinate order, or named with
#'   the coordinate names. Roles are inherited from `window` when omitted;
#'   otherwise all coordinates are spatial.
#' @param axis_units Optional character vector of coordinate units, unnamed or
#'   named with the coordinate axes. Units are inherited from `window` when
#'   omitted.
#' @param time_unit Optional temporal unit, one of `"days"`, `"hours"`,
#'   `"mins"`, `"secs"`, or `"weeks"`. It is inherited from `window` or the
#'   temporal entry of `axis_units` when omitted, and otherwise defaults to
#'   `"days"`.
#' @param time_origin Optional temporal origin. Numeric, `Date`, and
#'   `POSIXct`/`POSIXt` origins are supported according to the class of the
#'   temporal coordinate. It is inherited from `window` when available; without
#'   a window, numeric time defaults to origin zero and date-time input defaults
#'   to its minimum observed value.
#' @param window Optional [pp_window()] defining a reference domain and axis
#'   metadata. The covariate may use a subset of its axes, but matching roles,
#'   units, temporal origin, and coordinate reference system must be compatible.
#' @param crs Optional non-empty character string describing the coordinate
#'   reference system. It is inherited from `window` when available and must
#'   agree with it when supplied explicitly.
#' @param name Optional non-empty covariate name. By default, a field uses the
#'   name of its value column and a geometric support uses `"support"`.
#'
#' @return An object inheriting from `"pp_covariate"` and from either
#' `"pp_covariate_field"` or `"pp_covariate_support"`. Its principal components
#' are `coords`, `values`, `axis_info`, `time_info`, `window`, `crs`, `name`, and
#' `metadata`. The `values` component is `NULL` for a geometric support.
#'
#' @examples
#' w <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#'
#' # A continuous field observed on a regular support
#' field_data <- expand.grid(
#'   x = seq(0, 1, length.out = 5),
#'   y = seq(0, 1, length.out = 5)
#' )
#' field_data$cov1 <- with(field_data, 15 + 3 * x - 2 * y)
#'
#' cov1 <- pp_covariate(
#'   data = field_data,
#'   coords = c("x", "y"),
#'   value = "cov1",
#'   window = w,
#'   name = "cov1"
#' )
#' cov1
#'
#' # A geometric support used to calculate nearest-source distances
#' stations <- pp_covariate(
#'   data = data.frame(
#'     x = c(0.2, 0.8, 0.5),
#'     y = c(0.2, 0.3, 0.9)
#'   ),
#'   coords = c("x", "y"),
#'   value = NULL,
#'   window = w,
#'   name = "station_distance"
#' )
#'
#' # A space-time covariate with dates
#' origin <- as.Date("2024-01-01")
#' wst <- pp_window(
#'   lower = c(x = 0, time = 0),
#'   upper = c(x = 1, time = 10),
#'   axis_roles = c(x = "space", time = "time"),
#'   axis_units = c(x = "km", time = "days"),
#'   time_origin = origin
#' )
#' st_field <- pp_covariate(
#'   data = data.frame(
#'     x = rep(c(0, 0.5, 1), each = 3),
#'     time = origin + rep(c(0, 5, 10), times = 3),
#'     value = seq_len(9)
#'   ),
#'   coords = c("x", "time"),
#'   value = "value",
#'   window = wst,
#'   name = "space_time_field"
#' )
#'
#' @seealso
#' [evaluate_covariate()] and [extract_covariates()] for evaluation,
#' [pp_grid()] for regular targets, and [pp_window()] for reference domains.
#' Covariate objects can be visualised with `plot()`.
#'
#' @family covariate tools
#' @export
pp_covariate <- function(
    data,
    coords = NULL,
    d = NULL,
    value = NULL,
    axis_roles = NULL,
    axis_units = NULL,
    time_unit = NULL,
    time_origin = NULL,
    window = NULL,
    crs = NULL,
    name = NULL
) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data frame or matrix.", call. = FALSE)
  }

  data <- as.data.frame(data, stringsAsFactors = FALSE)

  if (ncol(data) == 0L) {
    stop("`data` must contain at least one column.", call. = FALSE)
  }

  if (is.null(names(data)) || any(!nzchar(names(data)))) {
    names(data) <- paste0("V", seq_len(ncol(data)))
  }

  if (anyDuplicated(names(data))) {
    stop("`data` must have unique column names.", call. = FALSE)
  }

  if (!is.null(window) && !inherits(window, "pp_window")) {
    stop(
      "`window` must be NULL or an object of class 'pp_window'.",
      call. = FALSE
    )
  }

  if (!is.null(coords) && !is.null(d)) {
    stop("Supply either `coords` or `d`, not both.", call. = FALSE)
  }

  if (is.null(coords) && is.null(d)) {
    stop(
      "Specify covariate coordinate columns through `coords` or `d`.",
      call. = FALSE
    )
  }

  if (!is.null(d)) {
    if (!is.numeric(d) ||
        length(d) != 1L ||
        !is.finite(d) ||
        d != as.integer(d) ||
        d < 1L ||
        d > ncol(data)) {
      stop(
        "`d` must be a positive integer not larger than the number of data columns.",
        call. = FALSE
      )
    }

    coord_idx <- seq_len(as.integer(d))
  } else {
    coord_idx <- .pp_resolve_columns(
      data = data,
      columns = coords,
      arg_name = "coords"
    )
  }

  if (anyDuplicated(coord_idx)) {
    stop("`coords` cannot contain duplicated columns.", call. = FALSE)
  }

  coord_names <- names(data)[coord_idx]

  if (!is.null(window)) {
    missing_window_axes <- setdiff(coord_names, window$axis_info$name)

    if (length(missing_window_axes) > 0L) {
      stop(
        sprintf(
          "Covariate coordinate axis/axes are absent from `window`: %s.",
          paste(missing_window_axes, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    window_order <- window$axis_info$name[
      window$axis_info$name %in% coord_names
    ]

    coord_idx <- coord_idx[match(window_order, coord_names)]
    coord_names <- names(data)[coord_idx]
  }

  d <- length(coord_idx)

  axis_roles <- .pp_align_axis_argument(
    x = axis_roles,
    axis_names = coord_names,
    arg_name = "axis_roles"
  )

  if (is.null(axis_roles) && !is.null(window)) {
    axis_roles <- window$axis_info$role[
      match(coord_names, window$axis_info$name)
    ]
  }

  if (is.null(axis_roles)) {
    axis_roles <- rep("space", d)
  }

  axis_roles <- tolower(trimws(axis_roles))

  if (!all(axis_roles %in% c("space", "time"))) {
    stop(
      "`axis_roles` may contain only 'space' and 'time'.",
      call. = FALSE
    )
  }

  if (!is.null(window)) {
    window_roles <- window$axis_info$role[
      match(coord_names, window$axis_info$name)
    ]

    if (!identical(axis_roles, window_roles)) {
      stop(
        "`axis_roles` is not consistent with the roles stored in `window`.",
        call. = FALSE
      )
    }
  }

  n_time <- sum(axis_roles == "time")

  if (n_time > 1L) {
    stop(
      "A covariate may currently contain at most one temporal axis.",
      call. = FALSE
    )
  }

  if (n_time == 1L && sum(axis_roles == "space") == 0L) {
    stop(
      "A spatio-temporal covariate must contain at least one spatial axis.",
      call. = FALSE
    )
  }

  if (n_time == 0L &&
      (!is.null(time_unit) || !is.null(time_origin))) {
    stop(
      "`time_unit` and `time_origin` require one temporal covariate axis.",
      call. = FALSE
    )
  }

  axis_units <- .pp_align_axis_argument(
    x = axis_units,
    axis_names = coord_names,
    arg_name = "axis_units"
  )

  if (is.null(axis_units) && !is.null(window)) {
    axis_units <- window$axis_info$unit[
      match(coord_names, window$axis_info$name)
    ]
  }

  if (is.null(axis_units)) {
    axis_units <- rep(NA_character_, d)
  }

  if (n_time == 1L) {
    time_idx <- which(axis_roles == "time")
    window_time_info <- if (!is.null(window)) window$time_info else NULL

    if (is.null(time_unit) && !is.null(window_time_info)) {
      time_unit <- window_time_info$unit
    }

    if (is.null(time_unit) &&
        !is.na(axis_units[time_idx]) &&
        nzchar(axis_units[time_idx])) {
      time_unit <- axis_units[time_idx]
    }

    if (is.null(time_unit)) {
      time_unit <- "days"
    }

    time_unit <- .pp_validate_time_unit(time_unit)

    if (!is.na(axis_units[time_idx]) &&
        nzchar(axis_units[time_idx]) &&
        !identical(axis_units[time_idx], time_unit)) {
      stop(
        "`axis_units` for the time axis must agree with `time_unit`.",
        call. = FALSE
      )
    }

    axis_units[time_idx] <- time_unit

    if (!is.null(window_time_info)) {
      if (is.null(time_origin)) {
        time_origin <- window_time_info$origin
      } else if (!.pp_same_time_origin(
        time_origin,
        window_time_info$origin
      )) {
        stop(
          "`time_origin` is not consistent with the time origin stored in `window`.",
          call. = FALSE
        )
      }
    }
  }

  if (!is.null(window)) {
    window_units <- window$axis_info$unit[
      match(coord_names, window$axis_info$name)
    ]

    if (!identical(axis_units, window_units)) {
      stop(
        "`axis_units` is not consistent with the units stored in `window`.",
        call. = FALSE
      )
    }
  }

  crs <- .pp_cov_validate_crs(crs)

  if (!is.null(window)) {
    if (is.null(window$crs)) {
      if (!is.null(crs)) {
        stop(
          "The supplied `window` has no CRS. Define its CRS in `pp_window()` first.",
          call. = FALSE
        )
      }
    } else if (is.null(crs)) {
      crs <- window$crs
    } else if (!identical(crs, window$crs)) {
      stop(
        "`crs` is not consistent with the CRS stored in `window`.",
        call. = FALSE
      )
    }
  }

  coords_matrix <- matrix(
    NA_real_,
    nrow = nrow(data),
    ncol = d,
    dimnames = list(NULL, coord_names)
  )

  time_info <- NULL

  for (j in seq_len(d)) {
    variable <- data[[coord_idx[j]]]

    if (axis_roles[j] == "time") {
      converted <- .pp_convert_time_coordinate(
        x = variable,
        time_unit = time_unit,
        time_origin = time_origin,
        axis_name = coord_names[j]
      )

      coords_matrix[, j] <- converted$values

      time_info <- list(
        axis = coord_names[j],
        origin = converted$origin,
        original_class = converted$original_class,
        unit = converted$unit
      )
    } else {
      if (!is.numeric(variable)) {
        stop(
          sprintf(
            "Spatial coordinate `%s` must be numeric.",
            coord_names[j]
          ),
          call. = FALSE
        )
      }

      variable <- as.numeric(variable)

      if (any(!is.finite(variable))) {
        stop(
          sprintf(
            "Spatial coordinate `%s` contains missing or non-finite values.",
            coord_names[j]
          ),
          call. = FALSE
        )
      }

      coords_matrix[, j] <- variable
    }
  }

  if (!is.null(window) && nrow(coords_matrix) > 0L) {
    bounds <- window$bounds[, coord_names, drop = FALSE]

    outside <- sweep(
      coords_matrix,
      2L,
      bounds["lower", ],
      FUN = "<"
    ) |
      sweep(
        coords_matrix,
        2L,
        bounds["upper", ],
        FUN = ">"
      )

    if (any(outside)) {
      bad_rows <- which(rowSums(outside) > 0L)

      stop(
        sprintf(
          "Some covariate support locations fall outside `window`. First affected row(s): %s.",
          paste(utils::head(bad_rows, 5L), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  values <- NULL
  value_name <- NULL

  if (!is.null(value)) {
    value_idx <- .pp_resolve_columns(
      data = data,
      columns = value,
      arg_name = "value"
    )

    if (length(value_idx) != 1L) {
      stop(
        "`value` must identify exactly one numerical covariate column.",
        call. = FALSE
      )
    }

    if (value_idx %in% coord_idx) {
      stop(
        "`value` cannot also be a coordinate column.",
        call. = FALSE
      )
    }

    values <- data[[value_idx]]

    if (!is.numeric(values)) {
      stop(
        "`value` must be a continuous numerical variable.",
        call. = FALSE
      )
    }

    values <- as.numeric(values)

    if (any(!is.finite(values))) {
      stop(
        "`value` contains missing or non-finite values.",
        call. = FALSE
      )
    }

    value_name <- names(data)[value_idx]
  }

  default_name <- if (is.null(values)) {
    "support"
  } else {
    value_name
  }

  name <- .pp_cov_validate_name(
    name = name,
    default = default_name
  )

  axis_info <- data.frame(
    name = coord_names,
    role = axis_roles,
    unit = axis_units,
    stringsAsFactors = FALSE
  )

  ignored_columns <- setdiff(
    names(data),
    c(coord_names, value_name)
  )

  out <- list(
    coords = coords_matrix,
    values = values,
    axis_info = axis_info,
    time_info = time_info,
    window = window,
    crs = crs,
    name = name,
    metadata = list(
      original_column_names = names(data),
      ignored_columns = ignored_columns,
      value_column = value_name
    )
  )

  class(out) <- if (is.null(values)) {
    c("pp_covariate_support", "pp_covariate")
  } else {
    c("pp_covariate_field", "pp_covariate")
  }

  out
}


#' Methods for multidimensional covariate objects
#'
#' Print, summarise, and plot objects created by [pp_covariate()].
#'
#' @details
#' `print.pp_covariate()` gives a concise description of the covariate type,
#' dimensions, support size, axes, and value range when present.
#' `summary.pp_covariate()` returns a structured summary, which is formatted by
#' `print.summary_pp_covariate()`.
#'
#' `plot.pp_covariate()` may display the original support or first evaluate the
#' covariate at `at`. A regular two-dimensional target is displayed as a raster
#' by default; irregular or higher-dimensional targets are displayed as coloured
#' points. Three-dimensional displays use `plot3D::scatter3D()`. A point pattern
#' may be overlaid to compare event locations with the covariate surface.
#'
#' @param x A `pp_covariate` object for `print()` and `plot()`, or a
#'   `summary_pp_covariate` object for `print.summary_pp_covariate()`.
#' @param object A `pp_covariate` object to be summarised.
#' @param digits Number of significant digits used when printing numeric values.
#' @param at Optional target locations supplied as a `pp`, `pp_grid`, data frame,
#'   or matrix. If omitted, the original covariate support is displayed. If
#'   supplied, `x` is first evaluated at these locations.
#' @param pattern Optional `pp` object overlaid on the covariate display as a
#'   second set of plotting symbols. It must be compatible with the displayed
#'   axes.
#' @param axes Optional character or numeric vector selecting two or three axes.
#'   It is required when the displayed target has more than three dimensions.
#' @param method Evaluation method passed to [evaluate_covariate()] when `at` is
#'   supplied. `NULL` uses the default method for the covariate type.
#' @param p,k,axis_scale Inverse-distance and distance-scaling arguments passed
#'   to [evaluate_covariate()]. See that function for their interpretation.
#' @param display Display type, matched to one of `"auto"`, `"points"`, or
#'   `"raster"`. `"auto"` uses a raster for a regular two-dimensional target and
#'   points otherwise. Raster display is available only for a regular
#'   two-dimensional target.
#' @param show_support Logical value indicating whether the original covariate
#'   support locations should be overlaid on an evaluated display.
#' @param show_window Logical value indicating whether the associated
#'   observation-window boundary should be drawn when available.
#' @param palette Name of an HCL palette used for covariate values.
#' @param missing_col Colour assigned to missing evaluated values.
#' @param point_size,point_pch Size multiplier and plotting symbol for the main
#'   covariate locations.
#' @param support_size,support_pch Size multiplier and optional plotting symbol
#'   for the original support overlay. If `support_pch = NULL`, a suitable
#'   default symbol is selected.
#' @param pattern_size,pattern_pch Size multiplier and plotting symbol for the
#'   optional point-pattern overlay.
#' @param theta,phi Viewing angles, in degrees, for three-dimensional displays.
#' @param max_display Positive integer limiting the number of locations drawn in
#'   point-based displays. A representative subset is used when necessary.
#' @param cex_axis,cex_lab,cex_main,cex_legend Positive scale factors for axis
#'   annotations, axis labels, titles, and colour-bar text.
#' @param colorbar_ticks Approximate positive number of ticks displayed on the
#'   colour bar.
#' @param colorbar_digits Optional non-negative integer fixing the number of
#'   decimal places in colour-bar labels. If `NULL`, a suitable value is chosen
#'   automatically.
#' @param main Optional character string used as the plot title. `NULL` creates
#'   a title from the covariate name and evaluation method.
#' @param ... Additional arguments. They are currently ignored.
#'
#' @return The print and plot methods invisibly return `x`.
#' `summary.pp_covariate()` returns an object of class
#' `"summary_pp_covariate"`.
#'
#' @examples
#' w <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' dat <- expand.grid(
#'   x = seq(0, 1, length.out = 6),
#'   y = seq(0, 1, length.out = 6)
#' )
#' dat$value <- with(dat, x + y)
#'
#' covariate <- pp_covariate(
#'   dat,
#'   coords = c("x", "y"),
#'   value = "value",
#'   window = w,
#'   name = "linear_field"
#' )
#'
#' print(covariate)
#' sc <- summary(covariate)
#' print(sc)
#' plot(covariate)
#'
#' @seealso
#' [pp_covariate()], [evaluate_covariate()], [extract_covariates()], and
#' [pp_grid()].
#' @family covariate tools
#' @rdname pp_covariate-methods
#' @method print pp_covariate
#' @export
print.pp_covariate <- function(x, digits = 4, ...) {
  .pp_cov_validate_object(x)

  is_field <- .pp_cov_is_field(x)

  cat("<", class(x)[1], ">\n", sep = "")

  if (is_field) {
    cat("Continuous covariate field: ", x$name, "\n\n", sep = "")
  } else {
    cat("Geometric covariate support: ", x$name, "\n\n", sep = "")
  }

  cat("Support locations: ", nrow(x$coords), "\n", sep = "")
  cat("Dimensions: ", ncol(x$coords), "\n", sep = "")

  axis_labels <- vapply(
    seq_len(nrow(x$axis_info)),
    function(j) {
      .pp_axis_label(
        x$axis_info$name[j],
        x$axis_info$unit[j]
      )
    },
    character(1)
  )

  cat("Axes: ", paste(axis_labels, collapse = ", "), "\n", sep = "")

  if (any(x$axis_info$role == "time")) {
    time_idx <- which(x$axis_info$role == "time")

    time_label <- .pp_axis_label(
      x$axis_info$name[time_idx],
      x$axis_info$unit[time_idx]
    )

    cat("Temporal axis: ", time_label, "\n", sep = "")
  }

  if (!is.null(x$crs)) {
    cat("CRS: ", x$crs, "\n", sep = "")
  }

  if (!is.null(x$window)) {
    cat(
      "Reference window: supplied (",
      x$window$dimension,
      " dimensions)\n",
      sep = ""
    )
  } else {
    cat("Reference window: not supplied\n")
  }

  if (is_field && length(x$values) > 0L) {
    cat(
      "Value range: [",
      .pp_format_number(min(x$values), digits),
      ", ",
      .pp_format_number(max(x$values), digits),
      "]\n",
      sep = ""
    )
  }

  cat("\nSupport bounding box:\n")

  .pp_print_bounds(
    bounds = .pp_cov_observed_bounds(x$coords),
    axis_info = x$axis_info,
    digits = digits
  )

  invisible(x)
}


#' @rdname pp_covariate-methods
#' @method summary pp_covariate
#' @export
summary.pp_covariate <- function(object, ...) {
  .pp_cov_validate_object(object)

  value_summary <- NULL

  if (.pp_cov_is_field(object)) {
    values <- object$values

    if (length(values) == 0L) {
      value_summary <- c(
        min = NA_real_,
        q1 = NA_real_,
        median = NA_real_,
        mean = NA_real_,
        sd = NA_real_,
        q3 = NA_real_,
        max = NA_real_
      )
    } else {
      q <- stats::quantile(
        values,
        probs = c(0.25, 0.50, 0.75),
        names = FALSE
      )

      value_summary <- c(
        min = min(values),
        q1 = q[1],
        median = q[2],
        mean = mean(values),
        sd = stats::sd(values),
        q3 = q[3],
        max = max(values)
      )
    }
  }

  time_summary <- NULL

  if (any(object$axis_info$role == "time")) {
    time_idx <- which(object$axis_info$role == "time")
    time_values <- object$coords[, time_idx]

    restored_range <- NULL

    if (!is.null(object$time_info) &&
        length(time_values) > 0L) {
      restored_range <- .pp_restore_time(
        values = range(time_values),
        time_info = object$time_info
      )
    }

    time_summary <- list(
      axis = object$axis_info$name[time_idx],
      unit = object$axis_info$unit[time_idx],
      numeric_range = if (length(time_values) == 0L) {
        c(NA_real_, NA_real_)
      } else {
        range(time_values)
      },
      original_range = restored_range,
      origin = object$time_info$origin
    )
  }

  out <- list(
    subclass = class(object)[1],
    name = object$name,
    type = if (.pp_cov_is_field(object)) "field" else "support",
    n_support = nrow(object$coords),
    dimension = ncol(object$coords),
    axis_info = object$axis_info,
    coordinate_summary = .pp_coordinate_summary(
      coords = object$coords,
      axis_info = object$axis_info
    ),
    support_bounds = .pp_cov_observed_bounds(object$coords),
    values = value_summary,
    time = time_summary,
    has_window = !is.null(object$window),
    crs = object$crs
  )

  class(out) <- "summary_pp_covariate"

  out
}


#' @rdname pp_covariate-methods
#' @method print summary_pp_covariate
#' @export
print.summary_pp_covariate <- function(x, digits = 4, ...) {
  cat("Covariate: ", x$name, "\n", sep = "")
  cat("Type: ", x$type, "\n", sep = "")
  cat("Support locations: ", x$n_support, "\n", sep = "")
  cat("Dimensions: ", x$dimension, "\n", sep = "")

  if (!is.null(x$crs)) {
    cat("CRS: ", x$crs, "\n", sep = "")
  }

  cat(
    "Reference window: ",
    if (x$has_window) "supplied" else "not supplied",
    "\n\n",
    sep = ""
  )

  cat("Coordinate summaries:\n")

  coordinate_summary <- x$coordinate_summary
  numeric_columns <- c("min", "q1", "median", "mean", "sd", "q3", "max")

  for (column in numeric_columns) {
    coordinate_summary[[column]] <- vapply(
      coordinate_summary[[column]],
      .pp_format_number,
      character(1),
      digits = digits
    )
  }

  if (!is.null(coordinate_summary$unit)) {
    coordinate_summary$unit[is.na(coordinate_summary$unit)] <- ""
  }
  print(coordinate_summary, row.names = FALSE)

  cat("\nSupport bounding box:\n")

  .pp_print_bounds(
    bounds = x$support_bounds,
    axis_info = x$axis_info,
    digits = digits
  )

  if (!is.null(x$values)) {
    cat("\nCovariate-value summary:\n")

    value_summary <- vapply(
      x$values,
      .pp_format_number,
      character(1),
      digits = digits
    )

    print(
      as.data.frame(
        as.list(value_summary),
        check.names = FALSE
      ),
      row.names = FALSE
    )
  }

  if (!is.null(x$time)) {
    cat("\nTemporal information:\n")
    cat("  Axis: ", x$time$axis, "\n", sep = "")
    cat("  Numeric range: [",
        .pp_format_number(x$time$numeric_range[1], digits),
        ", ",
        .pp_format_number(x$time$numeric_range[2], digits),
        "]\n",
        sep = "")

    if (!is.null(x$time$original_range)) {
      cat(
        "  Original range: [",
        format(x$time$original_range[1]),
        ", ",
        format(x$time$original_range[2]),
        "]\n",
        sep = ""
      )
    }
  }

  invisible(x)
}


#' Evaluate one covariate at target locations
#'
#' Evaluates a continuous covariate field by inverse-distance weighting or
#' nearest-neighbour assignment, or evaluates a geometric support as distance to
#' the nearest support location.
#'
#' @details
#' For a field and `method = "idw"`, the interpolated value at target location
#' \eqn{u} is
#'
#' \deqn{\widehat z(u) =
#' \frac{\sum_{i \in N_k(u)} d(u,u_i)^{-p} z_i}
#' {\sum_{i \in N_k(u)} d(u,u_i)^{-p}},}
#'
#' where \eqn{z_i} is the value at support location \eqn{u_i}, \eqn{p>0}, and
#' \eqn{N_k(u)} contains all support locations or the `k` nearest locations.
#' If a target coincides with one or more support locations, the mean of the
#' coincident values is returned directly.
#'
#' With `method = "nearest"`, the value at a nearest support location is
#' assigned; tied nearest values are averaged. For a geometric support,
#' `method = "distance_min"` returns the minimum scaled Euclidean distance to
#' the support.
#'
#' Distances are computed as
#'
#' \deqn{d(u,v) = \left\{\sum_j [(u_j-v_j)/s_j]^2\right\}^{1/2},}
#'
#' where \eqn{s_j} is the corresponding entry of `axis_scale`. The default is
#' one on every axis. Because spatial and temporal coordinates usually have
#' different units and scales, a warning is issued when a spatio-temporal
#' covariate is evaluated without an explicit `axis_scale`.
#'
#' Target objects must contain every axis used by the covariate and must have
#' compatible roles, units, temporal origin, and coordinate reference system.
#' Evaluation is performed in chunks to limit the size of the distance matrix.
#'
#' @param covariate A `pp_covariate` object created by [pp_covariate()]. Fields
#'   support `"idw"` and `"nearest"`; geometric supports use
#'   `"distance_min"`.
#' @param at Target locations. Supply a `pp` object, a `pp_grid` object, or a
#'   data frame/matrix containing all coordinate axes required by `covariate`.
#' @param method Optional evaluation method. Allowed values are `"idw"` and
#'   `"nearest"` for a covariate field and `"distance_min"` for a geometric
#'   support. The default is `"idw"` for a field and `"distance_min"` for a
#'   support.
#' @param p Strictly positive finite power used by inverse-distance weighting.
#'   Larger values give more influence to nearby support points. Ignored for
#'   methods other than `"idw"`.
#' @param k Optional positive integer giving the number of nearest support
#'   locations used by inverse-distance weighting. `NULL` uses the full support.
#'   Ignored for methods other than `"idw"`.
#' @param axis_scale Optional positive numeric vector of axis-specific scales.
#'   It may be unnamed, in covariate-axis order, or named exactly with the
#'   covariate axes. Coordinate differences are divided by these scales before
#'   Euclidean distances are calculated. `NULL` uses one for every axis.
#' @param chunk_size Positive integer giving the maximum number of target
#'   locations processed in each distance-matrix chunk. Smaller values reduce
#'   peak memory use at the cost of additional loop overhead.
#'
#' @return A numeric vector with one evaluated value for each target location,
#' in the same row order as `at`.
#'
#' @references
#' Shepard, D. (1968). A two-dimensional interpolation function for
#' irregularly-spaced data. *Proceedings of the 1968 ACM National Conference*,
#' 517--524. \doi{10.1145/800186.810616}.
#'
#' @examples
#' w <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' support <- expand.grid(
#'   x = seq(0, 1, length.out = 4),
#'   y = seq(0, 1, length.out = 4)
#' )
#' support$value <- with(support, sin(pi * x) + y)
#'
#' field <- pp_covariate(
#'   support,
#'   coords = c("x", "y"),
#'   value = "value",
#'   window = w,
#'   name = "field"
#' )
#' targets <- data.frame(
#'   x = c(0.15, 0.50, 0.85),
#'   y = c(0.20, 0.50, 0.80)
#' )
#'
#' evaluate_covariate(field, targets, method = "idw", p = 2)
#' evaluate_covariate(field, targets, method = "nearest")
#'
#' sources <- pp_covariate(
#'   data.frame(x = c(0.2, 0.8), y = c(0.2, 0.7)),
#'   coords = c("x", "y"),
#'   value = NULL,
#'   window = w,
#'   name = "source_distance"
#' )
#' evaluate_covariate(sources, targets)
#'
#' @seealso
#' [pp_covariate()] to construct covariates, [extract_covariates()] to evaluate
#' several covariates together, and [pp_grid()] to create regular targets.
#'
#' @family covariate tools
#' @export
evaluate_covariate <- function(
    covariate,
    at,
    method = NULL,
    p = 2,
    k = NULL,
    axis_scale = NULL,
    chunk_size = 1000L
) {
  .pp_cov_validate_object(covariate)

  method <- .pp_cov_validate_method(
    covariate = covariate,
    method = method
  )

  if (!is.numeric(chunk_size) ||
      length(chunk_size) != 1L ||
      !is.finite(chunk_size) ||
      chunk_size < 1L ||
      chunk_size != as.integer(chunk_size)) {
    stop(
      "`chunk_size` must be one positive integer.",
      call. = FALSE
    )
  }

  chunk_size <- as.integer(chunk_size)

  if (method == "idw") {
    if (!is.numeric(p) ||
        length(p) != 1L ||
        !is.finite(p) ||
        p <= 0) {
      stop(
        "`p` must be one strictly positive finite number.",
        call. = FALSE
      )
    }

    if (!is.null(k)) {
      if (!is.numeric(k) ||
          length(k) != 1L ||
          !is.finite(k) ||
          k < 1L ||
          k != as.integer(k)) {
        stop(
          "`k` must be NULL or one positive integer.",
          call. = FALSE
        )
      }

      k <- as.integer(k)
    }
  }

  if (nrow(covariate$coords) == 0L) {
    stop(
      "The covariate has no support locations to evaluate.",
      call. = FALSE
    )
  }

  target <- .pp_cov_target_coordinates(
    covariate = covariate,
    at = at
  )

  axis_scale <- .pp_cov_validate_axis_scale(
    covariate = covariate,
    axis_scale = axis_scale
  )

  n_target <- nrow(target$coords)

  if (n_target == 0L) {
    return(numeric(0))
  }

  result <- numeric(n_target)

  chunks <- .pp_cov_chunk_indices(
    n = n_target,
    chunk_size = chunk_size
  )

  for (index in chunks) {
    distance_matrix <- .pp_cov_distance_matrix(
      target = target$coords[index, , drop = FALSE],
      source = covariate$coords,
      axis_scale = axis_scale
    )

    if (method == "distance_min") {
      result[index] <- apply(
        distance_matrix,
        1L,
        min
      )
    } else if (method == "nearest") {
      result[index] <- .pp_cov_nearest_values(
        distance_matrix = distance_matrix,
        values = covariate$values
      )
    } else {
      result[index] <- .pp_cov_idw_values(
        distance_matrix = distance_matrix,
        values = covariate$values,
        p = p,
        k = k
      )
    }
  }

  names(result) <- NULL

  result
}


#' Extract multiple covariates at common target locations
#'
#' Evaluates several covariate fields or geometric supports at the same target
#' locations and combines the results in a single data frame.
#'
#' @details
#' Each covariate is evaluated independently by [evaluate_covariate()]. Methods
#' may be supplied by covariate name or in list order. Additional arguments such
#' as the inverse-distance power, nearest-neighbour count, axis scales, or chunk
#' size are supplied separately for each covariate through `method_args`.
#'
#' Covariates may involve different subsets of the target axes. If
#' `keep_coords = TRUE`, the output retains every target coordinate required by
#' at least one covariate. Covariate names must be unique and must not conflict
#' with retained coordinate names.
#'
#' @param covariates A non-empty list of `pp_covariate` objects. The list should
#'   have unique names. If it is unnamed, the internal `name` of each covariate
#'   is used, and those names must be unique.
#' @param at Common target locations supplied as a `pp` object, `pp_grid`
#'   object, data frame, or matrix. It must contain all axes required by the
#'   covariates.
#' @param methods Optional character vector with one evaluation method per
#'   covariate. It may be named exactly with the covariate-list names or unnamed
#'   in list order. `NULL` uses the default method for each object: `"idw"` for
#'   fields and `"distance_min"` for geometric supports.
#' @param method_args Optional named list of per-covariate argument lists passed
#'   to [evaluate_covariate()]. Names must be a subset of the covariate-list
#'   names. The protected arguments `covariate`, `at`, and `method` cannot be
#'   included.
#' @param keep_coords Logical value indicating whether the required target
#'   coordinate columns should be included before the evaluated covariates.
#'
#' @return A data frame with one row per target location. It contains the
#' required target coordinates when `keep_coords = TRUE`, followed by one
#' numeric column for each covariate.
#'
#' @examples
#' w <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' support <- expand.grid(
#'   x = seq(0, 1, length.out = 5),
#'   y = seq(0, 1, length.out = 5)
#' )
#' support$cov1 <- with(support, 10 + x - 2 * y)
#'
#' cov1 <- pp_covariate(
#'   support,
#'   coords = c("x", "y"),
#'   value = "cov1",
#'   window = w,
#'   name = "cov1"
#' )
#' stations <- pp_covariate(
#'   data.frame(x = c(0.2, 0.8), y = c(0.25, 0.75)),
#'   coords = c("x", "y"),
#'   window = w,
#'   name = "station_distance"
#' )
#'
#' target_grid <- pp_grid(w, n = 4)
#' values <- extract_covariates(
#'   covariates = list(
#'     cov1 = cov1,
#'     station_distance = stations
#'   ),
#'   at = target_grid,
#'   methods = c(
#'     cov1 = "idw",
#'     station_distance = "distance_min"
#'   ),
#'   method_args = list(
#'     cov1 = list(p = 2, k = 4)
#'   )
#' )
#' head(values)
#'
#' @seealso
#' [evaluate_covariate()] for single-covariate evaluation, [pp_covariate()] for
#' construction, and [pp_grid()] for regular targets.
#'
#' @family covariate tools
#' @export
extract_covariates <- function(
    covariates,
    at,
    methods = NULL,
    method_args = NULL,
    keep_coords = TRUE
) {
  if (!is.list(covariates) || length(covariates) == 0L) {
    stop(
      "`covariates` must be a non-empty list of pp_covariate objects.",
      call. = FALSE
    )
  }

  for (covariate in covariates) {
    .pp_cov_validate_object(covariate)
  }

  covariate_names <- names(covariates)

  if (is.null(covariate_names) ||
      any(!nzchar(covariate_names))) {
    covariate_names <- vapply(
      covariates,
      function(z) z$name,
      character(1)
    )
  }

  if (anyDuplicated(covariate_names)) {
    stop(
      "`covariates` must have unique names, either in the list or in covariate objects.",
      call. = FALSE
    )
  }

  names(covariates) <- covariate_names

  if (is.null(methods)) {
    methods <- vapply(
      covariates,
      .pp_cov_default_method,
      character(1)
    )
  } else {
    if (!is.character(methods)) {
      stop(
        "`methods` must be NULL or a character vector.",
        call. = FALSE
      )
    }

    if (!is.null(names(methods))) {
      missing_methods <- setdiff(covariate_names, names(methods))
      extra_methods <- setdiff(names(methods), covariate_names)

      if (length(missing_methods) > 0L ||
          length(extra_methods) > 0L) {
        stop(
          "`methods` names must match the covariate-list names exactly.",
          call. = FALSE
        )
      }

      methods <- methods[covariate_names]
    } else if (length(methods) != length(covariates)) {
      stop(
        "`methods` must have one entry per covariate.",
        call. = FALSE
      )
    }
  }

  if (is.null(method_args)) {
    method_args <- list()
  }

  if (!is.list(method_args)) {
    stop(
      "`method_args` must be NULL or a named list.",
      call. = FALSE
    )
  }

  if (length(method_args) > 0L) {
    if (is.null(names(method_args)) ||
        any(!nzchar(names(method_args)))) {
      stop(
        "`method_args` must be a named list.",
        call. = FALSE
      )
    }

    unknown_args <- setdiff(names(method_args), covariate_names)

    if (length(unknown_args) > 0L) {
      stop(
        sprintf(
          "`method_args` contains unknown covariate name(s): %s.",
          paste(unknown_args, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  required_axes <- unique(
    unlist(
      lapply(
        covariates,
        function(z) z$axis_info$name
      ),
      use.names = FALSE
    )
  )

  output <- if (isTRUE(keep_coords)) {
    .pp_cov_output_coordinates(
      at = at,
      required_axes = required_axes
    )
  } else {
    data.frame()
  }

  coordinate_names <- names(output)

  duplicated_names <- intersect(
    coordinate_names,
    covariate_names
  )

  if (length(duplicated_names) > 0L) {
    stop(
      sprintf(
        "Covariate name(s) conflict with coordinate name(s): %s.",
        paste(duplicated_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  for (i in seq_along(covariates)) {
    name_i <- covariate_names[i]
    args_i <- method_args[[name_i]]

    if (is.null(args_i)) {
      args_i <- list()
    }

    if (!is.list(args_i)) {
      stop(
        sprintf(
          "`method_args[[\"%s\"]]` must be a list.",
          name_i
        ),
        call. = FALSE
      )
    }

    protected_names <- c("covariate", "at", "method")

    if (any(names(args_i) %in% protected_names)) {
      stop(
        sprintf(
          "`method_args[[\"%s\"]]` cannot contain: %s.",
          name_i,
          paste(protected_names, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    evaluation_args <- c(
      list(
        covariate = covariates[[i]],
        at = at,
        method = methods[i]
      ),
      args_i
    )

    output[[name_i]] <- do.call(
      evaluate_covariate,
      evaluation_args
    )
  }

  rownames(output) <- NULL

  output
}


#' Create a regular grid inside a point-process window
#'
#' Creates the Cartesian product of regular one-dimensional sequences spanning
#' every axis of a [pp_window()]. The resulting object is intended for
#' covariate evaluation, prediction, visualisation, and numerical diagnostics.
#'
#' @details
#' For axis-specific resolutions \eqn{n_1,\ldots,n_d}, the grid contains
#'
#' \deqn{N_{\mathrm{grid}} = \prod_{j=1}^{d} n_j}
#'
#' points, including the lower and upper boundary on every axis. Grid size can
#' therefore grow rapidly with dimension. `max_points` is checked before the
#' Cartesian product is constructed and prevents accidental creation of an
#' excessively large object.
#'
#' The grid inherits all axis metadata, temporal information, and the coordinate
#' reference system from `window`. A covariate already observed on a suitable
#' support does not require conversion to a `pp_grid`; this constructor is for
#' new regular target locations.
#'
#' @param window A `pp_window` object defining the grid bounds and axis
#'   metadata.
#' @param n Positive integer resolution per axis. Supply one value to use the
#'   same resolution on every axis, an unnamed vector with one value per axis,
#'   or a vector named exactly with the window axes. The default is 50 points
#'   per axis.
#' @param max_points Positive finite number giving the maximum permitted total
#'   number of grid points. An error is raised before allocation if the product
#'   of the requested resolutions exceeds this limit.
#'
#' @return An object of class `"pp_grid"` with components `coords`, `window`,
#' `axis_info`, `time_info`, `crs`, `resolution`, `dimension`, and `n_points`.
#'
#' @examples
#' # A regular two-dimensional grid
#' w <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 2)
#' )
#' g <- pp_grid(w, n = c(x = 5, y = 4))
#' g
#' head(g$coords)
#'
#' # A space-time grid with different resolutions by axis
#' wst <- pp_window(
#'   lower = c(x = 0, y = 0, time = 0),
#'   upper = c(x = 1, y = 1, time = 10),
#'   axis_roles = c(x = "space", y = "space", time = "time"),
#'   axis_units = c(x = "km", y = "km", time = "days"),
#'   time_origin = as.Date("2024-01-01")
#' )
#' gst <- pp_grid(wst, n = c(x = 4, y = 4, time = 3))
#' summary(gst)
#'
#' @seealso
#' [pp_window()] for window construction, [evaluate_covariate()] and
#' [extract_covariates()] for grid evaluation, and [predict_pp_fit()] for model
#' prediction at grid locations.
#'
#' @family point-pattern data
#' @family covariate tools
#' @export
pp_grid <- function(
    window,
    n = 50L,
    max_points = 1e6
) {
  if (!inherits(window, "pp_window")) {
    stop(
      "`window` must be an object of class 'pp_window'.",
      call. = FALSE
    )
  }

  if (!is.numeric(max_points) ||
      length(max_points) != 1L ||
      !is.finite(max_points) ||
      max_points < 1) {
    stop(
      "`max_points` must be one positive finite number.",
      call. = FALSE
    )
  }

  axis_names <- window$axis_info$name
  d <- length(axis_names)

  if (!is.numeric(n) ||
      any(!is.finite(n)) ||
      any(n != as.integer(n))) {
    stop(
      "`n` must contain finite integer grid resolutions.",
      call. = FALSE
    )
  }

  if (!is.null(names(n))) {
    if (anyDuplicated(names(n))) {
      stop(
        "`n` cannot contain duplicated axis names.",
        call. = FALSE
      )
    }

    missing_axes <- setdiff(axis_names, names(n))
    extra_axes <- setdiff(names(n), axis_names)

    if (length(missing_axes) > 0L ||
        length(extra_axes) > 0L) {
      stop(
        sprintf(
          "`n` must be named exactly with these axes: %s.",
          paste(axis_names, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    n <- n[axis_names]
  } else if (length(n) == 1L) {
    n <- rep(n, d)
  } else if (length(n) != d) {
    stop(
      sprintf(
        "`n` must have length 1 or length %d.",
        d
      ),
      call. = FALSE
    )
  }

  n <- as.integer(n)

  if (any(n < 2L)) {
    stop(
      "Every grid axis must contain at least two points.",
      call. = FALSE
    )
  }

  total_points <- prod(as.double(n))

  if (total_points > max_points) {
    stop(
      sprintf(
        paste0(
          "The requested grid contains %.0f points, exceeding `max_points` ",
          "(%.0f). Reduce `n` or increase `max_points` deliberately."
        ),
        total_points,
        max_points
      ),
      call. = FALSE
    )
  }

  axis_values <- lapply(
    seq_len(d),
    function(j) {
      seq(
        from = window$bounds["lower", axis_names[j]],
        to = window$bounds["upper", axis_names[j]],
        length.out = n[j]
      )
    }
  )

  names(axis_values) <- axis_names

  grid_data <- do.call(
    expand.grid,
    c(
      axis_values,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  )

  coords <- as.matrix(grid_data)
  storage.mode(coords) <- "double"
  colnames(coords) <- axis_names

  out <- list(
    coords = coords,
    window = window,
    axis_info = window$axis_info,
    time_info = window$time_info,
    crs = window$crs,
    resolution = stats::setNames(n, axis_names),
    dimension = d,
    n_points = nrow(coords)
  )

  class(out) <- "pp_grid"

  out
}


#' Methods for regular point-process grids
#'
#' Print and summarise objects created by [pp_grid()].
#'
#' @param x A `pp_grid` object for `print.pp_grid()` or a `summary_pp_grid`
#'   object for `print.summary_pp_grid()`.
#' @param object A `pp_grid` object to be summarised.
#' @param digits Number of significant digits used when printing numeric values.
#' @param ... Additional arguments. They are currently ignored.
#'
#' @return `print.pp_grid()` and `print.summary_pp_grid()` invisibly return `x`.
#' `summary.pp_grid()` returns an object of class `"summary_pp_grid"` containing
#' the dimension, number of points, resolution, axis information, bounds, and
#' coordinate reference system.
#'
#' @examples
#' w <- pp_window(
#'   lower = c(x = 0, y = 0),
#'   upper = c(x = 1, y = 1)
#' )
#' g <- pp_grid(w, n = c(x = 3, y = 4))
#' print(g)
#' sg <- summary(g)
#' print(sg)
#'
#' @seealso [pp_grid()] and [pp_window()].
#' @family point-pattern data
#' @family covariate tools
#' @rdname pp_grid-methods
#' @method print pp_grid
#' @export
print.pp_grid <- function(x, ...) {
  cat("<pp_grid>\n")
  cat("Regular evaluation grid\n\n")

  cat("Dimensions: ", x$dimension, "\n", sep = "")
  cat("Grid points: ", x$n_points, "\n", sep = "")

  cat(
    "Resolution: ",
    paste(
      paste0(names(x$resolution), " = ", x$resolution),
      collapse = ", "
    ),
    "\n",
    sep = ""
  )

  axis_labels <- vapply(
    seq_len(nrow(x$axis_info)),
    function(j) {
      .pp_axis_label(
        x$axis_info$name[j],
        x$axis_info$unit[j]
      )
    },
    character(1)
  )

  cat("Axes: ", paste(axis_labels, collapse = ", "), "\n", sep = "")

  if (!is.null(x$crs)) {
    cat("CRS: ", x$crs, "\n", sep = "")
  }

  invisible(x)
}


#' @rdname pp_grid-methods
#' @method summary pp_grid
#' @export
summary.pp_grid <- function(object, ...) {
  out <- list(
    dimension = object$dimension,
    n_points = object$n_points,
    resolution = object$resolution,
    axis_info = object$axis_info,
    bounds = object$window$bounds,
    crs = object$crs
  )

  class(out) <- "summary_pp_grid"

  out
}


#' @rdname pp_grid-methods
#' @method print summary_pp_grid
#' @export
print.summary_pp_grid <- function(x, digits = 4, ...) {
  cat("Regular evaluation grid\n\n")

  cat("Dimensions: ", x$dimension, "\n", sep = "")
  cat("Grid points: ", x$n_points, "\n", sep = "")

  cat(
    "Resolution: ",
    paste(
      paste0(names(x$resolution), " = ", x$resolution),
      collapse = ", "
    ),
    "\n\n",
    sep = ""
  )

  cat("Window bounds:\n")

  .pp_print_bounds(
    bounds = x$bounds,
    axis_info = x$axis_info,
    digits = digits
  )

  invisible(x)
}



# Covariate plotting helpers ------------------------------------------------

.pp_cov_plot_limits <- function(x) {
  x <- sort(as.numeric(x))

  if (length(x) != 2L || any(!is.finite(x))) {
    stop("Plot limits must contain two finite values.", call. = FALSE)
  }

  if (x[1] == x[2]) {
    padding <- max(abs(x[1]) * 0.05, 1)
    x <- x + c(-padding, padding)
  }

  x
}


.pp_cov_plot_bounds <- function(target, axes) {
  if (!is.null(target$window) &&
      all(axes %in% target$window$axis_info$name)) {
    bounds <- target$window$bounds[, axes, drop = FALSE]
  } else {
    coords <- target$coords[, axes, drop = FALSE]

    if (nrow(coords) == 0L) {
      stop(
        "Cannot determine plot bounds for an empty target without a window.",
        call. = FALSE
      )
    }

    bounds <- rbind(
      lower = apply(coords, 2L, min),
      upper = apply(coords, 2L, max)
    )
  }

  for (j in seq_len(ncol(bounds))) {
    bounds[, j] <- .pp_cov_plot_limits(bounds[, j])
  }

  bounds
}


.pp_cov_plot_axis_labels <- function(axis_info, axes) {
  index <- match(axes, axis_info$name)

  vapply(
    index,
    function(j) {
      .pp_axis_label(
        name = axis_info$name[j],
        unit = axis_info$unit[j]
      )
    },
    character(1)
  )
}


.pp_cov_plot_target_info <- function(covariate, at = NULL) {
  if (is.null(at)) {
    return(list(
      coords = covariate$coords,
      axis_info = covariate$axis_info,
      window = covariate$window,
      crs = covariate$crs,
      time_info = covariate$time_info,
      is_grid = FALSE,
      source = "covariate_support"
    ))
  }

  if (inherits(at, "pp") || inherits(at, "pp_grid")) {
    .pp_cov_check_target_metadata(
      covariate = covariate,
      target = at
    )

    return(list(
      coords = at$coords,
      axis_info = at$axis_info,
      window = at$window,
      crs = at$crs,
      time_info = at$time_info,
      is_grid = inherits(at, "pp_grid"),
      source = "target_object"
    ))
  }

  target_data <- .pp_cov_target_from_data(
    covariate = covariate,
    at = at
  )

  list(
    coords = target_data$coords,
    axis_info = covariate$axis_info,
    window = NULL,
    crs = covariate$crs,
    time_info = covariate$time_info,
    is_grid = FALSE,
    source = "target_data"
  )
}


.pp_cov_plot_resolve_axes <- function(target, axes = NULL) {
  available_axes <- target$axis_info$name
  dimension <- length(available_axes)

  if (is.null(axes)) {
    if (dimension > 3L) {
      stop(
        paste0(
          "The target has more than three dimensions. ",
          "Specify exactly two or three axes through `axes`."
        ),
        call. = FALSE
      )
    }

    if (dimension < 2L) {
      stop(
        "Covariate plots currently require at least two axes.",
        call. = FALSE
      )
    }

    return(available_axes)
  }

  if (is.character(axes)) {
    unknown_axes <- setdiff(axes, available_axes)

    if (length(unknown_axes) > 0L) {
      stop(
        sprintf(
          "Unknown plotting axis/axes: %s.",
          paste(unknown_axes, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    selected <- axes
  } else if (is.numeric(axes)) {
    if (any(!is.finite(axes)) ||
        any(axes != as.integer(axes))) {
      stop(
        "Numeric `axes` must contain integer positions.",
        call. = FALSE
      )
    }

    axes <- as.integer(axes)

    if (any(axes < 1L | axes > dimension)) {
      stop(
        "`axes` contains positions outside the available axes.",
        call. = FALSE
      )
    }

    selected <- available_axes[axes]
  } else {
    stop(
      "`axes` must be NULL, a character vector, or numeric positions.",
      call. = FALSE
    )
  }

  if (anyDuplicated(selected)) {
    stop("`axes` cannot contain duplicates.", call. = FALSE)
  }

  if (!length(selected) %in% c(2L, 3L)) {
    stop(
      "`axes` must select exactly two or three axes.",
      call. = FALSE
    )
  }

  selected
}


.pp_cov_plot_pretty_ticks <- function(limits, n = 5L) {
  limits <- .pp_cov_plot_limits(limits)

  ticks <- pretty(limits, n = n)

  tolerance <- diff(limits) * 1e-9

  ticks <- ticks[
    ticks >= limits[1] - tolerance &
      ticks <= limits[2] + tolerance
  ]

  if (length(ticks) < 2L) {
    ticks <- seq(limits[1], limits[2], length.out = n)
  }

  unique(ticks)
}


.pp_cov_plot_tick_digits <- function(ticks) {
  ticks <- sort(unique(as.numeric(ticks)))

  if (length(ticks) < 2L) {
    return(2L)
  }

  step <- min(diff(ticks))

  if (!is.finite(step) || step <= 0) {
    return(2L)
  }

  if (step >= 10) {
    return(0L)
  }

  if (step >= 1) {
    return(1L)
  }

  max(1L, ceiling(-log10(step)) + 1L)
}


.pp_cov_plot_tick_labels <- function(ticks, digits = NULL) {
  if (is.null(digits)) {
    digits <- .pp_cov_plot_tick_digits(ticks)
  }

  labels <- formatC(
    ticks,
    format = "f",
    digits = digits
  )

  labels <- sub("0+$", "", labels)
  labels <- sub("\\.$", "", labels)

  labels
}


.pp_cov_plot_colour_spec <- function(
    values,
    palette = "YlOrRd",
    missing_col = "grey75",
    colorbar_ticks = 5L,
    colorbar_digits = NULL
) {
  values <- as.numeric(values)
  observed <- is.finite(values)

  if (!any(observed)) {
    return(list(
      colours = rep(missing_col, length(values)),
      limits = NULL,
      ticks = NULL,
      tick_labels = NULL,
      palette = grDevices::hcl.colors(80L, palette),
      missing_col = missing_col,
      has_missing = TRUE
    ))
  }

  limits <- range(values[observed])

  if (limits[1] == limits[2]) {
    padding <- max(abs(limits[1]) * 0.05, 1)
    limits <- limits + c(-padding, padding)
  }

  colours_palette <- grDevices::hcl.colors(
    n = 80L,
    palette = palette
  )

  scaled <- (values[observed] - limits[1]) / diff(limits)

  colour_index <- floor(
    scaled * (length(colours_palette) - 1L)
  ) + 1L

  colour_index <- pmax(
    1L,
    pmin(length(colours_palette), colour_index)
  )

  colours <- rep(missing_col, length(values))
  colours[observed] <- colours_palette[colour_index]

  ticks <- .pp_cov_plot_pretty_ticks(
    limits = limits,
    n = colorbar_ticks
  )

  list(
    colours = colours,
    limits = limits,
    ticks = ticks,
    tick_labels = .pp_cov_plot_tick_labels(
      ticks = ticks,
      digits = colorbar_digits
    ),
    palette = colours_palette,
    missing_col = missing_col,
    has_missing = any(!observed)
  )
}


.pp_cov_plot_colourbar_panel <- function(
    spec,
    title,
    cex_axis = 1.2,
    cex_main = 1.35
) {
  graphics::par(
    mar = c(3.0, 0.2, 4.0, 5.2)
  )

  graphics::plot.new()

  if (is.null(spec$limits)) {
    graphics::plot.window(
      xlim = c(0, 1),
      ylim = c(0, 1)
    )

    graphics::title(
      main = title,
      cex.main = cex_main,
      font.main = 2
    )

    graphics::text(
      x = 0.5,
      y = 0.5,
      labels = "All values\nmissing",
      cex = cex_axis
    )

    return(invisible(NULL))
  }

  limits <- spec$limits
  span <- diff(limits)

  missing_space <- if (spec$has_missing) 0.18 * span else 0

  graphics::plot.window(
    xlim = c(0, 1),
    ylim = c(limits[1] - missing_space, limits[2]),
    xaxs = "i",
    yaxs = "i"
  )

  y_values <- seq(
    from = limits[1],
    to = limits[2],
    length.out = length(spec$palette) + 1L
  )

  for (i in seq_along(spec$palette)) {
    graphics::rect(
      xleft = 0.20,
      ybottom = y_values[i],
      xright = 0.56,
      ytop = y_values[i + 1L],
      col = spec$palette[i],
      border = NA
    )
  }

  graphics::rect(
    xleft = 0.20,
    ybottom = limits[1],
    xright = 0.56,
    ytop = limits[2],
    border = "grey25",
    lwd = 1
  )

  graphics::axis(
    side = 4,
    at = spec$ticks,
    labels = spec$tick_labels,
    las = 1,
    cex.axis = cex_axis,
    lwd = 1,
    lwd.ticks = 1
  )

  graphics::title(
    main = title,
    cex.main = cex_main,
    font.main = 2
  )

  if (spec$has_missing) {
    y_missing <- limits[1] - 0.09 * span

    graphics::points(
      x = 0.38,
      y = y_missing,
      pch = 15,
      cex = 1.35,
      col = spec$missing_col
    )

    graphics::text(
      x = 0.64,
      y = y_missing,
      labels = "Missing",
      adj = c(0, 0.5),
      cex = cex_axis
    )
  }

  invisible(NULL)
}


.pp_cov_plot_is_regular_grid_2d <- function(coords) {
  if (ncol(coords) != 2L || nrow(coords) == 0L) {
    return(FALSE)
  }

  x_values <- sort(unique(coords[, 1]))
  y_values <- sort(unique(coords[, 2]))

  if (length(x_values) < 2L || length(y_values) < 2L) {
    return(FALSE)
  }

  if (nrow(coords) != length(x_values) * length(y_values)) {
    return(FALSE)
  }

  keys <- paste(coords[, 1], coords[, 2], sep = "\r")

  !anyDuplicated(keys)
}


.pp_cov_plot_raster_matrix <- function(coords, values) {
  x_values <- sort(unique(coords[, 1]))
  y_values <- sort(unique(coords[, 2]))

  x_index <- match(coords[, 1], x_values)
  y_index <- match(coords[, 2], y_values)

  z_values <- matrix(
    NA_real_,
    nrow = length(x_values),
    ncol = length(y_values)
  )

  z_values[cbind(x_index, y_index)] <- values

  list(
    x = x_values,
    y = y_values,
    z = z_values
  )
}


.pp_cov_plot_downsample <- function(
    coords,
    values,
    max_display
) {
  n <- nrow(coords)

  if (n <= max_display) {
    return(list(
      coords = coords,
      values = values,
      sampled = FALSE
    ))
  }

  index <- unique(
    round(
      seq(
        from = 1,
        to = n,
        length.out = max_display
      )
    )
  )

  warning(
    sprintf(
      "Only %d of %d target locations are displayed. Use a smaller grid or increase `max_display` to show all points.",
      length(index),
      n
    ),
    call. = FALSE
  )

  list(
    coords = coords[index, , drop = FALSE],
    values = values[index],
    sampled = TRUE
  )
}


.pp_cov_plot_support_overlay <- function(
    covariate,
    axes,
    show_support,
    at_is_null
) {
  if (!isTRUE(show_support) || at_is_null) {
    return(NULL)
  }

  if (!all(axes %in% covariate$axis_info$name)) {
    warning(
      paste0(
        "Covariate support locations are not shown because they do not contain ",
        "all displayed axes."
      ),
      call. = FALSE
    )

    return(NULL)
  }

  covariate$coords[, axes, drop = FALSE]
}


.pp_cov_plot_pattern_overlay <- function(
    pattern,
    target,
    axes
) {
  if (is.null(pattern)) {
    return(NULL)
  }

  if (!inherits(pattern, "pp")) {
    stop(
      "`pattern` must be NULL or an object inheriting from class 'pp'.",
      call. = FALSE
    )
  }

  missing_axes <- setdiff(axes, pattern$axis_info$name)

  if (length(missing_axes) > 0L) {
    stop(
      sprintf(
        "The supplied pattern is missing plotting axis/axes: %s.",
        paste(missing_axes, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  target_info <- target$axis_info[
    match(axes, target$axis_info$name),
    ,
    drop = FALSE
  ]

  pattern_info <- pattern$axis_info[
    match(axes, pattern$axis_info$name),
    ,
    drop = FALSE
  ]

  if (!identical(target_info$role, pattern_info$role)) {
    stop(
      "The supplied pattern has incompatible plotting-axis roles.",
      call. = FALSE
    )
  }

  if (!.pp_cov_compatible_units(
    cov_units = target_info$unit,
    target_units = pattern_info$unit
  )) {
    stop(
      "The supplied pattern has incompatible plotting-axis units.",
      call. = FALSE
    )
  }

  if (!is.null(target$crs) &&
      !is.null(pattern$crs) &&
      !identical(target$crs, pattern$crs)) {
    stop(
      "The supplied pattern has an incompatible CRS.",
      call. = FALSE
    )
  }

  pattern$coords[, axes, drop = FALSE]
}


.pp_cov_plot_2d_points_panel <- function(
    coords,
    bounds,
    axis_labels,
    axis_roles,
    main,
    colours,
    support_coords = NULL,
    pattern_coords = NULL,
    point_size = 1.2,
    point_pch = 19,
    support_size = 1.1,
    support_pch = 1,
    pattern_size = 1.2,
    pattern_pch = 4,
    show_window = TRUE,
    cex_axis = 1.25,
    cex_lab = 1.45,
    cex_main = 1.35
) {
  xlim <- .pp_cov_plot_limits(bounds[, 1])
  ylim <- .pp_cov_plot_limits(bounds[, 2])

  graphics::par(
    mar = c(5.5, 5.7, 4.1, 1.0)
  )

  graphics::plot(
    NA_real_,
    NA_real_,
    type = "n",
    xlim = xlim,
    ylim = ylim,
    axes = FALSE,
    ann = FALSE,
    bty = "n",
    xaxs = "i",
    yaxs = "i",
    asp = if (all(axis_roles == "space")) 1 else NA
  )

  graphics::axis(
    side = 1,
    cex.axis = cex_axis,
    lwd = 1,
    lwd.ticks = 1
  )

  graphics::axis(
    side = 2,
    cex.axis = cex_axis,
    las = 1,
    lwd = 1,
    lwd.ticks = 1
  )

  graphics::mtext(
    text = axis_labels[1],
    side = 1,
    line = 3.4,
    cex = cex_lab
  )

  graphics::mtext(
    text = axis_labels[2],
    side = 2,
    line = 3.8,
    cex = cex_lab
  )

  graphics::title(
    main = main,
    cex.main = cex_main,
    font.main = 2
  )

  if (isTRUE(show_window)) {
    graphics::rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = ylim[2],
      border = "grey35",
      lwd = 1.3
    )
  }

  if (nrow(coords) > 0L) {
    graphics::points(
      x = coords[, 1],
      y = coords[, 2],
      pch = point_pch,
      cex = point_size,
      col = colours
    )
  }

  if (!is.null(support_coords) &&
      nrow(support_coords) > 0L) {
    graphics::points(
      x = support_coords[, 1],
      y = support_coords[, 2],
      pch = support_pch,
      cex = support_size,
      col = "black",
      lwd = 1.2
    )
  }

  if (!is.null(pattern_coords) &&
      nrow(pattern_coords) > 0L) {
    graphics::points(
      x = pattern_coords[, 1],
      y = pattern_coords[, 2],
      pch = pattern_pch,
      cex = pattern_size,
      col = "black",
      lwd = 1.3
    )
  }

  invisible(NULL)
}


.pp_cov_plot_2d_raster_panel <- function(
    coords,
    values,
    bounds,
    axis_labels,
    axis_roles,
    main,
    palette,
    support_coords = NULL,
    pattern_coords = NULL,
    support_size = 1.1,
    support_pch = 1,
    pattern_size = 1.2,
    pattern_pch = 4,
    show_window = TRUE,
    cex_axis = 1.25,
    cex_lab = 1.45,
    cex_main = 1.35
) {
  raster_data <- .pp_cov_plot_raster_matrix(
    coords = coords,
    values = values
  )

  xlim <- .pp_cov_plot_limits(bounds[, 1])
  ylim <- .pp_cov_plot_limits(bounds[, 2])

  graphics::par(
    mar = c(5.5, 5.7, 4.1, 1.0)
  )

  graphics::image(
    x = raster_data$x,
    y = raster_data$y,
    z = raster_data$z,
    col = palette,
    xlim = xlim,
    ylim = ylim,
    axes = FALSE,
    ann = FALSE,
    xaxs = "i",
    yaxs = "i",
    useRaster = TRUE,
    asp = if (all(axis_roles == "space")) 1 else NA
  )

  graphics::axis(
    side = 1,
    cex.axis = cex_axis,
    lwd = 1,
    lwd.ticks = 1
  )

  graphics::axis(
    side = 2,
    cex.axis = cex_axis,
    las = 1,
    lwd = 1,
    lwd.ticks = 1
  )

  graphics::mtext(
    text = axis_labels[1],
    side = 1,
    line = 3.4,
    cex = cex_lab
  )

  graphics::mtext(
    text = axis_labels[2],
    side = 2,
    line = 3.8,
    cex = cex_lab
  )

  graphics::title(
    main = main,
    cex.main = cex_main,
    font.main = 2
  )

  if (isTRUE(show_window)) {
    graphics::rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = ylim[2],
      border = "grey35",
      lwd = 1.3
    )
  }

  if (!is.null(support_coords) &&
      nrow(support_coords) > 0L) {
    graphics::points(
      x = support_coords[, 1],
      y = support_coords[, 2],
      pch = support_pch,
      cex = support_size,
      col = "black",
      lwd = 1.2
    )
  }

  if (!is.null(pattern_coords) &&
      nrow(pattern_coords) > 0L) {
    graphics::points(
      x = pattern_coords[, 1],
      y = pattern_coords[, 2],
      pch = pattern_pch,
      cex = pattern_size,
      col = "black",
      lwd = 1.3
    )
  }

  invisible(NULL)
}


.pp_cov_plot_3d_panel <- function(
    coords,
    bounds,
    axis_labels,
    main,
    colours,
    support_coords = NULL,
    pattern_coords = NULL,
    point_size = 1.2,
    point_pch = 19,
    support_size = 1.0,
    support_pch = 1,
    pattern_size = 1.2,
    pattern_pch = 4,
    theta = -45,
    phi = 20,
    show_window = TRUE,
    cex_axis = 1.15,
    cex_lab = 1.3,
    cex_main = 1.3
) {
  xlim <- .pp_cov_plot_limits(bounds[, 1])
  ylim <- .pp_cov_plot_limits(bounds[, 2])
  zlim <- .pp_cov_plot_limits(bounds[, 3])

  empty_target <- nrow(coords) == 0L

  if (empty_target) {
    coords <- matrix(
      c(mean(xlim), mean(ylim), mean(zlim)),
      nrow = 1L
    )

    colours <- NA
  }

  graphics::par(
    mar = c(4.0, 3.7, 4.0, 1.4)
  )

  plot3D::scatter3D(
    x = coords[, 1],
    y = coords[, 2],
    z = coords[, 3],
    xlim = xlim,
    ylim = ylim,
    zlim = zlim,
    xlab = axis_labels[1],
    ylab = axis_labels[2],
    zlab = axis_labels[3],
    main = main,
    col = colours,
    pch = if (empty_target) NA else point_pch,
    cex = if (empty_target) 0 else point_size,
    theta = theta,
    phi = phi,
    ticktype = "detailed",
    bty = if (show_window) "b" else "n",
    colkey = FALSE,
    cex.axis = cex_axis,
    cex.lab = cex_lab,
    cex.main = cex_main
  )

  if (!is.null(support_coords) &&
      nrow(support_coords) > 0L) {
    plot3D::points3D(
      x = support_coords[, 1],
      y = support_coords[, 2],
      z = support_coords[, 3],
      col = "black",
      pch = support_pch,
      cex = support_size,
      add = TRUE
    )
  }

  if (!is.null(pattern_coords) &&
      nrow(pattern_coords) > 0L) {
    plot3D::points3D(
      x = pattern_coords[, 1],
      y = pattern_coords[, 2],
      z = pattern_coords[, 3],
      col = "black",
      pch = pattern_pch,
      cex = pattern_size,
      add = TRUE
    )
  }

  invisible(NULL)
}


#' @rdname pp_covariate-methods
#' @method plot pp_covariate
#' @export
plot.pp_covariate <- function(
    x,
    at = NULL,
    pattern = NULL,
    axes = NULL,
    method = NULL,
    p = 2,
    k = NULL,
    axis_scale = NULL,
    display = c("auto", "points", "raster"),
    show_support = FALSE,
    show_window = TRUE,
    palette = "YlOrRd",
    missing_col = "grey75",
    point_size = 1.2,
    point_pch = 19,
    support_size = 1.1,
    support_pch = NULL,
    pattern_size = 1.2,
    pattern_pch = 4,
    theta = -45,
    phi = 20,
    max_display = 25000L,
    cex_axis = 1.25,
    cex_lab = 1.45,
    cex_main = 1.35,
    cex_legend = 1.2,
    colorbar_ticks = 5L,
    colorbar_digits = NULL,
    main = NULL,
    ...
) {
  .pp_cov_validate_object(x)

  display <- match.arg(display)

  if (!is.numeric(max_display) ||
      length(max_display) != 1L ||
      !is.finite(max_display) ||
      max_display < 1L ||
      max_display != as.integer(max_display)) {
    stop(
      "`max_display` must be one positive integer.",
      call. = FALSE
    )
  }

  max_display <- as.integer(max_display)

  if (is.null(support_pch)) {
    support_pch <- if (.pp_cov_is_support(x)) 4 else 1
  }

  target <- .pp_cov_plot_target_info(
    covariate = x,
    at = at
  )

  display_axes <- .pp_cov_plot_resolve_axes(
    target = target,
    axes = axes
  )

  n_display_axes <- length(display_axes)

  display_axis_info <- target$axis_info[
    match(display_axes, target$axis_info$name),
    ,
    drop = FALSE
  ]

  axis_labels <- .pp_cov_plot_axis_labels(
    axis_info = target$axis_info,
    axes = display_axes
  )

  bounds <- .pp_cov_plot_bounds(
    target = target,
    axes = display_axes
  )

  at_is_null <- is.null(at)

  values <- NULL
  method_used <- NULL

  if (at_is_null && .pp_cov_is_field(x)) {
    values <- x$values
  }

  if (!at_is_null) {
    method_used <- .pp_cov_validate_method(
      covariate = x,
      method = method
    )

    values <- evaluate_covariate(
      covariate = x,
      at = at,
      method = method_used,
      p = p,
      k = k,
      axis_scale = axis_scale
    )
  }

  if (is.null(main)) {
    if (at_is_null && .pp_cov_is_support(x)) {
      main <- paste0("Geometric support: ", x$name)
    } else if (at_is_null && .pp_cov_is_field(x)) {
      main <- paste0("Observed covariate field: ", x$name)
    } else if (identical(method_used, "distance_min")) {
      main <- paste0("Distance to nearest ", x$name)
    } else {
      main <- paste0(
        x$name,
        " evaluated by ",
        toupper(method_used)
      )
    }

    if (length(display_axes) < ncol(target$coords)) {
      main <- paste0(
        main,
        " | projection: ",
        paste(display_axes, collapse = ", ")
      )
    }
  }

  support_coords <- .pp_cov_plot_support_overlay(
    covariate = x,
    axes = display_axes,
    show_support = show_support,
    at_is_null = at_is_null
  )

  pattern_coords <- .pp_cov_plot_pattern_overlay(
    pattern = pattern,
    target = target,
    axes = display_axes
  )

  old_par <- graphics::par(no.readonly = TRUE)

  on.exit({
    graphics::layout(matrix(1L, nrow = 1L, ncol = 1L))
    graphics::par(old_par)
  }, add = TRUE)

  # Geometric support without an evaluation target --------------------------

  if (is.null(values)) {
    if (n_display_axes == 2L) {
      graphics::layout(matrix(1L, nrow = 1L, ncol = 1L))

      .pp_cov_plot_2d_points_panel(
        coords = target$coords[, display_axes, drop = FALSE],
        bounds = bounds,
        axis_labels = axis_labels,
        axis_roles = display_axis_info$role,
        main = main,
        colours = "black",
        pattern_coords = pattern_coords,
        point_size = point_size,
        point_pch = support_pch,
        pattern_size = pattern_size,
        pattern_pch = pattern_pch,
        show_window = show_window,
        cex_axis = cex_axis,
        cex_lab = cex_lab,
        cex_main = cex_main
      )
    } else {
      graphics::layout(matrix(1L, nrow = 1L, ncol = 1L))

      .pp_cov_plot_3d_panel(
        coords = target$coords[, display_axes, drop = FALSE],
        bounds = bounds,
        axis_labels = axis_labels,
        main = main,
        colours = "black",
        pattern_coords = pattern_coords,
        point_size = point_size,
        point_pch = support_pch,
        pattern_size = pattern_size,
        pattern_pch = pattern_pch,
        theta = theta,
        phi = phi,
        show_window = show_window,
        cex_axis = cex_axis,
        cex_lab = cex_lab,
        cex_main = cex_main
      )
    }

    return(invisible(x))
  }

  colour_spec <- .pp_cov_plot_colour_spec(
    values = values,
    palette = palette,
    missing_col = missing_col,
    colorbar_ticks = colorbar_ticks,
    colorbar_digits = colorbar_digits
  )

  target_coords <- target$coords[, display_axes, drop = FALSE]

  raster_eligible <-
    n_display_axes == 2L &&
    ncol(target$coords) == 2L &&
    all(display_axes == target$axis_info$name) &&
    .pp_cov_plot_is_regular_grid_2d(target_coords)

  if (identical(display, "auto")) {
    display <- if (raster_eligible) "raster" else "points"
  }

  if (identical(display, "raster") && !raster_eligible) {
    stop(
      paste0(
        "`display = 'raster'` requires a complete regular two-dimensional ",
        "grid containing exactly the displayed axes."
      ),
      call. = FALSE
    )
  }

  if (n_display_axes == 3L && identical(display, "raster")) {
    stop(
      "`display = 'raster'` is available only for two-dimensional displays.",
      call. = FALSE
    )
  }

  graphics::layout(
    mat = matrix(1:2, nrow = 1L),
    widths = c(5.6, 1.7)
  )

  if (n_display_axes == 2L && identical(display, "raster")) {
    .pp_cov_plot_2d_raster_panel(
      coords = target_coords,
      values = values,
      bounds = bounds,
      axis_labels = axis_labels,
      axis_roles = display_axis_info$role,
      main = main,
      palette = colour_spec$palette,
      support_coords = support_coords,
      pattern_coords = pattern_coords,
      support_size = support_size,
      support_pch = support_pch,
      pattern_size = pattern_size,
      pattern_pch = pattern_pch,
      show_window = show_window,
      cex_axis = cex_axis,
      cex_lab = cex_lab,
      cex_main = cex_main
    )
  } else {
    display_data <- .pp_cov_plot_downsample(
      coords = target_coords,
      values = values,
      max_display = max_display
    )

    display_colours <- .pp_cov_plot_colour_spec(
      values = display_data$values,
      palette = palette,
      missing_col = missing_col,
      colorbar_ticks = colorbar_ticks,
      colorbar_digits = colorbar_digits
    )$colours

    if (n_display_axes == 2L) {
      .pp_cov_plot_2d_points_panel(
        coords = display_data$coords,
        bounds = bounds,
        axis_labels = axis_labels,
        axis_roles = display_axis_info$role,
        main = main,
        colours = display_colours,
        support_coords = support_coords,
        pattern_coords = pattern_coords,
        point_size = point_size,
        point_pch = point_pch,
        support_size = support_size,
        support_pch = support_pch,
        pattern_size = pattern_size,
        pattern_pch = pattern_pch,
        show_window = show_window,
        cex_axis = cex_axis,
        cex_lab = cex_lab,
        cex_main = cex_main
      )
    } else {
      .pp_cov_plot_3d_panel(
        coords = display_data$coords,
        bounds = bounds,
        axis_labels = axis_labels,
        main = main,
        colours = display_colours,
        support_coords = support_coords,
        pattern_coords = pattern_coords,
        point_size = point_size,
        point_pch = point_pch,
        support_size = support_size,
        support_pch = support_pch,
        pattern_size = pattern_size,
        pattern_pch = pattern_pch,
        theta = theta,
        phi = phi,
        show_window = show_window,
        cex_axis = cex_axis,
        cex_lab = cex_lab,
        cex_main = cex_main
      )
    }
  }

  .pp_cov_plot_colourbar_panel(
    spec = colour_spec,
    title = x$name,
    cex_axis = cex_legend,
    cex_main = cex_main
  )

  invisible(x)
}





#' Print an observation window
#'
#' Print a compact description of a hyperrectangular observation window
#' created by [pp_window()]: its axes with their roles and units, their
#' bounds, the hypervolume, and any temporal or coordinate-reference metadata.
#'
#' @param x A `pp_window` object.
#' @param digits Number of significant digits used for bounds and hypervolume.
#' @param ... Additional arguments. They are currently ignored.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' w <- pp_window(
#'   lower = c(x = 0, y = 0, time = 0),
#'   upper = c(x = 10, y = 5, time = 30),
#'   axis_roles = c(x = "space", y = "space", time = "time"),
#'   axis_units = c(x = "km", y = "km", time = "days")
#' )
#' w
#'
#' @seealso [pp_window()].
#' @family point-pattern data
#' @method print pp_window
#' @export
print.pp_window <- function(x, digits = 4L, ...) {
  cat("<pp_window>\n")
  cat("Hyperrectangular observation window\n\n")
  cat("Dimension: ", x$dimension, "\n", sep = "")

  info <- x$axis_info
  for (j in seq_len(nrow(info))) {
    label <- .pp_axis_label(info$name[j], info$unit[j])
    cat(
      "  ", label, " (", info$role[j], "): [",
      format(signif(x$bounds["lower", j], digits)), ", ",
      format(signif(x$bounds["upper", j], digits)), "]\n",
      sep = ""
    )
  }

  cat("Hypervolume: ", format(signif(x$hypervolume, digits)), "\n", sep = "")

  origin <- if (is.null(x$time_info)) NULL else x$time_info$origin
  if (!is.null(origin) && !(is.numeric(origin) && isTRUE(all(origin == 0)))) {
    cat("Time origin: ", format(x$time_info$origin), "\n", sep = "")
  }

  if (!is.null(x$crs)) {
    cat("CRS: ", x$crs, "\n", sep = "")
  }

  invisible(x)
}
