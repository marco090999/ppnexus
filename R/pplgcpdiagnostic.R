# LGCP diagnostics ------------------------------------------------------------
#
# This module provides LGCP-specific second-order and first-order diagnostics.
# It deliberately reuses the package's existing K-function, residual, refitting,
# simulation, and GET global-envelope infrastructure, while calibrating the
# envelopes under the fitted log-Gaussian Cox process rather than under a
# Poisson process.
#
# Required package modules:
#   ppdata.R
#   ppfit.R
#   ppdiagnostic.R
#   ppresiduals.R
#   pplgcp.R
#
# Required package:
#   GET

.pp_lgcp_diagnostic_version <- "0.1.0"


# Dependencies and validation ------------------------------------------------

.pp_lgcp_diag_require <- function(functions, context = "LGCP diagnostics") {
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
          "Source ppdata.R, ppfit.R, ppdiagnostic.R, ppresiduals.R, ",
          "and pplgcp.R before pplgcpdiagnostic.R."
        ),
        context,
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.pp_lgcp_diag_validate_fit <- function(fit) {
  .pp_lgcp_diag_require(
    c(
      ".pp_lgcp_validate_fit",
      ".pp_diag_validate_fit",
      "K_fit_pp",
      "refit_pp",
      "residuals_pp",
      "residual_profile_pp"
    ),
    context = "LGCP diagnostics"
  )

  .pp_lgcp_validate_fit(fit)
  .pp_diag_validate_fit(fit$first_order, require_model_spec = TRUE)

  if (is.null(fit$pattern) || !inherits(fit$pattern, "pp")) {
    stop(
      "The fitted LGCP does not retain a valid observed `pp` pattern.",
      call. = FALSE
    )
  }

  if (is.null(fit$theta) || !is.numeric(fit$theta) ||
      length(fit$theta) == 0L || any(!is.finite(fit$theta))) {
    stop(
      "The fitted LGCP does not contain valid finite covariance estimates.",
      call. = FALSE
    )
  }

  if (is.null(fit$second_order) || !is.list(fit$second_order)) {
    stop(
      "The fitted LGCP does not retain a valid second-order fit.",
      call. = FALSE
    )
  }

  invisible(fit)
}


.pp_lgcp_diag_validate_logical <- function(x, arg_name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", arg_name), call. = FALSE)
  }

  isTRUE(x)
}


.pp_lgcp_diag_validate_named_list <- function(x, arg_name) {
  if (is.null(x)) {
    return(list())
  }

  if (!is.list(x)) {
    stop(sprintf("`%s` must be a named list or NULL.", arg_name), call. = FALSE)
  }

  if (length(x) > 0L &&
      (is.null(names(x)) || any(is.na(names(x))) || any(!nzchar(names(x))))) {
    stop(sprintf("`%s` must be a named list.", arg_name), call. = FALSE)
  }

  x
}


.pp_lgcp_diag_validate_simulation_args <- function(simulation_args) {
  simulation_args <- .pp_lgcp_diag_validate_named_list(
    simulation_args,
    "simulation_args"
  )

  prohibited <- intersect(
    names(simulation_args),
    c("object", "fit", "nsim", "seed", "output")
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

  simulation_args
}


.pp_lgcp_diag_validate_optim_control <- function(optim_control) {
  if (is.null(optim_control)) {
    return(list(maxit = 200L))
  }

  if (!is.list(optim_control)) {
    stop("`optim_control` must be a list or NULL.", call. = FALSE)
  }

  optim_control
}


.pp_lgcp_diag_validate_refit_flag <- function(refit) {
  .pp_lgcp_diag_validate_logical(refit, "refit")
}


.pp_lgcp_diag_warn_fractional_weights <- function(fit) {
  if (isTRUE(fit$first_order$is_fractionally_weighted)) {
    warning(
      paste0(
        "The retained first-order fit used fractional event weights. ",
        "These data-specific weights are not automatically reusable for simulated ",
        "patterns, so LGCP diagnostic refits use unweighted simulated events. ",
        "Interpret the resulting envelopes with caution."
      ),
      call. = FALSE
    )
  }

  invisible(NULL)
}


.pp_lgcp_diag_covariance_model <- function(fit) {
  model <- fit$model_spec$covariance_model

  if (is.null(model)) {
    model <- fit$covariance_model
  }

  if (!is.character(model) || length(model) != 1L ||
      is.na(model) || !nzchar(model)) {
    stop(
      "The fitted LGCP does not retain a valid covariance-model name.",
      call. = FALSE
    )
  }

  model
}


.pp_lgcp_diag_covariance_function <- function(fit) {
  # Current pplgcp.R objects use `covariance_function`. The final fallback via
  # `[[]]` also allows diagnostics to read early development objects that used
  # the reserved word "function" as a list key.
  covariance_function <- fit$model_spec$covariance_function

  if (is.null(covariance_function)) {
    covariance_function <- fit$covariance_function
  }

  if (is.null(covariance_function) &&
      !is.null(fit$model_spec[["function"]])) {
    covariance_function <- fit$model_spec[["function"]]
  }

  if (is.null(covariance_function) &&
      !is.null(fit[["function"]])) {
    covariance_function <- fit[["function"]]
  }

  if (!is.null(covariance_function) && !is.function(covariance_function)) {
    stop(
      "The covariance function retained by `fit` is not a function.",
      call. = FALSE
    )
  }

  covariance_function
}


.pp_lgcp_diag_resolve_correction <- function(fit, correction = NULL) {
  if (is.null(correction)) {
    correction <- fit$settings$correction
  }

  if (is.null(correction)) {
    correction <- fit$model_spec$correction
  }

  if (is.null(correction)) {
    correction <- "translate"
  }

  match.arg(correction, c("translate", "none"))
}


.pp_lgcp_diag_resolve_integration_n <- function(fit) {
  integration_n <- fit$settings$integration_n

  if (is.null(integration_n)) {
    integration_n <- fit$model_spec$integration_n
  }

  if (is.null(integration_n)) {
    axes <- .pp_lgcp_axis_sets(fit$pattern$axis_info)
    integration_n <- if (axes$is_st) c(100L, 100L) else 400L
  }

  integration_n
}


.pp_lgcp_diag_resolve_contrast_q <- function(fit) {
  q <- fit$second_order$contrast_q

  if (is.null(q)) {
    q <- fit$settings$contrast_q
  }

  if (is.null(q)) {
    q <- 1 / 4
  }

  q
}


.pp_lgcp_diag_resolve_contrast_p <- function(fit) {
  p <- fit$second_order$contrast_p

  if (is.null(p)) {
    p <- fit$settings$contrast_p
  }

  if (is.null(p)) {
    p <- 2
  }

  p
}


.pp_lgcp_diag_recover_user_contrast_weight <- function(fit, empirical_K) {
  resolved <- fit$second_order$weight

  if (is.null(resolved)) {
    return(NULL)
  }

  same_r <- !is.null(fit$empirical_K$r) &&
    isTRUE(all.equal(as.numeric(empirical_K$r), as.numeric(fit$empirical_K$r)))

  same_u <- if (is.null(empirical_K$u) && is.null(fit$empirical_K$u)) {
    TRUE
  } else if (!is.null(empirical_K$u) && !is.null(fit$empirical_K$u)) {
    isTRUE(all.equal(as.numeric(empirical_K$u), as.numeric(fit$empirical_K$u)))
  } else {
    FALSE
  }

  if (!same_r || !same_u) {
    return(NULL)
  }

  if (!is.numeric(resolved) || length(resolved) != length(empirical_K$Khat)) {
    return(NULL)
  }

  geometry <- if (is.matrix(empirical_K$Khat)) {
    outer(
      .pp_lgcp_trapezoid_weights(empirical_K$r),
      .pp_lgcp_trapezoid_weights(empirical_K$u),
      FUN = "*"
    )
  } else {
    .pp_lgcp_trapezoid_weights(empirical_K$r)
  }

  if (length(geometry) != length(resolved) || any(geometry <= 0)) {
    return(NULL)
  }

  user_weight <- resolved / geometry

  if (is.matrix(empirical_K$Khat)) {
    dim(user_weight) <- dim(empirical_K$Khat)
  }

  if (any(!is.finite(user_weight)) || any(user_weight < 0)) {
    return(NULL)
  }

  user_weight
}


.pp_lgcp_diag_default_r <- function(fit, r) {
  if (!is.null(r)) {
    return(r)
  }

  if (!is.null(fit$settings$r)) {
    return(fit$settings$r)
  }

  if (!is.null(fit$empirical_K$r)) {
    return(fit$empirical_K$r)
  }

  NULL
}


.pp_lgcp_diag_default_u <- function(fit, u) {
  if (!is.null(u)) {
    return(u)
  }

  if (!is.null(fit$settings$u)) {
    return(fit$settings$u)
  }

  if (!is.null(fit$empirical_K$u)) {
    return(fit$empirical_K$u)
  }

  NULL
}


.pp_lgcp_diag_validate_empty_pattern <- function(pattern, simulation_id) {
  if (nrow(pattern$coords) == 0L) {
    stop(
      paste0(
        "LGCP simulation ", simulation_id, " was empty, so the model cannot ",
        "be refitted by the current `fit_pp()`. The diagnostic is unsuitable ",
        "for fitted models having substantial probability of zero events."
      ),
      call. = FALSE
    )
  }

  invisible(pattern)
}


# Shared simulation and refitting helpers ------------------------------------

.pp_lgcp_diag_simulate_one <- function(fit, simulation_args, simulation_id) {
  simulation <- tryCatch(
    do.call(
      stats::simulate,
      c(
        list(
          object = fit,
          nsim = 1L,
          output = "pattern"
        ),
        simulation_args
      )
    ),
    error = function(error) {
      stop(
        sprintf(
          "LGCP simulation %d failed: %s",
          simulation_id,
          conditionMessage(error)
        ),
        call. = FALSE
      )
    }
  )

  if (!inherits(simulation, "pp")) {
    stop(
      sprintf(
        "LGCP simulation %d did not return a `pp` pattern.",
        simulation_id
      ),
      call. = FALSE
    )
  }

  simulation
}


.pp_lgcp_diag_refit_first_order <- function(
    fit,
    pattern,
    simulation_id
) {
  tryCatch(
    refit_pp(
      fit = fit$first_order,
      pattern = pattern
    ),
    error = function(error) {
      stop(
        sprintf(
          "First-order refit %d failed: %s",
          simulation_id,
          conditionMessage(error)
        ),
        call. = FALSE
      )
    }
  )
}


.pp_lgcp_diag_refit_full <- function(
    fit,
    pattern,
    r,
    u,
    correction,
    mark_n,
    k_chunk_size,
    max_mark_combinations,
    optim_control,
    simulation_id,
    verbose = FALSE
) {
  first_order <- .pp_lgcp_diag_refit_first_order(
    fit = fit,
    pattern = pattern,
    simulation_id = simulation_id
  )

  empirical_K <- tryCatch(
    K_fit_pp(
      fit = first_order,
      pattern = pattern,
      r = r,
      u = u,
      correction = correction,
      local = FALSE,
      chunk_size = k_chunk_size,
      mark_n = mark_n,
      max_mark_combinations = max_mark_combinations
    ),
    error = function(error) {
      stop(
        sprintf(
          "K-function calculation for full LGCP refit %d failed: %s",
          simulation_id,
          conditionMessage(error)
        ),
        call. = FALSE
      )
    }
  )

  cov_model <- .pp_lgcp_diag_covariance_model(fit)
  covariance_function <- .pp_lgcp_diag_covariance_function(fit)
  contrast_weight <- .pp_lgcp_diag_recover_user_contrast_weight(
    fit = fit,
    empirical_K = empirical_K
  )

  second_order <- tryCatch(
    .pp_lgcp_mincontrast(
      empirical_K = empirical_K,
      window = pattern$window,
      cov_model = cov_model,
      covariance_function = covariance_function,
      start = fit$theta,
      lower = fit$second_order$lower,
      upper = fit$second_order$upper,
      contrast_q = .pp_lgcp_diag_resolve_contrast_q(fit),
      contrast_p = .pp_lgcp_diag_resolve_contrast_p(fit),
      contrast_weight = contrast_weight,
      integration_n = .pp_lgcp_diag_resolve_integration_n(fit),
      optim_control = optim_control,
      verbose = verbose
    ),
    error = function(error) {
      stop(
        sprintf(
          "Second-order LGCP refit %d failed: %s",
          simulation_id,
          conditionMessage(error)
        ),
        call. = FALSE
      )
    }
  )

  theoretical_K <- empirical_K
  theoretical_K$Khat <- second_order$theoretical_K
  theoretical_K$residual <- empirical_K$Khat - second_order$theoretical_K
  theoretical_K$Ktheo <- second_order$theoretical_K
  theoretical_K$source <- "refitted_lgcp"

  model_spec <- fit$model_spec
  model_spec$covariance_model <- cov_model
  model_spec$covariance_function <- covariance_function
  model_spec$covariance_params <- second_order$theta
  model_spec$correction <- correction

  settings <- fit$settings
  settings$r <- empirical_K$r
  settings$u <- empirical_K$u
  settings$correction <- correction
  settings$k_chunk_size <- k_chunk_size
  settings$mark_n <- mark_n
  settings$max_ground_mark_combinations <- max_mark_combinations

  out <- list(
    call = match.call(),
    first_order = first_order,
    theta = second_order$theta,
    covariance_model = cov_model,
    covariance_function = covariance_function,
    empirical_K = empirical_K,
    theoretical_K = theoretical_K,
    second_order = second_order,
    pattern = pattern,
    model_spec = model_spec,
    settings = settings
  )

  class(out) <- "pp_lgcp_fit"
  out
}


.pp_lgcp_diag_theoretical_k <- function(fit, empirical_K, params = fit$theta) {
  .pp_lgcp_theoretical_k(
    empirical_K = empirical_K,
    params = params,
    cov_model = .pp_lgcp_diag_covariance_model(fit),
    covariance_function = .pp_lgcp_diag_covariance_function(fit),
    integration_n = .pp_lgcp_diag_resolve_integration_n(fit)
  )
}


.pp_lgcp_diag_k_residual_object <- function(
    empirical_K,
    theoretical_K,
    source
) {
  if (!is.numeric(theoretical_K) ||
      length(theoretical_K) != length(empirical_K$Khat)) {
    stop(
      "Internal error: theoretical and empirical K-functions are not aligned.",
      call. = FALSE
    )
  }

  if (is.matrix(empirical_K$Khat)) {
    dim(theoretical_K) <- dim(empirical_K$Khat)
  }

  residual <- empirical_K$Khat - theoretical_K
  zero_reference <- residual * 0

  out <- empirical_K
  out$Khat <- residual
  out$residual <- residual
  out$Ktheo <- zero_reference
  out$source <- source
  out$diagnostic_statistic <- "Khat_minus_fitted_LGCP_K"
  out$empirical_Khat <- empirical_K$Khat
  out$fitted_Khat <- theoretical_K

  out
}


.pp_lgcp_diag_get_p_interval <- function(envelope) {
  if (exists(".pp_res_get_p_interval", mode = "function", inherits = TRUE)) {
    return(.pp_res_get_p_interval(envelope))
  }

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

  if (is.null(p)) p <- envelope$p
  if (is.null(p)) p <- envelope$p.value

  if (!is.null(p)) {
    p <- as.numeric(p)
    if (length(p) == 1L && is.finite(p)) {
      return(stats::setNames(rep(p, 2L), c("lower", "upper")))
    }
  }

  stats::setNames(c(NA_real_, NA_real_), c("lower", "upper"))
}


# Second-order LGCP diagnostic -----------------------------------------------

#' Diagnose the second-order structure of a fitted LGCP
#'
#' Performs a global-envelope diagnostic for a fitted log-Gaussian Cox process
#' using the difference between the empirical inhomogeneous ground-process
#' \eqn{K}-function and the theoretical \eqn{K}-function implied by the fitted
#' latent Gaussian covariance. Simulated envelopes are generated under the
#' fitted LGCP, not under a Poisson process.
#'
#' @details
#' Let \eqn{\widehat K_{\mathrm{obs}}} denote the empirical inhomogeneous
#' ground-process \eqn{K}-function and let
#' \eqn{K_{\widehat\theta}} denote the theoretical \eqn{K}-function implied by
#' the fitted LGCP covariance parameters. The observed diagnostic curve is
#' \deqn{D_{\mathrm{obs}}=\widehat K_{\mathrm{obs}}-K_{\widehat\theta}.}
#'
#' With `refit = TRUE`, each primary simulation is drawn from the fitted LGCP,
#' the marginal first-order model and the second-order covariance parameters are
#' refitted, and the corresponding curve
#' \deqn{D_b=\widehat K_b-K_{\widehat\theta_b}}
#' is calculated. A global envelope is then obtained with
#' `GET::global_envelope_test()`.
#'
#' With `refit = FALSE`, simulated patterns are still drawn from the fitted
#' LGCP, but the original first-order fit and covariance parameters are kept
#' fixed. This provides a faster conditional plug-in envelope and does not
#' account for model refitting.
#'
#' The procedure is a practical plug-in parametric-bootstrap diagnostic. Even
#' with `refit = TRUE`, it does not perform the nested second simulation level
#' required for a fully adjusted test of a composite hypothesis; Monte Carlo
#' p-values should therefore be interpreted as plug-in calibration.
#'
#' For a purely spatial process, the diagnostic curve is indexed by spatial lag
#' \eqn{r}. For a spatio-temporal process, the residual surface
#' \eqn{D(r,u)} is flattened for the GET calculation, consistently with the
#' package's existing spatio-temporal envelope machinery; the original
#' \eqn{(r,u)} grid is retained in the returned object.
#'
#' @param fit A `pp_lgcp_fit` object returned by `fit_lgcp_pp()`.
#' @param nsim Number of primary LGCP simulations. Values of at least 199 are
#'   recommended for exploratory diagnostics and 499 or more for more stable
#'   final Monte Carlo inference.
#' @param r Optional spatial lags. By default, the lags retained in `fit` are
#'   reused.
#' @param u Optional temporal lags for a spatio-temporal LGCP. By default, the
#'   lags retained in `fit` are reused.
#' @param correction Optional edge correction, either `"translate"` or
#'   `"none"`. The fitted model's correction is reused when `NULL`.
#' @param refit Logical. If `TRUE` (default), refit both first- and second-order
#'   components to every simulated pattern. If `FALSE`, construct a conditional
#'   plug-in envelope with the original fitted parameters held fixed.
#' @param alpha Global-envelope significance level.
#' @param type Global-envelope ordering passed to `GET`; the default is
#'   `"erl"`.
#' @param alternative Alternative passed to `GET::global_envelope_test()`.
#' @param mark_n Numerical integration resolution for continuous joint marks.
#'   When `NULL`, the package default used by `K_fit_pp()` is retained.
#' @param chunk_size Pairwise chunk size used by `K_fit_pp()`.
#' @param max_mark_combinations Maximum number of combined mark states permitted
#'   during ground-intensity marginalisation.
#' @param simulation_args Optional named list forwarded to `simulate()` for the
#'   fitted LGCP. Arguments `object`, `fit`, `nsim`, `seed`, and `output` are
#'   reserved and cannot be supplied here. This is where latent-field settings
#'   such as `field_n`, `field_method`, or thinning controls can be changed.
#' @param optim_control Control list passed to the minimum-contrast optimisation
#'   used for each full LGCP refit when `refit = TRUE`.
#' @param keep_simulations Logical; retain the simulated `pp` patterns.
#' @param keep_refits Logical; retain the refitted `pp_lgcp_fit` objects. This
#'   option is ignored when `refit = FALSE`.
#' @param seed Optional random seed controlling the sequence of LGCP simulations
#'   and any random refitting steps.
#' @param verbose Logical; report simulation progress.
#' @param ... Additional named arguments forwarded to
#'   `GET::global_envelope_test()`.
#'
#' @return An object of class `pp_lgcp_diagnostic` and `pp_envelope`. Important
#'   components are:
#'   \itemize{
#'     \item `observed`: the observed \eqn{K}-residual object;
#'     \item `simulated`: the simulated \eqn{K}-residual objects;
#'     \item `curve_set`: the `GET` curve set;
#'     \item `envelope`: the global-envelope result;
#'     \item `curve_index`: the original spatial or spatio-temporal lag grid;
#'     \item `simulations` and `refits`: optionally retained Monte Carlo objects.
#'   }
#'
#' @references
#' Møller, J., Syversveen, A. R., and Waagepetersen, R. P. (1998).
#' Log Gaussian Cox processes. *Scandinavian Journal of Statistics*, 25,
#' 451--482. doi:10.1111/1467-9469.00115.
#'
#' Baddeley, A. J., Møller, J., and Waagepetersen, R. (2000).
#' Non- and semi-parametric estimation of interaction in inhomogeneous point
#' patterns. *Statistica Neerlandica*, 54, 329--350.
#'
#' Gabriel, E. and Diggle, P. J. (2009). Second-order analysis of
#' inhomogeneous spatio-temporal point process data. *Statistica Neerlandica*,
#' 63, 43--51.
#'
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H., and Hahn, U. (2017).
#' Global envelope tests for spatial processes. *Journal of the Royal
#' Statistical Society: Series B*, 79, 381--404.
#'
#' @seealso `fit_lgcp_pp()`, `simulate.pp_lgcp_fit()`, `K_fit_pp()`,
#'   `diagnose_first_order_lgcp_pp()`, `diagnose_pp()`
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
#'   expected_n = 100,
#'   cov_model = "exponential",
#'   covariance_params = c(sigma2 = 0.4, alpha = 0.12),
#'   field_n = 12,
#'   field_method = "nearest",
#'   output = "pattern",
#'   seed = 101
#' )
#'
#' fit <- fit_lgcp_pp(
#'   pattern = x,
#'   formula = ~ 1,
#'   cov_model = "exponential",
#'   r = seq(0.01, 0.25, length.out = 15),
#'   nc = 10,
#'   control = lgcp_control(integration_n = 150)
#' )
#'
#' # Small nsim only for an executable example. Use many more in analysis.
#' diag2 <- diagnose_lgcp_pp(
#'   fit = fit,
#'   nsim = 9,
#'   alpha = 0.10,
#'   simulation_args = list(
#'     field_n = 12,
#'     field_method = "nearest",
#'     keep_field = FALSE
#'   ),
#'   seed = 103,
#'   verbose = FALSE
#' )
#'
#' print(diag2)
#' summary(diag2)
#' plot(diag2)
#' }
#'
#' @export
diagnose_lgcp_pp <- function(
    fit,
    nsim = 199L,
    r = NULL,
    u = NULL,
    correction = NULL,
    refit = TRUE,
    alpha = 0.05,
    type = "erl",
    alternative = c("two.sided", "less", "greater"),
    mark_n = NULL,
    chunk_size = 2000L,
    max_mark_combinations = 100000L,
    simulation_args = list(),
    optim_control = list(maxit = 200L),
    keep_simulations = FALSE,
    keep_refits = FALSE,
    seed = NULL,
    verbose = TRUE,
    ...
) {
  .pp_lgcp_diag_validate_fit(fit)
  .pp_lgcp_diag_warn_fractional_weights(fit)
  .pp_lgcp_diag_require(
    c(
      ".pp_lgcp_theoretical_k",
      ".pp_lgcp_mincontrast",
      ".pp_lgcp_trapezoid_weights",
      ".pp_diag_curve_set"
    ),
    context = "diagnose_lgcp_pp()"
  )

  if (!requireNamespace("GET", quietly = TRUE)) {
    stop(
      paste0(
        "Package 'GET' is required for LGCP global-envelope diagnostics; ",
        "add it to Imports and install it."
      ),
      call. = FALSE
    )
  }

  nsim <- .pp_diag_validate_positive_integer(nsim, "nsim", min_value = 2L)
  alpha <- .pp_diag_validate_probability(alpha, "alpha")
  seed <- .pp_diag_validate_seed(seed)
  alternative <- match.arg(alternative)
  refit <- .pp_lgcp_diag_validate_refit_flag(refit)
  correction <- .pp_lgcp_diag_resolve_correction(fit, correction)
  simulation_args <- .pp_lgcp_diag_validate_simulation_args(simulation_args)
  optim_control <- .pp_lgcp_diag_validate_optim_control(optim_control)
  keep_simulations <- .pp_lgcp_diag_validate_logical(
    keep_simulations,
    "keep_simulations"
  )
  keep_refits <- .pp_lgcp_diag_validate_logical(keep_refits, "keep_refits")
  verbose <- .pp_lgcp_diag_validate_logical(verbose, "verbose")

  if (!isTRUE(refit) && isTRUE(keep_refits)) {
    warning(
      "`keep_refits = TRUE` is ignored because `refit = FALSE`.",
      call. = FALSE
    )
    keep_refits <- FALSE
  }

  chunk_size <- .pp_diag_validate_positive_integer(
    chunk_size,
    "chunk_size",
    min_value = 1L
  )
  max_mark_combinations <- .pp_diag_validate_positive_integer(
    max_mark_combinations,
    "max_mark_combinations",
    min_value = 1L
  )

  r <- .pp_lgcp_diag_default_r(fit, r)
  u <- .pp_lgcp_diag_default_u(fit, u)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  observed_empirical <- K_fit_pp(
    fit = fit$first_order,
    pattern = fit$pattern,
    r = r,
    u = u,
    correction = correction,
    local = FALSE,
    chunk_size = chunk_size,
    mark_n = mark_n,
    max_mark_combinations = max_mark_combinations
  )

  observed_theoretical <- .pp_lgcp_diag_theoretical_k(
    fit = fit,
    empirical_K = observed_empirical,
    params = fit$theta
  )

  observed <- .pp_lgcp_diag_k_residual_object(
    empirical_K = observed_empirical,
    theoretical_K = observed_theoretical,
    source = "observed_fitted_lgcp"
  )

  simulations <- if (isTRUE(keep_simulations)) vector("list", nsim) else NULL
  refits <- if (isTRUE(keep_refits) && isTRUE(refit)) {
    vector("list", nsim)
  } else {
    NULL
  }
  residual_simulated <- vector("list", nsim)

  for (i in seq_len(nsim)) {
    if (isTRUE(verbose)) {
      message(sprintf("LGCP second-order diagnostic simulation %d of %d", i, nsim))
    }

    simulation <- .pp_lgcp_diag_simulate_one(
      fit = fit,
      simulation_args = simulation_args,
      simulation_id = i
    )

    .pp_lgcp_diag_validate_empty_pattern(simulation, i)

    if (isTRUE(refit)) {
      refit_i <- .pp_lgcp_diag_refit_full(
        fit = fit,
        pattern = simulation,
        r = observed_empirical$r,
        u = observed_empirical$u,
        correction = correction,
        mark_n = mark_n,
        k_chunk_size = chunk_size,
        max_mark_combinations = max_mark_combinations,
        optim_control = optim_control,
        simulation_id = i,
        verbose = FALSE
      )

      empirical_i <- refit_i$empirical_K
      theoretical_i <- refit_i$second_order$theoretical_K
    } else {
      refit_i <- NULL

      empirical_i <- tryCatch(
        K_fit_pp(
          fit = fit$first_order,
          pattern = simulation,
          r = observed_empirical$r,
          u = observed_empirical$u,
          correction = correction,
          local = FALSE,
          chunk_size = chunk_size,
          mark_n = mark_n,
          max_mark_combinations = max_mark_combinations
        ),
        error = function(error) {
          stop(
            sprintf(
              "K-function calculation for simulation %d failed: %s",
              i,
              conditionMessage(error)
            ),
            call. = FALSE
          )
        }
      )

      theoretical_i <- .pp_lgcp_diag_theoretical_k(
        fit = fit,
        empirical_K = empirical_i,
        params = fit$theta
      )
    }

    residual_simulated[[i]] <- .pp_lgcp_diag_k_residual_object(
      empirical_K = empirical_i,
      theoretical_K = theoretical_i,
      source = if (isTRUE(refit)) {
        "simulated_refitted_lgcp"
      } else {
        "simulated_conditional_lgcp"
      }
    )

    if (isTRUE(keep_simulations)) {
      simulations[[i]] <- simulation
    }

    if (isTRUE(keep_refits) && isTRUE(refit)) {
      refits[[i]] <- refit_i
    }
  }

  curves <- .pp_diag_curve_set(observed, residual_simulated)

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
    method = if (isTRUE(refit)) {
      "lgcp_plug_in_full_refit"
    } else {
      "lgcp_conditional_plugin"
    },
    fit = fit,
    observed = observed,
    observed_empirical = observed_empirical,
    observed_theoretical = observed_theoretical,
    simulated = residual_simulated,
    curve_set = curve_set,
    envelope = envelope,
    curve_index = curves$grid_index,
    settings = list(
      nsim = nsim,
      alpha = alpha,
      type = type,
      alternative = alternative,
      correction = correction,
      refit = refit,
      mark_n = mark_n,
      chunk_size = chunk_size,
      max_mark_combinations = max_mark_combinations,
      simulation_args = simulation_args,
      optim_control = optim_control,
      seed = seed
    ),
    simulations = simulations,
    refits = refits
  )

  class(out) <- c("pp_lgcp_diagnostic", "pp_envelope")
  out
}


#' @rdname diagnose_lgcp_pp
#' @param x An object of class `pp_lgcp_diagnostic`.
#' @param ... Additional arguments passed to methods.
#' @method print pp_lgcp_diagnostic
#' @export
print.pp_lgcp_diagnostic <- function(x, ...) {
  if (!inherits(x, "pp_lgcp_diagnostic")) {
    stop("`x` must inherit from class 'pp_lgcp_diagnostic'.", call. = FALSE)
  }

  p_interval <- .pp_lgcp_diag_get_p_interval(x$envelope)

  cat("<pp_lgcp_diagnostic>\n")
  cat("LGCP second-order global-envelope diagnostic\n\n")
  cat("Primary simulations: ", x$settings$nsim, "\n", sep = "")
  cat(
    "Refitting: ",
    if (isTRUE(x$settings$refit)) "full LGCP" else "none (conditional plug-in)",
    "\n",
    sep = ""
  )
  cat("Envelope ordering: ", x$settings$type, "\n", sep = "")
  cat("Alternative: ", x$settings$alternative, "\n", sep = "")
  cat("Alpha: ", x$settings$alpha, "\n", sep = "")
  cat("Correction: ", x$settings$correction, "\n", sep = "")
  cat("K type: ", x$observed$type, "\n", sep = "")

  if (all(is.finite(p_interval))) {
    cat(
      "Monte Carlo p interval: [",
      format(p_interval[["lower"]]),
      ", ",
      format(p_interval[["upper"]]),
      "]\n",
      sep = ""
    )
  }

  cat(
    "Note: calibration is plug-in and does not include nested adjustment for a composite hypothesis.\n"
  )

  invisible(x)
}


#' @rdname diagnose_lgcp_pp
#' @param object An object of class `pp_lgcp_diagnostic`.
#' @method summary pp_lgcp_diagnostic
#' @export
summary.pp_lgcp_diagnostic <- function(object, ...) {
  if (!inherits(object, "pp_lgcp_diagnostic")) {
    stop("`object` must inherit from class 'pp_lgcp_diagnostic'.", call. = FALSE)
  }

  out <- list(
    method = object$method,
    nsim = object$settings$nsim,
    refit = object$settings$refit,
    alpha = object$settings$alpha,
    type = object$settings$type,
    alternative = object$settings$alternative,
    correction = object$settings$correction,
    K_type = object$observed$type,
    covariance_model = .pp_lgcp_diag_covariance_model(object$fit),
    theta = object$fit$theta,
    p_interval = .pp_lgcp_diag_get_p_interval(object$envelope),
    maximum_absolute_observed_K_residual = max(
      abs(object$observed$Khat),
      na.rm = TRUE
    ),
    GET_summary = tryCatch(
      summary(object$envelope),
      error = function(error) NULL
    )
  )

  class(out) <- "summary_pp_lgcp_diagnostic"
  out
}


#' @rdname diagnose_lgcp_pp
#' @param digits Number of significant digits used for printing.
#' @method print summary_pp_lgcp_diagnostic
#' @export
print.summary_pp_lgcp_diagnostic <- function(x, digits = 4, ...) {
  cat("LGCP second-order diagnostic summary\n\n")
  cat("Method: ", x$method, "\n", sep = "")
  cat("Primary simulations: ", x$nsim, "\n", sep = "")
  cat(
    "Refitting: ",
    if (isTRUE(x$refit)) "full LGCP" else "none (conditional plug-in)",
    "\n",
    sep = ""
  )
  cat("Covariance model: ", x$covariance_model, "\n", sep = "")
  cat("Alpha: ", x$alpha, "\n", sep = "")
  cat("Ordering: ", x$type, "\n", sep = "")
  cat("Alternative: ", x$alternative, "\n", sep = "")
  cat("K type: ", x$K_type, "\n", sep = "")
  cat("Correction: ", x$correction, "\n\n", sep = "")

  cat("Fitted covariance parameters:\n")
  print(round(x$theta, digits = digits))

  cat(
    "\nMaximum absolute observed K residual: ",
    formatC(
      x$maximum_absolute_observed_K_residual,
      format = "fg",
      digits = digits
    ),
    "\n",
    sep = ""
  )

  if (all(is.finite(x$p_interval))) {
    cat(
      "Monte Carlo p interval: [",
      formatC(x$p_interval[["lower"]], format = "fg", digits = digits),
      ", ",
      formatC(x$p_interval[["upper"]], format = "fg", digits = digits),
      "]\n",
      sep = ""
    )
  }

  if (!is.null(x$GET_summary)) {
    cat("\nGET summary:\n")
    print(x$GET_summary)
  }

  invisible(x)
}


#' @rdname diagnose_lgcp_pp
#' @param what Either `"envelope"` or `"residual"`.
#' @method plot pp_lgcp_diagnostic
#' @export
plot.pp_lgcp_diagnostic <- function(
    x,
    what = c("envelope", "residual"),
    ...
) {
  if (!inherits(x, "pp_lgcp_diagnostic")) {
    stop("`x` must inherit from class 'pp_lgcp_diagnostic'.", call. = FALSE)
  }

  what <- match.arg(what)

  if (identical(what, "residual")) {
    plot(x$observed, what = "residual", ...)
    return(invisible(x))
  }

  if (identical(x$observed$type, "spatial")) {
    plotted <- plot(x$envelope, ...)

    if (inherits(plotted, "ggplot")) {
      plotted <- plotted +
        ggplot2::labs(x = "Spatial lag r", y = "Khat(r) - K_LGCP(r)")
      print(plotted)
    }

    return(invisible(x))
  }

  warning(
    paste0(
      "GET represents the spatio-temporal K-residual surface as a flattened ",
      "vector. Showing the observed residual surface instead; the complete ",
      "GET result is available in `$envelope`."
    ),
    call. = FALSE
  )

  plot(x$observed, what = "residual", ...)
  invisible(x)
}


# First-order LGCP diagnostic -------------------------------------------------

#' Diagnose the first-order component of a fitted LGCP
#'
#' Constructs first-order residual profiles from the fitted marginal
#' ground-process intensity and calibrates them with global envelopes obtained
#' from simulations of the fitted log-Gaussian Cox process. Unlike
#' `diagnose_first_order_pp()`, the Monte Carlo reference process is the fitted
#' LGCP rather than a Poisson process with the same marginal intensity.
#'
#' @details
#' The observed cell innovations are computed with `residuals_pp()` from the
#' marginal first-order fit stored in `fit$first_order`. For a diagnostic cell
#' \eqn{B_j}, the raw innovation is
#' \deqn{N(B_j)-\int_{B_j}\widehat\rho(u)\,du,}
#' where \eqn{\widehat\rho} is the fitted marginal ground-process intensity.
#' These innovations are aggregated by `residual_profile_pp()` along coordinate
#' axes or fitted numeric covariates.
#'
#' The default `residual = "raw"` is recommended for LGCP diagnostics. The
#' existing `"pearson"` and `"deviance"` options are also available because
#' they can be useful as diagnostic transformations, but they remain
#' working-Poisson residual contributions: in particular, the Pearson
#' denominator uses the Poisson variance rather than the larger LGCP count
#' variance. Their global envelopes are nevertheless calibrated empirically by
#' simulation from the fitted LGCP.
#'
#' With `refit = "first_order"` (the default), every simulated LGCP pattern is
#' refitted only with the retained marginal first-order model before residuals
#' are recalculated. This directly targets first-order misspecification and is
#' substantially faster than refitting the whole LGCP. With `refit = "full"`,
#' both the marginal first-order model and the second-order covariance are
#' refitted in each primary simulation, although only the refitted marginal
#' intensity enters the residual calculation.
#'
#' As for `diagnose_lgcp_pp()`, this is a plug-in parametric-bootstrap
#' procedure. It does not perform nested simulation adjustment for a composite
#' hypothesis, so Monte Carlo p-values should be interpreted accordingly.
#'
#' @param fit A `pp_lgcp_fit` object returned by `fit_lgcp_pp()`.
#' @param variable Profile variable(s), or `"all"`, in
#'   `diagnose_first_order_lgcp_pp()`. Coordinate axes and fitted numeric
#'   covariates supported by `residual_profile_pp()` can be used. In the plot
#'   method, supply retained profile variable(s); by default all are plotted.
#' @param nsim Number of primary LGCP simulations.
#' @param grid Optional `pp_residual_grid` or `pp_grid` supplied to
#'   `residuals_pp()` for the observed residual calculation.
#' @param nc Optional regular diagnostic-grid resolution when `grid = NULL`.
#' @param form Either `"cumulative"` or `"binned"`.
#' @param residual Residual contribution to aggregate. The default is `"raw"`;
#'   `"pearson"` and `"deviance"` are available as working-Poisson diagnostic
#'   transformations.
#' @param bins Number of quantile-based bins when `breaks = NULL` and
#'   `form = "binned"`.
#' @param breaks Optional profile breaks. With multiple variables, provide a
#'   named list containing one numeric break vector per variable.
#' @param refit Refit strategy: `"first_order"` (default) refits only the
#'   marginal `pp_fit`; `"full"` refits the complete LGCP in each simulation.
#' @param alpha Global-envelope significance level.
#' @param type Global-envelope ordering passed to `GET`; the default is
#'   `"erl"`.
#' @param alternative Alternative passed to `GET::global_envelope_test()`.
#' @param mark_n Numerical integration resolution for continuous joint marks.
#' @param chunk_size Maximum number of joint-intensity evaluations per block in
#'   first-order residual calculations.
#' @param max_mark_combinations Maximum permitted combined mark states during
#'   ground-intensity marginalisation.
#' @param max_cells Maximum permitted number of first-order diagnostic cells.
#' @param simulation_args Optional named list forwarded to `simulate()` for the
#'   fitted LGCP. Arguments `object`, `fit`, `nsim`, `seed`, and `output` are
#'   reserved and cannot be supplied here.
#' @param optim_control Control list used only when `refit = "full"`, and passed
#'   to each second-order minimum-contrast optimisation.
#' @param keep_simulations Logical; retain simulated `pp` patterns.
#' @param keep_refits Logical; retain refitted models. With
#'   `refit = "first_order"`, retained objects are `pp_fit`; with
#'   `refit = "full"`, they are `pp_lgcp_fit`.
#' @param seed Optional random seed controlling the sequence of LGCP simulations
#'   and any random refitting steps.
#' @param verbose Logical; report simulation progress.
#' @param ... Additional named arguments forwarded to
#'   `GET::global_envelope_test()`.
#'
#' @return An object of class `pp_lgcp_first_order_diagnostic` and
#'   `pp_residual_envelope`, containing the observed residuals and profiles,
#'   one global envelope per selected variable, and optionally retained
#'   simulations and refits.
#'
#' @references
#' Møller, J., Syversveen, A. R., and Waagepetersen, R. P. (1998).
#' Log Gaussian Cox processes. *Scandinavian Journal of Statistics*, 25,
#' 451--482. doi:10.1111/1467-9469.00115.
#'
#' Baddeley, A. J., Rubak, E., and Møller, J. (2005). Residual analysis for
#' spatial point processes (with discussion). *Journal of the Royal Statistical
#' Society: Series B*, 67, 617--666.
#'
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H., and Hahn, U. (2017).
#' Global envelope tests for spatial processes. *Journal of the Royal
#' Statistical Society: Series B*, 79, 381--404.
#'
#' Baddeley, A. J., Rubak, E., and Turner, R. (2015).
#' *Spatial Point Patterns: Methodology and Applications with R*.
#' Chapman and Hall/CRC.
#'
#' @seealso `diagnose_lgcp_pp()`, `residuals_pp()`, `residual_profile_pp()`,
#'   `diagnose_first_order_pp()`, `fit_lgcp_pp()`
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
#'   log_intensity = function(coords, params) {
#'     params$beta0 + params$beta_x * coords$x
#'   },
#'   params = list(beta0 = log(80), beta_x = 0.5),
#'   cov_model = "exponential",
#'   covariance_params = c(sigma2 = 0.3, alpha = 0.12),
#'   field_n = 12,
#'   field_method = "nearest",
#'   output = "pattern",
#'   seed = 201
#' )
#'
#' fit <- fit_lgcp_pp(
#'   pattern = x,
#'   formula = ~ x,
#'   cov_model = "exponential",
#'   r = seq(0.01, 0.25, length.out = 15),
#'   nc = 10,
#'   control = lgcp_control(integration_n = 150)
#' )
#'
#' # Small nsim only for an executable example. Use many more in analysis.
#' diag1 <- diagnose_first_order_lgcp_pp(
#'   fit = fit,
#'   variable = "x",
#'   nsim = 9,
#'   alpha = 0.10,
#'   nc = 8,
#'   residual = "raw",
#'   simulation_args = list(
#'     field_n = 12,
#'     field_method = "nearest",
#'     keep_field = FALSE
#'   ),
#'   seed = 203,
#'   verbose = FALSE
#' )
#'
#' print(diag1)
#' summary(diag1)
#' plot(diag1, variable = "x")
#' }
#'
#' @export
diagnose_first_order_lgcp_pp <- function(
    fit,
    variable = "all",
    nsim = 199L,
    grid = NULL,
    nc = NULL,
    form = c("cumulative", "binned"),
    residual = c("raw", "pearson", "deviance"),
    bins = NULL,
    breaks = NULL,
    refit = c("first_order", "full"),
    alpha = 0.05,
    type = "erl",
    alternative = c("two.sided", "less", "greater"),
    mark_n = NULL,
    chunk_size = 100000L,
    max_mark_combinations = 100000L,
    max_cells = 1e6,
    simulation_args = list(),
    optim_control = list(maxit = 200L),
    keep_simulations = FALSE,
    keep_refits = FALSE,
    seed = NULL,
    verbose = TRUE,
    ...
) {
  .pp_lgcp_diag_validate_fit(fit)
  .pp_lgcp_diag_warn_fractional_weights(fit)
  .pp_lgcp_diag_require(
    c(
      ".pp_res_profile_list",
      ".pp_lgcp_mincontrast",
      ".pp_lgcp_trapezoid_weights"
    ),
    context = "diagnose_first_order_lgcp_pp()"
  )

  if (!requireNamespace("GET", quietly = TRUE)) {
    stop(
      paste0(
        "Package 'GET' is required for LGCP first-order global envelopes; ",
        "add it to Imports and install it."
      ),
      call. = FALSE
    )
  }

  nsim <- .pp_diag_validate_positive_integer(nsim, "nsim", min_value = 2L)
  alpha <- .pp_diag_validate_probability(alpha, "alpha")
  seed <- .pp_diag_validate_seed(seed)
  form <- match.arg(form)
  residual <- match.arg(residual)
  refit <- match.arg(refit)
  alternative <- match.arg(alternative)
  simulation_args <- .pp_lgcp_diag_validate_simulation_args(simulation_args)
  optim_control <- .pp_lgcp_diag_validate_optim_control(optim_control)
  keep_simulations <- .pp_lgcp_diag_validate_logical(
    keep_simulations,
    "keep_simulations"
  )
  keep_refits <- .pp_lgcp_diag_validate_logical(keep_refits, "keep_refits")
  verbose <- .pp_lgcp_diag_validate_logical(verbose, "verbose")

  chunk_size <- .pp_diag_validate_positive_integer(
    chunk_size,
    "chunk_size",
    min_value = 1L
  )
  max_mark_combinations <- .pp_diag_validate_positive_integer(
    max_mark_combinations,
    "max_mark_combinations",
    min_value = 1L
  )
  max_cells <- .pp_diag_validate_positive_integer(
    max_cells,
    "max_cells",
    min_value = 1L
  )

  if (!is.null(seed)) {
    set.seed(seed)
  }

  observed_residuals <- residuals_pp(
    fit = fit$first_order,
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

  full_refit_r <- .pp_lgcp_diag_default_r(fit, NULL)
  full_refit_u <- .pp_lgcp_diag_default_u(fit, NULL)
  full_refit_correction <- .pp_lgcp_diag_resolve_correction(fit, NULL)
  full_k_chunk_size <- fit$settings$k_chunk_size

  if (is.null(full_k_chunk_size)) {
    full_k_chunk_size <- 2000L
  }

  full_k_max_marks <- fit$settings$max_ground_mark_combinations

  if (is.null(full_k_max_marks)) {
    full_k_max_marks <- 100000L
  }

  full_k_mark_n <- fit$settings$mark_n

  for (i in seq_len(nsim)) {
    if (isTRUE(verbose)) {
      message(sprintf("LGCP first-order diagnostic simulation %d of %d", i, nsim))
    }

    simulation <- .pp_lgcp_diag_simulate_one(
      fit = fit,
      simulation_args = simulation_args,
      simulation_id = i
    )

    .pp_lgcp_diag_validate_empty_pattern(simulation, i)

    if (identical(refit, "first_order")) {
      refit_i <- .pp_lgcp_diag_refit_first_order(
        fit = fit,
        pattern = simulation,
        simulation_id = i
      )

      first_order_i <- refit_i
    } else {
      refit_i <- .pp_lgcp_diag_refit_full(
        fit = fit,
        pattern = simulation,
        r = full_refit_r,
        u = full_refit_u,
        correction = full_refit_correction,
        mark_n = full_k_mark_n,
        k_chunk_size = full_k_chunk_size,
        max_mark_combinations = full_k_max_marks,
        optim_control = optim_control,
        simulation_id = i,
        verbose = FALSE
      )

      first_order_i <- refit_i$first_order
    }

    residuals_i <- tryCatch(
      residuals_pp(
        fit = first_order_i,
        type = residual,
        grid = observed_residuals$grid,
        mark_n = mark_n,
        chunk_size = chunk_size,
        max_mark_combinations = max_mark_combinations,
        max_cells = max_cells
      ),
      error = function(error) {
        stop(
          sprintf(
            "Residual calculation for LGCP simulation %d failed: %s",
            i,
            conditionMessage(error)
          ),
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
            paste0(
              "LGCP simulation %d produced a residual-profile grid ",
              "incompatible with the observed `%s` profile."
            ),
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
    method = paste0("lgcp_plug_in_", refit, "_refit"),
    fit = fit,
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
      refit = refit,
      mark_n = mark_n,
      chunk_size = chunk_size,
      max_mark_combinations = max_mark_combinations,
      max_cells = max_cells,
      simulation_args = simulation_args,
      optim_control = optim_control,
      seed = seed
    ),
    simulations = simulations,
    refits = refits
  )

  class(out) <- c(
    "pp_lgcp_first_order_diagnostic",
    "pp_residual_envelope"
  )

  out
}


#' @rdname diagnose_first_order_lgcp_pp
#' @param x An object of class `pp_lgcp_first_order_diagnostic`.
#' @param ... Additional arguments passed to methods.
#' @method print pp_lgcp_first_order_diagnostic
#' @export
print.pp_lgcp_first_order_diagnostic <- function(x, ...) {
  if (!inherits(x, "pp_lgcp_first_order_diagnostic")) {
    stop(
      "`x` must inherit from class 'pp_lgcp_first_order_diagnostic'.",
      call. = FALSE
    )
  }

  cat("<pp_lgcp_first_order_diagnostic>\n")
  cat("LGCP-calibrated first-order residual-profile diagnostic\n\n")
  cat("Primary simulations: ", x$settings$nsim, "\n", sep = "")
  cat("Refit strategy: ", x$settings$refit, "\n", sep = "")
  cat("Profile form: ", x$settings$form, "\n", sep = "")
  cat("Residual contribution: ", x$settings$residual, "\n", sep = "")
  cat("Variables: ", paste(names(x$envelopes), collapse = ", "), "\n", sep = "")
  cat("Alpha: ", x$settings$alpha, "\n", sep = "")
  cat(
    "Note: envelopes are simulated under the fitted LGCP; p-values use plug-in calibration without nested composite adjustment.\n"
  )

  invisible(x)
}


#' @rdname diagnose_first_order_lgcp_pp
#' @param object An object of class `pp_lgcp_first_order_diagnostic`.
#' @method summary pp_lgcp_first_order_diagnostic
#' @export
summary.pp_lgcp_first_order_diagnostic <- function(object, ...) {
  if (!inherits(object, "pp_lgcp_first_order_diagnostic")) {
    stop(
      "`object` must inherit from class 'pp_lgcp_first_order_diagnostic'.",
      call. = FALSE
    )
  }

  rows <- lapply(
    names(object$envelopes),
    function(name_i) {
      p_interval <- .pp_lgcp_diag_get_p_interval(
        object$envelopes[[name_i]]$envelope
      )

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
    refit = object$settings$refit,
    alpha = object$settings$alpha,
    form = object$settings$form,
    residual = object$settings$residual,
    p_values = do.call(rbind, rows)
  )

  class(out) <- "summary_pp_lgcp_first_order_diagnostic"
  out
}


#' @rdname diagnose_first_order_lgcp_pp
#' @param digits Number of significant digits used for printing.
#' @method print summary_pp_lgcp_first_order_diagnostic
#' @export
print.summary_pp_lgcp_first_order_diagnostic <- function(
    x,
    digits = 4,
    ...
) {
  cat("LGCP first-order diagnostic summary\n\n")
  cat("Method: ", x$method, "\n", sep = "")
  cat("Primary simulations: ", x$nsim, "\n", sep = "")
  cat("Refit strategy: ", x$refit, "\n", sep = "")
  cat("Alpha: ", x$alpha, "\n", sep = "")
  cat("Profile form: ", x$form, "\n", sep = "")
  cat("Residual contribution: ", x$residual, "\n\n", sep = "")

  print(x$p_values, row.names = FALSE, digits = digits)
  invisible(x)
}


#' @rdname diagnose_first_order_lgcp_pp
#' @param what Either `"envelope"` or `"profile"`.
#' @method plot pp_lgcp_first_order_diagnostic
#' @export
plot.pp_lgcp_first_order_diagnostic <- function(
    x,
    variable = NULL,
    what = c("envelope", "profile"),
    ...
) {
  .pp_plot_profile_envelopes(x, variable = variable, what = match.arg(what), ...)
}
