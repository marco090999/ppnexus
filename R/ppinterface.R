# User interfaces for the LGCP and ETAS extensions ----------------------------
#
# fit_lgcp_pp() and fit_etas_pp() expose only the arguments that define the
# model. Settings of the first-order (background) Poisson fit are forwarded to
# fit_pp(), and numerical settings of the second-order or triggering stage are
# collected by lgcp_control() and etas_control(). The computations are carried
# out by the internal engines .pp_fit_lgcp_engine() (pplgcp.R) and
# .pp_fit_etas_engine() (ppetas2.R).


# Shared helpers --------------------------------------------------------------

.pp_as_control <- function(control, constructor, class_name, arg_name) {
  if (is.null(control)) {
    return(constructor())
  }

  if (inherits(control, class_name)) {
    return(control)
  }

  if (!is.list(control)) {
    stop(
      sprintf("`%s` must be created with %s or be a named list.", arg_name,
              paste0(sub("^pp_", "", class_name), "()")),
      call. = FALSE
    )
  }

  allowed <- names(formals(constructor))
  unknown <- setdiff(names(control), allowed)

  if (length(unknown) > 0L || (length(control) > 0L && is.null(names(control)))) {
    stop(
      sprintf(
        "Unknown `%s` setting(s): %s. Allowed settings: %s.",
        arg_name,
        paste(unknown, collapse = ", "),
        paste(allowed, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  do.call(constructor, control)
}


.pp_check_named_dots <- function(dots, fun_name) {
  if (length(dots) > 0L &&
      (is.null(names(dots)) || any(!nzchar(names(dots))))) {
    stop(
      sprintf("All additional arguments of %s must be named.", fun_name),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.pp_interface_warning <- function(fun_name, legacy, advice) {
  warning(
    sprintf(
      paste0(
        "%s: argument(s) %s belong to the previous interface. They are still ",
        "accepted, but please use %s instead."
      ),
      fun_name,
      paste0("`", legacy, "`", collapse = ", "),
      advice
    ),
    call. = FALSE
  )
}


# LGCP --------------------------------------------------------------------------

#' Control settings for LGCP fitting
#'
#' Collects the numerical settings of the second-order (minimum-contrast) stage
#' of [fit_lgcp_pp()]. All arguments have defaults, so `lgcp_control()` only
#' needs the settings that should be changed.
#'
#' @param start,lower,upper Optional named vectors or named lists giving the
#'   starting values and box constraints of the covariance parameters. Partial
#'   named overrides are allowed for the built-in covariance models; a custom
#'   covariance function requires explicit named `start` values.
#' @param correction Edge correction used for the empirical inhomogeneous
#'   K-function: `"translate"` (default) or `"none"`.
#' @param contrast_q Positive power applied to the empirical and theoretical
#'   K-functions before they are compared. Default `1/4`.
#' @param contrast_p Positive power of the absolute transformed discrepancy.
#'   Default `2`.
#' @param contrast_weight Optional numeric weights aligned with the empirical
#'   K-function, or a function of `r`, `u`, `Khat` and `Kmodel`. They are
#'   multiplied by trapezoidal integration weights.
#' @param integration_n Resolution of the numerical integration used for the
#'   theoretical K-function: one integer for spatial models, one or two
#'   integers (spatial, temporal) for spatio-temporal models. The default is
#'   400 for spatial and `c(100, 100)` for spatio-temporal models.
#' @param optim_control Control list passed to [stats::optim()]. The method is
#'   `"L-BFGS-B"` when `lower` or `upper` contain finite bounds (always the case
#'   for the built-in covariance models) and `"Nelder-Mead"` otherwise.
#' @param k_chunk_size Chunk size for the pairwise computations of the
#'   empirical K-function.
#' @param mark_n Integration resolution used to marginalise continuous joint
#'   marks when the empirical K-function is computed.
#' @param max_ground_mark_combinations Maximum number of combined mark states
#'   allowed in that marginalisation.
#' @param verbose Logical; print progress messages.
#'
#' @return An object of class `"pp_lgcp_control"`, a named list of settings.
#'
#' @examples
#' lgcp_control(integration_n = 200, optim_control = list(maxit = 500))
#'
#' @seealso [fit_lgcp_pp()].
#' @export
lgcp_control <- function(
    start = NULL,
    lower = NULL,
    upper = NULL,
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
  out <- list(
    start = start,
    lower = lower,
    upper = upper,
    correction = match.arg(correction),
    contrast_q = contrast_q,
    contrast_p = contrast_p,
    contrast_weight = contrast_weight,
    integration_n = integration_n,
    optim_control = optim_control,
    k_chunk_size = k_chunk_size,
    mark_n = mark_n,
    max_ground_mark_combinations = max_ground_mark_combinations,
    verbose = verbose
  )

  class(out) <- "pp_lgcp_control"
  out
}


#' Fit a log-Gaussian Cox process with an enriched first-order intensity
#'
#' Fits a spatial or spatio-temporal log-Gaussian Cox process (LGCP) in two
#' stages. The first-order intensity is the multidimensional Poisson model of
#' [fit_pp()], so it can contain coordinates, spatial and spatio-temporal
#' covariates, functional covariates, categorical and continuous marks, and
#' smooth terms. The covariance parameters of the latent Gaussian field are
#' then estimated by minimum contrast on the inhomogeneous K-function weighted
#' by the fitted first-order intensity.
#'
#' @details
#' The fitted LGCP is parameterised as
#' \deqn{\Lambda(u)=\rho(u;\beta)\exp\{Z(u)-C_\theta(0)/2\},}
#' so that \eqn{\rho(u;\beta)} is the first-order intensity estimated in the
#' first stage. For a stationary isotropic Gaussian field the pair correlation
#' function is \eqn{g_\theta(h)=\exp\{C_\theta(h)\}} in the spatial case and
#' \eqn{g_\theta(h,v)=\exp\{C_\theta(h,v)\}} in the spatio-temporal case.
#'
#' For `d` spatial dimensions the theoretical K-function is
#' \deqn{K_\theta(r)=a_d\int_0^r h^{d-1}\exp\{C_\theta(h)\}\,dh,}
#' where `a_d = 2*pi^(d/2)/Gamma(d/2)` is the surface area of the unit sphere;
#' with a temporal axis,
#' \deqn{K_\theta(r,u)=2a_d\int_0^r\int_0^u
#' h^{d-1}\exp\{C_\theta(h,v)\}\,dv\,dh.}
#'
#' The minimum-contrast criterion compares transformed empirical and
#' theoretical K-functions (fourth-root transformation and squared discrepancy
#' by default; see [lgcp_control()]).
#'
#' Because structure that the first-order intensity does not capture is
#' absorbed by the latent field, an accurate first-order model matters: a
#' missing trend is interpreted as clustering and inflates the estimated
#' variance and range of the field.
#'
#' The first-stage [fit_pp()] is a working Poisson fit: its model-based
#' covariance matrix does not include the extra-Poisson variation induced by
#' the latent field.
#'
#' @param pattern A `pp` point pattern with a finite `pp_window`, at least one
#'   spatial axis, and at most one temporal axis.
#' @param formula One-sided formula of the first-order intensity, as in
#'   [fit_pp()]: coordinates, covariates, marks and smooth terms such as
#'   `s(x)` can be combined.
#' @param covariates Optional named list of `pp_covariate` objects (including
#'   score fields from [as_pp_covariates()]) used in `formula`.
#' @param continuous_marks,categorical_marks Optional mark specifications for
#'   the first-order intensity, as in [fit_pp()].
#' @param cov_model Covariance model of the latent field. Spatial models:
#'   `"exponential"`, `"gaussian"`, `"matern"`; spatio-temporal models:
#'   `"separable_exponential"`, `"gneiting"`. Alternatively, a custom
#'   stationary isotropic covariance function of any subset of `h`, `u`,
#'   `params` and `d_space`; its validity as a covariance function is the
#'   user's responsibility, and named starting values must then be given in
#'   `lgcp_control(start = ...)`.
#' @param r Spatial lags at which the K-function is compared. If `NULL`,
#'   package defaults are used.
#' @param u Temporal lags, for spatio-temporal patterns only.
#' @param ... Further arguments passed to [fit_pp()] for the first-order fit,
#'   for example `nc`, `q`, `covariate_methods`, `covariate_args`,
#'   `categorical_levels`, `event_weights`, `gam_args` or `seed`.
#' @param control Settings of the second-order stage, created with
#'   [lgcp_control()] (or a named list of its arguments).
#'
#' @return An object of class `"pp_lgcp_fit"` containing the first-order
#'   `pp_fit` (`first_order`), the estimated covariance parameters (`theta`),
#'   the empirical and theoretical K-functions, minimum-contrast details, and
#'   the model specification required for simulation.
#'
#' @references
#' Møller, J., Syversveen, A. R., and Waagepetersen, R. P. (1998).
#' Log Gaussian Cox processes. *Scandinavian Journal of Statistics*, 25,
#' 451--482. doi:10.1111/1467-9469.00115.
#'
#' Baddeley, A. J., Møller, J., and Waagepetersen, R. (2000).
#' Non- and semi-parametric estimation of interaction in inhomogeneous point
#' patterns. *Statistica Neerlandica*, 54, 329--350.
#' doi:10.1111/1467-9574.00144.
#'
#' Gabriel, E. and Diggle, P. J. (2009). Second-order analysis of
#' inhomogeneous spatio-temporal point process data. *Statistica Neerlandica*,
#' 63, 43--51. doi:10.1111/j.1467-9574.2008.00407.x.
#'
#' Gneiting, T. (2002). Nonseparable, stationary covariance functions for
#' space-time data. *Journal of the American Statistical Association*, 97,
#' 590--600. doi:10.1198/016214502760047113.
#'
#' Davies, T. M. and Hazelton, M. L. (2013). Assessing minimum contrast
#' parameter estimation for spatial and spatiotemporal log-Gaussian Cox
#' processes. *Statistica Neerlandica*, 67, 355--389.
#'
#' @seealso [lgcp_control()], [fit_pp()], [K_fit_pp()], [simulate_lgcp_pp()],
#'   [simulate.pp_lgcp_fit()], [diagnose_lgcp_pp()].
#'
#' @examples
#' w <- pp_window(lower = c(x = 0, y = 0), upper = c(x = 1, y = 1))
#'
#' # An LGCP whose first-order intensity increases with x
#' # (true values: slope 1.5, sigma2 = 0.5, alpha = 0.08)
#' x <- simulate_lgcp_pp(
#'   window = w,
#'   log_intensity = function(coords) log(150) + 1.5 * coords$x,
#'   cov_model = "exponential",
#'   covariance_params = c(sigma2 = 0.5, alpha = 0.08),
#'   field_n = 20,
#'   output = "pattern",
#'   seed = 5
#' )
#'
#' # First-order trend in x (fit_pp() settings such as nc are passed through),
#' # second-order settings in lgcp_control()
#' fit <- fit_lgcp_pp(
#'   x,
#'   formula = ~ x,
#'   cov_model = "exponential",
#'   r = seq(0.01, 0.25, length.out = 12),
#'   nc = 8,
#'   control = lgcp_control(integration_n = 150)
#' )
#' fit
#' coef(fit$first_order$gam)
#'
#' \donttest{
#' # A custom covariance function (here equal to the exponential model)
#' fit_custom <- fit_lgcp_pp(
#'   x,
#'   formula = ~ x,
#'   cov_model = function(h, params) params[["sigma2"]] * exp(-h / params[["alpha"]]),
#'   r = seq(0.01, 0.25, length.out = 12),
#'   nc = 8,
#'   control = lgcp_control(start = c(sigma2 = 1, alpha = 0.1), integration_n = 150)
#' )
#' fit_custom$theta
#' }
#'
#' @export
fit_lgcp_pp <- function(
    pattern,
    formula = ~ 1,
    covariates = NULL,
    continuous_marks = NULL,
    categorical_marks = NULL,
    cov_model = NULL,
    r = NULL,
    u = NULL,
    ...,
    control = lgcp_control()
) {
  call <- match.call()
  dots <- list(...)
  .pp_check_named_dots(dots, "fit_lgcp_pp()")

  control <- .pp_as_control(control, lgcp_control, "pp_lgcp_control", "control")

  # Previous interface: second-order settings and `covariance_function` were
  # direct arguments.
  legacy <- intersect(names(dots), c(names(formals(lgcp_control)), "covariance_function"))

  covariance_function <- NULL

  if (is.function(cov_model)) {
    covariance_function <- cov_model
    cov_model <- NULL
  }

  if (length(legacy) > 0L) {
    .pp_interface_warning(
      "fit_lgcp_pp()",
      legacy,
      "`control = lgcp_control(...)` (and a function in `cov_model` for a custom covariance)"
    )

    if ("covariance_function" %in% legacy) {
      covariance_function <- dots$covariance_function
      dots$covariance_function <- NULL
    }

    moved <- setdiff(legacy, "covariance_function")
    if (length(moved) > 0L) {
      settings <- unclass(control)
      settings[moved] <- dots[moved]
      control <- do.call(lgcp_control, settings)
      dots[moved] <- NULL
    }
  }

  first_order_args <- setdiff(
    names(formals(fit_pp)),
    c("pattern", "formula", "covariates", "continuous_marks", "categorical_marks")
  )
  unknown <- setdiff(names(dots), first_order_args)

  if (length(unknown) > 0L) {
    stop(
      sprintf(
        paste0(
          "Unknown argument(s) of fit_lgcp_pp(): %s. Additional arguments are ",
          "passed to fit_pp(); second-order settings belong in `control = lgcp_control()`."
        ),
        paste(unknown, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  fit <- do.call(
    .pp_fit_lgcp_engine,
    c(
      list(
        pattern = pattern,
        formula = formula,
        covariates = covariates,
        continuous_marks = continuous_marks,
        categorical_marks = categorical_marks,
        cov_model = cov_model,
        covariance_function = covariance_function,
        r = r,
        u = u
      ),
      dots,
      unclass(control)
    )
  )

  fit$call <- call
  fit
}


# ETAS --------------------------------------------------------------------------

.pp_etas_parameter_names <- c("mu", "k0", "c", "p", "gamma", "d", "q")

.pp_etas_background_settings <- c(
  "covariate_methods", "covariate_args", "continuous_marks",
  "categorical_marks", "categorical_levels", "q", "dummy_layout", "nc",
  "max_tiles", "max_quadrature_points", "max_mark_combinations", "seed",
  "gam_args", "mark_n"
)

# Arguments of the previous interface that the algorithm never used, kept only
# for compatibility with etasFLP::etasclass(). They are passed to the engine,
# which ignores them or stops for unsupported values (e.g. `flp = TRUE`).
.pp_etas_inactive_arguments <- c(
  "time.update", "magn.threshold.back", "hdef", "hvarx", "hvary", "thinning",
  "flp", "m1", "cat.back", "back.smooth", "sectoday", "longlat.to.km",
  "compsqm", "epsmax", "is.backconstant"
)


#' Control settings for ETAS fitting
#'
#' Collects the numerical and algorithmic settings of [fit_etas_pp()]. All
#' arguments have defaults, so `etas_control()` only needs the settings that
#' should be changed.
#'
#' @param ndeclust Maximum number of alternating iterations between
#'   maximum-likelihood estimation of the triggering parameters and the
#'   weighted Poisson fit of the background.
#' @param min.declust.iter Minimum number of alternating iterations completed
#'   before convergence can stop the fit.
#' @param tol.rho Convergence tolerance on the maximum absolute change of the
#'   eventwise posterior background probabilities.
#' @param tol.par Convergence tolerance on the maximum relative change of the
#'   estimated parameters.
#' @param stop.on.increasing.AIC Logical; stop the alternation when the AIC
#'   increases with respect to the previous iteration.
#' @param n.iterweight Number of weighted background updates within each
#'   alternating iteration.
#' @param declustering Logical; if `FALSE`, the background is fitted once with
#'   unit weights and no alternation is performed.
#' @param iterlim Maximum number of iterations of the optimiser in each
#'   maximum-likelihood step.
#' @param usenlm Logical; optimise with [stats::nlm()] (default) or, if
#'   `FALSE`, with [stats::optim()] using `method`.
#' @param method Optimisation method passed to [stats::optim()] when
#'   `usenlm = FALSE`.
#' @param ntheta Number of angular directions (at least four) used to
#'   approximate the spatial integral of the triggering function over the
#'   rectangular study region.
#' @param wp Optional initial background probabilities, one per selected
#'   event, in `[0, 1]`.
#' @param offset Numeric scalar or eventwise vector added to the triggering
#'   linear predictor defined by `formula1`.
#' @param tmax Optional upper time limit (original catalogue scale) of the
#'   fitting period.
#' @param long.range,lat.range Optional length-two vectors restricting the
#'   spatial fitting rectangle (original catalogue scale).
#' @param description Optional character string stored in the fitted object.
#' @param verbose Logical; print progress messages.
#' @param verbose.convergence Logical; print convergence measures at each
#'   alternating iteration.
#' @param trace Logical; print the optimiser trace.
#'
#' @return An object of class `"pp_etas_control"`, a named list of settings.
#'
#' @examples
#' etas_control(ndeclust = 10, verbose = TRUE)
#'
#' @seealso [fit_etas_pp()].
#' @export
etas_control <- function(
    ndeclust = 5,
    min.declust.iter = 2,
    tol.rho = 1e-3,
    tol.par = 1e-3,
    stop.on.increasing.AIC = FALSE,
    n.iterweight = 1,
    declustering = TRUE,
    iterlim = 50,
    usenlm = TRUE,
    method = "BFGS",
    ntheta = 36,
    wp = NULL,
    offset = 0,
    tmax = NULL,
    long.range = NULL,
    lat.range = NULL,
    description = "",
    verbose = FALSE,
    verbose.convergence = FALSE,
    trace = FALSE
) {
  out <- list(
    ndeclust = ndeclust,
    min.declust.iter = min.declust.iter,
    tol.rho = tol.rho,
    tol.par = tol.par,
    stop.on.increasing.AIC = stop.on.increasing.AIC,
    n.iterweight = n.iterweight,
    declustering = declustering,
    iterlim = iterlim,
    usenlm = usenlm,
    method = method,
    ntheta = ntheta,
    wp = wp,
    offset = offset,
    tmax = tmax,
    long.range = long.range,
    lat.range = lat.range,
    description = description,
    verbose = verbose,
    verbose.convergence = verbose.convergence,
    trace = trace
  )

  class(out) <- "pp_etas_control"
  out
}


#' Fit an ETAS model with an enriched parametric background
#'
#' Fits a spatio-temporal (or purely temporal) epidemic-type aftershock
#' sequence (ETAS) model whose background rate is the multidimensional Poisson
#' model of [fit_pp()]. The background can therefore contain coordinates,
#' spatial covariates, functional covariates, marks and smooth terms, with
#' interpretable effects and standard errors. The triggering component follows
#' the parametrisation of `etasFLP`.
#'
#' @details
#' For an event history \eqn{\mathcal{H}_t} the conditional intensity is
#'
#' \deqn{\lambda(t, \mathbf{s} \mid \mathcal{H}_t)
#' = \mu b(\mathbf{s};\boldsymbol{\beta})
#' + \sum_{t_j<t} h(t-t_j, \mathbf{s}-\mathbf{s}_j,
#' m_j;\boldsymbol{\theta}),}
#'
#' where \eqn{b(\mathbf{s};\boldsymbol{\beta})} is the background spatial
#' density, given by the fitted Poisson model normalised to integrate to one
#' over the window, and \eqn{h} is the triggering contribution of past events.
#' The seven ETAS parameters are `mu` (background rate), `k0` (productivity),
#' `c` and `p` (temporal decay), `d` and `q` (spatial decay) and `gamma`
#' (magnitude dependence of the spatial scale). Further productivity effects
#' enter through the linear predictor defined by `formula1`.
#'
#' Estimation alternates between maximum likelihood for the triggering
#' parameters and a weighted fit of the background with [fit_pp()], in which
#' each event is weighted by its posterior probability of being a background
#' event (stochastic declustering),
#' \deqn{\rho_i = \frac{\mu b(\mathbf{s}_i)}{\lambda(t_i,\mathbf{s}_i)}.}
#' Iteration stops when the changes of these probabilities and of the
#' parameters fall below the tolerances of [etas_control()], or after
#' `ndeclust` iterations. Triggering contributions are computed by compiled
#' Fortran code adapted from `etasFLP`.
#'
#' The spatial triggering parameters are identifiable only when the spatial
#' scale of the triggering kernel (about \eqn{\sqrt{d}}) is small compared
#' with the observation window, and when the catalogue contains enough
#' triggered events. When the Hessian of the log-likelihood is singular, the
#' standard errors are `NA` and a warning is issued.
#'
#' @section Previous interface:
#' Earlier versions accepted the starting values (`mu`, `k0`, ...),
#' `params.ind`, the background settings with suffix `.bg` (`nc.bg`, ...) and
#' the control settings as separate arguments, together with several inactive
#' arguments inherited from `etasFLP::etasclass()`. These are still accepted
#' through `...`, with a warning, and are translated to `start`, `fixed`,
#' `background` and `control`.
#'
#' @param catalog A `"pp_etas_catalog"` created by [etas_catalog()] or
#'   [as_etas_catalog()].
#' @param formula.bg One-sided formula of the background Poisson model on the
#'   working spatial coordinates `x` and `y`, for example `~ 1`, `~ x + y` or
#'   `~ s(x, y) + elevation`.
#' @param covariates.bg Optional named list of `pp_covariate` objects used in
#'   `formula.bg`.
#' @param formula1 Triggering formula with `time` on the left-hand side; its
#'   right-hand side defines the productivity linear predictor, e.g.
#'   `"time ~ magnitude - 1"`.
#' @param magn.threshold Completeness magnitude: events below it are removed
#'   and magnitudes are centred at it.
#' @param start Optional named list (or named numeric vector) of starting
#'   values for any of `mu`, `k0`, `c`, `p`, `gamma`, `d`, `q`, and `betacov`
#'   (the coefficients of `formula1`). Omitted values of `mu`, `k0`, `c` and
#'   `d` are computed from the data; the defaults of the others are
#'   `p = 1.01`, `gamma = 0`, `q = 1.5`, `betacov = 0.7`.
#' @param fixed Character vector naming the ETAS parameters that are kept at
#'   their starting values instead of being estimated. The default keeps
#'   `gamma` fixed (at zero unless set in `start`).
#' @param onlytime Logical; if `TRUE`, fit the purely temporal ETAS model
#'   (constant background; `gamma`, `d` and `q` fixed at zero).
#' @param background Optional named list of further settings of the background
#'   [fit_pp()] fit: `nc`, `q`, `covariate_methods`, `covariate_args`,
#'   `continuous_marks`, `categorical_marks`, `categorical_levels`,
#'   `dummy_layout`, `max_tiles`, `max_quadrature_points`,
#'   `max_mark_combinations`, `seed`, `gam_args`, `mark_n`.
#' @param control Algorithmic settings created with [etas_control()] (or a
#'   named list of its arguments).
#' @param ... Arguments of the previous interface only (see the section
#'   *Previous interface*).
#'
#' @return An object of class `"pp_etas_fit"`. Important components include
#'   `parameters` and `standard_errors` (ETAS and triggering coefficients),
#'   `background_fit` (the final weighted `pp_fit`), the eventwise
#'   `background_density`, `triggering_intensity`, `total_intensity`,
#'   `background_probability` and `triggered_probability`, the likelihood
#'   components `loglik`, `integral` and `trigger_integral`, the iteration
#'   `history`, and `convergence`, `optimisation` and `settings`.
#'
#' @references
#' Ogata, Y. (1988). Statistical models for earthquake occurrences and residual
#' analysis for point processes. *Journal of the American Statistical
#' Association*, 83(401), 9--27.
#'
#' Zhuang, J., Ogata, Y. and Vere-Jones, D. (2002). Stochastic declustering of
#' space-time earthquake occurrences. *Journal of the American Statistical
#' Association*, 97(458), 369--380.
#'
#' Chiodi, M. and Adelfio, G. (2017). Mixed non-parametric and parametric
#' estimation techniques in R package etasFLP for earthquakes' description.
#' *Journal of Statistical Software*, 76(3), 1--29.
#'
#' @seealso [etas_control()], [simulate_etas_pp()], [etas_catalog()],
#'   [decluster_etas_pp()], [fit_pp()].
#'
#' @examples
#' \donttest{
#' # Simulate a catalogue whose background rate increases with x
#' sim <- simulate_etas_pp(
#'   params = c(mu = 0.5, k0 = 1.5e-5, c = 0.05, p = 1.2,
#'              gamma = 0, d = 0.001, q = 2),
#'   m0 = 2.5, b = 1, tmin = 0, t.lag = 500,
#'   long.range = c(0, 1), lat.range = c(0, 1), longlat.to.km = FALSE,
#'   n_support_obs = 50, mult_support = 1,
#'   bg_lp_type = "linear", bg_lp_coefs = c(x_std = 0.6, y_std = 0),
#'   trig_lp_coefs = c(m_rel = 0.3),
#'   seed = 19
#' )
#'
#' fit <- fit_etas_pp(sim$catalog, formula.bg = ~ x)
#' fit
#' summary(fit)
#'
#' # Background effect, with its standard error
#' summary(fit$background_fit$gam)$p.table
#' }
#'
#' @family ETAS modelling
#' @export
fit_etas_pp <- function(
    catalog,
    formula.bg = ~ 1,
    covariates.bg = NULL,
    formula1 = "time ~ magnitude - 1",
    magn.threshold = 2.5,
    start = NULL,
    fixed = "gamma",
    onlytime = FALSE,
    background = list(),
    control = etas_control(),
    ...
) {
  call <- match.call()
  dots <- list(...)
  .pp_check_named_dots(dots, "fit_etas_pp()")

  control <- .pp_as_control(control, etas_control, "pp_etas_control", "control")

  if (is.null(start)) {
    start <- list()
  }

  if (!is.list(start)) {
    if (!is.numeric(start) || is.null(names(start))) {
      stop("`start` must be a named list or a named numeric vector.", call. = FALSE)
    }
    start <- as.list(start)
  }

  if (is.null(background)) {
    background <- list()
  }

  if (!is.list(background) ||
      (length(background) > 0L && is.null(names(background)))) {
    stop("`background` must be a named list of fit_pp() settings.", call. = FALSE)
  }

  # Translate arguments of the previous interface ------------------------------

  legacy_start <- intersect(names(dots), c(.pp_etas_parameter_names, "betacov"))
  legacy_bg <- intersect(names(dots), paste0(.pp_etas_background_settings, ".bg"))
  legacy_control <- intersect(names(dots), names(formals(etas_control)))
  legacy_inactive <- intersect(names(dots), .pp_etas_inactive_arguments)
  legacy_fixed <- intersect(names(dots), "params.ind")

  unknown <- setdiff(
    names(dots),
    c(legacy_start, legacy_bg, legacy_control, legacy_inactive, legacy_fixed)
  )

  if (length(unknown) > 0L) {
    stop(
      sprintf(
        paste0(
          "Unknown argument(s) of fit_etas_pp(): %s. Use `start`, `fixed`, ",
          "`background = list(...)` or `control = etas_control(...)`."
        ),
        paste(unknown, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  legacy <- c(legacy_start, legacy_fixed, legacy_bg, legacy_control, legacy_inactive)

  if (length(legacy) > 0L) {
    .pp_interface_warning(
      "fit_etas_pp()",
      legacy,
      "`start`, `fixed`, `background = list(...)` and `control = etas_control(...)`"
    )
  }

  start[legacy_start] <- dots[legacy_start]

  if (length(legacy_bg) > 0L) {
    background[sub("\\.bg$", "", legacy_bg)] <- dots[legacy_bg]
  }

  if (length(legacy_control) > 0L) {
    settings <- unclass(control)
    settings[legacy_control] <- dots[legacy_control]
    control <- do.call(etas_control, settings)
  }

  # Starting values and fixed parameters -----------------------------------------

  bad_start <- setdiff(names(start), c(.pp_etas_parameter_names, "betacov"))

  if (length(bad_start) > 0L) {
    stop(
      sprintf(
        "Unknown name(s) in `start`: %s. Allowed: %s.",
        paste(bad_start, collapse = ", "),
        paste(c(.pp_etas_parameter_names, "betacov"), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (length(legacy_fixed) > 0L) {
    params.ind <- dots$params.ind
  } else {
    if (is.null(fixed)) {
      fixed <- character(0)
    }

    if (!is.character(fixed) || !all(fixed %in% .pp_etas_parameter_names)) {
      stop(
        sprintf(
          "`fixed` must name ETAS parameters among: %s.",
          paste(.pp_etas_parameter_names, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    params.ind <- !(.pp_etas_parameter_names %in% fixed)
    names(params.ind) <- .pp_etas_parameter_names
  }

  # Background settings -----------------------------------------------------------

  bad_bg <- setdiff(names(background), .pp_etas_background_settings)

  if (length(bad_bg) > 0L) {
    stop(
      sprintf(
        "Unknown `background` setting(s): %s. Allowed: %s.",
        paste(bad_bg, collapse = ", "),
        paste(.pp_etas_background_settings, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  background_args <- background
  if (length(background_args) > 0L) {
    names(background_args) <- paste0(names(background_args), ".bg")
  }

  fit <- do.call(
    .pp_fit_etas_engine,
    c(
      list(
        catalog = catalog,
        formula.bg = formula.bg,
        covariates.bg = covariates.bg,
        formula1 = formula1,
        magn.threshold = magn.threshold,
        params.ind = params.ind,
        onlytime = onlytime
      ),
      start,
      background_args,
      unclass(control),
      dots[legacy_inactive]
    )
  )

  fit$call <- call
  fit
}
