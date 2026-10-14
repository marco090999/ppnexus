# ppnexus 0.1.0

## Simpler interfaces for the LGCP and ETAS extensions

* `fit_lgcp_pp()` now has nine arguments. Settings of the first-order fit
  (`nc`, `q`, `covariate_methods`, ...) are passed to `fit_pp()` through `...`,
  and the numerical settings of the minimum-contrast stage are collected by
  the new `lgcp_control()`. A custom covariance is given as a function in
  `cov_model`.
* `fit_etas_pp()` now has ten arguments: starting values in `start`, fixed
  parameters in `fixed`, further background settings in `background`, and
  algorithmic settings in the new `etas_control()`. Arguments inherited from
  `etasFLP::etasclass()` that had no effect have been removed from the
  interface. Calls using the previous arguments still work, with a warning
  that shows the new syntax. Progress messages are now off by default.
* `fit_etas_pp()` warns when standard errors cannot be computed because the
  Hessian is singular (typically a catalogue with little triggering).

## Further bug fixes

* `fit_lgcp_pp()` with a custom covariance function returned the starting
  values: with unbounded parameters L-BFGS-B stopped at the first step.
  Nelder-Mead is now used when no bound is finite.
* `simulate_lgcp_pp()`, `simulate()` for fitted LGCPs and the LGCP
  diagnostics no longer emit the grid-bound warning of `simulate_pp()`, which
  is not actionable for a random latent field (violations still trigger an
  automatic restart).

## Graphics and output

* Global envelopes of residual profiles are drawn as bands, one figure per
  variable, with the variable and residual type on the axes. K-function
  envelopes are labelled accordingly.
* Shorter panel titles in `plot()` for marked patterns; ETAS catalogue
  legends have a background box; missing axis units are printed as blanks.

## Examples

* All examples run without warnings. The `fit_etas_pp()`, `decluster_etas_pp()`
  and `pp_etas_fit` examples use an informative simulated catalogue, and the
  generic covariate name `cov1` replaces `temperature`.

## Testing and documentation

* Added a `testthat` (3rd edition) test suite. Tests check numerical results
  against closed-form values (intercept-only and categorical-mark fits,
  translation-corrected K-function, FPCA identities, LGCP K-function),
  parameter recovery from simulated data, and the compiled ETAS triggering
  routine against a plain R implementation.
* Added the vignette "Getting started with ppnexus".
* Added a `print()` method for `pp_window` objects.

## Compiled code

* The ETAS triggering intensity is now computed by `src/etas_trigger.f90`,
  adapted from the Fortran routine `etasfull8newserial` of etasFLP
  (GPL >= 2). ppnexus no longer depends on etasFLP and is now licensed
  under GPL (>= 2).

## Bug fixes

* `fit_pp()` no longer emits one `non-integer x` warning per quadrature row
  when cubature weights have non-integer reciprocals (any window other than
  special cases such as the unit square). The Poisson family used internally
  now evaluates `aic()` and `ls()` with the continuous extension of the
  Poisson log-probability; estimates are unchanged.
* `fit_pp()` no longer warns that `method` and `args` are ignored for joint
  continuous marks when neither was supplied.
* `simulate_etas_pp()` now returns a valid catalogue with `exploded = TRUE`
  when `max_events` is exceeded, instead of failing an internal validation
  check.
* `plot()` for `pp_envelope` objects now draws the global envelope; the
  ggplot object returned by GET was previously not printed.
* The colour bar of `plot()` for `pp_residuals` objects no longer fails with
  "figure margins too large" on ordinary device sizes.
