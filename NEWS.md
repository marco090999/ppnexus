# ppnexus 0.1.0

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
