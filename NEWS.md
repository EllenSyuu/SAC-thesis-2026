# covarchimax 1.0.0

First public release. It accompanies the M.Sc. thesis *Semiparametric
Analysis of Precipitation Extremes in Southern Québec using Archimax
Copulas* (Ziyi Zhou, McGill University, 2026).

## Estimation

* `estimate_theta()`, `estimate_theta_pairs()`: moment estimator of the
  Clayton generator parameter (Chatelain, Fougères and Nešlehová, 2020),
  with a delete-one jackknife standard error computed in O(n^2).
* `pickands_at()`, `pickands_curve()`, `extremal_coefficient()`:
  Pickands-type, CFG-type and maximum-likelihood estimators of the Pickands
  dependence function in any dimension, given theta. The MLE always
  satisfies the Pickands bounds.
* `pickands_cov_mle()`: covariate-dependent Pickands dependence function,
  fitted pointwise on the simplex, with an identity link (bounds enforced
  for every observed covariate profile) or a logit link (bounds by
  construction).

## Inference and diagnostics

* `boot_years()`, `boot_ci()`: block bootstrap with per-replicate seeds and
  percentile intervals.
* `gev_coef_table()`, `gev_sandwich_vcov()`: Godambe sandwich standard
  errors for pooled (composite-likelihood) GEV fits.
* `tail_dependence_archimax()`, `tail_dependence_empirical()`: model-implied
  and empirical tail-dependence coefficients.

## Failures are reported, not hidden

* `pickands_cov_mle()` returns per-grid-point diagnostics (status, objective,
  number of feasible starts, boundary flag). Infeasible problems give `NA`
  plus a warning, not a starting value.
* `boot_years()` counts failed replicates and reports the first error.
* `gev_backward()` records failed fits in its selection path.
* `validate_daily()`, `validate_monthly()` stop with actionable messages.

## Simulation

* `rclayton_archimax()`: exact sampler for bivariate Clayton-Archimax
  copulas with a logistic Pickands function.
