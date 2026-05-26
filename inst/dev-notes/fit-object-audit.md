# Fit Object Audit

This note records the shared and family-specific fields for the package fit
objects after the v0.3.0 interface cleanup.

## Shared Structure

All fitted model objects inherit from the common parent class `ad_fit` in
addition to their family-specific class:

- `c("gau_fit", "ad_fit")`
- `c("cat_fit", "ad_fit")`
- `c("inad_fit", "ad_fit")`

The standard methods rely on these shared fields:

- `log_l`: maximized log-likelihood
- `n_params`: number of free parameters
- `aic`: Akaike information criterion when stored by the fitter
- `bic`: Bayesian information criterion when stored by the fitter
- `settings`: list containing at least `order`, `n_subjects`, and `n_time`
- `convergence`: convergence metadata when available

## Family-Specific Fields

`gau_fit` stores Gaussian parameters:

- `mu`: marginal mean profile
- `sigma`: innovation standard deviations
- `phi`: AD coefficients
- `tau`: optional block effects

`cat_fit` stores categorical probabilities:

- `marginal`: marginal or initial conditional probabilities
- `transition`: transition probability arrays
- `cell_counts`: observed or expected cell counts used for standard errors

`inad_fit` stores count-model parameters:

- `alpha`: thinning parameters
- `theta`: innovation parameters
- `tau`: optional block effects
- `nb_inno_size`: optional negative-binomial innovation size parameters

## Standard Method Decisions

- `coef()` returns a named numeric vector by default and accepts
  `type = "list"` for family-specific inspection.
- `confint()` returns a standard matrix with `lower` and `upper` columns.
  The structured `ci_gau()`, `ci_cat()`, and `ci_inad()` helpers are retained.
- `confint.inad_fit()` requires the original data matrix through `y` because
  INAD confidence intervals currently recompute data-dependent observed
  information/profile quantities that are not stored on the fit object.
- INAD block-effect `tau` intervals are profile-likelihood intervals, and the
  displayed `tau` standard errors are derived from those profile interval
  widths. The remaining INAD confidence-interval rows use their existing
  Wald-style standard errors and intervals.
- `vcov()` returns a named diagonal matrix using available standard errors.
  Entries are `NA` when a standard error is not available. Off-diagonal
  covariances are not estimated, so downstream Wald summaries should be
  interpreted as diagonal/independence approximations unless the relevant
  model block is known to have negligible covariance.
- `vcov.inad_fit()` also requires the original data matrix through `y` for the
  same reason as `confint.inad_fit()`.
- `anova()` compares nested fits by likelihood-ratio statistics in the order
  supplied by the user.

## AIC/BIC Helper Policy

The family-specific helpers `aic_gau()`, `aic_cat()`, `aic_inad()`,
`bic_gau()`, `bic_cat()`, and `bic_inad()` are retained as thin compatibility
wrappers for existing user code. They are not deprecated in v0.3.0, but new
examples and manuscript/replication material should prefer the standard
`AIC()` and `BIC()` generics because the fit classes now provide `logLik()`
methods.
