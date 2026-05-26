# antedep 0.3.0

## New features
- Fit objects now inherit from a shared `ad_fit` parent class while retaining
  their family-specific classes (`gau_fit`, `cat_fit`, and `inad_fit`).
- Added `fit_ad()` as a unified fitting entry point for Gaussian, categorical,
  and INAD workflows.
- Added `anova()` methods for likelihood-ratio comparison of nested
  antedependence fits.

## Improvements
- `ci_inad()` now reports block-effect `tau` standard errors from the same
  profile-likelihood interval width used for the existing `tau` confidence
  intervals, keeping tau SE and CI columns aligned in downstream tables.
- `coef()` now returns a named numeric vector by default for standard R
  compatibility; use `type = "list"` for the previous structured inspection
  format.
- `confint()` now returns a standard two-column matrix for all fit families;
  the structured `ci_gau()`, `ci_cat()`, and `ci_inad()` helpers remain
  available.
- Added broader `vcov()` support using confidence-interval standard errors
  where available, improving compatibility with tools such as
  `lmtest::coeftest()`.
- The family-specific `aic_*()` and `bic_*()` helpers are retained as thin
  compatibility wrappers for existing code, but examples and new documentation
  prefer the standard `AIC()` and `BIC()` generics.

# antedep 0.2.0

## New features
- New S3 methods `deviance()`, `confint()`, and `vcov()` for `gau_fit`,
  `cat_fit`, and `inad_fit` objects, improving compatibility with standard
  R model-fitting workflows.
- New `summary.partial_corr()` method prints a lag-by-lag table of
  intervenor-adjusted partial correlations (mean absolute value, range, and
  number of significant pairs at each lag).

## Improvements
- Matrix inversion in `ci_inad()`, `fit_gau_em()`, and `cat_test_stats.R`
  now uses `chol2inv(chol(.))` instead of `solve()` for symmetric
  positive-definite matrices, improving numerical stability and efficiency.

# antedep 0.1.0

## New features
- `fit_inad()` gains a `nb_inno_size_ub` argument (default 50) that caps the
  upper bound of the negative-binomial innovation size parameter during
  optimization, improving numerical stability for near-Poisson data.
- `test_order_gau()` accepts `order_null` and `order_alt` as convenience aliases
  for `p` and the absolute alternative order; both are also returned in the
  result object.

## Bug fixes
- `ci_inad()`: fixed a sign error in the observed Fisher information for the
  negative-binomial innovation size parameter; the Hessian term
  `(r + u) / (r + λ)²` was added instead of subtracted, producing confidence
  intervals that were too wide.
- `ci_inad()`: the numerical second derivative for `nb_inno_size` CIs now
  retries with progressively smaller step sizes (×0.1, ×0.01) before falling
  back to NA, avoiding spurious failures when the default step lands in a
  non-finite region.
- `test_homogeneity_inad()`: degrees of freedom for LRT tests involving
  `innovation = "nbinom"` are now computed from the actual number of NB size
  parameters in the fitted models rather than assuming a fixed count of 1.
  This corrects LRT statistics and p-values whenever `nb_inno_size` is fitted
  as a time-varying vector.
- `ci_inad()` tau profile CI: `nb_inno_size` (negative-binomial innovation
  dispersion) is now held fixed at its full-model MLE during profile refits,
  consistent with the constrained-fit paradigm used throughout the package.
  Previously it was re-optimised as a nuisance parameter, which could widen
  the interval to the point of crossing zero even when the LRT clearly rejects
  the null (Variant 1 vs Variant 2 fix).
- `ci_inad()` tau profile CI: the bracket search in `.ci_tau_profile_inad`
  no longer imposes an artificial upper cap (`max(|tau_mle| + 1, 1)`) on the
  search range. The maximum bracket iterations are increased from 20 to 50
  and the initial step size is set to `max(0.1, |tau_mle| * 0.2)`, preventing
  the search from stalling for large or near-zero MLEs.

## Initial release
- Initial CRAN submission candidate for Gaussian AD, categorical AD, and INAD workflows.
- Added/expanded examples for key user-facing modeling functions (`fit_*`, `em_*`, `simulate_*`, `logL_*`).
- Refreshed missing-data notes and harmonized documentation links and metadata.
- `logL_gau()` default missing-data behavior is now `na_action = "fail"` (previously
  marginalization-first in earlier drafts). For missing inputs, pass
  `na_action = "marginalize"` or `na_action = "complete"` explicitly.
- Added packaged datasets `labor_force_cat` (categorical labor-force sequences) and
  `race_100km` (continuous 100km race split times).
