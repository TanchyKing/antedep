# antedep 0.3.0

## New features
- `fit_inad()` gains a `constrain` argument for time-invariance constraints
  (`"alpha"`, `"alpha1"`, `"alpha2"`, `"theta"`, `"kappa"`, `"innovation"`,
  `"both"`, `"all"`). Constrained models are refitted by maximizing the
  observed-data log-likelihood jointly over every remaining free parameter,
  including negative binomial innovation sizes, and are finished with Newton
  steps using the exact score and observed information. `n_params`, `AIC`,
  and `BIC` count only the free parameters.
- `ci_inad()` now computes the observed information exactly by Louis' identity
  for all free parameters jointly: models of order 0, 1, and 2, with or without
  block effects, with or without time-invariance constraints, and including
  the negative binomial innovation size. `method = "hessian"` gives the
  numerical observed-data Hessian (fixed absolute steps) as a check. Intervals
  for positive parameters are formed on the log scale by default
  (`log_scale = TRUE`); the full covariance matrix is returned as `vcov`.
  Block effects receive Wald intervals in `tau`; for unconstrained fits,
  profile-likelihood intervals are also returned in `tau_profile`.
- `em_inad()` is now an exact EM algorithm for all nine thinning-innovation
  combinations. The E-step evaluates the full conditional distribution of the
  latent thinned counts, so the expected complete-data log-likelihood is exact
  even for negative binomial innovations, whose log-likelihood is not linear in
  the latent counts; negative binomial innovation sizes are now estimated.
  With block effects the innovation parameters are updated by conditional
  maximization (ECM). The returned fit includes `n_params`, `aic`, and `bic`.
- `test_stationarity_inad()` and `run_stationarity_tests_inad()` support
  constraints on negative binomial innovation sizes (`"kappa"`,
  `"innovation"`, `"all"`).
- Fit objects now inherit from a shared `ad_fit` parent class while retaining
  their family-specific classes (`gau_fit`, `cat_fit`, and `inad_fit`).
- Added `fit_ad()` as a unified fitting entry point for Gaussian, categorical,
  and INAD workflows.
- Added `anova()` methods for likelihood-ratio comparison of nested
  antedependence fits.

## Bug fixes
- Constrained fits in `test_stationarity_inad()` held negative binomial
  innovation sizes fixed at their unconstrained estimates instead of
  re-estimating them, giving incorrect likelihood ratio statistics and BIC
  values for negative binomial innovations.
- `.chi_bar_pvalue()` (used by `test_order_inad(use_chibar = TRUE)`) counted the
  point mass at zero as contributing 1 to the upper tail, inflating every
  chi-bar-square p-value by 2^-r.
- Parameter counts in the tables of `test_order_inad()`, `bic_order_inad()`,
  and `test_stationarity_inad()` now include one negative binomial innovation
  size per time point and the coefficient alpha_{2,1} of order-2 models, as in
  `fit_inad()`. For order 2, the `"alpha1"` constraint now includes
  alpha_{2,1} instead of fixing it at zero.
- `em_inad()` restricted the latent thinned counts to be no larger than the
  count being thinned, which is wrong for Poisson and negative binomial
  thinning, and capped their coefficients below 1.
- The Bell parameter is now obtained from the Bell mean to machine precision
  (previously `uniroot()` default tolerance), so Bell log-likelihoods with
  block effects are smooth in the parameters.
- Louis-based intervals in `ci_inad()` for order-2 models used only the lag-1
  coefficient, and Bell intervals with block effects omitted the chain rule
  from the innovation mean to the Bell parameter.

## Improvements
- Block-effect `tau` standard errors and confidence intervals are aligned in
  downstream tables (Wald in `tau`; profile likelihood in `tau_profile`).
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
