# antedep 0.1.0

## Bug fixes
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
