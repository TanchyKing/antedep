# antedep

`antedep` fits antedependence models for longitudinal data:

- Gaussian AD models for continuous outcomes
- INAD models for count outcomes
- Categorical AD (CAT) models for discrete state outcomes

## Installation

```r
# install.packages("remotes")
remotes::install_github("TanchyKing/antedep")
```

## Production-Readiness Matrix

| Model | Data type | Complete-data fit/logLik | Missing-data fit/logLik | Notes |
|---|---|---|---|---|
| AD | Continuous | Ready | Ready (`fit_ad`, `logL_ad`) | Missing-data fit uses EM or observed-data likelihood modes |
| INAD | Counts | Ready | Ready (`fit_inad`, `logL_inad`) | Missing-data fit supports `na_action = "marginalize"` |
| CAT | Categorical states | Ready | Ready (`fit_cat`, `logL_cat`) | Missing-data fit supports orders 0, 1, 2 |

## Quick Start

### AD example (continuous)

```r
library(antedep)
set.seed(1)

y <- simulate_ad(n_subjects = 80, n_time = 5, order = 1)
fit <- fit_ad(y, order = 1)
fit$log_l
```

### Missing-data workflow

```r
library(antedep)
set.seed(1)

y <- simulate_inad(n_subjects = 60, n_time = 5, order = 1)
y[sample(length(y), 20)] <- NA

# Fit observed-data likelihood under MAR
fit_miss <- fit_inad(y, order = 1, na_action = "marginalize")
fit_miss$log_l
```

## Known Limitations

- Missing-data confidence intervals and likelihood-ratio inference are currently not implemented (`ci_*`, `lrt_*`); functions fail with clear error messages on incomplete data.
- Current model order support is up to order 2 for AD, INAD, and CAT.
- CAT missing-data marginalization currently supports orders 0, 1, and 2.

## Vignette

Source vignette: `vignettes/antedep-intro.Rmd`.
Rendered site article: `docs/articles/antedep-intro.html`.
Integrated function reference: `docs/reference/index.html`.

## Function Reference Site (pkgdown)

- Local build: `Rscript -e "pkgdown::build_site(preview = FALSE)"`
- CI deployment: `.github/workflows/pkgdown.yml` (GitHub Pages via Actions)

## Local Check Guide

- See `LOCAL_CHECKS.md` for tarball-mode vs directory-mode checks and expected NOTE messages.
