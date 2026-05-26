# Man Page Title Audit

This note tracks the documentation title audit requested after the JSS
editorial comments.

## Procedure

Run the following after regenerating documentation:

```r
help(package = "antedep")
```

Then check that visible topic titles use consistent title case and that roxygen
`@title` fields match the generated `.Rd` titles.

## Current Snapshot

Captured from `help(package = "antedep")` after building and checking
`antedep_0.3.0.tar.gz`:

```text
Bell                    The Bell distribution
anova.antedep           Likelihood-Ratio Comparisons for Nested AD Fits
as.ts.antedep_sim       Convert Antedependence Simulation Output to a Time-Series Object
bic_cat                 Bayesian and Akaike Information Criteria for Antedependence Fits
bic_order_cat           BIC-Based Order Selection for Categorical AD Models
bic_order_gau           BIC-Based Order Selection for Gaussian AD Models
bic_order_inad          BIC-Based Order Selection for INAD Models
bolus_inad              Morphine Bolus Analgesia Counts
cattle_growth           Cattle Growth Data (Treatments A and B)
ci_cat                  Confidence Intervals for Fitted Categorical AD Models
ci_gau                  Confidence Intervals for Fitted Gaussian AD Models
ci_inad                 Confidence Intervals for Fitted INAD Models
coef.antedep            Extract Model Coefficients from Antedependence Fits
confint.antedep         Confidence Intervals for Antedependence Model Fits
deviance.antedep        Deviance for Antedependence Model Fits
fit_ad                  Unified Antedependence Model Fitting Interface
logLik.antedep          Log-Likelihood for Antedependence Model Fits
nobs.antedep            Number of Subjects for Antedependence Model Fits
vcov.antedep            Variance-Covariance Matrices for Antedependence Model Fits
```

The full help index was checked during `R CMD check --as-cran --no-manual`,
which completed with `Status: OK`.

## Current Decision

The v0.3.0 documentation pass keeps standard-method pages grouped under shared
titles such as:

- `Log-likelihood for antedependence model fits`
- `Extract model coefficients from antedependence fits`
- `Confidence intervals for antedependence model fits`
- `Variance-covariance matrices for antedependence model fits`
- `Likelihood-ratio comparisons for nested AD fits`

Family-specific helper pages are retained for detailed APIs such as
`ci_gau()`, `ci_cat()`, `ci_inad()`, `fit_gau()`, `fit_cat()`, and `fit_inad()`.

## Follow-up

The v0.3.0 standard-method pages were updated to title case. If a future venue
requires strict title case for every print and summary helper as well, audit
those helper pages separately because many R packages intentionally use
descriptive method titles for those entries.
