# Missing Data Support (Current State)

This note reflects the integrated package behavior on the `refactor/remove-ad-main` line.
It replaces earlier implementation-draft instructions that referenced temporary files.

## Current API Behavior

### Gaussian (`gau`)

- `logL_gau()` supports `na_action = c("marginalize", "complete", "fail")`.
- `logL_gau()` default is `"marginalize"`.
- `fit_gau()` supports `na_action = c("fail", "complete", "em")`.
- `fit_gau()` default is `"fail"`.
- `fit_gau(na_action = "em")` is the missing-data estimation path.
- For `order = 2`, `fit_gau(..., na_action = "em")` is available but documented as provisional.

### Categorical (`cat`)

- `logL_cat()` supports `na_action = c("fail", "complete", "marginalize")`.
- `fit_cat()` supports `na_action = c("fail", "complete", "marginalize", "em")`.
- `fit_cat(..., na_action = "em")` currently supports `order` 0 and 1.
- For `order = 2`, use `na_action = "marginalize"` (EM path errors by design).

### Integer-valued (`inad`)

- `logL_inad()` supports `na_action = c("fail", "complete", "marginalize")`.
- `fit_inad()` supports `na_action = c("fail", "complete", "marginalize")`.

## Usage

Use package functions directly after loading/installing the package.
Do not `source()` individual files from `R/`.

```r
library(antedep)

# Gaussian observed-data likelihood (default is marginalize)
ll <- logL_gau(y, order = 1, mu = mu, phi = phi, sigma = sigma)

# Gaussian fitting with EM for missing values
fit <- fit_gau(y, order = 1, na_action = "em")

# Categorical and INAD fitting with marginalization
fit_cat_obj <- fit_cat(y_cat, order = 1, na_action = "marginalize")
fit_inad_obj <- fit_inad(y_inad, order = 1, na_action = "marginalize")
```

## Tests

Targeted missing-data tests live under `tests/testthat/`, including:

- `test-logL_gau.R`
- `test-missing_gau.R`
- `test-loglik_inad.R`
- `test-inference-missing-data-guards.R`

## Notes

If behavior and this file diverge, treat function documentation and source as canonical and update this file.