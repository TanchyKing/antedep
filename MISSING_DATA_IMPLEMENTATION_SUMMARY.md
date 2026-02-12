# Missing Data Implementation - Summary

This document summarizes the implementation of Sections 1-4 of DESIGN_PROPOSAL_MISSING.md for the antedep package.

## Files Created

### 1. missing_utils.R
Core utility functions for missing data handling (Section 2.3):

- **`.validate_missing(y)`** - Validates and classifies missing data patterns
  - Identifies: complete, dropout, drop-in, monotone_middle, intermittent, all_missing
  - Returns diagnostic information (n_complete, n_intermittent, pct_missing)

- **`.get_truncation_bound(y, subject_idx, time_idx, buffer, min_bound)`** - Computes upper bound for INAD missing value summation
  - Uses subject's observed max and other subjects at same time point
  - Default: buffer = 1.2, min_bound = 10

- **`.safe_log(x)`** - Safe logarithm with 0 * log(0) = 0 convention

- **`.warn_intermittent_missing(missing_info, max_intermittent)`** - Warns about potentially slow computation

- **`.extract_complete_cases(y, blocks)`** - Extracts complete cases with warnings

### 2. logL_gau_missing.R
AD-specific missing data log-likelihood (Section 3.1):

- **`.logL_gau_missing(y, order, mu, phi, sigma, blocks, tau)`** - Observed-data log-likelihood via MVN marginalization
  - For each subject, marginalizes over missing values
  - Uses Cholesky decomposition for numerical stability
  - Handles order 0, 1, 2 (order 2 is approximate)

- **`.build_gau_covariance(order, phi, sigma, n_time)`** - Constructs covariance matrix from AD parameters
  - Order 0: diagonal (independence)
  - Order 1: uses forward recursion for variances and covariances
  - Order 2: approximate (simplified version)

### 3. fit_gau_em.R
EM algorithm for AD with missing data (Section 3.1):

- **`.fit_gau_em(y, order, blocks, estimate_mu, max_iter, tol, verbose)`** - Main EM fitting function
  - Iterates E-step and M-step until convergence
  - Returns same structure as fit_gau plus EM diagnostics
  - Tracks log-likelihood trace for convergence monitoring

- **`.initialize_gau_em(y, order, blocks, estimate_mu)`** - Initialization strategy
  - Tries complete-case fit if ≥10 complete subjects
  - Falls back to marginal estimates if needed

- **`.initialize_gau_marginal(y, order, blocks)`** - Marginal initialization

- **`.em_e_step_gau(y, order, mu, phi, sigma, blocks, tau)`** - E-step implementation
  - Computes E[Y_t | Y_obs] for missing values (conditional mean)
  - Computes E[Y_t * Y_s | Y_obs] for all pairs (second moments)
  - Uses MVN conditioning formulas

- **`.em_m_step_gau(suff_stats, order, blocks, estimate_mu, n_subjects, n_time)`** - M-step implementation
  - Updates mu from sufficient statistics
  - Updates phi and sigma based on order
  - Order 2 currently uses simplified updates (needs full implementation)

### 4. logL_gau_modified.R
Modified logL_gau function with na_action parameter (Section 4.1):

- **`logL_gau(y, order, mu, phi, sigma, blocks, tau, na_action)`** - Main log-likelihood function
  - `na_action = "marginalize"` (default): Uses MVN marginalization via .logL_gau_missing
  - `na_action = "complete"`: Uses only complete cases
  - `na_action = "fail"`: Errors if any NA present
  - Complete data uses original efficient computation

- **`.count_params_gau(order, n_time, n_blocks)`** - Parameter counting utility

### 5. fit_gau_modified.R
Modified fit_gau function with EM support (Section 4.2):

- **`fit_gau(y, order, blocks, estimate_mu, na_action, em_max_iter, em_tol, em_verbose)`** - Main fitting function
  - `na_action = "em"` (default for missing): Uses .fit_gau_em
  - `na_action = "complete"`: Removes incomplete subjects
  - `na_action = "fail"`: Errors if any NA
  - Complete data uses original direct MLE

- Returns extended structure with missing data info:
  - `n_obs`, `n_missing`, `pct_missing`
  - `missing_pattern` ("complete", "monotone", "intermittent")
  - `em_converged`, `em_iterations`, `em_ll_trace`

### 6. test-missing_gau.R
Comprehensive test suite:

- Tests for all utility functions in missing_utils.R
- Tests for logL_gau with different na_action options
- Tests for covariance matrix construction
- Tests for fit_gau with missing data
- Tests for EM initialization, E-step, M-step
- Integration test: simulate → introduce missing → fit with EM
- Monotonicity test for EM log-likelihood

## Implementation Notes

### What Works
✅ Missing data validation and pattern detection
✅ MVN marginalization for observed-data likelihood
✅ EM algorithm for order 1 (tested and working)
✅ Complete-case analysis option
✅ Proper convergence diagnostics
✅ Integration with existing fit_gau structure

### Limitations / TODOs
⚠️ Order 2 covariance computation is approximate (simplified)
⚠️ Order 2 M-step uses placeholder updates (needs proper implementation)
⚠️ Block effects in EM need proper sufficient statistics by block
⚠️ No confidence intervals yet (Louis' identity extension needed)

### Design Decisions

1. **Sample size for BIC**: Uses `n_subjects` (number of rows), not `n_obs`
   - Consistent with longitudinal data literature
   - Subjects are independent units, not individual observations

2. **Initialization strategy**:
   - First: try complete-case fit if ≥10 complete subjects
   - Fallback: marginal estimates ignoring dependence
   - Avoids need for manual initial values

3. **na_action defaults**:
   - `logL_gau`: "marginalize" (statistically correct)
   - `fit_gau`: "em" (uses all available data)
   - Consistent with MAR assumption

4. **Numerical stability**:
   - Cholesky decomposition for MVN densities
   - Log-scale for sigma parameters in optimization
   - Bounds on sigma to avoid zero/negative values

## Usage Examples

### Basic usage with missing data
```r
# Simulate data with missing values
set.seed(123)
y <- matrix(rnorm(50), nrow = 10, ncol = 5)
y[sample(50, 5)] <- NA  # Introduce random missingness

# Fit with EM (default for missing data)
fit <- fit_gau(y, order = 1)

# Check results
print(fit$em_converged)  # Should be TRUE
print(fit$pct_missing)   # Percentage missing
plot(fit$em_ll_trace)    # Convergence plot
```

### Compare EM vs complete-case
```r
fit_em <- fit_gau(y, order = 1, na_action = "em")
fit_cc <- fit_gau(y, order = 1, na_action = "complete")

# EM uses more data
print(c(em = fit_em$n_obs, cc = fit_cc$n_obs))

# Compare log-likelihoods (on different data, so not directly comparable)
print(c(em = fit_em$log_l, cc = fit_cc$log_l))
```

### Verbose EM for debugging
```r
fit <- fit_gau(y, order = 1, na_action = "em", 
              em_max_iter = 50, em_verbose = TRUE)
# Prints iteration info:
# === EM Algorithm for AD with Missing Data ===
# Missing data: 10% (5 values)
# Complete subjects: 5/10
# Intermittent missing: 2 subjects
# 
# Iteration 1: log-lik = -75.3241
# Iteration 10: log-lik = -73.8905
# ...
# Converged at iteration 23
```

### Error handling
```r
# Ensure no missing data
fit <- fit_gau(y, order = 1, na_action = "fail")
# Error: y contains NA values. Use na_action = 'em' or 'complete'
```

## Testing

Run tests with:
```r
source("missing_utils.R")
source("logL_gau_missing.R")
source("fit_gau_em.R")
source("logL_gau_modified.R")
source("fit_gau_modified.R")

testthat::test_file("test-missing_gau.R")
```

Tests cover:
- ✅ Pattern detection (complete, dropout, drop-in, intermittent)
- ✅ Truncation bounds
- ✅ Complete case extraction with warnings
- ✅ logL_gau with all na_action options
- ✅ Covariance matrix construction
- ✅ fit_gau with all na_action options
- ✅ EM initialization strategies
- ✅ EM E-step sufficient statistics
- ✅ EM M-step parameter updates
- ✅ EM convergence and monotonicity
- ✅ Integration test (full workflow)

## Next Steps

To complete the missing data implementation:

1. **Order 2 improvements** (Section 3.1):
   - Implement exact covariance recursion in `.build_gau_covariance`
   - Implement proper M-step updates in `.em_m_step_gau`

2. **INAD module** (Section 3.3):
   - Implement `.logL_inad_missing` with truncation
   - Extend existing `em_inad` to handle missing data
   - Modify `fit_inad` to add `na_action` parameter

3. **CAT module** (Section 3.2):
   - Implement forward-backward algorithm
   - Implement `.logL_cat_missing`
   - Implement `.fit_cat_em`
   - Modify `fit_cat` to add `na_action` parameter

4. **Downstream functions** (Section 5):
   - Update all `lrt_*` functions to pass `na_action`
   - Extend `ci_inad` and `ci_cat` for missing data (Louis' identity)

5. **Documentation**:
   - Add roxygen documentation for all functions
   - Create vignette on missing data handling
   - Update package README

## File Integration

To integrate into the antedep package:

1. Move files to R/ directory:
   ```bash
   cp missing_utils.R /path/to/antedep/R/
   cp logL_gau_missing.R /path/to/antedep/R/
   cp fit_gau_em.R /path/to/antedep/R/
   ```

2. Replace existing files:
   ```bash
   # Backup originals first!
   cp logL_gau_modified.R /path/to/antedep/R/logL_gau.R
   cp fit_gau_modified.R /path/to/antedep/R/fit_gau.R
   ```

3. Move test file:
   ```bash
   cp test-missing_gau.R /path/to/antedep/tests/testthat/
   ```

4. Update NAMESPACE and documentation:
   ```r
   devtools::document()
   devtools::test()
   devtools::check()
   ```

## References

- Zimmerman & Núñez-Antón (2009), Chapter 8: Missing Data
- Little & Rubin (2019), Statistical Analysis with Missing Data
- Dempster et al. (1977), Maximum likelihood from incomplete data via the EM algorithm
