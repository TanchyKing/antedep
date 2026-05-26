# INAD Validation Summary v2 - Prediction and Marginal Moments

This v2 note implements the revision plan in `inad-validation-revision-todo.md`.
The main change from v1 is scope discipline: the simulation section is
now INAD-only correct-specification validation, and the bolus predictive
section is restricted to rolling one-step forecasting. Recursive multi-step
tables are retained only as smoothing-sensitivity context.

## 1. Scope and setup

### 1.1 What this note covers

- Simulation: absolute prediction precision for constrained-alpha INAD under its own NBT-NBI-INADFE(1) DGP.
- Bolus in-sample fit: BIC comparison across INAD and non-INAD count models.
- Bolus prediction: rolling one-step leave-one-per-group-out prediction on all 30 x 35 = 1050 held-out pairs.
- Bolus marginal moments: Henderson-and-Shimakura-style time-by-time mean and variance comparison.
- Package implications: what already landed in `antedep`, and what remains.

CGFM is included only in the marginal-moment comparison, using published Henderson-Shimakura estimates. No R implementation is available for refitting CGFM inside each cross-validation fold, so CGFM is not included in the predictive comparison.

### 1.2 Evaluation framework

The predictive target is conditional new-subject forecasting. A model is fit on training subjects, early observations of held-out subjects are revealed, and later observations are predicted from the fitted conditional distribution. The v2 headline uses rolling one-step prediction at t = 2, ..., 12, conditioning on the realized history through t - 1.

All predictive comparisons are plug-in MLE conditional. INAD uses Monte Carlo forward simulation from the fitted process; NB and Poisson GLMs use analytic plug-in PMFs; NB GLMM uses the fitted random-intercept model for prediction; tscount uses a pooled log-linear count time-series fit. Lower RMSE, lower log score, and lower RPS are better.

### 1.3 Metric primer

- RMSE measures point prediction error of the predictive mean, on the original count scale.
- Log score is `-log p(Y_obs)`. It rewards sharp probability at the realized count but is sensitive to finite-support smoothing when the PMF is simulation-based.
- RPS compares the predictive CDF with the realized count indicator over the full integer support. It is a proper probabilistic score and is less sensitive than log score to a single tail probability.

## 2. Models compared - formulations

Notation: subject i has group g(i), time t = 1, ..., T, and group 1 is the reference.

### 2.1 NB GLM

`Y_it ~ NB(mu_it, theta)`, with `log mu_it = beta0 + beta_group[g(i)] + beta_time[t]`. This is the marginal no-lag form used for the BIC and marginal-moment comparison.

### 2.2 NB GLMM

Two forms are used. The predictive form is lag-aware: `log mu_it = beta0 + beta_group[g(i)] + beta_time[t] + gamma log(Y_i,t-1 + 1) + b_i`, with `b_i ~ N(0, sigma_b^2)`. The marginal-moment/BIC form drops the lag term so closed-form marginal moments can be obtained by integrating the random intercept.

### 2.3 Poisson GLM

The predictive pipeline uses the lag-aware form `Y_it ~ Poisson(mu_it)`, `log mu_it = beta0 + beta_group[g(i)] + beta_time[t] + gamma log(Y_i,t-1 + 1)`. Rows at t = 1 are dropped because no lag is available.

### 2.4 tscount

The implemented comparator is an exploratory pooled `tscount::tsglm` fit. Training subjects are concatenated into one artificial sequence, with `model = list(past_obs = 1)` and `xreg = model.matrix(~ group + time_fac, ...)`. This introduces subject-boundary lag artifacts, so the BIC row is labelled non-comparable and the model is treated as a surrogate count-time-series baseline rather than a clean panel likelihood.

### 2.5 INAD constrained-alpha

NBT-NBI-INADFE(1): `Y_it = alpha o Y_i,t-1 + epsilon_it`, with constant alpha, negative-binomial thinning, and NB innovations with time-varying mean/size plus group effect tau. This is the BIC-supported primary INAD variant for bolus.

### 2.6 INAD unconstrained

Same model class, but alpha is allowed to vary over time. It is reported as a sensitivity because the bolus LRT detects mild non-stationarity even though BIC prefers the constrained model.

### 2.7 CGFM independent/shared frailty

CGFM rows use Henderson-Shimakura published full-data estimates only. In the local notation, the gamma frailty has mean 1 and frailty variance psi. Independent frailty has time-specific independent frailties; shared frailty uses a single subject-level frailty across all times. These rows are used for marginal mean and variance only.

## 3. Simulation - INAD prediction precision

DGP: NBT-NBI-INADFE(1), constant alpha = 0.35, time-varying NB innovation parameters, tau = c(0, 1.25), 100 subjects, 12 time points, 1000 simulation replicates. Because the DGP is INAD, this section is correct-specification validation rather than a neutral benchmark; cross-model baseline comparisons are not retained here.

### 3.1 Parameter recovery

| Fit | alpha MAE | theta RMSE | NB size RMSE | |tau2 bias| |
| --- | --- | --- | --- | --- |
| INAD constrained-alpha | 0.026 | 0.419 | 1.598 | 0.197 |
| INAD unconstrained | 0.089 | 0.576 | 1.598 | 0.197 |

The constrained-alpha fit improves alpha and theta recovery, as expected under a constant-alpha DGP. NB innovation size is shared by construction in the current constrained-alpha POC.

### 3.2 Rolling one-step precision

| Metric | Mean | SE | 2 5% | 97 5% | Reps |
| --- | --- | --- | --- | --- | --- |
| log_score | 2.508 | 0.002 | 2.403 | 2.611 | 1000 |
| rps | 1.857 | 0.003 | 1.653 | 2.069 | 1000 |
| rmse | 3.516 | 0.007 | 3.068 | 4.032 | 1000 |
| pit | 0.494 | 0.001 | 0.456 | 0.531 | 1000 |
| cover80 | 0.853 | 0.001 | 0.812 | 0.891 | 1000 |
| cover95 | 0.967 | 0.000 | 0.942 | 0.985 | 1000 |

PIT is close to 0.5 and the 80%/95% empirical coverages are slightly conservative, which is acceptable for the simulation validation target.

## 4. Bolus in-sample fit - BIC comparison

BIC is computed on the full bolus data. The tscount row has the best numeric BIC but is not directly comparable because it is a pooled artificial sequence rather than a proper subject-level panel likelihood. Among comparable rows, constrained-alpha INAD has the lowest BIC, with NB GLMM close behind.

| Model | Specification | BIC | Comparable | Note |
| --- | --- | --- | --- | --- |
| tscount tsglm | exploratory pooled log-linear NB tsglm | 3894.7 | no | Pooled artificial sequence; subject-boundary lag artifacts; not directly comparable to a proper panel likelihood. |
| INAD constrained-alpha | NBT-NBI-INADFE(1), constant alpha | 4243.5 | yes | BIC uses n_subjects; constrained POC inherits nb_inno_size from the unconstrained fit. |
| NB GLMM | marginal-moment no-lag: y ~ group + time_fac + (1 | subject) | 4264.8 | yes | No-lag form used for in-sample marginal fit and closed-form marginal moments. |
| INAD unconstrained | NBT-NBI-INADFE(1), time-varying alpha | 4297.8 | yes | BIC uses n_subjects, matching antedep fit-object convention. |
| Poisson GLM | lag-aware predictive: y ~ group + factor(time) + log(y_lag1 + 1) | 4334.6 | yes | Lag-aware form used in the predictive pipeline; asymmetric with marginal NB GLM. |
| CGFM independent frailty | published Henderson-Shimakura independent frailty | 4442.0 | yes | Reconstructed from published full log-likelihood and 12 time effects + group + frailty variance. |
| NB GLM | marginal no-lag: y ~ group + factor(time) | 4476.8 | yes | Marginal no-lag form, aligned with the dissertation BIC comparison. |
| CGFM shared frailty | published Henderson-Shimakura shared frailty | 4556.8 | yes | Reconstructed from published full log-likelihood and 12 time effects + group + frailty variance. |

## 5. Bolus prediction - rolling one-step

Design: all 30 x 35 cross-group held-out patient pairs, each fold holding out one group-1 and one group-2 patient and fitting on the remaining 63 patients. Standard errors are patient-level clustered, so overlapping folds do not receive naive independent-fold uncertainty.

### 5.1 Stationarity and primary-model choice

The full-data LRT for constant alpha rejects (p = 0.000234), but BIC strongly prefers constrained-alpha INAD over unconstrained INAD (4244 vs 4298). The interpretation is mild detectable non-stationarity that is not worth 11 extra alpha parameters under BIC. Constrained-alpha INAD is therefore primary; unconstrained INAD is sensitivity.

### 5.2 Constrained-alpha INAD vs baselines

Deltas are INAD minus reference; negative is INAD better.

| Reference | Delta RPS (t) | Delta RMSE (t) | Delta log score (t) |
| --- | --- | --- | --- |
| NB GLM | -0.396 (-22.3) | -0.689 (-23.9) | -0.123 (-18.6) |
| NB GLMM | -0.153 (-12.3) | -0.217 (-10.7) | -0.056 (-10.8) |
| Poisson GLM | -0.078 (-16.3) | -0.012 (-1.9) | -0.175 (-30.6) |
| tscount | -0.041 (-8.6) | -0.105 (-14.0) | 0.002 (0.6) |

### 5.3 Unconstrained INAD vs baselines

| Reference | Delta RPS (t) | Delta RMSE (t) | Delta log score (t) |
| --- | --- | --- | --- |
| NB GLM | -0.386 (-23.0) | -0.659 (-24.6) | -0.118 (-18.9) |
| NB GLMM | -0.144 (-12.3) | -0.187 (-10.1) | -0.052 (-10.4) |
| Poisson GLM | -0.069 (-14.0) | 0.018 (2.9) | -0.170 (-29.3) |
| tscount | -0.031 (-7.2) | -0.075 (-8.8) | 0.007 (2.3) |

### 5.4 Smoothing sensitivity context

The v2 headline is rolling one-step, but the existing recursive h = 2 smoothing sweep is retained as methodological context. It shows why log score is treated more cautiously than RPS: changing epsilon changes the single-point tail mass more than it changes integrated scores.

| Epsilon | Delta log score vs NB GLM (t) | Delta RPS vs NB GLM (t) |
| --- | --- | --- |
| 1/101 | -0.000 (-0.0) | -0.024 (-0.5) |
| 1/1001 | 0.047 (1.7) | 0.062 (1.7) |
| 1/10001 | 0.068 (2.0) | 0.065 (1.9) |

### 5.5 Bolus prediction headline

On rolling one-step prediction, constrained-alpha INAD beats the marginal NB GLM, NB GLMM, Poisson GLM, and tscount on RPS. It also beats all four on RMSE except that the Poisson RMSE margin is small and only borderline by t-statistic. Log score is better than NB GLM, NB GLMM, and Poisson GLM, and essentially tied with tscount. The unconstrained sensitivity tells the same qualitative story.

## 6. Bolus marginal mean and variance comparison

This section compares full-data fitted marginal moments with empirical bolus moments. CGFM is included here only, using published Henderson-Shimakura estimates. Relative discrepancy is `100 * (fitted - empirical) / empirical`; MARD is the mean absolute relative discrepancy.

### 6.1 MARD summary

| Model | G1 mean | G1 variance | G2 mean | G2 variance |
| --- | --- | --- | --- | --- |
| INAD constrained alpha | 8.35 | 26.74 | 9.31 | 22.79 |
| INAD unconstrained | 5.92 | 29.05 | 5.83 | 19.72 |
| NB GLM marginal | 8.89 | 42.36 | 6.94 | 21.01 |
| NB GLMM marginal | 8.35 | 58.35 | 8.24 | 24.25 |
| tscount tsglm MC | 12.95 | 61.03 | 10.12 | 27.82 |
| CGFM independent frailty | 8.94 | 42.44 | 6.87 | 20.75 |
| CGFM shared frailty | 9.86 | 29.01 | 6.01 | 41.83 |

### 6.2 G1 variance sensitivity

G1 time 11 has a very small empirical variance, so it can dominate relative-error summaries. The ranking is robust when it is removed.

| Model | MARD all 12 | MARD excl  t=11 |
| --- | --- | --- |
| INAD constrained alpha | 26.74 | 15.65 |
| INAD unconstrained | 29.05 | 20.37 |
| NB GLM marginal | 42.36 | 29.70 |
| NB GLMM marginal | 58.35 | 43.73 |
| tscount tsglm MC | 61.03 | 47.98 |
| CGFM independent frailty | 42.44 | 29.53 |
| CGFM shared frailty | 29.01 | 23.14 |

### 6.3 Group 1 mean relative discrepancy (%)

| Time | Empirical | C-alpha INAD | U-INAD | NB GLM | NB GLMM | tscount | CGFM ind  | CGFM shared |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
|  1 | 9.30 | -1.73 | -1.57 | -12.28 | -4.37 | 0.95 | -12.19 | -13.07 |
|  2 | 5.50 | 6.34 | -4.45 | -9.03 | -4.43 | 37.18 | -9.04 | -9.04 |
|  3 | 5.43 | 12.34 | 4.46 | 1.60 | 2.83 | 20.31 | 1.76 | 2.78 |
|  4 | 5.13 | -0.28 | 11.44 | 15.98 | 18.44 | 27.46 | 15.52 | 17.85 |
|  5 | 7.57 | -5.78 | -2.53 | -5.26 | -1.62 | -2.33 | -6.18 | -5.23 |
|  6 | 5.33 | -1.18 | -1.98 | -0.94 | 3.13 | 0.58 | -1.39 | -0.40 |
|  7 | 3.93 | 12.22 | 7.52 | 10.89 | 12.22 | 9.70 | 10.57 | 12.81 |
|  8 | 3.73 | 5.83 | 0.46 | 5.32 | 8.97 | 5.11 | 5.41 | 6.47 |
|  9 | 4.60 | -11.77 | -15.28 | -13.54 | -10.44 | -16.46 | -13.59 | -14.45 |
| 10 | 4.93 | 2.20 | -9.66 | -4.79 | -3.01 | -6.15 | -4.50 | -4.50 |
| 11 | 3.50 | 17.93 | 6.60 | 15.32 | 18.67 | 15.20 | 15.86 | 18.20 |
| 12 | 3.47 | -22.63 | -5.11 | 11.72 | 12.10 | 13.91 | 11.27 | 13.52 |

### 6.4 Group 1 variance relative discrepancy (%)

| Time | Empirical | C-alpha INAD | U-INAD | NB GLM | NB GLMM | tscount | CGFM ind  | CGFM shared |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
|  1 | 33.46 | 24.05 | 24.41 | 21.58 | 55.95 | -2.43 | 22.07 | -28.95 |
|  2 | 22.05 | -8.25 | -20.91 | -21.83 | -7.59 | 50.54 | -21.70 | -50.07 |
|  3 | 15.98 | 43.98 | 35.55 | 27.76 | 41.89 | 91.07 | 28.35 | -18.20 |
|  4 | 17.02 | 9.67 | 37.65 | 36.80 | 54.62 | 80.48 | 36.11 | -12.83 |
|  5 | 25.08 | -2.95 | 2.17 | 28.71 | 50.60 | 51.04 | 26.77 | -22.21 |
|  6 | 11.26 | 31.56 | 27.58 | 68.00 | 95.81 | 88.71 | 67.01 | 7.28 |
|  7 | 17.58 | -29.57 | -31.94 | -22.31 | -14.25 | -12.16 | -22.54 | -47.89 |
|  8 | 9.24 | 10.39 | 2.04 | 24.36 | 42.11 | 36.79 | 24.76 | -15.91 |
|  9 | 9.97 | 6.80 | 1.08 | 17.40 | 34.46 | 23.45 | 17.49 | -23.27 |
| 10 | 13.17 | -0.30 | 0.03 | 17.54 | 31.48 | 21.18 | 18.38 | -23.76 |
| 11 | 4.26 | 148.75 | 124.57 | 181.71 | 219.24 | 204.65 | 184.44 | 93.61 |
| 12 | 7.98 | -4.61 | 40.73 | 40.38 | 52.23 | 69.89 | 39.68 | -4.13 |

### 6.5 Group 2 mean relative discrepancy (%)

| Time | Empirical | C-alpha INAD | U-INAD | NB GLM | NB GLMM | tscount | CGFM ind  | CGFM shared |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
|  1 | 10.17 | 1.08 | 1.53 | 12.75 | 16.01 | -0.30 | 12.80 | 10.56 |
|  2 | 6.49 | 16.07 | 5.07 | 8.45 | 7.53 | 45.12 | 8.37 | 7.29 |
|  3 | 7.86 | 2.27 | -4.58 | -1.23 | -5.65 | 16.01 | -1.14 | -1.14 |
|  4 | 9.29 | -22.79 | -11.99 | -9.86 | -13.13 | -0.06 | -10.28 | -9.38 |
|  5 | 9.63 | -4.08 | 0.28 | 4.66 | 2.58 | 11.23 | 3.59 | 3.59 |
|  6 | 7.37 | 0.43 | 0.38 | 0.75 | -1.00 | 0.27 | 0.24 | 0.24 |
|  7 | 6.60 | -0.63 | -6.34 | -7.10 | -11.27 | -6.98 | -7.42 | -6.49 |
|  8 | 5.74 | 6.25 | 0.24 | -3.75 | -6.01 | -3.24 | -3.72 | -3.72 |
|  9 | 4.91 | 26.40 | 20.27 | 13.77 | 11.23 | 14.97 | 13.64 | 11.39 |
| 10 | 6.34 | 13.46 | 11.91 | 4.10 | 0.09 | 5.24 | 4.36 | 3.32 |
| 11 | 6.26 | 0.41 | -4.82 | -9.32 | -11.93 | -8.30 | -8.95 | -8.03 |
| 12 | 5.89 | -17.81 | 2.57 | -7.49 | -12.40 | -9.68 | -7.92 | -7.00 |

### 6.6 Group 2 variance relative discrepancy (%)

| Time | Empirical | C-alpha INAD | U-INAD | NB GLM | NB GLMM | tscount | CGFM ind  | CGFM shared |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
|  1 | 59.21 | -13.43 | -12.74 | 27.93 | 48.63 | -3.37 | 28.32 | -29.74 |
|  2 | 31.08 | -1.03 | -6.58 | 0.43 | 8.09 | 81.10 | 0.50 | -40.22 |
|  3 | 60.01 | -39.15 | -41.30 | -38.02 | -37.53 | -14.40 | -37.78 | -62.92 |
|  4 | 50.74 | -34.60 | -18.97 | -16.03 | -13.91 | 4.97 | -16.55 | -49.92 |
|  5 | 43.06 | -15.02 | -9.95 | 38.65 | 46.91 | 73.54 | 36.36 | -21.39 |
|  6 | 30.71 | -21.40 | -22.56 | 11.96 | 18.67 | 12.29 | 11.17 | -33.27 |
|  7 | 24.95 | -10.02 | -12.22 | -1.80 | -1.31 | 15.24 | -2.19 | -38.63 |
|  8 | 22.84 | -16.50 | -21.75 | -10.44 | -6.57 | -6.62 | -10.22 | -43.68 |
|  9 | 16.08 | 22.96 | 17.02 | 29.77 | 35.68 | 51.37 | 29.76 | -21.24 |
| 10 | 30.35 | -27.70 | -20.68 | -8.04 | -6.41 | 13.40 | -7.45 | -44.44 |
| 11 | 29.73 | -34.75 | -37.19 | -27.99 | -25.55 | -14.88 | -27.33 | -53.90 |
| 12 | 33.81 | -36.86 | -15.69 | -41.04 | -41.72 | -42.62 | -41.40 | -62.54 |

### 6.7 Marginal-moment headline

The INAD model class gives the closest overall reproduction of the empirical marginal moments. Unconstrained INAD is best for G1 means, G2 means, and G2 variances; constrained-alpha INAD is best for G1 variances. This supports reporting both INAD variants in the marginal-moment table rather than forcing one bolus variant to carry every criterion.

## 7. Cross-cutting findings

- Predictions are plug-in MLE conditional, not posterior predictive.
- INAD predictive samples are integer-valued by construction; predictive means are real-valued and are the right point forecasts for RMSE.
- Tau confidence intervals are profile-likelihood intervals. Tau SEs are derived from profile-CI width for internal consistency. A 3+ group unit test remains pending.
- The engineering pipeline completed the R = 1000 simulation, 1050 bolus all-pairs prediction, and full-data marginal-moment drivers without fit failures.

## 8. Implications for antedep

Already landed: `.simulate_inad_forward()` and public `predict.inad_fit()` v1 in `R/predict_inad.R`, with `type = c("mean", "sample")` and tests in `tests/testthat/test-predict-inad.R`.

Remaining: a shared cross-family `type = "distribution"` API after the gau/cat POCs; a proper `alpha_constraint = "constant"` path in `fit_inad()` with joint re-estimation; exact one-step PMFs as future methodological cleanup; and optional scoring helpers only if distribution-valued prediction becomes public.

Next sequence: start the gau and cat script-first prediction POCs, stop at smoke-test validation, then decide the shared prediction API before moving gau/cat prediction methods into the package.

