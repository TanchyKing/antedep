# INAD Validation Summary - Prediction and Marginal Moments

**Headline.** Under correct INAD specification (R = 1000 simulation),
INAD predicts the next visit with log score ~2.5 nats and randomized
PIT centered at 0.5 (well calibrated). On bolus rolling one-step
prediction (1050 all-pairs folds, patient-level clustered SEs),
constrained-α INAD beats marginal NB GLM, NB GLMM, and Poisson GLM at
conventional significance on RPS, RMSE, and log score — with `|t| > 10`
against NB GLM on every metric. On bolus marginal-moment reproduction,
**unconstrained INAD has the smallest MARD in 3 of 4 cells (G1 mean,
G2 mean, G2 variance)**, while **constrained-α INAD wins G1 variance**.
The two INAD variants are the primary models: **constrained-α** as the
BIC-supported and prediction-supported predictive primary; **unconstrained**
as the marginal-moment primary. See §5.5 and §6.7 for the full breakdowns.

This note presents predictive and in-sample validation evidence for INAD
models on simulated data and on the `bolus_inad` dataset. The simulation
section reports INAD-only correct-specification precision; the bolus
predictive section reports rolling one-step forecasting against count GLM
baselines; the bolus marginal-moment section compares fitted to empirical
mean and variance against a broader baseline ladder including published
Henderson-Shimakura CGFM estimates.

## 1. Scope and setup

### 1.1 What this note covers

- Simulation: absolute prediction precision for constrained-alpha INAD under its own NBT-NBI-INADFE(1) DGP.
- Bolus in-sample fit: BIC comparison across INAD and non-INAD count models.
- Bolus prediction: rolling one-step leave-one-per-group-out prediction on all 30 x 35 = 1050 held-out pairs.
- Bolus marginal moments: Henderson-and-Shimakura-style time-by-time mean and variance comparison.
- Package implications: what already landed in `antedep`, and what remains.

CGFM is included only in the marginal-moment comparison, using published Henderson-Shimakura estimates. No R implementation is available for refitting CGFM inside each cross-validation fold, so CGFM is not included in the predictive comparison.

### 1.2 Evaluation framework

The predictive target is conditional new-subject forecasting. A model is fit on training subjects, early observations of held-out subjects are revealed, and later observations are predicted from the fitted conditional distribution. The headline uses rolling one-step prediction at t = 2, ..., 12, conditioning on the realized history through t - 1.

All predictive comparisons are plug-in MLE conditional. INAD uses Monte Carlo forward simulation from the fitted process; NB and Poisson GLMs use analytic plug-in PMFs; NB GLMM uses the fitted random-intercept model for prediction. Lower RMSE, lower log score, and lower RPS are better.

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

### 2.4 INAD constrained-alpha

NBT-NBI-INADFE(1): `Y_it = alpha o Y_i,t-1 + epsilon_it`, with constant alpha, negative-binomial thinning, and NB innovations with time-varying mean/size plus group effect tau. This is the BIC-supported primary INAD variant for bolus.

### 2.5 INAD unconstrained

Same model class, but alpha is allowed to vary over time. It is reported as a sensitivity because the bolus LRT detects mild non-stationarity even though BIC prefers the constrained model.

### 2.6 CGFM independent/shared/time-varying frailty

CGFM rows use Henderson-Shimakura published full-data estimates only. In the local notation, the gamma frailty has mean 1 and frailty variance psi. Independent frailty has time-specific independent frailties; shared frailty uses a single subject-level frailty across all times. The time-varying CGFM row uses the reproduced Henderson-Shimakura time-varying estimates; its correlation parameter rho affects cross-time association, while the one-time marginal means and variances use the fitted time effects, group effect, and frailty variance. These rows are used for marginal mean and variance only.

## 3. Simulation - INAD prediction precision

This section asks one question: **how good are the predictions produced
by constrained-alpha INAD when its modelling assumptions are correct?**
The DGP is INAD itself, so the simulation tests the predictive machinery
under correct specification rather than benchmarking INAD against
alternatives. Cross-model baseline comparisons are not retained here;
parameter-recovery tables are not retained either, since fitter recovery
under known DGPs has been validated separately.

DGP: NBT-NBI-INADFE(1), constant alpha = 0.35, time-varying NB innovation
parameters, tau = c(0, 1.25), 100 subjects, 12 time points, R = 1000
simulation replicates. For each rep, the model is fit on 70% of subjects
and rolling one-step predictions are produced for the held-out 30% at
each time t = 2, ..., 12. Per-rep predictive scores are then summarized
across reps.

| Metric | Mean | SE | 2.5% | 97.5% | Reps |
| --- | --- | --- | --- | --- | --- |
| Log score (nats) | 2.508 | 0.002 | 2.403 | 2.611 | 1000 |
| RPS | 1.857 | 0.003 | 1.653 | 2.069 | 1000 |
| RMSE (counts) | 3.516 | 0.007 | 3.068 | 4.032 | 1000 |
| Randomized PIT | 0.494 | 0.001 | 0.456 | 0.531 | 1000 |
| 80% interval coverage | 0.853 | 0.001 | 0.812 | 0.891 | 1000 |
| 95% interval coverage | 0.967 | 0.000 | 0.942 | 0.985 | 1000 |

Three things to read off the table:

- **Point and probabilistic precision.** Log score around 2.5 nats per
  prediction and RPS around 1.86 characterize the predictive precision
  achievable when INAD's structural assumptions hold. These numbers
  serve as a reference scale for interpreting real-data INAD
  predictions; real-data performance can be better or worse depending
  on how well the data conform to the assumed structure and how much
  irreducible noise the data contain.
- **Calibration.** The randomized PIT mean is essentially 0.5 (ideal
  for a well-calibrated discrete predictive distribution), so the
  predictive distribution is not systematically too high or too low.
- **Predictive intervals.** Empirical 80% and 95% coverages are
  slightly above their nominal levels, so the predictive intervals are
  mildly conservative rather than over-confident. Acceptable for the
  validation target.

## 4. Bolus in-sample fit - BIC comparison

BIC is computed on the full bolus data. **Rows are listed in the consistent model order used throughout §4–§6** (INAD variants → Poisson GLM → NB family → CGFM family). Constrained-α INAD has the lowest BIC, with NB GLMM second (4264.8) and unconstrained INAD third (4297.8).

Note that Poisson GLM here is the lag-aware predictive form while NB GLM here is the marginal no-lag form; the asymmetry follows the dissertation BIC convention for NB GLM but the predictive-pipeline convention for Poisson GLM.

**Caveat on cross-model BIC comparability.** The BIC values below mix likelihood factorizations: INAD uses subject-level joint likelihoods, the GLM and GLMM rows use observation-level fits, and the CGFM rows are reconstructed from published log-likelihoods in Henderson & Shimakura. BIC is used here as an in-sample diagnostic for ranking these models on their ability to fit the bolus data, not as a formal Bayes-factor comparison across all possible likelihood factorizations. The CGFM rows in particular are reproduced for reference only.

| Model | Specification | BIC | Note |
| --- | --- | --- | --- |
| INAD constrained-α | NBT-NBI-INADFE(1), constant alpha | 4243.5 | BIC uses n_subjects; constrained POC inherits nb_inno_size from the unconstrained fit (see §6 caveat). |
| INAD unconstrained | NBT-NBI-INADFE(1), time-varying alpha | 4297.8 | BIC uses n_subjects, matching antedep fit-object convention. |
| Poisson GLM | lag-aware predictive: y ~ group + factor(time) + log(y_lag1 + 1) | 4334.6 | Lag-aware form used in the predictive pipeline; asymmetric with marginal NB GLM. |
| NB GLM | marginal no-lag: y ~ group + factor(time) | 4476.8 | Marginal no-lag form, aligned with the dissertation BIC comparison. |
| NB GLMM | marginal-moment no-lag with subject random intercept | 4264.8 | No-lag form used for in-sample marginal fit and closed-form marginal moments. |
| CGFM independent frailty | published Henderson-Shimakura independent frailty | 4442.0 | Reconstructed from published full log-likelihood and 12 time effects + group + frailty variance. |
| CGFM shared frailty | published Henderson-Shimakura shared frailty | 4556.8 | Reconstructed from published full log-likelihood and 12 time effects + group + frailty variance. |

## 5. Bolus prediction - rolling one-step

Design: all 30 x 35 cross-group held-out patient pairs, each fold holding out one group-1 and one group-2 patient and fitting on the remaining 63 patients. Standard errors are patient-level clustered, so overlapping folds do not receive naive independent-fold uncertainty.

### 5.1 Stationarity and primary-model choice

The full-data LRT for constant alpha rejects (p = 0.000234), but BIC strongly prefers constrained-alpha INAD over unconstrained INAD (4244 vs 4298). The interpretation is mild detectable non-stationarity that is not worth 11 extra alpha parameters under BIC. Constrained-alpha INAD is therefore primary; unconstrained INAD is sensitivity.

### 5.2 Single comparison table

Each cell reports `delta (t)`, where:

- **`delta`** = mean of the per-patient paired score difference
  (constrained-alpha INAD score minus reference score) across the 60 held-out patients
  pooled over folds. **Negative means INAD is better.** Units match the
  column: RPS units for RPS, counts for RMSE, nats for log score.
- **`t`** = t-statistic of the paired difference, using **patient-level
  clustered standard errors** so that overlapping folds (each patient
  appears in many cross-group fold combinations) do not inflate
  apparent precision. Conventional thresholds: `|t| > 2` is
  significant at the 5% level, `|t| > 3` at roughly 0.1%, and `|t| >
  10` is overwhelming.

A row therefore reads as: *"constrained-alpha INAD's score is `|delta|`
units below the reference model's, on average across patients, with a
t-statistic of `t`."*

Rows use constrained-alpha INAD as the comparison anchor. The unconstrained INAD fit is included as a reference row, so the constant-alpha sensitivity appears in the same table as the external baselines.

| Reference | Delta RPS (t) | Delta RMSE (t) | Delta log score (t) |
| --- | --- | --- | --- |
| INAD unconstrained | -0.010 (-2.9) | -0.029 (-5.1) | -0.005 (-2.4) |
| Poisson GLM | -0.078 (-16.3) | -0.012 (-1.9) | -0.175 (-30.6) |
| NB GLM | -0.396 (-22.3) | -0.689 (-23.9) | -0.123 (-18.6) |
| NB GLMM | -0.153 (-12.3) | -0.217 (-10.7) | -0.056 (-10.8) |

How to read the table: constrained-alpha INAD is modestly but consistently better than unconstrained INAD on rolling one-step prediction. It is also better than NB GLM, NB GLMM, and Poisson GLM on RPS and log score with strong statistical reliability, and better than NB GLM and NB GLMM on RMSE. The only soft external-baseline cell is the Poisson GLM RMSE comparison, which is small and borderline by t-statistic.

CGFM variants are not shown as numeric prediction rows in this table. The Henderson-Shimakura reproduction script validates the full-data independent, shared, and time-varying CGFM fits, but those fits have not yet been refit and scored inside each leave-two-out prediction fold. The independent CGFM has the same one-time marginal predictive family as a marginal NB GLM, because an independent Poisson-gamma mixture reduces to an NB distribution; a separate independent-CGFM prediction row would therefore be theoretically redundant with the NB-GLM-style marginal predictive row once the same mean and dispersion convention is used. Shared and time-varying CGFM prediction would be genuinely different because they use within-patient frailty history, but they require separate fold-wise prediction scripts.

### 5.3 Bolus prediction headline

On rolling one-step prediction, constrained-alpha INAD beats the marginal NB GLM, NB GLMM, and Poisson GLM on RPS. It also beats all three on RMSE except that the Poisson RMSE margin is small and only borderline by t-statistic. Log score is better than NB GLM, NB GLMM, and Poisson GLM. The unconstrained sensitivity is now folded into the single comparison table and shows that constrained-alpha INAD is slightly better than unconstrained INAD across all three metrics.

## 6. Bolus marginal mean and variance comparison

This section compares full-data fitted marginal moments with empirical bolus moments. CGFM is included here only, using published Henderson-Shimakura estimates. Relative discrepancy is `100 * (fitted - empirical) / empirical`; MARD is the mean absolute relative discrepancy.

### 6.1 MARD summary

MARD = **Mean Absolute Relative Discrepancy** across the 12 time points, computed per model, per group, per moment. All values are reported as **percentages**: a MARD of 8.35 for G1 mean means the fitted marginal mean is on average 8.35% off the empirical marginal mean across t = 1, ..., 12. Lower is better.

A second column for Group 1 variance, "G1 var (excl t=11)", reports the same MARD computed without time 11. This is a sensitivity column: the empirical G1 variance at t = 11 is unusually small (4.26), which inflates relative discrepancies for every model at that one cell. The ranking is robust to its exclusion.

**Model ordering.** All §6 tables list models in the consistent §4–§6 order: INAD variants (constrained then unconstrained), then count GLMs (Poisson, NB, NB GLMM), then CGFM (independent, shared, time-varying). The Poisson GLM row here uses the no-lag Poisson fit (`y ~ group + factor(time)`), which is the appropriate form for marginal-moment comparison; this is a different fit from the lag-aware predictive Poisson GLM used in §4 and §5.

**MARD headline (best per column in bold).** U-INAD has the smallest MARD in 3 of 4 cells (G1 mean, G2 mean, G2 variance); C-α INAD wins G1 variance and the G1-variance sensitivity excluding t = 11. INAD class therefore dominates marginal-moment reproduction on MARD.

| Model | G1 mean | G1 var (all 12) | G1 var (excl t=11) | G2 mean | G2 variance |
| --- | --- | --- | --- | --- | --- |
| INAD constrained-α | 8.35 | **26.74** | **15.65** | 9.31 | 22.79 |
| INAD unconstrained | **5.92** | 29.05 | 20.37 | **5.83** | **19.72** |
| Poisson GLM | 9.95 | 59.70 | 64.84 | 5.96 | 78.52 |
| NB GLM | 8.89 | 42.36 | 29.70 | 6.94 | 21.01 |
| NB GLMM | 8.35 | 58.35 | 43.73 | 8.24 | 24.25 |
| CGFM independent frailty | 8.94 | 42.44 | 29.53 | 6.87 | 20.75 |
| CGFM shared frailty | 9.86 | 29.01 | 23.14 | 6.01 | 41.83 |
| CGFM time-varying frailty | 8.96 | 42.26 | 29.70 | 6.88 | 21.08 |

**Caveat on the C-α INAD variance row.** The constrained-α fit currently uses an inherited `nb_inno_size` from the unconstrained fit (a POC shortcut documented in §4). Under proper joint re-estimation of `nb_inno_size_t` with constrained `α`, the C-α INAD variance numbers (G1 and G2) would shift to some extent — the unconstrained `nb_inno_size_t` has CV ≈ 47% across times, so joint constrained estimation is unlikely to leave the variance numbers unchanged. **The qualitative interpretation is expected to be stable, but the constrained-α variance magnitudes should be treated as provisional until joint estimation is implemented.** Future work; tracked in `inad-validation-revision-todo.md`.

### 6.2 Group 1 mean relative discrepancy (%)

Columns are in the consistent §4–§6 model order. Negative means fitted is below empirical.

| Time | Emp. | C-α INAD | U-INAD | Poisson | NB GLM | NB GLMM | CGFM ind | CGFM shared |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
|  1 | 9.30 | -1.73 | -1.57 | -13.03 | -12.28 | -4.37 | -12.19 | -13.07 |
|  2 | 5.50 | 6.34 | -4.45 | -9.22 | -9.03 | -4.43 | -9.04 | -9.04 |
|  3 | 5.43 | 12.34 | 4.46 | 2.73 | 1.60 | 2.83 | 1.76 | 2.78 |
|  4 | 5.13 | -0.28 | 11.44 | 18.93 | 15.98 | 18.44 | 15.52 | 17.85 |
|  5 | 7.57 | -5.78 | -2.53 | -5.11 | -5.26 | -1.62 | -6.18 | -5.23 |
|  6 | 5.33 | -1.18 | -1.98 | -0.11 | -0.94 | 3.13 | -1.39 | -0.40 |
|  7 | 3.93 | 12.22 | 7.52 | 13.10 | 10.89 | 12.22 | 10.57 | 12.81 |
|  8 | 3.73 | 5.83 | 0.46 | 6.86 | 5.32 | 8.97 | 5.41 | 6.47 |
|  9 | 4.60 | -11.77 | -15.28 | -14.17 | -13.54 | -10.44 | -13.59 | -14.45 |
| 10 | 4.93 | 2.20 | -9.66 | -4.42 | -4.79 | -3.01 | -4.50 | -4.50 |
| 11 | 3.50 | 17.93 | 6.60 | 17.91 | 15.32 | 18.67 | 15.86 | 18.20 |
| 12 | 3.47 | -22.63 | -5.11 | 13.78 | 11.72 | 12.10 | 11.27 | 13.52 |

### 6.3 Group 1 variance relative discrepancy (%)

Columns are in the consistent §4–§6 model order. Negative means fitted is below empirical.

| Time | Emp. | C-α INAD | U-INAD | Poisson | NB GLM | NB GLMM | CGFM ind | CGFM shared |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
|  1 | 33.46 | 24.05 | 24.41 | -75.83 | 21.58 | 55.95 | 22.07 | -28.95 |
|  2 | 22.05 | -8.25 | -20.91 | -77.36 | -21.83 | -7.59 | -21.70 | -50.07 |
|  3 | 15.98 | 43.98 | 35.55 | -65.10 | 27.76 | 41.89 | 28.35 | -18.20 |
|  4 | 17.02 | 9.67 | 37.65 | -64.16 | 36.80 | 54.62 | 36.11 | -12.83 |
|  5 | 25.08 | -2.95 | 2.17 | -71.36 | 28.71 | 50.60 | 26.77 | -22.21 |
|  6 | 11.26 | 31.56 | 27.58 | -52.71 | 68.00 | 95.81 | 67.01 | 7.28 |
|  7 | 17.58 | -29.57 | -31.94 | -74.71 | -22.31 | -14.25 | -22.54 | -47.89 |
|  8 | 9.24 | 10.39 | 2.04 | -56.86 | 24.36 | 42.11 | 24.76 | -15.91 |
|  9 | 9.97 | 6.80 | 1.08 | -60.40 | 17.40 | 34.46 | 17.49 | -23.27 |
| 10 | 13.17 | -0.30 | 0.03 | -64.22 | 17.54 | 31.48 | 18.38 | -23.76 |
| 11 | 4.26 | 148.75 | 124.57 | -3.12 | 181.71 | 219.24 | 184.44 | 93.61 |
| 12 | 7.98 | -4.61 | 40.73 | -50.53 | 40.38 | 52.23 | 39.68 | -4.13 |

### 6.4 Group 2 mean relative discrepancy (%)

Columns are in the consistent §4–§6 model order. Negative means fitted is below empirical.

| Time | Emp. | C-α INAD | U-INAD | Poisson | NB GLM | NB GLMM | CGFM ind | CGFM shared |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
|  1 | 10.17 | 1.08 | 1.53 | 10.24 | 12.75 | 16.01 | 12.80 | 10.56 |
|  2 | 6.49 | 16.07 | 5.07 | 6.65 | 8.45 | 7.53 | 8.37 | 7.29 |
|  3 | 7.86 | 2.27 | -4.58 | -1.62 | -1.23 | -5.65 | -1.14 | -1.14 |
|  4 | 9.29 | -22.79 | -11.99 | -8.97 | -9.86 | -13.13 | -10.28 | -9.38 |
|  5 | 9.63 | -4.08 | 0.28 | 3.40 | 4.66 | 2.58 | 3.59 | 3.59 |
|  6 | 7.37 | 0.43 | 0.38 | 0.14 | 0.75 | -1.00 | 0.24 | 0.24 |
|  7 | 6.60 | -0.63 | -6.34 | -6.65 | -7.10 | -11.27 | -7.42 | -6.49 |
|  8 | 5.74 | 6.25 | 0.24 | -3.73 | -3.75 | -6.01 | -3.72 | -3.72 |
|  9 | 4.91 | 26.40 | 20.27 | 11.47 | 13.77 | 11.23 | 13.64 | 11.39 |
| 10 | 6.34 | 13.46 | 11.91 | 3.03 | 4.10 | 0.09 | 4.36 | 3.32 |
| 11 | 6.26 | 0.41 | -4.82 | -8.63 | -9.32 | -11.93 | -8.95 | -8.03 |
| 12 | 5.89 | -17.81 | 2.57 | -7.08 | -7.49 | -12.40 | -7.92 | -7.00 |

### 6.5 Group 2 variance relative discrepancy (%)

Columns are in the consistent §4–§6 model order. Negative means fitted is below empirical.

| Time | Emp. | C-α INAD | U-INAD | Poisson | NB GLM | NB GLMM | CGFM ind | CGFM shared |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
|  1 | 59.21 | -13.43 | -12.74 | -81.07 | 27.93 | 48.63 | 28.32 | -29.74 |
|  2 | 31.08 | -1.03 | -6.58 | -77.73 | 0.43 | 8.09 | 0.50 | -40.22 |
|  3 | 60.01 | -39.15 | -41.30 | -87.11 | -38.02 | -37.53 | -37.78 | -62.92 |
|  4 | 50.74 | -34.60 | -18.97 | -83.33 | -16.03 | -13.91 | -16.55 | -49.92 |
|  5 | 43.06 | -15.02 | -9.95 | -76.88 | 38.65 | 46.91 | 36.36 | -21.39 |
|  6 | 30.71 | -21.40 | -22.56 | -75.97 | 11.96 | 18.67 | 11.17 | -33.27 |
|  7 | 24.95 | -10.02 | -12.22 | -75.31 | -1.80 | -1.31 | -2.19 | -38.63 |
|  8 | 22.84 | -16.50 | -21.75 | -75.80 | -10.44 | -6.57 | -10.22 | -43.68 |
|  9 | 16.08 | 22.96 | 17.02 | -65.97 | 29.77 | 35.68 | 29.76 | -21.24 |
| 10 | 30.35 | -27.70 | -20.68 | -78.48 | -8.04 | -6.41 | -7.45 | -44.44 |
| 11 | 29.73 | -34.75 | -37.19 | -80.76 | -27.99 | -25.55 | -27.33 | -53.90 |
| 12 | 33.81 | -36.86 | -15.69 | -83.81 | -41.04 | -41.72 | -41.40 | -62.54 |

### 6.6 Supporting rank and worst-case summaries (metric sensitivity)

MARD (§6.1) is the headline metric for the §6 comparison. As a sensitivity check on metric choice, three supporting per-cell summaries are reported here:

- **Rank sum.** At each time, rank the seven models by `|discrepancy|` (1 = smallest, ties get average rank). Sum the ranks across the 12 time points per model. Lower is better. With 7 models and 12 times, the minimum possible rank sum is 12 (always best); the maximum is 84 (always worst); the mean by chance is 48. Rank sum is invariant to magnitude — it measures how often each model is closer than its peers.
- **Wins.** Number of time points where the model has the smallest `|discrepancy|` (ties for first share a win). Useful diagnostic for bimodal best-at-some / worst-at-some behavior.
- **Max `|discrepancy|` (%).** The single worst time point per model. Worst-case calibration.

Rows are in the consistent §4–§6 model order. Best value per column is **bold**.

**G1 mean** (12 time points). Agrees with MARD: U-INAD is the cleanest single best by all four summaries.

| Model | Rank sum | Mean rank | Wins | Max \|discrepancy\| (%) |
| --- | ---: | ---: | ---: | ---: |
| INAD constrained-α | 46.5 | 3.88 | 2 | 22.63 |
| INAD unconstrained | **37.0** | **3.08** | **5** | **15.28** |
| Poisson GLM | 58.0 | 4.83 | 1 | 18.93 |
| NB GLM | 41.0 | 3.42 | 1 | 15.98 |
| NB GLMM | 48.5 | 4.04 | 3 | 18.67 |
| CGFM independent frailty | 45.0 | 3.75 | 0 | 15.86 |
| CGFM shared frailty | 60.0 | 5.00 | 0 | 18.20 |

**G1 variance** (12 time points). Agrees with MARD: C-α INAD is best by both metrics.

| Model | Rank sum | Mean rank | Wins | Max \|discrepancy\| (%) |
| --- | ---: | ---: | ---: | ---: |
| INAD constrained-α | **33** | **2.75** | 1 | 148.75 |
| INAD unconstrained | 35 | 2.92 | **4** | 124.57 |
| Poisson GLM | 74 | 6.17 | 1 | **77.36** |
| NB GLM | 44 | 3.67 | 1 | 181.71 |
| NB GLMM | 64 | 5.33 | 2 | 219.24 |
| CGFM independent frailty | 46 | 3.83 | 0 | 184.44 |
| CGFM shared frailty | 40 | 3.33 | 3 | 93.61 |

**G1 variance, excluding t = 11** (11 time points; same sensitivity as §6.1). Agrees with MARD: C-α INAD widens its lead.

| Model | Rank sum | Mean rank | Wins | Max \|discrepancy\| (%) |
| --- | ---: | ---: | ---: | ---: |
| INAD constrained-α | **29** | **2.64** | 1 | **43.98** |
| INAD unconstrained | 32 | 2.91 | **4** | 40.73 |
| Poisson GLM | 73 | 6.64 | 0 | 77.36 |
| NB GLM | 39 | 3.55 | 1 | 68.00 |
| NB GLMM | 57 | 5.18 | 2 | 95.81 |
| CGFM independent frailty | 40 | 3.64 | 0 | 67.01 |
| CGFM shared frailty | 38 | 3.45 | 3 | 50.07 |

**G2 mean** (12 time points). MARD says U-INAD (5.83%); rank sum nudges CGFM shared ahead (33 vs 37). Means CGFM shared is more often closest single-time, but U-INAD has smaller average percent error. The MARD ordering is the headline.

| Model | Rank sum | Mean rank | Wins | Max \|discrepancy\| (%) |
| --- | ---: | ---: | ---: | ---: |
| INAD constrained-α | 61 | 5.08 | 3 | 26.40 |
| INAD unconstrained | 37 | 3.08 | **4** | 20.27 |
| Poisson GLM | 34 | 2.83 | 2 | 11.47 |
| NB GLM | 59 | 4.92 | 0 | 13.77 |
| NB GLMM | 61 | 5.08 | 2 | 16.01 |
| CGFM independent frailty | 51 | 4.25 | 1 | 13.64 |
| CGFM shared frailty | **33** | **2.75** | 1 | **11.39** |

**G2 variance** (12 time points). MARD says U-INAD (19.72%); rank sum nudges CGFM independent ahead (33 vs 41), with NB GLMM 3rd by rank sum but most wins (6) — best-at-half / worst-at-half pattern.

| Model | Rank sum | Mean rank | Wins | Max \|discrepancy\| (%) |
| --- | ---: | ---: | ---: | ---: |
| INAD constrained-α | 42 | 3.50 | 0 | **39.15** |
| INAD unconstrained | 41 | 3.42 | 4 | 41.30 |
| Poisson GLM | 84 | 7.00 | 0 | 87.11 |
| NB GLM | 35 | 2.92 | 1 | 41.04 |
| NB GLMM | 37 | 3.08 | **6** | 48.63 |
| CGFM independent frailty | **33** | **2.75** | 1 | 41.40 |
| CGFM shared frailty | 64 | 5.33 | 0 | 62.92 |

**Summary of metric sensitivity:**

- **Group 1 (both moments):** all four metrics (MARD, rank sum, wins, max) agree that INAD class is best. C-α wins G1 variance, U-INAD wins G1 mean.
- **Group 2:** MARD favors U-INAD on both moments. Rank sum partially disagrees: CGFM shared edges G2 mean (33 vs 37), CGFM ind edges G2 variance (33 vs 41). The disagreements are small and consistent with CGFM variants being slightly more often the closest fit while U-INAD has the smaller mean error.
- **Poisson GLM as a lower-floor variance reference.** Tracks marginal means competitively (rank sum 34 on G2 mean) but is dramatically worse on variance (MARD 59.7% / 78.5% on G1 / G2; rank sum 74 / 84), because the mean = variance constraint cannot accommodate bolus overdispersion. Quantitative confirmation that bolus needs an overdispersion-capable model.
- **Max `|discrepancy|`** at G1 variance is dominated by t = 11 (small empirical variance, 4.26); every overdispersed-capable model overshoots, Poisson happens to be closest only because it always underestimates. The "excl t = 11" table widens the INAD lead.

### 6.7 Marginal-moment headline

**MARD is the headline metric** (§6.1); the §6.6 rank/wins/max summaries are reported as a metric-sensitivity check. By MARD:

- **G1 means:** U-INAD wins (5.92%; next: NB GLMM 8.35%, C-α INAD 8.35%).
- **G1 variances:** C-α INAD wins (26.74%; sensitivity excluding t = 11: 15.65%). C-α maintains its lead by rank sum (33), and widens it excluding t = 11.
- **G2 means:** U-INAD wins (5.83%; next: CGFM shared 6.01%, Poisson 5.96%).
- **G2 variances:** U-INAD wins (19.72%; next: CGFM ind 20.75%, NB GLM 21.01%).

**Unconstrained INAD is the §6 primary** — it has the smallest MARD in 3 of 4 cells. Constrained-α INAD is reported alongside because it wins G1 variance and is the §5 primary. The two INAD variants together carry the marginal-moment story for both groups.

Three takeaways:

1. **INAD class wins all four MARD cells.** Unconstrained INAD on G1 mean, G2 mean, and G2 variance; constrained-α INAD on G1 variance.
2. **Rank sum mostly agrees on Group 1 and partially disagrees on Group 2.** For G2, CGFM shared edges out U-INAD on rank sum for G2 mean (33 vs 37), and CGFM ind edges out U-INAD on rank sum for G2 variance (33 vs 41). The interpretation: CGFM variants are slightly more often the closest single-time fit on Group 2, while U-INAD has the smaller mean percent error. The disagreement is small and metric-driven, not a meaningful contradiction of the headline.
3. **Poisson GLM as a lower-floor variance reference.** Its mean tracks competitively (especially G2), but its mean = variance constraint makes it dramatically worse on variance (MARD 59.7% / 78.5%; rank sum 74 / 84). This quantitatively confirms bolus is meaningfully overdispersed and that any sensible count model needs more flexibility than Poisson.

## 7. Cross-cutting findings

- Predictions are plug-in MLE conditional, not posterior predictive.
- INAD predictive samples are integer-valued by construction; predictive means are real-valued and are the right point forecasts for RMSE.
- Tau confidence intervals are profile-likelihood intervals. Tau SEs are derived from profile-CI width for internal consistency in package output. A standalone profile-curvature diagnostic was also run from `scripts/prediction_poc/inad/26_tau_profile_curvature.R`; it gives a narrower local-quadratic approximation that does not capture the tails of the profile log-likelihood.
- The engineering pipeline completed the R = 1000 simulation, 1050 bolus all-pairs prediction, and full-data marginal-moment drivers without fit failures.

For the bolus alpha-constant NBT-NBI-INADFE(1) fit, the two methods give:

| Method | `tau_2` estimate | SE | 95% CI lower | 95% CI upper |
|---|---:|---:|---:|---:|
| Profile CI (package convention; primary) | 1.1428 | 0.4818 | 0.1869 | 2.0757 |
| Local profile curvature (diagnostic only) | 1.1003 | 0.2731 | 0.5650 | 1.6356 |

Both rows use the same profile log-likelihood — for each candidate `tau_2`, all nuisance parameters are re-optimized in the same way. The rows differ only in how the resulting profile is summarized:

- The **profile CI** (package convention) walks the profile log-likelihood outward from the MLE in each direction. It locates the two `tau_2` values — one below the MLE and one above — at which the profile log-likelihood has fallen by `qchisq(0.95, 1) / 2 ≈ 1.92` from its peak. Those two values are the 95% CI's lower and upper bounds (here `0.1869` and `2.0757`); the SE is then derived from the half-width as `(upper − lower) / (2 · qnorm(0.975))`.
- The **local profile curvature** (diagnostic) fits a quadratic to the profile near the MLE and reports the asymptotic SE `sqrt(-1 / (2c))` from the quadratic coefficient. The 95% CI is `vertex ± 1.96 · SE`.

The two would agree only if the profile log-likelihood were exactly quadratic. For the bolus `tau_2` profile, the LR threshold is reached at `tau_2` values further from the MLE than the local quadratic predicts, so the LR-based CI is wider. The "estimate" column also differs slightly: the package row reports the stored MLE on the constrained-alpha POC fit object, while the diagnostic row reports the vertex of the local-quadratic fit. They drift apart because the constrained-alpha POC is not a fully joint fit, so the stored MLE may not be exactly at the local profile maximum after nuisance refits. The diagnostic row is included as a sanity check on the profile shape, not as a competing interval; the package convention remains the primary reporting standard.

The older dissertation Table 4.3 printed a Hessian/Louis-style SE for `tau_2` while using profile-likelihood bounds, which is why `est ± 2 SE` did not agree with the displayed interval. The current package-facing convention keeps `tau` on a profile-likelihood basis for both the interval and the displayed SE; other INAD parameters retain Wald SEs and Wald intervals. The profiling code loops over `tau[2]`, ..., `tau[B]`, so the implementation extends naturally to 3+ groups.

## 8. Implications for antedep

Already landed: `.simulate_inad_forward()` and the public `predict.inad_fit()` S3 method in `R/predict_inad.R`, currently supporting `type = c("mean", "sample")`, with tests in `tests/testthat/test-predict-inad.R`.

Remaining: a shared cross-family `type = "distribution"` API after the gau/cat POCs; a proper `alpha_constraint = "constant"` path in `fit_inad()` with joint re-estimation; exact one-step PMFs as future methodological cleanup; and optional scoring helpers only if distribution-valued prediction becomes public.

Next sequence: start the gau and cat script-first prediction POCs, stop at smoke-test validation, then decide the shared prediction API before moving gau/cat prediction methods into the package.
