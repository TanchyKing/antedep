# INAD Validation Summary — Prediction and Marginal Moments

This note consolidates two streams of empirical validation for INAD models:
out-of-sample predictive performance against alternative count time-series
models, and in-sample marginal mean and variance reproduction on the
`bolus_inad` dataset. It covers the controlled simulation study (R = 1000),
the bolus all-pairs predictive comparison (1050 folds with patient-level
clustering), and the bolus marginal moment comparison reported in a
Henderson-and-Shimakura-style time-by-time format. The full design is in
`prediction-plan.md`; scripts and output sit under `scripts/prediction_poc/`.

## 1. Scope and Setup

### 1.1 What this note covers

- Simulation predictive comparison (R = 1000) — §2.
- Bolus predictive comparison (1050 all-pairs folds) — §3.
- Bolus marginal mean and variance comparison — §4.
- Cross-cutting findings on plug-in prediction, integer-valued outputs,
  log-score sensitivity, and tau interval reporting — §5.
- Implications for the `antedep` package — §6.
- Bridge to Gaussian and categorical AD prediction work — §7.

### 1.2 Evaluation framework

**Predictive design.** Conditional new-subject forecasting throughout. Fit
on training subjects using all 12 time points; reveal a held-out subject's
history through `T_obs` and predict the suffix from the model's conditional
distribution. Temporal holdout was rejected because the time-indexed
parameters do not extend past `T_train` without additional assumptions.

**Two prediction tasks.**

- *Rolling one-step* at `t = 2, …, 12`, conditioning on the realized
  history through `t-1`. The "given everything observed so far, predict
  the next visit" task.
- *Recursive multi-step* from split time `T_split = 8` for horizons
  `h ∈ {1, 2, 3, 4}`, propagated through the predictive distribution.

**Headline metrics.** RPS and RMSE. Log score reported alongside but
treated as secondary because it depends on the empirical PMF at a single
point and is sensitive to the smoothing constant on multi-step horizons
(see §5.3).

**INAD predictive object.** Monte Carlo estimate of the **plug-in
MLE-conditional predictive distribution** at `n_sims = 1000` paths,
marginalized over latent thinning and innovation noise. Not a posterior
predictive — no integration over parameter uncertainty. NB GLM and
Poisson GLM baselines use **analytic plug-in PMFs** at the MLE; NB GLMM
and tscount use analytic marginalization or Monte Carlo recursion as
appropriate. All models share the same plug-in convention; only the
predictive-PMF construction differs across families, so the comparison
is internally consistent.

**Marginal-moment comparison.** For each candidate model fit on the full
bolus data, the model-implied marginal `E[Y_t | group]` and
`Var[Y_t | group]` are computed per time per group and compared against the
empirical bolus moments in a Henderson-and-Shimakura Table-3-style format.
Closed-form recursion for INAD and NB GLM; analytic random-intercept
marginalization for NB GLMM; Monte Carlo panel recursion seeded from
empirical time-1 counts for tscount; published estimates with closed-form
marginal-moment formulas for the CGFM rows.

### 1.3 Models compared

| Model | Role |
|---|---|
| INAD constrained-α | Primary INAD fit — α held constant across time |
| INAD unconstrained | Sensitivity INAD fit — `α_t` estimated freely |
| NB GLM | Practitioner NB-family baseline |
| NB GLMM | NB GLM with subject random intercept (fairness check) |
| tscount | Generalized-linear autoregressive count time-series comparator |
| Poisson GLM | Lighter secondary comparator |
| CGFM independent frailty | Henderson–Shimakura model, marginal-moment row |
| CGFM shared frailty | Henderson–Shimakura model, marginal-moment row |

The constrained-α fit adapts the constrained machinery from
`test_stationarity_inad(..., constrain = "alpha")` into an
`inad_fit`-shaped object. `nb_inno_size` is inherited from the
unconstrained fit; a fully joint re-estimation under the constraint is
documented as a remaining package item (§6.2).

## 2. Simulation Results (R = 1000)

**DGP:** NBT-NBI-INADFE(1) with constant `α = 0.35`, time-varying NB
innovation parameters, two-group fixed effects (`τ = c(0, 1.25)`).
Subjects: 100 (50 per group). Time points: 12. Train/test: 70/30 stratified
by group. Output: `scripts/prediction_poc/output/inad/final_r1000/` (NB GLM
comparison, R = 1000) and `scripts/prediction_poc/output/inad/extra_baselines_r200/`
(extra baseline ladder, R = 200).

### 2.1 Parameter recovery (R = 1000)

| Metric | Constrained-α | Unconstrained |
|---|---|---|
| α MAE vs 0.35 | 0.026 | 0.089 |
| θ RMSE | 0.419 | 0.576 |
| nb_inno_size RMSE | 1.598 | 1.598 (shared) |
| \|τ₂ bias\| | 0.197 | 0.197 |

Constrained-α substantially improves both α and θ recovery on the
simulated DGP. nb_inno_size is shared between the two INAD fits by
construction.

### 2.2 Constrained-α INAD vs NB GLM (R = 1000)

All deltas reported as `INAD − NB GLM`; negative is INAD better. Every
cell on RPS and RMSE is statistically reliable at `|t| ≥ 2` on across-rep
paired differences.

| Task | h | ΔRPS (t) | ΔRMSE (t) | Δlog-S (t) |
|---|---|---|---|---|
| rolling | 1 | −0.0155 (−27.5) | −0.0247 (−24.8) | −0.0031 (−6.3) |
| recursive | 1 | −0.0144 (−7.8) | −0.0215 (−6.8) | −0.0030 (−1.9) |
| recursive | 2 | −0.0070 (−4.9) | −0.0142 (−6.3) | +0.0053 (+4.0) |
| recursive | 3 | −0.0067 (−4.2) | −0.0113 (−5.1) | +0.0034 (+2.3) |
| recursive | 4 | −0.0070 (−4.6) | −0.0094 (−4.5) | ≈ 0 |

Constrained-α INAD beats NB GLM at every cell on RPS and RMSE. Log score
flips sign at recursive `h = 2, 3`, consistent with empirical-PMF log
score being sensitive to tail smoothing at multi-step horizons (see §5.3).

### 2.3 vs NB GLMM, Poisson GLM, tscount (R = 200, same DGP)

Constrained-α INAD beats all three baselines on RPS and RMSE at every
(task, horizon) cell, with the single soft cell being recursive h=3 vs
NB GLMM at `t ≈ −1.5 to −1.9` on RPS and RMSE — directionally INAD-better
but at the conventional significance edge.

### 2.4 Simulation headline

Under the correctly specified NBT-NBI-INADFE(1) DGP, constrained-α INAD
gives uniform RPS and RMSE gains over the evaluated baselines across both
prediction tasks and all horizons. Evidence level by baseline:

- **NB GLM comparison confirmed at R = 1000.**
- **Extra baseline ladder (NB GLMM, Poisson GLM, tscount) checked at R = 200.**

This establishes that the INAD structural assumptions translate into
predictive advantage when they hold.

## 3. Bolus Predictive Comparison (1050 all-pairs folds)

**Data:** `bolus_inad`, 65 subjects (30 group 1, 35 group 2), 12 time
points. **Design:** every cross-group pair held out — 30 × 35 = 1050
folds, each holding out one patient from group 1 and one from group 2,
fit on remaining 63. SEs computed by **patient-level clustering** to
respect the dependence between overlapping folds. **Output:**
`scripts/prediction_poc/output/inad/bolus_leave_one_per_group_all_pairs_marginal_nb_glm/`.

### 3.1 Stationarity check on full data

- LRT for constant α: `p = 0.000234` — rejects constancy.
- BIC: constrained 4244, unconstrained 4298 (Δ = 54) — prefers constrained.

LRT and BIC test different things. LRT detects any time-variation in α;
BIC penalizes the 11 extra parameters. Real bolus α has detectable but
mild time-variation that does not survive parsimony adjustment.
**Constrained-α is the BIC-supported primary; unconstrained is reported
as sensitivity.**

### 3.2 Constrained-α INAD vs NB GLM (patient-clustered)

| Task | h | ΔRPS (t) | ΔRMSE (t) | Δlog-S (t) |
|---|---|---|---|---|
| rolling | 1 | −0.396 (−22.3) | −0.689 (−23.9) | −0.123 (−18.6) |
| recursive | 1 | −0.222 (−8.2) | −0.361 (−8.3) | −0.103 (−8.1) |
| recursive | 2 | −0.086 (−5.7) | −0.165 (−6.8) | −0.011 (−1.4) |
| recursive | 3 | −0.073 (−8.5) | −0.113 (−8.2) | −0.014 (−2.8) |
| recursive | 4 | −0.007 (−0.6) | −0.081 (−4.5) | −0.019 (−2.1) |

Constrained-α INAD beats NB GLM on RPS and RMSE at every task/horizon
cell, with rolling one-step the largest and most reliable gap. The single
non-reliable cell on RPS is recursive h=4 (`t = −0.6`). Log score is
INAD-better at every cell, with recursive h=2 the only non-significant
log-score cell.

### 3.3 Unconstrained INAD vs NB GLM (patient-clustered)

| Task | h | ΔRPS (t) | ΔRMSE (t) | Δlog-S (t) |
|---|---|---|---|---|
| rolling | 1 | −0.386 (−23.0) | −0.659 (−24.6) | −0.118 (−18.9) |
| recursive | 1 | −0.199 (−8.3) | −0.342 (−8.8) | −0.082 (−7.0) |
| recursive | 2 | −0.091 (−4.6) | −0.144 (−4.5) | −0.008 (−0.9) |
| recursive | 3 | −0.095 (−12.8) | −0.167 (−13.3) | −0.026 (−6.1) |
| recursive | 4 | −0.082 (−8.4) | −0.091 (−5.8) | −0.054 (−7.7) |

Unconstrained INAD also beats NB GLM at every cell on RPS and RMSE, with
the cleanest gaps on the rolling task and on the longer recursive
horizons (h=3, h=4). The two INAD variants are close in headline metrics;
neither is strictly dominated.

### 3.4 Constrained-α INAD vs broader baseline ladder — rolling one-step

| Reference | ΔRPS (t) | ΔRMSE (t) | Δlog-S (t) |
|---|---|---|---|
| NB GLM | −0.396 (−22.3) | −0.689 (−23.9) | −0.123 (−18.6) |
| NB GLMM | −0.153 (−12.3) | −0.217 (−10.7) | −0.056 (−10.8) |
| Poisson GLM | −0.078 (−16.3) | −0.012 (−1.9) | −0.175 (−30.6) |
| tscount | −0.041 (−8.6) | −0.105 (−14.0) | +0.002 (+0.6) |

Constrained-α INAD beats NB GLM, NB GLMM, and Poisson GLM at conventional
significance on rolling one-step across RPS, RMSE, and log score (the
sole exception being Poisson RMSE at `t = −1.9`, just at the edge). Beats
tscount on RPS and RMSE; log score against tscount is essentially tied.

### 3.5 Smoothing sensitivity for the recursive log-score behaviour

For a sanity check on the recursive log-score gap, ε was varied on the
35-fold balanced rotation against NB GLM (constrained-α INAD,
recursive h=2):

| ε | Δlog-S vs NB GLM (t) |
|---|---|
| 1/10001 | +0.068 (+2.04) |
| 1/1001 (baseline) | +0.047 (+1.68) |
| 1/101 | −0.0005 (−0.02) |

The recursive `h = 2` log-score gap on the 35-fold pre-rotation was
**epsilon-sensitive**: at a more generous smoothing constant it
disappeared. This was the diagnostic that led to the smoothing-sensitivity
documentation in §5.3. The 1050-fold all-pairs result above (§3.2) shows
the all-pairs log-score gap is itself smaller and INAD-favourable at
recursive h=2, consistent with the smoothing diagnosis.

### 3.6 BIC-selected INAD on bolus

Per-fold training BIC selection between constrained and unconstrained gives
essentially the same predictive result as the constrained-only row on
rolling one-step. Training-BIC selection rates across the 1050 folds:

| Selected fit | Folds | Proportion |
|---|---:|---:|
| Constrained-α | 1034 | 98.48% |
| Unconstrained | 16 | 1.52% |

This pre-empts the "you picked the model after seeing test results"
reviewer concern. Provenance note: the BIC-selected per-fold files live
under the earlier all-pairs output directory
(`bolus_leave_one_per_group_all_pairs/`), while the new headline bolus
results in §3.2–§3.4 are from the marginal-NB-GLM all-pairs output
(`bolus_leave_one_per_group_all_pairs_marginal_nb_glm/`).

### 3.7 Bolus prediction headline

- **Rolling one-step:** constrained-α INAD beats NB GLM cleanly on all
  three metrics (RPS `t = −22.3`, RMSE `t = −23.9`, log score `t = −18.6`),
  and beats the wider baseline ladder consistently.
- **Recursive multi-step:** constrained-α INAD beats NB GLM on RPS and
  RMSE at every horizon (with one non-significant cell at h=4 on RPS),
  and on log score at h=1, h=3, h=4.
- **Constrained vs unconstrained on bolus:** both are INAD-favourable
  against NB GLM at every cell; differences between the two are small
  relative to their gap from the GLM baselines.

## 4. Bolus Marginal Mean and Variance Comparison

### 4.1 Method

For each model fit on the full bolus data, the model-implied marginal
mean and variance are computed per group per time and compared against
empirical sample mean and variance, in a Henderson-and-Shimakura
Table-3-style format. INAD and NB GLM use closed-form recursion; NB GLMM
uses analytic random-intercept marginalization (`λ_mean = μ_cond_zero ·
exp(σ²_b / 2)`, with the corresponding second-moment expression for the
variance); tscount uses Monte Carlo panel recursion seeded from empirical
time-1 counts; the CGFM rows use published parameter estimates with
their closed-form marginal-moment formulas. All values are evaluated on
the empirical bolus group/time grid. Output:
`scripts/prediction_poc/output/inad/marginal_moments/`.

Relative discrepancy is reported as `100 × (fitted − empirical) /
empirical`. Negative means under-fit; positive means over-fit.

### 4.2 Empirical bolus moments (Henderson Table-3 style)

| Time | G1 mean | G1 variance | G2 mean | G2 variance |
|---:|---:|---:|---:|---:|
| 1 | 9.300 | 33.459 | 10.171 | 59.205 |
| 2 | 5.500 | 22.052 | 6.486 | 31.081 |
| 3 | 5.433 | 15.978 | 7.857 | 60.008 |
| 4 | 5.133 | 17.016 | 9.286 | 50.739 |
| 5 | 7.567 | 25.082 | 9.629 | 43.064 |
| 6 | 5.333 | 11.264 | 7.371 | 30.711 |
| 7 | 3.933 | 17.582 | 6.600 | 24.953 |
| 8 | 3.733 | 9.237 | 5.743 | 22.844 |
| 9 | 4.600 | 9.972 | 4.914 | 16.081 |
| 10 | 4.933 | 13.168 | 6.343 | 30.350 |
| 11 | 3.500 | 4.259 | 6.257 | 29.726 |
| 12 | 3.467 | 7.982 | 5.886 | 33.810 |

### 4.3 Group 1 means — fitted by model

| Time | Empirical | C-α INAD | U-INAD | NB GLM | NB GLMM | tscount | CGFM ind. | CGFM shared |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 9.30 | 9.14 | 9.15 | 8.16 | 8.89 | 9.39 | 8.17 | 8.08 |
| 2 | 5.50 | 5.85 | 5.26 | 5.00 | 5.26 | 7.55 | 5.00 | 5.00 |
| 3 | 5.43 | 6.10 | 5.68 | 5.52 | 5.59 | 6.54 | 5.53 | 5.58 |
| 4 | 5.13 | 5.12 | 5.72 | 5.95 | 6.08 | 6.54 | 5.93 | 6.05 |
| 5 | 7.57 | 7.13 | 7.38 | 7.17 | 7.44 | 7.39 | 7.10 | 7.17 |
| 6 | 5.33 | 5.27 | 5.23 | 5.28 | 5.50 | 5.36 | 5.26 | 5.31 |
| 7 | 3.93 | 4.41 | 4.23 | 4.36 | 4.41 | 4.32 | 4.35 | 4.44 |
| 8 | 3.73 | 3.95 | 3.75 | 3.93 | 4.07 | 3.92 | 3.94 | 3.97 |
| 9 | 4.60 | 4.06 | 3.90 | 3.98 | 4.12 | 3.84 | 3.97 | 3.94 |
| 10 | 4.93 | 5.04 | 4.46 | 4.70 | 4.78 | 4.63 | 4.71 | 4.71 |
| 11 | 3.50 | 4.13 | 3.73 | 4.04 | 4.15 | 4.03 | 4.06 | 4.14 |
| 12 | 3.47 | 2.68 | 3.29 | 3.87 | 3.89 | 3.95 | 3.86 | 3.94 |

### 4.4 Group 1 means — relative discrepancy (%)

| Time | C-α INAD | U-INAD | NB GLM | NB GLMM | tscount | CGFM ind. | CGFM shared |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | −1.73 | −1.57 | −12.28 | −4.37 | 0.95 | −12.19 | −13.07 |
| 2 | 6.34 | −4.45 | −9.03 | −4.43 | 37.18 | −9.04 | −9.04 |
| 3 | 12.34 | 4.46 | 1.60 | 2.83 | 20.31 | 1.76 | 2.78 |
| 4 | −0.28 | 11.44 | 15.98 | 18.44 | 27.46 | 15.52 | 17.85 |
| 5 | −5.78 | −2.53 | −5.26 | −1.62 | −2.33 | −6.18 | −5.23 |
| 6 | −1.18 | −1.98 | −0.94 | 3.13 | 0.58 | −1.39 | −0.40 |
| 7 | 12.22 | 7.52 | 10.89 | 12.22 | 9.70 | 10.57 | 12.81 |
| 8 | 5.83 | 0.46 | 5.32 | 8.97 | 5.11 | 5.41 | 6.47 |
| 9 | −11.77 | −15.28 | −13.54 | −10.44 | −16.46 | −13.59 | −14.45 |
| 10 | 2.20 | −9.66 | −4.79 | −3.01 | −6.15 | −4.50 | −4.50 |
| 11 | 17.93 | 6.60 | 15.32 | 18.67 | 15.20 | 15.86 | 18.20 |
| 12 | −22.63 | −5.11 | 11.72 | 12.10 | 13.91 | 11.27 | 13.52 |
| **MARD** | **8.35** | **5.92** | 8.89 | 8.35 | 12.95 | 8.94 | 9.86 |

Unconstrained INAD gives the lowest mean absolute relative discrepancy
(MARD = 5.92%) for the Group 1 mean trajectory.

### 4.5 Group 1 variances — fitted by model

| Time | Empirical | C-α INAD | U-INAD | NB GLM | NB GLMM | tscount | CGFM ind. | CGFM shared |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 33.46 | 41.51 | 41.63 | 40.68 | 52.18 | 32.65 | 40.84 | 23.77 |
| 2 | 22.05 | 20.23 | 17.44 | 17.24 | 20.38 | 33.20 | 17.27 | 11.01 |
| 3 | 15.98 | 23.00 | 21.66 | 20.41 | 22.67 | 30.53 | 20.51 | 13.07 |
| 4 | 17.02 | 18.66 | 23.42 | 23.28 | 26.31 | 30.71 | 23.16 | 14.83 |
| 5 | 25.08 | 24.34 | 25.63 | 32.28 | 37.77 | 37.88 | 31.80 | 19.51 |
| 6 | 11.26 | 14.82 | 14.37 | 18.92 | 22.06 | 21.26 | 18.81 | 12.08 |
| 7 | 17.58 | 12.38 | 11.97 | 13.66 | 15.08 | 15.44 | 13.62 | 9.16 |
| 8 | 9.24 | 10.20 | 9.43 | 11.49 | 13.13 | 12.63 | 11.52 | 7.77 |
| 9 | 9.97 | 10.65 | 10.08 | 11.71 | 13.41 | 12.31 | 11.72 | 7.65 |
| 10 | 13.17 | 13.13 | 13.17 | 15.48 | 17.31 | 15.96 | 15.59 | 10.04 |
| 11 | 4.26 | 10.59 | 9.56 | 12.00 | 13.60 | 12.97 | 12.11 | 8.24 |
| 12 | 7.98 | 7.61 | 11.23 | 11.20 | 12.15 | 13.56 | 11.15 | 7.65 |

### 4.6 Group 1 variances — relative discrepancy (%)

| Time | C-α INAD | U-INAD | NB GLM | NB GLMM | tscount | CGFM ind. | CGFM shared |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 24.05 | 24.41 | 21.58 | 55.95 | −2.43 | 22.07 | −28.95 |
| 2 | −8.25 | −20.91 | −21.83 | −7.59 | 50.54 | −21.70 | −50.07 |
| 3 | 43.98 | 35.55 | 27.76 | 41.89 | 91.07 | 28.35 | −18.20 |
| 4 | 9.67 | 37.65 | 36.80 | 54.62 | 80.48 | 36.11 | −12.83 |
| 5 | −2.95 | 2.17 | 28.71 | 50.60 | 51.04 | 26.77 | −22.21 |
| 6 | 31.56 | 27.58 | 68.00 | 95.81 | 88.71 | 67.01 | 7.28 |
| 7 | −29.57 | −31.94 | −22.31 | −14.25 | −12.16 | −22.54 | −47.89 |
| 8 | 10.39 | 2.04 | 24.36 | 42.11 | 36.79 | 24.76 | −15.91 |
| 9 | 6.80 | 1.08 | 17.40 | 34.46 | 23.45 | 17.49 | −23.27 |
| 10 | −0.30 | 0.03 | 17.54 | 31.48 | 21.18 | 18.38 | −23.76 |
| 11 | 148.75 | 124.57 | 181.71 | 219.24 | 204.65 | 184.44 | 93.61 |
| 12 | −4.61 | 40.73 | 40.38 | 52.23 | 69.89 | 39.68 | −4.13 |
| **MARD (all 12)** | **26.74** | 29.06 | 42.37 | 58.35 | 61.03 | 42.44 | 29.01 |
| **MARD (excl. t=11)** | **15.65** | 20.37 | 29.69 | 43.69 | 47.91 | 29.55 | 23.14 |

Constrained-α INAD gives the closest Group 1 variance trajectory.
CGFM shared frailty is competitive (29.01% all 12 cells; 23.14% excluding
t=11), but constrained-α pulls farther ahead on the sensitivity excluding
the small-empirical-variance cell at t=11.

### 4.7 Group 2 means — fitted by model

| Time | Empirical | C-α INAD | U-INAD | NB GLM | NB GLMM | tscount | CGFM ind. | CGFM shared |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 10.17 | 10.28 | 10.33 | 11.47 | 11.80 | 10.14 | 11.47 | 11.25 |
| 2 | 6.49 | 7.53 | 6.81 | 7.03 | 6.97 | 9.41 | 7.03 | 6.96 |
| 3 | 7.86 | 8.04 | 7.50 | 7.76 | 7.41 | 9.12 | 7.77 | 7.77 |
| 4 | 9.29 | 7.17 | 8.17 | 8.37 | 8.07 | 9.28 | 8.33 | 8.41 |
| 5 | 9.63 | 9.24 | 9.66 | 10.08 | 9.88 | 10.71 | 9.97 | 9.97 |
| 6 | 7.37 | 7.40 | 7.40 | 7.43 | 7.30 | 7.39 | 7.39 | 7.39 |
| 7 | 6.60 | 6.56 | 6.18 | 6.13 | 5.86 | 6.14 | 6.11 | 6.17 |
| 8 | 5.74 | 6.10 | 5.76 | 5.53 | 5.40 | 5.56 | 5.53 | 5.53 |
| 9 | 4.91 | 6.21 | 5.91 | 5.59 | 5.47 | 5.65 | 5.58 | 5.47 |
| 10 | 6.34 | 7.20 | 7.10 | 6.60 | 6.35 | 6.68 | 6.62 | 6.55 |
| 11 | 6.26 | 6.28 | 5.96 | 5.67 | 5.51 | 5.74 | 5.70 | 5.75 |
| 12 | 5.89 | 4.84 | 6.04 | 5.44 | 5.16 | 5.32 | 5.42 | 5.47 |

### 4.8 Group 2 means — relative discrepancy (%)

| Time | C-α INAD | U-INAD | NB GLM | NB GLMM | tscount | CGFM ind. | CGFM shared |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1.08 | 1.53 | 12.75 | 16.01 | −0.30 | 12.80 | 10.56 |
| 2 | 16.07 | 5.07 | 8.45 | 7.53 | 45.12 | 8.37 | 7.29 |
| 3 | 2.27 | −4.58 | −1.23 | −5.65 | 16.01 | −1.14 | −1.14 |
| 4 | −22.79 | −11.99 | −9.86 | −13.13 | −0.06 | −10.28 | −9.38 |
| 5 | −4.08 | 0.28 | 4.66 | 2.58 | 11.23 | 3.59 | 3.59 |
| 6 | 0.43 | 0.38 | 0.75 | −1.00 | 0.27 | 0.24 | 0.24 |
| 7 | −0.63 | −6.34 | −7.10 | −11.27 | −6.98 | −7.42 | −6.49 |
| 8 | 6.25 | 0.24 | −3.75 | −6.01 | −3.24 | −3.72 | −3.72 |
| 9 | 26.40 | 20.27 | 13.77 | 11.23 | 14.97 | 13.64 | 11.39 |
| 10 | 13.46 | 11.91 | 4.10 | 0.09 | 5.24 | 4.36 | 3.32 |
| 11 | 0.41 | −4.82 | −9.32 | −11.93 | −8.30 | −8.95 | −8.03 |
| 12 | −17.81 | 2.57 | −7.49 | −12.40 | −9.68 | −7.92 | −7.00 |
| **MARD** | 9.31 | **5.83** | 6.94 | 8.24 | 10.12 | 6.87 | 5.93 |

Unconstrained INAD gives the lowest Group 2 mean MARD (5.83%).
CGFM shared frailty is close (5.93%).

### 4.9 Group 2 variances — fitted by model

| Time | Empirical | C-α INAD | U-INAD | NB GLM | NB GLMM | tscount | CGFM ind. | CGFM shared |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 59.21 | 51.25 | 51.66 | 75.74 | 88.00 | 57.21 | 75.97 | 41.60 |
| 2 | 31.08 | 30.76 | 29.03 | 31.21 | 33.59 | 56.29 | 31.24 | 18.58 |
| 3 | 60.01 | 36.51 | 35.23 | 37.19 | 37.49 | 51.37 | 37.33 | 22.25 |
| 4 | 50.74 | 33.18 | 41.12 | 42.61 | 43.68 | 53.26 | 42.34 | 25.41 |
| 5 | 43.06 | 36.59 | 38.78 | 59.71 | 63.27 | 74.73 | 58.72 | 33.85 |
| 6 | 30.71 | 24.14 | 23.78 | 34.38 | 36.44 | 34.48 | 34.14 | 20.49 |
| 7 | 24.95 | 22.45 | 21.90 | 24.50 | 24.63 | 28.76 | 24.41 | 15.31 |
| 8 | 22.84 | 19.07 | 17.88 | 20.46 | 21.34 | 21.33 | 20.51 | 12.87 |
| 9 | 16.08 | 19.77 | 18.82 | 20.87 | 21.82 | 24.34 | 20.87 | 12.67 |
| 10 | 30.35 | 21.94 | 24.07 | 27.91 | 28.40 | 34.42 | 28.09 | 16.86 |
| 11 | 29.73 | 19.40 | 18.67 | 21.41 | 22.13 | 25.30 | 21.60 | 13.70 |
| 12 | 33.81 | 21.35 | 28.51 | 19.93 | 19.70 | 19.40 | 19.81 | 12.67 |

### 4.10 Group 2 variances — relative discrepancy (%)

| Time | C-α INAD | U-INAD | NB GLM | NB GLMM | tscount | CGFM ind. | CGFM shared |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | −13.43 | −12.74 | 27.93 | 48.63 | −3.37 | 28.32 | −29.74 |
| 2 | −1.03 | −6.58 | 0.43 | 8.09 | 81.10 | 0.50 | −40.22 |
| 3 | −39.15 | −41.30 | −38.02 | −37.53 | −14.40 | −37.78 | −62.92 |
| 4 | −34.60 | −18.97 | −16.03 | −13.91 | 4.97 | −16.55 | −49.92 |
| 5 | −15.02 | −9.95 | 38.65 | 46.91 | 73.54 | 36.36 | −21.39 |
| 6 | −21.40 | −22.56 | 11.96 | 18.67 | 12.29 | 11.17 | −33.27 |
| 7 | −10.02 | −12.22 | −1.80 | −1.31 | 15.24 | −2.19 | −38.63 |
| 8 | −16.50 | −21.75 | −10.44 | −6.57 | −6.62 | −10.22 | −43.68 |
| 9 | 22.96 | 17.02 | 29.77 | 35.68 | 51.37 | 29.76 | −21.24 |
| 10 | −27.70 | −20.68 | −8.04 | −6.41 | 13.40 | −7.45 | −44.44 |
| 11 | −34.75 | −37.19 | −27.99 | −25.55 | −14.88 | −27.33 | −53.90 |
| 12 | −36.86 | −15.69 | −41.04 | −41.72 | −42.62 | −41.40 | −62.54 |
| **MARD** | 22.79 | **19.72** | 21.01 | 24.25 | 27.82 | 20.75 | 41.83 |

Unconstrained INAD gives the lowest Group 2 variance MARD (19.72%);
CGFM independent frailty (20.75%) and NB GLM (21.01%) are close. CGFM
shared frailty is markedly worse here (41.83%) — the shared-frailty
restriction overshrinks the marginal variance for the larger-mean group.

### 4.11 Best-model summary across the four moment cells

| Cell | Best model | MARD |
|---|---|---|
| G1 means | Unconstrained INAD | 5.92% |
| G1 variances | Constrained-α INAD | 26.74% |
| G2 means | Unconstrained INAD | 5.83% |
| G2 variances | Unconstrained INAD | 19.72% |

Pooled across groups:

| Pooled metric | Best model | MARD |
|---|---|---|
| Means (pooled across both groups) | Unconstrained INAD | 5.88% |
| Variances (pooled across both groups) | Unconstrained INAD | 24.39% |

Sensitivity excluding the small-denominator cell at G1 t=11:

| Cell (sensitivity) | Best model | MARD | Note |
|---|---|---|---|
| G1 variances excl. t=11 | Constrained-α INAD | 15.65% | Ranking unchanged; INAD still leads |

### 4.12 Marginal-moment headline

The two INAD variants give the closest overall reproduction of the
empirical marginal mean and variance trajectories. Unconstrained INAD has
the smallest relative discrepancy for the mean trajectories in both
groups and for Group 2 variance, whereas constrained-α INAD gives the
closest Group 1 variance trajectory. The parsimonious constant-α fit
remains competitive on all four cells, consistent with its BIC-supported
role as the primary bolus model; allowing mild time variation in α
through the unconstrained fit improves marginal-mean tracking,
consistent with the LRT result for non-constancy of α. The advantage is
aggregate rather than pointwise, since individual visit-specific moments
are occasionally closer under the comparator models.

## 5. Cross-Cutting Findings

### 5.1 Plug-in vs posterior predictive

INAD predictions throughout are Monte Carlo estimates of the
**conditional predictive PMF at the MLE plug-in**, not posterior
predictives. NB GLM and Poisson GLM baselines use **analytic plug-in
PMFs** at the MLE; NB GLMM and tscount use analytic marginalization or
Monte Carlo recursion respectively. All models share the same plug-in
convention — what differs is only the predictive-PMF construction.
This is the standard convention for `predict()` in maximum-likelihood
frameworks. A future extension would propagate parameter uncertainty
(parametric bootstrap or Hessian-based draws) — useful on small-n
datasets like bolus but not required for the current scope.

### 5.2 Integer-valued property of predictions

The predictive **distribution** is integer-valued by construction (every
MC path is integer). The predictive **mean** is real-valued by definition
(the expectation of an integer random variable is not generally an
integer) and is the right object for RMSE/MAE. Rounding `μ̂` would
introduce bias at small means and would not change probabilistic scores.
Integer point forecasts, if needed, are the predictive mode or median of
the MC sample.

### 5.3 Why log score behaves differently from RPS and RMSE

Log score depends on the predictive PMF at a single point (the observed
value), so it is sensitive to the smoothing constant `ε` applied to the
empirical PMF and to tail behavior of the simulated support. RPS and
RMSE integrate across the predictive support and are more robust. This
explains why log score can show mixed signs across horizons (especially
at recursive multi-step) even when RPS and RMSE consistently favor INAD.
The §3.5 sensitivity sweep on the bolus 35-fold rotation made this
behaviour explicit; the all-pairs result in §3.2 is INAD-favourable on
log score across all cells, consistent with the smoothing diagnosis at
the larger sample size.

### 5.4 Tau SE and CI consistency

Tau confidence intervals are computed by profile likelihood; for
consistency, tau SEs are now derived from the profile-CI half-width
(`SE = (upper − lower) / (2 · qnorm(0.975))`). The `est ± 2·SE ≈ (lb, ub)`
relationship therefore holds for tau by construction. Other INAD
parameters (α, θ, nb_inno_size) continue to use mutually consistent Wald
SE and Wald CI.

Test-coverage gap: the per-component tau SE-from-CI logic is currently
exercised only on the two-group bolus case (`tau[2]` alone). A simulated
3+ group unit test is still pending so the logic is verified beyond
`tau[2]`.

### 5.5 Engineering validation

- All R scripts parse cleanly; all simulation and real-data runs completed
  without fit failures across R = 1000 simulated reps, 1050 bolus folds
  for the predictive comparison, and the full-data fits used for the
  marginal-moment comparison.
- Run times: simulation R = 1000 ≈ 5.6 h; bolus all-pairs ≈ 10.7 h
  (38591 s with 10 workers); marginal-moment driver minutes.
- Patient-level clustered SEs correctly account for fold dependence in
  the all-pairs design; fold-level SEs would have understated uncertainty.

## 6. Implications for the `antedep` Package

### 6.1 Landed

- **Internal forward simulator.** `.simulate_inad_forward()` is in
  `R/predict_inad.R` as an internal helper. It seeds the recursion from
  supplied history, applies `tau[blocks[i]]` per subject, and errors
  when `tau` is non-trivial but `blocks` is `NULL`. Covered by
  `tests/testthat/test-predict-inad.R`.
- **`predict.inad_fit()` v1 — public registered S3 method.** Defined
  in `R/predict_inad.R`, exported via `@export`, and registered through
  `S3method(predict, inad_fit)` in `NAMESPACE`. Supports
  `type = c("mean", "sample")` for order-1 INAD fits.

  ```r
  predict(object, newdata = NULL, blocks = NULL, h = 1L,
          type = c("mean", "sample"),
          n_sims = 1000L, seed = NULL, ...)
  ```

  - `newdata = NULL` returns in-sample fitted predictions over the last
    `h` training occasions. Non-`NULL` `newdata` returns conditional
    forecasts for new subjects; `blocks` is required iff the fit has
    non-zero block effects.
  - `type = "mean"` returns an `n_test × h` matrix of MC predictive
    means; `type = "sample"` returns an `n_test × h × n_sims` integer
    array of raw paths.
  - Documented as plug-in MLE-conditional — same convention as
    `stats::predict.glm`.

### 6.2 Remaining

- **Distribution API (`predict.inad_fit()` v1.1).** Adds `type =
  "distribution"` returning a smoothed empirical PMF over a finite
  support derived from the MC samples. Held until the gau/cat POCs
  settle the shared cross-family contract for distribution-valued
  returns, so the same API shape lands consistently across all three
  families.
- **Proper `alpha_constraint = "constant"` in `fit_inad()`.** The
  constrained-α fit currently exposed only through the POC-internal
  wrapper inherits `nb_inno_size` from the unconstrained fit. The
  user-facing version should be a proper path in `fit_inad()` with
  joint re-estimation of `(α, θ_t, nb_inno_size_t, τ)` under the
  constraint. Timing: after the gau/cat POCs settle the cross-family
  API shape for structural constraints. Skip if the constrained-α
  story stops being central.
- **Exact one-step PMF helpers — future work, not urgent.** Rolling
  one-step bolus is already a clean INAD win without exact PMFs (§3.2),
  so an exact one-step PMF would provide methodological cleanup rather
  than a result-changing improvement. Multi-step exact PMFs are harder
  to derive in closed form and are deferred.
- **Scoring helpers — decide later.** `log_score()`, `rps()`,
  `randomized_pit()` for count predictions on finite support. Useful if
  `predict()` returns distributions; not strictly required if users
  depend on `scoringRules` directly. Decide together with the v1.1
  `type = "distribution"` work.

### 6.3 Sequencing for the remaining items

1. **Start gau and cat prediction POCs in parallel** while INAD
   prediction machinery is fresh. Reuse the existing
   `scripts/prediction_poc/common/` scoring infrastructure (see §7).
2. **After the gau/cat POCs**, settle the shared cross-family
   `type = "distribution"` API. Then add INAD v1.1, gau, and cat
   distribution returns under a single consistent contract.
3. **Implement proper `alpha_constraint = "constant"` in `fit_inad()`**
   — only after step 2 reveals whether parallel constraint mechanisms
   are needed for gau and cat.
4. **Defer exact one-step PMFs** as documented future work.

## 7. Bridge to Gaussian and Categorical Modules

The infrastructure built for INAD generalizes cleanly to `gau_fit` and
`cat_fit`. The forward simulator is the family-specific part; the
scoring driver, finite-support smoothing, patient-level clustering, and
report shape are reusable across families.

For Gaussian:

- One-step and h-step predictive distributions are closed-form
  (linear-Gaussian state-space recursion). No Monte Carlo needed for
  the mean or variance.
- Validation should compare the analytic recursion against MC from
  `simulate_gau()` once, then use analytic.
- Baselines: pooled lag regression, time-marginal mean, LMM with AR
  errors.

For categorical:

- One-step is a transition lookup; h-step is exact recursive
  matrix-vector multiplication on the lag `p`-tuple distribution. No
  Monte Carlo.
- Smoothing/backoff policy needed for unseen lag tuples in the test set.
- Baselines: multinomial logistic regression with lag dummies, Markov
  chain, time-marginal proportion.

Both should reuse the existing `scripts/prediction_poc/common/` scoring
infrastructure (the `ppc_score_pmf`, `ppc_score_samples` helpers handle
the discrete-support case; Gaussian needs a continuous analogue but the
clustered-bootstrap and per-fold patient aggregation generalize
unchanged).

The same evaluation design — conditional new-subject forecasting,
rolling one-step + recursive multi-step, stratified train/test split,
R = 1000 simulation + a real-data application with patient-level
clustering — should carry across.
