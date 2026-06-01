---
title: INAD Prediction Metrics and CGFM Comparison Summary
---

This note gives a compact, presentation-oriented summary of prediction metrics,
the Henderson-Shimakura Table 4 reproduction, and the bolus rolling one-step
prediction comparison. The main comparison is between integer-valued
antedependence (INAD), negative-binomial generalized linear model (NB GLM),
negative-binomial generalized linear mixed model (NB GLMM), and correlated
gamma frailty model (CGFM) predictors.

# Supplement To Last Week's Validation Summary

## Marginal-Moment MARD Table

This section adds the CGFM time-varying frailty row to the bolus marginal
moment MARD table from last week's validation summary. MARD means mean
absolute relative discrepancy, reported as a percentage; lower is better.

The NB GLM and CGFM independent frailty models have the same one-time
negative-binomial predictive family under the same time and group linear
predictor. The table therefore uses a single NB GLM / CGFM independent row
based on the reproduced Henderson-Shimakura independent-frailty moments.

\begin{center}
\small
\begin{tabular}{lrrrrr}
\hline
Model & G1 mean & G1 var & G1 var excl. & G2 mean & G2 var \\
\hline
C-alpha INAD & 8.35 & 26.74 & 15.65 & 9.31 & 22.79 \\
U-INAD & 5.92 & 29.05 & 20.37 & 5.83 & 19.72 \\
Poisson GLM & 9.95 & 59.70 & 64.84 & 5.96 & 78.52 \\
NB GLM / CGFM ind. & 8.94 & 42.44 & 29.53 & 6.87 & 20.75 \\
NB GLMM & 8.35 & 58.35 & 43.73 & 8.24 & 24.25 \\
CGFM shared & 9.86 & 29.01 & 23.14 & 6.01 & 41.83 \\
CGFM time-varying & 8.96 & 42.26 & 29.70 & 6.88 & 21.08 \\
\hline
\end{tabular}
\end{center}

Adding the CGFM time-varying row does not change the table's main reading:
the INAD class still gives the best MARD in the table, with unconstrained
INAD best for G1 mean, G2 mean, and G2 variance, and constrained-alpha INAD
best for G1 variance.

# Henderson-Shimakura Table 4 Reproduction

Henderson and Shimakura (2003), Table 4, was reproduced from the bolus data to
check the CGFM reference models.

## Model-Level Checks

\begin{center}
\small
\begin{tabular}{p{0.22\textwidth}rrp{0.43\textwidth}}
\hline
Model & Reproduced LL & Paper LL & Key parameter check \\
\hline
Independent frailty & -2191.77 & -2191.77 &
\(\beta = 0.34\); \(\theta = 0.49\); information SE and robust SE match printed values. \\
Shared frailty & -2249.16 & -2249.16 &
\(\beta = 0.33\); \(\theta = 0.24\); SEs match printed values. \\
Time-varying frailty & -2159.65 & -2159.5 &
\(\beta = 0.3428\), \(\theta = 0.4907\), \(\rho = 0.8487\); parameter estimates match Table 4 to printed precision. \\
\hline
\end{tabular}
\end{center}

The full time-varying fit converged from data and recovered \(\rho \approx
0.85\), so the paper's reported correlation parameter is reproducible rather
than merely paper-calibrated.

The likelihood column verifies reproduction within each model row; it is not a
cross-model likelihood comparison. The independent and shared rows are full
log-likelihoods, whereas the time-varying row is the normalized pairwise
composite log-likelihood used by Henderson and Shimakura. The time-varying
objective value is close to, but not an exact one-decimal reproduction of, the
paper value (-2159.65 from data versus -2159.5 reported).

## Time-Varying Row Reproduced From Data

| Parameter | Script estimate | Script SE | Paper Table 4 |
| --- | ---: | ---: | ---: |
| alpha1 | 2.1043 | 0.10 | 2.10 (0.10) |
| alpha2 | 1.6164 | 0.13 | 1.62 (0.13) |
| alpha3 | 1.7118 | 0.12 | 1.71 (0.12) |
| alpha4 | 1.7898 | 0.11 | 1.79 (0.11) |
| alpha5 | 1.9664 | 0.10 | 1.97 (0.10) |
| alpha6 | 1.6621 | 0.10 | 1.66 (0.10) |
| alpha7 | 1.4741 | 0.13 | 1.47 (0.13) |
| alpha8 | 1.3701 | 0.12 | 1.37 (0.12) |
| alpha9 | 1.3743 | 0.11 | 1.37 (0.11) |
| alpha10 | 1.5398 | 0.11 | 1.54 (0.11) |
| alpha11 | 1.3910 | 0.10 | 1.39 (0.10) |
| alpha12 | 1.3479 | 0.12 | 1.35 (0.12) |
| beta | 0.3428 | 0.13 | 0.34 (0.13) |
| theta | 0.4907 | 0.06 | 0.49 (0.06) |
| rho | 0.8487 | 0.03 | 0.85 (0.03) |

# Prediction Metrics

Let \(Y_i\) be the observed count for prediction case \(i\), let
\(\hat m_i = E(Y_i \mid \mathcal H_i)\) be the predictive mean, and let
\(\hat p_i(y) = P(Y_i = y \mid \mathcal H_i)\) be the predictive probability
mass function conditional on the available history \(\mathcal H_i\). Let
\(\hat F_i(k) = \sum_{y \le k} \hat p_i(y)\) be the predictive CDF.
RPS is the primary metric for the bolus comparison below. Log score is the
second distributional score, and RMSE is kept as a point-forecast check.

## RPS

**In words.** RPS, the ranked probability score, measures the full predictive
CDF against the realized count. For count data, it sums the squared distance
between the predictive CDF and the step function induced by the observed
count. Lower RPS is better. Compared with log score, RPS is less driven by a
single tail probability and is often more stable for discrete empirical PMFs.

**Formula.**

For prediction case \(i\), let the finite scoring support be
\(k = 0, 1, \ldots, K_i\). Then

\[
\mathrm{RPS}_i
= \sum_{k=0}^{K_i}
\left\{
  \hat F_i(k) - \mathbf 1(Y_i \le k)
\right\}^2.
\]

The reported RPS is the average over prediction cases:

\[
\mathrm{RPS}
= \frac{1}{n} \sum_{i=1}^n \mathrm{RPS}_i.
\]

In this validation work, \(K_i\) is not a single global maximum. It is a
case-specific finite support chosen to include the observed value and the
predictive support used for that model's scoring calculation. For
simulation-based INAD forecasts this means the support extends to the maximum
of the observed value and the simulated predictive samples. This finite-support
choice is why the support/smoothing caveat matters for RPS and log score.

## Log score

**In words.** Log score measures how much probability the predictive
distribution assigns to the value that actually occurs. It is a proper
probabilistic scoring rule: a model is rewarded for assigning high probability
to the observed count and penalized for putting too little mass there. Lower
log score is better: the score is defined as \(-\log \hat p\), so smaller
values mean the model placed more probability on the observed count. For
simulation-based empirical PMFs, log score can be sensitive to tail smoothing
because it evaluates the PMF at one point.

**Formula.**

\[
\mathrm{LogScore}
= - \frac{1}{n} \sum_{i=1}^n \log \hat p_i(Y_i).
\]

For an individual prediction case,

\[
\mathrm{LogScore}_i = -\log \hat p_i(Y_i).
\]

## RMSE

**In words.** RMSE measures point-forecast accuracy. It compares the predictive
mean with the observed count, squares the error, averages over prediction
cases, and then returns to the original count scale by taking the square root.
Lower RMSE is better. Because errors are squared, RMSE penalizes large misses
more heavily than MAE.

**Formula.**

\[
\mathrm{RMSE}
= \left\{ \frac{1}{n} \sum_{i=1}^n (\hat m_i - Y_i)^2 \right\}^{1/2}.
\]

# Additional Prediction Metrics

Prediction comparisons usually use a mix of point, distributional, and
calibration metrics. The following metrics are useful context but are not the
main bolus summary metrics reported below.

| Metric | Type | What it measures | Useful here? |
| --- | --- | --- | --- |
| MAE | Point forecast | Mean absolute error, \(\frac{1}{n}\sum_i \lvert \hat m_i - Y_i \rvert\). | Redundant with RMSE for the present summary; useful as a robustness check. |
| PIT / randomized PIT | Calibration | Whether predictive distributions are calibrated; well-calibrated PIT values should be approximately uniform on \([0,1]\). | Useful as a diagnostic plot or calibration check, but not used as a one-number headline metric here. |
| Brier score | Probabilistic | Squared error for binary or categorical event probabilities. | More natural for categorical outcomes than count outcomes. |
| Interval coverage | Calibration | Whether nominal prediction intervals cover at the advertised rate. | High-value future supplement once coverage is computed consistently for all displayed models. |
| Interval score | Probabilistic / interval | Rewards narrow intervals with correct coverage and penalizes misses. | High-value future supplement if prediction intervals are added. |
| CRPS | Probabilistic | Continuous analogue of RPS. | Useful robustness check if a continuous approximation is desired. |
| Dawid-Sebastiani score | Probabilistic moment score | Uses predictive mean and variance. | Possible supplement, but less informative than full-PMF scores for counts. |

In plain language, PIT asks whether the predictive distribution is calibrated:
after randomization for discreteness, observations should fall throughout their
predictive distributions as if they were uniform draws on \([0,1]\). A PIT
histogram or interval-coverage table would be the preferred calibration display;
the present summary does not use PIT as a one-number comparison.

For the bolus comparison, the reported metric set is **RPS (primary), log
score, and RMSE**. RPS and log score evaluate the predictive distribution;
RMSE records whether the distributional gains come with comparable point
accuracy. Ad hoc folded diagnostics are omitted because they are not proper
scores and can be hard to interpret as prediction evidence.

# Bolus Rolling One-Step Prediction Comparison

Design: the bolus data comprise 65 patients in two treatment groups of size 30
and 35. Each fold holds out one patient from each group and fits on the
remaining 63, so the all-pairs leave-one-per-group-out design gives
\(30 \times 35 = 1050\) folds. Prediction is rolling one-step at
\(t = 2,\ldots,12\).

The NB GLM row is aligned with the independent-CGFM reference: under the same
time and group linear predictor, the independent Poisson-gamma frailty model
integrates to the same negative-binomial predictive family as the NB GLM.

The \(t\)-statistics are computed after aggregating to held-out patients. This
is more conservative than treating all 1050 all-pairs folds as independent,
because each patient appears in many folds.

For this comparison, RPS is treated as the primary prediction metric: it is a
proper distributional score for count-valued predictive CDFs and is less
sensitive than log score to one observed tail probability. Log score is the
second distributional score, and RMSE is the point-forecast check.

The table reports paired differences:

\[
\Delta = \text{score(constrained-alpha INAD)} - \text{score(reference)}.
\]

Negative values mean constrained-alpha INAD has the lower, better score. The
number in parentheses is the patient-level paired \(t\)-statistic.
Here, constrained-alpha INAD denotes the INAD fit that keeps the thinning
parameter \(\alpha\) constant across time.

| Model | Delta RPS (t) | Delta RMSE (t) | Delta log score (t) |
| --- | ---: | ---: | ---: |
| Unconstrained INAD | -0.012 (-0.6) | -0.028 (-0.9) | -0.006 (-0.7) |
| Poisson GLM | -0.081 (-2.9) | -0.013 (-0.4) | -0.181 (-5.4) |
| NB GLM | -0.405 (-3.9) | -0.627 (-4.2) | -0.124 (-3.3) |
| NB GLMM | -0.158 (-2.2) | -0.164 (-1.7) | -0.059 (-1.99) |
| CGFM shared | -0.072 (-1.1) | 0.015 (0.2) | -0.144 (-3.2) |
| CGFM time-varying | -0.066 (-2.3) | 0.020 (0.5) | -0.017 (-1.3) |

The absolute score scale is shown below for context. Lower values are better;
RMSE is averaged at the held-out-patient level, matching the paired comparison
above.

| Model | RPS | RMSE | Log score |
| --- | ---: | ---: | ---: |
| Constrained-alpha INAD | 2.207 | 3.870 | 2.673 |
| Unconstrained INAD | 2.219 | 3.898 | 2.679 |
| Poisson GLM | 2.288 | 3.883 | 2.854 |
| NB GLM | 2.612 | 4.497 | 2.797 |
| NB GLMM | 2.365 | 4.034 | 2.732 |
| CGFM shared | 2.279 | 3.855 | 2.817 |
| CGFM time-varying | 2.273 | 3.851 | 2.690 |

The NB GLM row also represents the independent-CGFM predictive family under
the same time and group linear predictor. RPS and log-score values are
sensitive to finite-support and smoothing conventions; this sensitivity is
most relevant for the closest RPS comparisons, especially CGFM time-varying.
The much larger RPS gain over NB GLM is least exposed to this issue; the
Poisson and NB GLMM RPS gains should still be read with the support-convention
caveat in mind.

## Reading The Table By Metric

**RPS.** This is INAD's strongest broad metric. Constrained-alpha INAD is
indistinguishable from unconstrained INAD and better than Poisson GLM, NB GLM,
NB GLMM, and CGFM time-varying by the conventional
\(|t| > 2\) threshold. It is directionally but not significantly better than
CGFM shared. If the primary RPS family of reference comparisons is adjusted
conservatively for multiplicity, the strongest RPS wins, especially NB GLM
and Poisson, remain the most secure; the NB GLMM and time-varying CGFM RPS
gains should be read as supportive rather than standalone proof.

**Log score.** Constrained-alpha INAD is comparable to unconstrained INAD. It
clearly beats Poisson GLM, NB GLM, and CGFM shared. Against NB GLMM, the
log-score advantage is borderline (t = -1.99, approximately p = 0.05). It is
essentially tied with CGFM time-varying.

**RMSE.** Constrained-alpha INAD is comparable to unconstrained INAD, Poisson
GLM, CGFM shared, and CGFM time-varying. It clearly beats NB GLM and is
directionally better than NB GLMM but not quite at the conventional threshold.
The defensible point-forecast statement is therefore comparable point
accuracy, not uniformly better point prediction.

# Summary Conclusion

The current evidence supports the following prediction statement: INAD's
advantage is in probabilistic distributional one-step prediction, not in
uniformly better point forecasts.

1. **INAD's clearest real-data advantage is probabilistic distributional
   prediction, especially RPS.** It improves RPS over Poisson GLM, NB GLM,
   NB GLMM, and time-varying CGFM, with the strongest support against NB GLM;
   the CGFM shared comparison is directionally favorable but not significant.
2. **The distributional evidence is coherent against the simpler
   count-regression baselines.** Against Poisson and NB GLM, INAD improves
   both RPS and log score. Against NB GLMM, the RPS result is clearer than the
   borderline log-score result.
3. **Constrained-alpha and unconstrained INAD are practically tied.** The
   constrained-alpha fit is slightly better on RPS, RMSE, and log score, but
   none of those differences is statistically meaningful; it remains the
   primary row because it is the validation parameterization used throughout.
4. **CGFM time-varying is the toughest comparator.** INAD has the better RPS,
   but ties on RMSE and log score. CGFM shared is also close on RPS. Until a
   common-support sensitivity check is fully documented, the closest RPS
   margins should be framed as supportive rather than overclaimed.
5. **Point-error metrics do not support a universal prediction win.** RMSE
   shows comparable point accuracy outside the NB GLM comparison. This result
   should be acknowledged rather than hidden.

This is a one-dataset, rolling-one-step validation on the bolus data; recursive
multi-step prediction and fully standardized support/smoothing sensitivity for
RPS and log score remain follow-up checks rather than settled evidence here.
