# Prediction for Antedependence Models — Design and Script-First Plan

This note plans a prediction capability for the three model families
(`gau_fit`, `cat_fit`, `inad_fit`). The plan is staged: a proof-of-concept
(POC) script set first, then in-package `predict.*` methods. The promotion
criteria are about *correctness, calibration, stability, and speed* — not
about beating every baseline on real data.

## 1. Motivation

`predict()` is missing from the standard-methods story the package now tells.
Adding it would:

- mirror `lm`, `glm`, `arima`, and friends with a familiar entry point;
- expose the model's own forecast distribution so users can compare against
  GLM, AR, and naive baselines honestly; and
- support natural longitudinal workflows (next-visit prediction, conditional
  forecasting for new subjects).

## 2. Evaluation framework

### 2.1 Evaluation designs

Name and define each design honestly.

- **Conditional new-subject forecasting (primary).** Fit parameters on
  training subjects using all `T` time points. Reveal observations
  `Y_{i,1:T_split}` for held-out test subjects. Predict
  `Y_{i,(T_split+1):T}` from the model's conditional distribution. This is
  the natural design for the package because it uses the time-indexed
  parameter estimates without extrapolation, and it isolates "does the lag
  structure carry signal for new subjects?" from any stationarity
  assumption.
- **Temporal holdout (only under explicit assumptions).** Fit on times
  `1..T_train`; predict times `T_train+1..T`. The current fitters
  produce time-indexed parameters, so this design is not generally
  feasible. It is valid only when one of the following is explicit and
  documented per experiment:
  1. a stationary / time-invariant fit is used (note: in this package
     "homogeneous" means shared-across-blocks, not time-invariant — the
     relevant property here is time-invariance, tested by
     `test_stationarity_*`);
  2. the last estimated parameter is frozen forward;
  3. a parametric time trend is fit and extrapolated; or
  4. it is labelled a heuristic, not a standard AD forecast.
- **Walk-forward / expanding window.** Repeated temporal holdout. Inherits
  the same caveats; useful only when one of the four assumptions above
  holds.

### 2.2 Horizon and predictive distribution

- Rolling one-step ahead: at each held-out time, condition on the observed
  history through the previous time point. Gaussian — closed-form Gaussian predictive;
  categorical — exact transition lookup; INAD — closed-form conditional
  mean, predictive distribution by Monte Carlo for the POC (exact
  one-step PMFs are tractable for some thinning/innovation combinations
  but deferred per §5.4).
- Recursive multi-step ahead: condition on observed history through a split
  time, then propagate future lags from the predictive distribution.
  Gaussian uses state-space recursion; categorical uses exact discrete
  recursion; INAD uses Monte Carlo.
- Always produce a *predictive distribution*, not just a mean. Probabilistic
  scores (log-score, CRPS, RPS, Brier) are the headline metric; point
  scores (RMSE, MAE) are reported alongside.

### 2.3 Uncertainty quantification

Independent unit = subject. Aggregate `(subject, time)` scores within the
held-out set; report paired differences vs the best baseline with a block
bootstrap CI over subjects.

## 3. Gaussian module

### 3.1 Model

For order `p`,
$$Y_{i,t} = \mu_t + \tau_{b(i)} + \sum_{k=1}^p \phi_{k,t}\,(Y_{i,t-k} - \mu_{t-k} - \tau_{b(i)}) + \varepsilon_{i,t},
\quad \varepsilon_{i,t} \sim \mathcal{N}(0, \sigma_t^2).$$

### 3.2 One-step predictive distribution (closed form)

$$\hat\mu_{i,T+1} = \mu_{T+1} + \tau_{b(i)} + \sum_{k=1}^p \phi_{k,T+1}\,(y_{i,T+1-k} - \mu_{T+1-k} - \tau_{b(i)}),
\qquad \hat\sigma^2_{i,T+1} = \sigma^2_{T+1}.$$

Start the POC here.

### 3.3 Multi-step predictive distribution (state-space recursion)

The mean *can* be iterated by plugging predicted values into future lag
slots, but the predictive *variance* must propagate through the lag state.
Stack the last `p` values into a state vector
$\mathbf{X}_t = (Y_t, Y_{t-1}, \ldots, Y_{t-p+1})^\top$ and write the AD
recursion as a linear-Gaussian state-space update:
$$\mathbf{X}_{t+1} = F_{t+1} \mathbf{X}_t + G_{t+1} + \mathbf{w}_{t+1},
\quad \mathbf{w}_{t+1} \sim \mathcal{N}(0, Q_{t+1}).$$

Predictive moments at horizon `h` come from iterating `F`, `G`, and `Q`.
Validate the analytic recursion against a Monte Carlo run from
`simulate_gau()` driven by the fitted parameters before reporting any
multi-step result.

### 3.4 Baselines

Keep the starting set light. Heavier baselines (per-subject AR, mixed-effects
with AR errors) only after the POC is producing scores.

| Baseline | Purpose |
|---|---|
| Last observation carried forward (`ŷ = y_T`) | Floor |
| Time-marginal mean (`ŷ = μ̂_t` from training) | Tests temporal-structure value |
| Pooled lag regression (single `φ` across subjects, OLS) | Practitioner default |

### 3.5 Metrics

RMSE, MAE, Gaussian log-score (analytic), CRPS for Gaussian
(`scoringRules::crps_norm`), 80% / 95% predictive-interval coverage, PIT
histogram.

### 3.6 Gotchas

- Block effects must be carried into the predictor for held-out subjects.
  Prefer stratifying the conditional new-subject split by block so every
  block appears in both train and test; drop any unavoidable
  unseen-block subjects and report how many.
- Temporal holdout requires one of the §2.1 assumptions to be made explicit.

## 4. Categorical module

### 4.1 Model

For order `p` and `K` categories, the conditional law
$$P(Y_{i,t} = j \mid Y_{i,t-1}, \ldots, Y_{i,t-p}) = \pi_{j \mid y_{i,t-1:t-p},\, t}$$
is encoded in `cat_fit$transition`.

### 4.2 Transition extraction across fit variants

`cat_fit` stores parameters differently across the cases the fitter handles.
The POC needs one helper that returns transition arrays in a canonical shape
for all four:

| Case | Storage |
|---|---|
| homogeneous, complete data | single `marginal`, single `transition` array |
| heterogeneous (blocks), complete data | per-block lists `marginal[[g]]`, `transition[[g]]` |
| homogeneous, marginalized (EM-based) | same shape as homogeneous complete; settings flag differs |
| heterogeneous, marginalized | same shape as heterogeneous complete; settings flag differs |

Add `.extract_cat_transitions(fit)` in the script `common/` layer; if it
proves correct, it becomes a candidate internal helper in `R/`.

### 4.3 One-step and multi-step (exact)

One-step: lookup `π[\cdot, j, T+1]` indexed by the realized lag `p`-tuple.

Multi-step (exact, no Monte Carlo): maintain a length-`K^p` predictive
distribution over the current lag `p`-tuple of states. At each forecast
step, left-multiply by the time-`(T+s)` transition tensor and marginalize
the oldest slot. Cost `O(h \cdot K^{p+1})` per subject, fine for `K`
small and `p ≤ 2`.

### 4.4 Fallback for unseen lag tuples

The test set can present lag tuples that have zero estimated
probability — log score becomes `-Inf` if not handled. Choose one policy
per experiment and document it:

1. **Laplace smoothing** with a small `α` added to each transition cell.
2. **Marginal backoff**: substitute the marginal `P(Y_t = j \mid t)` when
   the conditioning tuple is unseen.
3. **Previous-order backoff**: fall back to the AD(`p-1`) transition
   conditioned on the available tail of the tuple.

Recommend (1) for the POC with `α = 1 / n_train_obs` and report sensitivity
in the report.

### 4.5 Baselines

| Baseline | Purpose |
|---|---|
| Last-state (`ŷ = Y_{i,T}`) | Floor |
| Time-marginal proportion | Tests lag value |
| Multinomial logistic with lag dummies (`nnet::multinom`) | Practitioner default |

### 4.6 Metrics

Multiclass log-loss, multiclass Brier score, per-class reliability diagram,
RPS if the categories are ordinal. Top-1 accuracy is reported but not
treated as the headline.

## 5. INAD module

### 5.1 Model

For order `p`,
$$Y_{i,t} = \sum_{k=1}^p \alpha_{k,t} \circ Y_{i,t-k} + \varepsilon_{i,t},$$
where `∘` is binomial / Poisson / NB thinning and `ε_{i,t}` is
Poisson / Bell / NB.

The INAD POC starts with one fixed data-generating structure:
NBT-NBI-INADFE(1), meaning negative-binomial thinning, negative-binomial
innovation, an order-1 lag, and fixed block/treatment effects. Because the
simulation structure is known, this POC does **not** do model selection.

`fit_inad()` estimates a time-varying `α_t` by default; the DGP uses a
constant `α = 0.35`. The POC therefore runs the comparison **two ways**
so that "fit the known structure" is an honest claim and not a
hand-wave:

- **Unconstrained fit (primary).** Use `fit_inad()` as-is; `α_t` is
  estimated freely. Measures realistic user-workflow performance.
  Whether `α̂_t` recovers a near-constant pattern is itself a diagnostic.
- **Constrained-α fit (supplementary).** A POC-local helper in
  `scripts/prediction_poc/common/` either wraps a fitter that holds
  `α` constant or adapts the constrained machinery from
  `test_stationarity_inad(..., constrain = "alpha")` into an
  `inad_fit`-shaped object. **Implementation caveat:** the currently
  available constrained machinery fixes `nb_inno_size` from the
  unconstrained fit rather than jointly re-estimating it under the
  constrained-`α` model, so the supplementary fit **matches the DGP's
  constant-α structure** but does not jointly re-estimate every
  parameter under that structure unless the helper extends to a joint
  fit over `(α, θ_t, n̂b_inno_size_t, τ)`. The constancy of `α` is the
  load-bearing comparison; full joint re-estimation is an optional
  extension worth adding if early scores suggest the
  `nb_inno_size`-from-unconstrained step is doing real damage.

The two fits share the forward simulator, baselines, and scoring layer.
The difference between their scores isolates the cost of having to
estimate `α`'s time-variation when in truth it is constant.

Both fits are then evaluated on the two prediction tasks from §2.2:

1. rolling one-step conditional prediction at `t = 2, 3, ..., 12`,
   conditioning on the realized history `Y_{i, 1:(t-1)}` — independent
   of any split time, used as the across-time calibration diagnostic; and
2. recursive multi-step conditional prediction from split time
   `T_split = 8`, evaluating `Y_{i, 9}` (`h = 1`), `Y_{i, 10}`
   (`h = 2`), `Y_{i, 11}` (`h = 3`), and `Y_{i, 12}` (`h = 4`) drawn
   recursively from the predictive distribution — the applied
   "forecast the next four visits" task.

Both tasks live inside the conditional new-subject design from §2.1:
parameter estimates at every `t` come from training subjects fit on the
full 12-point grid, so no temporal extrapolation is required at any
horizon.

Only after both tasks behave correctly in this simulation should the same
two prediction tasks be run on `bolus_inad` and compared with GLM baselines.

Simulation setting for the first INAD POC:

- Subjects and time: use two groups and 12 time points to mirror
  `bolus_inad`; start with balanced groups, e.g. `n_per_group = 50`
  (`N = 100` total), then vary `N` only in sensitivity checks.
- Thinning: use a time-homogeneous negative-binomial thinning parameter,
  because the bolus analysis suggests a constant alpha is adequate. Start
  with `alpha = 0.35`; this is large enough to create visible lag
  dependence while staying comfortably below the package's `alpha < 1`
  constraint for NB thinning.
- Innovation mean (`theta`, the NB mean parameter): allow time variation.
  A reasonable starting vector for 12 time points is
  `c(2.2, 2.4, 2.8, 3.2, 3.6, 3.3, 3.0, 2.7, 2.4, 2.1, 1.9, 1.7)`.
- Innovation size (`nb_inno_size`, the NB dispersion/size parameter):
  allow time variation, e.g.
  `c(1.4, 1.4, 1.6, 1.8, 2.0, 2.0, 1.9, 1.8, 1.7, 1.6, 1.5, 1.5)`.
  These values keep moderate overdispersion without making the simulation
  dominated by extreme counts.
- Group effect: use two groups with a distinguishable additive innovation
  effect, e.g. `tau = c(0, 1.25)`. This should create a visible treatment
  difference while keeping effective NB means positive at all times.
- Replication: run `R = 1000` simulation replicates. Use parallel computing
  with reproducible seeds and a worker count derived from the local device,
  e.g. `n_workers = max(1, floor(2 * parallel::detectCores(logical = TRUE) / 3))`.
  Record the detected core count, worker count, RNG seed, and elapsed time
  in the POC output.

**Subject split.** Conditional new-subject design per §2.1, stratified by
group so both treatment groups appear in train and test. Default 70/30
stratified split, drawn fresh per rep with a rep-indexed seed.

**Sanity check before scoring (parameter recovery).** On each rep,
record per-fit recovery against the DGP:

- Unconstrained fit: compare `α̂_t` (length-12 vector) to the scalar
  `0.35` elementwise; report MAE and the across-`t` shape (mean,
  min, max) to surface whether the unconstrained fitter recovers a
  near-constant pattern. Compare `θ̂_t`, `n̂b_inno_size_t`, and `τ̂` to
  their DGP values (bias and RMSE per time / scalar bias and RMSE for
  `τ`).
- Constrained-α fit: compare scalar `α̂` to `0.35` (bias and RMSE).
  Same diagnostics for `θ̂_t`, `n̂b_inno_size_t`, and `τ̂`.

Reported alongside predictive scores so a fitter problem and a predictor
problem can be told apart: if recovery is fine but predictive scores are
bad, the predictor is at fault; if recovery is bad, fix the fitter
before interpreting any predictive comparison.

Early smoke diagnostics found two useful side results. First, increasing
`n_per_group` from 20 to 50 removed the spurious high-`nb_inno_size` basin,
so NB-size recovery is treated as a sample-size diagnostic rather than an
iteration-budget artifact. Second, the constrained-α fit reduced `θ_t`
RMSE relative to the unconstrained fit, suggesting that fixing α at a
consistent constant estimate can concentrate information on the innovation
mean even if predictive scores are not yet clearly better. A one-rep
default-initialization check at `n_per_group = 50` converged to the same
solution as truth-based initialization, so initialization is not expected
to be a major confounder; mention this only briefly in the report.

**Aggregation across `R = 1000` reps:**

- Per-rep paired score differences (INAD minus NB GLM) for log score,
  RPS, and RMSE, computed separately for rolling one-step and for each
  recursive horizon `h ∈ {1, 2, 3, 4}`. Run separately for the
  unconstrained fit and the constrained-α fit.
- Across-rep summary: mean, SD, and `2.5% / 97.5%` interval of each
  paired difference.
- Per-rep randomized PIT; stack across reps into a single PIT histogram
  per task and per horizon (§5.6).

**Smoke setting (run before any full launch).** Validate the pipeline
end-to-end with `R = 5`, `n_sims = 100`, `n_per_group = 20`, same seeds
and same parallel infrastructure. Confirms that fit → forward simulator
→ scoring → aggregation produces sensible output and that both the
unconstrained and constrained-α fits run to completion. Only after the
smoke run is clean should the full `R = 1000` configuration be launched.

**Runtime expectation (full run).** Rough order of magnitude is
`R = 1000` reps × (two INAD fits + one NB-GLM fit + `n_sims = 1000`
forward paths × two prediction tasks × ~100 subjects × ≤ 12 steps).
Expect tens of minutes to a few hours wall-clock on commodity hardware
at the `n_workers` setting above. Record measured elapsed time in the
POC output; if a sub-step accidentally goes quadratic in `R`, `n_sims`,
or `N`, it shows up here first.

Before the full launch, run an intermediate staged check with `R = 200`,
`n_sims = 1000`, `n_per_group = 50`, and `fit_max_iter = 50`. This gives
a real runtime estimate for the full run and reduces Monte Carlo/replicate
noise enough to see whether INAD is materially above, at, or below NB GLM.
If INAD clearly beats NB GLM on the headline probabilistic scores, use the
remaining 800 reps as a confirmation run; if it remains indistinguishable,
report parity rather than spending a full-run budget prematurely.

### 5.2 One-step mean (closed form)

All three thinning operators are linear in their input, so
$$\hat\mu_{i,T+1} = \sum_{k=1}^p \alpha_{k,T+1}\, y_{i,T+1-k} + m(\theta_{T+1}, r_{T+1}).$$
This gives a point forecast but not a predictive distribution.

**Integer-valued property: do not round `μ̂`.** Every *realization* of an
INAD process is integer-valued by construction (thinning produces
integers, innovations are integer-valued draws). That integer property
is preserved at the distribution level — every Monte Carlo path returned
by `simulate_inad_forward()` (§5.3) is integer-valued, and the empirical
predictive PMF (§5.4) has integer support. The *conditional mean* `μ̂`
above is real-valued by definition (the expectation of an integer
variable is generally not an integer) and is the right object to compare
against `y_{i,t}` for RMSE/MAE — this is the standard convention for
count-data forecasting (`tscount`, `INGARCH`, classical Poisson
regression). Rounding `μ̂` would introduce bias at small means and would
not change probabilistic scores, which consume the integer-valued PMF
directly. The NB-GLM and Poisson-GLM baselines (§5.5) also report
real-valued conditional means, so the comparison is fair. Integer point
forecasts, if ever needed, are the predictive **mode** or **median**
of the MC sample — both are byproducts of §5.4 and not required for the
POC scoring.

**Predictions vs single realizations (no, scores are not noisy in
`n_sims`).** A single Monte Carlo path drawn from the predictive
distribution is itself a random integer realization and can vary
substantially between draws — by design, because `Y_t` is a random
variable. The *prediction* used for scoring is never a single path:

- Point forecast (RMSE / MAE) = predictive **mean**. For one-step this
  is the closed-form `α·y_{t-1} + m(θ_t, r_t)`. For recursive `h > 1`
  the mean is also analytically iterable, because thinning is linear in
  `Y` (so `E[α ∘ Y] = α · E[Y]`); in practice we read the mean off the
  MC ensemble as `mean(samples)`, which is a Monte Carlo estimate of
  that deterministic quantity with Monte Carlo SE
  `sd(samples) / sqrt(n_sims)`, recorded in the smoke/full output.
- Probabilistic forecast (log score, RPS) = the **entire empirical PMF**
  over all `n_sims = 1000` paths.

Both summaries aggregate over the ensemble, so what enters each score
is stable to MC noise even though individual paths are not. The
intrinsic variability of `Y_t` shows up identically in every competing
model's predictive distribution; proper scoring rules (log score, RPS)
handle that variability correctly by rewarding distributions that put
mass in the right place with the right spread. RMSE alone would
penalize every model for `Y_t`'s irreducible randomness, which is one
of the reasons §5.6 keeps RMSE only as a secondary metric.

### 5.3 Conditional forward simulator (new helper)

The existing `simulate_inad()` always starts at `t = 1` with no conditioning
argument (verified at [simulate_inad.R:59](../../R/simulate_inad.R)). For
prediction, the POC needs a new helper:

```r
simulate_inad_forward(
  fit,                     # an inad_fit
  history,                 # n_test x T_obs matrix; full observed prefix
                           # (helper validates dimensions against fit
                           #  and uses only the last p columns for the
                           #  recursion)
  blocks     = NULL,       # length-n_test integer; required iff fit$tau
                           # is non-trivial (e.g. INADFE designs)
  start_time = ncol(history) + 1L,
  h          = 1L,
  n_sims     = 1000L,
  seed       = NULL
)
# returns: array[n_test, h, n_sims] of simulated forward paths
```

Implementation reuses the thinning and innovation samplers from
`simulate_inad()` but seeds the recursion from the supplied `history`
instead of from `t = 1`, and applies `tau[blocks[i]]` per subject when
the fit has block effects. The helper should error if `fit$tau` is
non-trivial and `blocks` is `NULL`. Live in the POC `common/` directory
until proven; then promote to `R/` as an internal helper used by
`predict.inad_fit()`.

### 5.4 Predictive distribution

- h-step ahead: Monte Carlo via `simulate_inad_forward()`. Build an
  empirical PMF over the `n_sims` paths at each `(subject, t)`.
- Smoothing for log score (finite support, not infinite). For each
  `(subject, t)`, form the empirical PMF on the **union** of the
  simulated values and the observed held-out value; add
  `ε = 1 / (n_sims + 1)` to every cell in that finite support; then
  renormalize. This guarantees the observed value has nonzero predictive
  probability without committing to an infinite support. Report
  sensitivity to `ε` when scores between INAD and the NB-GLM baseline
  are close.
- Smoothing sensitivity: when log-score differences drive a conclusion,
  evaluate at three epsilon values: the main setting
  `1 / (n_sims + 1)`, a smaller value `1 / (10 * n_sims + 1)`, and a
  larger value `1 / (n_sims / 10 + 1)` when `n_sims >= 1000`.
  The current all-pairs bolus output stores scores and diagnostics, not raw
  MC paths or empirical PMF counts, so this cannot be recovered from the
  saved CSV alone. A targeted rerun/rescore must store per-case predictive
  counts or enough MC output to recompute log score and RPS at multiple
  epsilon values.
- RPS: evaluate over the finite support `0:max(simulated, observed)`,
  applying the same smoothing before computing the cumulative PMF. This
  bounds the support explicitly and matches the log-score construction.
- One-step exact PMFs are tractable in special cases (binomial thinning
  with Poisson innovation gives a sum-of-binomials-plus-Poisson, etc.).
  Defer to a later iteration; not required for the POC.

### 5.5 Baselines

Pair innovation families fairly. For the starting NBT-NBI-INADFE(1)
simulation and the subsequent `bolus_inad` analysis, compare against a
negative-binomial GLM with lag, time, and block/treatment terms. A Poisson
GLM can be reported as a lighter secondary baseline, but it should not be
the only comparator for an NB-innovation DGP.

| Baseline | Purpose |
|---|---|
| Last observation carried forward | Floor |
| Time-marginal mean count | Tests lag value |
| NB GLM with lag + time + block/treatment | Family-matched practitioner default |
| Pooled Poisson GLM with lag + time + block/treatment | Lighter secondary comparator |
| NB GLMM with lag + time + block/treatment + subject random intercept | Fairness check for subject-level heterogeneity |
| `tscount::tsglm` with lag and time/block regressors | Exploratory INAR-style comparator |

Exact formula used for both GLM baselines (long-format training data
with one row per `(subject, time)` and `y_lag1` carrying the previous
count; rows at `t = 1` are dropped):

```r
MASS::glm.nb(
  y ~ group + factor(time) + log(y_lag1 + 1),
  data = train_long
)
stats::glm(
  y ~ group + factor(time) + log(y_lag1 + 1),
  family = poisson(link = "log"),
  data = train_long
)
```

The lag covariate is `log(y_lag1 + 1)` rather than raw `y_lag1` to match
the log link. Report the alternative `y ~ group + factor(time) + y_lag1`
specification in a sensitivity row if the conclusion is sensitive to
this choice. For rolling one-step prediction, GLM forecasts use the
realized `y_lag1`. For recursive multi-step prediction, forecasts at
`t > T_split` use the GLM's own iterated-mean recursion with `y_lag1`
replaced by the previous step's predicted mean.

**GLM predictive distribution.** Log score and RPS need a full
predictive distribution, not just a mean.

- NB GLM: mean `μ̂` from `predict(fit_nb, type = "response")`; dispersion
  from `fit_nb$theta` (MLE of the NB size parameter returned by
  `MASS::glm.nb()`). Predictive PMF at row `(subject, t)` is
  `dnbinom(y, size = fit_nb$theta, mu = μ̂)`.
- Poisson GLM: rate `λ̂` from `predict(fit_pois, type = "response")`.
  Predictive PMF is `dpois(y, lambda = λ̂)`.

Evaluate the GLM PMFs on the same finite support used for the INAD PMF
(§5.4), with the same `ε` smoothing, so log scores and RPS are
comparable across models. Smoothing a *known* analytic NB or Poisson
PMF slightly perturbs it, which is an acceptable cost for cross-model
comparability; flag this convention in the report so a reader is not
surprised that the GLM is evaluated under the same finite-support
scoring rule as the empirical INAD PMF.

`tscount::tsglm` (INAR-style) is a reasonable second-tier addition but
not required for the first POC. If it is used, label it exploratory:
`tsglm` is a univariate count-time-series model, not a native panel model.
The script therefore pools the subject-by-time training sequence and does
not reset the internal lag at subject boundaries. This still gives a useful
literature-facing sensitivity check, but it should not be described as a
fully fair panel INAR fit.

The NB GLMM is the cleaner added baseline. Fit it with
`glmmTMB::glmmTMB(y ~ group + factor(time) + log(y_lag1 + 1) +
(1 | subject), family = nbinom2(link = "log"))`. For conditional
new-subject prediction, use population-level predictions (`re.form = NA`,
`allow.new.levels = TRUE`) so the baseline does not leak held-out outcomes
through estimated test-subject random effects.

After the current R=1000 confirmation run completes, rerun the R=200
simulation with the extra baselines before deciding whether the final
writeup needs both of them. Report all three metrics, but keep RPS and RMSE
as the headline if log score remains small, mixed, and smoothing-sensitive.

### 5.6 Metrics

- RMSE / MAE on the predictive mean.
- Log score on the smoothed empirical PMF (§5.4).
- Ranked probability score on the finite support
  `0:max(simulated, observed)`, with the same smoothing.
- Randomized PIT histogram, reported **per horizon** for the recursive
  multi-step task — each of `h = 1, 2, 3, 4` separately, because
  calibration typically degrades with horizon and pooling masks that.
- Predictive-interval coverage at nominal 80% and 95% levels, also
  reported per horizon for recursive prediction.

Report each metric for rolling one-step and for each recursive horizon
separately; do not pool across tasks.

### 5.7 Bolus real-data prediction design

Run the bolus analysis only after the simulated-DGP pipeline is stable. Use
the same two prediction tasks as the simulation:

- rolling one-step: for each held-out patient, predict `t = 2, ..., 12`
  conditional on the patient's realized history through `t - 1`;
- recursive multi-step: reveal times `1:T_split` and predict the remaining
  horizons using the model's own forward recursion. Start with
  `T_split = 8`, hence `h = 1, 2, 3, 4` for times `9:12`.

The `bolus_inad` data contain 65 patients over 12 time points, with 30 in
group 1 and 35 in group 2. Because the sample is small, use a
balanced leave-one-per-group-out design rather than a large test split.
The primary design is the full all-pairs version:

1. In each fold, hold out exactly two patients, one from each group.
2. Fit on the remaining 63 patients.
3. Predict both held-out patients under the rolling one-step and recursive
   multi-step tasks.
4. Repeat over all `30 * 35 = 1050` cross-group held-out pairs.

This design keeps every fold treatment-balanced and uses almost all patients
for fitting, which is more appropriate for `n = 65` than a conventional
70/30 split. The lighter 35-fold rotation is useful only as a smoke check;
it should not be the final bolus analysis.

Because all-pairs folds are dependent (each group-1 patient appears in
35 folds and each group-2 patient appears in 30 folds), fold-level standard
errors are too optimistic. Report all-pairs means, but compute uncertainty
from patient-level clustering or a patient-level bootstrap/permutation
interval. Summaries should be reported overall, by prediction task, by
horizon, and by treatment group.

Primary model choice: constrained-alpha NBT-NBI-INADFE(1), because both the
simulation and the bolus stationarity analysis support time-homogeneous
alpha. Also report the unconstrained fit as a sensitivity check. Primary
baselines are NB GLM and Poisson GLM; NB GLMM is the preferred extra
baseline if `glmmTMB` is available. Treat `tscount::tsglm` as exploratory
for the same panel-data reasons described in section 5.5.

Bolus robustness checks:

- **Training-fold BIC-selected INAD.** In each fold, select constrained-alpha
  or unconstrained INAD using BIC computed from training-fold fits only,
  then use that selected model's already-computed predictions as a new
  `inad_bic_selected` row. This is a post-processing step for the current
  all-pairs output because both INAD fits and their training log-likelihoods
  are saved. Report the constrained/unconstrained selection rate.
- **Smoothing sensitivity for log score.** Target the recursive bolus
  horizons where INAD and NB GLM differ on log score. Since the saved
  all-pairs output does not retain raw MC paths or empirical PMF counts,
  this requires a targeted rerun/rescore that stores enough predictive
  distribution detail to evaluate multiple epsilon values.
- **Larger `n_sims` recursive sensitivity.** Run only if the smoothing
  sweep is inconclusive. Keep it narrow: bolus only, recursive horizons
  only, INAD fits only (especially the BIC-selected or constrained-alpha
  primary), with larger `n_sims` to assess whether MC tail noise drives
  the log-score result. Do not rerun the full baseline suite unless the
  targeted check changes the substantive conclusion.

## 6. Cross-module scoring infrastructure

All three families plug into the same scoring layer because the interface is
*predictive distribution per (subject, time)*:

- Predictive object: list with `family`, plus
  `mean` + `variance` (Gaussian), `pmf` (categorical), or
  `samples` (INAD), keyed by `(subject, time)`.
- Scoring driver: dispatch on `family` to the right rule. Depend on
  `scoringRules` where it implements the rule (Gaussian log-score and
  CRPS for `gau_fit`; sample-based log-score and CRPS for `inad_fit`;
  ranked probability score for ordinal `cat_fit` / discrete `inad_fit`).
  Confirm exact exported function names at implementation time —
  `scoringRules` is not currently a package dependency. Write Brier and
  multiclass log-loss directly.
- Block bootstrap helper over subjects for paired score differences.

## 7. Script-first plan

```
scripts/prediction_poc/
  README.md                  # decisions: design, fallback policies, M, seeds
  common/
    holdout.R                # conditional-new-subject and temporal splitters
    extract_cat.R            # canonical transitions across fit variants
    sim_inad_forward.R       # conditional forward simulator for INAD
    predictive.R             # shared predictive-object constructor
    score.R                  # rmse, mae, logS, CRPS, RPS, Brier, PIT, bootstrap
    plot.R                   # PIT and reliability plots
  gau/
    01_fit.R                 # fit_gau() + pooled-lag regression baseline
    02_predict_1step.R       # closed-form one-step
    02_predict_hstep.R       # state-space recursion + MC validation
    03_score.R
  cat/
    01_fit.R                 # fit_cat() + multinom baseline
    02_predict.R             # exact one-step + h-step recursion + fallback
    03_score.R
  inad/
    00_simulate_nbt_nbi.R    # known NBT-NBI-INADFE(1) DGP; no model selection
    01_fit.R                 # fit known INADFE(1) + NB/Poisson GLM baselines
    02_predict.R             # rolling one-step + recursive h-step MC prediction
    03_score.R               # simulation first, then bolus_inad if simulation passes
  report/
    report.Rmd               # three families side-by-side: scores + PIT plots
```

Datasets:

- Gaussian: `cattle_growth`, `race_100km`.
- Categorical: `labor_force_cat`.
- INAD: start with a simulated NBT-NBI-INADFE(1) DGP. If prediction behaves
  correctly there, repeat the same two prediction tasks on `bolus_inad`.
- Plus a simulated DGP per family using `simulate_*()` so the predictive
  log-score has a reference optimum — analytic DGP entropy for Gaussian
  (`½ log(2πeσ²)`) and categorical (`−Σ p log p`); Monte Carlo
  estimate of the DGP entropy for INAD (large sample from the true
  model, mean log-score under the true predictive distribution).
  Important as a sanity check on the scoring driver itself.

Per family, every experiment follows the same five steps:

1. Split (conditional new-subject by default).
2. Fit AD model(s) and baseline(s). For the INAD POC, fit the known
   NBT-NBI-INADFE(1) structure; do not run order/thinning/innovation
   selection.
3. Produce predictive objects for every held-out `(subject, time)` under
   both rolling one-step and recursive multi-step tasks.
4. Score with the shared driver.
5. Bootstrap paired differences vs the best baseline.

## 8. Promotion criteria

Promote a family to a `predict.*_fit()` method when **all** of these hold
for that family:

1. **Correctness under known DGPs.** Predictive moments match Monte Carlo
   to within Monte Carlo error on simulated data; log-score on a held-out
   sample from the DGP converges to the analytic or Monte Carlo reference
   described in §7 as `n` grows.
2. **Reasonable calibration.** Randomized PIT roughly uniform for
   Gaussian/INAD; reliability diagram near the diagonal for categorical
   on at least one real dataset.
3. **Stable API.** Predictive-object shape settled; transition-extraction
   helper (cat) and forward simulator (INAD) settled.
4. **Acceptable speed.** Initial target: <5 s for 100 subjects ×
   horizon 5 × `M = 1000` (INAD); effectively instant for Gaussian and
   categorical. The INAD number is a target to revisit after prototype
   timings — Bell and NB innovation paths in pure R may need
   vectorization or a compiled inner loop before they hit it. Treat the
   target as a guide, not a strict promotion blocker, until measured.

Real-data performance vs baselines is reported in the report and the
release notes, but is *not* a pass/fail gate. A model is allowed to expose
its model-implied predictions even when a simpler baseline wins on some
datasets.

## 9. Future package API (sketch — do not build yet)

```r
predict(object,
        newdata = NULL,        # NULL: in-sample fitted predictions on training data
                               # non-NULL: conditional forecasts for new subjects,
                               #   supplying their observed history
        h       = 1,
        type    = c("mean", "quantile", "sample", "distribution"),
        n_sims  = 1000,        # INAD only
        level   = 0.95,        # for type = "quantile"
        ...)
```

- `gau_fit`: state-space recursion; `type = "sample"` from the analytic
  Gaussian unless the user requests parameter-uncertainty propagation.
- `cat_fit`: exact recursion; `type = "distribution"` returns the
  predictive PMF over `K` categories with the chosen fallback policy
  exposed as an argument.
- `inad_fit`: Monte Carlo via the (then internal) forward simulator;
  `type = "distribution"` returns the smoothed empirical PMF.

A separate `ad_score()` helper is optional and probably better delegated
to `scoringRules` plus a thin adapter.

## 10. Status

- [x] Stand up `scripts/prediction_poc/` skeleton and `common/` utilities.
- [ ] Gaussian POC: one-step on `cattle_growth` and simulated DGP.
      Smoke scaffold for the simulated DGP is implemented in
      `scripts/prediction_poc/gau/01_smoke.R` and passes.
- [ ] Gaussian POC: multi-step state-space recursion, validated against MC.
      The smoke run validates the order-1 analytic recursion against Monte
      Carlo samples; the full simulation and real-data analysis are still
      pending.
- [ ] Categorical POC: transition extraction helper + exact recursion +
      fallback policy on `labor_force_cat` and simulated DGP.
- [x] INAD POC: simulate from known NBT-NBI-INADFE(1), fit that known
      structure without model selection, and evaluate rolling one-step plus
      recursive multi-step prediction.
      Smoke mode, the staged `R = 200` run, and the combined `R = 1000`
      confirmation run have completed successfully.
- [x] INAD extra-baseline check: after the current R=1000 job finishes,
      install/check optional `tscount` and `glmmTMB`, then rerun `R = 200`
      with NB GLMM and exploratory `tscount::tsglm` comparators.
- [x] Bolus all-pairs prediction: run the same two prediction tasks on
      `bolus_inad` under all `30 * 35 = 1050` leave-one-per-group folds;
      compare with NB GLM, Poisson GLM, NB GLMM, and exploratory
      `tscount::tsglm`; compute patient-level clustered summaries.
- [x] Bolus robustness: add a training-fold BIC-selected INAD row from the
      saved all-pairs output and report constrained/unconstrained selection
      rate. Implemented as a post-processing script; training-fold BIC chose
      constrained-alpha in 1034/1050 folds.
- [ ] Bolus robustness: run targeted smoothing sensitivity for recursive
      log-score cells, storing enough INAD predictive-distribution detail to
      rescore at multiple epsilon values. Smoke test is implemented and
      passing.
- [ ] Optional bolus robustness: if smoothing sensitivity is inconclusive,
      run a narrow larger-`n_sims` recursive INAD-only sensitivity check.
      Smoke test is implemented and passing.
- [ ] `report.Rmd` with the three families side-by-side.
- [ ] Review the §8 promotion criteria; decide per family whether to
      promote, report real-data results honestly either way.
- [ ] If promoted: design doc for `predict.*_fit()` API, then implement.
