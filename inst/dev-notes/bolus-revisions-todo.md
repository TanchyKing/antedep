# Bolus Comparison Revisions — TODO

Follow-up work from advisor feedback on the bolus predictive comparison and
Table 4.3 in the dissertation. Three discrete issues, each with its own
scripts. Sequencing matters because some refits are shared across issues.

## Issue 1 — Switch the predictive NB GLM baseline from conditional to marginal

### Decision

Use **marginal NB GLM** (`y ~ group + factor(time)`) as the practitioner
baseline in the predictive comparison, matching the dissertation's external BIC
table. Do **not** include conditional NB GLM in this revision pass; that
lag-aware sensitivity can be added later if it becomes necessary.

### Rationale

The earlier conditional NB GLM in `pipeline_inad.R` is a different lag-aware
comparator because `log(y_lag1+1)` lets it use the same history INAD uses. The
marginal NB GLM is the apples-to-apples comparator for the BIC table and is the
only NB GLM carried into this revision pass.

### Tasks

- [x] Add marginal NB GLM fitter to `scripts/prediction_poc/common/pipeline_inad.R`.
- [x] Build prediction helpers that propagate marginal-mean forecasts forward
      (no lag covariate to update).
- [x] Add `21_marginal_glm_baseline_bolus.R` and run the 35-fold rotation
      pass as the fast predictive diagnostic.
- [ ] Rerun the bolus all-pairs predictive comparison with marginal NB GLM as
      the NB GLM comparator.
- [ ] Rerun the R=1000 simulation comparison with marginal NB GLM as primary,
      to keep the baseline choice consistent across analyses.
- [ ] Update `prediction-results-inad.md` headline tables and §3 narrative
      to reflect the new primary baseline.
- [ ] Update paper section narrative so the prediction baseline and the
      dissertation BIC table both refer to marginal NB GLM.

### Files to touch

- `scripts/prediction_poc/common/pipeline_inad.R` — add marginal-GLM fitter
  and predictor.
- New: `scripts/prediction_poc/inad/21_marginal_glm_baseline_bolus.R`.

### Effort

Implementation ~1 day. Bolus all-pairs rerun ~12 h unattended; simulation
rerun ~5–6 h unattended.

---

## Issue 2 — Make tau SE consistent with tau CI (both profile likelihood)

### Decision

Compute tau SE from the profile likelihood, matching the existing
profile-likelihood CI. The two will then satisfy `est ± qnorm(0.975) · SE ≈
(lower, upper)` (exactly if SE is derived from the CI; approximately if SE
is derived from profile-curvature). **Do not extend profile likelihood to
alpha, theta, or nb_inno_size** — their current Wald SE and Wald CI are
mutually consistent, and profile likelihood would add computational cost
without a presentation problem to solve.

### Rationale

The current asymmetry — Wald SE from the Louis observed information at
[ci_inad.R:303](../../R/ci_inad.R), profile CI from `.ci_tau_profile_inad` at
[ci_inad.R:574](../../R/ci_inad.R) — produces the `est ± 2·SE ≠ (lb, ub)`
relationship the advisor flagged in Table 4.3. Profile likelihood is the
more accurate method near the tau boundary, so switching the SE rather than
the CI is the right direction.

### Implementation choice

Two reasonable implementations:

1. **SE derived from CI:** `SE = (upper − lower) / (2 · qnorm(0.975))`.
   Trivial to implement. The relationship `est ± 2·SE = CI` holds by
   construction (modulo asymmetry of the profile CI; report symmetric SE
   anyway and note the CI may be asymmetric).
2. **SE from profile-likelihood curvature:** Evaluate the profile
   log-likelihood at MLE ± small step, fit a parabola, derive SE from the
   curvature. More principled but adds optimization cost.

Recommend option (1) for simplicity unless the profile CI is markedly
asymmetric — in which case symmetric-SE-from-CI mid-width is no longer a
clean summary and option (2) is preferred.

### Tasks

- [x] In `R/ci_inad.R`, modify the tau row construction so that
      `tau$se` is computed from the profile CI (or profile curvature),
      not from Louis information.
- [x] Update unit tests in `tests/testthat/test-ci-inad.R` to reflect
      the new tau SE values.
- [x] Add a NEWS.md entry under the next version describing the change.
- [ ] Regenerate dissertation Table 4.3 with the new SE column.
- [ ] One-sentence footnote on Table 4.3 noting that tau intervals are
      profile-likelihood-based.
- [x] Document in `inst/dev-notes/fit-object-audit.md` that tau alone uses
      profile-likelihood SE/CI while other INAD parameters use Wald.

### Files to touch

- `R/ci_inad.R`
- `tests/testthat/test-ci-inad.R`
- `NEWS.md`
- `inst/dev-notes/fit-object-audit.md`
- New diagnostic: `scripts/prediction_poc/inad/22_tau_ci_se_diagnostic.R`.

### Effort

~Half day implementation + tests; ~30 minutes verification on `bolus_inad`.

---

## Issue 3 — Fitted marginal mean and variance, time by time

### Decision

For each candidate model fit to the full bolus data, compute the fitted
marginal `E[Y_t | group]` and `Var[Y_t | group]` per time per group. The
primary comparison is a Henderson--Shimakura Table-3-style set of time-by-time
mean and variance tables with empirical bolus moments next to each model.
Plots and per-group L1 distances summarize those tables secondarily.

### Rationale

Direct answer to the advisor's BIC-vs-prediction tension. If INAD's fitted
marginal moments track empirical curves much more tightly than NB GLM's do,
that demonstrates a real unconditional-fit advantage that a predictive-score
table alone may not expose. Provides the missing piece between "INAD has much
smaller BIC" and the revised marginal-NB-GLM prediction comparison.

### Models to include

| Model | Marginal moment source |
|---|---|
| Empirical bolus | Sample mean and variance of `y[, t]` per group |
| INAD constrained-α | Analytic recursion |
| INAD unconstrained | Analytic recursion |
| NB GLM marginal | Closed form |
| NB GLMM | Analytic for marginal NB-with-random-intercept (or MC over random effect) |
| tscount (INAR-style) | `simulate.tsglm()` Monte Carlo |
| CGFM independent frailty (Henderson & Shimakura) | Published bolus Table 4 fit + closed-form marginal moments |
| CGFM shared frailty (Henderson & Shimakura) | Published bolus Table 4 fit + closed-form marginal moments |

### INAD marginal-moment recursion (order 1)

Mean:
$$E[Y_{i,t}] = \alpha \cdot E[Y_{i,t-1}] + \theta_t + \tau_{b(i)}, \qquad E[Y_{i,1}] = \theta_1 + \tau_{b(i)}.$$

Variance (NB thinning + NB innovation):
$$\text{Var}[Y_{i,t}] = \alpha(1+\alpha)\,E[Y_{i,t-1}] + \alpha^2 \text{Var}[Y_{i,t-1}] + \theta_t + \theta_t^2 / r_t.$$

(Confirm the thinning-variance term against the package's NB-thinning code
before use; the leading coefficient depends on the parameterization.)

### CGFM moment-formula extraction

The repo and the installed comparison packages do not currently expose a
Henderson--Shimakura CGFM fitter. For this revision pass, use the published
bolus Table 4 estimates for the independent- and shared-frailty CGFM rows.
Their marginal means and variances are read directly from the time-specific
log means, group coefficient, and gamma-frailty variance in the paper. Keep
those rows labeled as published fits rather than local refits.

### Tasks

- [x] Refit locally available candidate models on full bolus data (no
      cross-validation needed):
  - INAD constrained-α and unconstrained.
  - Marginal NB GLM.
  - NB GLMM.
  - tscount.
  - Use published Henderson--Shimakura Table 4 parameters for the CGFM
    independent- and shared-frailty rows.
- [x] Implement `marginal_moments(fit, group_levels, n_time)` returning
      mean and variance vectors per group per time, per model.
- [x] Henderson--Shimakura Table-3-style empirical mean/variance table by
      time and per-group fitted mean/variance comparison tables.
- [x] Companion signed relative-discrepancy tables,
      `100 * (fitted - empirical) / empirical`, for the four
      group-by-moment combinations.
- [x] Per-group figure: empirical mean curve + each model's fitted mean curve
      overlaid. Same for variance.
- [x] L1-distance table: per (model, group), `sum_t |fitted_mean_t -
      empirical_mean_t|` and the variance analogue.
- [ ] Add figures and table to the bolus section of the paper, framed as
      "Marginal-moment fit on bolus."
- [ ] One-paragraph interpretation in the paper: if INAD wins this comparison,
      it bridges the BIC gap to the predictive tie. If the comparison is also
      a tie, it's still informative and the BIC gap moves to "explained by
      likelihood-density evaluation rather than marginal-moment fit."

### Files to touch

- New: `scripts/prediction_poc/inad/23_marginal_moments_bolus.R` — driver.
- New: `scripts/prediction_poc/common/marginal_moments.R` — per-model
  marginal-moment computation.
- New: `scripts/prediction_poc/common/cgfm_marginal.R` — CGFM moments
  per Henderson & Shimakura.
- New figures: `scripts/prediction_poc/output/inad/marginal_moments/`.

### Effort

~3 days. CGFM implementation is the main work; INAD recursion is short;
NB GLM marginal is one line; NB GLMM and tscount are MC wrappers.

---

## Sequencing

1. **Issue 2 first** (~half day, lowest risk). Resolves Table 4.3
   immediately and is independent of the others.
2. **Issue 1 baseline rerun and Issue 3 model refits in parallel.** The
   full-data fits needed for marginal moments overlap with Issue 1's
   refitting, so do both in one scripting session.
3. **Bolus all-pairs rerun (Issue 1, ~12 h)** can run overnight while
   Issue 3's marginal-moment figures are produced from the same fits.
4. **Simulation rerun (Issue 1, ~6 h)** runs after the bolus is finished.
5. **Update paper sections** with revised numbers, the new marginal-moment
   figure, and the consistent baseline choice.

Estimated total: ~5 days of script work + ~18 h unattended compute. The
critical-path output is the marginal-moment figure (Issue 3).

---

## Open methodological questions

- **Simulation baseline consistency.** Switching the simulation comparison
  to marginal NB GLM keeps the baseline choice consistent across analyses,
  but the simulation DGP has time-varying innovation parameters that the
  marginal NB GLM cannot capture. The writeup should state that the comparator
  is *consistently* marginal NB GLM because the revision aligns prediction with
  the dissertation BIC comparison.
- **Profile likelihood for non-tau parameters.** Deferred. Current
  Wald SE and Wald CI for alpha/theta/nb_inno_size are mutually consistent
  and computationally cheap. Revisit only if a user-visible inconsistency
  surfaces.
- **Marginal-moment comparison: paper figure or supplementary?** Recommend
  headline figure in the bolus section. It directly resolves the
  BIC-vs-prediction tension and is more informative than a long footnote.
- **CGFM Type I vs Type II naming.** Cross-check the names used in the
  dissertation Table 4.x vs the labels in Henderson & Shimakura before
  publishing so readers can map them.
