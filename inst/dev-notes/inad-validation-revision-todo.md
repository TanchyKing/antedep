# INAD Validation Summary - Revision TODO

Status: complete as of 2026-05-23.

This file tracks the v2 rewrite of `inad-validation-summary.md`. The v1
summary is preserved at `inad-validation-summary-v1.md`; the revised draft
was generated at `inad-validation-summary-v2.md`, reviewed during this pass,
and then copied to the canonical `inad-validation-summary.md`. The HTML render
was regenerated at `inad-validation-summary.html`.

## Direction Decisions

- [x] Simulation reports INAD-only prediction precision. Because the
      simulation DGP is INAD, the simulation section is treated as
      correct-specification validation rather than a neutral cross-model
      benchmark.
- [x] CGFM is dropped from the predictive comparison. It remains in the
      marginal-moment comparison only, using published Henderson-Shimakura
      estimates.
- [x] Model formulations are stated explicitly before the revised results.
- [x] Bolus prediction is restricted to rolling one-step prediction in the
      headline section; recursive results are retained only as smoothing
      sensitivity context.

## Phase 1 - Supporting Script Work

- [x] Added `scripts/prediction_poc/inad/24_validation_revision_support.R`.
- [x] Computed full-bolus BIC comparison outputs:
      `scripts/prediction_poc/output/inad/validation_revision/bolus_bic_comparison.csv`.
- [x] Reconstructed CGFM BIC rows from Henderson-Shimakura Table 4
      log-likelihoods where full log-likelihood was available.
- [x] Extracted INAD-only rolling one-step predictive precision from the
      existing R = 1000 simulation:
      `simulation_inad_precision_r1000.csv`.
- [x] Kept the existing recursive smoothing sensitivity as methodological
      context rather than rerunning rolling smoothing, since the revised
      headline no longer rests on recursive log-score behavior.
- [x] Wrote run provenance to
      `run_info_validation_revision.csv`.

## Phase 2 - Document Rewrite

- [x] Added early CGFM scope statement.
- [x] Added metric primer for RMSE, log score, and RPS.
- [x] Added model formulation section.
- [x] Rewrote simulation section to INAD-only precision.
- [x] Added bolus BIC comparison table.
- [x] Rewrote bolus prediction section to rolling one-step only.
- [x] Dropped BIC-selected INAD subsection from the headline document.
- [x] Expanded smoothing sensitivity context.
- [x] Retained marginal-moment comparison with CGFM scope clarified.
- [x] Kept plug-in prediction, integer-valued prediction, tau SE/CI, and
      engineering-validation notes.
- [x] Updated package implications to landed-vs-remaining structure.
- [x] Updated internal cross-references under the new section order.
- [x] Kept the gau/cat bridge as forward-looking next work.

## Phase 3 - Review Checkpoints

- [x] Confirmed model formulations sufficiently for the dev-note revision.
- [x] Confirmed section ordering by implementing it in v2.
- [x] Confirmed CGFM scope wording.
- [x] Confirmed "INAD predictive precision under DGP" content:
      absolute log score, RPS, RMSE, PIT, and coverage summaries.
- [x] Drafted `inad-validation-summary-v2.md`.
- [x] Reviewed generated v2 for table/output consistency.
- [x] Replaced `inad-validation-summary.md` with the v2 content.
- [x] Regenerated `inad-validation-summary.html`.

## Remaining Non-Blocking Items

These are intentionally left as future work rather than TODO blockers for the
v2 summary:

- Add a simulated 3+ group unit test for tau SE-from-CI behavior.
- Decide the shared cross-family `type = "distribution"` prediction API after
  the gau/cat prediction POCs.
- Implement a proper joint `alpha_constraint = "constant"` path in `fit_inad()`
  only if the constrained-alpha story remains central after the cross-family
  API work.
- Treat exact one-step PMFs as future methodological cleanup, not a blocker for
  the current validation summary.
- **Misspecification simulation studies for §3.** The current §3 simulation
  uses INAD as both DGP and model, which is correct-specification validation
  only. A more complete picture would simulate from each baseline DGP — NB
  GLM, NB GLMM, Poisson GLM, and CGFM (independent and shared frailty) — and
  ask how INAD predictions hold up against the true model under
  misspecification. Substantial implementation: requires DGP samplers for
  each baseline, a refit-and-score driver, and honest interpretation that
  does not pretend INAD's structural advantage carries over universally.
  Out of scope for the v2 summary; defer until cross-DGP robustness becomes
  a reviewer-driven requirement.
