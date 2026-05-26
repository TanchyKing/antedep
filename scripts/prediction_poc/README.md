# Prediction POC Scripts

This directory contains script-first prediction experiments. The INAD track is
implemented first; Gaussian and categorical tracks are placeholders for later
work.

## INAD Smoke Run

From `D:/UI/MyPackage`, run:

```sh
Rscript antedep/scripts/prediction_poc/inad/03_score.R smoke
```

From the package root (`D:/UI/MyPackage/antedep`), run:

```sh
Rscript scripts/prediction_poc/inad/03_score.R smoke
```

Smoke settings are `R = 5`, `n_sims = 100`, and `n_per_group = 20`. Outputs are
written to:

```text
scripts/prediction_poc/output/inad/smoke/
```

## INAD Full Run

After the smoke output is clean:

```sh
Rscript scripts/prediction_poc/inad/06_intermediate_r200.R
```

The intermediate run uses `R = 200`, `n_sims = 1000`, `n_per_group = 50`,
and `fit_max_iter = 50`. Outputs are written to:

```text
scripts/prediction_poc/output/inad/intermediate_r200/
```

Use the intermediate output to decide whether the remaining 800 confirmation
replicates are warranted.

To run the remaining reps `201:1000` and combine them with the existing
`R = 200` output:

```sh
Rscript scripts/prediction_poc/inad/07_remaining_r800.R
```

This runner saves 20-rep chunks under:

```text
scripts/prediction_poc/output/inad/remaining_r800/chunks/
```

If interrupted, rerun the same command; completed chunks are skipped. After all
chunks finish, combined `R = 1000` outputs are written to:

```text
scripts/prediction_poc/output/inad/final_r1000/
```

## Extra INAD Baselines

The optional extra-baseline runner adds:

- `glmmTMB` negative-binomial GLMM with a subject random intercept.
- `tscount::tsglm` as an exploratory INAR-style comparator.

These packages are not package dependencies. After any long-running INAD job
has finished, install/check them with:

```sh
Rscript scripts/prediction_poc/inad/08_install_extra_baselines.R
```

Then rerun the `R = 200` experiment with all baselines scored on the same INAD
finite support:

```sh
Rscript scripts/prediction_poc/inad/09_extra_baselines_r200.R
```

Outputs are written to:

```text
scripts/prediction_poc/output/inad/extra_baselines_r200/
```

The `tscount` baseline is marked exploratory because `tsglm` is a univariate
time-series model, not a native panel model. The script pools the subject-by-time
training sequence, so subject-boundary lags are not reset inside `tscount`.

For the full scripted run:

```sh
Rscript scripts/prediction_poc/inad/03_score.R full
```

Full settings are `R = 1000`, `n_sims = 1000`, and `n_per_group = 50`. Worker
count is `max(1, floor(2 * parallel::detectCores(logical = TRUE) / 3))`.

## Files

- `common/config.R`: shared configuration, package loading, output paths.
- `common/holdout.R`: stratified subject splits and long lag-1 data.
- `common/sim_inad_forward.R`: conditional INAD forward simulator.
- `common/constrained_inad.R`: POC constrained-alpha fit adapter.
- `common/score.R`: finite-support log score, RPS, PIT, coverage helpers.
- `common/pipeline_inad.R`: INAD simulation, fitting, prediction, scoring.
- `common/extra_baselines.R`: optional `glmmTMB` and `tscount` baselines.
- `inad/00_simulate_nbt_nbi.R`: single known-DGP simulation.
- `inad/01_fit.R`: one-rep fitting check.
- `inad/02_predict.R`: one-rep prediction/scoring check.
- `inad/03_score.R`: smoke/full replicated runner.
- `inad/04_smoke_iter50.R`: smoke diagnostic with full iteration budget.
- `inad/05_smoke_n50_iter50.R`: sample-size diagnostic at full planned `n`.
- `inad/06_intermediate_r200.R`: staged `R = 200` intermediate run.
- `inad/07_remaining_r800.R`: chunked reps `201:1000` plus final combine.
- `inad/08_install_extra_baselines.R`: optional baseline package installer.
- `inad/09_extra_baselines_r200.R`: R=200 rerun with extra baselines.
