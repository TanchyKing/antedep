# antedep TODO Roadmap (as of 2026-02-08)

## Immediate focus - Reviewer feedback sprint (effort ascending)

- [x] (XS) Unify BIC documentation style across model families (Comment 2)
  - Scope: `bic_cat`, `bic_ad`, `bic_inad` docs.
  - Actions:
    - Add consistent `@details` text to `bic_ad` and `bic_inad` (same formula-style level as `bic_cat`).
    - Standardize `@return` wording to "A numeric scalar BIC value." for all three.
  - Done when:
    - `man/bic_*.Rd` pages show matching `Value` phrasing and parallel `Details` sections.

- [x] (XS) Standardize maximum-likelihood wording in `fit_*` descriptions (Comment 6)
  - Scope: `fit_ad`, `fit_cat`, `fit_inad` docs.
  - Actions:
    - Ensure each `fit_*` help page explicitly states that estimation is likelihood-based / MLE.
    - Keep wording parallel so users can compare model families quickly.
  - Done when:
    - All three `man/fit_*.Rd` descriptions use consistent estimation language.

- [x] (S) Clarify Gaussian-specific scope and naming strategy for AD helpers (Comment 1)
  - Scope: `fit_ad`, `bic_ad`, `bic_order_ad` and related docs.
  - Actions:
    - Explicitly document that `*_ad` currently refers to Gaussian AD functions.
    - Decide non-breaking naming plan: keep existing names and optionally add alias wrappers (`fit_gau`, `bic_gau`, `bic_order_gau`).
  - Done when:
    - No ambiguity remains about overlap with CAT/INAD function families.

- [x] (S) Improve EM interface uniformity/discoverability across families (Comment 5)
  - Scope: `em_inad`, AD EM pathway, docs/examples.
  - Actions:
    - Document clearly that Gaussian EM is available via `fit_ad(..., na_action = "em")`.
    - Decide whether to export `em_ad` wrapper for symmetry (without changing AD internals).
    - Record explicit rationale for CAT (no separate EM helper yet, or planned roadmap).
  - Done when:
    - Users can immediately see where EM is available for each data type.

- [x] (M) Add AIC counterparts for AD and INAD (Comment 3)
  - Scope: `aic_ad`, `aic_inad`, exports, docs, tests.
  - Actions:
    - Implement and export `aic_ad` and `aic_inad` to mirror `aic_cat`.
    - Add unit tests against manual formulas and include examples in docs.
  - Done when:
    - AIC API is symmetric across Gaussian/CAT/INAD model families.

- [x] (L) Add Gaussian confidence-interval helper (`ci_ad`) (Comment 4)
  - Scope: new `ci_ad` API, docs, tests.
  - Actions:
    - Define CI target parameters (e.g., `mu`, `phi`, `sigma`, optional `tau`) and method (Wald/profile/bootstrap where feasible).
    - Implement `ci_ad` with guardrails for unsupported modes and missing-data behavior.
    - Add documentation/examples and tests for standard usage.
  - Done when:
    - CI coverage exists for Gaussian/CAT/INAD families with clearly documented scope.

- [x] (XL) Add additional datasets from book/paper sources (Comment 7)
  - Scope: package data inventory, `data-raw/`, docs, vignettes.
  - Progress (2026-02-11):
    - Added book companion datasets: `cattle_growth` and `cochlear_implant` (source: `stat.uiowa.edu/~dzimmer/Data-for-AD/`).
    - Added categorical cochlear variant `cochlear_implant_cat` for CAT workflows tied to the Xie-Zimmerman application context.
    - Added reproducible raw-data scripts in `data-raw/` and packaged `.rda` objects with documentation.
  - Actions:
    - Add candidate datasets (e.g., cattle, cochlear implant, 2013 Xie-Zimmerman categorical data), subject to redistribution rights.
    - Add dataset documentation (`.Rd`), provenance, preprocessing scripts, and reproducible loaders.
    - Integrate examples into README/vignette so each model family has real-data workflows.
  - Done when:
    - `data/` includes multiple documented datasets beyond `bolus_inad` with clear licensing/source notes.

## P0 - Release blockers

- [x] Clean Git workspace hygiene (post-checkpoint)
  - Issue: commit `65b86dc` included generated build/check artifacts and showed many LF/CRLF warnings.
  - Clean now:
    - Untrack generated artifacts: `antedep.Rcheck/`, `..Rcheck/`, `antedep_0.1.0.tar.gz`, `tests/testthat/Rplots.pdf`.
    - Update `.gitignore` so build/check outputs and transient test plots are not staged.
    - Add `.gitattributes` to normalize line endings and reduce Windows LF/CRLF warning noise.
  - Done when:
    - Running local build/check does not introduce new tracked artifacts.
    - `git status` remains clean after checks (except intentional source/doc changes).
    - Commit/push no longer emits mass LF/CRLF warnings for package text files.

- [x] Add documentation for all exported functions and datasets
  - Scope: all exports in `NAMESPACE` + `bolus_inad`.
  - Done when:
    - `man/` contains `.Rd` entries for all exported objects.
    - `R CMD check` no longer reports "Undocumented code objects" or "Undocumented data sets".

- [x] Run full package checks (not skipped)
  - Scope: examples, tests, and vignettes.
  - Done when:
    - `R CMD check` passes with no ERROR/WARNING on local machine.
    - Any unavoidable NOTE/WARNING is documented in a short decision log.

- [x] Fix stale missing-data example script
  - File: `examples/example_missing_data.R`.
  - Issue: script sources non-existing root files (`missing_utils.R`, `fit_ad_em.R`, etc.).
  - Done when:
    - Script runs from package root using package functions only.
    - No `source("...")` paths that do not exist.

- [x] Clean test layout inconsistency
  - File: `tests/testthat/lrt_homogeneity_inad.R`.
  - Issue: appears to be implementation code, not a `test_that` file.
  - Done when:
    - Non-test code is moved to `R/` (if needed) or file is replaced with real tests.
    - `tests/testthat` only contains test files.

## P1 - Statistical completeness

- [x] Complete AD order-2 missing-data implementation
  - Files: `R/fit_ad_em.R`, `R/logL_ad_missing.R`.
  - Current gap: simplified order-2 updates.
  - Done when:
    - Exact order-2 covariance and M-step are implemented.
    - Added tests prove stable convergence and agreement with simulated truth.

- [x] Add missing-data support for INAD
  - Scope: `fit_inad`, `logL_inad`, related inference paths.
  - Done when:
    - `na_action` strategy is defined and implemented.
    - Missing-data likelihood/EM path is tested on monotone + intermittent patterns.

- [x] Add missing-data support for CAT
  - Scope: `fit_cat`, `logL_cat`, and downstream CAT tests.
  - Done when:
    - Missingness handling is implemented and documented.
    - Unit tests cover MAR-style incomplete categorical sequences.

## P2 - Inference and workflow hardening

- [x] Extend CI/inference utilities for missing-data fits
  - Scope: `ci_*` and `lrt_*` behavior with incomplete data.
  - Done when:
    - Missing-data compatible inference is supported or explicitly rejected with clear errors.

- [x] Add automated check workflow
  - Scope: local script and/or GitHub Actions.
  - Done when:
    - One command runs document + test + check.
    - CI runs on push/PR and publishes results.

- [x] Expand user-facing docs
  - Scope: main README + vignette usage matrix (AD/CAT/INAD, complete/missing data).
  - Done when:
    - New users can identify which methods are production-ready.
    - Known limitations are documented in one place.

## P3 - Workflow polish

- [x] Document local check modes and expected notes/warnings
  - Scope: add a short local-check script/guide for directory-mode vs tarball-mode checks.
  - Done when:
    - A one-command local check is available with cleanup of `*.Rcheck` directories.
    - The team has documented which NOTE/WARNING messages are expected in directory-mode checks.
  - Implemented in: `LOCAL_CHECKS.md`

- [ ] Release memo: regenerate local PDF reference manual when needed
  - Command: `R CMD Rd2pdf --output=antedep-reference-manual.pdf .`
  - Policy: keep local only (ignored by git; excluded from package build).

- [ ] Implement missing-data CI/LRT inference for AD/INAD/CAT
  - Scope: support `ci_*` and `lrt_*` on incomplete data (or model-specific documented approximations).
  - Progress (2026-02-09):
    - Implemented: INAD missing-data LRT wrappers (`lrt_order_inad`, `lrt_homogeneity_inad`, `lrt_stationarity_inad`, and `run_*` helpers) via `na_action = "marginalize"`.
    - Implemented: CAT missing-data order/homogeneity LRT (`lrt_order_cat`, `lrt_homogeneity_cat`, `run_order_tests_cat`) via `na_action = "marginalize"`.
    - Remaining: AD missing-data LRT/mean/covariance tests, CAT missing-data time-invariance/stationarity tests, and missing-data CI (`ci_inad`, `ci_cat`).
  - Done when:
    - Missing-data inference paths are implemented and tested for AD, INAD, and CAT.
    - README/vignette limitations section is updated accordingly.

## Suggested execution order

1. Immediate focus: reviewer feedback sprint (7 comments, XS -> XL).
2. P0 documentation and check cleanup.
3. P0 example/test layout fixes.
4. P1 AD order-2 completion.
5. P1 INAD missing-data support.
6. P1 CAT missing-data support.
7. P2 inference + CI hardening.
