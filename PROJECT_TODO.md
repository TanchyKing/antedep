# antedep TODO Roadmap (as of 2026-02-08)

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

- [ ] Add automated check workflow
  - Scope: local script and/or GitHub Actions.
  - Done when:
    - One command runs document + test + check.
    - CI runs on push/PR and publishes results.

- [ ] Expand user-facing docs
  - Scope: main README + vignette usage matrix (AD/CAT/INAD, complete/missing data).
  - Done when:
    - New users can identify which methods are production-ready.
    - Known limitations are documented in one place.

## P3 - Workflow polish

- [ ] Document local check modes and expected notes/warnings
  - Scope: add a short local-check script/guide for directory-mode vs tarball-mode checks.
  - Done when:
    - A one-command local check is available with cleanup of `*.Rcheck` directories.
    - The team has documented which NOTE/WARNING messages are expected in directory-mode checks.

## Suggested execution order

1. P0 documentation and check cleanup.
2. P0 example/test layout fixes.
3. P1 AD order-2 completion.
4. P1 INAD missing-data support.
5. P1 CAT missing-data support.
6. P2 inference + CI hardening.
