# antedep Project Context (2026-02-08)

## Current state

- Package: `antedep` (`DESCRIPTION` version `0.1.0`).
- Main modules implemented: AD (Gaussian), INAD (count), CAT (categorical).
- Missing-data support is integrated for AD (`fit_ad`, `logL_ad`, EM + marginal likelihood helpers).
- INAD/CAT missing-data integration is not complete yet.

## Key findings from audit

- Public API is broad (`43` exports in `NAMESPACE`).
- Tests are substantial (`tests/testthat`: `18` files, `141` test blocks, `432` expectations).
- `R CMD check` log shows warning due missing documentation entries.
- `man/` is currently empty.
- `examples/example_missing_data.R` is stale (sources root-level files that do not exist).
- `tests/testthat/lrt_homogeneity_inad.R` appears to contain implementation code, not tests.

## Primary roadmap

- See `PROJECT_TODO.md` for prioritized TODO list (P0/P1/P2).

## How to resume in a new chat

Use this prompt at the start:

`Read PROJECT_CONTEXT.md and PROJECT_TODO.md in d:\\UI\\MyPackage\\antedep, then continue from the top P0 items.`
