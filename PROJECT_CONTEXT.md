# antedep Project Context (updated 2026-02-09)

## Current state

- Package: `antedep` (`DESCRIPTION` version `0.1.0`).
- Main modules implemented: AD (Gaussian), INAD (count), CAT (categorical).
- Missing-data support is integrated for AD (`fit_ad`, `logL_ad`, EM + marginal likelihood helpers),
  INAD, and CAT (`fit_cat`, `logL_cat` with `na_action`).

## Key findings from audit

- Public API is broad (`43` exports in `NAMESPACE`).
- Tests are substantial (`tests/testthat`: `18` files, `141` test blocks, `432` expectations).
- Prior P0 documentation/check issues were resolved in subsequent commits.
- Current focus has shifted to P2 inference/workflow hardening tasks.

## Git hygiene issue (2026-02-08)

- Checkpoint commit `65b86dc` successfully pushed, but included generated artifacts.
- Push/commit showed mass LF/CRLF conversion warnings on Windows (`core.autocrlf=true`).
- Cleanup is needed before more feature work:
  - Untrack generated artifacts: `antedep.Rcheck/`, `..Rcheck/`, `antedep_0.1.0.tar.gz`, `tests/testthat/Rplots.pdf`.
  - Add/adjust ignore rules so check/build outputs are never staged again.
  - Add a `.gitattributes` policy so text files stay normalized and warning noise is reduced.

## Primary roadmap

- See `PROJECT_TODO.md` for prioritized TODO list (P0/P1/P2).

## How to resume in a new chat

Use this prompt at the start:

`Read PROJECT_CONTEXT.md and PROJECT_TODO.md in d:\\UI\\MyPackage\\antedep, then continue from the top P0 items.`
