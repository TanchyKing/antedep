# Local Check Modes

This guide documents the two local check modes used in `antedep` and the expected NOTE messages.

## 1) Tarball mode (recommended before push / release)

One command (document + tests + build + check, with cleanup):

```bash
Rscript scripts/check-package.R --as-cran --no-manual
```

Useful variants:

- Skip roxygen docs refresh:

```bash
Rscript scripts/check-package.R --skip-document --as-cran --no-manual
```

- Keep old artifacts (debugging):

```bash
Rscript scripts/check-package.R --as-cran --no-manual --no-clean
```

## 2) Directory mode (direct `R CMD check` on source tree)

PowerShell (Windows):

```powershell
& "$env:R_HOME\bin\R.exe" CMD check --as-cran --no-manual .
```

General shell:

```bash
R CMD check --as-cran --no-manual .
```

## Expected check outcomes

As of 2026-02-09:

1. Tarball mode (`scripts/check-package.R`) should pass with no ERROR/WARNING.
2. The expected NOTE is:
   - `New submission` (CRAN incoming check for first submission).
3. `unable to verify current time` can appear on some systems/environments.

Directory mode (`R CMD check .`) is expected to be noisier in a normal git working tree.
Typical NOTE/WARNING in directory mode includes:

- source-package warning (`Checking should be performed on sources prepared by 'R CMD build'`)
- hidden/non-standard top-level files (`.git`, `.github`, `docs`, project notes)
- vignette/`inst/doc` warnings
- local check-directory warnings (`..Rcheck`, `antedep.Rcheck`)

Use tarball mode as the release/CRAN decision gate.

## Optional local artifact: PDF reference manual

Generate the integrated function manual from `man/*.Rd`:

```powershell
& "$env:R_HOME\bin\R.exe" CMD Rd2pdf --output=antedep-reference-manual.pdf .
```

This PDF is intentionally local-only (`antedep-reference-manual.pdf` is ignored in git and excluded from package build).
