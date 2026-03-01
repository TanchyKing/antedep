## Test environments

- Local: Windows 11, R 4.5.2
- CI: GitHub Actions `ubuntu-latest` (`R CMD check --as-cran --no-manual`)

## R CMD check results

- Local targeted checks completed:
  - `testthat::test_dir("tests/testthat", load_package = "source")`
  - `pkgdown::build_site_github_pages(new_process = FALSE, install = FALSE)`
- CI full package check runs on push to `master` via `.github/workflows/r-check.yml`.

## Notes

- This is a first CRAN submission for `antedep` (version 0.1.0).
- The package includes no compiled code.
- Reverse dependency checks are not applicable (no known downstream packages).

