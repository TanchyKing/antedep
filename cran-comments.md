## Test environments

- Local: Windows 11, R 4.5.2 (`R CMD check --as-cran`)

## R CMD check results

0 errors | 1 warning | 0 notes

- The warning, "'qpdf' is needed for checks on size reduction of PDFs", is
  local: qpdf is not installed on the test machine.

## Notes

- This is an update from 0.2.0 to 0.3.1. Version 0.3.0 was not released on
  CRAN; 0.3.1 differs from it only in examples and documentation.
- Main changes (see NEWS.md): an exact EM algorithm for all nine
  thinning-innovation combinations; observed information computed exactly by
  Louis' identity for confidence intervals; time-invariance constrained fits;
  a new `kkt_inad()` for checking first-order conditions; and corrections to
  direct maximization, the missing-data likelihood, the chi-bar-square
  p-value, and the parameter counts in test tables.
- `test_order_inad()` no longer reports an ordinary chi-square p-value for
  comparisons between orders, which are boundary problems.
- The package includes no compiled code.
- There are no reverse dependencies.
