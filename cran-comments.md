## Test environments

- Local: Windows 11, R 4.5.2
- CI: GitHub Actions `ubuntu-latest` (`R CMD check --as-cran --no-manual`)

## R CMD check results

Status: OK — 0 ERRORs, 0 WARNINGs, 0 NOTEs (local Windows 11, R 4.5.2).

## Notes

- This is a patch update (0.1.0 → 0.2.0) adding new S3 methods
  (`deviance`, `confint`, `vcov`, `summary.partial_corr`) and improving
  internal matrix inversion efficiency (`chol2inv` replaces `solve`).
- The package includes no compiled code.
- Reverse dependency checks are not applicable (no known downstream packages).
