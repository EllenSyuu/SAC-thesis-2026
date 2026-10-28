## Submission

First submission of covarchimax (1.0.0).

## Test environments

* Local: Ubuntu 24.04, R 4.3.3
* GitHub Actions: ubuntu-24.04, R 4.3.3

## R CMD check results

`R CMD check --as-cran` (local, with the PDF manual): 0 errors | 0 warnings | 1 note

* "checking for future file timestamps ... unable to verify current time".
  This comes from the build machine having no access to the time server and
  does not concern the package.

On CRAN, the "New submission" note is expected as well.

## Notes for reviewers

* Every exported function has an executable example. All examples together
  run in a few seconds.
* The vignette uses simulated data only.
* References in the Description field carry DOIs.
* The repository also contains a reproducible thesis analysis (`analysis/`,
  `data/`, `results/`). It is excluded from the package tarball via
  `.Rbuildignore`.
