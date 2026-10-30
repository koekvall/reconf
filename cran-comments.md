## Submission notes

This update follows the first release (2026-10-21) closely because the
RcppEigen maintainers asked for it: RcppEigen will move to Eigen 5, which
removes `Eigen::MappedSparseMatrix`, and reconf 0.2.0 would then fail to
compile (https://github.com/RcppCore/RcppEigen/issues/151). Version 0.2.1
uses `Eigen::Map<Eigen::SparseMatrix<double>>` instead, which compiles with
both the current RcppEigen and its release candidate. The update also adds
warnings when an internal maximization stops away from a root of the score.

## Test environments

- Local: macOS 27 (Apple Silicon), R 4.6.1, with RcppEigen 0.3.4.0.2 from
  CRAN and with the release candidate 0.4.9.9.5 from
  https://rcppcore.r-universe.dev

## R CMD check results

0 errors | 0 warnings | 0 notes

## Downstream dependencies

None.
