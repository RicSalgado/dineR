# Changelog

## dineR 2.0.0

### New features

- [`estimation()`](https://ricsalgado.github.io/dineR/reference/estimation.md)
  now supports parallel execution across the regularisation path via the
  `cores` argument. Setting `cores > 1` distributes each lambda value
  across worker processes using `doSNOW`, yielding 2–3× speed-ups on
  medium-to-large problems (p ≥ 100 or nlambda ≥ 15). Sequential
  execution (`cores = 1`) remains the default.

- [`data_generator()`](https://ricsalgado.github.io/dineR/reference/data_generator.md)
  now accepts asymmetric sample sizes: `n_Y` can be specified
  independently of `n_X`, allowing the two samples to have different
  numbers of observations.

- All differential network matrices returned by
  [`estimation()`](https://ricsalgado.github.io/dineR/reference/estimation.md)
  are now stored as sparse matrices (`dgCMatrix` class via the `Matrix`
  package), reducing memory usage for high-dimensional problems.

### New vignettes

- **Parallelisation** — covers how to switch between sequential and
  parallel modes, documents benchmark results across five problem sizes,
  and provides guidance on choosing the number of cores.

- **Estimation** — step-by-step walkthrough of data generation and the
  estimation workflow.

- **Data Generator** — documents the
  [`data_generator()`](https://ricsalgado.github.io/dineR/reference/data_generator.md)
  function and its outputs in detail.

- **Differential Networks** — end-to-end tutorial on generating data and
  estimating a differential network.

### Bug fixes and improvements

- Fixed `summary.estimation()` S3 method signature to match the
  `summary` generic (`object, ...`), resolving an R CMD check warning.

- Fixed partial argument matching ambiguity in
  [`data_generator()`](https://ricsalgado.github.io/dineR/reference/data_generator.md)
  where `n` matched both `n_X` and `n_Y`.

- Added missing `@importFrom foreach foreach %dopar%` directive,
  resolving undefined global variable notes in R CMD check.

- Added `Matrix` and `foreach` to `Imports` and `doParallel` to
  `Suggests` in `DESCRIPTION`.

## dineR 1.0.1

CRAN release: 2021-11-15

- Initial CRAN release.
