# News
All notable changes will be documented in this file.

## [0.3.3] - UNRELEASED

### Added
- `vst.flavor` argument to  `vst()` to allow for invoking running updated regularization (sctransform v2, proposed in [Satija and Choudhary, 2021](https://doi.org/10.1101/2021.07.07.451498). See paper for details.
- `scale_factor` to `correct()` to allow for a custom library size when correcting counts


## [0.3.2.9008] - 2021-07-28
### Added
- Add future.seed = TRUE to all `future_lapply()` calls

### Changed
- Wrap MASS::theta.ml() in suppressWarnings()

### Fixed
- Fix logical comparison of vectors of length one in `diff_mean_test()`

## [0.3.2.9003] - 2020-02-11
### Added
- `compare` argument to the nonparametric differential expression test `diff_mean_test()` to allow for multiple comparisons and various ways to specify which groups to compare
- Input checking at various places in `vst()` and `diff_mean_test()`

### Changed
- Major speed improvements for `diff_mean_test()`
- Changed the `labels` argument to `group_labels` in `diff_mean_test()`

### Fixed
- Fix bug where factors in cell attributes gave error when checking for NA, NaN, inf


## [0.3.2] - 2020-12-16
### Added
- Ability to control the values of latent variables when calculating corrected counts
- Offset model as method, including the ability to use a single estimated theta for all genes
- Nonparametric differential expression test for sparse non-negative data

### Changed
- Improve poor coefficient initialization in quasi poisson regression
- When plotting model, do not show density by default; change bandwidth to `bw.nrd0`
- Updates to C++ code to use sparse matrices as S4 objects
- Add check for NA, NaN, Inf values in cell attributes

### Fixed
- Remove biocViews from DESCRIPTION - not needed and was causing problems with deploying shiny apps
- Fix bug where a coefficient was given the wrong name when using `glmGamPoi` (only affected runs with a batch variable set)


## [0.3.1] - 2020-10-08
### Added
- Add a `qpoisson` method for parameter estimation that uses fast Rcpp quasi poisson regression where possible (based on `Rfast` package); this adds `RcppArmadillo` dependency

### Changed
- Remove `poisson_fast` method (replaced by `qpoisson`)
- Use `matrixStats` package and remove `RcppEigen` dependency
- Use quasi poisson regression where possible
- Define cell detection event as counts >= 0.01 (instead of > 0) - this only matters to people playing around with fractional counts (see [issue #65](https://github.com/satijalab/sctransform/issues/65))
- Internal code restructuring and improvements

### Fixed
- Fix inefficiency of using `match.call()` in `vst()` when called via `do.call`

## [0.3] - 2020-09-19
### Added
- Add support for `glmGamPoi` as method to estimate the model parameters; thanks @yuhanH for his pull request
- Add option to use `theta.mm` or`theta.ml` to estimate theta when `method = 'poisson'` or `method = 'nb_fast'`
- Add a `poisson_fast` method for parameter estimation that uses the `speedglm` package and `theta.mm` by default
- Add ability to plot overdispersion factor in `plot_model_pars`
- Add and return time stamps at various steps in the `vst` function
- Add functions to calculate grouped arithmetic and geometric mean per row for sparse matrices (`dgCMatrix`)	- might come in handy some time

### Changed
- Default theta regularization is now based on overdispersion factor (`1 + m / theta` where m is the geometric mean of the observed counts) not `log10(theta)`; old behavior available via `theta_regularization` parameter
- Refactored model fitting code - is now more efficient when using parallel processing
- Changed how message and progress bar output is controlled; integer `verbosity` parameter controls all output: 0 for no output, 1 for only messages, 2 for messages and progress bars
- Increased default bin size (genes being processed simultaneously) from 256 to 500
- Better input checking for cell attributes; more efficient calculation of missing ones

### Fixed
- Some non-regularized model parameters were not plotted

## [0.2.1] - 2019-12-17
### Added
- Add function to generate data given the output of a vst run
- Add cpp support for dense integer matrices
- Minimum variance parameter added to vst function

## [0.2.0] - 2019-04-12
### Added
- Rcpp versions of utility functions
- Helper functions to get corrected UMI and variance of pearson residuals for large UMI matrices

### Changed
- lots of things
