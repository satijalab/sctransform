# News
All notable changes will be documented in this file.

## [0.3] - 2020-09-19
### Added
- Add support for `glmGamPoi` as method to estimate the model parameters; thanks @yuhanH for his pull request
- Add option to use `theta.md`, `theta.mm` or`theta.ml` to estimate theta when `method = 'poisson'`
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
