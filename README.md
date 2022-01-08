# sctransform
## R package for normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression

This package was developed by Christoph Hafemeister in [Rahul Satija's lab](https://satijalab.org/) at the New York Genome Center. Core functionality of this package has been integrated into [Seurat](https://satijalab.org/seurat/), an R package designed for QC, analysis, and exploration of single cell RNA-seq data.

## Quick start

```r
# Install sctransform from CRAN
# install.packages("sctransform")

# Or the development version from GitHub:
# install.packages("remotes")
remotes::install_github("satijalab/sctransform", ref="develop")

normalized_data <- sctransform::vst(umi_count_matrix)$y
```

To invoke the `v2` flavor:

```r
normalized_data <- sctransform::vst(umi_count_matrix, vst.flavor="v2")$y

# Using Seurat
seurat_object <- Seurat::SCTransform(seurat_object, vst.flavor="v2")
```

## Help

For usage examples see vignettes in inst/doc or use the built-in help after installation  
`?sctransform::vst`  

Available vignettes:  

- [Variance stabilizing transformation](https://rawgit.com/satijalab/sctransform/supp_html/supplement/variance_stabilizing_transformation.html)  
- [Using sctransform in Seurat](https://rawgit.com/satijalab/sctransform/supp_html/supplement/seurat.html)  

## Known Issues

* `node stack overflow` error when Rfast package is loaded. The Rfast package does not play nicely with the future.apply package. Try to avoid loading the Rfast package. See discussions: https://github.com/RfastOfficial/Rfast/issues/5 https://github.com/satijalab/sctransform/issues/108

Please use [the issue tracker](https://github.com/satijalab/sctransform/issues) if you encounter a problem

## News
For a detailed change log have a look at the file [NEWS.md](https://github.com/satijalab/sctransform/blob/master/NEWS.md)

### v0.3.2
This release improves the coefficient initialization in quasi poisson regression that sometimes led to errors. There are also some minor bug fixes and a new non-parametric differential expression test for sparse non-negative data (`diff_mean_test`, [this vignette](https://rawgit.com/satijalab/sctransform/supp_html/supplement/np_diff_mean_test.html) gives some details).

### v0.3.1
This release fixes a performance regression when `sctransform::vst` was called via `do.call`, as is the case in the Seurat wrapper. 

Additionally, model fitting is significantly faster now, because we use a fast Rcpp quasi poisson regression implementation (based on `Rfast` package). This applies to methods `poisson`, `qpoisson` and `nb_fast`.

The `qpoisson` method is new and uses the dispersion parameter from the quasi poisson regression directly to estimate `theta` for the NB model. This can speed up the model fitting step considerably, while giving similar results to the other methods. [This vignette](https://rawgit.com/satijalab/sctransform/supp_html/supplement/method_comparison.html) compares the methods.

### v0.3
The latest version of `sctransform` now supports the [glmGamPoi](https://github.com/const-ae/glmGamPoi) package to speed up the model fitting step. You can see more about the different methods supported and how they compare in terms of results and speed [in this new vignette](https://rawgit.com/satijalab/sctransform/supp_html/supplement/method_comparison.html).

Also note that default theta regularization is now based on overdispersion factor (`1 + m / theta` where m is the geometric mean of the observed counts) not `log10(theta)`. The old behavior is still available via `theta_regularization` parameter. You can see how this changes (or doesn't change) the results [in this new vignette](https://rawgit.com/satijalab/sctransform/supp_html/supplement/theta_regularization.html).


## References

- Hafemeister, C. & Satija, R. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. Genome Biol 20, 296 (December 23, 2019).  [https://doi.org/10.1186/s13059-019-1874-1](https://doi.org/10.1186/s13059-019-1874-1). An early version of this work was used in the paper [Developmental diversification of cortical inhibitory interneurons, Nature 555, 2018](https://github.com/ChristophH/in-lineage).

- Choudhary, S. & Satija, R. Comparison and evaluation of statistical error models for scRNA-seq. bioRxiv (2021). [https://doi.org/10.1101/2021.07.07.451498](https://doi.org/10.1101/2021.07.07.451498)

