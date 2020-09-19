# sctransform
## R package for normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression

This package was developed by Christoph Hafemeister in [Rahul Satija's lab](https://satijalab.org/) at the New York Genome Center. Core functionality of this package has been integrated into [Seurat](https://satijalab.org/seurat/), an R package designed for QC, analysis, and exploration of single cell RNA-seq data.

## Quick start
`devtools::install_github(repo = 'ChristophH/sctransform')`  
`normalized_data <- sctransform::vst(umi_count_matrix)$y`

(you can also install from CRAN: `install.packages('sctransform'))`)

## Help
For usage examples see vignettes in inst/doc or use the built-in help after installation  
`?sctransform::vst`  

Available vignettes:  
[Variance stabilizing transformation](https://rawgit.com/ChristophH/sctransform/develop/supplement/variance_stabilizing_transformation.html)  
[Using sctransform in Seurat](https://rawgit.com/ChristophH/sctransform/develop/supplement/seurat.html)  

## News
The latest version of `sctransform` now supports the [glmGamPoi](https://github.com/const-ae/glmGamPoi) package to speed up the model fitting step. You can see more about the different methods supported and how they compare in terms of results and speed [in this new vignette](https://rawgit.com/ChristophH/sctransform/develop/supplement/method_comparison.html).

Also note that default theta regularization is now based on overdispersion factor (`1 + m / theta` where m is the geometric mean of the observed counts) not `log10(theta)`. The old behavior is still available via `theta_regularization` parameter. You can see how this changes (or doesn't change) the results [in this new vignette](https://rawgit.com/ChristophH/sctransform/develop/supplement/theta_regularization.html).

For a detailed change log have a look at the file [NEWS.md](https://github.com/ChristophH/sctransform/blob/develop/NEWS.md)


## Reference
Hafemeister, C. & Satija, R. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. Genome Biol 20, 296 (December 23, 2019). [https://doi.org/10.1186/s13059-019-1874-1](https://doi.org/10.1186/s13059-019-1874-1)

An early version of this work was used in the paper [Developmental diversification of cortical inhibitory interneurons, Nature 555, 2018](https://github.com/ChristophH/in-lineage).
