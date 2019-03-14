# sctransform
## R package for normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression

This packaged was developed by Christoph Hafemeister in [Rahul Satija's lab](https://satijalab.org/) at the New York Genome Center. Core functionality of this package has been integrated into [Seurat](https://satijalab.org/seurat/), an R package designed for QC, analysis, and exploration of single cell RNA-seq data.

## Quick start
`devtools::install_github(repo = 'ChristophH/sctransform')`  
`normalized_data <- sctransform::vst(umi_count_matrix)$y`

## Help
For usage examples see vignettes in inst/doc or use the built-in help after installation  
`?sctransform::vst`  

Available vignettes:  
[Variance stabilizing transformation](https://rawgit.com/ChristophH/sctransform/master/inst/doc/variance_stabilizing_transformation.html)  
[Using sctransform in Seurat](https://rawgit.com/ChristophH/sctransform/master/inst/doc/seurat.html)  

## Reference
[Hafemeister, C. & Satija, R. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. bioRxiv 576827 (2019). doi:10.1101/576827](https://www.biorxiv.org/content/10.1101/576827v1)

An early version of this work was used in the paper [Developmental diversification of cortical inhibitory interneurons, Nature 555, 2018](https://github.com/ChristophH/in-lineage).
