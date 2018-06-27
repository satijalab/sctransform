# sctransform
R package for modeling single cell UMI expression data using regularized negative binomial regression

This packaged was developed by Christoph Hafemeister in [Rahul Satija's lab](https://satijalab.org/) at the New York Genome Center. A previous version of this work was used in the paper [Developmental diversification of cortical inhibitory interneurons, Nature 555, 2018](https://github.com/ChristophH/in-lineage) and is, in part, also implemented in [Seurat](https://satijalab.org/seurat/), an R package designed for QC, analysis, and exploration of single cell RNA-seq data.

This package is in beta status, please sanity check any results, and notify me of any issues you find.

## Quick start
`devtools::install_github(repo = 'ChristophH/sctransform')`  
`normalized_data <- sctransform::vst(umi_count_matrix)$y`

## Help
For usage examples see vignettes in inst/doc or use the built-in help after installation  
`?sctransform::vst`  
`browseVignettes(package = 'sctransform')`

Available vignettes:  
[Variance stabilizing transformation](https://rawgit.com/ChristophH/sctransform/master/inst/doc/variance_stabilizing_transformation.html)  
[Differential expression](https://rawgit.com/ChristophH/sctransform/master/inst/doc/differential_expression.html)  
[Batch correction](https://rawgit.com/ChristophH/sctransform/master/inst/doc/batch_correction.html)  
[Denoising](https://rawgit.com/ChristophH/sctransform/master/inst/doc/denoising.html)  

