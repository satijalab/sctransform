# sctransform
R package for single cell expression data transformation and normalization

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

