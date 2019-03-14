## ----setup, include = FALSE----------------------------------------------
library('Matrix')
library('ggplot2')
library('reshape2')
library('sctransform')
library('knitr')
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  digits = 2,
  fig.width=8, fig.height=5, dpi=100, out.width = '70%'
)
library('Seurat')
#old_theme <- theme_set(theme_classic(base_size=8))

## ----eval=FALSE----------------------------------------------------------
#  devtools::install_github(repo = 'ChristophH/sctransform', ref = 'develop')
#  devtools::install_github(repo = 'satijalab/seurat', ref = 'release/3.0')
#  library(Seurat)
#  library(sctransform)

## ----load_data, warning=FALSE, message=FALSE, cache = T------------------
pbmc_data <- Read10X(data.dir = "~/Downloads/pbmc3k_filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc_data)

## ----apply_sct, warning=FALSE, message=FALSE, cache = T------------------
# Note that this single command replaces NormalizeData, ScaleData, and FindVariableFeatures.
# Transformed data will be available in the SCT assay, which is set as the default after running sctransform
pbmc <- SCTransform(object = pbmc, verbose = FALSE)

## ----pca, fig.width=5, fig.height=5, cache = T---------------------------
# These are now standard steps in the Seurat workflow for visualization and clustering
pbmc <- RunPCA(object = pbmc, verbose = FALSE)
pbmc <- RunUMAP(object = pbmc, dims = 1:20, verbose = FALSE)
pbmc <- FindNeighbors(object = pbmc, dims = 1:20, verbose = FALSE)
pbmc <- FindClusters(object = pbmc, verbose = FALSE)
DimPlot(object = pbmc, label = TRUE) + NoLegend()

## ----fplot, fig.width = 10, fig.height=10, cache = F---------------------
# These are now standard steps in the Seurat workflow for visualization and clustering
FeaturePlot(object = pbmc, features = c("CD8A","GZMK","CCL5","S100A4"), pt.size = 0.3)
FeaturePlot(object = pbmc, features = c("S100A4","CCR7","CD4","ISG15"), pt.size = 0.3)
FeaturePlot(object = pbmc, features = c("TCL1A","FCER2","XCL1","FCGR3A"), pt.size = 0.3)

