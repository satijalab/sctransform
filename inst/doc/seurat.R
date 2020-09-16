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
  fig.width=13, fig.height=6.5, dpi=100, out.width = '95%'
)
library('Seurat')
#old_theme <- theme_set(theme_classic(base_size=8))
future::plan(strategy = 'multicore', workers = 4)
options(future.globals.maxSize = 10 * 1024 ^ 3)

## ----eval=FALSE----------------------------------------------------------
#  devtools::install_github(repo = 'ChristophH/sctransform', ref = 'develop')
#  devtools::install_github(repo = 'satijalab/seurat', ref = 'release/3.0')
#  library(Seurat)
#  library(sctransform)
#  library(ggplot2)

## ----load_data, warning=FALSE, message=FALSE, cache = T------------------
pbmc_data <- Read10X(data.dir = "~/Downloads/pbmc3k_filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc_data)

## ----logtrans, warning=FALSE, message=FALSE, cache = T-------------------
pbmc_logtransform <- pbmc
pbmc_logtransform <- NormalizeData(pbmc_logtransform, verbose = FALSE)
pbmc_logtransform <- FindVariableFeatures(pbmc_logtransform, verbose = FALSE) 
pbmc_logtransform <- ScaleData(pbmc_logtransform, verbose = FALSE) 
pbmc_logtransform <- RunPCA(pbmc_logtransform, verbose = FALSE) 
pbmc_logtransform <- RunUMAP(pbmc_logtransform,dims = 1:20, verbose = FALSE)

## ----apply_sct, warning=FALSE, message=FALSE, cache = T------------------
# Note that this single command replaces NormalizeData, ScaleData, and FindVariableFeatures.
# Transformed data will be available in the SCT assay, which is set as the default after running sctransform
pbmc <- SCTransform(object = pbmc, verbose = FALSE)

## ----pca, cache = T------------------------------------------------------
# These are now standard steps in the Seurat workflow for visualization and clustering
pbmc <- RunPCA(object = pbmc, verbose = FALSE)
pbmc <- RunUMAP(object = pbmc, dims = 1:20, verbose = FALSE)
pbmc <- FindNeighbors(object = pbmc, dims = 1:20, verbose = FALSE)
pbmc <- FindClusters(object = pbmc, verbose = FALSE)

## ----viz, cache = T------------------------------------------------------
pbmc_logtransform$clusterID <- Idents(pbmc)
Idents(pbmc_logtransform) <- 'clusterID'
plot1 <- DimPlot(object = pbmc, label = TRUE) + NoLegend() + ggtitle('sctransform') 
plot2 <- DimPlot(object = pbmc_logtransform, label = TRUE)
plot2 <- plot2 + NoLegend() + ggtitle('Log-normalization') 
CombinePlots(list(plot1,plot2))

## ----fplot, cache = T----------------------------------------------------
VlnPlot(object = pbmc, features = c("CD8A","GZMK","CCL5","S100A4","S100A4","CCR7","ISG15","CD3D"), pt.size = 0.2,ncol = 4)

## ----fplot2, cache = T---------------------------------------------------
FeaturePlot(object = pbmc, features = c("CD8A","GZMK","CCL5","S100A4","S100A4","CCR7"), pt.size = 0.2,ncol = 3)

## ----fplot3, cache = T---------------------------------------------------
FeaturePlot(object = pbmc, features = c("CD3D","ISG15","TCL1A","FCER2","XCL1","FCGR3A"), pt.size = 0.2,ncol = 3)

