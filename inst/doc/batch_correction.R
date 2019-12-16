## ----setup, include = FALSE----------------------------------------------
library('Matrix')
library('ggplot2')
library('reshape2')
library('sctransform')
library('knitr')
knit_hooks$set(optipng = hook_optipng)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  digits = 2,
  optipng = '-o7',
  fig.width=5, fig.height=3, dpi=100, out.width = '70%'
)
old_theme <- theme_set(theme_classic(base_size=8))

## ---- warning=FALSE------------------------------------------------------
# some of the vst steps can use multiple cores
# We use the Future API for parallel processing; set parameters here
future::plan(strategy = 'multicore', workers = 4)
options(future.globals.maxSize = 10 * 1024 ^ 3)

# load data and cluster identities, and calculate some cell attributes
cm <- readRDS(file = '~/Projects/data/BipolarCell2016_GSE81904_sparse_matrix.Rds')
cell_attr <- read.table(file = '~/Projects/data/BipolarCell2016_GSE81904_ClustAssignFile.txt', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
cell_attr <- cell_attr[colnames(cm), ]
cell_attr$CLUSTER <- factor(cell_attr$CLUSTER)
cell_attr$replicate <- as.integer(sapply(strsplit(rownames(cell_attr), '_'), function(x) gsub('Bipolar', '', x[1])))
cell_attr$batch <- factor(as.integer(cell_attr$replicate %in% c(5, 6)) + 1)
cell_attr$gene <- Matrix::colSums(cm > 0)
cell_attr$umi <- as.integer(Matrix::colSums(cm))
cell_attr$log_umi <- log10(cell_attr$umi)
cell_attr$umi_per_gene <- cell_attr$umi / cell_attr$gene
cell_attr$log_umi_per_gene <- log10(cell_attr$umi_per_gene)

# remove 50 weird 'cells' 
weird <- cell_attr$log_umi_per_gene < 0.05 | cell_attr$log_umi_per_gene > 0.35
cm <- cm[, !weird]
cell_attr <- cell_attr[!weird, ]

# remove genes as in the paper; goes from 24904 down to 13181
bad_genes <- Matrix::rowSums(cm > 0) < 30 | Matrix::rowSums(cm) < 60
cm <- cm[!bad_genes, ]

# use fewer cells to speed up vignette compilation
set.seed(42)
keep <- sample(1:ncol(cm), 10000)
cm <- cm[, keep]
cell_attr <- cell_attr[keep, ]

# apply vst
set.seed(42)
vst_out <- vst(cm, cell_attr = cell_attr, latent_var = 'log_umi_per_gene', bin_size = 128, show_progress = FALSE)

## ---- fig.width=4.5, fig.height=5.5, out.width='49%', fig.show='hold'----
pca <- irlba::prcomp_irlba(t(vst_out$y), n = 35)
tsne <- Rtsne::Rtsne(scale(pca$x), pca = FALSE)
cell_attr$tSNE1 <- tsne$Y[, 1]
cell_attr$tSNE2 <- tsne$Y[, 2]
ggplot(cell_attr, aes(tSNE1, tSNE2, color = batch)) + geom_point(alpha=0.5, shape=16) + 
  theme(legend.position="bottom")
ggplot(cell_attr, aes(tSNE1, tSNE2, color = factor(CLUSTER))) + geom_point(alpha=0.5, shape=16) +
  theme(legend.position="bottom")

## ---- fig.width=4.5, fig.height=5.5, out.width='49%', fig.show='hold', warning=FALSE----
set.seed(42)
vst_out2 <- vst(cm, cell_attr = cell_attr, latent_var = 'log_umi_per_gene', batch_var = 'batch', bin_size = 128, show_progress = FALSE)

# dimensionality reduction via PCA
pca2 <- irlba::prcomp_irlba(t(vst_out2$y), n = 35)

tsne2 <- Rtsne::Rtsne(scale(pca2$x), pca = FALSE)
cell_attr$tSNE1 <- tsne2$Y[, 1]
cell_attr$tSNE2 <- tsne2$Y[, 2]
ggplot(cell_attr, aes(tSNE1, tSNE2, color = batch)) + geom_point(alpha=0.5, shape=16) + 
  theme(legend.position="bottom")
ggplot(cell_attr, aes(tSNE1, tSNE2, color = factor(CLUSTER))) + geom_point(alpha=0.5, shape=16) +
  theme(legend.position="bottom")

## ---- fig.width=4.5, fig.height=7, out.width='49%', fig.show='hold'------
# plot expression and model for one gene as example
cell_attr$idx <- NA
cell_attr$idx[order(cell_attr$batch, -cell_attr$log_umi_per_gene)] <- 1:ncol(cm)
goi <- 'mt-Rnr2'
plot_model(vst_out, cm, goi, x_var = 'idx', cell_attr = cell_attr,
           batches = cell_attr$batch, plot_residual = TRUE, show_density = FALSE)
plot_model(vst_out2, cm, goi, x_var = 'idx', cell_attr = cell_attr,
           batches = cell_attr$batch, plot_residual = TRUE, show_density = FALSE)

