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
  fig.width=4, fig.height=2.5, dpi=100, out.width = '70%'
)
old_theme <- theme_set(theme_classic(base_size=10))

## ----load_data-----------------------------------------------------------
pbmc_clusters <- readRDS(file = "~/Projects/data/pbmc3k_celltypes.rds")
pbmc_data <- readRDS(file = "~/Projects/data/pbmc3k_umi_counts.Rds")
pbmc_data <- pbmc_data[, names(pbmc_clusters)]

class(pbmc_data)
dim(pbmc_data)

## ---- fig.width=4, fig.height=2.5, warning=FALSE-------------------------
# some of the vst steps can use multiple cores
# We use the Future API for parallel processing; set parameters here
future::plan(strategy = 'multicore', workers = 4)
options(future.globals.maxSize = 10 * 1024 ^ 3)

set.seed(43)
vst_out <- sctransform::vst(pbmc_data, latent_var = c('log_umi_per_gene'), return_gene_attr = TRUE, return_cell_attr = TRUE, show_progress = FALSE)

## ------------------------------------------------------------------------
res1 <- sctransform:::compare_expression(x = vst_out, umi = pbmc_data, group = pbmc_clusters, 
                                        val1 = 'CD14+ Monocytes', 
                                        val2 = 'FCGR3A+ Monocytes', 
                                        show_progress = FALSE)

## ------------------------------------------------------------------------
head(subset(res1, log_fc < 0), 10)
head(subset(res1, log_fc > 0), 10)

## ---- fig.height=3.5-----------------------------------------------------
ggplot(res1, aes(log_fc, -log10(p_value))) + geom_point(alpha=0.4, shape=16)

## ------------------------------------------------------------------------
res2 <- sctransform:::compare_expression(x = vst_out, umi = pbmc_data, group = pbmc_clusters, val1 = setdiff(pbmc_clusters, 'CD8 T cells'), val2 = 'CD8 T cells', show_progress = FALSE)
head(subset(res2, log_fc > 0), 10)

## ---- fig.width=4, fig.height=4, out.width='49%', fig.show='hold'--------
goi <- rownames(subset(res2, log_fc > 0))[1:3]
df <- melt(t(as.matrix(pbmc_data[goi, ])), varnames = c('cell', 'gene'), value.name = 'UMI')
df$cluster <- pbmc_clusters
ggplot(df, aes(x = gene, y = log10(UMI + 1), color = cluster == 'CD8 T cells')) + geom_violin(scale = 'width')
df <- melt(t(as.matrix(vst_out$y[goi, ])), varnames = c('cell', 'gene'), value.name = 'Pearson_residual')
df$cluster <- pbmc_clusters
ggplot(df, aes(x = gene, y = Pearson_residual, color = cluster == 'CD8 T cells')) + geom_violin(scale = 'width')


