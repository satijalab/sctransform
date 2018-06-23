## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  digits = 2,
  fig.width=4, fig.height=2.5, dpi=300, out.width = '70%'
)

## ----load_packages, include=FALSE----------------------------------------
library('Matrix')
library('ggplot2')
library('reshape2')
library('sctransform')

## ----load_data-----------------------------------------------------------
pbmc_data <- readRDS(file = "~/Projects/data/pbmc3k_umi_counts.Rds")
class(pbmc_data)
dim(pbmc_data)

## ---- fig.width=4, fig.height=2.5----------------------------------------
options(mc.cores = 7)
set.seed(43)
vst_out <- sctransform::vst(pbmc_data, latent_var = c('log_umi_per_gene'), return_gene_attr = TRUE, return_cell_attr = TRUE, show_progress = FALSE)

## ------------------------------------------------------------------------
var_genes <- rownames(vst_out$gene_attr)[scale(vst_out$gene_attr$residual_variance)[, 1] > 1]
pca <- irlba::prcomp_irlba(t(vst_out$y[var_genes, ]), n = 20)
ce <- pca$x[, 1:10]
hcl <- hclust(dist(ce), method='ward.D2')
clustering <- cutree(hcl, 10)

## ------------------------------------------------------------------------
res1 <- sctransform::compare_expression(x = vst_out, umi = pbmc_data, group = clustering, val1 = 1, val2 = 2, show_progress = FALSE)

## ------------------------------------------------------------------------
head(subset(res1, log_fc < 0), 10)
head(subset(res1, log_fc > 0), 10)

## ---- fig.height=3.5-----------------------------------------------------
ggplot(res1, aes(log_fc, -log10(p_value))) + geom_point(alpha=0.4, shape=16)

## ------------------------------------------------------------------------
res2 <- sctransform::compare_expression(x = vst_out, umi = pbmc_data, group = clustering, val1 = setdiff(clustering, 3), val2 = 3, show_progress = FALSE)

## ---- fig.width=4, fig.height=4, out.width='49%', fig.show='hold'--------
goi <- rownames(subset(res2, log_fc > 0))[1:3]
df <- melt(t(as.matrix(pbmc_data[goi, ])), varnames = c('cell', 'gene'), value.name = 'UMI')
df$cluster <- clustering
ggplot(df, aes(x = gene, y = log10(UMI + 1), color = cluster == 3)) + geom_violin(scale = 'width')
df <- melt(t(as.matrix(vst_out$y[goi, ])), varnames = c('cell', 'gene'), value.name = 'Pearson_residual')
df$cluster <- clustering
ggplot(df, aes(x = gene, y = Pearson_residual, color = cluster == 3)) + geom_violin(scale = 'width')


