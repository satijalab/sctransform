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
old_theme <- theme_set(theme_classic(base_size=8))

## ----load_data-----------------------------------------------------------
pbmc_data <- readRDS(file = "~/Projects/data/pbmc3k_umi_counts.Rds")
class(pbmc_data)
dim(pbmc_data)

## ----calc_attributes-----------------------------------------------------
gene_attr <- data.frame(mean = rowMeans(pbmc_data), 
                        detection_rate = rowMeans(pbmc_data > 0),
                        var = apply(pbmc_data, 1, var))
gene_attr$log_mean <- log10(gene_attr$mean)
gene_attr$log_var <- log10(gene_attr$var)
rownames(gene_attr) <- rownames(pbmc_data)
cell_attr <- data.frame(n_umi = colSums(pbmc_data),
                        n_gene = colSums(pbmc_data > 0))
rownames(cell_attr) <- colnames(pbmc_data)

## ----mean_var_rel, warning=FALSE, fig.cap='Mean-variance relationship'----
ggplot(gene_attr, aes(log_mean, log_var)) + 
  geom_point(alpha=0.3, shape=16) + 
  geom_density_2d(size = 0.3) +
  geom_abline(intercept = 0, slope = 1, color='red')

## ----mean_dr_rel, warning=FALSE, fig.cap='Mean-detection-rate relationship'----
# add the expected detection rate under Poisson model
x = seq(from = -3, to = 2, length.out = 1000)
poisson_model <- data.frame(log_mean = x, detection_rate = 1 - dpois(0, lambda = 10^x))
ggplot(gene_attr, aes(log_mean, detection_rate)) + 
  geom_point(alpha=0.3, shape=16) + 
  geom_line(data=poisson_model, color='red') +
  theme_gray(base_size = 8)

## ----umi_gene_rel, warning=FALSE, fig.cap='UMI detected genes relationship'----
ggplot(cell_attr, aes(n_umi, n_gene)) + 
  geom_point(alpha=0.3, shape=16) + 
  geom_density_2d(size = 0.3)

## ---- fig.width=4, fig.height=2.5, warning=FALSE-------------------------
# We use the Future API for parallel processing; set parameters here
future::plan(strategy = 'multicore', workers = 4)
options(future.globals.maxSize = 10 * 1024 ^ 3)

set.seed(44)
vst_out <- sctransform::vst(pbmc_data, latent_var = c('log_umi'), return_gene_attr = TRUE, return_cell_attr = TRUE, show_progress = FALSE)
sctransform::plot_model_pars(vst_out)

## ---- fig.width=5, fig.height=3.5, warning=FALSE-------------------------
sctransform::plot_model(vst_out, pbmc_data, c('MALAT1', 'RPL10', 'FTL'), plot_residual = TRUE)

## ---- fig.width=5, fig.height=2, warning=FALSE---------------------------
sctransform::plot_model(vst_out, pbmc_data, c('FTL'), plot_residual = TRUE, show_nr = TRUE, arrange_vertical = FALSE)

## ---- fig.width=5, fig.height=4, warning=FALSE---------------------------
sctransform::plot_model(vst_out, pbmc_data, c('GNLY', 'S100A9'), plot_residual = TRUE, show_nr = TRUE)

## ---- fig.keep=TRUE------------------------------------------------------
ggplot(vst_out$gene_attr, aes(residual_mean)) + geom_histogram(binwidth=0.01)
ggplot(vst_out$gene_attr, aes(residual_variance)) + geom_histogram(binwidth=0.1) + geom_vline(xintercept=1, color='red') + xlim(0, 10)

## ------------------------------------------------------------------------
ggplot(vst_out$gene_attr, aes(log10(gmean), residual_variance)) + geom_point(alpha=0.3, shape=16) +
  geom_density_2d(size = 0.3)

## ------------------------------------------------------------------------
head(round(vst_out$gene_attr[order(-vst_out$gene_attr$residual_variance), ], 2), 22)

