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
  fig.width=6, fig.height=4, dpi=100, out.width = '70%'
)
old_theme <- theme_set(theme_classic(base_size=8))

## ---- warning=FALSE------------------------------------------------------
# some of the vst steps can use multiple cores
# We use the Future API for parallel processing; set parameters here
future::plan(strategy = 'multicore', workers = 4)
options(future.globals.maxSize = 10 * 1024 ^ 3)

cm <- readRDS('~/Projects/data/in-lineage_dropseq_CGE3_digitial_expression.Rds')

set.seed(42)
vst_out <- sctransform::vst(cm, latent_var = 'log_umi_per_gene', bin_size = 128, return_cell_attr = TRUE, show_progress = FALSE)

## ------------------------------------------------------------------------
pca <- irlba::prcomp_irlba(t(vst_out$y), n = 2)

# fit principal curve through first two PCs
pricu <- princurve::principal.curve(pca$x, smoother='lowess', f=0.5, stretch=333)
# cell projection onto curve is maturation score
maturation_score <- pricu$lambda/max(pricu$lambda)

## ---- fig.width=8, fig.height=4, out.width = '100%'----------------------
y_smooth <- sctransform::smooth_via_pca(vst_out$y, do_plot = TRUE)

## ------------------------------------------------------------------------
cm_corrected <- sctransform::correct(vst_out, data = y_smooth, show_progress = FALSE)

## ---- fig.width=7, fig.height=7, out.width='100%'------------------------
goi <- c('Nes', 'Ccnd2', 'Tuba1a')
df <- list()
df[[1]] <- melt(t(as.matrix(cm[goi, ])), varnames = c('cell', 'gene'), value.name = 'value')
df[[1]]$type <- 'UMI'
df[[1]]$maturation_rank <- rank(maturation_score)
df[[2]] <- melt(t(as.matrix(vst_out$y[goi, ])), varnames = c('cell', 'gene'), value.name = 'value')
df[[2]]$type <- 'Pearson residual'
df[[2]]$maturation_rank <- rank(maturation_score)
df[[3]] <- melt(t(as.matrix(y_smooth[goi, ])), varnames = c('cell', 'gene'), value.name = 'value')
df[[3]]$type <- 'corrected Pearson residual'
df[[3]]$maturation_rank <- rank(maturation_score)
df[[4]] <- melt(t(as.matrix(cm_corrected[goi, ])), varnames = c('cell', 'gene'), value.name = 'value')
df[[4]]$type <- 'corrected UMI'
df[[4]]$maturation_rank <- rank(maturation_score)
df <- do.call(rbind, df)
df$gene <- factor(df$gene, ordered=TRUE, levels=unique(df$gene))
df$type <- factor(df$type, ordered=TRUE, levels=unique(df$type))
ggplot(df, aes(maturation_rank, value)) + geom_point(alpha=0.5, shape=16) + 
  geom_density_2d(size=0.5) + facet_grid(type ~ gene, scales = 'free')

