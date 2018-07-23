#' Compare gene expression between two groups
#'
#' @param x A list that provides model parameters and optionally meta data; use output of vst function
#' @param umi A matrix of UMI counts with genes as rows and cells as columns
#' @param group A vector indicating the groups
#' @param val1 A vector indicating the values of the group vector to treat as group 1
#' @param val2 A vector indicating the values of the group vector to treat as group 2
#' @param method Either 'LRT' for likelihood ratio test, or 't_test' for t-test
#' @param bin_size Number of genes that are processed between updates of progress bar
#' @param cell_attr Data frame of cell meta data
#' @param min_cells A gene has to be detected in at least this many cells in at least one of the groups being compared to be tested
#' @param show_progress Show progress bar
#'
#' @return Data frame of results
#'
#' @import Matrix
#' @import parallel
#' @importFrom stats model.matrix
#'
#' @export
#'
#' @examples
#' \dontrun{
#' vst_out <- vst(pbmc, return_cell_attr = TRUE)
#' # create fake clusters
#' clustering <- 1:ncol(pbmc) %/% 100
#' res <- compare_expression(x = vst_out, umi = pbmc, group = clustering, val1 = 0, val2 = 3)
#' }
#'
compare_expression <- function(x, umi, group, val1, val2, method = 'LRT', bin_size = 256,
                               cell_attr = x$cell_attr, min_cells = 5, show_progress = TRUE) {
  if (! method %in% c('LRT', 't_test')) {
    stop('method needs to be either \'LRT\' or \'t_test\'')
  }
  do_fast = FALSE
  if (method == 't_test') {
    do_fast = TRUE
  }
  sel1 <- which(group %in% val1)
  sel2 <- which(group %in% val2)
  use_cells <- c(sel1, sel2)
  cell_attr <- cell_attr[use_cells, ]
  group <- factor(c(rep(0, length(sel1)), rep(1, length(sel2))))
  weights <- c(rep(1/length(sel1), length(sel1)), rep(1/length(sel2), length(sel2)))
  genes <- rownames(x$model_pars_fit)[rownames(x$model_pars_fit) %in% rownames(umi)]
  cells_group1 <- rowSums(umi[genes, sel1] > 0)
  cells_group2 <- rowSums(umi[genes, sel2] > 0)
  genes <- genes[cells_group1 >= min_cells | cells_group2 >= min_cells]
  message('Testing for differential gene expression between two groups')
  message('Cells in group 1: ', length(sel1))
  message('Cells in group 2: ', length(sel2))
  message('Testing ', length(genes), ' genes')

  regressor_data <- model.matrix(as.formula(gsub('^y', '', x$model_str)), cell_attr)
  # process genes in batches
  bin_ind <- ceiling(x = 1:length(x = genes) / bin_size)
  max_bin <- max(bin_ind)
  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  res <- list()
  for (i in 1:max_bin) {
    genes_bin <- genes[bin_ind == i]
    if (do_fast) {
      bin_res <- mclapply(genes_bin, function(gene) {
        model_comparison_ttest(x$y[gene, use_cells], group)
      })
    } else {
      mu <- x$model_pars_fit[genes_bin, -1, drop=FALSE] %*% t(regressor_data)  # in log space
      y <- as.matrix(umi[genes_bin, use_cells])
      bin_res <- mclapply(genes_bin, function(gene) {
        model_comparison_lrt(y[gene, ], mu[gene, ], x$model_pars_fit[gene, 'theta'], group, weights)
      })
    }
    res[[i]] <- do.call(rbind, bin_res)
    if (show_progress) {
      setTxtProgressBar(pb, i)
    }
  }
  if (show_progress) {
    close(pb)
  }
  res <- do.call(rbind, res)
  rownames(res) <- genes
  colnames(res) <- c('p_value', 'log_fc')
  res <- as.data.frame(res)
  res <- res[order(res$p_value), ]
  return(res)
}

#' @importFrom stats glm offset anova
#' @importFrom MASS negative.binomial
model_comparison_lrt <- function(y, offs, theta, group, weights = NULL) {
  fam <- negative.binomial(theta = theta)
  mod0 <- glm(y ~ 1 + offset(offs), family = fam, weights = weights)
  mod1 <- glm(y ~ 1 + offset(offs) + group, family = fam, weights = weights)
  p_val <- anova(mod0, mod1, test = 'LRT')$'Pr(>Chi)'[2]
  fold_change <- log2(exp(mod1$coefficients[2]))
  return(c(p_val, fold_change))
}

#' @importFrom stats t.test
model_comparison_ttest <- function(y, group) {
  tt <- t.test(y ~ group)
  return(c(tt$p.value, diff(tt$estimate)))
}


