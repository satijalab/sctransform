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
#' @param y Only used if methtod = 't_test', this is the residual matrix; default is x$y
#' @param min_cells A gene has to be detected in at least this many cells in at least one of the groups being compared to be tested
#' @param weighted Balance the groups by using the appropriate weights
#' @param randomize Boolean indicating whether to shuffle group labels - only set to TRUE when testing methods
#' @param verbose Whether to show messages; default is TRUE
#' @param show_progress Show progress bar; default is same as verbose
#'
#' @return Data frame of results
#'
#' @import Matrix
#' @importFrom future.apply future_lapply
#' @importFrom stats model.matrix p.adjust pchisq
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
                               cell_attr = x$cell_attr, y = x$y, min_cells = 5,
                               weighted = TRUE, randomize = FALSE, verbose = TRUE,
                               show_progress = verbose) {
  if (! method %in% c('LRT', 'LRT_free', 'LRT_reg', 't_test')) {
    stop('method needs to be either \'LRT\', \'LRT_free\', \'LRT_reg\' or \'t_test\'')
  }
  if ('DE_test_group' %in% colnames(cell_attr)) {
    stop('DE_test_group cannot be a column name in cell attributes')
  }
  sel1 <- which(group %in% val1)
  sel2 <- which(group %in% val2)
  # randomize
  # if (randomize) {
  #   sel.rnd <- sample(x = c(sel1, sel2), replace = FALSE)
  #   sel1 <- sel.rnd[1:length(sel1)]
  #   sel2 <- sel.rnd[(length(sel1)+1):length(sel.rnd)]
  # }
  use_cells <- c(sel1, sel2)
  group <- factor(c(rep(0, length(sel1)), rep(1, length(sel2))))
  cell_attr <- cell_attr[use_cells, ]
  cell_attr$DE_test_group <- group
  if (weighted) {
    weights <- c(rep(1/length(sel1), length(sel1)), rep(1/length(sel2), length(sel2)))
    #weights <- c(rep(1/length(sel2), length(sel1)), rep(1/length(sel1), length(sel2)))
    weights <- weights / sum(weights) * length(use_cells)
  } else {
    weights <- rep(1, length(use_cells))
  }
  print(table(weights))
  genes <- rownames(x$model_pars_fit)[rownames(x$model_pars_fit) %in% rownames(umi)]
  cells_group1 <- rowSums(umi[genes, sel1] > 0)
  cells_group2 <- rowSums(umi[genes, sel2] > 0)
  genes <- genes[cells_group1 >= min_cells | cells_group2 >= min_cells]
  if (verbose) {
    message('Testing for differential gene expression between two groups')
    message('Cells in group 1: ', length(sel1))
    message('Cells in group 2: ', length(sel2))
    message('Testing ', length(genes), ' genes')
  }

  regressor_data <- model.matrix(as.formula(gsub('^y', '', x$model_str)), cell_attr)
  if (!is.null(dim(x$model_pars_nonreg))) {
    regressor_data_nonreg <- model.matrix(as.formula(gsub('^y', '', x$model_str_nonreg)), cell_attr)
    regressor_data <- cbind(regressor_data, regressor_data_nonreg)
  }

  # process genes in batches
  bin_ind <- ceiling(x = 1:length(x = genes) / bin_size)
  max_bin <- max(bin_ind)
  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  res <- list()
  for (i in 1:max_bin) {
    genes_bin <- genes[bin_ind == i]
    if (method == 't_test') {
      bin_res <- future_lapply(genes_bin, function(gene) {
        model_comparison_ttest(y[gene, use_cells], group)
      })
    }
    if (method == 'LRT') {
      mu <- x$model_pars_fit[genes_bin, -1, drop=FALSE] %*% t(regressor_data)  # in log space
      y <- as.matrix(umi[genes_bin, use_cells])
      bin_res <- future_lapply(genes_bin, function(gene) {
        model_comparison_lrt(y[gene, ], mu[gene, ], x$model_pars_fit[gene, 'theta'], group, weights)
      })
    }
    if (method == 'LRT_reg') {
      LB <- min(x$genes_log_mean_step1)
      UB <- max(x$genes_log_mean_step1)

      y <- as.matrix(umi[genes_bin, use_cells, drop=FALSE])
      if (randomize) {
        y <- t(apply(y, 1, sample))
        #y <- t(apply(y, 1, function(x) ceiling(pmax(0, rnorm(n = length(x), mean = 0, sd = 2)))))
      }
      # get estimated model parameters and expected counts for all cells combined
      #y_log_mean <- log10(base::rowMeans(y))
      y_log_mean <- log10(apply(y, 1, function(x) mean(x * weights)))
      y_log_mean <- pmax(LB, pmin(y_log_mean, UB))
      names(y_log_mean) <- rownames(y)
      mp <- reg_pars(x$genes_log_mean_step1, x$model_pars, y_log_mean, x$arguments$bw_adjust)
      if (!is.null(dim(x$model_pars_nonreg))) {
        mp <- cbind(mp, x$model_pars_nonreg[genes_bin, ])
      }
      mu <- exp(tcrossprod(mp[, -1, drop=FALSE], regressor_data))
      sq_dev <- sapply(1:nrow(mu), function(i) sq_deviance_residual(y[i, ], mu[i, ], mp[i, 'theta']))

      # same per group
      y0 <- y[, group==0]
      y_log_mean0 <- log10(base::rowMeans(y0))
      y_log_mean0 <- pmax(LB, pmin(y_log_mean0, UB))
      names(y_log_mean0) <- rownames(y)
      mp0 <- reg_pars(x$genes_log_mean_step1, x$model_pars, y_log_mean0, x$arguments$bw_adjust)
      if (!is.null(dim(x$model_pars_nonreg))) {
        mp0 <- cbind(mp0, x$model_pars_nonreg[genes_bin, ])
      }
      mu0 <- exp(tcrossprod(mp0[, -1, drop=FALSE], regressor_data[group==0, ]))
      sq_dev0 <- sapply(1:nrow(mu0), function(i) sq_deviance_residual(y0[i, ], mu0[i, ], mp0[i, 'theta']))

      y1 <- y[, group==1]
      y_log_mean1 <- log10(base::rowMeans(y1))
      y_log_mean1 <- pmax(LB, pmin(y_log_mean1, UB))
      names(y_log_mean1) <- rownames(y)
      mp1 <- reg_pars(x$genes_log_mean_step1, x$model_pars, y_log_mean1, x$arguments$bw_adjust)
      if (!is.null(dim(x$model_pars_nonreg))) {
        mp1 <- cbind(mp1, x$model_pars_nonreg[genes_bin, ])
      }
      mu1 <- exp(tcrossprod(mp1[, -1, drop=FALSE], regressor_data[group==1, ]))
      sq_dev1 <- sapply(1:nrow(mu1), function(i) sq_deviance_residual(y1[i, ], mu1[i, ], mp1[i, 'theta']))

      #pvals <- pchisq(base::rowSums(cbind(sq_dev0, sq_dev1)) - base::rowSums(sq_dev), df = 1, lower.tail = FALSE)
      pvals <- pchisq(base::colSums(sq_dev * weights) -
                      base::colSums(rbind(sq_dev0, sq_dev1) * weights),
                      df = 3, lower.tail = FALSE)

      #fold_change <- log2(10 ^ (y_log_mean1 - y_log_mean0))
      # tmp stuff for fold change
      mu0 <- tcrossprod(mp0[, -1, drop=FALSE], regressor_data)
      mu1 <- tcrossprod(mp1[, -1, drop=FALSE], regressor_data)
      fold_change <- apply(log2(exp(mu1 - mu0)), 1, mean)
      #if (max(fold_change) > 0.4) browser()
      if ('SON' %in% genes_bin) browser()
      bin_res <- list(cbind(pvals, fold_change))
    }
    if (method == 'LRT_free') {
      y <- as.matrix(umi[genes_bin, use_cells])
      # get estimated theta
      bw <- bw.SJ(x$genes_log_mean_step1)
      y_log_mean <- log10(base::rowMeans(y))
      o <- order(y_log_mean)
      y_theta <- rep(NA_real_, nrow(y))
      y_theta[o] <- 10 ^ ksmooth(x = x$genes_log_mean_step1, y = log10(x$model_pars[, 'theta']),
                                 x.points = y_log_mean, bandwidth = bw, kernel='normal')$y
      names(y_theta) <- genes_bin
      bin_res <- future_lapply(genes_bin, function(gene) {
        return(model_comparison_lrt_free3(gene, y[gene, ], y_theta[gene], x$model_str, cell_attr, group, weights, randomize))
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
  res$fdr <- p.adjust(res$p_value, method='fdr')
  res <- res[order(res$p_value, -abs(res$log_fc)), ]
  res$mean1 <- rowMeans(umi[rownames(res), sel1])
  res$mean2 <- rowMeans(umi[rownames(res), sel2])
  res$mean <- rowMeans(umi[rownames(res), use_cells])
  res$mean_weighted <- (res$mean1 + res$mean2) / 2
  return(res)
}

compare_expression_full <- function(umi, cell_attr, group, val1, val2,
                                    latent_var = c('log_umi'),
                                    batch_var = NULL,
                                    latent_var_nonreg = NULL,
                                    n_genes = 2000,
                                    method = 'poisson',
                                    bin_size = 256,
                                    min_cells = 3,
                                    bw_adjust = 2,
                                    min_frac = 0,
                                    verbose = TRUE,
                                    show_progress = verbose) {
  sel1 <- which(group %in% val1)
  sel2 <- which(group %in% val2)

  det1 <- rowMeans(umi[, sel1] > 0)
  det2 <- rowMeans(umi[, sel2] > 0)
  umi <- umi[det1 >= min_frac | det2 >= min_frac, ]

  cells1 <- rowSums(umi[, sel1] > 0)
  cells2 <- rowSums(umi[, sel2] > 0)
  umi <- umi[cells1 >= min_cells | cells2 >= min_cells, ]

  vst.out0 <- vst(umi = umi[, c(sel1, sel2)],
                  cell_attr = cell_attr[c(sel1, sel2), ],
                  latent_var = latent_var,
                  batch_var = batch_var,
                  latent_var_nonreg = latent_var_nonreg,
                  n_genes = n_genes,
                  n_cells = NULL,
                  method = method,
                  do_regularize = TRUE,
                  res_clip_range = c(-Inf, Inf),
                  bin_size = bin_size,
                  min_cells = min_cells,
                  return_cell_attr = FALSE,
                  return_gene_attr = FALSE,
                  residual_type = 'deviance',
                  bw_adjust = bw_adjust,
                  verbose = verbose,
                  show_progress = show_progress)

  vst.out1 <- vst(umi = umi[, sel1],
                  cell_attr = cell_attr[sel1, ],
                  latent_var = latent_var,
                  batch_var = batch_var,
                  latent_var_nonreg = latent_var_nonreg,
                  n_genes = n_genes,
                  n_cells = NULL,
                  method = 'nb_theta_given', #method,
                  do_regularize = TRUE,
                  res_clip_range = c(-Inf, Inf),
                  bin_size = bin_size,
                  min_cells = min_cells,
                  return_cell_attr = FALSE,
                  return_gene_attr = FALSE,
                  residual_type = 'deviance',
                  bw_adjust = bw_adjust,
                  theta_given = vst.out0$model_pars_fit[, 'theta'],
                  verbose = verbose,
                  show_progress = show_progress)

  vst.out2 <- vst(umi = umi[, sel2],
                  cell_attr = cell_attr[sel2, ],
                  latent_var = latent_var,
                  batch_var = batch_var,
                  latent_var_nonreg = latent_var_nonreg,
                  n_genes = n_genes,
                  n_cells = NULL,
                  method = 'nb_theta_given', #method
                  do_regularize = TRUE,
                  res_clip_range = c(-Inf, Inf),
                  bin_size = bin_size,
                  min_cells = min_cells,
                  return_cell_attr = FALSE,
                  return_gene_attr = FALSE,
                  residual_type = 'deviance',
                  bw_adjust = bw_adjust,
                  theta_given = vst.out0$model_pars_fit[, 'theta'],
                  verbose = verbose,
                  show_progress = show_progress)

  genes <- union(rownames(vst.out1$y), rownames(vst.out2$y))
  genes_both <- intersect(rownames(vst.out1$y), rownames(vst.out2$y))
  genes1 <- setdiff(rownames(vst.out1$y), genes_both)
  genes2 <- setdiff(rownames(vst.out2$y), genes_both)

  sq_dev_one <- base::rowSums(vst.out0$y[genes, ]^2 * 1)
  sq_dev_two <- rep(0, length(sq_dev_one))
  names(sq_dev_two) <- genes
  sq_dev_two[rownames(vst.out1$y)] <- base::rowSums(vst.out1$y^2 * 1)
  sq_dev_two[rownames(vst.out2$y)] <- sq_dev_two[rownames(vst.out2$y)] + base::rowSums(vst.out2$y^2 * 1)
  pvals <- pchisq(sq_dev_one - sq_dev_two,
                  df = 3, lower.tail = FALSE)
  # get log-fold change
  log_fc <- rep(NA_real_, length(sq_dev_one))
  names(log_fc) <- genes

  regressor_data <- model.matrix(as.formula(gsub('^y', '', vst.out0$model_str)), cell_attr[c(sel1, sel2), ])
  if (!is.null(dim(vst.out0$model_pars_nonreg))) {
    regressor_data_nonreg <- model.matrix(as.formula(gsub('^y', '', vst.out0$model_str_nonreg)), cell_attr[c(sel1, sel2), ])
    regressor_data <- cbind(regressor_data, regressor_data_nonreg)
  }
  mp1 <- cbind(vst.out1$model_pars_fit, vst.out1$model_pars_nonreg)
  mp2 <- cbind(vst.out2$model_pars_fit, vst.out2$model_pars_nonreg)
  mu1 <- tcrossprod(mp1[genes_both, -1, drop=FALSE], regressor_data)
  mu2 <- tcrossprod(mp2[genes_both, -1, drop=FALSE], regressor_data)
  log_fc[genes_both] <- apply(log2(exp(mu2 - mu1)), 1, mean)
  log_fc[genes1] <- -Inf
  log_fc[genes2] <- Inf

  res <- data.frame(p_value = pvals, log_fc = log_fc)
  res$fdr <- p.adjust(res$p_value, method='fdr')
  res <- res[order(res$p_value, -abs(res$log_fc)), ]
  res$mean1 <- rowMeans(umi[rownames(res), sel1])
  res$mean2 <- rowMeans(umi[rownames(res), sel2])
  res$det1 <- rowMeans(umi[rownames(res), sel1] > 0)
  res$det2 <- rowMeans(umi[rownames(res), sel2] > 0)

  # tmp stuff
  # goi <- 'MALAT1'
  # y <- umi[goi, c(sel1, sel2)]
  # grp <- c(rep('A', length(sel1)), rep('B', length(sel2)))
  # df <- data.frame(y=y, log_umi=cell_attr[c(sel1, sel2), 'log_umi'], grp=grp)
  # mod0 <- glm.nb(y ~ log_umi, data = df)
  # mod1 <- glm.nb(y ~ log_umi + grp, data = df)
  # mod1 <- glm(y ~ log_umi + grp, data = df, family = negative.binomial(theta=mod0$theta))
  # mod1 <- glm(y ~ log_umi:grp, data = df, family = negative.binomial(theta=mod0$theta))

  return(res)
}


# function to get regularized model parameters
reg_pars <- function(x, y.mat, x.points, bw.adjust) {
  bw <- bw.SJ(x) * bw.adjust
  o <- order(x.points)
  y.mat.out <- matrix(NA_real_, length(x.points), ncol(y.mat))
  y.mat.out[o, 1] <- 10 ^ ksmooth(x = x, y = log10(y.mat[, 1]), x.points = x.points,
                                  bandwidth = bw*3, kernel='normal')$y
  for (i in 2:ncol(y.mat)) {
    y.mat.out[o, i] <- ksmooth(x = x, y = y.mat[, i], x.points = x.points,
                               bandwidth = bw, kernel='normal')$y
  }
  colnames(y.mat.out) <- colnames(y.mat)
  rownames(y.mat.out) <- names(x.points)
  if (any(apply(is.na(y.mat.out), 1, any))) {
    browser()
  }
  return(y.mat.out)
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

# fixed overdispersion (theta)
# different slopes
model_comparison_lrt_free1 <- function(gene, y, theta, model_str, cell_attr, group, weights = NULL, randomize = FALSE) {
  if (randomize) {
    y <- sample(y)
  }
  mod0 <- glm(as.formula(model_str), data = cell_attr, family = negative.binomial(theta=theta), weights = weights)
  mod1_str <- paste0(model_str, ' + DE_test_group')
  mod1 <- glm(as.formula(mod1_str), data = cell_attr, family = negative.binomial(theta=theta), weights = weights)
  p_val <- anova(mod0, mod1, test = 'Chisq', dispersion = 1)$'Pr(>Chi)'[2]
  fold_change <- log2(exp(rev(mod1$coefficients)[1]))
  return(c(p_val, fold_change))
}

# fixed overdispersion (theta)
# fixed slopes
#' @importFrom stats pchisq
model_comparison_lrt_free2 <- function(gene, y, theta, model_str, cell_attr, group, weights = NULL, randomize = FALSE) {
  if (randomize) {
    y <- sample(y)
  }
  mod0 <- glm(as.formula(model_str), data = cell_attr, family = negative.binomial(theta=theta), weights = weights)
  offs <- log(mod0$fitted.values) - mod0$coefficients[1]
  mod1 <- glm(y ~ 1 + offset(offs) + group, family = negative.binomial(theta=theta), weights = weights)
  deviance_diff <- mod0$deviance - mod1$deviance
  p_val <- pchisq(q = deviance_diff, df = 1, lower.tail = FALSE)
  fold_change <- log2(exp(rev(mod1$coefficients)[1]))
  return(c(p_val, fold_change))
}

# fixed overdispersion (theta)
# different per-group slopes
#' @importFrom stats predict
model_comparison_lrt_free3 <- function(gene, y, theta, model_str, cell_attr, group, weights = NULL, randomize = FALSE) {
  if (randomize) {
    y <- sample(y)
  }
  mod0 <- glm(as.formula(model_str), data = cell_attr, family = negative.binomial(theta=theta), weights = weights)
  mod1_str <- paste(c('y ~', '(', gsub('^y ~ ', '', model_str), ') : DE_test_group + DE_test_group'), collapse=' ')
  mod1 <- glm(as.formula(mod1_str), data = cell_attr, family = negative.binomial(theta=theta), weights = weights)
  p_val <- anova(mod0, mod1, test = 'Chisq', dispersion = 1)$'Pr(>Chi)'[2]
  # to get fold change, predict data
  tmp.ca0 <- cell_attr
  tmp.ca0$DE_test_group <- factor(0)
  tmp.ca1 <- cell_attr
  tmp.ca1$DE_test_group <- factor(1)
  fold_change <- log2(median(predict(mod1, newdata = tmp.ca1, type = 'response')/predict(mod0, newdata = tmp.ca0, type = 'response')))
  return(c(p_val, fold_change))
}

model_comparison_lrt_free <- function(gene, y, theta, model_str, cell_attr, group, weights = NULL) {
  #print(gene)
  # model 0
  #mod0 <- MASS::glm.nb(as.formula(model_str), data = cell_attr, weights = weights)
  #fit1 <- glm(as.formula(model_str), data = cell_attr, family = poisson, weights = weights)
  #theta1 <- as.numeric(x = theta.ml(y = y, mu = fit1$fitted, weights = weights))
  #theta1b <- max(0.1, theta1)
  mod0 <- 0
  try(mod0 <- glm(as.formula(model_str), data = cell_attr, family = negative.binomial(theta=theta), weights = weights), silent = TRUE)
  if (class(mod0)[1] == 'numeric') {
    print('mod0 failed')
    browser()
  }

  # model 1
  #mod1_str <- paste(c('y ~', '(', gsub('^y ~ ', '', model_str), ') : DE_test_group'), collapse=' ')
  #mod1_str <- paste0(model_str, ' + DE_test_group')
  #mod1 <- MASS::glm.nb(as.formula(mod1_str), data = cell_attr, weights = weights)
  #fit2 <- glm(as.formula(mod1_str), data = cell_attr, family = poisson, weights = weights)
  #theta2 <- as.numeric(x = theta.ml(y = y, mu = fit2$fitted, weights = weights))
  #theta2b <- max(0.1, theta2)
  #mod1 <- 0
  #try(mod1 <- glm(as.formula(mod1_str), data = cell_attr, family = negative.binomial(theta=theta), weights = weights), silent = TRUE)
  #if (class(mod1)[1] == 'numeric') {
  #  print('mod1 failed')
  #  browser()
  #}

  #if (sum(y[group==0]) == 0 | sum(y[group==1]) == 0) {
  # print(theta1)
  # print(theta1b)
  # print(theta2)
  # print(theta2b)
  #print(anova(mod0, mod1, test = 'Chisq'))
  #browser()
  #}

  #print(mod0)
  #print(mod1)
  #print(anova(mod0, mod1, test = 'Chisq'))

  #p_val <- anova(mod0, mod1, test = 'Chisq')$'Pr(Chi)'[2]
  #p_val <- anova(mod0, mod1, test = 'Chisq')$'Pr(>Chi)'[2]
  #fold_change <- log2(exp(rev(mod1$coefficients)[1]))

  # alternative model 1 and p-value calculation
  #mod0.o <- glm(y ~ 1 + offset(log(mod0$fitted.values)), family = negative.binomial(theta=theta), weights = weights)
  #offs <- predict(mod0, newdata = cell_attr) - mod0$coefficients[1]
  offs <- log(mod0$fitted.values) - mod0$coefficients[1]
  mod1.o <- glm(y ~ 1 + offset(offs) + group, family = negative.binomial(theta=theta), weights = weights)


  grp.intercept <- mod1.o$coefficients
  if (grp.intercept[1] > grp.intercept[2]) {
    p_val <- summary(mod1.o)$coefficients[1, 4]
    fold_change <- log2(exp(diff(grp.intercept)))
  } else {
    p_val <- summary(mod1.o)$coefficients[2, 4]
    fold_change <- log2(exp(diff(grp.intercept)))
  }
  #mod1.o <- glm(y ~ 1 + offset(log(mod0$fitted.values)) + group, family = negative.binomial(theta=theta), weights = weights)
  #p_val <- anova(mod0.o, mod1.o, test = 'Chisq')$'Pr(>Chi)'[2]
  #deviance_diff <- sum(residuals(mod0, type='deviance')^2) - sum(residuals(mod1.o, type='deviance')^2)
  #p_val <- pchisq(q = deviance_diff, df = 1, lower.tail = FALSE)
  #fold_change <- log2(exp(rev(mod1.o$coefficients)[1]))

  # if (mean(y[group==0]) > 0.03 & mean(y[group==1]) == 0) {
  if (gene == 'OGFOD1') {
    browser()
  }

  return(c(p_val, fold_change))
}

#' @importFrom stats t.test
model_comparison_ttest <- function(y, group) {
  tt <- t.test(y ~ group)
  return(c(tt$p.value, diff(tt$estimate)))
}


