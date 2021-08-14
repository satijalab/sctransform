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
#' @param verbosity An integer specifying whether to show only messages (1), messages and progress bars (2) or nothing (0) while the function is running; default is 2
#' @param verbose Deprecated; use verbosity instead
#' @param show_progress Deprecated; use verbosity instead
#'
#' @return Data frame of results
#'
#' @import Matrix
#' @importFrom future.apply future_lapply
#' @importFrom stats model.matrix p.adjust pchisq
#'
compare_expression <- function(x, umi, group, val1, val2, method = 'LRT', bin_size = 256,
                               cell_attr = x$cell_attr, y = x$y, min_cells = 5,
                               weighted = TRUE, randomize = FALSE, verbosity = 2,
                               verbose = NULL, show_progress = NULL) {
  # Take care of deprecated arguments
  if (!is.null(verbose)) {
    warning("The 'verbose' argument is deprecated as of v0.3. Use 'verbosity' instead. (in sctransform::vst)", immediate. = TRUE, call. = FALSE)
    verbosity <- as.numeric(verbose)
  }
  if (!is.null(show_progress)) {
    warning("The 'show_progress' argument is deprecated as of v0.3. Use 'verbosity' instead. (in sctransform::vst)", immediate. = TRUE, call. = FALSE)
    if (show_progress) {
      verbosity <- 2
    } else {
      verbosity <- min(verbosity, 1)
    }
  }

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
  if (verbosity > 0) {
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
  if (verbosity > 1) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  res <- list()
  for (i in 1:max_bin) {
    genes_bin <- genes[bin_ind == i]
    if (method == 't_test') {
      bin_res <- future_lapply(
        X = genes_bin, 
        FUN = function(gene) {model_comparison_ttest(y[gene, use_cells], group)},
        future.seed = TRUE)
    }
    if (method == 'LRT') {
      mu <- x$model_pars_fit[genes_bin, -1, drop=FALSE] %*% t(regressor_data)  # in log space
      y <- as.matrix(umi[genes_bin, use_cells])
      bin_res <- future_lapply(
        X = genes_bin, 
        FUN = function(gene) {
          model_comparison_lrt(y[gene, ], mu[gene, ], x$model_pars_fit[gene, 'theta'], group, weights)},
        future.seed = TRUE)
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
      bin_res <- future_lapply(
        X = genes_bin, 
        FUN = function(gene) {
          return(model_comparison_lrt_free3(gene, y[gene, ], y_theta[gene], x$model_str, cell_attr, group, weights, randomize))
        }, 
        future.seed = TRUE)
    }
    res[[i]] <- do.call(rbind, bin_res)
    if (verbosity > 1) {
      setTxtProgressBar(pb, i)
    }
  }
  if (verbosity > 1) {
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
                                    verbosity = 2) {
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
                  verbosity = verbosity)

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
                  verbosity = verbosity)

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
                  verbosity = verbosity)

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

#' Non-parametric differential expression test for sparse non-negative data
#'
#' @param y A matrix of counts; must be (or inherit from) class dgCMatrix; genes are row,
#' cells are columns
#' @param group_labels The group labels (e.g. cluster identities); 
#' will be converted to factor
#' @param compare Specifies which groups to compare, see details; default is 'each_vs_rest'
#' @param R The number of random permutations used to derive the p-values; default is 99
#' @param log2FC_th Threshold to remove genes from testing; absolute log2FC must be at least
#' this large for a gene to be tested; default is \code{log2(1.2)}
#' @param mean_th Threshold to remove genes from testing; gene mean must be at least this
#' large for a gene to be tested; default is 0.05
#' @param cells_th Threshold to remove genes from testing; gene must be detected (non-zero count)
#' in at least this many cells in the group with higher mean; default is 5
#' @param only_pos Test only genes with positive fold change (mean in group 1 > mean in group2); 
#' default is FALSE
#' @param only_top_n Test only the this number of genes from both ends of the log2FC spectrum
#' after all of the above filters have been applied; useful to get only the top markers; 
#' only used if set to a numeric value; default is NULL
#' @param mean_type Which type of mean to use; if \code{'geometric'} (default) the geometric mean is
#' used; to avoid \code{log(0)} we use \code{log1p} to add 1 to all counts and log-transform, 
#' calculate the arithmetic mean, and then back-transform and subtract 1 using \code{exp1m}; if
#' this parameter is set to \code{'arithmetic'} the data is used as is
#' @param verbosity Integer controlling how many messages the function prints; 
#' 0 is silent, 1 (default) is not
#'
#' @return Data frame of results
#' 
#' @section Details:
#' This model-free test is applied to each gene (row) individually but is
#' optimized to make use of the efficient sparse data representation of
#' the input. A permutation null distribution us used to assess the 
#' significance of the observed difference in mean between two groups.
#' 
#' The observed difference in mean is compared against a distribution
#' obtained by random shuffling of the group labels. For each gene every 
#' random permutation yields a difference in mean and from the population of
#' these background differences we estimate a mean and standard
#' deviation for the null distribution. 
#' This mean and standard deviation are used to turn the observed
#' difference in mean into a z-score and then into a p-value. Finally,
#' all p-values (for the tested genes) are adjusted using the Benjamini & Hochberg
#' method (fdr). The log2FC values in the output are \code{log2(mean1 / mean2)}.
#' Empirical p-values are also calculated: \code{emp_pval = (b + 1) / (R + 1)}
#' where b is the number of times the absolute difference in mean from a random 
#' permutation is at least as large as the absolute value of the observed difference
#' in mean, R is the number of random permutations. This is an upper bound of
#' the real empirical p-value that would be obtained by enumerating all possible
#' group label permutations.
#' 
#' There are multiple ways the group comparisons can be specified based on the compare
#' parameter. The default, \code{'each_vs_rest'}, does multiple comparisons, one per 
#' group vs all remaining cells. \code{'all_vs_all'}, also does multiple comparisons, 
#' covering all groups pairs. If compare is set to a length two character vector, e.g.
#' \code{c('T-cells', 'B-cells')}, one comparison between those two groups is done.
#' To put multiple groups on either side of a single comparison, use a list of length two. 
#' E.g. \code{compare = list(c('cluster1', 'cluster5'), c('cluster3'))}.
#' 
#' @import Matrix
#' @importFrom matrixStats rowMeans2 rowSds
#' @importFrom stats p.adjust pnorm
#' 
#' @export
#'
#' @examples
#' \donttest{
#' clustering <- 1:ncol(pbmc) %% 2
#' vst_out <- vst(pbmc, return_corrected_umi = TRUE)
#' de_res <- diff_mean_test(y = vst_out$umi_corrected, group_labels = clustering)
#' }
#'
diff_mean_test <- function(y, group_labels, 
                           compare = 'each_vs_rest', 
                           R = 99, log2FC_th = log2(1.2), 
                           mean_th = 0.05, cells_th = 5, only_pos = FALSE,
                           only_top_n = NULL,
                           mean_type = 'geometric', 
                           verbosity = 1) {
  if (is.na(match(x = mean_type, table = c('geometric', 'arithmetic')))) {
    stop('mean_type must be geometric or arithmetic')
  }
  if (!inherits(x = y, what = 'dgCMatrix')) {
    stop('y must be a dgCMatrix')
  }
  if (R < 13) {
    stop('R must be at least 13')
  }
  if (!is.null(only_top_n) & (!is.numeric(only_top_n) | length(only_top_n) > 1)) {
    stop('only_top_n must be NULL or a single numeric value')
  }
  group_labels <- droplevels(as.factor(group_labels))
  lab_tab <- table(group_labels)
  group_levels <- levels(group_labels)
  G <- length(group_levels)
  if (length(group_labels) != ncol(y)) {
    stop('length of group labels must be equal to the number of columns in y')
  }
  
  if (verbosity > 0) {
    message('Non-parametric DE test for count data')
    message(sprintf('Using %s mean and %d random permutations', mean_type, R))
    message('Input: ', nrow(y), ' genes, ', ncol(y), ' cells; ', G, ' groups')
  }
  
  # Set up the comparisons we want to do; each comparison is a list
  # name1, name2, labels grp1, labels grp2
  if (compare[1] == 'each_vs_rest' && G == 2) {
    compare <- group_levels
    if (verbosity > 0) {
      message('There are only two groups in the data. Changing compare argument from "each_vs_rest" to group levels')
    }
  }
  if (compare[1] == 'each_vs_rest') {
    comparisons <- lapply(group_levels, function(x) list(x, 'rest', x, setdiff(group_levels, x)))
  } else if (compare[1] == 'all_vs_all') {
    comparisons <- list()
    for (i in 1:(G-1)) {
      for (j in (i+1):G) {
        comparisons[[length(comparisons) + 1]] <- list(group_levels[i], group_levels[j], group_levels[i], group_levels[j])
      }
    }
  } else if (inherits(x = compare, what = 'character') &&
      length(compare) == 2 &&
      all(compare %in% group_levels))  {
    if (compare[1] == compare[2]) {
      stop('Group 1 and 2 need to be different - please check your compare argument')
    }
    comparisons <- list(list(compare[1], compare[2], compare[1], compare[2]))
  } else if (inherits(x = compare, what = 'list') &&
             length(compare) == 2 &&
             all(unlist(lapply(compare, inherits, what = 'character')))) {
    compare <- lapply(compare, unique)
    if (length(intersect(compare[[1]], compare[[2]])) > 0) {
      stop('Intersection between group 1 and 2 - please check your compare argument')
    }
    comparisons <- list(list('group1', 'group2', compare[[1]], compare[[2]]))
  } else {
    stop("Make sure the compare argument is 'each_vs_rest' or 'all_vs_all' or a length 2 
          character vector with both entries present in the group_labels argument or 
          a list of length 2 with each entry being a character vector of group labels")
  }
  
  # for all the genes, get the number of non-zero observations per group
  cells <- row_nonzero_count_grouped_dgcmatrix(matrix = y, group = group_labels)
  # if we want to use the geometric mean, it's fastest to convert all counts to
  # log1p upfront, then use expm1 of arithmetic mean later on
  if (mean_type == 'geometric') {
    y@x <- log(y@x + 1)
    means <- row_mean_grouped_dgcmatrix(matrix = y, group = group_labels, shuffle = FALSE)
  } else {
    means <- row_mean_grouped_dgcmatrix(matrix = y, group = group_labels, shuffle = FALSE)
  }
  
  # Run the test for each comparison
  res_lst <- lapply(comparisons, function(comp) {
    # we might only be using a subset of the input cells; set up here
    sel_columns1 <- group_labels %in% comp[[3]]
    sel_columns2 <- group_labels %in% comp[[4]]
    sel_columns <- sel_columns1 | sel_columns2
    comp_group_labels <- factor(sel_columns2[sel_columns])
    
    if (verbosity > 0) {
      message(sprintf('Comparing %s (group1, N = %d) to %s (group2, N = %d)',
                     comp[[1]], sum(sel_columns1), comp[[2]], sum(sel_columns2)))
    }
    if (sum(sel_columns1) == 0 || sum(sel_columns2) == 0) {
      return()
    }
    comp_cells <- do.call(cbind, lapply(comp[3:4], function(x) combine_counts(cells, x)))
    comp_means <- do.call(cbind, lapply(comp[3:4], function(x) combine_means(means, lab_tab, x, mean_type)))
    
    res <- data.frame(gene = rownames(means),
                      group1 = comp[[1]],
                      mean1 = comp_means[, 1],
                      cells1 = comp_cells[, 1],
                      group2 = comp[[2]],
                      mean2 = comp_means[, 2],
                      cells2 = comp_cells[, 2])
    res$mean_diff <- res$mean1 - res$mean2
    res$log2FC <- log2(res$mean1 / res$mean2)
    
    # remove genes according to the filters
    if (log2FC_th > 0 || mean_th > 0 || cells_th > 0 || only_pos || !is.null(only_top_n)) {
      sel1 <- abs(res$log2FC) >= log2FC_th
      sel2 <- res$mean1 >= mean_th | res$mean2 >= mean_th
      sel3 <- (res$log2FC >= 0 & res$cells1 >= cells_th) | (res$log2FC <= 0 & res$cells2 >= cells_th)
      if (only_pos) {
        sel4 <- res$log2FC > 0
      } else {
        sel4 <- TRUE
      }
      res <- res[sel1 & sel2 & sel3 & sel4, , drop = FALSE]
      if (!is.null(only_top_n)) {
        sel0 <- rank(-res$log2FC) <= only_top_n
        if (!only_pos) {
          sel0 <- sel0 | rank(res$log2FC) <= only_top_n
        }
        res <- res[sel0, , drop = FALSE]
      }
      if (verbosity > 0) {
        message(sprintf('Keeping %d genes after initial filtering', nrow(res)))
      }
      # handle the case where no genes remain after filtering
      if (nrow(res) == 0) {
        return(res)
      }
    }
    
    # now get the empirical null distribution for mean_diff
    y_ss <- y[rownames(res), sel_columns, drop = FALSE]
    if (mean_type == 'geometric') {
      mean_diff_rnd <- do.call(cbind, lapply(1:R, function(i) {
        means_r <- expm1(row_mean_grouped_dgcmatrix(matrix = y_ss, group = comp_group_labels, shuffle = TRUE))
        means_r[, 1, drop = FALSE] - means_r[, 2, drop = FALSE]
      }))
    } else {
      mean_diff_rnd <- do.call(cbind, lapply(1:R, function(i) {
        means_r <- row_mean_grouped_dgcmatrix(matrix = y_ss, group = comp_group_labels, shuffle = TRUE)
        means_r[, 1, drop = FALSE] - means_r[, 2, drop = FALSE]
      }))
    }
    
    # use null distribution to get empirical p-values
    # also approximate null with normal and derive z-scores and p-values
    res$emp_pval <- (rowSums((abs(mean_diff_rnd) - abs(res$mean_diff)) >= 0) + 1) / (R + 1)
    res$emp_pval_adj <- p.adjust(res$emp_pval, method = 'BH')
    #res$zscore <- (res$mean_diff - rowMeans2(mean_diff_rnd)) / rowSds(mean_diff_rnd)
    sds <- sqrt(rowSums(mean_diff_rnd^2)/(R-1))
    res$zscore <- (res$mean_diff - rowMeans2(mean_diff_rnd)) / sds
    res$pval <- 2 * pnorm(-abs(res$zscore))
    res$pval_adj <- p.adjust(res$pval, method = 'BH')
    
    if (length(comparisons) > 1) {
      rownames(res) <- NULL
    }
    return(res)
  })
  res <- Reduce(rbind, res_lst)
  if (length(compare) == 1 && compare == 'each_vs_rest' && !is.null(res)) {
    res$group1 <- factor(res$group1, levels = group_levels)
    res$group2 <- factor(res$group2)
  } 
  if (length(compare) == 1 && compare == 'all_vs_all' && !is.null(res)) {
    res$group1 <- factor(res$group1, levels = group_levels)
    res$group2 <- factor(res$group2, levels = group_levels)
  }
  return(res)
}

# helper functions

combine_counts <- function(group_counts, columns) {
 as.matrix(rowSums(group_counts[, columns, drop = FALSE]))
}

# combine per-group-mean to get the mean spanning multiple groups
# in an act of irrational premature optimization, we pass the 
# log-space mean when mean_type is geometric - need to make sure to 
# transform with exp1m before returning
combine_means <- function(means, n_items, columns, mean_type) {
  if (length(columns) == 1) {
    if (mean_type == 'arithmetic') {
      return(means[, columns, drop = FALSE])
    }
    if (mean_type == 'geometric') {
      return(expm1(means[, columns, drop = FALSE]))
    }
  }
  means <- means[, columns]
  n_items <- n_items[columns]
  tmp <- sweep(x = means, MARGIN = 2, STATS = n_items, FUN = '*')
  
  if (mean_type == 'arithmetic') {
    return(as.matrix(rowSums(tmp) / sum(n_items)))
  }
  if (mean_type == 'geometric') {
    return(as.matrix(expm1(rowSums(tmp) / sum(n_items))))
  }
}

#' Find differentially expressed genes that are conserved across samples
#'
#' @param y A matrix of counts; must be (or inherit from) class dgCMatrix; genes are rows,
#' cells are columns
#' @param group_labels The group labels (i.e. clusters or time points); 
#' will be converted to factor
#' @param sample_labels The sample labels; will be converted to factor
#' @param balanced Boolean, see details for explanation; default is TRUE
#' @param compare Specifies which groups to compare, see details; currently only 'each_vs_rest' 
#' (the default) is supported
#' @param pval_th P-value threshold used to call a gene differentially expressed when summarizing 
#' the tests per gene
#' @param ... Parameters passed to diff_mean_test
#' 
#' @return Data frame of results
#' 
#' @section Details:
#' This function calls diff_mean_test repeatedly and aggregates the results per group and gene.
#' 
#' If balanced is TRUE (the default), it is assumed that each sample spans multiple groups, 
#' as would be the case when merging or integrating samples from the same tissue followed by 
#' clustering. Here the group labels would be the clusters and cluster markers would have support
#' in each sample.
#' 
#' If balanced is FALSE, an unbalanced design is assumed where each sample contributes to one
#' group. An example is a time series experiment where some samples are taken from time point 
#' 1 while other samples are taken from time point 2. The time point would be the group label
#' and the goal would be to identify differentially expressed genes between time points that
#' are supported by many between-sample comparisons.
#' 
#' Output columns:
#' \describe{
#' \item{group1}{Group label of the frist group of cells}
#' \item{group2}{Group label of the second group of cells; currently fixed to 'rest'}
#' \item{gene}{Gene name (from rownames of input matrix)}
#' \item{n_tests}{The number of tests this gene participated in for this group}
#' \item{log2FC_min,median,max}{Summary statistics for log2FC across the tests}
#' \item{mean1,2_median}{Median of group mean across the tests}
#' \item{pval_max}{Maximum of p-values across tests}
#' \item{de_tests}{Number of tests that showed this gene having a log2FC going in the same
#' direction as log2FC_median and having a p-value <= pval_th}
#' }
#' 
#' The output is ordered by group1, -de_tests, -abs(log2FC_median), pval_max
#' 
#' @import Matrix
#' @importFrom dplyr n group_by summarise arrange
#' @importFrom rlang .data
#' 
#' @export
#'
#' @examples
#' \donttest{
#' clustering <- 1:ncol(pbmc) %% 2
#' sample_id <- 1:ncol(pbmc) %% 3
#' vst_out <- vst(pbmc, return_corrected_umi = TRUE)
#' de_res <- diff_mean_test_conserved(y = vst_out$umi_corrected, 
#' group_labels = clustering, sample_labels = sample_id)
#' }
#'
diff_mean_test_conserved <- function(y, group_labels, sample_labels, balanced = TRUE, 
                                     compare = 'each_vs_rest', pval_th = 1e-4, ...) {
  if (!inherits(x = y, what = 'dgCMatrix')) {
    stop('y must be a dgCMatrix')
  }
  group_labels <- droplevels(as.factor(group_labels))
  sample_labels <- droplevels(as.factor(sample_labels))
  
  res <- NULL
  if (compare[1] == 'each_vs_rest') {
    if (balanced) {
      res_lst <- lapply(levels(sample_labels), function(sl) {
        sel <- sample_labels == sl
        res <- diff_mean_test(y = y[, sel], group_labels = group_labels[sel], 
                              compare = compare, ...)
        if (!is.null(res)) {
          res$sample <- sl
        }
        res
      })
    } else {
      # fix special case when there are only two groups
      if (length(levels(group_labels)) == 2) {
        group_labels_to_do <- levels(group_labels)[1]
        gl_rest <- levels(group_labels)[2]
      } else {
        group_labels_to_do <- levels(group_labels)
        gl_rest <- 'rest'
      }
      # for each group, compare each sample against all samples that are not in group individually
      res_lst <- lapply(group_labels_to_do, function(gl) {
        gl_sel <- group_labels == gl
        samples_in_group <- sample_labels[gl_sel]
        res_lst <- lapply(unique(samples_in_group), function(sl_in_group) {
          sl_in_group_sel <- sample_labels == sl_in_group
          other_samples_not_in_group <- sample_labels[!(gl_sel | sl_in_group_sel)]
          res_lst <- lapply(unique(other_samples_not_in_group), function(osl_not_in_group) {
            sel <- (gl_sel & sl_in_group_sel) | (!gl_sel & sample_labels == osl_not_in_group)
            tmp_group <- c(gl, gl_rest)[as.numeric(!gl_sel & sample_labels == osl_not_in_group) + 1]
            res <- diff_mean_test(y = y[, sel], group_labels = tmp_group[sel], 
                                  compare = c(gl, gl_rest), ...)
            if (!is.null(res)) {
              res$sample1 <- sl_in_group
              res$sample2 <- osl_not_in_group
            }
            res
          })
          do.call(rbind, res_lst)
        })
        do.call(rbind, res_lst)
      })
    }
    res <- do.call(rbind, res_lst)
    levels(res$group1) <- levels(group_labels)
  }
  if (!is.null(res)) {
    res <- group_by(res, .data$group1, .data$group2, .data$gene) %>%
      summarise(n_tests = n(),
                log2FC_min = min(.data$log2FC), 
                log2FC_median = median(.data$log2FC),
                log2FC_max = max(.data$log2FC),
                mean1_median = median(.data$mean1),
                mean2_median = median(.data$mean2),
                pval_max = max(.data$pval),
                de_tests = sum((sign(.data$log2FC) == sign(.data$log2FC_median)) &
                               (.data$pval <= pval_th)),
                .groups = 'drop') %>%
      arrange(.data$group1, -.data$de_tests, -abs(.data$log2FC_median), .data$pval_max)
  }
  return(res)
}



# non-parametric differential expression test
np_de_test <- function(y, labels, N = 100, S = 100, randomize = FALSE) {
  if (!inherits(x = y, what = 'matrix')) {
    stop('y must be a matrix')
  }
  labels <- droplevels(as.factor(labels))
  if (length(levels(labels)) > 2) {
    stop('only two groups can be compared')
  }
  if (N < 50 || S < 50) {
    stop('N and S must both be at least 50')
  }
  if (ncol(y) != length(labels)) {
    stop('number of columns in y and length of label vector must match')
  }
  labels <- as.integer(labels)-1L
  if (randomize) {
    res <- apply(y, 1, function(x) distribution_shift(mean_boot_grouped(x, sample(labels), N = 100, S = 100)))
  } else {
    res <- apply(y, 1, function(x) distribution_shift(mean_boot_grouped(x, labels, N = 100, S = 100)))
  }
  #res <- data.frame(t(res))
  #colnames(res) <- c('q16_a', 'q50_a', 'q84_a', 'q16_b', 'q50_b', 'q84_b', 'div', 'z')
  res <- data.frame(t(res[c(2, 5, 7, 8), ]))
  colnames(res) <- c('mu1', 'mu2', 'div', 'z')
  gene <- rownames(res)
  res <- cbind(gene, res)
  rownames(res) <- NULL
  return(res)
}

