
#' Geometric mean per row
#'
#' @param x matrix of class \code{matrix} or \code{dgCMatrix}
#' @param eps small value to add to x to avoid log(0); default is 1
#'
#' @return geometric means
row_gmean <- function(x, eps = 1) {
  if (inherits(x = x, what = 'matrix')) {
    return(exp(rowMeans(log(x + eps))) - eps)
  }
  if (inherits(x = x, what = 'dgCMatrix')) {
    ret <- row_gmean_dgcmatrix(x = x@x, i = x@i, rows = nrow(x), cols = ncol(x), eps = eps)
    names(ret) <- rownames(x)
    return(ret)
  }
  stop('matrix x needs to be of class matrix or dgCMatrix')
}

#' Geometric mean per row grouped by a factor
#'
#' @param x matrix of class \code{dgCMatrix}
#' @param group factor to group the columns by (will be converted using \code{as.factor} and \code{droplevels})
#' @param eps small value to add to x to avoid log(0); default is 1
#'
#' @return matrix of geometric means
row_gmean_grouped <- function(x, group, eps = 1) {
  group <- droplevels(as.factor(group))
  if (inherits(x = x, what = 'dgCMatrix')) {
    ret <- row_gmean_grouped_dgcmatrix(x = x@x, i = x@i, p = x@p,
                                       group = as.integer(group) - 1,
                                       groups = length(levels(group)),
                                       rows = nrow(x),
                                       eps = eps)
    rownames(ret) <- rownames(x)
    colnames(ret) <- levels(group)
    return(ret)
  }
  stop('matrix x needs to be of class dgCMatrix')
}

#' Arithmetic mean per row grouped by a factor
#'
#' @param x matrix of class \code{dgCMatrix}
#' @param group factor to group the columns by (will be converted using \code{as.factor} and \code{droplevels})
#'
#' @return matrix of arithmetic means
row_mean_grouped <- function(x, group) {
  group <- droplevels(as.factor(group))
  if (inherits(x = x, what = 'dgCMatrix')) {
    ret <- row_mean_grouped_dgcmatrix(x = x@x, i = x@i, p = x@p,
                                      group = as.integer(group) - 1,
                                      groups = length(levels(group)),
                                      rows = nrow(x))
    rownames(ret) <- rownames(x)
    colnames(ret) <- levels(group)
    return(ret)
  }
  stop('matrix x needs to be of class dgCMatrix')
}

#' Variance per row
#'
#' @param x matrix of class \code{matrix} or \code{dgCMatrix}
#'
#' @return variances
row_var <- function(x) {
  if (inherits(x = x, what = 'matrix')) {
    ret <- switch(storage.mode(x),
                  'double' = row_var_dense_d(x),
                  'integer' = row_var_dense_i(x),
                  stop('Unknown matrix storage mode'))
    names(ret) <- rownames(x)
    return(ret)
  }
  if (inherits(x = x, what = 'dgCMatrix')) {
    ret <- row_var_dgcmatrix(x = x@x, i = x@i, rows = nrow(x), cols = ncol(x))
    names(ret) <- rownames(x)
    return(ret)
  }
  stop('matrix x needs to be of class matrix or dgCMatrix')
}



#' Identify outliers
#'
#' @param y Dependent variable
#' @param x Independent variable
#' @param th Outlier score threshold
#'
#' @return Boolean vector
#'
#' @importFrom stats aggregate
#'
is_outlier <- function(y, x, th = 10) {
  #bin.width <- var(x) * bw.SJ(x)
  bin.width <- (max(x) - min(x)) * bw.SJ(x) / 2
  eps <- .Machine$double.eps * 10
  breaks1 <- seq(from = min(x) - eps, to = max(x) + bin.width, by = bin.width)
  breaks2 <- seq(from = min(x) - eps - bin.width/2, to = max(x) + bin.width, by = bin.width)
  score1 <- robust_scale_binned(y, x, breaks1)
  score2 <- robust_scale_binned(y, x, breaks2)
  return(pmin(abs(score1), abs(score2)) > th)
}

#' Robust scale using median and mad per bin
#'
#' @param y Numeric vector
#' @param x Numeric vector
#' @param breaks Numeric vector of breaks
#'
#' @return Numeric vector of scaled score
#'
#' @importFrom stats aggregate
#'
robust_scale_binned <- function(y, x, breaks) {
  bins <- cut(x = x, breaks = breaks, ordered_result = TRUE)
  tmp <- aggregate(x = y, by = list(bin=bins), FUN = robust_scale)
  score <- rep(0, length(x))
  o <- order(bins)
  if (inherits(x = tmp$x, what = 'list')) {
    score[o] <- unlist(tmp$x)
  } else {
    score[o] <- as.numeric(t(tmp$x))
  }
  return(score)
}

#' Robust scale using median and mad
#'
#' @param x Numeric
#'
#' @return Numeric
#'
#' @importFrom stats median mad
#'
robust_scale <- function(x) {
  return((x - median(x)) / (mad(x) + .Machine$double.eps))
}

pearson_residual <- function(y, mu, theta, min_var = -Inf) {
  model_var <- mu + mu^2 / theta
  model_var[model_var < min_var] <- min_var
  return((y - mu) / sqrt(model_var))
}

sq_deviance_residual <- function(y, mu, theta, wt=1) {
  2 * wt * (y * log(pmax(1, y)/mu) - (y + theta) * log((y + theta)/(mu + theta)))
}

deviance_residual <- function(y, mu, theta, wt=1) {
  r <- 2 * wt * (y * log(pmax(1, y)/mu) - (y + theta) * log((y + theta)/(mu + theta)))
  sqrt(r) * sign(y - mu)
}

#' Return Pearson or deviance residuals of regularized models
#'
#' @param vst_out The output of a vst run
#' @param umi The UMI count matrix that will be used
#' @param residual_type What type of residuals to return; can be 'pearson' or 'deviance'; default is 'pearson'
#' @param res_clip_range Numeric of length two specifying the min and max values the results will be clipped to; default is c(-sqrt(ncol(umi)), sqrt(ncol(umi)))
#' @param min_variance Lower bound for the estimated variance for any gene in any cell when calculating pearson residual; default is vst_out$arguments$min_variance
#' @param cell_attr Data frame of cell meta data
#' @param bin_size Number of genes to put in each bin (to show progress)
#' @param verbose Whether to show messages; default is TRUE
#' @param show_progress Whether to print progress bar; default is verbose
#'
#' @return A matrix of residuals
#'
#' @export
#'
#' @examples
#' \dontrun{
#' vst_out <- vst(pbmc)
#' pearson_res <- get_residuals(vst_out, pbmc)
#' deviance_res <- get_residuals(vst_out, pbmc, residual_type = 'deviance')
#' }
#'
get_residuals <- function(vst_out, umi, residual_type = 'pearson',
                          res_clip_range = c(-sqrt(ncol(umi)), sqrt(ncol(umi))),
                          min_variance = vst_out$arguments$min_variance,
                          cell_attr = vst_out$cell_attr, bin_size = 256,
                          verbose = TRUE, show_progress = verbose) {
  regressor_data <- model.matrix(as.formula(gsub('^y', '', vst_out$model_str)), cell_attr)
  model_pars <- vst_out$model_pars_fit
  if (!is.null(dim(vst_out$model_pars_nonreg))) {
    regressor_data_nonreg <- model.matrix(as.formula(gsub('^y', '', vst_out$model_str_nonreg)), cell_attr)
    regressor_data <- cbind(regressor_data, regressor_data_nonreg)
    model_pars <- cbind(vst_out$model_pars_fit, vst_out$model_pars_nonreg)
  }

  genes <- rownames(umi)[rownames(umi) %in% rownames(model_pars)]
  if (verbose) {
    message('Calculating residuals of type ', residual_type, ' for ', length(genes), ' genes')
  }
  bin_ind <- ceiling(x = 1:length(x = genes) / bin_size)
  max_bin <- max(bin_ind)
  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  res <- matrix(NA_real_, length(genes), nrow(regressor_data), dimnames = list(genes, rownames(regressor_data)))
  for (i in 1:max_bin) {
    genes_bin <- genes[bin_ind == i]
    mu <- exp(tcrossprod(model_pars[genes_bin, -1, drop=FALSE], regressor_data))
    y <- as.matrix(umi[genes_bin, , drop=FALSE])
    res[genes_bin, ] <- switch(residual_type,
      'pearson' = pearson_residual(y, mu, model_pars[genes_bin, 'theta'], min_var = min_variance),
      'deviance' = deviance_residual(y, mu, model_pars[genes_bin, 'theta'])
    )
    if (show_progress) {
      setTxtProgressBar(pb, i)
    }
  }
  if (show_progress) {
    close(pb)
  }
  res[res < res_clip_range[1]] <- res_clip_range[1]
  res[res > res_clip_range[2]] <- res_clip_range[2]
  return(res)
}

#' Return variance of residuals of regularized models
#'
#' This never creates the full residual matrix and can be used to determine highly variable genes.
#'
#' @param vst_out The output of a vst run
#' @param umi The UMI count matrix that will be used
#' @param residual_type What type of residuals to return; can be 'pearson' or 'deviance'; default is 'pearson'
#' @param res_clip_range Numeric of length two specifying the min and max values the residuals will be clipped to; default is c(-sqrt(ncol(umi)), sqrt(ncol(umi)))
#' @param min_variance Lower bound for the estimated variance for any gene in any cell when calculating pearson residual; default is vst_out$arguments$min_variance
#' @param cell_attr Data frame of cell meta data
#' @param bin_size Number of genes to put in each bin (to show progress)
#' @param verbose Whether to show messages; default is TRUE
#' @param show_progress Whether to print progress bar; default is same as verbose
#'
#' @return A vector of residual variances (after clipping)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' vst_out <- vst(pbmc)
#' res_var <- get_residual_var(vst_out, pbmc)
#' }
#'
get_residual_var <- function(vst_out, umi, residual_type = 'pearson',
                             res_clip_range = c(-sqrt(ncol(umi)), sqrt(ncol(umi))),
                             min_variance = vst_out$arguments$min_variance,
                             cell_attr = vst_out$cell_attr, bin_size = 256,
                             verbose = TRUE, show_progress = verbose) {
  regressor_data <- model.matrix(as.formula(gsub('^y', '', vst_out$model_str)), cell_attr)
  model_pars <- vst_out$model_pars_fit
  if (!is.null(dim(vst_out$model_pars_nonreg))) {
    regressor_data_nonreg <- model.matrix(as.formula(gsub('^y', '', vst_out$model_str_nonreg)), cell_attr)
    regressor_data <- cbind(regressor_data, regressor_data_nonreg)
    model_pars <- cbind(vst_out$model_pars_fit, vst_out$model_pars_nonreg)
  }

  genes <- rownames(umi)[rownames(umi) %in% rownames(model_pars)]
  if (verbose) {
    message('Calculating variance for residuals of type ', residual_type, ' for ', length(genes), ' genes')
  }
  bin_ind <- ceiling(x = 1:length(x = genes) / bin_size)
  max_bin <- max(bin_ind)
  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  res <- matrix(NA_real_, length(genes))
  names(res) <- genes
  for (i in 1:max_bin) {
    genes_bin <- genes[bin_ind == i]
    mu <- exp(tcrossprod(model_pars[genes_bin, -1, drop=FALSE], regressor_data))
    y <- as.matrix(umi[genes_bin, , drop=FALSE])
    res_mat <- switch(residual_type,
                      'pearson' = pearson_residual(y, mu, model_pars[genes_bin, 'theta'], min_var = min_variance),
                      'deviance' = deviance_residual(y, mu, model_pars[genes_bin, 'theta']))
    res_mat[res_mat < res_clip_range[1]] <- res_clip_range[1]
    res_mat[res_mat > res_clip_range[2]] <- res_clip_range[2]
    res[genes_bin] <- row_var(res_mat)
    if (show_progress) {
      setTxtProgressBar(pb, i)
    }
  }
  if (show_progress) {
    close(pb)
  }
  return(res)
}

#' Return average variance under negative binomial model
#'
#' This is based on the formula var = mu + mu^2 / theta
#'
#' @param vst_out The output of a vst run
#' @param cell_attr Data frame of cell meta data
#' @param use_nonreg Use the non-regularized parameter estimates; boolean; default is FALSE
#' @param bin_size Number of genes to put in each bin (to show progress)
#' @param verbose Whether to show messages; default is TRUE
#' @param show_progress Whether to print progress bar; default is same as verbose
#'
#' @return A named vector of variances (the average across all cells), one entry per gene.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' vst_out <- vst(pbmc)
#' res_var <- get_model_var(vst_out)
#' }
#'
get_model_var <- function(vst_out, cell_attr = vst_out$cell_attr, use_nonreg = FALSE,
                          bin_size = 256, verbose = TRUE, show_progress = verbose) {
  regressor_data <- model.matrix(as.formula(gsub('^y', '', vst_out$model_str)), cell_attr)
  if (use_nonreg) {
    model_pars <- vst_out$model_pars
  } else {
    model_pars <- vst_out$model_pars_fit
  }
  if (!is.null(dim(vst_out$model_pars_nonreg))) {
    regressor_data_nonreg <- model.matrix(as.formula(gsub('^y', '', vst_out$model_str_nonreg)), cell_attr)
    regressor_data <- cbind(regressor_data, regressor_data_nonreg)
    model_pars <- cbind(vst_out$model_pars_fit, vst_out$model_pars_nonreg)
  }

  genes <- rownames(model_pars)
  if (verbose) {
    message('Calculating model variance for ', length(genes), ' genes')
  }
  bin_ind <- ceiling(x = 1:length(x = genes) / bin_size)
  max_bin <- max(bin_ind)
  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  res <- matrix(NA_real_, length(genes))
  names(res) <- genes
  for (i in 1:max_bin) {
    genes_bin <- genes[bin_ind == i]
    mu <- exp(tcrossprod(model_pars[genes_bin, -1, drop=FALSE], regressor_data))
    model_var = mu + mu^2 / model_pars[genes_bin, 'theta']
    res[genes_bin] <- rowMeans(model_var)
    if (show_progress) {
      setTxtProgressBar(pb, i)
    }
  }
  if (show_progress) {
    close(pb)
  }
  return(res)
}
