
# Check cell attributes; add missing ones
make_cell_attr <- function(umi, cell_attr, latent_var, batch_var, latent_var_nonreg, verbosity) {
  if (is.null(cell_attr)) {
    cell_attr <- data.frame(row.names = colnames(umi))
  }

  # Make sure count matrix has row and column names
  if (is.null(rownames(umi)) || is.null(colnames(umi))) {
    stop('count matrix must have row and column names')
  }

  # Make sure rownames of cell attributes match cell names in count matrix
  if (!identical(rownames(cell_attr), colnames(umi))) {
    stop('cell attribute row names must match column names of count matrix')
  }

  # Do not allow certain variable names
  no_good <- c('(Intercept)', 'Intercept')
  if (any(no_good %in% c(latent_var, batch_var, latent_var_nonreg))) {
    stop('Do not use the following variable names for a latent variable or batch variable: ', paste(no_good, collapse = ', '))
  }

  # these are the cell attributes that we know how to calculate given the count matrix
  known_attr <- c('umi', 'gene', 'log_umi', 'log_gene', 'umi_per_gene', 'log_umi_per_gene')
  # these are the missing cell attributes specified in latent_var
  missing_attr <- setdiff(c(latent_var, batch_var, latent_var_nonreg), colnames(cell_attr))
  if (length(missing_attr) > 0) {
    if (verbosity > 0) {
      message('Calculating cell attributes from input UMI matrix: ', paste(missing_attr, collapse = ', '))
    }
    unknown_attr <- setdiff(missing_attr, known_attr)
    if (length(unknown_attr) > 0) {
      stop(sprintf('Unknown cell attributes: %s. Check latent_var, batch_var and latent_var_nonreg and make sure the variables are in cell_attr', paste(unknown_attr, collapse = ', ')))
    }
    new_attr <- list()
    if (any(c('umi', 'log_umi', 'umi_per_gene', 'log_umi_per_gene') %in% missing_attr)) {
      new_attr$umi <- colSums(umi)
      new_attr$log_umi <- log10(new_attr$umi)
    }
    if (any(c('gene', 'log_gene', 'umi_per_gene', 'log_umi_per_gene') %in% missing_attr)) {
      new_attr$gene <- colSums(umi > 0)
      new_attr$log_gene <- log10(new_attr$gene)
    }
    if (any(c('umi_per_gene', 'log_umi_per_gene') %in% missing_attr)) {
      new_attr$umi_per_gene <- new_attr$umi / new_attr$gene
      new_attr$log_umi_per_gene <- log10(new_attr$umi_per_gene)
    }
    new_attr <- do.call(cbind, new_attr)
    cell_attr <- cbind(cell_attr, new_attr[, setdiff(colnames(new_attr), colnames(cell_attr)), drop = FALSE])
  }

  # make sure no NA, NaN, Inf values are in cell attributes - they would cause
  # problems later on
  for (ca in c(latent_var, batch_var, latent_var_nonreg)) {
    ca_values <- cell_attr[, ca]
    if (any(is.na(ca_values)) ||
        any(is.nan(ca_values)) ||
        any(is.infinite(ca_values))) {
      stop('cell attribute "', ca, '" contains NA, NaN, or infinite value')
    }
  }

  return(cell_attr)
}

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
    ret <- row_gmean_dgcmatrix(matrix = x, eps = eps)
    names(ret) <- rownames(x)
    return(ret)
  }
  stop('matrix x needs to be of class matrix or dgCMatrix')
}



#' Variance per row
#'
#' @param x matrix of class \code{matrix} or \code{dgCMatrix}
#'
#' @return variances
#'
#' @importFrom matrixStats rowVars
row_var <- function(x) {
  if (inherits(x = x, what = 'matrix')) {
    ret <- rowVars(x)
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


pearson_residual2 <- function(y, mu, theta, min_vars) {
  model_var <- mu + mu^2 / theta
  for (row in 1:nrow(model_var)){
    var_row <- model_var[row,]
    min_var <- min_vars[row]
    var_row[var_row < min_var] <- min_var
    model_var[row,] <- var_row
  }
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
#' @param verbosity An integer specifying whether to show only messages (1), messages and progress bars (2) or nothing (0) while the function is running; default is 2
#' @param verbose Deprecated; use verbosity instead
#' @param show_progress Deprecated; use verbosity instead
#'
#' @return A matrix of residuals
#'
#' @export
#'
#' @examples
#' \donttest{
#' vst_out <- vst(pbmc, return_cell_attr = TRUE)
#' pearson_res <- get_residuals(vst_out, pbmc)
#' deviance_res <- get_residuals(vst_out, pbmc, residual_type = 'deviance')
#' }
#'
get_residuals <- function(vst_out, umi, residual_type = 'pearson',
                          res_clip_range = c(-sqrt(ncol(umi)), sqrt(ncol(umi))),
                          min_variance = vst_out$arguments$min_variance,
                          cell_attr = vst_out$cell_attr, bin_size = 256,
                          verbosity = vst_out$arguments$verbosity,
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

  # min_variance estimated using median umi
  if (min_variance == "umi_median"){
    # Maximum pearson residual for non-zero median UMI is 5
    min_var <- (get_nz_median(umi) / 5)^2
    if (verbosity > 0) {
      message(paste("Setting min_variance based on median UMI: ", min_var))
    }
  } else {
    if (verbosity > 0) {
      message(paste("Setting min_variance to: ", min_variance))
    }
    min_var <- min_variance
  }
  regressor_data <- model.matrix(as.formula(gsub('^y', '', vst_out$model_str)), cell_attr)
  model_pars <- vst_out$model_pars_fit
  if (!is.null(dim(vst_out$model_pars_nonreg))) {
    regressor_data_nonreg <- model.matrix(as.formula(gsub('^y', '', vst_out$model_str_nonreg)), cell_attr)
    regressor_data <- cbind(regressor_data, regressor_data_nonreg)
    model_pars <- cbind(vst_out$model_pars_fit, vst_out$model_pars_nonreg)
  }

  genes <- rownames(umi)[rownames(umi) %in% rownames(model_pars)]
  if (verbosity > 0) {
    message('Calculating residuals of type ', residual_type, ' for ', length(genes), ' genes')
  }
  bin_ind <- ceiling(x = 1:length(x = genes) / bin_size)
  max_bin <- max(bin_ind)
  if (verbosity > 1) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  res <- matrix(NA_real_, length(genes), nrow(regressor_data), dimnames = list(genes, rownames(regressor_data)))
  for (i in 1:max_bin) {
    genes_bin <- genes[bin_ind == i]
    mu <- exp(tcrossprod(model_pars[genes_bin, -1, drop=FALSE], regressor_data))

    y <- as.matrix(umi[genes_bin, , drop=FALSE])
    if (min_variance == "model_mean") {
      mu_mean_var <- matrixStats::rowMeans2(mu)
      res[genes_bin, ] <- switch(residual_type,
                                 'pearson' = pearson_residual2(y, mu, model_pars[genes_bin, 'theta'], min_vars = mu_mean_var),
                                 'deviance' = deviance_residual(y, mu, model_pars[genes_bin, 'theta']),
                                 stop('residual_type ', residual_type, ' unknown - only pearson and deviance supported at the moment'))
    } else if (min_variance == "model_median") {
        mu_median_var <- matrixStats::rowMedians(mu)
        res[genes_bin, ] <- switch(residual_type,
                                   'pearson' = pearson_residual2(y, mu, model_pars[genes_bin, 'theta'], min_vars = mu_median_var),
                                   'deviance' = deviance_residual(y, mu, model_pars[genes_bin, 'theta']),
                                   stop('residual_type ', residual_type, ' unknown - only pearson and deviance supported at the moment'))

      } else {
        res[genes_bin, ] <- switch(residual_type,
                                   'pearson' = pearson_residual(y, mu, model_pars[genes_bin, 'theta'], min_var = min_var),
                                   'deviance' = deviance_residual(y, mu, model_pars[genes_bin, 'theta']),
                                   stop('residual_type ', residual_type, ' unknown - only pearson and deviance supported at the moment'))
    }

    if (verbosity > 1) {
      setTxtProgressBar(pb, i)
    }
  }
  if (verbosity > 1) {
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
#' @param verbosity An integer specifying whether to show only messages (1), messages and progress bars (2) or nothing (0) while the function is running; default is 2
#' @param verbose Deprecated; use verbosity instead
#' @param show_progress Deprecated; use verbosity instead
#'
#' @return A vector of residual variances (after clipping)
#'
#' @export
#'
#' @examples
#' \donttest{
#' vst_out <- vst(pbmc, return_cell_attr = TRUE)
#' res_var <- get_residual_var(vst_out, pbmc)
#' }
#'
get_residual_var <- function(vst_out, umi, residual_type = 'pearson',
                             res_clip_range = c(-sqrt(ncol(umi)), sqrt(ncol(umi))),
                             min_variance = vst_out$arguments$min_variance,
                             cell_attr = vst_out$cell_attr, bin_size = 256,
                             verbosity = vst_out$arguments$verbosity,
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

  regressor_data <- model.matrix(as.formula(gsub('^y', '', vst_out$model_str)), cell_attr)
  model_pars <- vst_out$model_pars_fit
  if (!is.null(dim(vst_out$model_pars_nonreg))) {
    regressor_data_nonreg <- model.matrix(as.formula(gsub('^y', '', vst_out$model_str_nonreg)), cell_attr)
    regressor_data <- cbind(regressor_data, regressor_data_nonreg)
    model_pars <- cbind(vst_out$model_pars_fit, vst_out$model_pars_nonreg)
  }

  genes <- rownames(umi)[rownames(umi) %in% rownames(model_pars)]
  # min_variance estimated using median umi
  if (min_variance == "umi_median"){
    # Maximum pearson residual for non-zero median UMI is 5
    min_var <- (get_nz_median(umi, genes) / 5)^2
    if (verbosity > 0) {
      message(paste("Setting min_variance based on median UMI: ", min_var))
    }
  } else {
    if (verbosity > 0) {
      message(paste("Setting min_variance to: ", min_variance))
    }
    min_var <- min_variance
  }
  if (verbosity > 0) {
    message('Calculating variance for residuals of type ', residual_type, ' for ', length(genes), ' genes')
  }
  bin_ind <- ceiling(x = 1:length(x = genes) / bin_size)
  max_bin <- max(bin_ind)
  if (verbosity > 1) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  res <- matrix(NA_real_, length(genes))
  names(res) <- genes
  for (i in 1:max_bin) {
    genes_bin <- genes[bin_ind == i]
    mu <- exp(tcrossprod(model_pars[genes_bin, -1, drop=FALSE], regressor_data))
    y <- as.matrix(umi[genes_bin, , drop=FALSE])
    if (min_variance == "model_mean") {
      mu_mean_var <- matrixStats::rowMeans2(mu)
      res_mat <- switch(residual_type,
                        'pearson' = pearson_residual2(y, mu, model_pars[genes_bin, 'theta'], min_vars = mu_mean_var),
                        'deviance' = deviance_residual(y, mu, model_pars[genes_bin, 'theta']),
                        stop('residual_type ', residual_type, ' unknown - only pearson and deviance supported at the moment'))
    } else if (min_variance == "model_median") {
      mu_median_var <- matrixStats::rowMedians(mu)
      res_mat <- switch(residual_type,
                        'pearson' = pearson_residual2(y, mu, model_pars[genes_bin, 'theta'], min_vars = mu_median_var),
                        'deviance' = deviance_residual(y, mu, model_pars[genes_bin, 'theta']),
                        stop('residual_type ', residual_type, ' unknown - only pearson and deviance supported at the moment'))

    } else {
      res_mat <- switch(residual_type,
                        'pearson' = pearson_residual(y, mu, model_pars[genes_bin, 'theta'], min_var = min_var),
                        'deviance' = deviance_residual(y, mu, model_pars[genes_bin, 'theta']),
                        stop('residual_type ', residual_type, ' unknown - only pearson and deviance supported at the moment'))
      }
    res_mat[res_mat < res_clip_range[1]] <- res_clip_range[1]
    res_mat[res_mat > res_clip_range[2]] <- res_clip_range[2]
    res[genes_bin] <- row_var(res_mat)
    if (verbosity > 1) {
      setTxtProgressBar(pb, i)
    }
  }
  if (verbosity > 1) {
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
#' @param verbosity An integer specifying whether to show only messages (1), messages and progress bars (2) or nothing (0) while the function is running; default is 2
#' @param verbose Deprecated; use verbosity instead
#' @param show_progress Deprecated; use verbosity instead
#'
#' @return A named vector of variances (the average across all cells), one entry per gene.
#'
#' @export
#'
#' @examples
#' \donttest{
#' vst_out <- vst(pbmc, return_cell_attr = TRUE)
#' res_var <- get_model_var(vst_out)
#' }
#'
get_model_var <- function(vst_out, cell_attr = vst_out$cell_attr, use_nonreg = FALSE,
                          bin_size = 256, verbosity = 2,
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
  if (verbosity > 0) {
    message('Calculating model variance for ', length(genes), ' genes')
  }
  bin_ind <- ceiling(x = 1:length(x = genes) / bin_size)
  max_bin <- max(bin_ind)
  if (verbosity > 1) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  res <- matrix(NA_real_, length(genes))
  names(res) <- genes
  for (i in 1:max_bin) {
    genes_bin <- genes[bin_ind == i]
    mu <- exp(tcrossprod(model_pars[genes_bin, -1, drop=FALSE], regressor_data))
    model_var = mu + mu^2 / model_pars[genes_bin, 'theta']
    res[genes_bin] <- rowMeans(model_var)
    if (verbosity > 1) {
      setTxtProgressBar(pb, i)
    }
  }
  if (verbosity > 1) {
    close(pb)
  }
  return(res)
}


#' Get median of non zero UMIs from a count matrix using a subset of genes (slow)
#'
#' @param umi Count matrix
#' @param genes List of genes to calculate statistics. Default is NULL which returns the non-zero median using all genes
#'
#' @return A numeric value representing the median of non-zero entries from the UMI matrix
get_nz_median <- function(umi, genes = NULL){
  cm.T <- Matrix::t(umi)
  n_g <- dim(umi)[1]
  allnonzero <- c()
  if (is.null(genes)) {
    gene_index <- seq(1, nrow(umi))
  } else {
    gene_index <- which(genes %in% rownames(umi))
  }
  for (g in gene_index) {
    m_i <- cm.T@x[(cm.T@p[g] + 1):cm.T@p[g + 1]]
    allnonzero <- c(allnonzero, m_i)
  }
  return (median(allnonzero, na.rm = TRUE))
}

#' Get median of non zero UMIs from a count matrix
#'
#' @param umi Count matrix
#'
#' @return A numeric value representing the median of non-zero entries from the UMI matrix
get_nz_median2 <- function(umi){
  return (median(umi@x))
}

