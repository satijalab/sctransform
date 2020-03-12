
#' Smooth data by PCA
#'
#' Perform PCA, identify significant dimensions, and reverse the rotation using only significant dimensions.
#'
#' @param x A data matrix with genes as rows and cells as columns
#' @param elbow_th The fraction of PC sdev drop that is considered significant; low values will lead to more PCs being used
#' @param dims_use Directly specify PCs to use, e.g. 1:10
#' @param max_pc Maximum number of PCs computed
#' @param do_plot Plot PC sdev and sdev drop
#' @param scale. Boolean indicating whether genes should be divided by standard deviation after centering and prior to PCA
#'
#' @return Smoothed data
#'
#' @importFrom graphics par plot abline
#'
#' @export
#'
#' @examples
#' \donttest{
#' vst_out <- vst(pbmc)
#' y_smooth <- smooth_via_pca(vst_out$y, do_plot = TRUE)
#' }
#'
smooth_via_pca <- function(x, elbow_th = 0.025, dims_use = NULL, max_pc = 100, do_plot = FALSE,
                           scale. = FALSE) {
  requireNamespace('irlba', quietly = TRUE)
  # perform pca
  if (scale.) {
    scale. <- apply(x, 1, 'sd')
  } else {
    scale. <- rep(1, nrow(x))
  }
  pca <- irlba::prcomp_irlba(t(x), n = max_pc, center = TRUE, scale. = scale.)

  if (is.null(dims_use)) {
    pca_sdev_drop <- c(diff(pca$sdev), 0) / -pca$sdev
    max_dim <- rev(which(pca_sdev_drop > elbow_th))[1]
    dims_use <- 1:max_dim

    if (do_plot) {
      par(mfrow=c(1,2))
      plot(pca$sdev)
      abline(v = max_dim + 0.5, col='red')
      plot(pca_sdev_drop)
      abline(h = elbow_th, col='red')
      abline(v = max_dim + 0.5, col='red')
      par(mfrow=c(1,1))
    }
  }

  new_x <- pca$rotation[, dims_use] %*% t(pca$x[, dims_use]) * pca$scale + pca$center
  dimnames(new_x) <- dimnames(x)
  return(new_x)
}


#' Correct data by setting all latent factors to their median values and reversing the regression model
#'
#' @param x A list that provides model parameters and optionally meta data; use output of vst function
#' @param data The name of the entry in x that holds the data
#' @param cell_attr Provide cell meta data holding latent data info
#' @param do_round Round the result to integers
#' @param do_pos Set negative values in the result to zero
#' @param verbose Whether to show messages; default is TRUE
#' @param show_progress Whether to print progress bar; default is same as verbose
#'
#' @return Corrected data as UMI counts
#'
#' @export
#'
#' @examples
#' \donttest{
#' vst_out <- vst(pbmc, return_cell_attr = TRUE)
#' umi_corrected <- correct(vst_out)
#' }
#'
correct <- function(x, data = 'y', cell_attr = x$cell_attr, do_round = TRUE,
                    do_pos = TRUE, verbose = TRUE, show_progress = verbose) {
  if (is.character(data)) {
    data <- x[[data]]
  }
  # when correcting, set all latent variables to median values
  cell_attr[, x$arguments$latent_var] <- apply(cell_attr[, x$arguments$latent_var, drop=FALSE], 2, function(x) rep(median(x), length(x)))
  regressor_data <- model.matrix(as.formula(gsub('^y', '', x$model_str)), cell_attr)

  genes <- rownames(data)
  bin_size <- x$arguments$bin_size
  bin_ind <- ceiling(x = 1:length(x = genes) / bin_size)
  max_bin <- max(bin_ind)
  if (verbose) {
    message('Computing corrected count matrix for ', length(genes), ' genes')
  }
  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  corrected_data <- matrix(NA_real_, length(genes), nrow(regressor_data), dimnames = list(genes, rownames(regressor_data)))
  for (i in 1:max_bin) {
    genes_bin <- genes[bin_ind == i]
    pearson_residual <- data[genes_bin, ]
    coefs <- x$model_pars_fit[genes_bin, -1]
    theta <- x$model_pars_fit[genes_bin, 1]
    mu <- exp(tcrossprod(coefs, regressor_data))
    variance <- mu + mu^2 / theta
    corrected_data[genes_bin, ] <- mu + pearson_residual * sqrt(variance)
    if (show_progress) {
      setTxtProgressBar(pb, i)
    }
  }
  if (show_progress) {
    close(pb)
  }

  if (do_round) {
    corrected_data <- round(corrected_data, 0)
  }
  if (do_pos) {
    corrected_data[corrected_data < 0] <- 0
  }
  return(corrected_data)
}

#' Correct data by setting all latent factors to their median values and reversing the regression model
#'
#' This version does not need a matrix of Pearson residuals. It takes the count matrix as input and
#' calculates the residuals on the fly. The corrected UMI counts will be rounded to the nearest
#' integer and negative values clipped to 0.
#'
#' @param x A list that provides model parameters and optionally meta data; use output of vst function
#' @param umi The count matrix
#' @param cell_attr Provide cell meta data holding latent data info
#' @param verbose Whether to show messages; default is TRUE
#' @param show_progress Whether to print progress bar; default is same as verbose
#'
#' @return Corrected data as UMI counts
#'
#' @importFrom methods as
#'
#' @export
#'
#' @examples
#' \donttest{
#' vst_out <- vst(pbmc, return_cell_attr = TRUE)
#' umi_corrected <- correct_counts(vst_out, pbmc)
#' }
#'
correct_counts <- function(x, umi, cell_attr = x$cell_attr, verbose = TRUE,
                           show_progress = verbose) {
  regressor_data_orig <- model.matrix(as.formula(gsub('^y', '', x$model_str)), cell_attr)
  # when correcting, set all latent variables to median values
  cell_attr[, x$arguments$latent_var] <- apply(cell_attr[, x$arguments$latent_var, drop=FALSE], 2, function(x) rep(median(x), length(x)))
  regressor_data <- model.matrix(as.formula(gsub('^y', '', x$model_str)), cell_attr)

  genes <- rownames(umi)[rownames(umi) %in% rownames(x$model_pars_fit)]
  bin_size <- x$arguments$bin_size
  bin_ind <- ceiling(x = 1:length(x = genes) / bin_size)
  max_bin <- max(bin_ind)
  if (verbose) {
    message('Computing corrected UMI count matrix')
  }
  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  #corrected_data <- matrix(NA_real_, length(genes), nrow(regressor_data), dimnames = list(genes, rownames(regressor_data)))
  corrected_data <- list()
  for (i in 1:max_bin) {
    genes_bin <- genes[bin_ind == i]
    coefs <- x$model_pars_fit[genes_bin, -1, drop=FALSE]
    theta <- x$model_pars_fit[genes_bin, 1]
    # get pearson residuals
    mu <- exp(tcrossprod(coefs, regressor_data_orig))
    variance <- mu + mu^2 / theta
    y <- as.matrix(umi[genes_bin, , drop=FALSE])
    pearson_residual <- (y - mu) / sqrt(variance)
    # generate output
    mu <- exp(tcrossprod(coefs, regressor_data))
    variance <- mu + mu^2 / theta
    y.res <- mu + pearson_residual * sqrt(variance)
    y.res <- round(y.res, 0)
    y.res[y.res < 0] <- 0
    corrected_data[[length(corrected_data) + 1]] <- as(y.res, Class = 'dgCMatrix')
    if (show_progress) {
      setTxtProgressBar(pb, i)
    }
  }
  if (show_progress) {
    close(pb)
  }
  corrected_data <- do.call(what = rbind, args = corrected_data)

  return(corrected_data)
}

reverse_regression <- function(pearson_residual, theta, coefs, data) {
  mu <- exp(data %*% coefs)[, 1]
  variance <- mu + mu^2 / theta
  return(mu + pearson_residual * sqrt(variance))
}
