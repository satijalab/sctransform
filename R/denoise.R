
#' Smooth data by PCA
#'
#' @param x A data matrix with genes as rows and cells as columns
#' @param elbow_th The fraction of PC sdev drop that is considered significant; low values will lead to more PCs being used
#' @param dims_use Directly specify PCs to use, e.g. 1:10
#' @param max_pc Maximum number of PCs computed
#' @param do_plot Plot PC sdev and sdev drop
#' @param scale. Boolean indicating whether genes should be divided by standard deviation after centering and prior to PCA
#'
#' @return De-noised data
#'
#' @importFrom graphics par plot abline
#'
#' @export
#'
smooth_via_pca <- function(x, elbow_th = 0.025, dims_use = NULL, max_pc = 100, do_plot = FALSE,
                           scale. = FALSE) {
  requireNamespace('irlba')
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


#' Denoise data by setting all latent factors to their median values and reversing the regression model
#'
#' @param x A list that provides model parameters and optionally meta data; use output of vst function
#' @param data The name of the entry in x that holds the data
#' @param cell_attr Provide cell meta data holding latent data info
#' @param do_round Round the result to integers
#' @param do_pos Set negative values in the result to zero
#' @param show_progress Whether to print progress bar
#'
#' @return De-noised data as UMI counts
#'
#' @export
#'
denoise <- function(x, data = 'y', cell_attr = x$cell_attr, do_round = TRUE, do_pos = TRUE,
                    show_progress = TRUE) {
  if (is.character(data)) {
    data <- x[[data]]
  }
  # when denoising, set all latent variables to median values
  cell_attr[, x$arguments$latent_var] <- apply(cell_attr[, x$arguments$latent_var, drop=FALSE], 2, function(x) rep(median(x), length(x)))
  regressor_data <- model.matrix(as.formula(gsub('^y', '', x$model_str)), cell_attr)

  genes <- rownames(data)
  bin_size <- x$arguments$bin_size
  bin_ind <- ceiling(x = 1:length(x = genes) / bin_size)
  max_bin <- max(bin_ind)
  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  denoised_data <- matrix(NA, length(genes), nrow(regressor_data), dimnames = list(genes, rownames(regressor_data)))
  for (i in 1:max_bin) {
    genes_bin <- genes[bin_ind == i]
    pearson_residual <- data[genes_bin, ]
    coefs <- x$model_pars_fit[genes_bin, -1]
    theta <- x$model_pars_fit[genes_bin, 1]
    mu <- exp(tcrossprod(coefs, regressor_data))
    variance <- mu + mu^2 / theta
    denoised_data[genes_bin, ] <- mu + pearson_residual * sqrt(variance)
    if (show_progress) {
      setTxtProgressBar(pb, i)
    }
  }
  if (show_progress) {
    close(pb)
  }

  if (do_round) {
    denoised_data <- round(denoised_data, 0)
  }
  if (do_pos) {
    denoised_data[denoised_data < 0] <- 0
  }
  return(denoised_data)
}

reverse_regression <- function(pearson_residual, theta, coefs, data) {
  mu <- exp(data %*% coefs)[, 1]
  variance <- mu + mu^2 / theta
  return(mu + pearson_residual * sqrt(variance))
}
