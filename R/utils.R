#' Geometric mean
#'
#' @param x array
#' @param eps small value to add to x to avoid log(0); default is 1
#'
#' @return geometric mean
gmean <- function(x, eps=1) {
  exp(sum(log(x+eps))/length(x))-eps
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
  if (class(tmp$x) == 'list') {
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

pearson_residual <- function(y, mu, theta) {
  (y - mu) / sqrt(mu + mu^2 / theta)
}

sq_deviance_residual <- function(y, mu, theta, wt=1) {
  2 * wt * (y * log(pmax(1, y)/mu) - (y + theta) * log((y + theta)/(mu + theta)))
}

deviance_residual <- function(y, mu, theta, wt=1) {
  r <- 2 * wt * (y * log(pmax(1, y)/mu) - (y + theta) * log((y + theta)/(mu + theta)))
  sqrt(r) * sign(y - mu)
}

#' Return deviance residuals of regularized models
#'
#' @param vst_out The output of a vst run
#' @param umi The UMI count matrix that will be converted
#' @param cell_attr Data frame of cell meta data
#' @param bin_size Number of genes to put in each bin (to show progress)
#' @param show_progress Whether to print progress bar
#'
#' @return A matrix of deviance residuals
#'
#' @export
#'
#' @examples
#' \dontrun{
#' vst_out <- vst(pbmc, return_gene_attr = TRUE)
#' dev_res <- get_deviance_residuals(vst_out, pbmc)
#' }
#'
get_deviance_residuals <- function(vst_out, umi, cell_attr = vst_out$cell_attr,
                                   bin_size = 256, show_progress = TRUE) {
  regressor_data <- model.matrix(as.formula(gsub('^y', '', vst_out$model_str)), cell_attr)
  model_pars <- vst_out$model_pars_fit
  if (!is.null(dim(vst_out$model_pars_nonreg))) {
    regressor_data_nonreg <- model.matrix(as.formula(gsub('^y', '', vst_out$model_str_nonreg)), cell_attr)
    regressor_data <- cbind(regressor_data, regressor_data_nonreg)
    model_pars <- cbind(vst_out$model_pars_fit, vst_out$model_pars_nonreg)
  }

  genes <- rownames(model_pars)[rownames(model_pars) %in% rownames(umi)]
  bin_ind <- ceiling(x = 1:length(x = genes) / bin_size)
  max_bin <- max(bin_ind)
  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  res <- matrix(NA, length(genes), nrow(regressor_data), dimnames = list(genes, rownames(regressor_data)))
  for (i in 1:max_bin) {
    genes_bin <- genes[bin_ind == i]
    mu <- exp(tcrossprod(model_pars[genes_bin, -1, drop=FALSE], regressor_data))
    y <- as.matrix(umi[genes_bin, , drop=FALSE])
    res[genes_bin, ] <- deviance_residual(y, mu, model_pars[genes_bin, 'theta'])
    if (show_progress) {
      setTxtProgressBar(pb, i)
    }
  }
  if (show_progress) {
    close(pb)
  }
  return(res)
}
