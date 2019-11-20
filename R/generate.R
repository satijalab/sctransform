#' Generate data from regularized models.
#'
#' Generate data from regularized models. If background_only is set to TRUE no residuals are
#' added to the simulated data.
#'
#' @param vst_out A list that provides model parameters and optionally meta data; use output of vst function
#' @param data The name of the entry in x that holds the data (Pearson residuals)
#' @param cell_attr Provide cell meta data holding latent data info
#' @param background_only Set to TRUE to only generate the background signal (do not add back residuals)
#' @param do_round Round the result to integers
#' @param do_pos Set negative values in the result to zero
#' @param make_sparse Make the output a dgCMatrix (if do_pos and do_round are also TRUE)
#'
#' @return Generated data
#'
#' @export
#'
#' @examples
#' vst_out <- vst(pbmc, return_cell_attr = TRUE, res_clip_range = c(-Inf, Inf))
#' generated_data <- generate(vst_out)
#'
generate <- function(vst_out, data = 'y', genes = rownames(vst_out$model_pars_fit),
                     cell_attr = x$cell_attr, background_only = FALSE,
                     do_round = TRUE, do_pos = TRUE,
                     make_sparse = FALSE) {
  genes <- genes[genes %in% rownames(vst_out$model_pars_fit)]
  # when we generate data, use the cell attributes from the original data
  if (is.character(data)) {
    data <- vst_out[[data]]
  }
  data <- data[genes, , drop = FALSE]
  # get model parameters
  mp <- vst_out$model_pars_fit[genes, , drop = FALSE]
  coefs <- mp[, -1, drop=FALSE]
  theta <- mp[, 1]
  ca <- vst_out$cell_attr
  regressor_data <- cbind(rep(1, nrow(ca)), ca[, colnames(coefs)[-1]])
  # calculate expected values
  mu <- exp(tcrossprod(coefs, regressor_data))
  x.sim <- t(sapply(rownames(mu), function(gene) {
    gene.mu <- mu[gene, ]
    x <- MASS::rnegbin(n = length(gene.mu), mu = gene.mu, theta = theta[gene])
    if (!background_only) {
      this.sd <- sqrt(gene.mu + gene.mu^2 / theta[gene])
      # add back residuals
      x <- x + this.sd * data[gene, ]
    }
    return(x)
  }))
  if (do_round) {
    x.sim <- round(x.sim, 0)
  }
  if (do_pos) {
    x.sim[x.sim < 0] <- 0
  }
  if (do_round & do_pos & make_sparse) {
    x.sim <- as(x.sim, Class = 'dgCMatrix')
  }
  return(x.sim)
}
