#' Generate data from regularized models.
#'
#' Generate data from regularized models. This generates data from the background,
#' i.e. no residuals are added to the simulated data. The cell attributes for the
#' generated cells are sampled from the input with replacement.
#'
#' @param vst_out A list that provides model parameters and optionally meta data; use output of vst function
#' @param genes The gene names for which to generate data; default is rownames(vst_out$model_pars_fit)
#' @param cell_attr Provide cell meta data holding latent data info; default is vst_out$cell_attr
#' @param n_cells Number of cells to generate; default is nrow(cell_attr)
#'
#' @return Generated data as dgCMatrix
#'
#' @importFrom methods as
#'
#' @export
#'
#' @examples
#' \donttest{
#' vst_out <- vst(pbmc, return_cell_attr = TRUE)
#' generated_data <- generate(vst_out)
#' }
#'
generate <- function(vst_out, genes = rownames(vst_out$model_pars_fit),
                     cell_attr = vst_out$cell_attr,
                     n_cells = nrow(cell_attr)) {
  genes <- genes[genes %in% rownames(vst_out$model_pars_fit)]
  # get model parameters
  mp <- vst_out$model_pars_fit[genes, , drop = FALSE]
  coefs <- mp[, -1, drop=FALSE]
  theta <- mp[, 1]
  # we sample from the original list of cell attributes when we generate data
  # choose cells here
  idx <- sample(x = nrow(cell_attr), size = n_cells, replace = TRUE)
  regressor_data <- cbind(rep(1, length(idx)), cell_attr[idx, colnames(coefs)[-1]])
  # calculate expected values
  mu <- exp(tcrossprod(coefs, regressor_data))
  x.sim <- t(sapply(rownames(mu), function(gene) {
    gene.mu <- mu[gene, ]
    x <- MASS::rnegbin(n = length(gene.mu), mu = gene.mu, theta = theta[gene])
    return(x)
  }))
  x.sim <- as(x.sim, Class = 'dgCMatrix')
  return(x.sim)
}
