#' Variance stabilizing transformation for UMI count data
#'
#' Apply variance stabilizing transformation to UMI count data using a regularized Negative Binomial regression model.
#' This will remove unwanted effects from UMI data and return Pearson residuals.
#' Uses mclapply; you can set the number of cores it will use to n with command options(mc.cores = n).
#' If n_genes is set, only a (somewhat-random) subset of genes is used for estimating the
#' initial model parameters.
#'
#' @param umi A matrix of UMI counts with genes as rows and cells as columns
#' @param cell_attr A data frame containing the dependent variables; if omitted a data frame with umi and gene will be generated
#' @param latent_var The dependent variables to regress out as a character vector; must match column names in cell_attr; default is c("log_umi_per_gene")
#' @param batch_var The dependent variables indicating which batch a cell belongs to; no batch interaction terms used if omiited
#' @param n_genes Number of genes to use when estimating parameters (default uses 2000 genes, set to NULL to use all genes)
#' @param method Method to use for initial parameter estimation; one of 'poisson', 'nb_fast', 'nb'
#' @param res_clip_range Numeric of length two specifying the min and max values the results will be clipped to
#' @param bin_size Number of genes to put in each bin (to show progress)
#' @param min_cells Only use genes that have been detected in at least this many cells
#' @param return_cell_attr Make cell attributes part of the output
#' @param return_gene_attr Calculate gene attributes and make part of output
#' @param show_progress Whether to print progress bar
#'
#' @return A list with components
#' \item{y}{Matrix of transformed data, i.e. Pearson residuals}
#' \item{model_str}{Character representation of the model formula}
#' \item{model_pars}{Matrix of estimated model parameters per gene (theta and regression coefficients)}
#' \item{model_pars_fit}{Matrix of fitted / regularized model parameters}
#' \item{arguments}{List of function call arguments}
#' \item{cell_attr}{Data frame of cell meta data (optional)}
#' \item{gene_attr}{Data frame with gene attributes such as mean, detection rate, etc. (optional)}
#'
#' @import Matrix
#' @import parallel
#' @importFrom MASS theta.ml glm.nb negative.binomial
#' @importFrom stats glm ksmooth model.matrix as.formula approx density poisson var
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
#' @examples
#'
#'
vst <- function(umi,
                cell_attr = NULL,
                latent_var = c('log_umi_per_gene'),
                batch_var = NULL,
                n_genes = 2000,
                method = 'poisson',
                res_clip_range = c(-50, 50),
                bin_size = 256,
                min_cells = 5,
                return_cell_attr = FALSE,
                return_gene_attr = FALSE,
                show_progress = TRUE) {
  arguments <- as.list(environment())[-c(1, 2)]
  if (is.null(cell_attr)) {
    message('Calculating cell meta data for input UMI matrix')
    cell_attr <- data.frame(umi = colSums(umi),
                            gene = colSums(umi > 0))
    cell_attr$log_umi <- log10(cell_attr$umi)
    cell_attr$log_gene <- log10(cell_attr$gene)
    cell_attr$log_umi_per_gene <- log10(cell_attr$umi / cell_attr$gene)
  }
  if (!all(latent_var %in% colnames(cell_attr))) {
    stop('Not all latent variables present in meta data')
  }
  if (!is.null(batch_var)) {
    if (!batch_var %in% colnames(cell_attr)) {
      stop('Batch variable not present in meta data')
    }
  }

  gene_cell_count <- rowSums(umi > 0)
  genes <- rownames(umi)[gene_cell_count >= min_cells]
  umi <- umi[genes, ]
  genes_log_mean <- log10(rowMeans(umi))

  if (!is.null(n_genes)) {
    # density-sample genes to speed up the first step
    log_mean_dens <- density(x = genes_log_mean, bw = 'nrd', adjust = 1)
    sampling_prob <- 1 / (approx(x = log_mean_dens$x, y = log_mean_dens$y, xout = genes_log_mean)$y + .Machine$double.eps)
    genes_step1 <- sample(x = genes, size = n_genes, prob = sampling_prob)
  } else {
    genes_step1 <- genes
  }

  if (!is.null(batch_var)) {
    model_str <- paste0('y ~ (', paste(latent_var, collapse = ' + '), ') : ', batch_var, ' + ', batch_var, ' + 0')
  } else {
    model_str <- paste0('y ~ ', paste(latent_var, collapse = ' + '))
  }

  bin_ind <- ceiling(x = 1:length(x = genes_step1) / bin_size)
  max_bin <- max(bin_ind)
  message(paste("Regressing out", paste(latent_var, collapse = ' ')))
  message('Model formula is ', model_str)
  message('First step: Poisson regression (to get initial model), and estimate theta per gene')
  message('Using ', length(x = genes_step1), ' genes')

  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  model_pars <- list()
  for (i in 1:max_bin) {
    genes_bin_regress <- genes_step1[bin_ind == i]
    umi_bin <- as.matrix(umi[genes_bin_regress, ])
    model_pars[[i]] <- do.call(rbind,
      mclapply(
        X = genes_bin_regress,
        FUN = function(j) {
          y <- umi_bin[j, ]
          if (method == 'poisson') {
            fit <- glm(as.formula(model_str), data = cell_attr, family = poisson)
            theta <- as.numeric(x = theta.ml(y = y, mu = fit$fitted))
            return(c(theta, fit$coefficients))
          }
          if (method == 'nb_fast') {
            fit <- glm(as.formula(model_str), data = cell_attr, family = poisson)
            theta <- as.numeric(x = theta.ml(y = y, mu = fit$fitted))
            fit2 <- 0
            try(fit2 <- glm(as.formula(model_str), data = cell_attr, family = negative.binomial(theta=theta)), silent=TRUE)
            if (class(fit2)[1] == 'numeric') {
              return(c(theta, fit$coefficients))
            } else {
              return(c(theta, fit2$coefficients))
            }
          }
          if (method == 'nb') {
            fit <- 0
            try(fit <- glm.nb(as.formula(model_str), data = cell_attr), silent=TRUE)
            if (class(fit)[1] == 'numeric') {
              fit <- glm(as.formula(model_str), data = cell_attr, family = poisson)
              fit$theta <- as.numeric(x = theta.ml(y = y, mu = fit$fitted))
            }
            return(c(fit$theta, fit$coefficients))
          }
        }
      )
    )
    if (show_progress) {
      setTxtProgressBar(pb, i)
    }
  }
  model_pars <- do.call(rbind, model_pars)
  if (show_progress) {
    close(pb)
  }
  rownames(model_pars) <- genes_step1
  colnames(model_pars)[1] <- 'theta'

  # look for outliers in the parameters
  # outliers are thos that do not fit the relationship with the mean
  outliers <- apply(model_pars, 2, function(y) is_outlier(y, genes_log_mean[genes_step1], 100))
  outliers <- apply(outliers, 1, any)
  if (sum(outliers) > 0) {
    message('Found ', sum(outliers), ' outliers - those will be ignored in fitting/regularization step\n')
    model_pars <- model_pars[!outliers, ]
    genes_step1 <- rownames(model_pars)
  }

  # take results from step 1 and fit/predict parameters to all genes
  o <- order(genes_log_mean)
  model_pars_fit <- matrix(NA, length(genes), ncol(model_pars),
                              dimnames = list(genes, colnames(model_pars)))
  # fit / regularize theta
  model_pars_fit[o, 'theta'] <- 10 ^ ksmooth(x = genes_log_mean[genes_step1], y = log10(model_pars[, 'theta']),
                                                x.points = genes_log_mean, bandwidth = 0.3)$y

  if (is.null(batch_var)){
    # global fit / regularization for all coefficients
    for (i in 2:ncol(model_pars)) {
      model_pars_fit[o, i] <- ksmooth(x = genes_log_mean[genes_step1], y = model_pars[, i],
                                         x.points = genes_log_mean, bandwidth = 0.3)$y
    }
  } else {
    # fit / regularize per batch
    batches <- unique(cell_attr[, batch_var])
    for (b in batches) {
      sel <- cell_attr[, batch_var] == b
      batch_genes_log_mean <- log10(rowMeans(umi[, sel]))
      # in case some genes have not been observed in this batch
      batch_genes_log_mean <- pmax(batch_genes_log_mean, min(genes_log_mean))
      batch_o <- order(batch_genes_log_mean)
      for (i in which(grepl(paste0(batch_var, b), colnames(model_pars)))) {
        model_pars_fit[batch_o, i] <- ksmooth(x = batch_genes_log_mean[genes_step1], y = model_pars[, i],
                                              x.points = batch_genes_log_mean, bandwidth = 0.3)$y
      }
    }
  }

  message('Second step: Pearson residuals using fitted parameters for ', length(x = genes), ' genes')
  # use all fitted values in NB model
  regressor_data <- model.matrix(as.formula(gsub('^y', '', model_str)), cell_attr)
  bin_ind <- ceiling(x = 1:length(x = genes) / bin_size)
  max_bin <- max(bin_ind)
  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  res <- matrix(NA, length(genes), nrow(regressor_data), dimnames = list(genes, rownames(regressor_data)))
  for (i in 1:max_bin) {
    genes_bin <- genes[bin_ind == i]
    mu <- exp(model_pars_fit[genes_bin, -1, drop=FALSE] %*% t(regressor_data))
    y <- as.matrix(umi[genes_bin, ])
    res[genes_bin, ] <- (y - mu) / sqrt(mu + mu^2 / model_pars_fit[genes_bin, 'theta'])
    if (show_progress) {
      setTxtProgressBar(pb, i)
    }
  }
  if (show_progress) {
    close(pb)
  }
  res[res < res_clip_range[1]] <- res_clip_range[1]
  res[res > res_clip_range[2]] <- res_clip_range[2]

  rv <- list(y = res,
             model_str = model_str,
             model_pars = model_pars,
             model_pars_fit = model_pars_fit,
             arguments = arguments)

  if (return_cell_attr) {
    rv[['cell_attr']] <- cell_attr
  }
  if (return_gene_attr) {
    message('Calculating gene attributes')
    gene_attr <- data.frame(
      detection_rate = gene_cell_count[genes] / ncol(umi),
      mean = 10 ^ genes_log_mean,
      variance = apply(umi, 1, var),
      residual_mean = apply(res, 1, mean),
      residual_variance = apply(res, 1, var)
    )
    rv[['gene_attr']] <- gene_attr
  }
  return(rv)
}
