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
#' @param latent_var_nonreg The non-regularized dependent variables to regress out as a character vector; must match column names in cell_attr; default is NULL
#' @param n_genes Number of genes to use when estimating parameters (default uses 2000 genes, set to NULL to use all genes)
#' @param n_cells Number of cells to use when estimating parameters (default uses all cells)
#' @param method Method to use for initial parameter estimation; one of 'poisson', 'nb_fast', 'nb'
#' @param do_regularize Boolean that, if set to FALSE, will bypass parameter regularization
#' @param res_clip_range Numeric of length two specifying the min and max values the results will be clipped to; default is c(-sqrt(ncol(umi)), sqrt(ncol(umi)))
#' @param bin_size Number of genes to put in each bin (to show progress)
#' @param min_cells Only use genes that have been detected in at least this many cells
#' @param return_cell_attr Make cell attributes part of the output
#' @param return_gene_attr Calculate gene attributes and make part of output
#' @param return_dev_residuals If set to TRUE output will be deviance residuals, NOT Pearson residuals; default is FALSE
#' @param return_corrected_umi If set to TRUE output will contain corrected UMI matrix; see \code{denoise} function
#' @param bw_adjust Kernel bandwidth adjustment factor used during regurlarization; factor will be applied to output of bw.SJ; default is 3
#' @param theta_given Named numeric vector of fixed theta values for the genes; will only be used if method is set to nb_theta_given; default is NULL
#' @param show_progress Whether to print progress bar
#'
#' @return A list with components
#' \item{y}{Matrix of transformed data, i.e. Pearson residuals}
#' \item{umi_corrected}{Matrix of corrected UMI counts (optional)}
#' \item{model_str}{Character representation of the model formula}
#' \item{model_pars}{Matrix of estimated model parameters per gene (theta and regression coefficients)}
#' \item{model_pars_outliers}{Vector indicating whether a gene was considered to be an outlier}
#' \item{model_pars_fit}{Matrix of fitted / regularized model parameters}
#' \item{model_str_nonreg}{Character representation of model for non-regularized variables}
#' \item{model_pars_nonreg}{Model parameters for non-regularized variables}
#' \item{genes_log_mean_step1}{log-mean of genes used in initial step of parameter estimation}
#' \item{cells_step1}{Cells used in initial step of parameter estimation}
#' \item{arguments}{List of function call arguments}
#' \item{cell_attr}{Data frame of cell meta data (optional)}
#' \item{gene_attr}{Data frame with gene attributes such as mean, detection rate, etc. (optional)}
#'
#' @section Details:
#' In the first step of the algorithm, per-gene glm model parameters are learned. This step can be done
#' on a subset of genes and/or cells to speed things up.
#' If \code{method} is set to 'poisson', glm will be called with \code{family = poisson} and
#' the negative binomial theta parameter will be estimated using the response residuals in
#' \code{MASS::theta.ml}.
#' If \code{method} is set to 'nb_fast', glm coefficients and theta are estimated as in the
#' 'poisson' method, but coefficients are then re-estimated using a proper negative binomial
#' model in a second call to glm with
#' \code{family = MASS::negative.binomial(theta = theta)}.
#' If \code{method} is set to 'nb', coefficients and theta are estimated by a single call to
#' \code{MASS::glm.nb}.
#'
#' @import Matrix
#' @import parallel
#' @importFrom MASS theta.ml glm.nb negative.binomial
#' @importFrom stats glm ksmooth model.matrix as.formula approx density poisson var bw.SJ
#' @importFrom utils txtProgressBar setTxtProgressBar capture.output
#'
#' @export
#'
#' @examples
#' \dontrun{
#' vst_out <- vst(pbmc)
#' }
#'
vst <- function(umi,
                cell_attr = NULL,
                latent_var = c('log_umi'),
                batch_var = NULL,
                latent_var_nonreg = NULL,
                n_genes = 2000,
                n_cells = NULL,
                method = 'poisson',
                do_regularize = TRUE,
                res_clip_range = c(-sqrt(ncol(umi)), sqrt(ncol(umi))),
                bin_size = 256,
                min_cells = 5,
                return_cell_attr = FALSE,
                return_gene_attr = FALSE,
                return_dev_residuals = FALSE,
                return_corrected_umi = FALSE,
                bw_adjust = 3,
                theta_given = NULL,
                show_progress = TRUE) {
  arguments <- as.list(environment())[-c(1, 2)]
  start_time <- Sys.time()
  if (is.null(cell_attr)) {
    message('Calculating cell meta data for input UMI matrix')
    cell_attr <- data.frame(umi = colSums(umi),
                            gene = colSums(umi > 0))
    cell_attr$log_umi <- log10(cell_attr$umi)
    cell_attr$log_gene <- log10(cell_attr$gene)
    cell_attr$umi_per_gene <- cell_attr$umi / cell_attr$gene
    cell_attr$log_umi_per_gene <- log10(cell_attr$umi_per_gene)
  }
  if (!all(latent_var %in% colnames(cell_attr))) {
    stop('Not all latent variables present in cell attributes')
  }
  if (!is.null(batch_var)) {
    if (!batch_var %in% colnames(cell_attr)) {
      stop('Batch variable not present in cell attributes')
    }
    batch_levels <- levels(cell_attr[, batch_var])
  }

  # we will generate output for all genes detected in at least min_cells cells
  # but for the first step of parameter estimation we might use only a subset of genes
  genes_cell_count <- rowSums(umi > 0)
  genes <- rownames(umi)[genes_cell_count >= min_cells]
  umi <- umi[genes, ]
  genes_log_mean <- log10(rowMeans(umi))

  if (!is.null(n_cells) && n_cells < ncol(umi)) {
    # downsample cells to speed up the first step
    cells_step1 <- sample(x = colnames(umi), size = n_cells)
    if (!is.null(batch_var)) {
      dropped_batch_levels <- setdiff(batch_levels, levels(droplevels(cell_attr[cells_step1, batch_var])))
      if (length(dropped_batch_levels) > 0) {
        stop('Dropped batch levels ', dropped_batch_levels, ', set n_cells higher')
      }
    }
    genes_cell_count_step1 <- rowSums(umi[, cells_step1] > 0)
    genes_step1 <- rownames(umi)[genes_cell_count_step1 >= min_cells]
    genes_log_mean_step1 <- log10(rowMeans(umi[genes_step1, cells_step1]))
  } else {
    cells_step1 <- colnames(umi)
    genes_step1 <- genes
    genes_log_mean_step1 <- genes_log_mean
  }

  data_step1 <- cell_attr[cells_step1, ]

  if (!is.null(n_genes) && n_genes < length(genes_step1)) {
    # density-sample genes to speed up the first step
    log_mean_dens <- density(x = genes_log_mean_step1, bw = 'nrd', adjust = 1)
    sampling_prob <- 1 / (approx(x = log_mean_dens$x, y = log_mean_dens$y, xout = genes_log_mean_step1)$y + .Machine$double.eps)
    genes_step1 <- sample(x = genes_step1, size = n_genes, prob = sampling_prob)
    genes_log_mean_step1 <- log10(rowMeans(umi[genes_step1, cells_step1]))
  }

  if (!is.null(batch_var)) {
    model_str <- paste0('y ~ (', paste(latent_var, collapse = ' + '), ') : ', batch_var, ' + ', batch_var, ' + 0')
  } else {
    model_str <- paste0('y ~ ', paste(latent_var, collapse = ' + '))
  }

  bin_ind <- ceiling(x = 1:length(x = genes_step1) / bin_size)
  max_bin <- max(bin_ind)
  message('Variance stabilizing transformation of count matrix of size ', nrow(umi), ' by ', ncol(umi))
  message('Model formula is ', model_str)

  model_pars <- get_model_pars(genes_step1, bin_size, umi, model_str, cells_step1, method, data_step1, theta_given, show_progress)

  if (do_regularize) {
    model_pars[, 'theta'] <- log10(model_pars[, 'theta'])
    model_pars_fit <- reg_model_pars(model_pars, genes_log_mean_step1, genes_log_mean, cell_attr,
                                     batch_var, cells_step1, genes_step1, umi, bw_adjust)
    model_pars[, 'theta'] <- 10^model_pars[, 'theta']
    model_pars_fit[, 'theta'] <- 10^model_pars_fit[, 'theta']
    model_pars_outliers <- attr(model_pars_fit, 'outliers')
  } else {
    model_pars_fit <- model_pars
    model_pars_outliers <- rep(FALSE, nrow(model_pars))
  }

  # use all fitted values in NB model
  regressor_data <- model.matrix(as.formula(gsub('^y', '', model_str)), cell_attr)

  if (!is.null(latent_var_nonreg)) {
    message('Estimating parameters for following non-regularized variables: ', latent_var_nonreg)
    if (!is.null(batch_var)) {
      model_str_nonreg <- paste0('y ~ (', paste(latent_var_nonreg, collapse = ' + '), ') : ', batch_var, ' + ', batch_var, ' + 0')
    } else {
      model_str_nonreg <- paste0('y ~ ', paste(latent_var_nonreg, collapse = ' + '))
    }

    model_pars_nonreg <- get_model_pars_nonreg(genes, bin_size, model_pars_fit, regressor_data, umi, model_str_nonreg, cell_attr, show_progress)

    regressor_data_nonreg <- model.matrix(as.formula(gsub('^y', '', model_str_nonreg)), cell_attr)
    model_pars_final <- cbind(model_pars_fit, model_pars_nonreg)
    regressor_data_final <- cbind(regressor_data, regressor_data_nonreg)
    #model_pars_final[, '(Intercept)'] <- model_pars_final[, '(Intercept)'] + model_pars_nonreg[, '(Intercept)']
    #model_pars_final <- cbind(model_pars_final, model_pars_nonreg[, -1, drop=FALSE])
    # model_str <- paste0(model_str, gsub('^y ~ 1', '', model_str2))
  } else {
    model_str_nonreg <- ''
    model_pars_nonreg <- c()
    model_pars_final <- model_pars_fit
    regressor_data_final <- regressor_data
  }

  message('Second step: Pearson residuals using fitted parameters for ', length(x = genes), ' genes')
  bin_ind <- ceiling(x = 1:length(x = genes) / bin_size)
  max_bin <- max(bin_ind)
  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  res <- matrix(NA, length(genes), nrow(regressor_data_final), dimnames = list(genes, rownames(regressor_data_final)))
  for (i in 1:max_bin) {
    genes_bin <- genes[bin_ind == i]
    mu <- exp(tcrossprod(model_pars_final[genes_bin, -1, drop=FALSE], regressor_data_final))
    y <- as.matrix(umi[genes_bin, , drop=FALSE])
    if (return_dev_residuals) {
      res[genes_bin, ] <- deviance_residual(y, mu, model_pars_final[genes_bin, 'theta'])
    } else {
      res[genes_bin, ] <- (y - mu) / sqrt(mu + mu^2 / model_pars_final[genes_bin, 'theta'])
    }
    if (show_progress) {
      setTxtProgressBar(pb, i)
    }
  }
  if (show_progress) {
    close(pb)
  }
  rv <- list(y = res,
             model_str = model_str,
             model_pars = model_pars,
             model_pars_outliers = model_pars_outliers,
             model_pars_fit = model_pars_fit,
             model_str_nonreg = model_str_nonreg,
             model_pars_nonreg = model_pars_nonreg,
             arguments = arguments,
             genes_log_mean_step1 = genes_log_mean_step1,
             cells_step1 = cells_step1,
             cell_attr = cell_attr)
  rm(res)
  gc(verbose = FALSE)

  if (return_corrected_umi) {
    rv$umi_corrected <- sctransform::denoise(rv, do_round = TRUE, do_pos = TRUE,
                                             show_progress = show_progress)
    rv$umi_corrected <- as(object = rv$umi_corrected, Class = 'dgCMatrix')
  }

  rv$y[rv$y < res_clip_range[1]] <- res_clip_range[1]
  rv$y[rv$y > res_clip_range[2]] <- res_clip_range[2]

  if (!return_cell_attr) {
    rv[['cell_attr']] <- NULL
  }

  if (return_gene_attr) {
    message('Calculating gene attributes')
    gene_attr <- data.frame(
      detection_rate = genes_cell_count[genes] / ncol(umi),
      mean = 10 ^ genes_log_mean,
      variance = apply(umi, 1, var),
      residual_mean = rowMeans(rv$y)
    )
    if (requireNamespace('matrixStats', quietly = TRUE)) {
      gene_attr$residual_variance = matrixStats::rowVars(rv$y)
    } else {
      message('Consider installing matrixStats package for faster gene attribute calculation.')
      gene_attr$residual_variance = apply(rv$y, 1, var)
    }
    rv[['gene_attr']] <- gene_attr
  }

  message('Wall clock passed: ', capture.output(print(Sys.time() - start_time)))
  return(rv)
}


get_model_pars <- function(genes_step1, bin_size, umi, model_str, cells_step1, method, data_step1, theta_given, show_progress) {
  bin_ind <- ceiling(x = 1:length(x = genes_step1) / bin_size)
  max_bin <- max(bin_ind)
  message('Get Negative Binomial regression parameters per gene')
  message('Using ', length(x = genes_step1), ' genes, ', length(x = cells_step1), ' cells')

  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  model_pars <- list()
  for (i in 1:max_bin) {
    genes_bin_regress <- genes_step1[bin_ind == i]
    umi_bin <- as.matrix(umi[genes_bin_regress, cells_step1, drop=FALSE])
    model_pars[[i]] <- do.call(rbind,
                               mclapply(
                                 X = genes_bin_regress,
                                 FUN = function(j) {
                                   y <- umi_bin[j, ]
                                   if (method == 'poisson') {
                                     fit <- glm(as.formula(model_str), data = data_step1, family = poisson)
                                     theta <- as.numeric(x = theta.ml(y = y, mu = fit$fitted))
                                     return(c(theta, fit$coefficients))
                                   }
                                   if (method == 'nb_theta_given') {
                                     theta <- theta_given[j]
                                     fit2 <- 0
                                     try(fit2 <- glm(as.formula(model_str), data = data_step1, family = negative.binomial(theta=theta)), silent=TRUE)
                                     if (class(fit2)[1] == 'numeric') {
                                       return(c(theta, glm(as.formula(model_str), data = data_step1, family = poisson)$coefficients))
                                     } else {
                                       return(c(theta, fit2$coefficients))
                                     }
                                   }
                                   if (method == 'nb_fast') {
                                     fit <- glm(as.formula(model_str), data = data_step1, family = poisson)
                                     theta <- as.numeric(x = theta.ml(y = y, mu = fit$fitted))
                                     fit2 <- 0
                                     try(fit2 <- glm(as.formula(model_str), data = data_step1, family = negative.binomial(theta=theta)), silent=TRUE)
                                     if (class(fit2)[1] == 'numeric') {
                                       return(c(theta, fit$coefficients))
                                     } else {
                                       return(c(theta, fit2$coefficients))
                                     }
                                   }
                                   if (method == 'nb') {
                                     fit <- 0
                                     try(fit <- glm.nb(as.formula(model_str), data = data_step1), silent=TRUE)
                                     if (class(fit)[1] == 'numeric') {
                                       fit <- glm(as.formula(model_str), data = data_step1, family = poisson)
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
  return(model_pars)
}

get_model_pars_nonreg <- function(genes, bin_size, model_pars_fit, regressor_data, umi, model_str_nonreg, cell_attr, show_progress) {
  bin_ind <- ceiling(x = 1:length(x = genes) / bin_size)
  max_bin <- max(bin_ind)
  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  model_pars_nonreg <- list()
  for (i in 1:max_bin) {
    genes_bin <- genes[bin_ind == i]
    mu <- tcrossprod(model_pars_fit[genes_bin, -1, drop=FALSE], regressor_data)
    umi_bin <- as.matrix(umi[genes_bin, ])
    model_pars_nonreg[[i]] <- do.call(rbind,
                                      mclapply(genes_bin, function(gene) {
                                        fam <- negative.binomial(theta = model_pars_fit[gene, 'theta'], link = 'log')
                                        y <- umi_bin[gene, ]
                                        offs <- mu[gene, ]
                                        fit <- glm(as.formula(model_str_nonreg), data = cell_attr, family = fam, offset=offs)
                                        return(fit$coefficients)
                                      }))
    if (show_progress) {
      setTxtProgressBar(pb, i)
    }
  }
  if (show_progress) {
    close(pb)
  }
  model_pars_nonreg <- do.call(rbind, model_pars_nonreg)
  rownames(model_pars_nonreg) <- genes
  return(model_pars_nonreg)
}

reg_model_pars <- function(model_pars, genes_log_mean_step1, genes_log_mean, cell_attr,
                           batch_var, cells_step1, genes_step1, umi, bw_adjust) {
  genes <- names(genes_log_mean)
  # look for outliers in the parameters
  # outliers are those that do not fit the overall relationship with the mean at all
  outliers <- apply(model_pars, 2, function(y) is_outlier(y, genes_log_mean_step1))
  outliers <- apply(outliers, 1, any)
  if (sum(outliers) > 0) {
    message('Found ', sum(outliers), ' outliers - those will be ignored in fitting/regularization step\n')
    model_pars <- model_pars[!outliers, ]
    genes_step1 <- rownames(model_pars)
    genes_log_mean_step1 <- genes_log_mean_step1[!outliers]
  }

  # select bandwidth to be used for smoothing
  bw <- bw.SJ(genes_log_mean_step1) * bw_adjust

  # for parameter predictions
  x_points <- pmax(genes_log_mean, min(genes_log_mean_step1))
  x_points <- pmin(x_points, max(genes_log_mean_step1))

  # take results from step 1 and fit/predict parameters to all genes
  o <- order(x_points)
  model_pars_fit <- matrix(NA, length(genes), ncol(model_pars),
                           dimnames = list(genes, colnames(model_pars)))

  # fit / regularize theta
  model_pars_fit[o, 'theta'] <- ksmooth(x = genes_log_mean_step1, y = model_pars[, 'theta'],
                                             x.points = x_points, bandwidth = bw, kernel='normal')$y

  if (is.null(batch_var)){
    # global fit / regularization for all coefficients
    for (i in 2:ncol(model_pars)) {
      model_pars_fit[o, i] <- ksmooth(x = genes_log_mean_step1, y = model_pars[, i],
                                      x.points = x_points, bandwidth = bw, kernel='normal')$y
    }
  } else {
    # fit / regularize per batch
    batches <- unique(cell_attr[, batch_var])
    for (b in batches) {
      sel <- cell_attr[, batch_var] == b & rownames(cell_attr) %in% cells_step1
      batch_genes_log_mean_step1 <- log10(rowMeans(umi[genes_step1, sel]))
      sel <- cell_attr[, batch_var] == b
      batch_genes_log_mean <- log10(rowMeans(umi[, sel]))
      # in case some genes have not been observed in this batch
      batch_genes_log_mean <- pmax(batch_genes_log_mean, min(genes_log_mean))
      batch_o <- order(batch_genes_log_mean)
      for (i in which(grepl(paste0(batch_var, b), colnames(model_pars)))) {
        model_pars_fit[batch_o, i] <- ksmooth(x = batch_genes_log_mean_step1, y = model_pars[, i],
                                              x.points = batch_genes_log_mean, bandwidth = bw, kernel='normal')$y
      }
    }
  }
  attr(model_pars_fit, 'outliers') <- outliers
  return(model_pars_fit)
}
