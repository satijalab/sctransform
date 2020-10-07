#' @useDynLib sctransform
NULL


#' Variance stabilizing transformation for UMI count data
#'
#' Apply variance stabilizing transformation to UMI count data using a regularized Negative Binomial regression model.
#' This will remove unwanted effects from UMI data and return Pearson residuals.
#' Uses future_lapply; you can set the number of cores it will use to n with plan(strategy = "multicore", workers = n).
#' If n_genes is set, only a (somewhat-random) subset of genes is used for estimating the
#' initial model parameters.
#'
#' @param umi A matrix of UMI counts with genes as rows and cells as columns
#' @param cell_attr A data frame containing the dependent variables; if omitted a data frame with umi and gene will be generated
#' @param latent_var The independent variables to regress out as a character vector; must match column names in cell_attr; default is c("log_umi")
#' @param batch_var The dependent variables indicating which batch a cell belongs to; no batch interaction terms used if omiited
#' @param latent_var_nonreg The non-regularized dependent variables to regress out as a character vector; must match column names in cell_attr; default is NULL
#' @param n_genes Number of genes to use when estimating parameters (default uses 2000 genes, set to NULL to use all genes)
#' @param n_cells Number of cells to use when estimating parameters (default uses all cells)
#' @param method Method to use for initial parameter estimation; one of 'poisson', 'qpoisson', 'nb_fast', 'nb', 'nb_theta_given', 'glmGamPoi'
#' @param do_regularize Boolean that, if set to FALSE, will bypass parameter regularization and use all genes in first step (ignoring n_genes); default is FALSE
#' @param theta_regularization Method to use to regularize theta; use 'log_theta' for the behavior prior to version 0.3; default is 'od_factor'
#' @param res_clip_range Numeric of length two specifying the min and max values the results will be clipped to; default is c(-sqrt(ncol(umi)), sqrt(ncol(umi)))
#' @param bin_size Number of genes to process simultaneously; this will determine how often the progress bars are updated and how much memory is being used; default is 500
#' @param min_cells Only use genes that have been detected in at least this many cells; default is 5
#' @param residual_type What type of residuals to return; can be 'pearson', 'deviance', or 'none'; default is 'pearson'
#' @param return_cell_attr Make cell attributes part of the output; default is FALSE
#' @param return_gene_attr Calculate gene attributes and make part of output; default is TRUE
#' @param return_corrected_umi If set to TRUE output will contain corrected UMI matrix; see \code{correct} function
#' @param min_variance Lower bound for the estimated variance for any gene in any cell when calculating pearson residual; default is -Inf
#' @param bw_adjust Kernel bandwidth adjustment factor used during regurlarization; factor will be applied to output of bw.SJ; default is 3
#' @param gmean_eps Small value added when calculating geometric mean of a gene to avoid log(0); default is 1
#' @param theta_estimation_fun Character string indicating which method to use to estimate theta (when method = poisson); default is 'theta.ml', but 'theta.mm' seems to be a good and fast alternative
#' @param theta_given Named numeric vector of fixed theta values for the genes; will only be used if method is set to nb_theta_given; default is NULL
#' @param verbosity An integer specifying whether to show only messages (1), messages and progress bars (2) or nothing (0) while the function is running; default is 2
#' @param verbose Deprecated; use verbosity instead
#' @param show_progress Deprecated; use verbosity instead
#'
#' @return A list with components
#' \item{y}{Matrix of transformed data, i.e. Pearson residuals, or deviance residuals; empty if \code{residual_type = 'none'}}
#' \item{umi_corrected}{Matrix of corrected UMI counts (optional)}
#' \item{model_str}{Character representation of the model formula}
#' \item{model_pars}{Matrix of estimated model parameters per gene (theta and regression coefficients)}
#' \item{model_pars_outliers}{Vector indicating whether a gene was considered to be an outlier}
#' \item{model_pars_fit}{Matrix of fitted / regularized model parameters}
#' \item{model_str_nonreg}{Character representation of model for non-regularized variables}
#' \item{model_pars_nonreg}{Model parameters for non-regularized variables}
#' \item{genes_log_gmean_step1}{log-geometric mean of genes used in initial step of parameter estimation}
#' \item{cells_step1}{Cells used in initial step of parameter estimation}
#' \item{arguments}{List of function call arguments}
#' \item{cell_attr}{Data frame of cell meta data (optional)}
#' \item{gene_attr}{Data frame with gene attributes such as mean, detection rate, etc. (optional)}
#' \item{times}{Time stamps at various points in the function}
#'
#' @section Details:
#' In the first step of the algorithm, per-gene glm model parameters are learned. This step can be done
#' on a subset of genes and/or cells to speed things up.
#' If \code{method} is set to 'poisson', a poisson regression is done and
#' the negative binomial theta parameter is estimated using the response residuals in
#' \code{theta_estimation_fun}.
#' If \code{method} is set to 'qpoisson', coefficients and overdispersion (phi) are estimated by quasi 
#' poisson regression and theta is estimated based on phi and the mean fitted value - this is currently 
#' the fastest method with results very similar to 'glmGamPoi'
#' If \code{method} is set to 'nb_fast', coefficients and theta are estimated as in the
#' 'poisson' method, but coefficients are then re-estimated using a proper negative binomial
#' model in a second call to glm with \code{family = MASS::negative.binomial(theta = theta)}.
#' If \code{method} is set to 'nb', coefficients and theta are estimated by a single call to
#' \code{MASS::glm.nb}.
#' If \code{method} is set to 'glmGamPoi', coefficients and theta are estimated by a single call to
#' \code{glmGamPoi::glm_gp}.
#'
#' @import Matrix
#' @importFrom future.apply future_lapply
#' @importFrom MASS theta.ml theta.mm glm.nb negative.binomial
#' @importFrom stats glm glm.fit df.residual ksmooth model.matrix as.formula approx density poisson var bw.SJ
#' @importFrom utils txtProgressBar setTxtProgressBar capture.output
#' @importFrom methods as
#'
#' @export
#'
#' @examples
#' \donttest{
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
                theta_regularization = 'od_factor',
                res_clip_range = c(-sqrt(ncol(umi)), sqrt(ncol(umi))),
                bin_size = 500,
                min_cells = 5,
                residual_type = 'pearson',
                return_cell_attr = FALSE,
                return_gene_attr = TRUE,
                return_corrected_umi = FALSE,
                min_variance = -Inf,
                bw_adjust = 3,
                gmean_eps = 1,
                theta_estimation_fun = 'theta.ml',
                theta_given = NULL,
                verbosity = 2,
                verbose = NULL,
                show_progress = NULL) {
  arguments <- as.list(environment())
  arguments <- arguments[!names(arguments) %in% c("umi", "cell_attr")]

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


  # Check for suggested package
  if (method == "glmGamPoi") {
    glmGamPoi_check <- requireNamespace("glmGamPoi", quietly = TRUE)
    if (!glmGamPoi_check){
      stop('Please install the glmGamPoi package. See https://github.com/const-ae/glmGamPoi for details.')
    }
  }
  

  times <- list(start_time = Sys.time())

  cell_attr <- make_cell_attr(umi, cell_attr, latent_var, batch_var, latent_var_nonreg, verbosity)
  if (!is.null(batch_var)) {
    cell_attr[, batch_var] <- as.factor(cell_attr[, batch_var])
    batch_levels <- levels(cell_attr[, batch_var])
  }

  # we will generate output for all genes detected in at least min_cells cells
  # but for the first step of parameter estimation we might use only a subset of genes
  genes_cell_count <- rowSums(umi >= 0.01)
  genes <- rownames(umi)[genes_cell_count >= min_cells]
  umi <- umi[genes, ]
  genes_log_gmean <- log10(row_gmean(umi, eps = gmean_eps))

  if (!do_regularize) {
    if (verbosity > 0) {
      message('do_regularize is set to FALSE, will use all genes')
    }
    n_genes <- NULL
  }

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
    genes_log_gmean_step1 <- log10(row_gmean(umi[genes_step1, cells_step1], eps = gmean_eps))
  } else {
    cells_step1 <- colnames(umi)
    genes_step1 <- genes
    genes_log_gmean_step1 <- genes_log_gmean
  }

  data_step1 <- cell_attr[cells_step1, , drop = FALSE]

  if (!is.null(n_genes) && n_genes < length(genes_step1)) {
    # density-sample genes to speed up the first step
    log_gmean_dens <- density(x = genes_log_gmean_step1, bw = 'nrd', adjust = 1)
    sampling_prob <- 1 / (approx(x = log_gmean_dens$x, y = log_gmean_dens$y, xout = genes_log_gmean_step1)$y + .Machine$double.eps)
    genes_step1 <- sample(x = genes_step1, size = n_genes, prob = sampling_prob)
    genes_log_gmean_step1 <- log10(row_gmean(umi[genes_step1, cells_step1], eps = gmean_eps))
  }

  if (!is.null(batch_var)) {
    model_str <- paste0('y ~ (', paste(latent_var, collapse = ' + '), ') : ', batch_var, ' + ', batch_var, ' + 0')
  } else {
    model_str <- paste0('y ~ ', paste(latent_var, collapse = ' + '))
  }

  bin_ind <- ceiling(x = 1:length(x = genes_step1) / bin_size)
  max_bin <- max(bin_ind)
  if (verbosity > 0) {
    message('Variance stabilizing transformation of count matrix of size ', nrow(umi), ' by ', ncol(umi))
    message('Model formula is ', model_str)
  }

  times$get_model_pars = Sys.time()
  model_pars <- get_model_pars(genes_step1, bin_size, umi, model_str, cells_step1,
                               method, data_step1, theta_given, theta_estimation_fun,
                               verbosity)
  # make sure theta is not too small
  min_theta <- 1e-7
  if (any(model_pars[, 'theta'] < min_theta)) {
    if (verbosity > 0) {
      msg <- sprintf('There are %d estimated thetas smaller than %g - will be set to %g', sum(model_pars[, 'theta'] < min_theta), min_theta, min_theta)
      message(msg)
    }
    model_pars[, 'theta'] <- pmax(model_pars[, 'theta'], min_theta)
  }


  times$reg_model_pars = Sys.time()
  if (do_regularize) {
    model_pars_fit <- reg_model_pars(model_pars, genes_log_gmean_step1, genes_log_gmean, cell_attr,
                                     batch_var, cells_step1, genes_step1, umi, bw_adjust, gmean_eps,
                                     theta_regularization, verbosity)
    model_pars_outliers <- attr(model_pars_fit, 'outliers')
  } else {
    model_pars_fit <- model_pars
    model_pars_outliers <- rep(FALSE, nrow(model_pars))
  }

  # use all fitted values in NB model
  regressor_data <- model.matrix(as.formula(gsub('^y', '', model_str)), cell_attr)

  if (!is.null(latent_var_nonreg)) {
    if (verbosity > 0) {
      message('Estimating parameters for following non-regularized variables: ', latent_var_nonreg)
    }
    if (!is.null(batch_var)) {
      model_str_nonreg <- paste0('y ~ (', paste(latent_var_nonreg, collapse = ' + '), ') : ', batch_var, ' + ', batch_var, ' + 0')
    } else {
      model_str_nonreg <- paste0('y ~ ', paste(latent_var_nonreg, collapse = ' + '))
    }

    times$get_model_pars_nonreg = Sys.time()
    model_pars_nonreg <- get_model_pars_nonreg(genes, bin_size, model_pars_fit, regressor_data, umi, model_str_nonreg, cell_attr, verbosity)

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

  times$get_residuals = Sys.time()
  if (!residual_type == 'none') {
    if (verbosity > 0) {
      message('Second step: Get residuals using fitted parameters for ', length(x = genes), ' genes')
    }
    bin_ind <- ceiling(x = 1:length(x = genes) / bin_size)
    max_bin <- max(bin_ind)
    if (verbosity > 1) {
      pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
    }
    res <- matrix(NA_real_, length(genes), nrow(regressor_data_final), dimnames = list(genes, rownames(regressor_data_final)))
    for (i in 1:max_bin) {
      genes_bin <- genes[bin_ind == i]
      mu <- exp(tcrossprod(model_pars_final[genes_bin, -1, drop=FALSE], regressor_data_final))
      y <- as.matrix(umi[genes_bin, , drop=FALSE])
      res[genes_bin, ] <- switch(residual_type,
        'pearson' = pearson_residual(y, mu, model_pars_final[genes_bin, 'theta'], min_var = min_variance),
        'deviance' = deviance_residual(y, mu, model_pars_final[genes_bin, 'theta']),
        stop('residual_type ', residual_type, ' unknown - only pearson and deviance supported at the moment')
      )
      if (verbosity > 1) {
        setTxtProgressBar(pb, i)
      }
    }
    if (verbosity > 1) {
      close(pb)
    }
  } else {
    if (verbosity > 0) {
      message('Skip calculation of full residual matrix')
    }
    res <- matrix(data = NA, nrow = 0, ncol = 0)
  }

  rv <- list(y = res,
             model_str = model_str,
             model_pars = model_pars,
             model_pars_outliers = model_pars_outliers,
             model_pars_fit = model_pars_fit,
             model_str_nonreg = model_str_nonreg,
             model_pars_nonreg = model_pars_nonreg,
             arguments = arguments,
             genes_log_gmean_step1 = genes_log_gmean_step1,
             cells_step1 = cells_step1,
             cell_attr = cell_attr)
  rm(res)
  gc(verbose = FALSE)

  times$correct_umi = Sys.time()
  if (return_corrected_umi) {
    if (residual_type != 'pearson') {
      message("Will not return corrected UMI because residual type is not set to 'pearson'")
    } else {
      rv$umi_corrected <- sctransform::correct(rv, do_round = TRUE, do_pos = TRUE,
                                               verbosity = verbosity)
      rv$umi_corrected <- as(object = rv$umi_corrected, Class = 'dgCMatrix')
    }
  }

  rv$y[rv$y < res_clip_range[1]] <- res_clip_range[1]
  rv$y[rv$y > res_clip_range[2]] <- res_clip_range[2]

  if (!return_cell_attr) {
    rv[['cell_attr']] <- NULL
  }

  times$get_gene_attr = Sys.time()
  if (return_gene_attr) {
    if (verbosity > 0) {
      message('Calculating gene attributes')
    }
    gene_attr <- data.frame(
      detection_rate = genes_cell_count[genes] / ncol(umi),
      gmean = 10 ^ genes_log_gmean,
      variance = row_var(umi))
    if (ncol(rv$y) > 0) {
      gene_attr$residual_mean = rowMeans(rv$y)
      gene_attr$residual_variance = row_var(rv$y)
    }
    rv[['gene_attr']] <- gene_attr
  }

  if (verbosity > 0) {
    message('Wall clock passed: ', capture.output(print(Sys.time() - times$start_time)))
  }
  times$done = Sys.time()
  rv$times <- times
  return(rv)
}


get_model_pars <- function(genes_step1, bin_size, umi, model_str, cells_step1,
                           method, data_step1, theta_given, theta_estimation_fun,
                           verbosity) {
  bin_ind <- ceiling(x = 1:length(x = genes_step1) / bin_size)
  max_bin <- max(bin_ind)
  if (verbosity > 0) {
    message('Get Negative Binomial regression parameters per gene')
    message('Using ', length(x = genes_step1), ' genes, ', length(x = cells_step1), ' cells')
  }

  if (verbosity > 1) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  model_pars <- list()
  for (i in 1:max_bin) {
    genes_bin_regress <- genes_step1[bin_ind == i]
    umi_bin <- as.matrix(umi[genes_bin_regress, cells_step1, drop=FALSE])
    if (!is.null(theta_given)) {
      theta_given_bin <- theta_given[genes_bin_regress]
    }

    # umi_bin is a matrix of counts - we want a model per row
    # if there are multiple workers, split up the matrix in chunks of n rows
    # where n is the number of workers
    n_workers <- 1
    if (future::supportsMulticore()) {
      n_workers <- future::nbrOfWorkers()
    }
    genes_per_worker <- nrow(umi_bin) / n_workers + .Machine$double.eps
    index_vec <- 1:nrow(umi_bin)
    index_lst <- split(index_vec, ceiling(index_vec/genes_per_worker))

    # the index list will have at most n_workers entries, each one defining which genes to work on
    par_lst <- future_lapply(
      X = index_lst,
      FUN = function(indices) {
        umi_bin_worker <- umi_bin[indices, , drop = FALSE]
        if (method == 'poisson') {
          return(fit_poisson(umi = umi_bin_worker, model_str = model_str, data = data_step1, theta_estimation_fun = theta_estimation_fun))
        }
        if (method == 'qpoisson') {
          return(fit_qpoisson(umi = umi_bin_worker, model_str = model_str, data = data_step1))
        }
        if (method == 'nb_theta_given') {
          theta_given_bin_worker <- theta_given_bin[indices]
          return(fit_nb_theta_given(umi = umi_bin_worker, model_str = model_str, data = data_step1, theta_given = theta_given_bin_worker))
        }
        if (method == 'nb_fast') {
          return(fit_nb_fast(umi = umi_bin_worker, model_str = model_str, data = data_step1, theta_estimation_fun = theta_estimation_fun))
        }
        if (method == 'nb') {
          return(fit_nb(umi = umi_bin_worker, model_str = model_str, data = data_step1))
        }
        if (method == "glmGamPoi") {
          return(fit_glmGamPoi(umi = umi_bin_worker, model_str = model_str, data = data_step1))
        }
      }
    )
    model_pars[[i]] <- do.call(rbind, par_lst)

    if (verbosity > 1) {
      setTxtProgressBar(pb, i)
    }
  }
  model_pars <- do.call(rbind, model_pars)
  if (verbosity > 1) {
    close(pb)
  }
  rownames(model_pars) <- genes_step1
  colnames(model_pars)[1] <- 'theta'
  return(model_pars)
}

get_model_pars_nonreg <- function(genes, bin_size, model_pars_fit, regressor_data, umi, model_str_nonreg, cell_attr, verbosity) {
  bin_ind <- ceiling(x = 1:length(x = genes) / bin_size)
  max_bin <- max(bin_ind)
  if (verbosity > 1) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  model_pars_nonreg <- list()
  for (i in 1:max_bin) {
    genes_bin <- genes[bin_ind == i]
    mu <- tcrossprod(model_pars_fit[genes_bin, -1, drop=FALSE], regressor_data)
    umi_bin <- as.matrix(umi[genes_bin, ])
    model_pars_nonreg[[i]] <- do.call(rbind,
                                      future_lapply(genes_bin, function(gene) {
                                        fam <- negative.binomial(theta = model_pars_fit[gene, 'theta'], link = 'log')
                                        y <- umi_bin[gene, ]
                                        offs <- mu[gene, ]
                                        fit <- glm(as.formula(model_str_nonreg), data = cell_attr, family = fam, offset=offs)
                                        return(fit$coefficients)
                                      }))
    if (verbosity > 1) {
      setTxtProgressBar(pb, i)
    }
  }
  if (verbosity > 1) {
    close(pb)
  }
  model_pars_nonreg <- do.call(rbind, model_pars_nonreg)
  rownames(model_pars_nonreg) <- genes
  return(model_pars_nonreg)
}

reg_model_pars <- function(model_pars, genes_log_gmean_step1, genes_log_gmean, cell_attr,
                           batch_var, cells_step1, genes_step1, umi, bw_adjust, gmean_eps,
                           theta_regularization, verbosity) {
  genes <- names(genes_log_gmean)

  # we don't regularize theta directly
  # prior to v0.3 we regularized log10(theta)
  # now we transform to overdispersion factor
  # variance of NB is mu * (1 + mu / theta)
  # (1 + mu / theta) is what we call overdispersion factor here
  dispersion_par <- switch(theta_regularization,
    'log_theta' = log10(model_pars[, 'theta']),
    'od_factor' = log10(1 + 10^genes_log_gmean_step1 / model_pars[, 'theta']),
    stop('theta_regularization ', theta_regularization, ' unknown - only log_theta and od_factor supported at the moment')
  )

  model_pars <- model_pars[, colnames(model_pars) != 'theta']
  model_pars <- cbind(dispersion_par, model_pars)

  # look for outliers in the parameters
  # outliers are those that do not fit the overall relationship with the mean at all
  outliers <- apply(model_pars, 2, function(y) is_outlier(y, genes_log_gmean_step1))
  outliers <- apply(outliers, 1, any)
  if (sum(outliers) > 0) {
    if (verbosity > 0) {
      message('Found ', sum(outliers), ' outliers - those will be ignored in fitting/regularization step\n')
    }
    model_pars <- model_pars[!outliers, ]
    genes_step1 <- rownames(model_pars)
    genes_log_gmean_step1 <- genes_log_gmean_step1[!outliers]
  }

  # select bandwidth to be used for smoothing
  bw <- bw.SJ(genes_log_gmean_step1) * bw_adjust

  # for parameter predictions
  x_points <- pmax(genes_log_gmean, min(genes_log_gmean_step1))
  x_points <- pmin(x_points, max(genes_log_gmean_step1))

  # take results from step 1 and fit/predict parameters to all genes
  o <- order(x_points)
  model_pars_fit <- matrix(NA_real_, length(genes), ncol(model_pars),
                           dimnames = list(genes, colnames(model_pars)))

  # fit / regularize dispersion parameter
  model_pars_fit[o, 'dispersion_par'] <- ksmooth(x = genes_log_gmean_step1, y = model_pars[, 'dispersion_par'],
                                                 x.points = x_points, bandwidth = bw, kernel='normal')$y

  if (is.null(batch_var)){
    # global fit / regularization for all coefficients
    for (i in 2:ncol(model_pars)) {
      model_pars_fit[o, i] <- ksmooth(x = genes_log_gmean_step1, y = model_pars[, i],
                                      x.points = x_points, bandwidth = bw, kernel='normal')$y
    }
  } else {
    # fit / regularize per batch
    batches <- unique(cell_attr[, batch_var])
    for (b in batches) {
      sel <- cell_attr[, batch_var] == b & rownames(cell_attr) %in% cells_step1
      #batch_genes_log_gmean_step1 <- log10(rowMeans(umi[genes_step1, sel]))
      batch_genes_log_gmean_step1 <- log10(row_gmean(umi[genes_step1, sel], eps = gmean_eps))
      if (any(is.infinite(batch_genes_log_gmean_step1))) {
        if (verbosity > 0) {
          message('Some genes not detected in batch ', b, ' -- assuming a low mean.')
        }
        batch_genes_log_gmean_step1[is.infinite(batch_genes_log_gmean_step1) & batch_genes_log_gmean_step1 < 0] <- min(batch_genes_log_gmean_step1[!is.infinite(batch_genes_log_gmean_step1)])
      }
      sel <- cell_attr[, batch_var] == b
      #batch_genes_log_gmean <- log10(rowMeans(umi[, sel]))
      batch_genes_log_gmean <- log10(row_gmean(umi[, sel], eps = gmean_eps))
      # in case some genes have not been observed in this batch
      batch_genes_log_gmean <- pmax(batch_genes_log_gmean, min(batch_genes_log_gmean_step1))
      batch_o <- order(batch_genes_log_gmean)
      for (i in which(grepl(paste0(batch_var, b), colnames(model_pars)))) {
        model_pars_fit[batch_o, i] <- ksmooth(x = batch_genes_log_gmean_step1, y = model_pars[, i],
                                              x.points = batch_genes_log_gmean, bandwidth = bw, kernel='normal')$y
      }
    }
  }

  # back-transform dispersion parameter to theta
  theta <- switch(theta_regularization,
    'log_theta' = 10^model_pars_fit[, 'dispersion_par'],
    'od_factor' = 10^genes_log_gmean / (10^model_pars_fit[, 'dispersion_par'] - 1)
  )
  model_pars_fit <- model_pars_fit[, colnames(model_pars_fit) != 'dispersion_par']
  model_pars_fit <- cbind(theta, model_pars_fit)

  attr(model_pars_fit, 'outliers') <- outliers
  return(model_pars_fit)
}
