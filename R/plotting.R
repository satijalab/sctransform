#' Plot estimated and fitted model parameters
#'
#' @param vst_out The output of a vst run
#' @param xaxis Variable to plot on X axis; default is "gmean"
#' @param show_theta Whether to show the theta parameter; default is FALSE (only the overdispersion factor is shown)
#' @param show_var Whether to show the average model variance; default is FALSE
#' @param verbosity An integer specifying whether to show only messages (1), messages and progress bars (2) or nothing (0) while the function is running; default is 2
#' @param verbose Deprecated; use verbosity instead
#' @param show_progress Deprecated; use verbosity instead
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @import reshape2
#'
#' @export
#'
#' @examples
#' \donttest{
#' vst_out <- vst(pbmc, return_gene_attr = TRUE)
#' plot_model_pars(vst_out)
#' }
#'
plot_model_pars <- function(vst_out, xaxis="gmean", show_theta = FALSE, show_var = FALSE,
                            verbosity = 2, verbose = NULL, show_progress = NULL) {
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

  if (! 'gmean' %in% names(vst_out$gene_attr)) {
    stop('vst_out must contain a data frame named gene_attr with a column named gmean (perhaps call vst with return_gene_attr = TRUE)')
  }

  # first handle the per-gene estimates
  mp <- get_model_par_mat(vst_out, model_pars = vst_out$model_pars, use_nonreg = TRUE,
                          show_theta = show_theta, show_var = show_var,
                          verbosity = verbosity)
  ordered_par_names <- colnames(mp)

  # second the regularized estimates
  mp_fit <- get_model_par_mat(vst_out, model_pars = vst_out$model_pars_fit, use_nonreg = FALSE,
                              show_theta = show_theta, show_var = show_var,
                              verbosity = verbosity)

  mpnr <- vst_out$model_pars_nonreg
  if (!is.null(dim(mpnr))) {
    colnames(mpnr) <- paste0('nonreg:', colnames(mpnr))
    mp <- cbind(mp, mpnr)
    ordered_par_names <- c(ordered_par_names, colnames(mpnr))
  }
  # show estimated and regularized parameters
  df <- melt(mp, varnames = c('gene', 'parameter'), as.is = TRUE)
  df$parameter <- factor(df$parameter, levels = ordered_par_names)
  df_fit <- melt(mp_fit, varnames = c('gene', 'parameter'), as.is = TRUE)
  df_fit$parameter <- factor(df_fit$parameter, levels = ordered_par_names)
  df$is_outl <- vst_out$model_pars_outliers
  df$type <- 'single gene estimate'
  df_fit$type <- 'regularized'
  df_fit$is_outl <- FALSE

  if (startsWith(x = vst_out$arguments$method, prefix = 'offset') | (xaxis=="amean")) {
    df$x <- vst_out$gene_attr[df$gene, 'amean']
    df_fit$x <- vst_out$gene_attr[df_fit$gene, 'amean']
    xlab <- 'Arithmetic mean of gene [log10]'
  } else {
    df$x <- vst_out$gene_attr[df$gene, 'gmean']
    df_fit$x <- vst_out$gene_attr[df_fit$gene, 'gmean']
    xlab <- 'Geometric mean of gene [log10]'
  }

  df_plot <- rbind(df, df_fit)
  df_plot$parameter <- factor(df_plot$parameter, levels = ordered_par_names)

  if (!vst_out$arguments$do_regularize || startsWith(x = vst_out$arguments$method, prefix = 'offset')) {
    df_plot <- df_plot[df_plot$type == 'single gene estimate', ]
    legend_pos <- 'none'
  } else {
    legend_pos <- 'bottom'
  }

  g <- ggplot(df_plot, aes_(x=~log10(x), y=~value, color=~type)) +
    geom_point(data=df, aes_(shape=~is_outl), size=0.5, alpha=0.5) +
    scale_shape_manual(values=c(16, 4), guide = FALSE) +
    geom_point(data=df_fit, size=0.66, alpha=0.5, shape=16) +
    facet_wrap(~ parameter, scales = 'free_y', ncol = ncol(mp)) +
    theme(legend.position = legend_pos) +
    xlab(label = xlab)
  return(g)
}

# helper function to plot model parameters
get_model_par_mat <- function(vst_out, model_pars, use_nonreg, show_theta = FALSE, show_var = FALSE,
                              verbosity = 2) {
  mp <- model_pars
  # transform theta to overdispersion factor
  if (startsWith(x = vst_out$arguments$method, prefix = 'offset')) {
    mp[, 1] <- log10(1 + vst_out$gene_attr[rownames(mp), 'amean'] / mp[, 'theta'])
  } else {
    mp[, 1] <- log10(1 + vst_out$gene_attr[rownames(mp), 'gmean'] / mp[, 'theta'])
  }

  colnames(mp)[1] <- 'log10(od_factor)'
  ordered_par_names <- colnames(mp)[c(2:ncol(mp), 1)]
  if (show_theta) {
    mp <- cbind(mp, log10(model_pars[, 'theta']))
    colnames(mp)[ncol(mp)] <- 'log10(theta)'
    ordered_par_names <- c(ordered_par_names, 'log10(theta)')
  }
  if (show_var) {
    mp <- cbind(mp, log10(get_model_var(vst_out, use_nonreg = use_nonreg, verbosity = verbosity)))
    colnames(mp)[ncol(mp)] <- 'log10(model var)'
    ordered_par_names <- c(ordered_par_names, 'log10(model var)')
  }
  return(mp[, ordered_par_names])
}


# helper function to plot model fit for a single gene
# returns list with mean, sd, pearson residual
#' @importFrom stats model.matrix
get_nb_fit <- function(x, umi, gene, cell_attr, as_poisson = FALSE) {
  regressor_data <- model.matrix(as.formula(gsub('^y', '', x$model_str)), cell_attr)

  coefs <- x$model_pars_fit[gene, -1, drop=FALSE]
  theta <- x$model_pars_fit[gene, 1]
  if (as_poisson) {
    theta <- Inf
  }
  mu <- exp(coefs %*% t(regressor_data))[1, ]
  sd <- sqrt(mu + mu^2 / theta)
  res <- (umi[gene, ] - mu) / sd
  res <- pmin(res, x$arguments$res_clip_range[2])
  res <- pmax(res, x$arguments$res_clip_range[1])
  ret_df <- data.frame(mu = mu, sd = sd, res = res)
  # in case we have individual (non-regularized) parameters
  if (gene %in% rownames(x$model_pars)) {
    coefs <- x$model_pars[gene, -1, drop=FALSE]
    theta <- x$model_pars[gene, 1]
    ret_df$mu_nr <- exp(coefs %*% t(regressor_data))[1, ]
    ret_df$sd_nr <- sqrt(ret_df$mu_nr + ret_df$mu_nr^2 / theta)
    ret_df$res_nr <- (umi[gene, ] - ret_df$mu_nr) / ret_df$sd_nr
    ret_df$res_nr <- pmin(ret_df$res_nr, x$arguments$res_clip_range[2])
    ret_df$res_nr <- pmax(ret_df$res_nr, x$arguments$res_clip_range[1])
  } else {
    ret_df$mu_nr <- NA_real_
    ret_df$sd_nr <- NA_real_
    ret_df$res_nr <- NA_real_
  }
  return(ret_df)
}

#' Plot observed UMI counts and model
#'
#' @param x The output of a vst run
#' @param umi UMI count matrix
#' @param goi Vector of genes to plot
#' @param x_var Cell attribute to use on x axis; will be taken from x$arguments$latent_var[1] by default
#' @param cell_attr Cell attributes data frame; will be taken from x$cell_attr by default
#' @param do_log Log10 transform the UMI counts in plot
#' @param show_fit Show the model fit
#' @param show_nr Show the non-regularized model (if available)
#' @param plot_residual Add panels for the Pearson residuals
#' @param batches Manually specify a batch variable to break up the model plot in segments
#' @param as_poisson Fix model parameter theta to Inf, effectively showing a Poisson model
#' @param arrange_vertical Stack individual ggplot objects or place side by side
#' @param show_density Draw 2D density lines over points
#' @param gg_cmds Additional ggplot layer commands
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @import reshape2
#' @importFrom gridExtra grid.arrange
#' @importFrom stats bw.nrd0
#'
#' @export
#'
#' @examples
#' \donttest{
#' vst_out <- vst(pbmc, return_cell_attr = TRUE)
#' plot_model(vst_out, pbmc, 'EMC4')
#' }
#'
plot_model <- function(x, umi, goi, x_var = x$arguments$latent_var[1], cell_attr = x$cell_attr,
                       do_log = TRUE, show_fit = TRUE, show_nr = FALSE, plot_residual = FALSE,
                       batches = NULL, as_poisson = FALSE, arrange_vertical = TRUE,
                       show_density = FALSE, gg_cmds = NULL) {
  if (is.null(batches)) {
    if (!is.null(x$arguments$batch_var)) {
      batches <- cell_attr[, x$arguments$batch_var]
    } else {
      batches <- rep(1, nrow(cell_attr))
    }
  }
  df_list <- list()
  for (gene in goi) {
    nb_fit <- get_nb_fit(x, umi, gene, cell_attr, as_poisson)
    nb_fit$x <- cell_attr[, x_var]
    nb_fit$y <- umi[gene, ]
    nb_fit$batch <- batches
    nb_fit$gene <- gene
    nb_fit$ymin <- nb_fit$mu - nb_fit$sd
    nb_fit$ymax <- nb_fit$mu + nb_fit$sd
    nb_fit$ymin_nr <- nb_fit$mu_nr - nb_fit$sd_nr
    nb_fit$ymax_nr <- nb_fit$mu_nr + nb_fit$sd_nr
    if (do_log) {
      nb_fit$y <- log10(nb_fit$y + 1)
      nb_fit$mu <- log10(nb_fit$mu + 1)
      nb_fit$mu_nr <- log10(nb_fit$mu_nr + 1)
      nb_fit$ymin <- log10(pmax(nb_fit$ymin, 0) + 1)
      nb_fit$ymax <- log10(pmax(nb_fit$ymax, 0) + 1)
      nb_fit$ymin_nr <- log10(pmax(nb_fit$ymin_nr, 0) + 1)
      nb_fit$ymax_nr <- log10(pmax(nb_fit$ymax_nr, 0) + 1)
    }
    df_list[[length(df_list) + 1]] <- nb_fit[order(nb_fit$x), ]
  }
  df <- do.call(rbind, df_list)
  df$gene <- factor(df$gene, ordered=TRUE, levels=unique(df$gene))
  g <- ggplot(df, aes_(~x, ~y)) + geom_point(alpha=0.5, shape=16)
  if (show_density) {
    bandwidths <- c(bw.nrd0(g$data$x), bw.nrd0(g$data$y))
    g <- g + geom_density_2d(color = 'lightblue', size=0.5, h = bandwidths,
                             contour_var = "ndensity")
  }
  if (show_fit) {
    for (b in unique(df$batch)) {
      g <- g +
        geom_line(data = df[df$batch == b, ], aes_(~x, ~mu), color='deeppink', size = 1) +
        geom_ribbon(data = df[df$batch == b, ], aes_(x = ~x, ymin = ~ymin, ymax = ~ymax), alpha = 0.5, fill='deeppink')
    }
  }
  if (show_nr) {
    for (b in unique(df$batch)) {
      g <- g +
        geom_line(aes_(~x, ~mu_nr), color='blue', size = 1) +
        geom_ribbon(aes_(x = ~x, ymin = ~ymin_nr, ymax = ~ymax_nr), alpha = 0.5, fill='blue')
    }
  }
  g <- g + facet_grid(~gene) + xlab(paste('Cell', x_var)) + ylab('Gene UMI counts')
  if (do_log) {
    g <- g + ylab('Gene log10(UMI + 1)')
  }
  g <- g + gg_cmds
  if (length(goi) == 1) {
    g <- g + theme(strip.text = element_blank())
  }

  if (plot_residual) {
    ga_col = 1
    res_range <- range(df$res)
    g2 <- ggplot(df, aes_(~x, ~res)) + geom_point(alpha = 0.5, shape=16) +
      coord_cartesian(ylim = res_range) +
      facet_grid(~gene) + xlab(x) + ylab('Pearson residual') +
      xlab(paste('Cell', x_var)) + gg_cmds +
      theme(strip.text = element_blank()) # strip.background = element_blank(),
    if (show_density) {
      g2 <- g2 + geom_density_2d(color = 'lightblue', size=0.5)
    }
    if (show_nr) {
      g3 <- ggplot(df, aes_(~x, ~res_nr)) + geom_point(alpha = 0.5, shape=16) +
        coord_cartesian(ylim = res_range) +
        facet_grid(~gene) + xlab(x) + ylab('Pearson residual non-reg.') +
        xlab(paste('Cell', x_var)) + gg_cmds +
        theme(strip.text = element_blank()) # strip.background = element_blank(),
      if (show_density) {
        g3 <- g3 + geom_density_2d(color = 'lightblue', size=0.5)
      }
      if (!arrange_vertical) {
        ga_col = 3
      }
      return(grid.arrange(g, g2, g3, ncol=ga_col))
    } else {
      if (!arrange_vertical) {
        ga_col = 2
      }
      return(grid.arrange(g, g2, ncol=ga_col))
    }
  }

  return(g)
}

