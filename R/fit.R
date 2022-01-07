# Fir NB regression models using different approaches

fit_poisson <- function(umi, model_str, data, theta_estimation_fun) {
  regressor_data <- model.matrix(as.formula(gsub('^y', '', model_str)), data)
  dfr <- ncol(umi) - ncol(regressor_data)
  par_mat <- t(apply(umi, 1, function(y) {
    fit <- qpois_reg(regressor_data, y, 1e-9, 100, 1.0001, TRUE)
    theta <- switch(theta_estimation_fun,
      'theta.ml' = as.numeric(x = suppressWarnings(theta.ml(y = y, mu = fit$fitted))),
      'theta.mm' = as.numeric(x = theta.mm(y = y, mu = fit$fitted, dfr = dfr)),
      stop('theta_estimation_fun ', theta_estimation_fun, ' unknown - only theta.ml and theta.mm supported at the moment')
    )
    return(c(theta, fit$coefficients))
  }))
  return(par_mat)
}

fit_qpoisson <- function(umi, model_str, data) {
  regressor_data <- model.matrix(as.formula(gsub('^y', '', model_str)), data)
  par_mat <- t(apply(umi, 1, function(y) {
    fit <- qpois_reg(regressor_data, y, 1e-9, 100, 1.0001, FALSE)
    return(c(fit$theta.guesstimate, fit$coefficients))
  }))
  return(par_mat)
}

fit_nb_theta_given <- function(umi, model_str, data, theta_given) {
  par_lst <- lapply(1:nrow(umi), function(j) {
    y <- umi[j, ]
    theta <- theta_given[j]
    fit2 <- 0
    try(fit2 <- glm(as.formula(model_str), data = data, family = negative.binomial(theta=theta)), silent=TRUE)
    if (inherits(x = fit2, what = 'numeric')) {
      return(c(theta, glm(as.formula(model_str), data = data, family = poisson)$coefficients))
    } else {
      return(c(theta, fit2$coefficients))
    }
  })
  return(do.call(rbind, par_lst))
}

fit_nb_fast <- function(umi, model_str, data, theta_estimation_fun) {
  regressor_data <- model.matrix(as.formula(gsub('^y', '', model_str)), data)
  dfr <- ncol(umi) - ncol(regressor_data)
  par_mat <- apply(umi, 1, function(y) {
    fit <- qpois_reg(regressor_data, y, 1e-9, 100, 1.0001, TRUE)
    theta <- switch(theta_estimation_fun,
                    'theta.ml' = as.numeric(x = theta.ml(y = y, mu = fit$fitted)),
                    'theta.mm' = as.numeric(x = theta.mm(y = y, mu = fit$fitted, dfr = dfr)),
                    stop('theta_estimation_fun ', theta_estimation_fun, ' unknown - only theta.ml and theta.mm supported at the moment')
    )
    fit2 <- 0
    try(fit2 <- glm(as.formula(model_str), data = data, family = negative.binomial(theta=theta)), silent=TRUE)
    if (inherits(x = fit2, what = 'numeric')) {
      return(c(theta, fit$coefficients))
    } else {
      return(c(theta, fit2$coefficients))
    }
  })
  return(t(par_mat))
}

fit_nb <- function(umi, model_str, data) {
  par_mat <- apply(umi, 1, function(y) {
    fit <- 0
    try(fit <- glm.nb(as.formula(model_str), data = data), silent=TRUE)
    if (inherits(x = fit, what = 'numeric')) {
      fit <- glm(as.formula(model_str), data = data, family = poisson)
      fit$theta <- as.numeric(x = suppressWarnings(theta.ml(y = y, mu = fit$fitted)))
    }
    return(c(fit$theta, fit$coefficients))
  })
  return(t(par_mat))
}

# allow_inf_theta: if FALSE, replace theta by min(theta, rowmeans(mu)/1e-4)
#                  else allow theta = Inf (poisson)
fit_glmGamPoi <- function(umi, model_str, data, allow_inf_theta=FALSE) {
  fit <- glmGamPoi::glm_gp(data = umi,
                           design = as.formula(gsub("y", "", model_str)),
                           col_data = data,
                           size_factors = FALSE)

  fit$theta <- 1 / fit$overdispersions
  if (!allow_inf_theta){
    fit$theta <- pmin(1 / fit$overdispersions, rowMeans(fit$Mu) / 1e-4)
  }
  colnames(fit$Beta)[match(x = 'Intercept', colnames(fit$Beta))] <- "(Intercept)"
  return(cbind(fit$theta, fit$Beta))
}

fit_overdisp_mle <- function(umi, mu, intercept, slope){
  fit <- glmGamPoi::overdispersion_mle(umi,
                                       mu,
                                       model_matrix = NULL,
                                       do_cox_reid_adjustment = TRUE, #!is.null(model_matrix),
                                       global_estimate = FALSE,
                                       subsample = FALSE,
                                       max_iter = 200,
                                       verbose = FALSE)
  theta <- 1 / fit$estimate
  model_pars <- cbind(theta, intercept, slope)
  colnames(model_pars) <- c("theta", "(Intercept)", "log_umi")
  return (model_pars)
}

# Use log_umi as offset using glmGamPoi
fit_glmGamPoi_offset <- function(umi, model_str, data,  allow_inf_theta=FALSE) {
  # only intercept varies
  new_formula <- gsub("y", "", model_str)
  # remove log_umi from model formula if it is with batch variables
  new_formula <- gsub("\\+ log_umi", "", new_formula)
  # replace log_umi with 1 if it is the only formula
  new_formula <- gsub("log_umi", "1", new_formula)

  log10_umi <- data$log_umi
  stopifnot(!is.null(log10_umi))
  log_umi <- log(10^log10_umi)


  fit <- glmGamPoi::glm_gp(data = umi,
                           design = as.formula(new_formula),
                           col_data = data,
                           offset = log_umi,
                           size_factors = FALSE)
  fit$theta <- 1 / fit$overdispersions
  if (!allow_inf_theta){
    fit$theta <- pmin(1 / fit$overdispersions, rowMeans(fit$Mu) / 1e-4)
  }
  model_pars <- cbind(fit$theta,
                      fit$Beta[, "Intercept"],
                      rep(log(10), nrow(umi)))
  dimnames(model_pars) <- list(rownames(umi), c('theta', '(Intercept)', 'log_umi'))
  return(model_pars)
}