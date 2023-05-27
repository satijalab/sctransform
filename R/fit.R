# Fit NB regression models using different approaches

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
                                       do_cox_reid_adjustment = TRUE,
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
  log10_umi <- data$log_umi
  stopifnot(!is.null(log10_umi))
  log_umi <- log(10^log10_umi)

  # only intercept varies
  includes.batch_var <- FALSE
  if (grepl(pattern = "\\(log_umi\\) :", x = model_str)) {
    includes.batch_var <- TRUE
  }
  new_formula <- gsub("y", "", model_str)

  includes.log_umi <- grepl(pattern = "~ log_umi", x = new_formula)
  if (!includes.batch_var) {
    # if therse is no batch variable - remove log_umi and fix it
    # remove log_umi from model formula if it is with batch variables
    new_formula <- gsub(pattern = "\\+ log_umi", replacement = "", x = new_formula)
    # replace log_umi with 1 if it is the only formula

    new_formula <- gsub(pattern = "log_umi", replacement = "1", x = new_formula)
  } else {
    log_umi <- 0
  }
  fit <- glmGamPoi::glm_gp(data = umi,
                           design = as.formula(new_formula),
                           col_data = data,
                           offset = log_umi,
                           size_factors = FALSE)
  fit$theta <- 1 / fit$overdispersions
  if (!allow_inf_theta){
    fit$theta <- pmin(1 / fit$overdispersions, rowMeans(fit$Mu) / 1e-4)
  }
  if ("Intercept" %in% colnames(x = fit$Beta)) {
    if (includes.log_umi){
      model_pars <- cbind(fit$theta,
                          fit$Beta[, "Intercept"],
                          rep(log(10), nrow(umi)))
      dimnames(model_pars) <- list(rownames(umi), c('theta', '(Intercept)', 'log_umi'))

      n_coefficients <- ncol(fit$Beta)
      if (n_coefficients>1){
        model_pars <- cbind(model_pars, fit$Beta[, 2:n_coefficients])
        colnames(x = model_pars)[4:ncol(x = model_pars)] <- colnames(x = fit$Beta)[2:n_coefficients]
      }
    } else {
      model_pars <- cbind(fit$theta,
                          fit$Beta)
      dimnames(model_pars) <- list(rownames(umi), c('theta', colnames(x = fit$Beta)))
    }
  } else {
    if (!includes.batch_var){
      if (includes.log_umi){
        model_pars <- cbind(fit$theta,
                            rep(log(10), nrow(umi)),
                            fit$Beta)
        dimnames(model_pars) <- list(rownames(umi), c('theta', 'log_umi', colnames(x = fit$Beta)))
      } else {
        model_pars <- cbind(fit$theta,
                            fit$Beta)
        dimnames(model_pars) <- list(rownames(umi), c('theta', colnames(x = fit$Beta)))
      }
    } else {
      model_pars <- cbind(fit$theta,
                          fit$Beta)
      dimnames(model_pars) <- list(rownames(umi), c('theta', colnames(x = fit$Beta)))
    }
  }
  colnames(x = model_pars)[match(x = 'Intercept', table = colnames(x = model_pars))] <- "(Intercept)"
  return(model_pars)
}

fit_nb_offset <- function(umi, model_str, data, allow_inf_theta=FALSE) {
  # remove log_umi from model formula if it is with batch variables
  new_formula <- gsub("\\+ log_umi", "", model_str)
  # replace log_umi with 1 if it is the only formula
  new_formula <- gsub("log_umi", "1 + offset(log_umi)", new_formula)

  log10_umi <- data$log_umi
  stopifnot(!is.null(log10_umi))
  log_umi <- log(10^log10_umi)
  data$log_umi <- log_umi

  par_mat <- apply(umi, 1, function(y) {
    fit <- 0
    try(fit <- glm.nb(formula = as.formula(new_formula), data = data), silent=TRUE)
    if (inherits(x = fit, what = 'numeric')) {
      fit <- glm(formula = as.formula(new_formula), data = data, family = poisson)
      fit$theta <- Inf
      #fit$theta <- as.numeric(x = suppressWarnings(theta.ml(y = y, mu = fit$fitted)))
    }
    if (!allow_inf_theta){
      fit$theta <- pmin(fit$theta, mean(y) / 1e-4)
    }
    return(c(fit$theta, fit$coefficients))
  })
  model_pars <- t(par_mat)
  model_pars <- cbind(model_pars, rep(log(10), nrow(umi)))
  rownames(x = model_pars) <- rownames(x = umi)
  colnames(x = model_pars)[match(x = 'Intercept', table = colnames(x = model_pars))] <- "(Intercept)"
  return(model_pars)
}
