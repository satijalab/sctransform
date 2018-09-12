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
