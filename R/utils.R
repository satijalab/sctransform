#' Identify outliers
#'
#' @param y Dependent variable
#' @param x Independent variable
#' @param n Number of observations per bin
#' @param th Outlier score threshold
#'
#' @return Boolean vector
#'
#' @importFrom stats aggregate
#'
is_outlier <- function(y, x, n, th = 40) {
  o <- order(x)
  x <- x[o]
  y <- y[o]
  n.bins <- ceiling(length(x) / n)
  bins <- cut(rank(x, ties.method = 'first'), n.bins, ordered_result = TRUE)
  tmp <- aggregate(x = y, by = list(bin=bins), FUN = robust_scale)
  score <- rep(0, length(x))
  if (class(tmp$x) == 'list') {
    score[o] <- unlist(tmp$x)
  } else {
    score[o] <- as.numeric(t(tmp$x))
  }
  return(score > th)
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
