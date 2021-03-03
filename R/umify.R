data("umify_data", envir=environment())

#' Quantile normalization of cell-level data to match typical UMI count data
#' 
#' @param counts A matrix of class dgCMatrix with genes as rows and columns as cells
#' 
#' @return A UMI-fied count matrix
#' 
#' @section Details:
#' sctransform::vst operates under the assumption that gene counts approximately 
#' follow a Negative Binomial dristribution. For UMI-based data that seems to be 
#' the case, however, non-UMI data does not behave in the same way. 
#' In some cases it might be better to to apply a transformation to such data 
#' to make it look like UMI data. This function applies such a transformation function.
#'
#' Cells in the input matrix are processed independently. For each cell
#' the non-zero data is transformed to quantile values. Based on the number of genes 
#' detected a smooth function is used to predict the UMI-like counts.
#' 
#' The functions have be trained on various public data sets and come as part of the 
#' package (see umify_data data set in this package).
#' 
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise pull mutate case_when
#' @importFrom stats approxfun
#' @importFrom rlang .data
#' 
#' @export
#' 
#' @examples
#' \donttest{
#' silly_example <- umify(pbmc)
#' }
umify <- function(counts) {
  # check input
  if (!inherits(x = counts, what = 'dgCMatrix')) {
    stop('counts must be a dgCMatrix')
  }
  
  # load the group breaks needed to place cells into groups
  grp_breaks <- umify_data$grp_breaks
  K <- length(grp_breaks) - 1
  w <- mean(diff(grp_breaks))
  
  # given the umify data models, create functions that
  # predict the log count from the distribution quantile  
  # create one function per group membership (based on number of genes detected)
  apprx_funs <- group_by(umify_data$fit_df, .data$grp) %>%
    summarise(fun = list(approxfun(x = .data$q, y = .data$log_y, rule = 2:1)), .groups = 'drop') %>% 
    pull(.data$fun)
  names(apprx_funs) <- levels(umify_data$fit_df$grp)
  
  # for each cell in the input we need to know how many genes are detected
  # then determine the primary and secondary group and the weights
  # to be used for the linear interpolation
  ca <- data.frame(genes = diff(counts@p)) %>%
    mutate(log_genes = log10(.data$genes),
           grp1 = cut(.data$log_genes, breaks = grp_breaks, right = FALSE),
           weight = ((.data$log_genes - grp_breaks[1]) / w) %% 1,
           grp2 = case_when(.data$weight >= 0.5 ~ as.numeric(.data$grp1)+1, TRUE ~ as.numeric(.data$grp1)-1),
           grp2 = case_when(.data$grp2 < 1 ~ 1, .data$grp2 > K ~ K, TRUE ~ .data$grp2),
           grp2 = factor(levels(.data$grp1)[.data$grp2], levels = levels(.data$grp1)),
           weight = -abs(.data$weight-0.5)+1)
  
  if (any(is.na(ca$grp1))) {
    warning(sprintf('Cells with very few or too many genes detected. 
                     The lower limit is %d, the upper limit is %d. 
                     Setting non-zero values to NA for %d cells.', 
                    floor(10^min(grp_breaks)), 
                    floor(10^max(grp_breaks)), 
                    sum(is.na(ca$grp1))))
  }
  
  # predict UMIfied counts
  # per cell
  out_vec <- rep(NA_real_, length(counts@x))
  j <- 1
  while (j <= ncol(counts)) {
    if (!is.na(ca$grp1[j])) {
      i_start <- counts@p[j] + 1
      i_end <- counts@p[j+1]
      sel <- i_start:i_end
      q <- rank(counts@x[sel], ties.method = 'max') / length(sel)
      pred1 <- apprx_funs[[ca$grp1[j]]](q)
      pred2 <- apprx_funs[[ca$grp2[j]]](q)
      out_vec[sel] <- round(10^(pred1 * ca$weight[j] + pred2 * (1 - ca$weight[j])))
    }
    j <- j + 1
  }
  
  counts_new <- sparseMatrix(i = counts@i, 
                             p = counts@p,
                             x = out_vec,
                             dims = counts@Dim,
                             dimnames = counts@Dimnames,
                             giveCsparse = TRUE,
                             index1 = FALSE)
  return(counts_new)
}
