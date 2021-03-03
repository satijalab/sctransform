#' Peripheral Blood Mononuclear Cells (PBMCs)
#'
#' UMI counts for a subset of cells freely available from 10X Genomics
#'
#' @format A sparse matrix (dgCMatrix, see Matrix package) of molecule counts.
#' There are 914 rows (genes) and 283 columns (cells). This is a downsampled
#' version of a 3K PBMC dataset available from 10x Genomics.
#'
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k}
"pbmc"

#' Transformation functions for umify
#' 
#' The functions have been trained on various public data sets and relate quantile
#' values to log-counts. Here the expected values at various points are given.
#' 
#' @format A list of length two. The first element is a data frame with group, quantile and 
#' log-counts values. The second element is a vector of breaks to be used with cut to group
#' observations.
"umify_data"
