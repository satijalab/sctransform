% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_nz_median2}
\alias{get_nz_median2}
\title{Get median of non zero UMIs from a count matrix}
\usage{
get_nz_median2(umi, genes = NULL)
}
\arguments{
\item{umi}{Count matrix}

\item{genes}{A vector of genes to consider for calculating
the median. Default is NULL which uses all genes.}
}
\value{
A numeric value representing the median of non-zero entries from the UMI matrix
}
\description{
Get median of non zero UMIs from a count matrix
}
