#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
NumericVector row_mean_dgcmatrix(NumericVector x, IntegerVector i, int rows, int cols) {
  NumericVector ret(rows, 0.0);
  for (int k=0; k<x.length(); ++k) {
    ret[i[k]] += x[k];
  }
  for (int k=0; k<rows; ++k) {
    ret[k] /= cols;
  }
  return ret;
}

// [[Rcpp::export]]
NumericMatrix row_mean_grouped_dgcmatrix(NumericVector x, IntegerVector i, IntegerVector p,
                                         IntegerVector group, int groups, int rows) {
  NumericMatrix ret(rows, groups);
  IntegerVector groupsize(groups, 0);

  int col = 0;
  for (int k=0; k<x.length(); ++k) {
    while (k>=p[col]) {
      ++col;
    }
    ret(i[k], group[col-1]) += x[k];
  }

  for (int k=0; k<group.length(); ++k) {
    ++groupsize[group[k]];
  }

  for (int j=0; j<groups; ++j) {
    for (int k=0; k<rows; ++k) {
      ret(k, j) /= groupsize[j];
    }
  }
  return ret;
}

// [[Rcpp::export]]
NumericVector row_gmean_dgcmatrix(NumericVector x, IntegerVector i, int rows, int cols, double eps) {
  NumericVector ret(rows, 0.0);
  IntegerVector nzero(rows, cols);
  for (int k=0; k<x.length(); ++k) {
    ret[i[k]] += log(x[k] + eps);
    nzero[i[k]] -= 1;
  }
  for (int k=0; k<rows; ++k) {
    ret[k] = exp((ret[k] + log(eps) * nzero[k]) / cols) - eps;
  }
  return ret;
}

// [[Rcpp::export]]
NumericMatrix row_gmean_grouped_dgcmatrix(NumericVector x, IntegerVector i, IntegerVector p,
                                         IntegerVector group, int groups, int rows, double eps) {
  NumericMatrix ret(rows, groups);
  IntegerMatrix nzero(rows, groups);
  IntegerVector groupsize(groups, 0);

  for (int k=0; k<group.length(); ++k) {
    ++groupsize[group[k]];
  }

  for (int k=0; k<groups; ++k) {
    IntegerMatrix::Column col = nzero(_, k);
    col = col + groupsize[k];
  }

  int col = 0;
  for (int k=0; k<x.length(); ++k) {
    while (k>=p[col]) {
      ++col;
    }
    ret(i[k], group[col-1]) += log(x[k] + eps);
    nzero(i[k], group[col-1]) -= 1;
  }

  for (int j=0; j<groups; ++j) {
    for (int k=0; k<rows; ++k) {
      ret(k, j) = exp((ret(k, j) + log(eps) * nzero(k, j)) / groupsize[j]) - eps;
    }
  }
  return ret;
}

// [[Rcpp::export]]
NumericVector row_var_dgcmatrix(NumericVector x, IntegerVector i, int rows, int cols) {
  NumericVector rowmean(rows, 0.0);
  for (int k=0; k<x.length(); ++k) {
    rowmean[i[k]] += x[k];
  }
  for (int k=0; k<rows; ++k) {
    rowmean[k] /= cols;
  }
  NumericVector rowvar(rows, 0.0);
  IntegerVector nzero(rows, cols);
  for (int k=0; k<x.length(); ++k) {
    rowvar[i[k]] += pow(x[k] - rowmean[i[k]], 2);
    nzero[i[k]] -= 1;
  }
  for (int k=0; k<rows; ++k) {
    rowvar[k] = (rowvar[k] + (pow(rowmean[k], 2) * nzero[k])) / (cols - 1);
  }
  return rowvar;
}

// [[Rcpp::export]]
NumericVector row_var_dense_d(Eigen::Map<Eigen::MatrixXd> x) {
  NumericVector out(x.rows());
  for(int i=0; i < x.rows(); ++i){
    Eigen::ArrayXd r = x.row(i).array();
    double rowMean = r.mean();
    out[i] = (r - rowMean).square().sum() / (x.cols() - 1);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector row_var_dense_i(Eigen::Map<Eigen::MatrixXi> x) {
  NumericVector out(x.rows());
  for(int i=0; i < x.rows(); ++i){
    Eigen::ArrayXd r = (x.row(i).array()).cast<double>();
    double rowMean = r.mean();
    out[i] = (r - rowMean).square().sum() / (x.cols() - 1);
  }
  return out;
}
