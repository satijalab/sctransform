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
