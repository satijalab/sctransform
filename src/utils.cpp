// [[Rcpp::depends(RcppArmadillo)]]

#include "RcppArmadillo.h"

using namespace Rcpp;


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


// The following function was taken from the Rfast package
// with kind permission from the authors.
// It has been slightly adopted for our use case here.
// [[Rcpp::export]]
List qpois_reg(NumericMatrix X, NumericVector Y, const double tol, const int maxiters, 
               const double minphi, const bool returnfit){
  const unsigned int n=X.nrow(), pcols=X.ncol(), d=pcols;
  
  arma::colvec b_old(d, arma::fill::zeros), b_new(d), L1(d), yhat(n), y(Y.begin(), n, false), m(n), phi(n);
  arma::vec unique_vals;
  arma::mat L2, x(X.begin(), n, pcols, false), x_tr(n, pcols);
  double dif;
  
  // Identify the intercept term(s) and initialize the coefficients
  for(int i=0;i<pcols;++i){
    unique_vals = arma::unique(x.unsafe_col(i));
    
    if(unique_vals.n_elem==1){
      b_old(i)=log(mean(y));
      break;
    }
    if((unique_vals.n_elem==2) & (unique_vals[0] == 0 | unique_vals[1] == 0)){
      b_old(i)=arma::as_scalar(y.t()*x.unsafe_col(i));
      b_old(i)=b_old(i)/sum(x.unsafe_col(i));
      b_old(i)=log(std::max(1e-9, b_old(i)));
    }
  }
  
  x_tr=x.t();
  int ij=2;
  
  for(dif=1.0;dif>tol;){
    yhat=x*b_old;
    m=(exp(yhat));
    phi=y-m;
    L1=x_tr*phi;
    L2=x.each_col()%m;
    L2=x_tr*L2;
    b_new=b_old+solve(L2,L1,arma::solve_opts::fast);
    dif=sum(abs(b_new-b_old));
    b_old=b_new;
    if(++ij==maxiters)
      break;
  }
  double p=sum(arma::square(phi)/m)/(n-pcols);
  NumericVector coefs = NumericVector(b_new.begin(), b_new.end());
  coefs.names() = colnames(X);
  
  List l;
  l["coefficients"]=coefs;
  l["phi"]=p;
  l["theta.guesstimate"]=mean(m)/(std::max(p, minphi)-1);
  if(returnfit){
    l["fitted"]=NumericVector(m.begin(), m.end());
  }
  
  return l;
}
