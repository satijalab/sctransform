// [[Rcpp::depends(RcppArmadillo)]]

#include "RcppArmadillo.h"
#include "math.h"

using namespace Rcpp;

// from Rcpp gallery https://gallery.rcpp.org/articles/stl-random-shuffle/
// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
// inline int randWrapper(const int n) { return floor(unif_rand()*n); }
std::random_device rd;
std::mt19937 randWrapper(rd());


// [[Rcpp::export]]
NumericVector row_mean_dgcmatrix(S4 matrix) {
  NumericVector x = matrix.slot("x");
  IntegerVector i = matrix.slot("i");
  IntegerVector dim = matrix.slot("Dim");
  int rows = dim[0];
  int cols = dim[1];

  NumericVector ret(rows, 0.0);
  int x_length = x.length();
  for (int k=0; k<x_length; ++k) {
    ret[i[k]] += x[k];
  }
  for (int k=0; k<rows; ++k) {
    ret[k] /= cols;
  }
  List dn = matrix.slot("Dimnames");
  if (dn[0] != R_NilValue) {
    ret.attr("names") = as<CharacterVector>(dn[0]);
  }
  return ret;
}

// [[Rcpp::export]]
NumericMatrix row_mean_grouped_dgcmatrix(S4 matrix, IntegerVector group,
                                         bool shuffle) {
  NumericVector x = matrix.slot("x");
  IntegerVector i = matrix.slot("i");
  IntegerVector p = matrix.slot("p");
  IntegerVector dim = matrix.slot("Dim");
  int rows = dim[0];
  int cols = dim[1];
  CharacterVector levs = group.attr("levels");
  int groups = levs.length();
  NumericMatrix ret(rows, groups);
  IntegerVector groupsize(groups, 0);
  int x_length = x.length();

  if (shuffle) {
    group = clone(group);
    std::shuffle(group.begin(), group.end(), randWrapper);
  }

  int col = 0;
  for (int k=0; k<x_length; ++k) {
    while (k>=p[col]) {
      ++col;
      ++groupsize[group[col-1]-1];
    }
    ret(i[k], group[col-1]-1) += x[k];
  }
  while (col < cols) {
    ++col;
    ++groupsize[group[col-1]-1];
  }

  for (int j=0; j<groups; ++j) {
    if (groupsize[j] == 0) {
      ret(_, j) = rep(NumericVector::get_na(), rows);
    } else{
      ret(_, j) = ret(_, j) / groupsize[j];
    }
  }
  colnames(ret) = levs;
  List dn = matrix.slot("Dimnames");
  if (dn[0] != R_NilValue) {
    rownames(ret) = as<CharacterVector>(dn[0]);
  }
  return ret;
}

// [[Rcpp::export]]
NumericVector row_gmean_dgcmatrix(S4 matrix, double eps) {
  NumericVector x = matrix.slot("x");
  IntegerVector i = matrix.slot("i");
  IntegerVector dim = matrix.slot("Dim");
  int rows = dim[0];
  int cols = dim[1];

  NumericVector ret(rows, 0.0);
  IntegerVector nzero(rows, cols);
  int x_length = x.length();
  double log_eps = log(eps);

  for (int k=0; k<x_length; ++k) {
    ret[i[k]] += log(x[k] + eps);
    nzero[i[k]] -= 1;
  }
  for (int k=0; k<rows; ++k) {
    ret[k] = exp((ret[k] + log_eps * nzero[k]) / cols) - eps;
  }
  List dn = matrix.slot("Dimnames");
  if (dn[0] != R_NilValue) {
    ret.attr("names") = as<CharacterVector>(dn[0]);
  }
  return ret;
}

// [[Rcpp::export]]
NumericMatrix row_gmean_grouped_dgcmatrix(S4 matrix, IntegerVector group,
                                          double eps, bool shuffle) {
  NumericVector x = matrix.slot("x");
  IntegerVector i = matrix.slot("i");
  IntegerVector p = matrix.slot("p");
  IntegerVector dim = matrix.slot("Dim");
  int rows = dim[0];
  int cols = dim[1];
  CharacterVector levs = group.attr("levels");
  int groups = levs.length();
  NumericMatrix ret(rows, groups);
  IntegerVector groupsize(groups, 0);
  int x_length = x.length();
  IntegerMatrix nonzero(rows, groups);
  double log_eps = log(eps);

  if (shuffle) {
    group = clone(group);
    std::shuffle(group.begin(), group.end(), randWrapper);
  }

  int col = 0;
  for (int k=0; k<x_length; ++k) {
    while (k>=p[col]) {
      ++col;
      ++groupsize[group[col-1]-1];
    }
    ret(i[k], group[col-1]-1) += log(x[k] + eps);
    ++nonzero(i[k], group[col-1]-1);
  }
  while (col < cols) {
    ++col;
    ++groupsize[group[col-1]-1];
  }

  for (int j=0; j<groups; ++j) {
    for (int k=0; k<rows; ++k) {
      ret(k, j) = exp((ret(k, j) + log_eps * (groupsize[j] - nonzero(k, j))) / groupsize[j]) - eps;
    }
  }
  colnames(ret) = levs;
  List dn = matrix.slot("Dimnames");
  if (dn[0] != R_NilValue) {
    rownames(ret) = as<CharacterVector>(dn[0]);
  }
  return ret;
}

// [[Rcpp::export]]
IntegerVector row_nonzero_count_dgcmatrix(S4 matrix) {
  IntegerVector i = matrix.slot("i");
  IntegerVector dim = matrix.slot("Dim");
  int rows = dim[0];

  IntegerVector ret(rows, 0);
  int i_len = i.length();
  for(int k = 0; k < i_len; ++k) {
    ret[i[k]]++;
  }

  List dn = matrix.slot("Dimnames");
  if (dn[0] != R_NilValue) {
    ret.attr("names") = as<CharacterVector>(dn[0]);
  }
  return ret;
}

// [[Rcpp::export]]
IntegerMatrix row_nonzero_count_grouped_dgcmatrix(S4 matrix, IntegerVector group) {
  IntegerVector p = matrix.slot("p");
  IntegerVector i = matrix.slot("i");
  int i_length = i.length();
  IntegerVector dim = matrix.slot("Dim");
  int rows = dim[0];
  CharacterVector levs = group.attr("levels");
  int groups = levs.length();
  IntegerMatrix ret(rows, groups);

  int col = 0;
  for (int k=0; k<i_length; ++k) {
    while (k>=p[col]) {
      ++col;
    }
    ret(i[k], group[col-1]-1)++;
  }

  colnames(ret) = levs;
  List dn = matrix.slot("Dimnames");
  if (dn[0] != R_NilValue) {
    rownames(ret) = as<CharacterVector>(dn[0]);
  }
  return ret;
}

// [[Rcpp::export]]
NumericVector row_var_dgcmatrix(NumericVector x, IntegerVector i, int rows, int cols) {
  NumericVector rowmean(rows, 0.0);
  int x_length = x.length();
  for (int k=0; k<x_length; ++k) {
    rowmean[i[k]] += x[k];
  }
  for (int k=0; k<rows; ++k) {
    rowmean[k] /= cols;
  }
  NumericVector rowvar(rows, 0.0);
  IntegerVector nzero(rows, cols);
  for (int k=0; k<x_length; ++k) {
    rowvar[i[k]] += pow(x[k] - rowmean[i[k]], 2);
    nzero[i[k]] -= 1;
  }
  for (int k=0; k<rows; ++k) {
    rowvar[k] = (rowvar[k] + (pow(rowmean[k], 2) * nzero[k])) / (cols - 1);
  }
  return rowvar;
}

// assume two groups
// assume group integers to be 0 and 1
// [[Rcpp::export]]
NumericVector grouped_mean_diff_per_row(NumericMatrix x, IntegerVector group, bool shuffle) {
  int nrows = x.nrow();
  int ncols = x.ncol();
  NumericMatrix tmp(2, nrows);
  IntegerVector groupsize(2);
  NumericVector ret(nrows, 0.0);

  if (shuffle) {
    group = clone(group);
    std::shuffle(group.begin(), group.end(), randWrapper);
  }

  for (int i = 0; i < ncols; i++) {
    ++groupsize(group(i));
    for (int j = 0; j < nrows; j++) {
      tmp(group(i), j) += x(j,i);
    }
  }
  for (int j = 0; j < nrows; j++) {
    ret(j) = (tmp(0, j) / groupsize(0)) - (tmp(1, j) / groupsize(1));
  }
  return ret;
}

// Bootstrapped mean
// [[Rcpp::export]]
NumericVector mean_boot(NumericVector x, int N, int S) {
  NumericVector ret(N);
  for (int i = 0; i < N; i++) {
    ret(i) = mean(sample(x, S, true));
  }
  return ret;
}

// Bootstrapped mean per group
// assume group vector uses contiguous integers starting from 0
// [[Rcpp::export]]
NumericMatrix mean_boot_grouped(NumericVector x, IntegerVector group, int N, int S) {
  // how many groups are there
  int groups = max(group) + 1;
  // we need as many columns
  NumericMatrix ret(N, groups);

  for (int g = 0; g < groups; g++) {
    NumericVector xg = x[group == g];
    ret(_, g) = mean_boot(xg, N, S);
  }
  return ret;
}

// Based on bootstrapped means per group quantify differences
// x is the result of mean_boot_grouped; make sure there are only two columns
// [[Rcpp::export]]
NumericVector distribution_shift(NumericMatrix x) {
  int N = x.nrow();
  int cs = 0;
  int cs_sum = 0;
  double sd0, sd1;
  // we also want to return three quantile scores per group
  IntegerVector qidx = {
    (int) floor(0.16 * N - 1),
    (int) round(0.5 * N - 1),
    (int) ceil(0.84 * N - 1)
  };
  // the quantile results and current index
  NumericVector res(8);
  int q0i = 0;
  int q1i = 0;
  // current rank per group
  int r0 = 0;
  int r1 = 0;
  arma::uvec indices = arma::sort_index(as<arma::mat>(x));
  arma::uvec::const_iterator it = indices.begin();
  arma::uvec::const_iterator it_end = indices.end();
  for (; it != it_end; ++it) {
    if ((*it) < N) {
      cs += 1;
      if ((q0i < 3) & (r0 == qidx[q0i])) {
        res[q0i] = x[*it];
        q0i++;
      }
      r0++;
    } else {
      cs -= 1;
      if ((q1i < 3) & (r1 == qidx[q1i])) {
        res[q1i+3] = x[*it];
        q1i++;
      }
      r1++;
    }
    cs_sum += cs;
  }
  res[6] = (double) cs_sum / N / N;
  // add z-score-like score
  if (res[4] > res[1]) {
    // second group has higher mean
    sd0 = res[2] - res[1];
    sd1 = res[4] - res[3];
  } else {
    sd0 = res[1] - res[0];
    sd1 = res[5] - res[4];
  }
  //res[7] = (res[4] - res[1]) / sqrt(sd0 * sd1);
  //res[7] = (res[4] - res[1]) / ((sd0 + sd1) / 2);
  res[7] = (res[4] - res[1]) / sqrt((sd0*sd0 + sd1*sd1) / 2);
  return res;
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
    if((unique_vals.n_elem==2) && ((unique_vals[0] == 0) || (unique_vals[1] == 0))){
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
