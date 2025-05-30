// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double eucNorm(arma::vec x, int d){
  double res = 0;
  for(int i = 0; i < d; i++) res = res + pow(x(i),2);
  return(sqrt(res));  
}

// [[Rcpp::export]]
double MAD(arma::vec x, double med, int n){
  for(int i = 0; i < n; i++) x(i) = abs(x(i) - med);
  return(arma::median(x));  
}

// [[Rcpp::export]]
arma::mat RPD_outl_C(arma::mat x, arma::mat X, arma::vec beta, bool quant, int n, int d, int m, int ndir, int ndir2) {
  arma::vec res(m, arma::fill::zeros);
  arma::vec Xu(n, arma::fill::zeros), xu(m, arma::fill::zeros);
  double Xm, XMAD, temp;
  bool cont = true;
  int j = 0, jraw = 0, Jmax = pow(10,6);
  
  if(beta(0)<=arma::datum::eps){
    quant = false;
    beta(0) = 0;
    }
  if(quant){
    if(beta(0) > 1 - arma::datum::eps) Rcpp::stop("The quantile in beta must be between 0 and 1.");
    // only if beta is given as the quantile value
    // find the value of cutoff for regularization based on beta
    arma::vec MADs(ndir, arma::fill::zeros);
    for(int j = 0; j < ndir; j++){
      arma::vec u(d, arma::fill::randn);
      u = u/eucNorm(u, d);
      for(int i = 0; i < n; i++){
        Xu(i) = dot(u, X.row(i));
      }
      Xm = arma::median(Xu);
      MADs(j) = MAD(Xu, Xm, n);  
    }
    beta = arma::quantile(MADs, beta);
    }
  
  j = 0;
  while(cont){
    jraw = jraw + 1;
    arma::vec u(d, arma::fill::randn);
    u = u/eucNorm(u, d);
    for(int i = 0; i < n; i++){
      Xu(i) = dot(u, X.row(i));
    }
    Xm = arma::median(Xu);
    XMAD = MAD(Xu, Xm, n);
    if(XMAD > beta(0)){
      j = j + 1;
      for(int i = 0; i < m; i++){
        xu(i) = dot(u, x.row(i));
        temp = abs(xu(i) - Xm)/XMAD;
        if(res(i)<temp) res(i) = temp;
      }
    }
    if(j >= ndir) cont = false;
    if(jraw >= Jmax){
      cont = false;
      Rcpp::stop("Exceeded maximum number of loops.");
    }
  }
  return(res);  
}