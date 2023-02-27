#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// [[Rcpp::export]]

NumericVector Rcpp_KerMVN(arma::vec mu, arma::mat Sigma) {
  // Kernel routine
  // Author : Oswaldo Gressani
  // Draw a vector from multivariate normal distribution
  arma::vec out = arma::mvnrnd(mu, Sigma);
  return wrap(out.t());
}




