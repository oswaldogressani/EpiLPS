/* ---------------------------------------------------
 Optimization routine (Gradient ascent) in C++
 Copyright, Oswaldo Gressani. All rights reserved.
 ------------------------------------------------------*/

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// [[Rcpp::export]]

List Rcpp_Kerhyperoptim(NumericVector x, NumericMatrix BB, Function grad,
                        double step){

  int iter = 0;
  double tol = 0.01;
  double dist = 3;
  NumericVector xinit = x;
  int K = BB.ncol();
  NumericVector theta0 = rep(1.5,K);

  // Gradient ascent
  while(dist > tol){
    List gradeval = grad(xinit, theta0);
    NumericVector xnew = wrap(as<arma::vec>(xinit) + step * as<arma::vec>(gradeval[0]));
    dist = sqrt(sum(pow((xnew-xinit),2)));
    xinit = xnew;
    theta0 = gradeval[1];
    iter = iter + 1;
  }

  return Rcpp::List::create(Named("iterations") = iter,
                      Named("res") = xinit);
}
