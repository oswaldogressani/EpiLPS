/* ---------------------------------------------------
Laplace approximation in C++ for nowcasting routine
Copyright, Oswaldo Gressani. All rights reserved.
 ------------------------------------------------------*/

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// [[Rcpp::export]]

List Rcpp_KerLaplaceNowcast(NumericVector xi0, NumericVector v, int dimxi,
                  Function Dlogpxi, Function D2logpxi){
  // Kernel routine
  // Author : Oswaldo Gressani
  int iter = 0;
  double tol = 1e-3;
  double dist = 3;
  arma::vec xinew(dimxi);
  NumericVector gg;
  NumericMatrix H;
  arma::mat invH(dimxi,dimxi);

    //Newton-Raphson
    while(dist > tol){
      gg = Dlogpxi(xi0, v);
      H = D2logpxi(xi0, v);
      invH = arma::inv(as<arma::mat>(H));
      xinew = as<arma::vec>(xi0) - invH * (as<arma::vec>(gg));
      NumericVector xinew2 = wrap(xinew);
      dist = sqrt(sum(pow((xi0-xinew2),2)));
      xi0 = xinew2;
      iter = iter + 1;
    }

  return Rcpp::List::create(Named("xistar") = xi0);
}
