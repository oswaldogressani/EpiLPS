/* ---------------------------------------------------
 Laplace approximation in C++
 Copyright, Oswaldo Gressani. All rights reserved.
 ------------------------------------------------------*/

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// [[Rcpp::export]]

List Rcpp_Laplace(double rho, double lambda, int K,
                  Function Dlogptheta, Function D2logptheta){

  int iter = 0;
  double tol = 1e-5;
  double dist = 3;
  NumericVector theta0 = rep(1.5,K);
  arma::vec thetanew(K);
  NumericVector gg;
  NumericMatrix H;
  arma::mat invH(K,K);


    //Newton-Raphson
    while(dist > tol){
      gg = Dlogptheta(theta0, rho, lambda);
      H = D2logptheta(theta0, rho, lambda);
      invH = arma::inv(as<arma::mat>(H));
      thetanew = as<arma::vec>(theta0) - invH * (as<arma::vec>(gg));
      NumericVector thetanew2 = wrap(thetanew);
      dist = sqrt(sum(pow((theta0-thetanew2),2)));
      theta0 = thetanew2;
      iter = iter + 1;
    }

  arma::mat Sigmastar = (-1) * arma::inv(as<arma::mat>(D2logptheta(theta0,rho,lambda)));
  arma:: cx_double logdetSigmastar = arma::log_det(Sigmastar);

  return Rcpp::List::create(Named("Lapmode") = theta0,
                      Named("Lapvar") = Sigmastar,
                      Named("logdetSigma") = logdetSigmastar);
}


