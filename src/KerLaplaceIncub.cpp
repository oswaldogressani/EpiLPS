/* ---------------------------------------------------
Laplace approximation in C++
Copyright, Oswaldo Gressani. All rights reserved.
 ------------------------------------------------------*/

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// [[Rcpp::export]]

List Rcpp_KerLaplaceIncub(NumericVector theta0, double lambda, int K,
                  Function Dlogptheta, Function D2logptheta){
  // Kernel routine
  // Author : Oswaldo Gressani
  int iter = 0;
  double tol = 1e-5;
  double dist = 3;
  arma::vec thetanew(K);
  NumericVector gg;
  NumericMatrix H;
  arma::mat invH(K,K);

    //Newton-Raphson
    while(dist > tol){
      gg = Dlogptheta(theta0, lambda);
      H = D2logptheta(theta0, lambda);
      invH = arma::inv(as<arma::mat>(H));
      thetanew = as<arma::vec>(theta0) - invH * (as<arma::vec>(gg));
      NumericVector thetanew2 = wrap(thetanew);
      dist = sqrt(sum(pow((theta0-thetanew2),2)));
      theta0 = thetanew2;
      iter = iter + 1;
      if(iter > 50) break;
    }

  arma::mat Hess = as<arma::mat>(D2logptheta(theta0,lambda));
  arma::mat Sigmastar = (-1) * arma::inv(Hess);
  arma::mat Sigmainv = (-1) * Hess;
  arma::cx_double logdetSigmastar = arma::log_det(Sigmastar);
  arma::cx_double detSigmastar = exp(logdetSigmastar);
  arma::cx_double invdetSigmastar = pow(detSigmastar,(-1));

 // double epsil = 1e-5;
 // double lambdainc = lambda + epsil;
 // arma::mat Sigmainc = (-1) * arma::inv(as<arma::mat>(D2logptheta(theta0,lambdainc)));
 // arma::mat diff = (1/epsil) * (Sigmainc-Sigmastar);
 // arma::cx_double dSigma =  invdetSigmastar * (arma::trace(detSigmastar * Sigmainv * diff));

  return Rcpp::List::create(Named("Lapmode") = theta0,
                            Named("Lapvar") = Sigmastar,
                            Named("logdetSigma") = logdetSigmastar);
}
