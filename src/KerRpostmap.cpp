/* ---------------------------------------------------
 R posterior summary for estimR in C++
 Copyright, Oswaldo Gressani. All rights reserved.
 ------------------------------------------------------*/

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// [[Rcpp::export]]

List Rcpp_KerRpostmap(NumericMatrix BB, NumericVector theta, NumericMatrix Covar,
                      NumericVector sinter, NumericVector MVvec) {
  // Kernel routine
  // Author : Oswaldo Gressani
  int K = BB.ncol();
  NumericVector muhat = wrap(exp(as<arma::mat>(BB) * as<arma::vec>(theta)));
  int n = muhat.length();
  int simax = sinter.length();
  IntegerVector tt = Rcpp::seq_len(n);

  //------------------ R(t) point estimate
  NumericMatrix R(1, n);
  R(0,0) = muhat[0];
  NumericVector mutemp = muhat[Rcpp::seq(0,(simax - 2))];
  NumericVector sitemp = sinter[Rcpp::seq(0,(simax - 2))];
  for(int t = 1; t <= (simax - 1); t++){
    R(0,t) = muhat[t] *
    pow(Rcpp::sum(Rcpp::rev(mutemp[Rcpp::seq(0, (t - 1))]) *
    sitemp[Rcpp::seq(0, (t - 1))]), (-1));
  }
  for(int t = simax; t <= (n - 1); t++){
    R(0,t) = muhat[t] *
    pow(Rcpp::sum(Rcpp::rev(muhat[seq((t - simax), (t - 1))]) *
    sinter), (-1));
  }
 NumericVector meanlogNorm = Rcpp::log(R);

  //------------------ R(t) logNormal posterior parameters
  NumericMatrix sdlogNorm (1, n);
  NumericVector dhstar_t (K);
  NumericVector BBk;
  NumericVector BBtemp;

  for(int t = 0; t <= (simax - 1); t++){
     for(int k = 0; k <= (K - 1); k++){
       if(t == 0){
         dhstar_t(k) = BB(t,k);
       } else{
        BBk = BB( _ , k );
        BBtemp = Rcpp::rev(BBk[Rcpp::seq(0, (t - 1))]);
        double add = (-1) * (pow(Rcpp::sum(
          Rcpp::rev(mutemp[Rcpp::seq(0, (t - 1))]) *
          sitemp[Rcpp::seq(0, (t - 1))]), (-1))) *
          (Rcpp::sum(Rcpp::rev(mutemp[Rcpp::seq(0, (t - 1))]) *
          sitemp[Rcpp::seq(0, (t - 1))] * BBtemp));
      dhstar_t(k) = BB(t,k) + add;
        }
      }
     NumericVector vectemp (K);
     for(int l = 0; l <= (K - 1); l++){
       vectemp[l] = MVvec[t] * sum(dhstar_t * Covar(_ , l)) * dhstar_t[l];
     }
     sdlogNorm(0,t) = sqrt(sum(vectemp));
  }

  for(int t = simax; t <= (n - 1); t++){
    for(int k = 0; k <= (K - 1); k++){
      BBk = BB( _ , k );
      BBtemp = Rcpp::rev(BBk[Rcpp::seq(0, (t - 1))]);
      double add = (-1) *  pow(Rcpp::sum(
        Rcpp::rev(muhat[seq((t - simax), (t - 1))]) * sinter), (-1)) *
        (Rcpp::sum(Rcpp::rev(muhat[seq((t - simax), (t - 1))]) *
        sinter * BBtemp));
      dhstar_t(k) = BB(t,k) + add;
    }

    NumericVector vectemp (K);
    for(int l = 0; l <= (K - 1); l++){
      vectemp[l] = MVvec[t] * sum(dhstar_t * Covar(_ , l)) * dhstar_t[l];
    }
    sdlogNorm(0,t) = sqrt(sum(vectemp));
  }

  NumericVector Rsd = sqrt((exp(pow(sdlogNorm,2)) - 1) * exp(2 * meanlogNorm +
    pow(sdlogNorm,2)));

  return Rcpp::List::create(Named("R") = R,
                            Named("meanlogNorm") = meanlogNorm,
                            Named("sdlogNorm") = sdlogNorm,
                            Named("Rsd") = Rsd);
}

























