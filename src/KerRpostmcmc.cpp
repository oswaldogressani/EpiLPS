/* ---------------------------------------------------
 Computing  summary R values for MCMC output in C++
 Copyright, Oswaldo Gressani. All rights reserved.
 ------------------------------------------------------*/

#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]

NumericVector Rcpp_KerRpostmcmc(int t, NumericMatrix BB, NumericVector sinter, 
                       NumericMatrix thetasample) {
  // Kernel routine
  // Author : Oswaldo Gressani
  int K = BB.ncol();
  int n = BB.nrow();
  int M = thetasample.nrow();
  int simax = sinter.length();
  
  // Computing R sample at time point t
  NumericMatrix Rsample (M, 1);
  NumericVector musample (n);
  
  for(int m = 0; m < M; m++){
    NumericVector mutemp (K);
    for(int j = 0; j < n; j++){
      for(int k = 0; k < K; k++){
        mutemp(k) = BB(j, k) * thetasample(m , k);
    }
      musample(j) = exp(sum(mutemp));
    }
    
  // compute Rsample
  if(t == 1){
    Rsample(m, 0) = musample(0);
  } else if (t >= 2 && t <= simax){
      NumericVector mutemp = musample[Rcpp::seq(0,(simax - 2))];
      NumericVector sitemp = sinter[Rcpp::seq(0,(simax - 2))];
      Rsample(m, 0) = musample(t - 1) *
      pow(Rcpp::sum(Rcpp::rev(mutemp[Rcpp::seq(0, (t - 2))]) *
      sitemp[Rcpp::seq(0, (t - 2))]), (-1));
  } else if (t > simax && t <= n){
      Rsample(m,0) = musample(t - 1) *
      pow(Rcpp::sum(Rcpp::rev(musample[seq((t - simax - 1), (t - 2))]) *
      sinter), (-1)); 
  }
  

  }

  return(Rsample);
}