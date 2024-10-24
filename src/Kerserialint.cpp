
/* ---------------------------------------------------
 Serial interval routine in C++
 Copyright, Oswaldo Gressani. All rights reserved.
 ------------------------------------------------------*/

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

List Rcpp_Kerserialint(NumericMatrix x, int B = 5000, 
              NumericVector p = NumericVector::create(0.05,0.25,0.50,0.75,0.95)) {
  
  double n = x.nrow();        // Sample size
  double ninv = 1 / n;        // Inverse of sample size
  int plen = p.length();      // Length of probability vector for quantiles
  NumericVector sL = x(_,0);  // Left boundary of serial interval window
  NumericVector sR = x(_,1);  // Right boundary of serial interval window
  NumericVector ds = sR - sL; 
  
  // Create the data matrix
  // NumericMatrix data(n,2);
  // data(_,0) = sL;
  // data(_,1) = sR;
  // colnames(data) = CharacterVector::create("sL","sR");
  
  // Create the order statistics from serial interval boundary values
  NumericVector sorder(2 * n);
  for(int i = 0; i < n; i++){
    sorder[i] = sL[i];
    sorder[i + n] = sR[i];
  }
  sorder = sorder.sort(); 
  
  // Evaluate cdf on order statistics
  NumericVector Fsorder(2 * n);
  for(int j = 0; j < (2*n); j++){
    NumericVector Fsorder_temp(n);
    for(int i = 0; i < n; i++){
      Fsorder_temp[i] = ((sorder[j] - sL[i]) / ds[i]) * 
        (sorder[j] >= sL[i] and sorder[j] < sR[i]) +  (sorder[j] >= sR[i]);
    }
    Fsorder[j] = ninv * sum(Fsorder_temp);
  } 
  NumericVector diffF = diff(Fsorder);
  
// Bootstrap
 int nfeats = 2 + plen;
 NumericMatrix featboot(B, nfeats);

 for(int b = 0; b < B; b++){
 NumericVector uvec = runif(n);
 NumericVector bootsample(n);
 for(int i = 0; i < n; i++){
   double u = uvec[i];
   int idxl = sum((u-Fsorder)>0) - 1;
   int idxr = idxl + 1;
   double deltaF = diffF[idxl];
   if(deltaF > 0){
    bootsample[i] = ((Fsorder[idxr]-u)*sorder[idxl] + 
      (u-Fsorder[idxl])*sorder[idxr])/deltaF;
   } else{
    bootsample[i] = sorder[idxl];
  }
 }
 
 // Fill feature matrix based on bootstrap sample
 featboot(b,0) = mean(bootsample);
 featboot(b,1) = sd(bootsample);
 
 // Compute quantiles following the default type 7 quantile rule in R
 NumericVector bootsort = bootsample.sort(); 
 for(int i = 0; i < plen; i++){
   double prob = p[i];
   double m = 1-prob;
   double j = floor(n * prob + m);
   double gamma = n * prob + m - j;
   featboot(b, i + 2) = (1 - gamma) * bootsort[j - 1] + gamma * bootsort[j];
  }
 }
 
 // Compute summary statistics
 NumericMatrix estimres(nfeats, 5);
 colnames(estimres) = CharacterVector::create("Estim", "CI90p_l", "CI90p_r",
          "CI95p_l", "CI95p_r");
 rownames(estimres) = CharacterVector::create("mean", "sd", "q05", "q25",
          "q50", "q75", "q95");
 NumericVector probCI = NumericVector::create(0.05,0.95,0.025,0.975);
 for(int i = 0; i < nfeats; i++){
   // Results for the mean of the serial interval
   NumericVector featvec = featboot(_,i);
   estimres(i, 0) = mean(featvec);
   NumericVector featvecsort = featvec.sort();
   for(int l = 0; l < 4; l++){
     double prob = probCI[l];
     double m = 1-prob;
     double j = floor(B * prob + m);
     double gamma = B * prob + m - j;
     estimres(i, l + 1)= (1 - gamma) * featvecsort[j - 1] + 
       gamma * featvecsort[j];
   }
 }

 return Rcpp::List::create(Named("estimres") = estimres);
  
  
}

