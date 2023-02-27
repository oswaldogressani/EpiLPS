/* ---------------------------------------------------
Cubic B-spline matrix in C++
Copyright, Oswaldo Gressani. All rights reserved.
------------------------------------------------------*/

#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]

NumericMatrix Rcpp_KercubicBspline(NumericVector x, double lower,
                                   double upper, int K){
  // Kernel routine
  // Author : Oswaldo Gressani
  int nx = x.length();
  int ndx = K-3;
  double dx = (upper - lower) / ndx;
  int nknots = ndx + 2 * 3 + 1;

  NumericVector knots(nknots);
    knots[0] = lower-3*dx;
      for(int i=1; i<nknots; i++){
        knots[i] = knots[i-1] + dx;
      }

  NumericMatrix B(nx, K);
  double temp;
  double cub;

      // Fill matrix B
      for(int i=0; i<nx; i++){
        for(int j=0; j<(nknots-4); j++){
          temp = 0;
          cub = x[i] - knots[j];
          if (cub > 0) {
            temp += pow(cub, 3);
            cub = x[i]-knots[j+1];
            if(cub > 0){
              temp -= 4*pow(cub,3);
              cub=x[i]-knots[j+2];
              if(cub>0){
                temp += 6*pow(cub,3);
                cub = x[i]-knots[j+3];
                if(cub>0){
                  temp-=4*pow(cub,3);
                  cub=x[i]-knots[j+4];
                  if(cub>0){
                    temp+=pow(cub,3);
                  }
                }
              }
            }
          }
          B(i,j) = temp / (6*pow(dx,3));
            if(std::abs(B(i,j))<1e-10){
              B(i,j)=0;
              }
        }
      }

  return(B);
}


