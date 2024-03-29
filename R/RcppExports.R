# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

Rcpp_KerLaplace <- function(theta0, rho, lambda, K, Dlogptheta, D2logptheta) {
    .Call(`_EpiLPS_Rcpp_KerLaplace`, theta0, rho, lambda, K, Dlogptheta, D2logptheta)
}

Rcpp_KerLaplaceIncub <- function(theta0, lambda, K, Dlogptheta, D2logptheta) {
    .Call(`_EpiLPS_Rcpp_KerLaplaceIncub`, theta0, lambda, K, Dlogptheta, D2logptheta)
}

Rcpp_KerLaplaceNowcast <- function(xi0, v, dimxi, Dlogpxi, D2logpxi) {
    .Call(`_EpiLPS_Rcpp_KerLaplaceNowcast`, xi0, v, dimxi, Dlogpxi, D2logpxi)
}

Rcpp_KerMVN <- function(mu, Sigma) {
    .Call(`_EpiLPS_Rcpp_KerMVN`, mu, Sigma)
}

Rcpp_KerRpostmap <- function(BB, theta, Covar, sinter, MVvec) {
    .Call(`_EpiLPS_Rcpp_KerRpostmap`, BB, theta, Covar, sinter, MVvec)
}

Rcpp_KerRpostmcmc <- function(t, BB, sinter, thetasample) {
    .Call(`_EpiLPS_Rcpp_KerRpostmcmc`, t, BB, sinter, thetasample)
}

Rcpp_KercubicBspline <- function(x, lower, upper, K) {
    .Call(`_EpiLPS_Rcpp_KercubicBspline`, x, lower, upper, K)
}

Rcpp_Kerhyperoptim <- function(x, BB, grad, step) {
    .Call(`_EpiLPS_Rcpp_Kerhyperoptim`, x, BB, grad, step)
}

