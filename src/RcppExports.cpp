// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Rcpp_KerLaplace
List Rcpp_KerLaplace(NumericVector theta0, double rho, double lambda, int K, Function Dlogptheta, Function D2logptheta);
RcppExport SEXP _EpiLPS_Rcpp_KerLaplace(SEXP theta0SEXP, SEXP rhoSEXP, SEXP lambdaSEXP, SEXP KSEXP, SEXP DlogpthetaSEXP, SEXP D2logpthetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta0(theta0SEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< Function >::type Dlogptheta(DlogpthetaSEXP);
    Rcpp::traits::input_parameter< Function >::type D2logptheta(D2logpthetaSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_KerLaplace(theta0, rho, lambda, K, Dlogptheta, D2logptheta));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_KerLaplaceIncub
List Rcpp_KerLaplaceIncub(NumericVector theta0, double lambda, int K, Function Dlogptheta, Function D2logptheta);
RcppExport SEXP _EpiLPS_Rcpp_KerLaplaceIncub(SEXP theta0SEXP, SEXP lambdaSEXP, SEXP KSEXP, SEXP DlogpthetaSEXP, SEXP D2logpthetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta0(theta0SEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< Function >::type Dlogptheta(DlogpthetaSEXP);
    Rcpp::traits::input_parameter< Function >::type D2logptheta(D2logpthetaSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_KerLaplaceIncub(theta0, lambda, K, Dlogptheta, D2logptheta));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_KerLaplaceNowcast
List Rcpp_KerLaplaceNowcast(NumericVector xi0, NumericVector v, int dimxi, Function Dlogpxi, Function D2logpxi);
RcppExport SEXP _EpiLPS_Rcpp_KerLaplaceNowcast(SEXP xi0SEXP, SEXP vSEXP, SEXP dimxiSEXP, SEXP DlogpxiSEXP, SEXP D2logpxiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type xi0(xi0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type dimxi(dimxiSEXP);
    Rcpp::traits::input_parameter< Function >::type Dlogpxi(DlogpxiSEXP);
    Rcpp::traits::input_parameter< Function >::type D2logpxi(D2logpxiSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_KerLaplaceNowcast(xi0, v, dimxi, Dlogpxi, D2logpxi));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_KerMVN
NumericVector Rcpp_KerMVN(arma::vec mu, arma::mat Sigma);
RcppExport SEXP _EpiLPS_Rcpp_KerMVN(SEXP muSEXP, SEXP SigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_KerMVN(mu, Sigma));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_KerRpostmap
List Rcpp_KerRpostmap(NumericMatrix BB, NumericVector theta, NumericMatrix Covar, NumericVector sinter, NumericVector MVvec);
RcppExport SEXP _EpiLPS_Rcpp_KerRpostmap(SEXP BBSEXP, SEXP thetaSEXP, SEXP CovarSEXP, SEXP sinterSEXP, SEXP MVvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type BB(BBSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Covar(CovarSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sinter(sinterSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type MVvec(MVvecSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_KerRpostmap(BB, theta, Covar, sinter, MVvec));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_KerRpostmcmc
NumericVector Rcpp_KerRpostmcmc(int t, NumericMatrix BB, NumericVector sinter, NumericMatrix thetasample);
RcppExport SEXP _EpiLPS_Rcpp_KerRpostmcmc(SEXP tSEXP, SEXP BBSEXP, SEXP sinterSEXP, SEXP thetasampleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type BB(BBSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sinter(sinterSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type thetasample(thetasampleSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_KerRpostmcmc(t, BB, sinter, thetasample));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_KercubicBspline
NumericMatrix Rcpp_KercubicBspline(NumericVector x, double lower, double upper, int K);
RcppExport SEXP _EpiLPS_Rcpp_KercubicBspline(SEXP xSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< double >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_KercubicBspline(x, lower, upper, K));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_Kerhyperoptim
List Rcpp_Kerhyperoptim(NumericVector x, NumericMatrix BB, Function grad, double step);
RcppExport SEXP _EpiLPS_Rcpp_Kerhyperoptim(SEXP xSEXP, SEXP BBSEXP, SEXP gradSEXP, SEXP stepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type BB(BBSEXP);
    Rcpp::traits::input_parameter< Function >::type grad(gradSEXP);
    Rcpp::traits::input_parameter< double >::type step(stepSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_Kerhyperoptim(x, BB, grad, step));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_Kerserialint
List Rcpp_Kerserialint(NumericMatrix x, int B, NumericVector p);
RcppExport SEXP _EpiLPS_Rcpp_Kerserialint(SEXP xSEXP, SEXP BSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_Kerserialint(x, B, p));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_EpiLPS_Rcpp_KerLaplace", (DL_FUNC) &_EpiLPS_Rcpp_KerLaplace, 6},
    {"_EpiLPS_Rcpp_KerLaplaceIncub", (DL_FUNC) &_EpiLPS_Rcpp_KerLaplaceIncub, 5},
    {"_EpiLPS_Rcpp_KerLaplaceNowcast", (DL_FUNC) &_EpiLPS_Rcpp_KerLaplaceNowcast, 5},
    {"_EpiLPS_Rcpp_KerMVN", (DL_FUNC) &_EpiLPS_Rcpp_KerMVN, 2},
    {"_EpiLPS_Rcpp_KerRpostmap", (DL_FUNC) &_EpiLPS_Rcpp_KerRpostmap, 5},
    {"_EpiLPS_Rcpp_KerRpostmcmc", (DL_FUNC) &_EpiLPS_Rcpp_KerRpostmcmc, 4},
    {"_EpiLPS_Rcpp_KercubicBspline", (DL_FUNC) &_EpiLPS_Rcpp_KercubicBspline, 4},
    {"_EpiLPS_Rcpp_Kerhyperoptim", (DL_FUNC) &_EpiLPS_Rcpp_Kerhyperoptim, 4},
    {"_EpiLPS_Rcpp_Kerserialint", (DL_FUNC) &_EpiLPS_Rcpp_Kerserialint, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_EpiLPS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
