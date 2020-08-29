// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// Generalized_Pareto_calc_function
double Generalized_Pareto_calc_function(int m, double c, double b, double rho);
RcppExport SEXP _SkeweDF_Generalized_Pareto_calc_function(SEXP mSEXP, SEXP cSEXP, SEXP bSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(Generalized_Pareto_calc_function(m, c, b, rho));
    return rcpp_result_gen;
END_RCPP
}
// Generalized_Pareto_calc_pi
double Generalized_Pareto_calc_pi(double k, double theta, double c, double b, double rho);
RcppExport SEXP _SkeweDF_Generalized_Pareto_calc_pi(SEXP kSEXP, SEXP thetaSEXP, SEXP cSEXP, SEXP bSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(Generalized_Pareto_calc_pi(k, theta, c, b, rho));
    return rcpp_result_gen;
END_RCPP
}
// Generalized_Pareto_calc_P0_delta
double Generalized_Pareto_calc_P0_delta(double theta, double c, double b, double rho, double delta);
RcppExport SEXP _SkeweDF_Generalized_Pareto_calc_P0_delta(SEXP thetaSEXP, SEXP cSEXP, SEXP bSEXP, SEXP rhoSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    rcpp_result_gen = Rcpp::wrap(Generalized_Pareto_calc_P0_delta(theta, c, b, rho, delta));
    return rcpp_result_gen;
END_RCPP
}
// Generalized_Pareto_calc_P0_iter
double Generalized_Pareto_calc_P0_iter(double theta, double c, double b, double rho, int iter);
RcppExport SEXP _SkeweDF_Generalized_Pareto_calc_P0_iter(SEXP thetaSEXP, SEXP cSEXP, SEXP bSEXP, SEXP rhoSEXP, SEXP iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    rcpp_result_gen = Rcpp::wrap(Generalized_Pareto_calc_P0_iter(theta, c, b, rho, iter));
    return rcpp_result_gen;
END_RCPP
}
// Generalized_Pareto
std::vector<double> Generalized_Pareto(int k, double theta, double c, double b, double rho);
RcppExport SEXP _SkeweDF_Generalized_Pareto(SEXP kSEXP, SEXP thetaSEXP, SEXP cSEXP, SEXP bSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(Generalized_Pareto(k, theta, c, b, rho));
    return rcpp_result_gen;
END_RCPP
}
// Kolmogorov_Waring_P0_calc
double Kolmogorov_Waring_P0_calc(double a, double b, double t);
RcppExport SEXP _SkeweDF_Kolmogorov_Waring_P0_calc(SEXP aSEXP, SEXP bSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(Kolmogorov_Waring_P0_calc(a, b, t));
    return rcpp_result_gen;
END_RCPP
}
// Kolmogorov_Waring_P0
double Kolmogorov_Waring_P0(double a, double b, double theta);
RcppExport SEXP _SkeweDF_Kolmogorov_Waring_P0(SEXP aSEXP, SEXP bSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(Kolmogorov_Waring_P0(a, b, theta));
    return rcpp_result_gen;
END_RCPP
}
// Kolmogorov_Waring
std::vector<double> Kolmogorov_Waring(int n, double a, double b, double theta);
RcppExport SEXP _SkeweDF_Kolmogorov_Waring(SEXP nSEXP, SEXP aSEXP, SEXP bSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(Kolmogorov_Waring(n, a, b, theta));
    return rcpp_result_gen;
END_RCPP
}
// RGHD_P0_calc_function
double RGHD_P0_calc_function(int k, int m, Rcpp::NumericVector r, Rcpp::NumericVector q);
RcppExport SEXP _SkeweDF_RGHD_P0_calc_function(SEXP kSEXP, SEXP mSEXP, SEXP rSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(RGHD_P0_calc_function(k, m, r, q));
    return rcpp_result_gen;
END_RCPP
}
// RGHD_P0_calc_pi
double RGHD_P0_calc_pi(int y, int m, Rcpp::NumericVector r, Rcpp::NumericVector q);
RcppExport SEXP _SkeweDF_RGHD_P0_calc_pi(SEXP ySEXP, SEXP mSEXP, SEXP rSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(RGHD_P0_calc_pi(y, m, r, q));
    return rcpp_result_gen;
END_RCPP
}
// RGHD_P0_calc
double RGHD_P0_calc(int sigma_upper, int m, Rcpp::NumericVector r, Rcpp::NumericVector q);
RcppExport SEXP _SkeweDF_RGHD_P0_calc(SEXP sigma_upperSEXP, SEXP mSEXP, SEXP rSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type sigma_upper(sigma_upperSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(RGHD_P0_calc(sigma_upper, m, r, q));
    return rcpp_result_gen;
END_RCPP
}
// RGHD_P0
double RGHD_P0(int m, Rcpp::NumericVector r, Rcpp::NumericVector q);
RcppExport SEXP _SkeweDF_RGHD_P0(SEXP mSEXP, SEXP rSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(RGHD_P0(m, r, q));
    return rcpp_result_gen;
END_RCPP
}
// RGHD
std::vector<double> RGHD(int J, int m, Rcpp::NumericVector r, Rcpp::NumericVector q, int P0_iter, bool P0_included);
RcppExport SEXP _SkeweDF_RGHD(SEXP JSEXP, SEXP mSEXP, SEXP rSEXP, SEXP qSEXP, SEXP P0_iterSEXP, SEXP P0_includedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type J(JSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type P0_iter(P0_iterSEXP);
    Rcpp::traits::input_parameter< bool >::type P0_included(P0_includedSEXP);
    rcpp_result_gen = Rcpp::wrap(RGHD(J, m, r, q, P0_iter, P0_included));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SkeweDF_Generalized_Pareto_calc_function", (DL_FUNC) &_SkeweDF_Generalized_Pareto_calc_function, 4},
    {"_SkeweDF_Generalized_Pareto_calc_pi", (DL_FUNC) &_SkeweDF_Generalized_Pareto_calc_pi, 5},
    {"_SkeweDF_Generalized_Pareto_calc_P0_delta", (DL_FUNC) &_SkeweDF_Generalized_Pareto_calc_P0_delta, 5},
    {"_SkeweDF_Generalized_Pareto_calc_P0_iter", (DL_FUNC) &_SkeweDF_Generalized_Pareto_calc_P0_iter, 5},
    {"_SkeweDF_Generalized_Pareto", (DL_FUNC) &_SkeweDF_Generalized_Pareto, 5},
    {"_SkeweDF_Kolmogorov_Waring_P0_calc", (DL_FUNC) &_SkeweDF_Kolmogorov_Waring_P0_calc, 3},
    {"_SkeweDF_Kolmogorov_Waring_P0", (DL_FUNC) &_SkeweDF_Kolmogorov_Waring_P0, 3},
    {"_SkeweDF_Kolmogorov_Waring", (DL_FUNC) &_SkeweDF_Kolmogorov_Waring, 4},
    {"_SkeweDF_RGHD_P0_calc_function", (DL_FUNC) &_SkeweDF_RGHD_P0_calc_function, 4},
    {"_SkeweDF_RGHD_P0_calc_pi", (DL_FUNC) &_SkeweDF_RGHD_P0_calc_pi, 4},
    {"_SkeweDF_RGHD_P0_calc", (DL_FUNC) &_SkeweDF_RGHD_P0_calc, 4},
    {"_SkeweDF_RGHD_P0", (DL_FUNC) &_SkeweDF_RGHD_P0, 3},
    {"_SkeweDF_RGHD", (DL_FUNC) &_SkeweDF_RGHD, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_SkeweDF(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
