// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// epanechnikov
double epanechnikov(double x);
RcppExport SEXP _logKDE_epanechnikov(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(epanechnikov(x));
    return rcpp_result_gen;
END_RCPP
}
// gaussian
double gaussian(double x);
RcppExport SEXP _logKDE_gaussian(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(gaussian(x));
    return rcpp_result_gen;
END_RCPP
}
// laplace
double laplace(double x);
RcppExport SEXP _logKDE_laplace(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(laplace(x));
    return rcpp_result_gen;
END_RCPP
}
// logistic
double logistic(double x);
RcppExport SEXP _logKDE_logistic(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(logistic(x));
    return rcpp_result_gen;
END_RCPP
}
// triangular
double triangular(double x);
RcppExport SEXP _logKDE_triangular(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(triangular(x));
    return rcpp_result_gen;
END_RCPP
}
// uniform
double uniform(double x);
RcppExport SEXP _logKDE_uniform(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(uniform(x));
    return rcpp_result_gen;
END_RCPP
}
// KDE_e
double KDE_e(double x, std::vector<double> xi, double h);
RcppExport SEXP _logKDE_KDE_e(SEXP xSEXP, SEXP xiSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(KDE_e(x, xi, h));
    return rcpp_result_gen;
END_RCPP
}
// KDE_g
double KDE_g(double x, std::vector<double> xi, double h);
RcppExport SEXP _logKDE_KDE_g(SEXP xSEXP, SEXP xiSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(KDE_g(x, xi, h));
    return rcpp_result_gen;
END_RCPP
}
// KDE_lo
double KDE_lo(double x, std::vector<double> xi, double h);
RcppExport SEXP _logKDE_KDE_lo(SEXP xSEXP, SEXP xiSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(KDE_lo(x, xi, h));
    return rcpp_result_gen;
END_RCPP
}
// KDE_t
double KDE_t(double x, std::vector<double> xi, double h);
RcppExport SEXP _logKDE_KDE_t(SEXP xSEXP, SEXP xiSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(KDE_t(x, xi, h));
    return rcpp_result_gen;
END_RCPP
}
// KDE_la
double KDE_la(double x, std::vector<double> xi, double h);
RcppExport SEXP _logKDE_KDE_la(SEXP xSEXP, SEXP xiSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(KDE_la(x, xi, h));
    return rcpp_result_gen;
END_RCPP
}
// KDE_u
double KDE_u(double x, std::vector<double> xi, double h);
RcppExport SEXP _logKDE_KDE_u(SEXP xSEXP, SEXP xiSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(KDE_u(x, xi, h));
    return rcpp_result_gen;
END_RCPP
}
// silverman
double silverman(std::vector<double> x);
RcppExport SEXP _logKDE_silverman(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(silverman(x));
    return rcpp_result_gen;
END_RCPP
}
// logKDE
std::vector<double> logKDE(const std::vector<double>& input, const std::vector<double>& support, double h, std::string method);
RcppExport SEXP _logKDE_logKDE(SEXP inputSEXP, SEXP supportSEXP, SEXP hSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type input(inputSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type support(supportSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< std::string >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(logKDE(input, support, h, method));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_logKDE_epanechnikov", (DL_FUNC) &_logKDE_epanechnikov, 1},
    {"_logKDE_gaussian", (DL_FUNC) &_logKDE_gaussian, 1},
    {"_logKDE_laplace", (DL_FUNC) &_logKDE_laplace, 1},
    {"_logKDE_logistic", (DL_FUNC) &_logKDE_logistic, 1},
    {"_logKDE_triangular", (DL_FUNC) &_logKDE_triangular, 1},
    {"_logKDE_uniform", (DL_FUNC) &_logKDE_uniform, 1},
    {"_logKDE_KDE_e", (DL_FUNC) &_logKDE_KDE_e, 3},
    {"_logKDE_KDE_g", (DL_FUNC) &_logKDE_KDE_g, 3},
    {"_logKDE_KDE_lo", (DL_FUNC) &_logKDE_KDE_lo, 3},
    {"_logKDE_KDE_t", (DL_FUNC) &_logKDE_KDE_t, 3},
    {"_logKDE_KDE_la", (DL_FUNC) &_logKDE_KDE_la, 3},
    {"_logKDE_KDE_u", (DL_FUNC) &_logKDE_KDE_u, 3},
    {"_logKDE_silverman", (DL_FUNC) &_logKDE_silverman, 1},
    {"_logKDE_logKDE", (DL_FUNC) &_logKDE_logKDE, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_logKDE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
