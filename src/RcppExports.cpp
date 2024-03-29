// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// gradIPS
arma::vec gradIPS(arma::vec b, arma::vec d, arma::mat& X, arma::mat& w, double treated_flag, arma::vec whs);
RcppExport SEXP _IPS_gradIPS(SEXP bSEXP, SEXP dSEXP, SEXP XSEXP, SEXP wSEXP, SEXP treated_flagSEXP, SEXP whsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type treated_flag(treated_flagSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type whs(whsSEXP);
    rcpp_result_gen = Rcpp::wrap(gradIPS(b, d, X, w, treated_flag, whs));
    return rcpp_result_gen;
END_RCPP
}
// gradLIPS
arma::vec gradLIPS(arma::vec b, arma::vec d, arma::vec z, arma::mat& X, arma::mat& w, arma::vec whs);
RcppExport SEXP _IPS_gradLIPS(SEXP bSEXP, SEXP dSEXP, SEXP zSEXP, SEXP XSEXP, SEXP wSEXP, SEXP whsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type whs(whsSEXP);
    rcpp_result_gen = Rcpp::wrap(gradLIPS(b, d, z, X, w, whs));
    return rcpp_result_gen;
END_RCPP
}
// linIPS
arma::mat linIPS(arma::vec bhat, arma::vec d, arma::vec pshat, arma::mat& X, arma::mat& w, double treated_flag, arma::vec whs);
RcppExport SEXP _IPS_linIPS(SEXP bhatSEXP, SEXP dSEXP, SEXP pshatSEXP, SEXP XSEXP, SEXP wSEXP, SEXP treated_flagSEXP, SEXP whsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type bhat(bhatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pshat(pshatSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type treated_flag(treated_flagSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type whs(whsSEXP);
    rcpp_result_gen = Rcpp::wrap(linIPS(bhat, d, pshat, X, w, treated_flag, whs));
    return rcpp_result_gen;
END_RCPP
}
// linLIPS
arma::mat linLIPS(arma::vec bhat, arma::vec d, arma::vec z, arma::vec ipshat, arma::mat& X, arma::mat& w, arma::vec whs);
RcppExport SEXP _IPS_linLIPS(SEXP bhatSEXP, SEXP dSEXP, SEXP zSEXP, SEXP ipshatSEXP, SEXP XSEXP, SEXP wSEXP, SEXP whsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type bhat(bhatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ipshat(ipshatSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type whs(whsSEXP);
    rcpp_result_gen = Rcpp::wrap(linLIPS(bhat, d, z, ipshat, X, w, whs));
    return rcpp_result_gen;
END_RCPP
}
// objIPS
double objIPS(arma::vec b, arma::vec d, arma::mat X, arma::mat w, double treated_flag, arma::vec whs);
RcppExport SEXP _IPS_objIPS(SEXP bSEXP, SEXP dSEXP, SEXP XSEXP, SEXP wSEXP, SEXP treated_flagSEXP, SEXP whsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type treated_flag(treated_flagSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type whs(whsSEXP);
    rcpp_result_gen = Rcpp::wrap(objIPS(b, d, X, w, treated_flag, whs));
    return rcpp_result_gen;
END_RCPP
}
// objLIPS
double objLIPS(arma::vec b, arma::vec d, arma::vec z, arma::mat& X, arma::mat& w, arma::vec whs);
RcppExport SEXP _IPS_objLIPS(SEXP bSEXP, SEXP dSEXP, SEXP zSEXP, SEXP XSEXP, SEXP wSEXP, SEXP whsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type whs(whsSEXP);
    rcpp_result_gen = Rcpp::wrap(objLIPS(b, d, z, X, w, whs));
    return rcpp_result_gen;
END_RCPP
}
// weightIPSexp
arma::mat weightIPSexp(arma::mat X, String X_trans);
RcppExport SEXP _IPS_weightIPSexp(SEXP XSEXP, SEXP X_transSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< String >::type X_trans(X_transSEXP);
    rcpp_result_gen = Rcpp::wrap(weightIPSexp(X, X_trans));
    return rcpp_result_gen;
END_RCPP
}
// weightIPSind
arma::mat weightIPSind(arma::mat X);
RcppExport SEXP _IPS_weightIPSind(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(weightIPSind(X));
    return rcpp_result_gen;
END_RCPP
}
// weightIPSproj
NumericMatrix weightIPSproj(const NumericMatrix& X);
RcppExport SEXP _IPS_weightIPSproj(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(weightIPSproj(X));
    return rcpp_result_gen;
END_RCPP
}
// weightIPSproj_uniq
NumericMatrix weightIPSproj_uniq(const NumericMatrix& X, const NumericVector& wgt);
RcppExport SEXP _IPS_weightIPSproj_uniq(SEXP XSEXP, SEXP wgtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type wgt(wgtSEXP);
    rcpp_result_gen = Rcpp::wrap(weightIPSproj_uniq(X, wgt));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_IPS_gradIPS", (DL_FUNC) &_IPS_gradIPS, 6},
    {"_IPS_gradLIPS", (DL_FUNC) &_IPS_gradLIPS, 6},
    {"_IPS_linIPS", (DL_FUNC) &_IPS_linIPS, 7},
    {"_IPS_linLIPS", (DL_FUNC) &_IPS_linLIPS, 7},
    {"_IPS_objIPS", (DL_FUNC) &_IPS_objIPS, 6},
    {"_IPS_objLIPS", (DL_FUNC) &_IPS_objLIPS, 6},
    {"_IPS_weightIPSexp", (DL_FUNC) &_IPS_weightIPSexp, 2},
    {"_IPS_weightIPSind", (DL_FUNC) &_IPS_weightIPSind, 1},
    {"_IPS_weightIPSproj", (DL_FUNC) &_IPS_weightIPSproj, 1},
    {"_IPS_weightIPSproj_uniq", (DL_FUNC) &_IPS_weightIPSproj_uniq, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_IPS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
