// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// assign_branchnum_cpp
List assign_branchnum_cpp(NumericVector furcations, LogicalVector is_tip);
RcppExport SEXP _treestruct_assign_branchnum_cpp(SEXP furcationsSEXP, SEXP is_tipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type furcations(furcationsSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type is_tip(is_tipSEXP);
    rcpp_result_gen = Rcpp::wrap(assign_branchnum_cpp(furcations, is_tip));
    return rcpp_result_gen;
END_RCPP
}
// calc_pathlen_cpp
NumericVector calc_pathlen_cpp(NumericVector len, NumericVector parent_idx);
RcppExport SEXP _treestruct_calc_pathlen_cpp(SEXP lenSEXP, SEXP parent_idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type len(lenSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type parent_idx(parent_idxSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_pathlen_cpp(len, parent_idx));
    return rcpp_result_gen;
END_RCPP
}
// calc_sa_above_cpp
NumericVector calc_sa_above_cpp(NumericVector sa, NumericVector parentrow);
RcppExport SEXP _treestruct_calc_sa_above_cpp(SEXP saSEXP, SEXP parentrowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type sa(saSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type parentrow(parentrowSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_sa_above_cpp(sa, parentrow));
    return rcpp_result_gen;
END_RCPP
}
// calc_total_x_above_internode_cpp
NumericVector calc_total_x_above_internode_cpp(NumericVector x, NumericVector parentrow);
RcppExport SEXP _treestruct_calc_total_x_above_internode_cpp(SEXP xSEXP, SEXP parentrowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type parentrow(parentrowSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_total_x_above_internode_cpp(x, parentrow));
    return rcpp_result_gen;
END_RCPP
}
// get_children
std::set<int> get_children(int thisrow, IntegerVector parentrow);
RcppExport SEXP _treestruct_get_children(SEXP thisrowSEXP, SEXP parentrowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type thisrow(thisrowSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type parentrow(parentrowSEXP);
    rcpp_result_gen = Rcpp::wrap(get_children(thisrow, parentrow));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_treestruct_assign_branchnum_cpp", (DL_FUNC) &_treestruct_assign_branchnum_cpp, 2},
    {"_treestruct_calc_pathlen_cpp", (DL_FUNC) &_treestruct_calc_pathlen_cpp, 2},
    {"_treestruct_calc_sa_above_cpp", (DL_FUNC) &_treestruct_calc_sa_above_cpp, 2},
    {"_treestruct_calc_total_x_above_internode_cpp", (DL_FUNC) &_treestruct_calc_total_x_above_internode_cpp, 2},
    {"_treestruct_get_children", (DL_FUNC) &_treestruct_get_children, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_treestruct(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
