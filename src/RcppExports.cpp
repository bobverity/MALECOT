// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// run_mcmc_biallelic_cpp
Rcpp::List run_mcmc_biallelic_cpp(Rcpp::List args);
RcppExport SEXP _MALECOT_run_mcmc_biallelic_cpp(SEXP argsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type args(argsSEXP);
    rcpp_result_gen = Rcpp::wrap(run_mcmc_biallelic_cpp(args));
    return rcpp_result_gen;
END_RCPP
}
// run_mcmc_multiallelic_cpp
Rcpp::List run_mcmc_multiallelic_cpp(Rcpp::List args);
RcppExport SEXP _MALECOT_run_mcmc_multiallelic_cpp(SEXP argsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type args(argsSEXP);
    rcpp_result_gen = Rcpp::wrap(run_mcmc_multiallelic_cpp(args));
    return rcpp_result_gen;
END_RCPP
}
// GTI_posterior_K_sim_cpp
Rcpp::List GTI_posterior_K_sim_cpp(Rcpp::List args);
RcppExport SEXP _MALECOT_GTI_posterior_K_sim_cpp(SEXP argsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type args(argsSEXP);
    rcpp_result_gen = Rcpp::wrap(GTI_posterior_K_sim_cpp(args));
    return rcpp_result_gen;
END_RCPP
}
// GTI_integrated_K_sim_cpp
Rcpp::List GTI_integrated_K_sim_cpp(Rcpp::List args);
RcppExport SEXP _MALECOT_GTI_integrated_K_sim_cpp(SEXP argsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type args(argsSEXP);
    rcpp_result_gen = Rcpp::wrap(GTI_integrated_K_sim_cpp(args));
    return rcpp_result_gen;
END_RCPP
}
// fix_labels_cpp
Rcpp::List fix_labels_cpp(Rcpp::List args_model);
RcppExport SEXP _MALECOT_fix_labels_cpp(SEXP args_modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type args_model(args_modelSEXP);
    rcpp_result_gen = Rcpp::wrap(fix_labels_cpp(args_model));
    return rcpp_result_gen;
END_RCPP
}
// GTI_evidence_sim_cpp
Rcpp::List GTI_evidence_sim_cpp(Rcpp::List args);
RcppExport SEXP _MALECOT_GTI_evidence_sim_cpp(SEXP argsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type args(argsSEXP);
    rcpp_result_gen = Rcpp::wrap(GTI_evidence_sim_cpp(args));
    return rcpp_result_gen;
END_RCPP
}
// call_hungarian_cpp
Rcpp::List call_hungarian_cpp(Rcpp::List args);
RcppExport SEXP _MALECOT_call_hungarian_cpp(SEXP argsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type args(argsSEXP);
    rcpp_result_gen = Rcpp::wrap(call_hungarian_cpp(args));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MALECOT_run_mcmc_biallelic_cpp", (DL_FUNC) &_MALECOT_run_mcmc_biallelic_cpp, 1},
    {"_MALECOT_run_mcmc_multiallelic_cpp", (DL_FUNC) &_MALECOT_run_mcmc_multiallelic_cpp, 1},
    {"_MALECOT_GTI_posterior_K_sim_cpp", (DL_FUNC) &_MALECOT_GTI_posterior_K_sim_cpp, 1},
    {"_MALECOT_GTI_integrated_K_sim_cpp", (DL_FUNC) &_MALECOT_GTI_integrated_K_sim_cpp, 1},
    {"_MALECOT_fix_labels_cpp", (DL_FUNC) &_MALECOT_fix_labels_cpp, 1},
    {"_MALECOT_GTI_evidence_sim_cpp", (DL_FUNC) &_MALECOT_GTI_evidence_sim_cpp, 1},
    {"_MALECOT_call_hungarian_cpp", (DL_FUNC) &_MALECOT_call_hungarian_cpp, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_MALECOT(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
