
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// generate scaffolds under biallelic model
Rcpp::List generate_scaffolds_biallelic_cpp(Rcpp::List args);

//------------------------------------------------
// run biallelic MCMC
Rcpp::List run_mcmc_biallelic_cpp(Rcpp::List args);

//------------------------------------------------
// generate scaffolds under multiallelic model
Rcpp::List generate_scaffolds_multiallelic_cpp(Rcpp::List args);

//------------------------------------------------
// run multiallelic MCMC
Rcpp::List run_mcmc_multiallelic_cpp(Rcpp::List args);

//------------------------------------------------
// run Hungarian algorithm on given cost matrix
Rcpp::List fix_labels_cpp(Rcpp::List args_model);

//------------------------------------------------
// estimate evidence quantiles by simulation
Rcpp::List GTI_evidence_sim_cpp(Rcpp::List args);

