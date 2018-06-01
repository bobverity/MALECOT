
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// run biallelic MCMC
Rcpp::List run_mcmc_biallelic_cpp(Rcpp::List args);

//------------------------------------------------
// run multiallelic MCMC
Rcpp::List run_mcmc_multiallelic_cpp(Rcpp::List args);
