
#include <chrono>
#include "main.h"
#include "misc.h"
#include "probability.h"
#include "MCMC_biallelic.h"
#include "Hungarian.h"

using namespace std;

//------------------------------------------------ new
// run biallelic MCMC
// [[Rcpp::export]]
Rcpp::List run_mcmc_biallelic_cpp(Rcpp::List args) {

  // extract arguments
  int K = rcpp_to_int(args["K"]);
  bool scaffold_on = rcpp_to_bool(args["scaffold_on"]);
  bool silent = rcpp_to_bool(args["silent"]);

  // start program
  if (!silent) {
    print("Running MCMC with K =", K);
  }
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // create MCMC object
  MCMC_biallelic m(args);
  
  // generate scaffold groupings
  if (scaffold_on) {
    m.scaffold_mcmc(args);
  }
  
  // run MCMC
  m.burnin_mcmc(args);
  m.sampling_mcmc(args);
  
  // end timer
  chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
  chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
  if (!silent) {
    print("   MCMC completed in", time_span.count(), "seconds");
  }
  
  // convert inputs to base C++ format
  //vector<double> x = rcpp_to_vector_double(args["x"]);
  //bool scaf_on = rcpp_to_bool(args["scaffold_on"]);
  /*
  // create MCMC object
  MCMC m(x, args);
  
  // generate scaffold groupings
  if (scaf_on) {
    m.scaffold_mcmc(args);
  }

  // run MCMC
  m.main_mcmc(args);

  // create return object
  Rcpp::List ret;
  ret.push_back(Rcpp::wrap( m.loglike_store ));
  ret.push_back(Rcpp::wrap( m.mu_store ));
  ret.push_back(Rcpp::wrap( m.qmatrix_final ));
  ret.push_back(Rcpp::wrap( m.mc_accept ));
  ret.push_back(Rcpp::wrap( m.scaf_accept ));
  ret.push_back(Rcpp::wrap( m.splitmerge_accept ));

  Rcpp::StringVector ret_names;
  ret_names.push_back("loglike");
  ret_names.push_back("mu");
  ret_names.push_back("qmatrix");
  ret_names.push_back("mc_accept");
  ret_names.push_back("scaf_accept");
  ret_names.push_back("splitmerge_accept");

  ret.names() = ret_names;
  return ret;
  */

  // create return object
  Rcpp::List ret;
  ret.push_back(Rcpp::wrap( m.burnin_loglike ));
  ret.push_back(Rcpp::wrap( m.sampling_loglike ));
  ret.push_back(Rcpp::wrap( m.qmatrix_final ));
  ret.push_back(Rcpp::wrap( m.coupling_accept ));

  Rcpp::StringVector ret_names;
  ret_names.push_back("burnin_loglike");
  ret_names.push_back("sampling_loglike");
  ret_names.push_back("qmatrix_final");
  ret_names.push_back("coupling_accept");

  ret.names() = ret_names;
  return ret;
}

//------------------------------------------------
// run multiallelic MCMC
// [[Rcpp::export]]
Rcpp::List run_mcmc_multiallelic_cpp(Rcpp::List args) {

  // create return object
  Rcpp::List ret;
  ret.push_back(Rcpp::wrap( -9 ));

  Rcpp::StringVector ret_names;
  ret_names.push_back("foo");

  ret.names() = ret_names;
  return ret;
}
