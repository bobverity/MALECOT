
#include <chrono>
#include "main.h"
#include "misc.h"
#include "probability.h"
#include "MCMC_biallelic.h"
#include "Hungarian.h"

using namespace std;

//------------------------------------------------
// generate scaffolds under biallelic model
// [[Rcpp::export]]
Rcpp::List generate_scaffolds_biallelic_cpp(Rcpp::List args) {
  
  // extract arguments
  int K = rcpp_to_int(args["K"]);
  bool silent = rcpp_to_bool(args["silent"]);
  
  // start program
  if (!silent) {
    print("Generating scaffolds for K =", K);
  }
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // create MCMC object
  MCMC_biallelic m(args);
  
  // generate scaffold groupings
  m.scaffold_mcmc(args);
  
  // end timer
  chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
  chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
  if (!silent) {
    print("   completed in", time_span.count(), "seconds");
  }
  
  // create return object
  Rcpp::List ret;
  ret.push_back(Rcpp::wrap( m.scaffold_group ));
  ret.push_back(Rcpp::wrap( m.scaffold_loglike ));
  
  Rcpp::StringVector ret_names;
  ret_names.push_back("scaffold_group");
  ret_names.push_back("scaffold_loglike");
  
  ret.names() = ret_names;
  return ret;
}

//------------------------------------------------
// run biallelic MCMC
// [[Rcpp::export]]
Rcpp::List run_mcmc_biallelic_cpp(Rcpp::List args) {
  
  // extract arguments
  int K = rcpp_to_int(args["K"]);
  bool silent = rcpp_to_bool(args["silent"]);
  
  // start program
  if (!silent) {
    print("Running MCMC for K =", K);
  }
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // create MCMC object
  MCMC_biallelic m(args);
  
  // run MCMC
  m.burnin_mcmc(args);
  m.sampling_mcmc(args);
  
  // end timer
  chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
  chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
  if (!silent) {
    print("   completed in", time_span.count(), "seconds");
  }
  
  // create return object
  Rcpp::List ret;
  ret.push_back(Rcpp::wrap( m.burnin_loglike ));
  ret.push_back(Rcpp::wrap( m.sampling_loglike ));
  ret.push_back(Rcpp::wrap( m.m_store ));
  ret.push_back(Rcpp::wrap( m.p_store ));
  ret.push_back(Rcpp::wrap( m.e1_store ));
  ret.push_back(Rcpp::wrap( m.e2_store ));
  ret.push_back(Rcpp::wrap( m.COI_mean_store ));
  ret.push_back(Rcpp::wrap( m.qmatrix_final ));
  ret.push_back(Rcpp::wrap( m.coupling_accept ));
  ret.push_back(Rcpp::wrap( m.scaf_trials ));
  ret.push_back(Rcpp::wrap( m.scaf_accept ));
  
  Rcpp::StringVector ret_names;
  ret_names.push_back("burnin_loglike");
  ret_names.push_back("sampling_loglike");
  ret_names.push_back("m_store");
  ret_names.push_back("p_store");
  ret_names.push_back("e1_store");
  ret_names.push_back("e2_store");
  ret_names.push_back("COI_mean_store");
  ret_names.push_back("q_matrix");
  ret_names.push_back("coupling_accept");
  ret_names.push_back("scaf_trials");
  ret_names.push_back("scaf_accept");
  
  ret.names() = ret_names;
  return ret;
}

//------------------------------------------------
// generate scaffolds under multiallelic model
// [[Rcpp::export]]
Rcpp::List generate_scaffolds_multiallelic_cpp(Rcpp::List args) {
  
  // create return object
  Rcpp::List ret;
  ret.push_back(Rcpp::wrap( -9 ));
  
  Rcpp::StringVector ret_names;
  ret_names.push_back("foo");
  
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

//------------------------------------------------
// run Hungarian algorithm on given cost matrix
// [[Rcpp::export]]
Rcpp::List fix_labels_cpp(Rcpp::List args_model) {
  
  // extract cost matrix
  vector<vector<double>> cost_mat = rcpp_to_mat_double(args_model["cost_mat"]);
  
  // other objects needed for Hungarian algorithm
  int K = cost_mat.size();
  vector<int> best_perm(K);
  vector<int> edges_left(K);
  vector<int> edges_right(K);
  vector<int> blocked_left(K);
  vector<int> blocked_right(K);
  
  // find best permutation of current labels using Hungarian algorithm
  best_perm = hungarian(cost_mat, edges_left, edges_right, blocked_left, blocked_right);
  
  // return as Rcpp object
  return Rcpp::List::create(Rcpp::Named("best_perm")=best_perm);
}
