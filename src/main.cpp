
#include <chrono>
#include "main.h"
#include "Parameters.h"
#include "Data.h"
#include "Lookup.h"
#include "misc_v1.h"
#include "probability.h"
#include "MCMC_biallelic.h"
#include "MCMC_multiallelic.h"
#include "hungarian.h"

using namespace std;

//------------------------------------------------
// run biallelic MCMC
// [[Rcpp::export]]
Rcpp::List run_mcmc_biallelic_cpp(Rcpp::List args) {
  
  // split argument lists
  Rcpp::List args_data = args["args_data"];
  Rcpp::List args_model = args["args_model"];
  Rcpp::List args_functions = args["args_functions"];
  Rcpp::List args_progress = args["args_progress"];
  
  // read in parameters into separate class
  Parameters parameters(args_model);
  
  // read in data into separate class
  Data_biallelic data(args_data);
  
  // define look-up tables
  Lookup lookup;
  if (parameters.precision != 0) {
    lookup.init_homohet();
  }
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // create MCMC object
  MCMC_biallelic mcmc_biallelic;
  
  // run MCMC
  mcmc_biallelic.burnin_mcmc(args_functions, args_progress);
  mcmc_biallelic.sampling_mcmc(args_functions, args_progress);
  
  // end timer
  chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
  chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
  if (!parameters.silent) {
    print("   completed in", time_span.count(), "seconds\n");
  }
  
  // create return object
  Rcpp::List ret;
  ret.push_back(Rcpp::wrap( mcmc_biallelic.loglike_burnin ));
  ret.push_back(Rcpp::wrap( mcmc_biallelic.loglike_sampling ));
  ret.push_back(Rcpp::wrap( mcmc_biallelic.m_store ));
  ret.push_back(Rcpp::wrap( mcmc_biallelic.p_store ));
  ret.push_back(Rcpp::wrap( mcmc_biallelic.e1_store ));
  ret.push_back(Rcpp::wrap( mcmc_biallelic.e2_store ));
  ret.push_back(Rcpp::wrap( mcmc_biallelic.COI_mean_store ));
  ret.push_back(Rcpp::wrap( mcmc_biallelic.qmatrix_final ));
  ret.push_back(Rcpp::wrap( mcmc_biallelic.p_accept ));
  ret.push_back(Rcpp::wrap( mcmc_biallelic.m_accept ));
  ret.push_back(Rcpp::wrap( mcmc_biallelic.e_accept ));
  ret.push_back(Rcpp::wrap( mcmc_biallelic.coupling_accept ));
  ret.push_back(Rcpp::wrap( mcmc_biallelic.rung_converged ));
  
  Rcpp::StringVector ret_names;
  ret_names.push_back("loglike_burnin");
  ret_names.push_back("loglike_sampling");
  ret_names.push_back("m_store");
  ret_names.push_back("p_store");
  ret_names.push_back("e1_store");
  ret_names.push_back("e2_store");
  ret_names.push_back("COI_mean_store");
  ret_names.push_back("qmatrix");
  ret_names.push_back("p_accept");
  ret_names.push_back("m_accept");
  ret_names.push_back("e_accept");
  ret_names.push_back("coupling_accept");
  ret_names.push_back("rung_converged");
  
  ret.names() = ret_names;
  return ret;
}

//------------------------------------------------
// run multiallelic MCMC
// [[Rcpp::export]]
Rcpp::List run_mcmc_multiallelic_cpp(Rcpp::List args) {
  
  // split argument lists
  Rcpp::List args_data = args["args_data"];
  Rcpp::List args_model = args["args_model"];
  Rcpp::List args_functions = args["args_functions"];
  Rcpp::List args_progress = args["args_progress"];
  
  // read in parameters into separate class
  Parameters parameters(args_model);
  
  // read in data into separate class
  Data_multiallelic data(args_data);
  
  // define look-up tables
  Lookup lookup;
  lookup.init_lgamma();
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // create MCMC object
  MCMC_multiallelic mcmc_multiallelic;
  
  // run MCMC
  mcmc_multiallelic.burnin_mcmc(args_functions, args_progress);
  mcmc_multiallelic.sampling_mcmc(args_functions, args_progress);
  
  // end timer
  chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
  chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
  if (!parameters.silent) {
    print("   completed in", time_span.count(), "seconds\n");
  }
  
  // create return object
  Rcpp::List ret;
  ret.push_back(Rcpp::wrap( mcmc_multiallelic.loglike_burnin ));
  ret.push_back(Rcpp::wrap( mcmc_multiallelic.loglike_sampling ));
  ret.push_back(Rcpp::wrap( mcmc_multiallelic.m_store ));
  ret.push_back(Rcpp::wrap( mcmc_multiallelic.p_store ));
  ret.push_back(Rcpp::wrap( mcmc_multiallelic.COI_mean_store ));
  ret.push_back(Rcpp::wrap( mcmc_multiallelic.qmatrix_final ));
  ret.push_back(Rcpp::wrap( mcmc_multiallelic.p_accept ));
  ret.push_back(Rcpp::wrap( mcmc_multiallelic.m_accept ));
  ret.push_back(Rcpp::wrap( mcmc_multiallelic.coupling_accept ));
  ret.push_back(Rcpp::wrap( mcmc_multiallelic.rung_converged ));
  
  Rcpp::StringVector ret_names;
  ret_names.push_back("loglike_burnin");
  ret_names.push_back("loglike_sampling");
  ret_names.push_back("m_store");
  ret_names.push_back("p_store");
  ret_names.push_back("COI_mean_store");
  ret_names.push_back("qmatrix");
  ret_names.push_back("p_accept");
  ret_names.push_back("m_accept");
  ret_names.push_back("coupling_accept");
  ret_names.push_back("rung_converged");
  
  ret.names() = ret_names;
  return ret;
}

//------------------------------------------------
// estimate quantiles of posterior probability of K by simulation
// [[Rcpp::export]]
Rcpp::List GTI_posterior_K_sim_cpp(Rcpp::List args) {
  
  // extract arguments
  vector<double> m = rcpp_to_vector_double(args["mean"]);
  vector<double> s = rcpp_to_vector_double(args["SE"]);
  int reps = rcpp_to_int(args["reps"]);
  int K = m.size();
  
  // obtain normalised draws
  vector<double> y(K);
  double y_max = 0;
  vector<double> y_trans(K);
  double y_trans_sum = 0;
  double y_trans_inv_sum = 0;
  vector<vector<double>> ret(K, vector<double>(reps));
  for (int i=0; i<reps; i++) {
    for (int k=0; k<K; k++) {
      y[k] = rnorm1(m[k], s[k]);
      if (k==0 || y[k]>y_max) {
        y_max = y[k];
      }
    }
    y_trans_sum = 0;
    for (int k=0; k<K; k++) {
      y_trans[k] = exp(y[k]-y_max);
      y_trans_sum += y_trans[k];
    }
    y_trans_inv_sum = 1/y_trans_sum;
    for (int k=0; k<K; k++) {
      ret[k][i] = y_trans[k]*y_trans_inv_sum;
    }
  }
  
  // return as Rcpp object
  return Rcpp::List::create(Rcpp::Named("ret")=ret);
}

//------------------------------------------------
// integrate log-evidence over K by simulation
// [[Rcpp::export]]
Rcpp::List GTI_integrated_K_sim_cpp(Rcpp::List args) {
  
  // extract arguments
  vector<double> m = rcpp_to_vector_double(args["mean"]);
  vector<double> s = rcpp_to_vector_double(args["SE"]);
  int reps = rcpp_to_int(args["reps"]);
  int K = m.size();
  
  // obtain integrated draws
  vector<double> y(K);
  double y_max = 0;
  double y_trans_sum = 0;
  double x = 0;
  double x_sum = 0;
  double x_sum_squared = 0;
  for (int i=0; i<reps; i++) {
    for (int k=0; k<K; k++) {
      y[k] = rnorm1(m[k], s[k]);
      if (k==0 || y[k]>y_max) {
        y_max = y[k];
      }
    }
    y_trans_sum = 0;
    for (int k=0; k<K; k++) {
      y_trans_sum += exp(y[k]-y_max);
    }
    x = y_max + log(y_trans_sum) - log(K);
    x_sum += x;
    x_sum_squared += x*x;
  }
  double x_mean = x_sum/double(reps);
  double x_var = x_sum_squared/double(reps) - x_mean*x_mean;
  if (x_var<0) {
    x_var = 0;
  }
  
  // return as Rcpp object
  return Rcpp::List::create(Rcpp::Named("mean")=x_mean,
                            Rcpp::Named("SE")=x_var);
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

//------------------------------------------------
// estimate evidence quantiles by simulation
// [[Rcpp::export]]
Rcpp::List GTI_evidence_sim_cpp(Rcpp::List args) {
  
  // extract arguments
  vector<double> m = rcpp_to_vector_double(args["mean"]);
  vector<double> s = rcpp_to_vector_double(args["SE"]);
  int reps = rcpp_to_int(args["reps"]);
  int K = m.size();
  
  // obtain normalised draws
  vector<double> y(K);
  double y_max = 0;
  vector<double> y_trans(K);
  double y_trans_sum = 0;
  double y_trans_inv_sum = 0;
  vector<vector<double>> ret(K, vector<double>(reps));
  for (int i=0; i<reps; i++) {
    for (int k=0; k<K; k++) {
      y[k] = rnorm1(m[k], s[k]);
      if (k==0 || y[k]>y_max) {
        y_max = y[k];
      }
    }
    y_trans_sum = 0;
    for (int k=0; k<K; k++) {
      y_trans[k] = exp(y[k]-y_max);
      y_trans_sum += y_trans[k];
    }
    y_trans_inv_sum = 1/y_trans_sum;
    for (int k=0; k<K; k++) {
      ret[k][i] = y_trans[k]*y_trans_inv_sum;
    }
  }
  
  // return as Rcpp object
  return Rcpp::List::create(Rcpp::Named("ret")=ret);
}


