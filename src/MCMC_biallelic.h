
#pragma once

#include <Rcpp.h>
#include "particle_biallelic.h"

//------------------------------------------------
// class defining MCMC for biallelic case
class MCMC_biallelic {
  
public:
  
  // PUBLIC OBJECTS
  
  // extract data and parameters
  std::vector<std::vector<int>> data;
  int n;
  int L;
  int K;
  int COI_model;
  int COI_max;
  double e1;
  double e2;
  bool estimate_error;
  int burnin;
  int samples;
  int rungs;
  bool auto_converge;
  bool coupling_on;
  bool scaffold_on;
  int scaffold_n;
  int scaffold_group_n;
  std::vector<std::vector<int>> scaffold_group;
  bool split_merge_on;
  bool solve_label_switching_on;
  double precision;
  int precision_size;
  double GTI_pow;
  bool silent;
  //output_format = rcpp_to_int(args["output_format"]);
  
  // lookup tables
  std::vector< std::vector<double> > lookup_homo;
  std::vector< std::vector<double> > lookup_het;
  
  // thermodynamic parameters
  std::vector<double> beta_raised_vec;
  std::vector<int> rung_order;
  int cold_rung;
  
  // vector of particles
  std::vector<particle_biallelic> particle_vec;
  
  // scaffold objects
  std::vector<double> scaffold_loglike;
  
  // ordering of labels
  std::vector<int> label_order;
  std::vector<int> label_order_new;
  
  // Q-matrices
  std::vector<std::vector<double>> log_qmatrix_running;
  std::vector<std::vector<double>> qmatrix_final;
  
  // objects for storing results
  std::vector<std::vector<double>> burnin_loglike;
  std::vector<std::vector<double>> sampling_loglike;
  std::vector<std::vector<int>> m_store;
  std::vector<std::vector<std::vector<double>>> p_store;
  std::vector<double> e1_store;
  std::vector<double> e2_store;
  std::vector<std::vector<double>> COI_mean_store;
  
  // objects for storing acceptance rates
  std::vector<std::vector<int>> p_accept;
  int e1_accept;
  int e2_accept;
  std::vector<int> coupling_accept;
  int scaf_trials;
  int scaf_accept;
  std::vector<int> split_merge_accept;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  MCMC_biallelic(Rcpp::List &args);
  
  // other functions
  void scaffold_mcmc(Rcpp::List &args);
  void burnin_mcmc(Rcpp::List &args);
  void sampling_mcmc(Rcpp::List &args);
  void update_qmatrix_running();
  void update_qmatrix_final();
  void metropolis_coupling();
};
