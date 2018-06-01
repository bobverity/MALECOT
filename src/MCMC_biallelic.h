
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
  //lambda;
  int COI_model;
  int COI_max;
  //COI_dispersion;
  double e1;
  double e2;
  bool estimate_error;
  //e1_max;
  //e2_max;
  int burnin;
  int samples;
  int rungs;
  bool auto_converge;
  bool coupling_on;
  bool scaffold_on;
  int scaffold_n;
  bool split_merge_on;
  bool solve_label_switching;
  double precision;
  int precision_size;
  double GTI_pow;
  int report_iteration;
  bool flush_console;
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
  std::vector<std::vector<std::vector<int>>> scaf_group;
  
  // ordering of labels
  std::vector<int> label_order;
  std::vector<int> label_order_new;
  
  // objects for storing results
  std::vector<std::vector<double>> burnin_loglike;
  std::vector<std::vector<double>> sampling_loglike;
  std::vector<std::vector<double>> log_qmatrix_running;
  std::vector<std::vector<double>> qmatrix_final;
  
  // objects for storing acceptance rates
  std::vector<std::vector<int>> p_accept;
  int e1_accept;
  int e2_accept;
  std::vector<int> coupling_accept;
  std::vector<int> scaf_accept;
  std::vector<int> splitmerge_accept;
  
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
