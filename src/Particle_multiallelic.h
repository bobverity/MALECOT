
#pragma once

#include <Rcpp.h>
#include "Data.h"
#include "Lookup.h"

//------------------------------------------------
// class defining particle for multiallelic case
class Particle_multiallelic : public Data_multiallelic, public Lookup {

public:
  // PUBLIC OBJECTS
  
  // beta_raised stores values of beta (the thermodynamic power), raised to the 
  // power GTI_pow
  double beta_raised;
  
  // scalar lambda value
  double lambda0;
  
  // proposal standard deviations
  std::vector< std::vector<double> > p_propSD;
  
  // COI objects
  std::vector<double> COI_mean_vec;
  std::vector<double> COI_mean_shape;
  std::vector<double> COI_mean_rate;
  std::vector<double> COI_mean_propSD;
  
  // grouping
  std::vector<int> group;
  
  // likelihood
  std::vector<std::vector<double>> loglike_old;
  std::vector<std::vector<double>> loglike_new;
  std::vector<std::vector<double>> loglike_new_group;
  double loglike;
  
  // COI and allele frequencies
  std::vector<int> m;
  std::vector<std::vector<std::vector<double>>> p;
  std::vector<std::vector<std::vector<double>>> logp;
  
  // qmatrices
  std::vector<std::vector<double>> log_qmatrix;
  std::vector<std::vector<double>> qmatrix;
  
  // objects used when updating draws
  std::vector<double> v;
  std::vector<double> sum_loglike_old_vec;
  std::vector<double> sum_loglike_new_vec;
  std::vector<std::vector<std::vector<double>>> p_prop;
  std::vector<std::vector<std::vector<double>>> logp_prop;
  //std::vector<double> COI_mean_prop;
  
  // initialise ordering of labels
  std::vector<int> label_order;
  std::vector<int> label_order_new;
  
  // objects for solving label switching problem
  std::vector<std::vector<double>> cost_mat;
  std::vector<int> best_perm;
  std::vector<int> best_perm_order;
  std::vector<int> edges_left;
  std::vector<int> edges_right;
  std::vector<int> blocked_left;
  std::vector<int> blocked_right;
  
  // store acceptance rates
  std::vector<std::vector<int>> p_accept;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Particle_multiallelic() {};
  Particle_multiallelic(double beta_raised);
  
  // other functions
  double logprob_genotype(const std::vector<int> &x, const std::vector<double> &logp, int m);
  void reset();
  void update_p(bool robbins_monro_on, int iteration);
  void update_m();
  void update_group();
  void update_COI_mean(bool robbins_monro_on, int iteration);
  void calculate_loglike();
  void solve_label_switching(const std::vector<std::vector<double>> &log_qmatrix_running);
  
};
