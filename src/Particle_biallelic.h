
#pragma once

#include <Rcpp.h>
#include "Data.h"
#include "Lookup.h"

//------------------------------------------------
// class defining particle for biallelic case
class Particle_biallelic : public Data_biallelic, public Lookup {
  
public:
  // PUBLIC OBJECTS
  
  // beta_raised stores values of beta (the thermodynamic power), raised to the 
  // power GTI_pow
  double beta_raised;
  
  // proposal standard deviations
  double e_propSD;
  std::vector< std::vector<double> > p_propSD;
  std::vector<double> m_prop_mean;
  std::vector<double> COI_mean_propSD;
  std::vector<double> COI_mean_propSD_v2;
  
  // COI objects
  std::vector<double> COI_mean_vec;
  std::vector<double> COI_mean_shape;
  std::vector<double> COI_mean_rate;
  
  // grouping
  std::vector<int> group;
  
  // likelihood
  std::vector<std::vector<double>> loglike_old;
  std::vector<std::vector<double>> loglike_new;
  std::vector<std::vector<double>> loglike_new_group;
  double loglike;
  
  // COI and allele frequencies
  std::vector<int> m;
  std::vector<std::vector<double>> p;
  
  // qmatrices
  std::vector<std::vector<double>> log_qmatrix;
  std::vector<std::vector<double>> qmatrix;
  
  // probability vectors and matrices used when updating draws
  std::vector<double> sum_loglike_old_vec;
  std::vector<double> sum_loglike_new_vec;
  std::vector<std::vector<double>> p_prop;
  std::vector<int> m_prop;
  std::vector<double> COI_mean_prop;
  
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
  std::vector<int> m_accept;
  int e_accept;
  std::vector<int> COI_mean_accept;
  std::vector<int> COI_mean_accept_v2;
  
  // function pointers
  double (Particle_biallelic::*logprob_genotype_ptr) (int, double, int, double, double);
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Particle_biallelic() {};
  Particle_biallelic(double beta_raised);
  
  // other functions
  double logprob_genotype(int S, double p, int m, double e1, double e2);
  double logprob_genotype_exact(int S, double p, int m, double e1, double e2);
  double logprob_genotype_lookup(int S, double p, int m, double e1, double e2);
  double logprob_genotype_lookup_error(int S, double p, int m, double e1, double e2);
  
  void reset(double beta_raised);
  double get_lambda(int i, int j);
  void update_p(bool robbins_monro_on, int iteration);
  void update_m(bool robbins_monro_on, int iteration);
  void update_group();
  void update_e(bool robbins_monro_on, int iteration);
  void update_COI_mean(bool robbins_monro_on, int iteration);
  void update_COI_mean_v2(bool robbins_monro_on, int iteration);
  void calculate_loglike();
  void solve_label_switching(const std::vector<std::vector<double>> &log_qmatrix_running);
  
};
