
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// class defining particle for biallelic case
class particle_biallelic {

public:
  
  // PUBLIC OBJECTS
  
  // pointers to observed data and basic data properties
  std::vector<std::vector<int>> *data_ptr;
  int n;
  int L;
  
  // pointers to lookup tables
  std::vector<std::vector<double>> *lookup_homo_ptr;
  std::vector<std::vector<double>> *lookup_het_ptr;
  
  // model and MCMC parameters
  int K;
  double precision;
  int precision_size;
  bool estimate_error;
  int COI_model;
  int COI_max;
  double COI_dispersion;
  double beta_raised;
  int scaffold_n;
  double e1;
  double e2;
  double e1_max;
  double e2_max;
  
  // proposal standard deviations
  double e1_propSD;
  double e2_propSD;
  std::vector< std::vector<double> > p_propSD;
  
  // COI objects
  std::vector<double> COI_mean;
  std::vector<double> COI_mean_shape;
  std::vector<double> COI_mean_rate;
  std::vector<double> COI_mean_propSD;
  
  // grouping
  std::vector<int> group;
  std::vector<int> group_order;
  
  // likelihood
  std::vector<std::vector<double>> loglike_old;
  std::vector<std::vector<double>> loglike_new;
  std::vector<std::vector<double>> loglike_new_group;
  double loglike;
  
  // COI and allele frequencies
  std::vector<int> m;
  std::vector<std::vector<double>> p;
  
  // probability vectors and matrices used when updating draws
  std::vector<std::vector<double>> log_qmatrix;
  std::vector<std::vector<double>> qmatrix;
  
  std::vector<double> sum_loglike_old_vec;
  std::vector<double> sum_loglike_new_vec;
  std::vector<double> p_prop;
  std::vector<double> COI_mean_prop;
  
  
  // Q matrix
  //std::vector<std::vector<double>> qmatrix;
  //std::vector<std::vector<double>> log_qmatrix;
  
  // initialise ordering of labels
  std::vector<int> label_order;
  std::vector<int> label_order_new;
  /*
  // objects for solving label switching problem
  std::vector<std::vector<double>> cost_mat;
  std::vector<int> best_perm;
  std::vector<int> best_perm_order;
  std::vector<int> edges_left;
  std::vector<int> edges_right;
  std::vector<int> blocked_left;
  std::vector<int> blocked_right;
  */
  // scaffold objects
  std::vector<std::vector<int>> scaf_group;
  std::vector<int> scaf_count;
  /*
  // objects for split-merge
  std::vector<int> splitmerge_targets;
  std::vector<double> splitmerge_mu;
  std::vector<int> splitmerge_group;
  std::vector<int> splitmerge_counts;
  std::vector<double> splitmerge_x_sum;
  */
  
  // store acceptance rates
  std::vector<std::vector<int>> p_accept;
  int e1_accept;
  int e2_accept;
  
  // function pointers
  double (particle_biallelic::*logprob_genotype_ptr) (int, double, int, double, double);
  
  // PUBLIC FUNCTIONS
  
  // constructors
  particle_biallelic();
  particle_biallelic(std::vector<std::vector<int>> &data, Rcpp::List &args, std::vector<std::vector<double>> &lookup_homo, std::vector<std::vector<double>> &lookup_het, double beta_raised_);
  
  // other functions
  void reset();
  void update_e(int which_e, bool robbins_monro_on, int iteration);
  void update_p(bool robbins_monro_on, int iteration);
  void update_m();
  void update_group();
  void update_COI_mean(bool robbins_monro_on, int iteration);
  void calculate_loglike();
  
  double logprob_genotype(int S, double p, int m, double e1, double e2);
  double logprob_genotype_exact(int S, double p, int m, double e1, double e2);
  double logprob_genotype_lookup(int S, double p, int m, double e1, double e2);
  double logprob_genotype_lookup_varE(int S, double p, int m, double e1, double e2);
  
  void get_group_order();
  double scaf_prop_logprob(const std::vector<int> &prop_group);
  void scaf_propose(int &scaf_accept);
  void splitmerge_propose(int &splitmerge_accept);
  void solve_label_switching(const std::vector<std::vector<double>> &log_qmatrix_running);
  
};
