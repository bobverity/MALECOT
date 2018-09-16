
#pragma once

#include <Rcpp.h>
#include "Data.h"
#include "Parameters.h"
#include "Lookup.h"
#include "Particle_biallelic.h"

//------------------------------------------------
// class defining MCMC for biallelic case
class MCMC_biallelic : public Data_biallelic, public Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // thermodynamic parameters
  std::vector<double> beta_raised_vec;
  std::vector<int> rung_order;
  int cold_rung;
  
  // vector of particles
  std::vector<Particle_biallelic> particle_vec;
  
  // ordering of labels
  std::vector<int> label_order;
  
  // Q-matrices
  std::vector<std::vector<double>> log_qmatrix_running;
  std::vector<std::vector<double>> qmatrix_final;
  
  // objects for storing results
  std::vector<std::vector<double>> loglike_burnin;
  std::vector<std::vector<double>> loglike_sampling;
  std::vector<std::vector<int>> m_store;
  std::vector<std::vector<std::vector<double>>> p_store;
  std::vector<double> e1_store;
  std::vector<double> e2_store;
  std::vector<std::vector<double>> COI_mean_store;
  
  // objects for storing acceptance rates
  std::vector<std::vector<int>> p_accept;
  int e_accept;
  std::vector<int> coupling_accept;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  MCMC_biallelic();
  
  // other functions
  void burnin_mcmc(Rcpp::List &args_functions, Rcpp::List &args_progress);
  void sampling_mcmc(Rcpp::List &args_functions, Rcpp::List &args_progress);
  void update_qmatrix_running();
  void update_qmatrix_final();
  void metropolis_coupling();
};
