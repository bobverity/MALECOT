
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// class containing all parameters
class Parameters {
  
public:
  
  // MCMC parameters
  static int K;
  static double precision;
  static int precision_size;
  static int burnin;
  static int samples;
  static int rungs;
  static double GTI_pow;
  static bool auto_converge;
  static int converge_test;
  static bool solve_label_switching_on;
  static bool coupling_on;
  static bool pb_markdown;
  static bool silent;
  
  // model parameters
  static std::vector<std::vector<double>> lambda;
  static bool lambda_scalar;
  static int COI_model;
  static int COI_max;
  static std::vector<int> COI_manual;
  static bool estimate_COI_mean;
  static double COI_mean;
  static double COI_dispersion;
  static bool estimate_error;
  static double e1;
  static double e2;
  static double e1_max;
  static double e2_max;
  
  // constructors
  Parameters() {};
  Parameters(const Rcpp::List &args);
  
};
