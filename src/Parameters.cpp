
#include "Parameters.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// declare static member variables

// MCMC parameters
int Parameters::K;
double Parameters::precision;
int Parameters::precision_size;
int Parameters::burnin;
int Parameters::samples;
int Parameters::rungs;
double Parameters::GTI_pow;
bool Parameters::auto_converge;
int Parameters::converge_test;
bool Parameters::solve_label_switching_on;
bool Parameters::coupling_on;
bool Parameters::pb_markdown;
bool Parameters::silent;

// model parameters
std::vector<std::vector<double>> Parameters::lambda;
int Parameters::COI_model;
int Parameters::COI_max;
std::vector<int> Parameters::COI_manual;
bool Parameters::estimate_COI_mean;
double Parameters::COI_mean;
double Parameters::COI_dispersion;
bool Parameters::estimate_error;
double Parameters::e1;
double Parameters::e2;
double Parameters::e1_max;
double Parameters::e2_max;

//------------------------------------------------
// constructor for parameters class
Parameters::Parameters(const Rcpp::List &args) {
  
  // MCMC parameters
  K = rcpp_to_int(args["K"]);
  precision = rcpp_to_double(args["precision"]);
  precision_size = (int)(1/precision);
  burnin = rcpp_to_int(args["burnin"]);
  samples = rcpp_to_int(args["samples"]);
  rungs = rcpp_to_int(args["rungs"]);
  GTI_pow = rcpp_to_double(args["GTI_pow"]);
  auto_converge = rcpp_to_bool(args["auto_converge"]);
  converge_test = rcpp_to_int(args["converge_test"]);
  solve_label_switching_on = rcpp_to_bool(args["solve_label_switching_on"]);
  coupling_on = rcpp_to_bool(args["coupling_on"]);
  pb_markdown = rcpp_to_bool(args["pb_markdown"]);
  silent = rcpp_to_bool(args["silent"]);
  
  // model parameters
  lambda = rcpp_to_mat_double(args["lambda"]);
  COI_model = rcpp_to_int(args["COI_model_numeric"]);
  COI_max = rcpp_to_int(args["COI_max"]);
  COI_manual = rcpp_to_vector_int(args["COI_manual"]);
  estimate_COI_mean = rcpp_to_bool(args["estimate_COI_mean"]);
  COI_mean = rcpp_to_double(args["COI_mean"]);
  COI_dispersion = rcpp_to_double(args["COI_dispersion"]);
  estimate_error = rcpp_to_bool(args["estimate_error"]);
  e1 = rcpp_to_double(args["e1"]);
  e2 = rcpp_to_double(args["e2"]);
  e1_max = rcpp_to_double(args["e1_max"]);
  e2_max = rcpp_to_double(args["e2_max"]);
  
}
