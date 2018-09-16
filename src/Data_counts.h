
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// class containing data in count format
class Data_counts {
  
public:
  
  std::vector<int> counts;
  std::vector<double> logprob;
  int K;
  int L;
  int COI_max;
  
  // constructors
  Data_counts() {};
  
  // other functions
  void init(int K, int L, int COI_max);
  void populate(const std::vector<std::vector<int>> &data, const std::vector<int> &group, const std::vector<int> &m);
  int get_counts(int k, int l, int m, int S);
  void set_logprob(int k, int l, int m, int S, double x);
  double get_logprob(int k, int l, int m, int S);
  
};

