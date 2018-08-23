
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// class containing bi-allelic data
class Data_biallelic {
  
public:
  
  static std::vector<std::vector<int>> data;
  static int n;
  static int L;
  
  // constructors
  Data_biallelic() {};
  Data_biallelic(const Rcpp::List &args);
  
};

//------------------------------------------------
// class containing multi-allelic data
class Data_multiallelic {
  
public:
  
  static std::vector<std::vector<std::vector<int>>> data;
  static std::vector<int> alleles;
  static int n;
  static int L;
  
  // constructors
  Data_multiallelic() {};
  Data_multiallelic(const Rcpp::List &args);
  
};
