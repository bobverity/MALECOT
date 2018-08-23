
#pragma once

#include <Rcpp.h>
#include "Parameters.h"

//------------------------------------------------
// class containing lookup tables
class Lookup : public Parameters {
  
public:
  
  // lookup tables
  static std::vector< std::vector<double> > lookup_homo;
  static std::vector< std::vector<double> > lookup_het;
  static std::vector< std::vector<double> > lookup_lgamma;
  
  // constructor
  Lookup() {};
  
  // public methods
  void init_homohet();
  void init_lgamma();
  
};

