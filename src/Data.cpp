
#include "Data.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// declare static member variables for class Data_biallelic

std::vector<std::vector<int>> Data_biallelic::data;
int Data_biallelic::n;
int Data_biallelic::L;

//------------------------------------------------
// constructor for Data_biallelic class
Data_biallelic::Data_biallelic(const Rcpp::List &args) {
  
  data = rcpp_to_mat_int(args["data"]);
  n = rcpp_to_int(args["n"]);
  L = rcpp_to_int(args["L"]);
}


//------------------------------------------------
// declare static member variables for class Data_multiallelic

std::vector<std::vector<int>> Data_multiallelic::data;
int Data_multiallelic::n;
int Data_multiallelic::L;

//------------------------------------------------
// constructor for Data_multiallelic class
Data_multiallelic::Data_multiallelic(const Rcpp::List &args) {
  
  data = rcpp_to_mat_int(args["data"]);
  n = rcpp_to_int(args["n"]);
  L = rcpp_to_int(args["L"]);
}
