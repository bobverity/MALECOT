
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// define very large/small numbers for catching overflow/underflow problems
#define OVERFLO   1e100
#define UNDERFLO   1e-100

//------------------------------------------------
// basic sum over elements in a vector (templated for different data types).
template<class TYPE>
TYPE sum(std::vector<TYPE> &x) {
    TYPE output = 0;
    for (int i=0; i<int(x.size()); i++)
        output += x[i];
    return output;
}

//------------------------------------------------
// mean of vector (templated for different data types)
template<class TYPE>
double mean(std::vector<TYPE> &x) {
    return sum(x)/double(x.size());
}

//------------------------------------------------
// min of vector (templated for different data types)
template<class TYPE>
TYPE min(std::vector<TYPE> x) {
    return *min_element(x.begin(), x.end());
}

//------------------------------------------------
// max of vector (templated for different data types)
template<class TYPE>
TYPE max(std::vector<TYPE> x) {
    return *max_element(x.begin(), x.end());
}

//------------------------------------------------
// push back multiple values to vector
template<class TYPE>
void push_back_multiple(std::vector<TYPE> &lhs, std::vector<TYPE> &rhs) {
    lhs.insert(lhs.end(), rhs.begin(), rhs.end());
}

//------------------------------------------------
// add two numbers together in log space. One number (but not both) is allowed to be -inf.
double logSum(double logA, double logB);

//------------------------------------------------
// subtract B from A in log space, i.e. do log(exp(logA)-exp(logB)). Note that logB cannot be greater than logA, if it is the function returns 0 and prints a warning.
double logSubtract(double logA, double logB);

//------------------------------------------------
// helper function for printing a single value or series of values (templated for different data types)
template<class TYPE>
void print(TYPE x) {
    Rcpp::Rcout << x << "\n";
    R_FlushConsole();
}
template<class TYPE1, class TYPE2>
void print(TYPE1 x1, TYPE2 x2) {
    Rcpp::Rcout << x1 << " " << x2 << "\n";
    R_FlushConsole();
}
template<class TYPE1, class TYPE2, class TYPE3>
void print(TYPE1 x1, TYPE2 x2, TYPE3 x3) {
    Rcpp::Rcout << x1 << " " << x2 << " " << x3 << "\n";
    R_FlushConsole();
}
template<class TYPE1, class TYPE2, class TYPE3, class TYPE4>
void print(TYPE1 x1, TYPE2 x2, TYPE3 x3, TYPE4 x4) {
    Rcpp::Rcout << x1 << " " << x2 << " " << x3 << " " << x4 << "\n";
    R_FlushConsole();
}
template<class TYPE1, class TYPE2, class TYPE3, class TYPE4, class TYPE5>
void print(TYPE1 x1, TYPE2 x2, TYPE3 x3, TYPE4 x4, TYPE5 x5) {
    Rcpp::Rcout << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << "\n";
    R_FlushConsole();
}
template<class TYPE1, class TYPE2, class TYPE3, class TYPE4, class TYPE5, class TYPE6>
void print(TYPE1 x1, TYPE2 x2, TYPE3 x3, TYPE4 x4, TYPE5 x5, TYPE6 x6) {
    Rcpp::Rcout << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << " " << x6 << "\n";
    R_FlushConsole();
}

//------------------------------------------------
// helper function for printing contents of a vector (templated for different data types)
template<class TYPE>
void printVector(std::vector<TYPE> &x) {
    for (int i=0; i<x.size(); i++) {
        Rcpp::Rcout << x[i] << " ";
    }
    Rcpp::Rcout << "\n";
    R_FlushConsole();
}

//------------------------------------------------
// helper function for printing contents of a matrix (templated for different data types)
template<class TYPE>
void printMatrix(std::vector< std::vector<TYPE> > &M) {
    for (int i=0; i<M.size(); i++) {
        for (int j=0; j<M[i].size(); j++) {
            Rcpp::Rcout << M[i][j] << " ";
        }
        Rcpp::Rcout << "\n";
    }
    Rcpp::Rcout << "\n";
    R_FlushConsole();
}

//------------------------------------------------
// helper function for printing contents of a 3D array (templated for different data types)
template<class TYPE>
void printArray(std::vector< std::vector< std::vector<TYPE> > > &x) {
    for (int i=0; i<x.size(); i++) {
        Rcpp::Rcout << "--- slice " << i+1 << " ---\n";
        for (int j=0; j<x[i].size(); j++) {
            for (int k=0; k<x[i][j].size(); k++) {
                Rcpp::Rcout << x[i][j][k] << " ";
            }
            Rcpp::Rcout << "\n";
        }
        Rcpp::Rcout << "\n";
    }
    Rcpp::Rcout << "\n";
    R_FlushConsole();
}

//------------------------------------------------
// print simple bar-graph composed of title followed by n stars
void printStars(std::string title="", int n=10);

//------------------------------------------------
// print "foo", with option number e.g. "foo2"
void foo(int n=0);

//------------------------------------------------
// print "bar", with option number e.g. "bar2"
void bar(int n=0);

//------------------------------------------------
// print "foobar", with option number e.g. "foobar2"
void foobar(int n=0);

//------------------------------------------------
// analogue of R function seq() for integers
std::vector<int> seq_int(int from, int to, int by=1);

//------------------------------------------------
// converts input from Rcpp::List format to bool format.
int Rcpp_to_bool(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to int format.
int Rcpp_to_int(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to double format.
double Rcpp_to_double(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to string format.
std::string Rcpp_to_string(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<int> format.
std::vector<int> Rcpp_to_vector_int(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<double> format.
std::vector<double> Rcpp_to_vector_double(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<string> format.
std::vector<std::string> Rcpp_to_vector_string(SEXP x);

// converts input from Rcpp::List format to vector<vector<bool>> format.
std::vector< std::vector<bool> > Rcpp_to_mat_bool(Rcpp::List x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<int>> format.
std::vector< std::vector<int> > Rcpp_to_mat_int(Rcpp::List x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<double>> format.
std::vector< std::vector<double> > Rcpp_to_mat_double(Rcpp::List x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<vector<double>>> format.
std::vector< std::vector< std::vector<double> > > Rcpp_to_array_double(Rcpp::List x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<vector<int>>> format.
std::vector< std::vector< std::vector<int> > > Rcpp_to_array_int(Rcpp::List x);


