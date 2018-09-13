
#include "Lookup.h"
#include "misc_v1.h"

using namespace std;

//------------------------------------------------
// declare static member variables for Lookup class

vector<vector<double>> Lookup::lookup_homo;
vector<vector<double>> Lookup::lookup_het;
vector<double> Lookup::lookup_lgamma;

//------------------------------------------------
// initialise homo/het lookup tables
void Lookup::init_homohet(){
  
  // no need for lookup if using zero precision
  if (precision == 0) {
    return;
  }
    
  // initialise tables
  lookup_homo = vector<vector<double>>(precision_size + 1, vector<double>(COI_max));
  lookup_het = vector<vector<double>>(precision_size + 1, vector<double>(COI_max));
  
  // populate tables
  for (int i=0; i<(precision_size + 1); i++) {
    
    // this allele frequency
    double p = i/double(precision_size);
    
    // loop through COI
    for (int m=0; m<COI_max; m++) {
      
      // if not estimating error terms then incorporate error at this stage, and
      // store log probabilities. Otherwise store raw probability and error will
      // be incorporated later
      if (!estimate_error) {
        
        // homo lookup
        lookup_homo[i][m] = log( (1.0-e1)*((double)pow(p,m+1)) + 0.5*e2*(1.0-(double)pow(p,m+1)-(double)pow(1.0-p,m+1)) );
        lookup_homo[i][m] = (lookup_homo[i][m]<(-OVERFLO)) ? -OVERFLO : lookup_homo[i][m];
        
        // het lookup
        lookup_het[i][m] = log( e1*((double)pow(p,m+1)) + e1*((double)pow(1.0-p,m+1)) + (1.0-e2)*(1.0-pow(p,m+1)-pow(1.0-p,m+1)) );
        lookup_het[i][m] = (lookup_het[i][m]<(-OVERFLO)) ? -OVERFLO : lookup_het[i][m];
        
      } else {
        // homo lookup
        lookup_homo[i][m] = (double)pow(p,m+1);
        
        // het lookup
        lookup_het[i][m] = (1.0-(double)pow(p,m+1)-(double)pow(1.0-p,m+1));
      }
    }
  }
  
}

//------------------------------------------------
// initialise lgamma lookup tables
void Lookup::init_lgamma() {
  
  // populate lgamma table
  lookup_lgamma = vector<double>(COI_max+2);
  for (int i=0; i<int(lookup_lgamma.size()); ++i) {
    lookup_lgamma[i] = lgamma(i);
  }
  
}
