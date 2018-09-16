
#include "Data_counts.h"
#include "misc_v1.h"

using namespace std;

//------------------------------------------------
// initialise
void Data_counts::init(int K, int L, int COI_max) {
  
  counts = vector<int>(K*L*COI_max*3);
  logprob = vector<double>(K*L*COI_max*3);
  this->K = K;
  this->L = L;
  this->COI_max = COI_max;
  
}

void Data_counts::populate(const vector<vector<int>> &data, const vector<int> &group, const vector<int> &m) {
  
  // empty counts vector
  fill(counts.begin(), counts.end(), 0);
  
  // populate counts vector
  int n = data.size();
  for (int i=0; i<n; ++i) {
    int this_group = group[i];
    int this_m = m[i];
    for (int l=0; l<L; ++l) {
      counts[L*COI_max*3*this_group + COI_max*3*l + 3*(this_m-1) + data[i][l]-1]++;
    }
  }
  
}

int Data_counts::get_counts(int k, int l, int m, int S) {
  return counts[L*COI_max*3*k + COI_max*3*l + 3*(m-1) + S-1];
}

void Data_counts::set_logprob(int k, int l, int m, int S, double x) {
  logprob[L*COI_max*3*k + COI_max*3*l + 3*(m-1) + S-1] = x;
}

double Data_counts::get_logprob(int k, int l, int m, int S) {
  return logprob[L*COI_max*3*k + COI_max*3*l + 3*(m-1) + S-1];
}
