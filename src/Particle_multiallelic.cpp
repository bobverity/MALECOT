
#include "Particle_multiallelic.h"
#include "misc_v1.h"
#include "probability.h"
#include "hungarian.h"

using namespace std;

//------------------------------------------------
// constructor for Particle_multiallelic class
Particle_multiallelic::Particle_multiallelic(double beta_raised) {
  
  // beta_raised stores values of beta (the thermodynamic power), raised to the
  // power GTI_pow
  this->beta_raised = beta_raised;
  
  // scalar lambda value
  lambda0 = lambda[0][0];
  
  // initialise proposal standard deviations
  p_propSD = vector<vector<double>>(K, vector<double>(L,1));
  
  // initialise COI_mean
  COI_mean_vec = vector<double>(K, COI_mean);
  COI_mean_shape = vector<double>(K);
  COI_mean_rate = vector<double>(K);
  COI_mean_propSD = vector<double>(K,1);
  
  // grouping
  group = vector<int>(n);
  
  // likelihood
  loglike_old = vector<vector<double>>(n, vector<double>(L));
  loglike_new = vector<vector<double>>(n, vector<double>(L));
  loglike_new_group = vector<vector<double>>(K, vector<double>(L));
  loglike = 0;
  
  // COI and allele frequencies
  m = vector<int>(n,COI_max);
  p = vector<vector<vector<double>>>(K, vector<vector<double>>(L));
  logp = p;
  for (int k=0; k<K; k++) {
    for (int j=0; j<L; j++) {
      p[k][j] = vector<double>(alleles[j], 1/double(alleles[j]));
      logp[k][j] = vector<double>(alleles[j], -log(alleles[j]));
    }
  }
  
  // qmatrices
  log_qmatrix = vector<vector<double>>(n, vector<double>(K));
  qmatrix = vector<vector<double>>(n, vector<double>(K));
  
  // objects used when updating draws
  v = vector<double>(max(observed_COI));
  sum_loglike_old_vec = vector<double>(K);
  sum_loglike_new_vec = vector<double>(K);
  p_prop = p;
  logp_prop = logp;
  COI_mean_prop = vector<double>(K);
  
  // initialise ordering of labels
  label_order = seq_int(0,K-1);
  label_order_new = vector<int>(K);
  
  // objects for solving label switching problem
  cost_mat = vector<vector<double>>(K, vector<double>(K));
  best_perm = vector<int>(K);
  best_perm_order = vector<int>(K);
  edges_left = vector<int>(K);
  edges_right = vector<int>(K);
  blocked_left = vector<int>(K);
  blocked_right = vector<int>(K);
  
  // store acceptance rates
  p_accept = vector<vector<int>>(K, vector<int>(L));
  
}

//------------------------------------------------
// calculate exact log-probabilty of genotypes. Credit goes to Inna Gerlovina
// for the clever way of looping over configurations
double Particle_multiallelic::logprob_genotype(const vector<int> &x, const vector<double> &logp, int m) {
  
  // skip over missing data
  if (x[0] == 0) {
    return 1;
  }
  
  // special case if single allele observed at this locus - only one configuration possible
  int x_size = x.size();
  if (x_size == 1) {
    return m*logp[x[0]-1];
  }
  
  // create objects for looping through configurations
  int vs = x_size - 1;
  int vmax = m - vs;
  int no = vmax;
  fill(v.begin(), v.end(), 1);
  
  // likelihood of first configuration
  double const1 = lookup_lgamma[m+1];
  double sp = const1 + no*logp[x[vs]-1] - lookup_lgamma[no+1];
  for (int j=0; j<vs; ++j) {
    sp += logp[x[j]-1];
  }
  double likelihood = exp(sp);
  
  // loop through all remaining configurations
  int i = 0;
  while (v[vs-1] < vmax) {
    
    // update configuration
    if (v[0] < vmax && no > 1) {
      v[0]++;
      no--;
      i = 0;
    } else {
      if (no <= 1) {
        i++;
      }
      while (v[i] == vmax && i < vs) {
        i++;
      }
      v[i]++;
      no = m - i;
      for (int k=i; k<vs; k++) {
        no -= v[k];
      }
      if (no < 1) {
        continue;
      }
      for (int k=0; k<i; k++) {
        v[k] = 1;
      }
    }
    
    // update likelihood
    sp = const1 + no*logp[x[vs]-1] - lookup_lgamma[no+1];
    for (int j=0; j<vs; ++j) {
      sp += v[j]*logp[x[j]-1] - lookup_lgamma[v[j]+1];
    }
    likelihood += exp(sp);
    
  }   // end while loop over configurations
  
  // calculate and store log-likelihood
  double ret = log(likelihood);
  ret = (ret < -OVERFLO) ? -OVERFLO : ret;    // catch -inf values
  
  // return log-likelihood
  return ret;
}

//------------------------------------------------
// reset particle
void Particle_multiallelic::reset() {
  
  // reset qmatrices
  for (int i=0; i<n; i++) {
    fill(qmatrix[i].begin(), qmatrix[i].end(), 0);
    fill(log_qmatrix[i].begin(), log_qmatrix[i].end(), 0);
  }
  
  // reset allele frequencies
  for (int k=0; k<K; k++) {
    for (int j=0; j<L; j++) {
      fill(p[k][j].begin(), p[k][j].end(), 1.0/double(alleles[j]));
    }
    fill(p_propSD[k].begin(), p_propSD[k].end(), 1);
  }
  
  // reset COI
  fill(m.begin(), m.end(), COI_max);
  for (int i=0; i<n; i++) {
    if (COI_manual[i] != -1) {
      m[i] = COI_manual[i];
    }
  }
  
  // reset COI_mean
  fill(COI_mean_vec.begin(), COI_mean_vec.end(), COI_mean);
  
  // reset order of labels
  label_order = seq_int(0,K-1);
  
  // initialise with random grouping
  for (int i=0; i<n; i++) {
    group[i] = sample2(0,K-1);
  }
  
  // calculate initial likelihood
  loglike = 0;
  for (int i=0; i<n; i++) {
    int this_group = group[i];
    int this_m = m[i];
    for (int l=0; l<L; l++) {
      loglike_old[i][l] = logprob_genotype(data[i][l], logp[this_group][l], this_m);
      loglike += loglike_old[i][l];
    }
  }
  
}

//------------------------------------------------
// update allele frequencies p
void Particle_multiallelic::update_p(bool robbins_monro_on, int iteration) {
  
  // loop through loci
  for (int l=0; l<L; l++) {
    
    // clear vectors storing sum over old and new likelihood
    fill(sum_loglike_old_vec.begin(), sum_loglike_old_vec.end(), 0);
    fill(sum_loglike_new_vec.begin(), sum_loglike_new_vec.end(), 0);
    
    // draw p for all demes
    for (int k=0; k<K; k++) {
      p_prop[k][l] = rmlogitnorm2(p[k][l], p_propSD[k][l]);
      for (int j=0; j<alleles[l]; ++j) {
        logp_prop[k][l][j] = log(p_prop[k][l][j]);
      }
    }
    
    // calculate likelihood under proposed p
    for (int i=0; i<n; i++) {
      int this_group = group[i];
      loglike_new[i][l] = logprob_genotype(data[i][l], logp_prop[this_group][l], m[i]);
      sum_loglike_old_vec[this_group] += loglike_old[i][l];
      sum_loglike_new_vec[this_group] += loglike_new[i][l];
    }
    
    // Metropolis step for all K
    for (int k=0; k<K; k++) {
      
      // catch impossible proposed values
      if (sum_loglike_new_vec[k] <= -OVERFLO) {
        continue;
      }
      
      // incorporate prior
      double log_prior_old = dsym_dirichlet1(p[k][l], lambda0);
      double log_prior_new = dsym_dirichlet1(p_prop[k][l], lambda0);
      
      // adjust for proposal density
      double log_prop_forwards = dmlogitnorm2(p_prop[k][l], p[k][l], p_propSD[k][l]);
      double log_prop_backwards = dmlogitnorm2(p[k][l], p_prop[k][l], p_propSD[k][l]);
      
      // Metropolis step
      double MH = (beta_raised*sum_loglike_new_vec[k] + log_prior_new - log_prop_forwards) - (beta_raised*sum_loglike_old_vec[k] + log_prior_old - log_prop_backwards);
      if (log(runif_0_1())<MH) {
        
        // update p
        p[k][l] = p_prop[k][l];
        logp[k][l] = logp_prop[k][l];
        
        // update loglike in all individuals from this group
        for (int i=0; i<n; i++) {
          if (group[i] == k) {
            loglike_old[i][l] = loglike_new[i][l];
          }
        }
        
        // Robbins-Monro positive update
        if (robbins_monro_on) {
          p_propSD[k][l] = exp(log(p_propSD[k][l]) + (1-0.23)/sqrt(iteration));
        }
        
        // update acceptance rates
        p_accept[k][l]++;
        
      } else {
        
        // Robbins-Monro negative update
        if (robbins_monro_on) {
          p_propSD[k][l] = exp(log(p_propSD[k][l]) - 0.23/sqrt(iteration));
        }
        
      }  // end Metropolis-Hastings step
    }  // end loop through demes
  }  // end loop through loci
  
}

//------------------------------------------------
// linear update of m
void Particle_multiallelic::update_m() {
  
  // define sum over old and new likelihood
  double sum_loglike_old = 0;
  double sum_loglike_new = 0;
  
  // loop through individuals
  for (int i=0; i<n; i++) {
    int this_group = group[i];
    
    // skip update if defined manually
    if (COI_manual[i] != -1) {
      continue;
    }
    
    // propose new m
    int m_prop = rbernoulli1(0.5);
    m_prop = (m_prop == 0) ? m[i]-1 : m[i]+1;
    
    // skip if outside range
    if (m_prop < observed_COI[i] || m_prop > COI_max) {
      continue;
    }
    
    // reset sum over old and new likelihood
    sum_loglike_old = 0;
    sum_loglike_new = 0;
    
    // calculate likelihood
    for (int l=0; l<L; l++) {
      loglike_new[i][l] = logprob_genotype(data[i][l], logp[this_group][l], m_prop);
      sum_loglike_old += beta_raised*loglike_old[i][l];
      sum_loglike_new += beta_raised*loglike_new[i][l];
    }
    
    // catch impossible proposed values
    if (sum_loglike_new <= -OVERFLO) {
      continue;
    }
    
    // apply Poisson or negative Binomial prior
    if (COI_model == 2) {
      sum_loglike_old += dpois1(m[i]-1, COI_mean_vec[this_group]-1);
      sum_loglike_new += dpois1(m_prop-1, COI_mean_vec[this_group]-1);
    } else if (COI_model == 3) {
      sum_loglike_old += dnbinom1(m[i]-1, COI_mean_vec[this_group]-1, COI_dispersion);
      sum_loglike_new += dnbinom1(m_prop-1, COI_mean_vec[this_group]-1, COI_dispersion);
    }
    
    // Metropolis step
    if (log(runif_0_1()) < (sum_loglike_new - sum_loglike_old)) {
      m[i] = m_prop;
      for (int l=0; l<L; l++) {
        loglike_old[i][l] = loglike_new[i][l];
      }
    }
    
  }   // end loop through individuals
  
}

//------------------------------------------------
// update group allocation
void Particle_multiallelic::update_group() {
  
  // no change if K==1
  if (K == 1) {
    return;
  }
  
  // loop through individuals
  for (int i=0; i<n; i++) {
    int this_group = group[i];
    
    // reset log_qmatrix for this individual
    fill(log_qmatrix[i].begin(), log_qmatrix[i].end(), 0);
    
    // calculate likelihood
    for (int k=0; k<K; k++) {
      fill(loglike_new_group[k].begin(), loglike_new_group[k].end(), 0);
      if (k == this_group) {  // likelihood already calculated for current group
        for (int j=0; j<L; j++) {
          loglike_new_group[k][j] = loglike_old[i][j];
          if (loglike_new_group[k][j] > -OVERFLO) {
            log_qmatrix[i][k] += beta_raised*loglike_new_group[k][j];
          } else {
            log_qmatrix[i][k] = -OVERFLO;
          }
        }
      } else {  // calculate likelihood for group k
        for (int j=0; j<L; j++) {
          loglike_new_group[k][j] = logprob_genotype(data[i][j], logp[k][j], m[i]);
          if (loglike_new_group[k][j] > -OVERFLO) {
            log_qmatrix[i][k] += beta_raised*loglike_new_group[k][j];
          } else {
            log_qmatrix[i][k] = -OVERFLO;
          }
        }
      }
      
      // add prior information
      if (COI_model == 2) {
        log_qmatrix[i][k] += dpois1(m[i]-1, COI_mean_vec[k]-1);
      } else if (COI_model == 3) {
        log_qmatrix[i][k] += dnbinom1(m[i]-1, COI_mean_vec[k]-1, COI_dispersion);
      }
      
      // limit small (-inf) values
      if (log_qmatrix[i][k] <= -OVERFLO) {
        log_qmatrix[i][k] = -OVERFLO;
      }
    }
    
    // exponentiate without underflow
    double log_qmatrix_max = max(log_qmatrix[i]);
    double qmatrix_sum = 0;
    for (int k=0; k<K; k++) {
      log_qmatrix[i][k] -= log_qmatrix_max;
      qmatrix[i][k] = exp(log_qmatrix[i][k]);
      qmatrix_sum += qmatrix[i][k];
    }
    for (int k=0; k<K; k++) {
      qmatrix[i][k] /= qmatrix_sum;
    }
    
    // sample new group
    int new_group = sample1(qmatrix[i], 1.0)-1;
    
    // if group has changed update likelihood
    if (new_group!=this_group) {
      for (int j=0; j<L; j++) {
        loglike_old[i][j] = loglike_new_group[new_group][j];
      }
    }
    
    // update group
    group[i] = new_group;
  }
  
}

//------------------------------------------------
// update mean COI
void Particle_multiallelic::update_COI_mean(bool robbins_monro_on, int iteration) {
  
  // split method between poisson and negative binomial
  // poisson model
  if (COI_model == 2) {
    
    // reset shape and rate vectors
    fill(COI_mean_shape.begin(), COI_mean_shape.end(), 0.25);
    fill(COI_mean_rate.begin(), COI_mean_rate.end(), 0.25);
    
    // update shape and rate vectors to posterior values
    for (int i=0; i<n; i++) {
      int this_group = group[i];
      COI_mean_shape[this_group] += m[i]-1;
      COI_mean_rate[this_group] += 1;
    }
    
    // draw new COI means from conditional posterior
    for (int k=0; k<K; k++) {
      COI_mean_vec[k] = rgamma1(COI_mean_shape[k], COI_mean_rate[k]) + 1;
    }
  } // end poisson model
  
  // negative binomial model
  if (COI_model == 3) {
    
    // clear vectors storing sum over old and new likelihood
    fill(sum_loglike_old_vec.begin(), sum_loglike_old_vec.end(), 0);
    fill(sum_loglike_new_vec.begin(), sum_loglike_new_vec.end(), 0);
    
    // draw COI_mean for all demes
    for (int k=0; k<K; k++) {
      COI_mean_prop[k] = rnorm1_interval(COI_mean_vec[k], COI_mean_propSD[k], 1, COI_max);
    }
    
    // calculate likelihood
    for (int i=0; i<n; i++) {
      sum_loglike_old_vec[group[i]] += dnbinom1(m[i]-1, COI_mean_vec[group[i]]-1, COI_dispersion);
      sum_loglike_new_vec[group[i]] += dnbinom1(m[i]-1, COI_mean_prop[group[i]]-1, COI_dispersion);
    }
    
    // Metropolis step for all demes
    for (int k=0; k<K; k++) {
      
      // Metropolis step
      if (log(runif_0_1())<(sum_loglike_new_vec[k] - sum_loglike_old_vec[k])) {
        COI_mean_vec[k] = COI_mean_prop[k];
        
        // Robbins-Monro positive update
        if (robbins_monro_on) {
          COI_mean_propSD[k] += (1-0.23)/sqrt(double(iteration));
        }
        
      } else {
      
      // Robbins-Monro negative update
      if (robbins_monro_on) {
        COI_mean_propSD[k] -= 0.23/sqrt(double(iteration));
        if (COI_mean_propSD[k] < UNDERFLO) {
          COI_mean_propSD[k] = UNDERFLO;
        }
      }
      
      }
    }  // end Metropolis step
  } // end negative binomial model
  
}

//------------------------------------------------
// solve label switching problem
void Particle_multiallelic::solve_label_switching(const vector<vector<double>> &log_qmatrix_running) {
  
  for (int k1=0; k1<K; k1++) {
    fill(cost_mat[k1].begin(), cost_mat[k1].end(), 0);
    for (int k2=0; k2<K; k2++) {
      for (int i=0; i<n; i++) {
        cost_mat[k1][k2] += qmatrix[i][label_order[k1]]*(log_qmatrix[i][label_order[k1]] - log_qmatrix_running[i][k2]);
      }
    }
  }
  
  // find best permutation of current labels using Hungarian algorithm
  best_perm = hungarian(cost_mat, edges_left, edges_right, blocked_left, blocked_right);

  // define best_perm_order
  for (int k=0; k<K; k++) {
    best_perm_order[best_perm[k]] = k;
  }
  
  // replace old label order with new
  for (int k=0; k<K; k++) {
    label_order_new[k] = label_order[best_perm_order[k]];
  }
  label_order = label_order_new;
  
}

//------------------------------------------------
// calculate overall log-likelihood
void Particle_multiallelic::calculate_loglike() {
  loglike = 0;
  for (int i=0; i<n; i++) {
    for (int j=0; j<L; j++) {
      loglike += loglike_old[i][j];
    }
  }
}

