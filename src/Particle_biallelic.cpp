
#include "Particle_biallelic.h"
#include "misc_v1.h"
#include "probability.h"
#include "hungarian.h"

using namespace std;

//------------------------------------------------
// constructor for Particle_biallelic class
Particle_biallelic::Particle_biallelic(double beta_raised) {
  
  // beta_raised stores values of beta (the thermodynamic power), raised to the
  // power GTI_pow
  this->beta_raised = beta_raised;
  
  // initialise data frequencies
  data_counts.init(K, L, COI_max);
  
  // initialise proposal standard deviations
  e_propSD = 1;
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
  p = vector<vector<double>>(K, vector<double>(L,0.5));
  
  // qmatrices
  log_qmatrix = vector<vector<double>>(n, vector<double>(K));
  qmatrix = vector<vector<double>>(n, vector<double>(K));
  
  // probability vectors and matrices used when updating draws
  sum_loglike_old_vec = vector<double>(K);
  sum_loglike_new_vec = vector<double>(K);
  p_prop = 0;
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
  e_accept = 0;
  
  // vary core likelihood function depending on user choices
  if (precision == 0) {
    logprob_genotype_ptr = &Particle_biallelic::logprob_genotype_exact;
  } else {
    if (estimate_error) {
      logprob_genotype_ptr = &Particle_biallelic::logprob_genotype_lookup_varE;
    } else {
      logprob_genotype_ptr = &Particle_biallelic::logprob_genotype_lookup;
    }
  }
  
}

//########################################################################################################
// genotype probabilities

//------------------------------------------------
// switch between functions for calculating log-probability of genotypes
double Particle_biallelic::logprob_genotype(int S, double p, int m, double e1, double e2) {
  return (this->*logprob_genotype_ptr)(S, p, m, e1, e2);
}

//------------------------------------------------
// exact log-probabilty of genotypes
double Particle_biallelic::logprob_genotype_exact(int S, double p, int m, double e1, double e2) {
  double ret = 0; // return value for missing data
  if (S == 1) {
    ret = log( (1.0-e1)*((double)pow(p,m)) + 0.5*e2*(1.0-(double)pow(p,m)-(double)pow(1.0-p,m)) );
  } else if (S == 3) {
    ret = log( (1.0-e1)*((double)pow(1.0-p,m)) + 0.5*e2*(1.0-(double)pow(p,m)-(double)pow(1.0-p,m)) );
  } else if (S == 2) {
    ret = log( e1*((double)pow(p,m)) + e1*((double)pow(1.0-p,m)) + (1.0-e2)*(1.0-pow(p,m)-pow(1.0-p,m)) );
  }
  ret = (ret < -OVERFLO) ? -OVERFLO : ret;    // catch -inf values
  return ret;
}

//------------------------------------------------
// log-probabilty of genotypes using lookup tables
double Particle_biallelic::logprob_genotype_lookup(int S, double p, int m, double e1, double e2) {
  int p_index = round(p*precision_size);
  double ret = 0; // return value for missing data
  if (S == 1) {
    ret = lookup_homo[p_index][m-1];
  } else if (S == 3) {
    ret = lookup_homo[precision_size-p_index][m-1];
  } else if (S == 2) {
    ret = lookup_het[p_index][m-1];
  }
  ret = (ret < -OVERFLO) ? -OVERFLO : ret;    // catch -inf values
  return ret;
}

//------------------------------------------------
// log-probabilty of genotypes using lookup tables, with error terms specified
double Particle_biallelic::logprob_genotype_lookup_varE(int S, double p, int m, double e1, double e2) {
  int p_index = (int)(p*precision_size);
  double ret = 0; // return value for missing data
  if (S == 1) {
    ret = log( (1.0-e1)*lookup_homo[p_index][m-1] + 0.5*e2*lookup_het[p_index][m-1]);
  } else if (S == 3) {
    ret = log( (1.0-e1)*lookup_homo[precision_size-p_index][m-1] + 0.5*e2*lookup_het[p_index][m-1]);
  } else if (S == 2) {
    ret = log( e1*lookup_homo[p_index][m-1] + e1*lookup_homo[precision_size-p_index][m-1] + (1.0-e2)*lookup_het[p_index][m-1] );
  }
  ret = (ret < -OVERFLO) ? -OVERFLO : ret;    // catch -inf values
  return ret;
}

//########################################################################################################

//------------------------------------------------
// reset particle
void Particle_biallelic::reset() {
  
  // reset qmatrices
  for (int i=0; i<n; i++) {
    fill(qmatrix[i].begin(), qmatrix[i].end(), 0);
    fill(log_qmatrix[i].begin(), log_qmatrix[i].end(), 0);
  }
  
  // reset allele frequencies
  for (int k=0; k<K; k++) {
    fill(p[k].begin(), p[k].end(), 0.5);
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
    for (int j=0; j<L; j++) {
      loglike_old[i][j] = logprob_genotype(data[i][j], p[group[i]][j], m[i], e1, e2);
      loglike += loglike_old[i][j];
    }
  }
  
}

//------------------------------------------------
// get lambda value for locus i, allele j
double Particle_biallelic::get_lambda(int i, int j) {
  if (lambda_scalar) {
    return lambda[0][0];
  } else {
    return lambda[i][j];
  }
  return 1.0;
}

//------------------------------------------------
// update allele frequencies p
void Particle_biallelic::update_p(bool robbins_monro_on, int iteration) {
  /*
  // empty counts table
  for (int k=0; k<K; ++k) {
    for (int l=0; l<L; ++l) {
      for (int i=0; i<3; ++i) {
        fill(data_counts[k][l][i].begin(), data_counts[k][l][i].end(), 0);
      }
    }
  }
  
  // populate counts table and calculate current likelihood
  for (int i=0; i<n; i++) {
    int this_group = group[i];
    for (int l=0; l<L; ++l) {
      data_counts[this_group][l][data[i][l]-1][m[i]-1]++;
    }
  }
  */
  
  // clear and repopulate data_counts
  data_counts.populate(data, group, m);
  
  vector<double> p_prop(K);
  
  // loop through loci
  for (int l=0; l<L; l++) {
    
    // clear vectors storing sum over old and new likelihood
    fill(sum_loglike_old_vec.begin(), sum_loglike_old_vec.end(), 0);
    fill(sum_loglike_new_vec.begin(), sum_loglike_new_vec.end(), 0);
    
    for (int k=0; k<K; ++k) {
      p_prop[k] = rnorm1_interval(p[k][l], p_propSD[k][l], 0, 1.0);
    }
    
    /*
    // calculate current likelihood
    for (int i=0; i<n; i++) {
      sum_loglike_old_vec[group[i]] += loglike_old[i][l];
    }
    */
    // calculate likelihood
    for (int i=0; i<n; i++) {
      int this_group = group[i];
      loglike_new[i][l] = logprob_genotype(data[i][l], p_prop[this_group], m[i], e1, e2);
      sum_loglike_old_vec[this_group] += loglike_old[i][l];
      sum_loglike_new_vec[this_group] += loglike_new[i][l];
    }
    
    // Metropolis step for all K
    for (int k=0; k<K; k++) {
      
      // draw p
      //p_prop = rnorm1_interval(p[k][l], p_propSD[k][l], 0, 1.0);
      
      /*
      // calculate and store logprobs from counts
      double tmp = 0;
      for (int m=0; m<COI_max; ++m) {
        int count1 = data_counts.get_counts(k, l, m+1, 1);
        if (count1 > 0) {
          tmp = logprob_genotype(1, p_prop[k], m+1, e1, e2);
          data_counts.set_logprob(k, l, m+1, 1, tmp);
          sum_loglike_new_vec[k] += count1 * tmp;
        }
        int count2 = data_counts.get_counts(k, l, m+1, 2);
        if (count2 > 0) {
          tmp = logprob_genotype(2, p_prop[k], m+1, e1, e2);
          data_counts.set_logprob(k, l, m+1, 2, tmp);
          sum_loglike_new_vec[k] += count2 * tmp;
        }
        int count3 = data_counts.get_counts(k, l, m+1, 3);
        if (count3 > 0) {
          tmp = logprob_genotype(3, p_prop[k], m+1, e1, e2);
          data_counts.set_logprob(k, l, m+1, 3, tmp);
          sum_loglike_new_vec[k] += count3 * tmp;
        }
      }
      */
      // catch impossible proposed values
      if (sum_loglike_new_vec[k] <= -OVERFLO) {
        continue;
      }
      
      // incorporate prior
      double log_prior_old = dbeta1(p[k][l], get_lambda(l,0), get_lambda(l,1));
      //double log_prior_new = dbeta1(p_prop, get_lambda(l,0), get_lambda(l,1));
      double log_prior_new = dbeta1(p_prop[k], get_lambda(l,0), get_lambda(l,1));
      
      // Metropolis step
      double MH = (beta_raised*sum_loglike_new_vec[k] + log_prior_new) - (beta_raised*sum_loglike_old_vec[k] + log_prior_old);
      if (log(runif_0_1())<MH) {
        
        // update p
        //p[k][l] = p_prop;
        p[k][l] = p_prop[k];
        
        // update loglike in all individuals from this group
        for (int i=0; i<n; i++) {
          if (group[i] == k) {
            loglike_old[i][l] = loglike_new[i][l];
            //loglike_old[i][l] = data_counts.get_logprob(k, l, m[i], data[i][l]);
          }
        }
        
        // Robbins-Monro positive update
        if (robbins_monro_on) {
          p_propSD[k][l] = exp(log(p_propSD[k][l]) + (1-0.23)/sqrt(iteration));
        } else {
          p_accept[k][l]++;
        }
        
      } else {
        
        // Robbins-Monro negative update
        if (robbins_monro_on) {
          p_propSD[k][l] = exp(log(p_propSD[k][l]) - 0.23/sqrt(iteration));
        }
        
      }  // end Metropolis step
      
      // limit p_propSD[k][j]
      p_propSD[k][l] = (p_propSD[k][l] > 1) ? 1 : p_propSD[k][l];
      
    }  // end loop through k
  }   // end loop through loci
  
}

//------------------------------------------------
// update error parameters e1 and e2
void Particle_biallelic::update_e(bool robbins_monro_on, int iteration) {
  
  // define sum over old and new likelihood
  double sum_loglike_old = 0;
  double sum_loglike_new = 0;
  
  // propose new value
  double e1_prop = rnorm1_interval(e1, e_propSD, 0, e1_max);
  double e2_prop = rnorm1_interval(e2, e_propSD, 0, e2_max);
  
  // calculate likelihood
  for (int i=0; i<n; i++) {
    for (int j=0; j<L; j++) {
      //loglike_new[i][j] = logprob_genotype(data[i][j], p[group[i]][j], m[i], e1_prop, e2_prop);
      sum_loglike_old += loglike_old[i][j];
      //sum_loglike_new += loglike_new[i][j];
    }
  }
  ///*
  // loop through loci
  for (int l=0; l<L; ++l) {
    
    // add to current likelihood
    //for (int i=0; i<n; i++) {
    //  sum_loglike_old += loglike_old[i][l];
    //}
    
    // calculate new likelihood
    for (int k=0; k<K; ++k) {
      
      double tmp = 0;
      for (int m=0; m<COI_max; ++m) {
        int count1 = data_counts.get_counts(k, l, m+1, 1);
        if (count1 > 0) {
          tmp = logprob_genotype(1, p[k][l], m+1, e1_prop, e2_prop);
          data_counts.set_logprob(k, l, m+1, 1, tmp);
          sum_loglike_new += count1 * tmp;
        }
        int count2 = data_counts.get_counts(k, l, m+1, 2);
        if (count2 > 0) {
          tmp = logprob_genotype(2, p[k][l], m+1, e1_prop, e2_prop);
          data_counts.set_logprob(k, l, m+1, 2, tmp);
          sum_loglike_new += count2 * tmp;
        }
        int count3 = data_counts.get_counts(k, l, m+1, 3);
        if (count3 > 0) {
          tmp = logprob_genotype(3, p[k][l], m+1, e1_prop, e2_prop);
          data_counts.set_logprob(k, l, m+1, 3, tmp);
          sum_loglike_new += count3 * tmp;
        }
      }
      
    }  // end loop over K
    
  }  // end loop over L
  //*/
  // catch impossible proposed values
  if (sum_loglike_new <= -OVERFLO) {
    return;
  }
  
  // Metropolis step
  if (log(runif_0_1()) < beta_raised*(sum_loglike_new - sum_loglike_old)) {
    
    // update e1 and e2
    e1 = e1_prop;
    e2 = e2_prop;
    
    // update loglike
    for (int i=0; i<n; i++) {
      for (int l=0; l<L; l++) {
        loglike_old[i][l] = data_counts.get_logprob(group[i], l, m[i], data[i][l]);
        //loglike_old[i][l] = loglike_new[i][l];
      }
    }
    
    // Robbins-Monro positive update
    if (robbins_monro_on) {
      e_propSD = exp(log(e_propSD) + (1-0.23)/sqrt(iteration));
    } else {
      e_accept++;
    }
    
  } else {
    
    // Robbins-Monro negative update step
    if (robbins_monro_on) {
      e_propSD = exp(log(e_propSD) - 0.23/sqrt(iteration));
    }
    
  } // end Metropolis step
  
  // limit e_propSD
  e_propSD = (e_propSD > 1) ? 1 : e_propSD;
  
}

//------------------------------------------------
// linear update of m
void Particle_biallelic::update_m() {
  
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
    m_prop = (m_prop==0) ? m[i]-1 : m[i]+1;
    
    // skip if outside range
    if (m_prop == 0 || m_prop > COI_max) {
      continue;
    }
    
    // reset sum over old and new likelihood
    sum_loglike_old = 0;
    sum_loglike_new = 0;
    
    // calculate likelihood
    for (int j=0; j<L; j++) {
      loglike_new[i][j] = logprob_genotype(data[i][j], p[this_group][j], m_prop, e1, e2);
      sum_loglike_old += beta_raised*loglike_old[i][j];
      sum_loglike_new += beta_raised*loglike_new[i][j];
    }
    /*
    // calculate likelihood
    for (int l=0; l<L; ++l) {
      if (data_counts[this_group][l][data[i][l]-1][m_prop-1] == 0) {
        data_counts_logprob[this_group][l][data[i][l]-1][m_prop-1] = logprob_genotype(data[i][l], p[this_group][l], m_prop, e1, e2);
      }
      loglike_new[i][l] = data_counts_logprob[this_group][l][data[i][l]-1][m_prop-1];
      sum_loglike_old += loglike_old[i][l];
      sum_loglike_new += loglike_new[i][l];
    }
    */
    // raise to beta power
    sum_loglike_old *= beta_raised;
    sum_loglike_new *= beta_raised;
    
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
      
      // update m
      m[i] = m_prop;
      
      // update loglike
      for (int j=0; j<L; j++) {
        loglike_old[i][j] = loglike_new[i][j];
      }
      /*
      // update counts
      for (int l=0; l<L; ++l) {
        data_counts[this_group][l][data[i][l]-1][m[i]-1]--;
        data_counts[this_group][l][data[i][l]-1][m_prop-1]++;
      }
      */
    }  // end Metropolis-Hastings
    
  }   // end loop through individuals
  
}

//------------------------------------------------
// update group allocation
void Particle_biallelic::update_group() {
  
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
          //if (data_counts[k][j][data[i][j]-1][m[i]-1] == 0) {
          //  data_counts_logprob[k][j][data[i][j]-1][m[i]-1] = logprob_genotype(data[i][j], p[k][j], m[i], e1, e2);
          //}
          //loglike_new_group[k][j] = data_counts_logprob[k][j][data[i][j]-1][m[i]-1];
          loglike_new_group[k][j] = logprob_genotype(data[i][j], p[k][j], m[i], e1, e2);
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
    
    // if group has changed
    if (new_group!=this_group) {
      
      //  update likelihood
      for (int j=0; j<L; j++) {
        loglike_old[i][j] = loglike_new_group[new_group][j];
      }
      /*
      // update counts
      for (int l=0; l<L; ++l) {
        data_counts[this_group][l][data[i][l]-1][m[i]-1]--;
        data_counts[new_group][l][data[i][l]-1][m[i]-1]++;
      }
      */
    }
    
    // update group
    group[i] = new_group;
    
  }  // end loop through individuals
  
}



//------------------------------------------------
// update mean COI
void Particle_biallelic::update_COI_mean(bool robbins_monro_on, int iteration) {
  
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
void Particle_biallelic::solve_label_switching(const vector<vector<double>> &log_qmatrix_running) {
  
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
void Particle_biallelic::calculate_loglike() {
  loglike = 0;
  for (int i=0; i<n; i++) {
    for (int j=0; j<L; j++) {
      loglike += loglike_old[i][j];
    }
  }
}
