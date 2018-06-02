
#include "particle_biallelic.h"
#include "misc.h"
#include "probability.h"
#include "Hungarian.h"

using namespace std;

//------------------------------------------------
// default constructor for particle_biallelic class
particle_biallelic::particle_biallelic() {}

//------------------------------------------------
// constructor for particle_biallelic class
particle_biallelic::particle_biallelic(vector<vector<int>> &data, Rcpp::List &args, vector<vector<double>> &lookup_homo, vector<vector<double>> &lookup_het, double beta_raised_) {
  
  // pointers to observed data and basic data properties
  data_ptr = &data;
  n = rcpp_to_int(args["n"]);
  L = rcpp_to_int(args["L"]);
  
  // pointers to lookup tables
  lookup_homo_ptr = &lookup_homo;
  lookup_het_ptr = &lookup_het;
  
  // model and MCMC parameters
  K = rcpp_to_int(args["K"]);
  precision = rcpp_to_double(args["precision"]);
  precision_size = 1/precision;
  estimate_error = rcpp_to_bool(args["estimate_error"]);
  COI_model = rcpp_to_int(args["COI_model_numeric"]);
  COI_max = rcpp_to_int(args["COI_max"]);
  COI_dispersion = rcpp_to_double(args["COI_dispersion"]);
  beta_raised = beta_raised_;
  scaffold_n = rcpp_to_int(args["scaffold_n"]);
  e1 = rcpp_to_double(args["e1"]);
  e2 = rcpp_to_double(args["e2"]);
  e1_max = rcpp_to_double(args["e1_max"]);
  e2_max = rcpp_to_double(args["e2_max"]);
  
  // initialise proposal standard deviations
  e1_propSD = 1;
  e2_propSD = 1;
  p_propSD = vector<vector<double>>(K, vector<double>(L,1));

  // initialise COI_mean
  COI_mean = vector<double>(K, COI_max);
  COI_mean_shape = vector<double>(K);
  COI_mean_rate = vector<double>(K);
  COI_mean_propSD = vector<double>(K,1);
  
  // grouping
  group = vector<int>(n);
  group_order = seq_int(0,K-1);
  
  // likelihood
  loglike_old = vector<vector<double>>(n, vector<double>(L));
  loglike_new = vector<vector<double>>(n, vector<double>(L));
  loglike_new_group = vector<vector<double>>(K, vector<double>(L));
  loglike = 0;
  
  // COI and allele frequencies
  m = vector<int>(n,COI_max);
  p = vector<vector<double>>(K, vector<double>(L,0.5));
  
  // probability vectors and matrices used when updating draws
  log_qmatrix = vector<vector<double>>(n, vector<double>(K));
  qmatrix = vector<vector<double>>(n, vector<double>(K));
  
  sum_loglike_old_vec = vector<double>(K);
  sum_loglike_new_vec = vector<double>(K);
  p_prop = vector<double>(K);
  COI_mean_prop = vector<double>(K);
  
  
  // initialise Q matrix
  //qmatrix = vector<vector<double>>(n, vector<double>(K));
  //log_qmatrix = vector<vector<double>>(n, vector<double>(K));
  
  // initialise ordering of labels
  label_order = seq_int(0,K-1);
  label_order_new = vector<int>(K);
  /*
  // objects for solving label switching problem
  cost_mat = vector<vector<double>>(K, vector<double>(K));
  best_perm = vector<int>(K);
  best_perm_order = vector<int>(K);
  edges_left = vector<int>(K);
  edges_right = vector<int>(K);
  blocked_left = vector<int>(K);
  blocked_right = vector<int>(K);
  */
  // initialise scaffold objects
  scaf_group = vector<vector<int>>(scaffold_n, vector<int>(n));
  scaf_count = vector<int>(scaffold_n);
  /*
  // initialise objects for split-merge
  splitmerge_targets = seq_int(0,K-1);
  splitmerge_mu = vector<double>(K);
  splitmerge_group = vector<int>(n);
  splitmerge_counts = vector<int>(K);
  splitmerge_x_sum = vector<double>(K);
  */
  
  // store acceptance rates
  p_accept = vector<vector<int>>(K, vector<int>(L));
  e1_accept = 0;
  e2_accept = 0;
  
  // vary core likelihood function depending on user choices
  if (precision==0) {
    logprob_genotype_ptr = &particle_biallelic::logprob_genotype_exact;
  } else {
    if (estimate_error) {
      logprob_genotype_ptr = &particle_biallelic::logprob_genotype_lookup_varE;
    } else {
      logprob_genotype_ptr = &particle_biallelic::logprob_genotype_lookup;
    }
  }
  
}

//########################################################################################################
// genotype probabilities

//------------------------------------------------
// switch between functions for calculating log-probability of genotypes
double particle_biallelic::logprob_genotype(int S, double p, int m, double e1, double e2) {
  return (this->*logprob_genotype_ptr)(S, p, m, e1, e2);
}

//------------------------------------------------
// exact log-probabilty of genotypes
double particle_biallelic::logprob_genotype_exact(int S, double p, int m, double e1, double e2) {
  double ret = 0; // return value for missing data
  if (S==1) {
    ret = log( (1.0-e1)*((double)pow(p,m)) + 0.5*e2*(1.0-(double)pow(p,m)-(double)pow(1.0-p,m)) );
  } else if (S==3) {
    ret = log( (1.0-e1)*((double)pow(1.0-p,m)) + 0.5*e2*(1.0-(double)pow(p,m)-(double)pow(1.0-p,m)) );
  } else if (S==2) {
    ret = log( e1*((double)pow(p,m)) + e1*((double)pow(1.0-p,m)) + (1.0-e2)*(1.0-pow(p,m)-pow(1.0-p,m)) );
  }
  ret = (ret < -OVERFLO) ? -OVERFLO : ret;    // catch -inf values
  return ret;
}

//------------------------------------------------
// log-probabilty of genotypes using lookup tables
double particle_biallelic::logprob_genotype_lookup(int S, double p, int m, double e1, double e2) {
  int p_index = round(p*precision_size);
  double ret = 0; // return value for missing data
  if (S==1) {
    ret = (*lookup_homo_ptr)[p_index][m-1];
  } else if (S==3) {
    ret = (*lookup_homo_ptr)[precision_size-p_index][m-1];
  } else if (S==2) {
    ret = (*lookup_het_ptr)[p_index][m-1];
  }
  ret = (ret < -OVERFLO) ? -OVERFLO : ret;    // catch -inf values
  return ret;
}

//------------------------------------------------
// log-probabilty of genotypes using lookup tables, with error terms specified
double particle_biallelic::logprob_genotype_lookup_varE(int S, double p, int m, double e1, double e2) {
  int p_index = round(p*precision_size);
  double ret = 0; // return value for missing data
  if (S==1) {
    ret = log( (1.0-e1)*(*lookup_homo_ptr)[p_index][m-1] + 0.5*e2*(*lookup_het_ptr)[p_index][m-1]);
  } else if (S==3) {
    ret = log( (1.0-e1)*(*lookup_homo_ptr)[precision_size-p_index][m-1] + 0.5*e2*(*lookup_het_ptr)[p_index][m-1]);
  } else if (S==2) {
    ret = log( e1*(*lookup_homo_ptr)[p_index][m-1] + e1*(*lookup_homo_ptr)[precision_size-p_index][m-1] + (1.0-e2)*(*lookup_het_ptr)[p_index][m-1] );
  }
  ret = (ret < -OVERFLO) ? -OVERFLO : ret;    // catch -inf values
  return ret;
}

//########################################################################################################

//------------------------------------------------
// reset particle
void particle_biallelic::reset() {
  
  // reset qmatrices
  for (int i=0; i<n; i++) {
    fill(qmatrix[i].begin(), qmatrix[i].end(), 0);
    fill(log_qmatrix[i].begin(), log_qmatrix[i].end(), 0);
  }
  
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
      loglike_old[i][j] = logprob_genotype((*data_ptr)[i][j], p[group[i]][j], m[i], e1, e2);
      loglike += loglike_old[i][j];
    }
  }
  
}

//------------------------------------------------
// update error parameters e1 or e2
void particle_biallelic::update_e(int which_e, bool robbins_monro_on, int iteration) {

  // define sum over old and new likelihood
  double sum_loglike_old = 0;
  double sum_loglike_new = 0;

  // propose new value
  double e_prop;
  if (which_e==1) {
    e_prop = rnorm1_interval(e1, e1_propSD, 0, e1_max);
  } else {
    e_prop = rnorm1_interval(e2, e2_propSD, 0, e2_max);
  }

  // calculate likelihood
  for (int i=0; i<n; i++) {
    for (int j=0; j<L; j++) {
      if (which_e==1) {
        loglike_new[i][j] = logprob_genotype((*data_ptr)[i][j], p[group[i]][j], m[i], e_prop, e2);
      } else {
        loglike_new[i][j] = logprob_genotype((*data_ptr)[i][j], p[group[i]][j], m[i], e1, e_prop);
      }
      sum_loglike_old += loglike_old[i][j];
      sum_loglike_new += loglike_new[i][j];
    }
  }

  // catch impossible proposed values
  if (sum_loglike_new < -OVERFLO) {
    return;
  }

  // Metropolis step
  if (log(runif_0_1())<beta_raised*(sum_loglike_new - sum_loglike_old)) {

    // update loglike
    for (int i=0; i<n; i++) {
      for (int j=0; j<L; j++) {
        loglike_old[i][j] = loglike_new[i][j];
      }
    }

    // update selected parameter, apply Robbins-Monro positive update step, and
    // update acceptance rate.
    if (which_e==1) {
      e1 = e_prop;
      if (robbins_monro_on) {
        e1_propSD  += (1-0.23)/sqrt(double(iteration));
      }
      e1_accept++;
    } else {
      e2 = e_prop;
      if (robbins_monro_on) {
        e2_propSD  += (1-0.23)/sqrt(double(iteration));
      }
      e2_accept++;
    }

  } else {

    // Robbins-Monro negative update step
    if (which_e==1) {
      if (robbins_monro_on) {
        e1_propSD  -= 0.23/sqrt(double(iteration));
        if (e1_propSD < UNDERFLO) {
          e1_propSD = UNDERFLO;
        }
      }
    } else {
      if (robbins_monro_on) {
        e2_propSD  -= 0.23/sqrt(double(iteration));
        if (e2_propSD < UNDERFLO) {
          e2_propSD = UNDERFLO;
        }
      }
    }

  } // end Metropolis step

}

//------------------------------------------------
// update allele frequencies p
void particle_biallelic::update_p(bool robbins_monro_on, int iteration) {

  // loop through loci
  for (int j=0; j<L; j++) {

    // clear vectors storing sum over old and new likelihood
    fill(sum_loglike_old_vec.begin(), sum_loglike_old_vec.end(), 0);
    fill(sum_loglike_new_vec.begin(), sum_loglike_new_vec.end(), 0);

    // draw p for all demes
    for (int k=0; k<K; k++) {
      p_prop[k] = rnorm1_interval(p[k][j], p_propSD[k][j], 0, 1.0);
    }

    // calculate likelihood
    for (int i=0; i<n; i++) {
      loglike_new[i][j] = logprob_genotype((*data_ptr)[i][j], p_prop[group[i]], m[i], e1, e2);
      sum_loglike_old_vec[group[i]] += loglike_old[i][j];
      sum_loglike_new_vec[group[i]] += loglike_new[i][j];
    }

    // TODO - add lambda prior

    // Metropolis step for all K
    for (int k=0; k<K; k++) {

      // catch impossible proposed values
      if (sum_loglike_new_vec[k] < -OVERFLO) {
        continue;
      }

      // Metropolis step
      if (log(runif_0_1())<beta_raised*(sum_loglike_new_vec[k] - sum_loglike_old_vec[k])) {

        // update p
        p[k][j] = p_prop[k];

        // update loglike in all individuals from this group
        for (int i=0; i<n; i++) {
          if (group[i]==k) {
            loglike_old[i][j] = loglike_new[i][j];
          }
        }

        // Robbins-Monro positive update
        if (robbins_monro_on) {
          p_propSD[k][j]  += (1-0.23)/sqrt(double(iteration));
        }

        // update acceptance rates
        p_accept[k][j]++;

      } else {

        // Robbins-Monro negative update
        if (robbins_monro_on) {
          p_propSD[k][j]  -= 0.23/sqrt(double(iteration));
          if (p_propSD[k][j] < UNDERFLO) {
            p_propSD[k][j] = UNDERFLO;
          }
        }

      }
    }  // end Metropolis step
  }   // end loop through loci

}

//------------------------------------------------
// linear update of m
void particle_biallelic::update_m() {

  // define sum over old and new likelihood
  double sum_loglike_old = 0;
  double sum_loglike_new = 0;

  // loop through individuals
  for (int i=0; i<n; i++) {
    int this_group = group[i];

    // propose new m
    int m_prop = rbernoulli1(0.5);
    m_prop = (m_prop==0) ? m[i]-1 : m[i]+1;

    // TODO - skip if proposed m impossible?
    //if (m_prop==1 && (*anyHet_ptr)[i]) {
    //    continue;
    //}

    // calculate likelihood
    if (m_prop>0 && m_prop<=COI_max) {

      // reset sum over old and new likelihood
      sum_loglike_old = 0;
      sum_loglike_new = 0;

      // calculate likelihood
      for (int j=0; j<L; j++) {
        loglike_new[i][j] = logprob_genotype((*data_ptr)[i][j], p[this_group][j], m_prop, e1, e2);
        sum_loglike_old += beta_raised*loglike_old[i][j];
        sum_loglike_new += beta_raised*loglike_new[i][j];
      }

      // catch impossible proposed values
      if (sum_loglike_new < -OVERFLO) {
        continue;
      }

      // apply Poisson or negative Binomial prior
      if (COI_model==2) {
        sum_loglike_old += dpois1(m[i]-1, COI_mean[this_group]-1);
        sum_loglike_new += dpois1(m_prop-1, COI_mean[this_group]-1);
      } else if (COI_model==3) {
        sum_loglike_old += dnbinom1(m[i]-1, COI_mean[this_group]-1, COI_dispersion);
        sum_loglike_new += dnbinom1(m_prop-1, COI_mean[this_group]-1, COI_dispersion);
      }

      // Metropolis step
      if (log(runif_0_1()) < (sum_loglike_new - sum_loglike_old)) {
        m[i] = m_prop;
        for (int j=0; j<L; j++) {
          loglike_old[i][j] = loglike_new[i][j];
        }
      }
    }

  }   // end loop through individuals

}

//------------------------------------------------
// update group allocation
void particle_biallelic::update_group() {

  // no change if K==1
  if (K==1) {
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
      if (k==this_group) {  // likelihood already calculated for current group
        for (int j=0; j<L; j++) {
          loglike_new_group[k][j] = loglike_old[i][j];
          log_qmatrix[i][k] += beta_raised*loglike_new_group[k][j];
        }
      } else {  // calculate likelihood for group k
        for (int j=0; j<L; j++) {
          loglike_new_group[k][j] = logprob_genotype((*data_ptr)[i][j], p[k][j], m[i], e1, e2);
          log_qmatrix[i][k] += beta_raised*loglike_new_group[k][j];
        }
      }

      // add prior information
      if (COI_model==2) {
        log_qmatrix[i][k] += dpois1(m[i]-1, COI_mean[k]-1);
      } else if (COI_model==3) {
        log_qmatrix[i][k] += dnbinom1(m[i]-1, COI_mean[k]-1, COI_dispersion);
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
void particle_biallelic::update_COI_mean(bool robbins_monro_on, int iteration) {
  
  // split method between poisson and negative binomial
  // poisson model
  if (COI_model==2) {
    
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
      COI_mean[k] = rgamma1(COI_mean_shape[k], COI_mean_rate[k])+1;
    }
  } // end poisson model
  
  // negative binomial model
  if (COI_model==3) {
    
    // clear vectors storing sum over old and new likelihood
    fill(sum_loglike_old_vec.begin(), sum_loglike_old_vec.end(), 0);
    fill(sum_loglike_new_vec.begin(), sum_loglike_new_vec.end(), 0);
    
    // draw COI_mean for all demes
    for (int k=0; k<K; k++) {
      COI_mean_prop[k] = rnorm1_interval(COI_mean[k], COI_mean_propSD[k], 1, COI_max);
    }
    
    // calculate likelihood
    for (int i=0; i<n; i++) {
      sum_loglike_old_vec[group[i]] += dnbinom1(m[i]-1, COI_mean[group[i]]-1, COI_dispersion);
      sum_loglike_new_vec[group[i]] += dnbinom1(m[i]-1, COI_mean_prop[group[i]]-1, COI_dispersion);
    }
    
    // Metropolis step for all demes
    for (int k=0; k<K; k++) {
      
      // Metropolis step
      if (log(runif_0_1())<(sum_loglike_new_vec[k] - sum_loglike_old_vec[k])) {
        COI_mean[k] = COI_mean_prop[k];
        
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
// get group order such that group_order[group[i]] is always-increasing over i
void particle_biallelic::get_group_order() {
  
  // get order of unique values in group
  vector<int> group_uniques = unique_int(group, K);
  group_order = order_unique_int(group_uniques, K);
}

//------------------------------------------------
// return log-probability of a particular scaffold group being drawn
/*
double particle_biallelic::scaf_prop_logprob(const vector<int> &prop_group) {
  
  int matches = 0;
  for (int i=0; i<scaffold_n; i++) {
    bool exact_match = true;
    for (int j=0; j<n; j++) {
      if (prop_group[j] != scaf_group[i][j]) {
        exact_match = false;
        break;
      }
    }
    if (exact_match) {
      matches ++;
    }
  }
  double ret = log(matches/double(scaffold_n));
  return ret;
}
*/
//------------------------------------------------
// propose swap with scaffold grouping
void particle_biallelic::scaf_propose(int &scaf_accept) {
  
  // re-order group allocation to be always-increasing
  //vector<int> group_inc = group_increasing();
  
  // initialise backwards proposal probability
  double prop_backwards = 0;//scaf_prop_logprob();
  
  // break if no chance of backwards move
  if (std::isinf(prop_backwards)) {
    return;
  }
  //print(prop_backwards);
  /*
  // store old group and propose new group from scaffolds
  vector<int> group_old = group;
  int rand1 = sample2(0, scaf_n-1);
  double prop_forwards = scaf_prop_logprob(scaf_group[rand1]);

  // propose new mu
  for (int k=0; k<K; k++) {
    double mu_post_var = 1/(beta*scaf_counts[rand1][k]/sigma_sq + 1/mu_prior_var);
    double mu_post_mean = mu_post_var * (beta*scaf_x_sum[rand1][k]/sigma_sq + mu_prior_mean/mu_prior_var);
    double mu_post_sd = sqrt(mu_post_var);
    scaf_mu[k] = rnorm1(mu_post_mean, mu_post_sd);
    prop_forwards += dnorm1(scaf_mu[k], mu_post_mean, mu_post_sd);
  }

  // calculate backwards probability of drawing current mu
  for (int k=0; k<K; k++) {
    double mu_post_var = 1/(beta*counts[k]/sigma_sq + 1/mu_prior_var);
    double mu_post_mean = mu_post_var * (beta*x_sum[k]/sigma_sq + mu_prior_mean/mu_prior_var);
    double mu_post_sd = sqrt(mu_post_var);
    prop_backwards += dnorm1(mu[k], mu_post_mean, mu_post_sd);
  }

  // calculate new likelihood
  double loglike_new = 0;
  for (int i=0; i<n; i++) {
    loglike_new += dnorm1((*x_ptr)[i], scaf_mu[scaf_group[rand1][i]], sigma);
  }

  // calculate old and new priors
  double logprior_old = 0;
  double logprior_new = 0;
  for (int k=0; k<K; k++) {
    logprior_old += dnorm1(mu[k], mu_prior_mean, sqrt(mu_prior_var));
    logprior_new += dnorm1(scaf_mu[k], mu_prior_mean, sqrt(mu_prior_var));
  }

  // Metropolis-Hastings
  double MH = ((beta*loglike_new + logprior_new) - (beta*loglike + logprior_old)) - (prop_forwards - prop_backwards);
  if (log(runif1()) < MH) {

    // jump to scaffold
    group = scaf_group[rand1];
    counts = scaf_counts[rand1];
    x_sum = scaf_x_sum[rand1];
    mu = scaf_mu;
    loglike = loglike_new;

    // update acceptance rate
    scaf_accept ++;
  }
*/
}

//------------------------------------------------
// split-merge proposal
void particle_biallelic::splitmerge_propose(int &splitmerge_accept) {
/*
  // only for K>=3
  if (K<3) {
    return;
  }

  // zero counts etc.
  fill(splitmerge_counts.begin(), splitmerge_counts.end(), 0);
  fill(splitmerge_x_sum.begin(), splitmerge_x_sum.end(), 0);

  // sample 3 groups (without replacement). Merge first to second, split third
  // with first. Count the number of splits that take place in the forward
  // proposal, and that must take place to do the reverse move.
  sample3(splitmerge_targets, 3);
  int n_forwards = 0;
  int n_backwards = 0;
  for (int i=0; i<n; i++) {

    // split and merge groups
    if (group[i]==splitmerge_targets[0]) {
      splitmerge_group[i] = splitmerge_targets[1];
    } else if (group[i]==splitmerge_targets[2]) {
      n_forwards ++;
      if (runif_0_1()<0.5) {
        splitmerge_group[i] = splitmerge_targets[0];
      } else {
        splitmerge_group[i] = splitmerge_targets[2];
      }
    } else {
      splitmerge_group[i] = group[i];
    }
    if (splitmerge_group[i]==splitmerge_targets[1]) {
      n_backwards ++;
    }

    // update counts etc.
    splitmerge_counts[splitmerge_group[i]] ++;
    splitmerge_x_sum[splitmerge_group[i]] += (*x_ptr)[i];
  }

  // initialise forward and backward proposal probabilities
  double prop_forwards = n_forwards*log(0.5);
  double prop_backwards = n_backwards*log(0.5);

  // propose new mu
  for (int k=0; k<K; k++) {
    double mu_post_var = 1/(beta*splitmerge_counts[k]/sigma_sq + 1/mu_prior_var);
    double mu_post_mean = mu_post_var * (beta*splitmerge_x_sum[k]/sigma_sq + mu_prior_mean/mu_prior_var);
    double mu_post_sd = sqrt(mu_post_var);
    splitmerge_mu[k] = rnorm1(mu_post_mean, mu_post_sd);
    prop_forwards += dnorm1(splitmerge_mu[k], mu_post_mean, mu_post_sd);
  }

  // calculate new likelihood
  double loglike_new = 0;
  for (int i=0; i<n; i++) {
    loglike_new += dnorm1((*x_ptr)[i], splitmerge_mu[splitmerge_group[i]], sigma);
  }

  // calculate old and new priors
  double logprior_old = 0;
  double logprior_new = 0;
  for (int k=0; k<K; k++) {
    logprior_old += dnorm1(mu[k], mu_prior_mean, sqrt(mu_prior_var));
    logprior_new += dnorm1(splitmerge_mu[k], mu_prior_mean, sqrt(mu_prior_var));
  }

  // Metropolis-Hastings
  double MH = ((beta*loglike_new + logprior_new) - (beta*loglike + logprior_old)) - (prop_forwards - prop_backwards);
  if (log(runif1()) < MH) {

    // jump to split-merge proposal
    group = splitmerge_group;
    counts = splitmerge_counts;
    x_sum = splitmerge_x_sum;
    mu = splitmerge_mu;
    loglike = loglike_new;

    // update acceptance rate
    splitmerge_accept ++;
  }
*/
}

//------------------------------------------------
// fix label switching problem
void particle_biallelic::solve_label_switching(const vector<vector<double>> &log_qmatrix_running) {
/*
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
*/
}

//------------------------------------------------
// calculate overall log-likelihood
void particle_biallelic::calculate_loglike() {
  loglike = 0;
  for (int i=0; i<n; i++) {
    for (int j=0; j<L; j++) {
      loglike += loglike_old[i][j];
    }
  }
}
