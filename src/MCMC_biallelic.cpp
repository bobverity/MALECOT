
#include "MCMC_biallelic.h"
#include "misc.h"
#include "probability.h"

using namespace std;

//------------------------------------------------
// constructor for MCMC_biallelic class
MCMC_biallelic::MCMC_biallelic(Rcpp::List &args) {
  
  // extract data and parameters
  data = rcpp_to_mat_int(args["data"]);
  n = rcpp_to_int(args["n"]);
  L = rcpp_to_int(args["L"]);
  K = rcpp_to_int(args["K"]);
  COI_model = rcpp_to_int(args["COI_model_numeric"]);
  COI_max = rcpp_to_int(args["COI_max"]);
  e1 = rcpp_to_double(args["e1"]);
  e2 = rcpp_to_double(args["e2"]);
  estimate_error = rcpp_to_bool(args["estimate_error"]);
  burnin = rcpp_to_int(args["burnin"]);
  samples = rcpp_to_int(args["samples"]);
  rungs = rcpp_to_int(args["rungs"]);
  auto_converge = rcpp_to_bool(args["auto_converge"]);
  coupling_on = rcpp_to_bool(args["coupling_on"]);
  scaffold_on = rcpp_to_bool(args["scaffold_on"]);
  scaffold_n = rcpp_to_int(args["scaffold_n"]);
  scaffold_group_n = rcpp_to_int(args["scaffold_group_n"]);
  split_merge_on = rcpp_to_bool(args["split_merge_on"]);
  solve_label_switching_on = rcpp_to_bool(args["solve_label_switching_on"]);
  precision = rcpp_to_double(args["precision"]);
  GTI_pow = rcpp_to_double(args["GTI_pow"]);
  silent = rcpp_to_bool(args["silent"]);
  //output_format = rcpp_to_int(args["output_format"]);
  
  // cannot use scaffolds if none available
  if (scaffold_group_n==0) {
    scaffold_on = false;
  }
  
  // create lookup tables
  if (precision!=0) {

    // initialise tables
    precision_size = int(1/precision);
    lookup_homo = vector<vector<double>>(precision_size+1, vector<double>(COI_max));
    lookup_het = vector<vector<double>>(precision_size+1, vector<double>(COI_max));

    // populate tables
    for (int i=0; i<(precision_size+1); i++) {
      double p = i/double(precision_size);
      for (int m=0; m<COI_max; m++) {

        // if error terms are fixed then incorporate at this stage, and store
        // log probabilities. Otherwise store raw probability and error will be
        // incorporated later
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
  } // end populate lookup tables

  // thermodynamic parameters. The object beta_raised_vec stores values of beta
  // (the thermodynamic power), raised to the power GTI_pow
  beta_raised_vec = vector<double>(rungs);
  for (int rung=0; rung<rungs; rung++) {
    beta_raised_vec[rung] = (rungs==1) ? 1 : pow((rung+1)/double(rungs), GTI_pow);
  }
  rung_order = seq_int(0,rungs-1);
  cold_rung = rung_order[rungs-1];

  // vector of particles
  particle_vec = vector<particle_biallelic>(rungs);
  for (int rung=0; rung<rungs; rung++) {
    particle_vec[rung] = particle_biallelic(data, args, lookup_homo, lookup_het, beta_raised_vec[rung]);
  }
  
  // initialise ordering of labels
  label_order = seq_int(0,K-1);
  label_order_new = vector<int>(K);
  
  // Q-matrices
  log_qmatrix_running = vector<vector<double>>(n, vector<double>(K));
  if (K==1) {
    qmatrix_final = vector<vector<double>>(n, vector<double>(K, samples));
  } else {
    qmatrix_final = vector<vector<double>>(n, vector<double>(K));
  }
  
  // objects for storing results
  burnin_loglike = vector<vector<double>>(rungs, vector<double>(burnin));
  sampling_loglike = vector<vector<double>>(rungs, vector<double>(samples));
  m_store = vector<vector<int>>(samples, vector<int>(n));
  p_store = vector<vector<vector<double>>>(K, vector<vector<double>>(L, vector<double>(samples)));
  e1_store = vector<double>(samples, e1);
  e2_store = vector<double>(samples, e2);
  COI_mean_store = vector<vector<double>>(samples, vector<double>(K));
  
  // objects for storing acceptance rates
  p_accept = vector<vector<int>>(K, vector<int>(L));
  e1_accept = 0;
  e2_accept = 0;
  coupling_accept = vector<int>(rungs-1);
  scaf_trials = 0;
  scaf_accept = 0;
  split_merge_accept = vector<int>(rungs);
}

//------------------------------------------------
// generate scaffold groupings
void MCMC_biallelic::scaffold_mcmc(Rcpp::List &args) {
  
  // print header
  if (!silent) {
    print("Generating", scaffold_n, "scaffolds");
  }
  
  // extract R functions
  Rcpp::Function test_convergence = args["test_convergence"];
  Rcpp::Function update_progress = args["update_progress"];
  
  // define points at which convergence checked
  vector<int> convergence_checkpoint(10);
  for (int i=0; i<10; i++) {
    convergence_checkpoint[i] = (i+1)*double(burnin)/10;
  }
  int checkpoint_i = 0;
  
  // initialise scaffold progress bar
  if (!silent) {
    update_progress(args, 1, 0, scaffold_n);
  }
  
  // create scaffold particle
  particle_biallelic particle_scaf(data, args, lookup_homo, lookup_het, 1);
  
  // objects for storing results
  scaffold_group = vector<vector<int>>(scaffold_n, vector<int>(n));
  scaffold_loglike = vector<double>(scaffold_n);
  
  // generate multiple scaffolds
  int n_non_converge = 0;
  for (int scaf_rep=0; scaf_rep<scaffold_n; scaf_rep++) {
    
    // reset particle
    particle_scaf.reset();
    
    // store log-likelihoods for checking convergence
    vector<double> scaffold_loglike_conv(burnin);
    
    // loop through burn-in iterations
    bool has_converged = false;
    for (int rep=0; rep<burnin; rep++) {
      
      // update error estimates
      if (estimate_error) {
        particle_scaf.update_e(1, true, rep+1);
        particle_scaf.update_e(2, true, rep+1);
      }
      
      // update p
      particle_scaf.update_p(true, rep+1);
      
      // update m
      particle_scaf.update_m();
      
      // update COI_means
      if (COI_model==2 || COI_model==3) {
        particle_scaf.update_COI_mean(true, rep+1);
      }
      
      // update group
      particle_scaf.update_group();
      
      // split-merge step
      //if (splitmerge_on) {
      //  particle_scaf.splitmerge_propose(dummy_accept);
      //}
      
      // calculate log-likelihood
      particle_scaf.calculate_loglike();
      
      // store log-likelihood
      scaffold_loglike_conv[rep] = particle_scaf.loglike;
      
      // check for convergence
      if (auto_converge && (rep+1)==convergence_checkpoint[checkpoint_i]) {
        
        // break if converged
        has_converged = rcpp_to_bool(test_convergence(scaffold_loglike_conv, rep+1));
        if (has_converged) {
          break;
        }
        checkpoint_i++;
      }
      
    } // end scaffold burn-in iterations
    
    // store final likelihood
    scaffold_loglike[scaf_rep] = particle_scaf.loglike;
    
    // check if still not converged at end of batches
    if (!has_converged) {
      n_non_converge++;
    }
    
    // re-order group allocation to be always-increasing, and store this scaffold grouping
    particle_scaf.get_group_increasing();
    scaffold_group[scaf_rep] = particle_scaf.group_increasing;
    
    // update scaffold progress bar
    if (!silent) {
      update_progress(args, 1, scaf_rep+1, scaffold_n);
    }
    
  } // end loop over scaffold_n
  
  // throw warning if some scaffolds not converged
  if (n_non_converge>0) {
    print("   Warning:", n_non_converge, "out of", scaffold_n,"scaffolds did not converge"); 
  }
  
}

//------------------------------------------------
// run burn-in phase of MCMC
void MCMC_biallelic::burnin_mcmc(Rcpp::List &args) {
  
  // print header
  if (!silent) {
    print("Burn-in phase");
  }
  
  // extract R functions
  Rcpp::Function test_convergence = args["test_convergence"];
  Rcpp::Function update_progress = args["update_progress"];
  
  // define points at which convergence checked
  vector<int> convergence_checkpoint(10);
  for (int i=0; i<10; i++) {
    convergence_checkpoint[i] = (i+1)*double(burnin)/10;
  }
  int checkpoint_i = 0;
  
  // reset particles
  for (int r=0; r<rungs; r++) {
    particle_vec[r].reset();
    particle_vec[r].beta_raised = beta_raised_vec[r];
  }
  rung_order = seq_int(0,rungs-1);
  
  // loop through burnin iterations
  bool convergence_reached = false;
  for (int rep=0; rep<burnin; rep++) {
    
    // update particles
    for (int r=0; r<rungs; r++) {
      int rung = rung_order[r];
      
      // update error estimates
      if (estimate_error) {
        particle_vec[rung].update_e(1, true, rep+1);
        particle_vec[rung].update_e(2, true, rep+1);
      }
      
      // update p
      particle_vec[rung].update_p(true, rep+1);
      
      // update m
      particle_vec[rung].update_m();
      
      // update COI_means
      if (COI_model==2 || COI_model==3) {
        particle_vec[rung].update_COI_mean(true, rep+1);
      }
      
      // update group
      particle_vec[rung].update_group();
      
      // calculate log-likelihood
      particle_vec[rung].calculate_loglike();
      
      // split-merge step
      if (split_merge_on) {
        //particle_vec[rung].splitmerge_propose(splitmerge_accept[r]);
      }
      
    } // end loop over rungs
    
    // Metropolis-coupling
    if (coupling_on) {
      metropolis_coupling();
    }
    
    // focus on cold rung
    cold_rung = rung_order[rungs-1];
    
    // methods that only apply when K>1
    if (K>1) {
      
      // propose swap with scaffolds
      if (scaffold_on) {
        particle_vec[cold_rung].scaf_propose(scaf_trials, scaf_accept);
      }
      
      // fix labels
      if (solve_label_switching_on) {
        particle_vec[cold_rung].solve_label_switching(log_qmatrix_running);
        label_order = particle_vec[cold_rung].label_order;
      }
      
      // add particle log_qmatrix to log_qmatrix_running
      for (int i=0; i<n; i++) {
        for (int k=0; k<K; k++) {
          log_qmatrix_running[i][k] = log_sum(log_qmatrix_running[i][k], particle_vec[cold_rung].log_qmatrix[i][label_order[k]]);
        }
      }
      
    }
    
    // store loglikelihood
    for (int r=0; r<rungs; r++) {
      int rung = rung_order[r];
      burnin_loglike[r][rep] = particle_vec[rung].loglike;
    }
    
    // update progress bars
    if (!silent) {
      int remainder = rep % int(ceil(double(burnin)/100));
      if (remainder==0 || (rep+1)==burnin) {
        update_progress(args, 2, rep+1, burnin);
      }
    }
    
    // check for convergence
    if (auto_converge && (rep+1)==convergence_checkpoint[checkpoint_i]) {
      
      // break if convergence reached
      convergence_reached = true;
      for (int r=0; r<rungs; r++) {
        bool this_converged = rcpp_to_bool(test_convergence(burnin_loglike[r], rep+1));
        if (!this_converged) {
          convergence_reached = false;
          break;
        }
      }
      if (convergence_reached) {
        if (!silent) {
          update_progress(args, 2, burnin, burnin);
          print("   converged within", rep+1, "iterations");
        }
        for (int r=0; r<rungs; r++) {
          burnin_loglike[r].resize(rep+1);
        }
        break;
      }
      checkpoint_i++;
    }
    
  } // end burn-in iterations
  
  // warning if still not converged
  if (!convergence_reached) {
    print("   Warning: convergence still not reached within", burnin, "iterations");
  }
}

//------------------------------------------------
// run sampling phase of MCMC
void MCMC_biallelic::sampling_mcmc(Rcpp::List &args) {
  
  // print header
  if (!silent) {
    print("Sampling phase");
  }
  
  // extract R functions
  Rcpp::Function update_progress = args["update_progress"];
  
  // reset acceptance rates
  for (int r=0; r<rungs; r++) {
    particle_vec[r].p_accept = vector<vector<int>>(K, vector<int>(L));
    particle_vec[r].e1_accept = 0;
    particle_vec[r].e2_accept = 0;
  }
  coupling_accept = vector<int>(rungs-1);
  scaf_trials = 0;
  scaf_accept = 0;
  split_merge_accept = vector<int>(rungs);
  
  // loop through sampling iterations
  for (int rep=0; rep<samples; rep++) {
    
    // update particles
    for (int r=0; r<rungs; r++) {
      int rung = rung_order[r];
      
      // update error estimates
      if (estimate_error) {
        particle_vec[rung].update_e(1, false, rep+1);
        particle_vec[rung].update_e(2, false, rep+1);
      }
      
      // update p
      particle_vec[rung].update_p(false, rep+1);
      
      // update m
      particle_vec[rung].update_m();
      
      // update COI_means
      if (COI_model==2 || COI_model==3) {
        particle_vec[rung].update_COI_mean(false, rep+1);
      }
      
      // update group
      particle_vec[rung].update_group();
      
      // calculate log-likelihood
      particle_vec[rung].calculate_loglike();
      
      // split-merge step
      if (split_merge_on) {
        //particle_vec[rung].splitmerge_propose(splitmerge_accept[r]);
      }
      
    } // end loop over rungs
    
    // Metropolis-coupling
    if (coupling_on) {
      metropolis_coupling();
    }
    
    // focus on cold rung
    cold_rung = rung_order[rungs-1];
    
    // methods that only apply when K>1
    if (K>1) {
      
      // propose swap with scaffolds
      if (scaffold_on) {
        particle_vec[cold_rung].scaf_propose(scaf_trials, scaf_accept);
      }
      
      // fix labels
      if (solve_label_switching_on) {
        particle_vec[cold_rung].solve_label_switching(log_qmatrix_running);
        label_order = particle_vec[cold_rung].label_order;
      }
      
      // add particle log_qmatrix to log_qmatrix_running
      for (int i=0; i<n; i++) {
        for (int k=0; k<K; k++) {
          log_qmatrix_running[i][k] = log_sum(log_qmatrix_running[i][k], particle_vec[cold_rung].log_qmatrix[i][label_order[k]]);
        }
      }
      
      // add particle qmatrix to qmatrix_final
      for (int i=0; i<n; i++) {
        for (int k=0; k<K; k++) {
          qmatrix_final[i][k] += particle_vec[cold_rung].qmatrix[i][label_order[k]];
        }
      }
      
    }
    
    // store loglikelihood
    for (int r=0; r<rungs; r++) {
      int rung = rung_order[r];
      sampling_loglike[r][rep] = particle_vec[rung].loglike;
    }
    
    // store m and p
    m_store[rep] = particle_vec[cold_rung].m;
    for (int k=0; k<K; k++) {
      for (int j=0; j<L; j++) {
        p_store[k][j][rep] = particle_vec[cold_rung].p[label_order[k]][j];
      }
    }
    
    // store e1 and e2
    if (estimate_error) {
      e1_store[rep] = particle_vec[cold_rung].e1;
      e2_store[rep] = particle_vec[cold_rung].e2;
    }
    
    // store COI_mean
    if (COI_model==2 || COI_model==3) {
      COI_mean_store[rep] = particle_vec[cold_rung].COI_mean;
    }
    
    // update progress bars
    if (!silent) {
      int remainder = rep % int(ceil(double(samples)/100));
      if (remainder==0 || (rep+1)==samples) {
        update_progress(args, 3, rep+1, samples);
      }
    }
    
  } // end sampling iterations
  
  // finalise Q-matrix
  for (int i=0; i<n; i++) {
    for (int j=0; j<K; j++) {
      qmatrix_final[i][j] /= samples;
    }
  }
  
  // store acceptance rates
  p_accept = particle_vec[cold_rung].p_accept;
  e1_accept = particle_vec[cold_rung].e1_accept;
  e2_accept = particle_vec[cold_rung].e2_accept;
}

//------------------------------------------------
// add qmatrix to qmatrix_running
void MCMC_biallelic::update_qmatrix_running() {
/*
  for (int i=0; i<n; i++) {
    for (int k=0; k<K; k++) {
      qmatrix_running[i][k] += particle_vec[cold_rung].qmatrix[i][particle_vec[cold_rung].label_order[k]];
      double lq = log(qmatrix_running[i][k]);
      if (lq < -OVERFLO) {
        log_qmatrix_running[i][k] = -OVERFLO;
      } else {
        log_qmatrix_running[i][k] = lq;
      }
    }
  }
*/
}

//------------------------------------------------
// add qmatrix to qmatrix_final
void MCMC_biallelic::update_qmatrix_final() {
/*
  for (int i=0; i<n; i++) {
    for (int k=0; k<K; k++) {
      qmatrix_final[i][k] += particle_vec[cold_rung].qmatrix[i][particle_vec[cold_rung].label_order[k]];
    }
  }
*/
}

//------------------------------------------------
// Metropolis-coupling to propose swaps between temperature rungs
void MCMC_biallelic::metropolis_coupling() {
  
  // loop over rungs, starting with the hottest chain and moving to the cold
  // chain. Each time propose a swap with the next rung up
  for (int i=0; i<(rungs-1); i++) {
    
    // define rungs of interest
    int rung1 = rung_order[i];
    int rung2 = rung_order[i+1];
    
    // get log-likelihoods and beta values of two chains in the comparison
    double log_like1 = particle_vec[rung1].loglike;
    double log_like2 = particle_vec[rung2].loglike;
    
    double beta_raised1 = particle_vec[rung1].beta_raised;
    double beta_raised2 = particle_vec[rung2].beta_raised;
    
    // calculate acceptance ratio (still in log space)
    double acceptance = (log_like2*beta_raised1 + log_like1*beta_raised2) - (log_like1*beta_raised1 + log_like2*beta_raised2);

    // accept or reject move
    double rand1 = runif1();
    if (log(rand1)<acceptance) {
      
      // swap beta values
      particle_vec[rung1].beta_raised = beta_raised2;
      particle_vec[rung2].beta_raised = beta_raised1;
      
      // swap rung order
      rung_order[i] = rung2;
      rung_order[i+1] = rung1;
      
      // update acceptance rates
      coupling_accept[i]++;
    }
  }
  
}
