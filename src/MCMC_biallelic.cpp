
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
  //lambda;
  COI_model = rcpp_to_int(args["COI_model_numeric"]);
  COI_max = rcpp_to_int(args["COI_max"]);
  //COI_dispersion;
  e1 = rcpp_to_double(args["e1"]);
  e2 = rcpp_to_double(args["e2"]);
  estimate_error = rcpp_to_bool(args["estimate_error"]);
  //e1_max;
  //e2_max;
  burnin = rcpp_to_int(args["burnin"]);
  samples = rcpp_to_int(args["samples"]);
  rungs = rcpp_to_int(args["rungs"]);
  auto_converge = rcpp_to_bool(args["auto_converge"]);
  coupling_on = rcpp_to_bool(args["coupling_on"]);
  scaffold_on = rcpp_to_bool(args["scaffold_on"]);
  scaffold_n = rcpp_to_int(args["scaffold_n"]);
  split_merge_on = rcpp_to_bool(args["split_merge_on"]);
  solve_label_switching_on = rcpp_to_bool(args["solve_label_switching_on"]);
  precision = rcpp_to_double(args["precision"]);
  GTI_pow = rcpp_to_double(args["GTI_pow"]);
  silent = rcpp_to_bool(args["silent"]);
  //output_format = rcpp_to_int(args["output_format"]);
  
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
    beta_raised_vec[rung] = (rungs==1) ? 1 : pow(rung/double(rungs-1), GTI_pow);
  }
  rung_order = seq_int(0,rungs-1);
  cold_rung = rung_order[rungs-1];

  // vector of particles
  particle_vec = vector<particle_biallelic>(rungs);
  for (int rung=0; rung<rungs; rung++) {
    particle_vec[rung] = particle_biallelic(data, args, lookup_homo, lookup_het, beta_raised_vec[rung]);
  }
  
  // scaffold objects
  scaf_group = vector<vector<vector<int>>>(rungs, vector<vector<int>>(scaffold_n, vector<int>(n)));
  scaf_count = vector<vector<int>>(rungs, vector<int>(scaffold_n));
  
  // initialise ordering of labels
  label_order = seq_int(0,K-1);
  label_order_new = vector<int>(K);
  
  // objects for storing results
  burnin_loglike = vector<vector<double>>(rungs, vector<double>(burnin));
  sampling_loglike = vector<vector<double>>(rungs, vector<double>(samples));
  log_qmatrix_running = vector<vector<double>>(n, vector<double>(K));
  qmatrix_final = vector<vector<double>>(n, vector<double>(K));
  
  // objects for storing acceptance rates
  p_accept = vector<vector<int>>(K, vector<int>(L));
  e1_accept = 0;
  e2_accept = 0;
  coupling_accept = vector<int>(rungs-1);
  scaf_accept = vector<int>(rungs);
  splitmerge_accept = vector<int>(rungs);
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
  
  // generate multiple scaffolds
  int n_non_converge = 0;
  for (int scaf_rep=0; scaf_rep<scaffold_n; scaf_rep++) {
    
    // reset particles
    for (int r=0; r<rungs; r++) {
      particle_vec[r].reset();
      particle_vec[r].beta_raised = beta_raised_vec[r];
    }
    rung_order = seq_int(0,rungs-1);
    
    // store log-likelihoods for checking convergence
    vector<vector<double>> scaf_log_like(rungs, vector<double>(burnin));
    
    // loop through burn-in iterations
    bool all_converged = false;
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
        
        // split-merge step
        //if (splitmerge_on) {
        //  particle_vec[rung].splitmerge_propose(dummy_accept);
        //}
        
        // calculate log-likelihood
        particle_vec[rung].calculate_loglike();
        
      } // end loop over rungs
      
      // Metropolis-coupling
      if (coupling_on) {
        metropolis_coupling();
      }
      
      // store log-likelihoods
      for (int r=0; r<rungs; r++) {
        int rung = rung_order[r];
        scaf_log_like[r][rep] = particle_vec[rung].loglike;
      }
      
      // check for convergence
      if (auto_converge && (rep+1)==convergence_checkpoint[checkpoint_i]) {
        
        // break if all rungs have converged
        all_converged = true;
        for (int r=0; r<rungs; r++) {
          int rung = rung_order[r];
          bool rung_converged = rcpp_to_bool(test_convergence(scaf_log_like[rung], rep+1));
          if (!rung_converged) {
            all_converged = false;
            break;
          }
        }
        if (all_converged) {
          break;
        }
        checkpoint_i++;
      }
      
    } // end scaffold burn-in iterations
    
    // check if still not converged at end of batches
    if (!all_converged) {
      n_non_converge++;
    }
    
    // loop over rungs
    for (int r=0; r<rungs; r++) {
      int rung = rung_order[r];
      
      // re-order group allocation to be always-increasing, and store this scaffold grouping
      particle_vec[rung].get_group_order();
      for (int i=0; i<n; i++) {
        scaf_group[rung][scaf_rep][i] = particle_vec[rung].group_order[particle_vec[rung].group[i]];
      }
    }
    
    // update scaffold progress bar
    if (!silent) {
      update_progress(args, 1, scaf_rep+1, scaffold_n);
    }
    
  } // end loop over scaffold_n
  
  // throw warning if some scaffolds not converged
  if (n_non_converge>0) {
    print("   Warning:", n_non_converge, "out of", scaffold_n,"scaffolds did not converge"); 
  }
  
  // count number of times each scaffold is represented
  if (scaffold_n>1) {
    for (int r=0; r<rungs; r++) {
      
      // get unique scaffolds
      vector<vector<int>> scaf_unique;
      scaf_unique.push_back(scaf_group[r][0]);
      vector<int> unique_count(1,1);
      vector<int> match_unique(1);
      
      for (int i=1; i<scaffold_n; i++) {
        bool new_unique = true;
        for (int j=0; j<int(scaf_unique.size()); j++) {
          if (vectors_identical(scaf_group[r][i], scaf_unique[j])) {
            new_unique = false;
            unique_count[j]++;
            match_unique.push_back(j);
            break;
          }
        }
        if (new_unique) {
          scaf_unique.push_back(scaf_group[r][i]);
          unique_count.push_back(1);
          match_unique.push_back(scaf_unique.size()-1);
        }
      }
      
      // store number of times each scaf_group is represented in unique list
      for (int i=0; i<scaffold_n; i++) {
        scaf_count[r][i] = unique_count[match_unique[i]];
      }
      
    } // end loop over rungs
  }
  
  // load these scaffolds back into particles
  for (int r=0; r<rungs; r++) {
    particle_vec[r].scaf_group = scaf_group[r];
    particle_vec[r].scaf_count = scaf_count[r];
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
      
      // propose swap with scaffolds
      if (scaffold_on) {
        particle_vec[rung].scaf_propose(scaf_accept[r]);
      }
      
      // split-merge step
      if (split_merge_on) {
        //particle_vec[rung].splitmerge_propose(splitmerge_accept[r]);
      }
      
    } // end loop over rungs
    
    // Metropolis-coupling
    if (coupling_on) {
      metropolis_coupling();
    }
    
    // store loglikelihood
    for (int r=0; r<rungs; r++) {
      int rung = rung_order[r];
      burnin_loglike[r][rep] = particle_vec[rung].loglike;
    }
    
    // methods that only apply when K>1
    if (K>1) {
      
      // focus on cold rung
      cold_rung = rung_order[rungs-1];
      
      // fix labels
      if (solve_label_switching_on) {
        //fixLabels(particleMat[0][coldRung]);
      }
      
      // add particle log_qmatrix to log_qmatrix_running
      for (int i=0; i<n; i++) {
        for (int k=0; k<K; k++) {
          log_qmatrix_running[i][k] = log_sum(log_qmatrix_running[i][k], particle_vec[cold_rung].log_qmatrix[i][label_order[k]]);
        }
      }
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
      convergence_reached = rcpp_to_bool(test_convergence(burnin_loglike[cold_rung], rep+1));
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
      
      // propose swap with scaffolds
      if (scaffold_on) {
        //particle_vec[rung].scaf_propose(scaf_accept[r]);
      }
      
      // split-merge step
      if (split_merge_on) {
        //particle_vec[rung].splitmerge_propose(splitmerge_accept[r]);
      }
      
    } // end loop over rungs
    
    // Metropolis-coupling
    if (coupling_on) {
      metropolis_coupling();
    }
    
    // store loglikelihood
    for (int r=0; r<rungs; r++) {
      int rung = rung_order[r];
      sampling_loglike[r][rep] = particle_vec[rung].loglike;
    }
    
    // methods that only apply when K>1
    if (K>1) {
      
      // focus on cold rung
      cold_rung = rung_order[rungs-1];
      
      // fix labels
      if (solve_label_switching_on) {
        //fixLabels(particleMat[0][coldRung]);
      }
      
      // add particle log_qmatrix to log_qmatrix_running
      for (int i=0; i<n; i++) {
        for (int k=0; k<K; k++) {
          log_qmatrix_running[i][k] = log_sum(log_qmatrix_running[i][k], particle_vec[cold_rung].log_qmatrix[i][label_order[k]]);
        }
      }
      
      // add particle qmatrix to 
      for (int i=0; i<n; i++) {
        for (int k=0; k<K; k++) {
          qmatrix_final[i][k] += particle_vec[cold_rung].qmatrix[i][label_order[k]];
        }
      }
      
    }
    
    // update progress bars
    if (!silent) {
      int remainder = rep % int(ceil(samples/100));
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
