
#include "MCMC_biallelic.h"
#include "misc.h"
#include "probability.h"

using namespace std;

//------------------------------------------------
// constructor for MCMC_biallelic class
MCMC_biallelic::MCMC_biallelic() {
  
  // thermodynamic parameters. The object beta_raised_vec stores values of beta
  // (the thermodynamic power), raised to the power GTI_pow
  beta_raised_vec = vector<double>(rungs);
  for (int rung=0; rung<rungs; rung++) {
    beta_raised_vec[rung] = (rungs==1) ? 1 : pow((rung+1)/double(rungs), GTI_pow);
  }
  rung_order = seq_int(0,rungs-1);
  cold_rung = rung_order[rungs-1];
  
  // vector of particles
  particle_vec = vector<Particle_biallelic>(rungs);
  for (int rung=0; rung<rungs; rung++) {
    particle_vec[rung] = Particle_biallelic(beta_raised_vec[rung]);
  }
  
  // initialise ordering of labels
  label_order = seq_int(0,K-1);
  
  // Q-matrices
  log_qmatrix_running = vector<vector<double>>(n, vector<double>(K));
  if (K == 1) {
    qmatrix_final = vector<vector<double>>(n, vector<double>(K, samples));
  } else {
    qmatrix_final = vector<vector<double>>(n, vector<double>(K));
  }
  
  // objects for storing results
  loglike_burnin = vector<vector<double>>(rungs, vector<double>(burnin));
  loglike_sampling = vector<vector<double>>(rungs, vector<double>(samples));
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
  
}

//------------------------------------------------
// run burn-in phase of MCMC
void MCMC_biallelic::burnin_mcmc(Rcpp::List &args_functions, Rcpp::List &args_progress) {
  
  // print header
  if (!silent) {
    print("Running MCMC for K =", K);
    print("Burn-in phase");
  }
  
  // read in R functions
  Rcpp::Function test_convergence = args_functions["test_convergence"];
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // define points at which convergence checked
  vector<int> convergence_checkpoint(1,converge_test);
  while(convergence_checkpoint.back()<burnin) {
    convergence_checkpoint.push_back(convergence_checkpoint.back()+converge_test);
  }
  int checkpoint_i = 0;
  
  // reset particles
  for (int r=0; r<rungs; r++) {
    particle_vec[r].reset();
    particle_vec[r].beta_raised = beta_raised_vec[r];
  }
  rung_order = seq_int(0,rungs-1);
  
  // loop through burnin iterations
  vector<bool> convergence_reached(rungs, false);
  bool all_convergence_reached = false;
  for (int rep=0; rep<burnin; rep++) {
    
    // update particles
    for (int r=0; r<rungs; r++) {
      
      // skip over converged rungs
      if (convergence_reached[r]) {
        continue;
      }
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
      if (COI_model == 2 || COI_model == 3) {
        if (estimate_COI_mean) {
          particle_vec[rung].update_COI_mean(true, rep+1);
        }
      }
      
      // update group
      particle_vec[rung].update_group();
      
      // calculate log-likelihood
      particle_vec[rung].calculate_loglike();
      
    } // end loop over rungs
    
    // methods that only apply when K>1
    if (K>1) {
      
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
      if (convergence_reached[r]) {
        continue;
      }
      int rung = rung_order[r];
      loglike_burnin[r][rep] = particle_vec[rung].loglike;
    }
    
    // update progress bars
    if (!silent) {
      if ((rep+1)==burnin) {
        update_progress(args_progress, "pb_burnin", rep+1, burnin);
      } else {
        int remainder = rep % int(ceil(double(burnin)/100));
        if (remainder==0 && !pb_markdown) {
          update_progress(args_progress, "pb_burnin", rep+1, burnin);
        }
      }
    }
    
    // check for convergence
    if (auto_converge && (rep+1)==convergence_checkpoint[checkpoint_i]) {
      
      // check for convergence of all unconverged chains
      for (int r=0; r<rungs; r++) {
        if (!convergence_reached[r]) {
          convergence_reached[r] = rcpp_to_bool(test_convergence(loglike_burnin[r], rep+1));
          if (convergence_reached[r]) {
            loglike_burnin[r].resize(rep+1);
          }
        }
      }
      // break if convergence reached
      all_convergence_reached = true;
      for (int r=0; r<rungs; r++) {
        if (!convergence_reached[r]) {
          all_convergence_reached = false;
          break;
        }
      }
      // end if all reached convergence
      if (all_convergence_reached) {
        if (!silent) {
          update_progress(args_progress, "pb_burnin", burnin, burnin);
          print("   converged within", rep+1, "iterations");
        }
        break;
      }
      checkpoint_i++;
      
    }  // end if auto_converge
    
  }  // end burn-in iterations
  
  // warning if still not converged
  if (!all_convergence_reached && !silent) {
    print("   Warning: convergence still not reached within", burnin, "iterations");
  }
  
}

//------------------------------------------------
// run sampling phase of MCMC
void MCMC_biallelic::sampling_mcmc(Rcpp::List &args_functions, Rcpp::List &args_progress) {
  
  // print header
  if (!silent) {
    print("Sampling phase");
  }
  
  // read in R functions
  Rcpp::Function test_convergence = args_functions["test_convergence"];
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // reset acceptance rates
  for (int r=0; r<rungs; r++) {
    particle_vec[r].p_accept = vector<vector<int>>(K, vector<int>(L));
    particle_vec[r].e1_accept = 0;
    particle_vec[r].e2_accept = 0;
  }
  coupling_accept = vector<int>(rungs-1);
  
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
        if (estimate_COI_mean) {
          particle_vec[rung].update_COI_mean(false, rep+1);
        }
      }
      
      // update group
      particle_vec[rung].update_group();
      
      // calculate log-likelihood
      particle_vec[rung].calculate_loglike();
      
    } // end loop over rungs
    
    // Metropolis-coupling
    if (coupling_on) {
      metropolis_coupling();
    }
    
    // focus on cold rung
    cold_rung = rung_order[rungs-1];
    
    // methods that only apply when K>1
    if (K>1) {
      
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
      loglike_sampling[r][rep] = particle_vec[rung].loglike;
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
      COI_mean_store[rep] = particle_vec[cold_rung].COI_mean_vec;
    }
    
    // update progress bars
    if (!silent) {
      if ((rep+1) == samples) {
        update_progress(args_progress, "pb_samples", rep+1, samples);
      } else {
        int remainder = rep % int(ceil(double(samples)/100));
        if (remainder == 0 && !pb_markdown) {
          update_progress(args_progress, "pb_samples", rep+1, samples);
        }
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
    double loglike1 = particle_vec[rung1].loglike;
    double loglike2 = particle_vec[rung2].loglike;
    
    double beta_raised1 = particle_vec[rung1].beta_raised;
    double beta_raised2 = particle_vec[rung2].beta_raised;
    
    // calculate acceptance ratio (still in log space)
    double acceptance = (loglike2*beta_raised1 + loglike1*beta_raised2) - (loglike1*beta_raised1 + loglike2*beta_raised2);

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
