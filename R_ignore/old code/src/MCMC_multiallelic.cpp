
#include <Rcpp.h>
#include "MCMC_multiallelic.h"
#include "misc.h"
#include "probability.h"
#include "Hungarian.h"

using namespace std;

//------------------------------------------------
// MCMC_multiallelic::
// constructor for MCMC_multiallelic class
MCMC_multiallelic::MCMC_multiallelic(Rcpp::List args_data, Rcpp::List args_model) {
    
    // extract raw data and basic data properties
    obsData_raw = Rcpp_to_mat_int(args_data["data"]);
    n = Rcpp_to_int(args_data["n"]);
    L = Rcpp_to_int(args_data["L"]);
    Jl = Rcpp_to_vector_int(args_data["Jl"]);
    
    // process data
    obsData = vector< vector< vector<int> > >(n, vector< vector<int> >(L));
    for (int i=0; i<int(obsData_raw.size()); i++) {
        int samp = obsData_raw[i][0] - 1;
        int loc = obsData_raw[i][1] - 1;
        int hap = obsData_raw[i][2];
        obsData[samp][loc].push_back(hap);
    }
    
    // model and MCMC parameters
    burnin = Rcpp_to_int(args_model["burnin"]);
    autoConverge = Rcpp_to_bool(args_model["autoConverge"]);
    samples = Rcpp_to_int(args_model["samples"]);
    rungs = Rcpp_to_int(args_model["rungs"]);
    MetropolisCouple = Rcpp_to_bool(args_model["MetropolisCouple"]);
    disruptivePropose = Rcpp_to_bool(args_model["disruptivePropose"]);
    chains = Rcpp_to_int(args_model["chains"]);
    parallel = Rcpp_to_bool(args_model["parallel"]);
    precision = Rcpp_to_double(args_model["precision"]);
    reportIteration = Rcpp_to_int(args_model["reportIteration"]);
    flushConsole = Rcpp_to_bool(args_model["flushConsole"]);
    solveLabelSwitching = Rcpp_to_bool(args_model["solveLabelSwitching"]);
    estimateError = Rcpp_to_bool(args_model["estimateError"]);
    lambda = 1.0;
    
    K = Rcpp_to_int(args_model["K"]);
    COI_max = Rcpp_to_int(args_model["COI_max"]);
    COI_model = Rcpp_to_int(args_model["COI_model_numeric"]);
    GTI_pow = Rcpp_to_double(args_model["GTI_pow"]);
    p_propSD_raw = Rcpp_to_double(args_model["p_propSD"]);
    p_RobbinsMonro = (p_propSD_raw<0);
    
    // lookup table for lgamma function
    lgamma_lookup_dim1 = max(Jl)+1;
    lgamma_lookup_dim2 = COI_max*n+1;
    lgamma_lookup = vector< vector<double> >(lgamma_lookup_dim1, vector<double>(lgamma_lookup_dim2));
    for (int i=0; i<lgamma_lookup_dim1; i++) {
        for (int j=0; j<lgamma_lookup_dim2; j++) {
            lgamma_lookup[i][j] = lgamma(j + i*lambda);
        }
    }
    
    // objects over chains and rungs
    particleMat = vector< vector<particle_multiallelic> >(chains, vector<particle_multiallelic>(rungs));
    betaPowMat = vector< vector<double> >(chains, vector<double>(rungs));
    rungOrderMat = vector< vector<int> >(chains);
    for (int c=0; c<chains; c++) {
        rungOrderMat[c] = seq_int(0,rungs-1);
        for (int r=0; r<rungs; r++) {
            betaPowMat[c][r] = (rungs==1) ? 1 : pow(r/double(rungs-1), GTI_pow);
            particleMat[c][r].initialise(obsData, Jl, args_model, lgamma_lookup, betaPowMat[c][r]);
        }
    }
    
    // running estimate of Qmatrix
    logQmatrix_running = vector< vector<double> >(n, vector<double>(K));
    Qmatrix_sum = vector< vector<double> >(n, vector<double>(K));
    
    // initialise ordering of labels
    labelOrder = seq_int(0,K-1);
    labelOrder_new = vector<int>(K);
    
    // initialise objects for Hungarian algorithm
    costMat = vector< vector<double> >(K, vector<double>(K));
    bestPerm = vector<int>(K);
    bestPermOrder = vector<int>(K);
    edgesLeft = vector<int>(K);
    edgesRight = vector<int>(K);
    blockedLeft = vector<int>(K);
    blockedRight = vector<int>(K);
    
    // objects for storing results
    burnin_logLike = vector< vector<double> >(chains);
    logLike_store = vector< vector<double> >(samples, vector<double>(rungs));
    m_store = vector< vector<int> >(samples, vector<int>(n));
    p_store = vector< vector< vector< vector<double> > > >(samples, vector< vector< vector<double> > >(K));
    COI_mean_store = vector< vector<double> >(samples, vector<double>(K));
    
    p_accept = vector< vector<int> >(K, vector<int>(L));
    disruptive_accept = 0;
    MC_accept = 0;
}

//------------------------------------------------
// MCMC_multiallelic::
// run burn-in phase of MCMC
void MCMC_multiallelic::burninMCMC(Rcpp::List args_functions) {
    
    // extract R functions
    Rcpp::Function fortyWinks = args_functions["fortyWinks"];
    //Rcpp::Function getGR = args_functions["getGR"];   // TODO - remove?
    Rcpp::Function testConvergence = args_functions["testConvergence"];
    
    print("   begin burn-in phase");
    
    // loop through burnin iterations
    bool convergenceReached = false;
    int rep;    // keep track of rep even when loop finishes
    for (rep=0; rep<burnin; rep++) {
        
        // report progress
        if (reportIteration>0) {
            if (((rep+1) % reportIteration)==0) {
                print("      iteration",rep+1);
                if (flushConsole) {
                    fortyWinks();
                }
            }
        }
        
        // update all chains and store loglikelihood
        for (int c=0; c<chains; c++) {
            
            // loop through rungs
            for (int r=0; r<rungs; r++) {
                
                // update m and configuration
                particleMat[c][r].update_sample();
                
                // update p
                particleMat[c][r].update_p();
                
                // update group
                particleMat[c][r].update_group();
                
                // update COI_means
                if (COI_model==2 || COI_model==3) {
                    particleMat[c][r].update_COI_mean(true, rep+1);
                }
                
                // disruptive proposal
                //if (disruptivePropose) {
                //    particleMat[c][r].propose_disruptive();
                //}
                
                // calculate loglikelihood
                particleMat[c][r].calculateLogLike();
            }
            
            // Metropolis coupling
            if (rungs>1 && MetropolisCouple) {
                MetropolisCoupling(particleMat[c], rungOrderMat[c]);
            }
            
            // focus on cold rung
            int coldRung = rungOrderMat[c][rungs-1];
            
            // store loglikelihood of cold rung
            burnin_logLike[c].push_back(particleMat[c][coldRung].logLike_allData);
            
        }
        
        // methods that only apply when K>1
        if (K>1) {
            
            // focus on cold rung of particleMat row 0
            int coldRung = rungOrderMat[0][rungs-1];
            
            // fix labels
            if (solveLabelSwitching) {
                fixLabels(particleMat[0][coldRung]);
            }
            
            // add logQmatrix_update to logQmatrix_running
            for (int i=0; i<n; i++) {
                for (int k=0; k<K; k++) {
                    logQmatrix_running[i][k] = logSum(logQmatrix_running[i][k], particleMat[0][coldRung].logQmatrix_update[i][labelOrder[k]]);
                }
            }
        }
        
        // calculate Gelman-Rubin statistic every 100 iterations
        if (autoConverge && rep>0 && ((rep+1) % int(1e2))==0) {
            
            //convergenceReached = getConvergence(getGR);
            convergenceReached = Rcpp_to_bool(testConvergence(burnin_logLike));
            
            // break if convergence reached
            if (convergenceReached) {
                print("   converged within", rep+1, "iterations");
                break;
            }
        }
        
    }   // end of burnin phase
    burnin_actual = rep+1;
    
    // TODO - remove
    // print disruptive proposal acceptance
    //for (int c=0; c<chains; c++) {
    //    for (int r=0; r<rungs; r++) {
    //        print(c, r, particleMat[c][r].disruptive_accept);
    //    }
    //}
    
    // end burnin phase message
    if (autoConverge && !convergenceReached) {
        print("   WARNING: convergence not reached within", burnin_actual, "iterations");
    }
    print("   --------------------");
    
}


//------------------------------------------------
// MCMC_multiallelic::
// run sampling phase of MCMC
void MCMC_multiallelic::runMCMC(Rcpp::List args_functions) {
    
    // extract R functions
    Rcpp::Function fortyWinks = args_functions["fortyWinks"];
    
    // use first row of burn-in matrix
    particleVec = particleMat[0];
    betaPowVec = betaPowMat[0];
    rungOrderVec = rungOrderMat[0];
    
    // begin sampling phase
    print("   begin sampling phase");
    
    // loop through MCMC iterations
    int coldRung = 0;    // keep track of coldRung even when loop finishes
    for (int rep=0; rep<samples; rep++) {
        
        // report progress
        if (reportIteration>0) {
            if (((rep+1) % reportIteration)==0) {
                print("      iteration",rep+1);
                if (flushConsole) {
                    fortyWinks();
                }
            }
        }
        
        // loop through rungs
        for (int r=0; r<rungs; r++) {
            
            // update m
            particleVec[r].update_sample();
            
            // update p
            particleVec[r].update_p();
            
            // update group
            particleVec[r].update_group();
            
            // update COI_means
            if (COI_model==2 || COI_model==3) {
                particleVec[r].update_COI_mean(false, rep+1);
            }
            
            // disruptive proposal
            //if (disruptivePropose) {
            //    particleVec[r].propose_disruptive();
            //}
            
            // calculate loglikelihood
            particleVec[r].calculateLogLike();
        }
        
        // Metropolis coupling
        if (rungs>1 && MetropolisCouple) {
            MetropolisCoupling(particleVec, rungOrderVec);
        }
        
        // focus on cold rung
        coldRung = rungOrderVec[rungs-1];
        
        // methods that only apply when K>1
        if (K>1) {
            
            // fix labels
            if (solveLabelSwitching) {
                fixLabels(particleVec[coldRung]);
            }
            
            // add logQmatrix_update to logQmatrix_running
            for (int i=0; i<n; i++) {
                for (int k=0; k<K; k++) {
                    logQmatrix_running[i][k] = logSum(logQmatrix_running[i][k], particleVec[coldRung].logQmatrix_update[i][labelOrder[k]]);
                }
            }
            
            // add Qmatrix values to sum
            for (int i=0; i<n; i++) {
                for (int k=0; k<K; k++) {
                    Qmatrix_sum[i][k] += particleVec[coldRung].Qmatrix_update[i][labelOrder[k]];
                }
            }
            
        }
        
        // store likelihood all rungs
        for (int r=0; r<rungs; r++) {
            logLike_store[rep][r] = particleVec[rungOrderVec[r]].logLike_allData;
        }
        
        // store m and p
        m_store[rep] = particleVec[coldRung].m;
        for (int k=0; k<K; k++) {
            p_store[rep][k] = particleVec[coldRung].p[labelOrder[k]];
        }
        
        // store COI_mean
        if (COI_model==2 || COI_model==3) {
            COI_mean_store[rep] = particleVec[coldRung].COI_mean;
        }
        
    }   // end of MCMC
    
    // TODO - remove
    // print disruptive proposal acceptance
    //for (int r=0; r<rungs; r++) {
    //    print(r, particleVec[r].disruptive_accept);
    //}
    
    // store acceptance rates
    p_accept = particleVec[coldRung].p_accept;
    disruptive_accept = particleVec[coldRung].disruptive_accept;
    // (MC_accept already stored)
    
    print("   sampling phase complete");
    print("   --------------------");
    
}

//------------------------------------------------
// MCMC_multiallelic::
// text
void MCMC_multiallelic::fixLabels(particle_multiallelic particle) {
    /*
    // calculate cost matrix from old and new Qmatrices
    for (int k1=0; k1<K; k1++) {
        fill(costMat[k1].begin(), costMat[k1].end(), 0);
        for (int k2=0; k2<K; k2++) {
            for (int i=0; i<n; i++) {
                costMat[k1][k2] += particle.Qmatrix_update[i][labelOrder[k1]]*(particle.logQmatrix_update[i][labelOrder[k1]]-logQmatrix_running[i][k2]);
            }
        }
    }
    
    // find best permutation of current labels using Hungarian algorithm
    bestPerm = hungarian(costMat, edgesLeft, edgesRight, blockedLeft, blockedRight);
    
    // catch returned error
    if (bestPerm[0]<0) {
        return;
    }
    
    // define bestPermOrder
    for (int k=0; k<K; k++) {
        bestPermOrder[bestPerm[k]] = k;
    }
    
    // replace old label order with new
    for (int k=0; k<K; k++) {
        labelOrder_new[k] = labelOrder[bestPermOrder[k]];
    }
    labelOrder = labelOrder_new;
    */
}

//------------------------------------------------
// MCMC_multiallelic::
// text
void MCMC_multiallelic::MetropolisCoupling(vector<particle_multiallelic> &particleVec, vector<int> &rungOrderVec) {
    /*
    // loop over rungs, starting with the hottest chain and moving to the cold chain. Each time propose a swap with a randomly chosen chain.
    for (int i1=0; i1<rungs; i1++) {
        
        // draw random value i2!=i1
        int i2 = sample2(1,rungs)-1;
        if (i2==i1) {
            continue;
        }
        int rung1 = rungOrderVec[i1];
        int rung2 = rungOrderVec[i2];
        
        // get log-likelihoods and beta values of two chains in the comparison
        double logLike1 = particleVec[rung1].logLike_allData;
        double logLike2 = particleVec[rung2].logLike_allData;
        
        double betaPow1 = particleVec[rung1].betaPow;
        double betaPow2 = particleVec[rung2].betaPow;
        
        // calculate acceptance ratio (still in log space)
        double acceptance = (logLike2*betaPow1 + logLike1*betaPow2) - (logLike1*betaPow1 + logLike2*betaPow2);
        
        // accept or reject move
        double rand1 = runif1();
        if (log(rand1)<acceptance) {
            
            // swap beta values
            double spareBetaPow = particleVec[rung1].betaPow;
            particleVec[rung1].betaPow = particleVec[rung2].betaPow;
            particleVec[rung2].betaPow = spareBetaPow;
            
            // swap rung order
            rungOrderVec[i1] = rung2;
            rungOrderVec[i2] = rung1;
            
            // update acceptance rates
            MC_accept++;
        }
        
    }
    */
}

