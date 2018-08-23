
#include <Rcpp.h>
#include <chrono>
#include "misc.h"
#include "probability.h"
#include "MCMC_biallelic.h"
#include "particle_biallelic.h"
#include "MCMC_multiallelic.h"
#include "Hungarian.h"

using namespace std;

//------------------------------------------------
// TODO - text

// [[Rcpp::export]]
Rcpp::List runMCMC_biallelic_cpp(Rcpp::List args_data, Rcpp::List args_model, Rcpp::List args_functions) {
    
    // extract some model parameters
    int K = Rcpp_to_int(args_model["K"]);
    
    // start program
    print("Running MCMC with K =", K);
    
    // start timer
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    
    // create pointer to entire MCMC object
    Rcpp::XPtr<MCMC_biallelic> mainMCMC(new MCMC_biallelic{args_data, args_model}, true);
    
    // run MCMC
    mainMCMC->burninMCMC(args_functions);
    mainMCMC->runMCMC(args_functions);
    
    // end timer
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
    print("   MCMC completed in", time_span.count(), "seconds");
    
    // create ret object
    Rcpp::List ret;
    ret.push_back(Rcpp::wrap( mainMCMC->burnin_actual ));
    ret.push_back(Rcpp::wrap( mainMCMC->burnin_logLike ));
    ret.push_back(Rcpp::wrap( mainMCMC->logLike_store ));
    ret.push_back(Rcpp::wrap( mainMCMC->Qmatrix_sum ));
    ret.push_back(Rcpp::wrap( mainMCMC->m_store ));
    ret.push_back(Rcpp::wrap( mainMCMC->p_store ));
    ret.push_back(Rcpp::wrap( mainMCMC->COI_mean_store ));
    ret.push_back(Rcpp::wrap( mainMCMC->e1_store ));
    ret.push_back(Rcpp::wrap( mainMCMC->e2_store ));
    ret.push_back(Rcpp::wrap( mainMCMC->p_accept ));
    ret.push_back(Rcpp::wrap( mainMCMC->e1_accept ));
    ret.push_back(Rcpp::wrap( mainMCMC->e2_accept ));
    ret.push_back(Rcpp::wrap( mainMCMC->disruptive_accept ));
    ret.push_back(Rcpp::wrap( mainMCMC->MC_accept ));
    ret.push_back(Rcpp::wrap( time_span.count() ));
    //ret.push_back(Rcpp::wrap( mainMCMC ));
    
    Rcpp::StringVector ret_names;
    ret_names.push_back("burnin_actual");
    ret_names.push_back("logLike_burnin");
    ret_names.push_back("logLike_store");
    ret_names.push_back("Qmatrix_sum");
    ret_names.push_back("m_store");
    ret_names.push_back("p_store");
    ret_names.push_back("COI_mean_store");
    ret_names.push_back("e1_store");
    ret_names.push_back("e2_store");
    ret_names.push_back("p_accept");
    ret_names.push_back("e1_accept");
    ret_names.push_back("e2_accept");
    ret_names.push_back("disruptive_accept");
    ret_names.push_back("MC_accept");
    ret_names.push_back("runTime");
    //ret_names.push_back("MCMC_ptr");
    
    ret.names() = ret_names;
    return ret;
}

//------------------------------------------------
// TODO - text

// [[Rcpp::export]]
Rcpp::List runMCMC_multiallelic_cpp(Rcpp::List args_data, Rcpp::List args_model, Rcpp::List args_functions) {
    
    // extract some model parameters
    int K = Rcpp_to_int(args_model["K"]);
    
    // start program
    print("Running MCMC with K =", K);
    
    // start timer
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    
    // create pointer to entire MCMC object
    Rcpp::XPtr<MCMC_multiallelic> mainMCMC(new MCMC_multiallelic{args_data, args_model}, true);
    
    // run MCMC
    mainMCMC->burninMCMC(args_functions);
    mainMCMC->runMCMC(args_functions);
    
    // end timer
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
    print("   MCMC completed in", time_span.count(), "seconds");
    
    // create ret object
    Rcpp::List ret;
    ret.push_back(Rcpp::wrap( mainMCMC->burnin_actual ));
    ret.push_back(Rcpp::wrap( mainMCMC->burnin_logLike ));
    ret.push_back(Rcpp::wrap( mainMCMC->logLike_store ));
    ret.push_back(Rcpp::wrap( mainMCMC->Qmatrix_sum ));
    ret.push_back(Rcpp::wrap( mainMCMC->m_store ));
    ret.push_back(Rcpp::wrap( mainMCMC->p_store ));
    ret.push_back(Rcpp::wrap( mainMCMC->COI_mean_store ));
    ret.push_back(Rcpp::wrap( mainMCMC->p_accept ));
    ret.push_back(Rcpp::wrap( mainMCMC->disruptive_accept ));
    ret.push_back(Rcpp::wrap( mainMCMC->MC_accept ));
    ret.push_back(Rcpp::wrap( time_span.count() ));
    //ret.push_back(Rcpp::wrap( mainMCMC ));
    
    Rcpp::StringVector ret_names;
    ret_names.push_back("burnin_actual");
    ret_names.push_back("logLike_burnin");
    ret_names.push_back("logLike_store");
    ret_names.push_back("Qmatrix_sum");
    ret_names.push_back("m_store");
    ret_names.push_back("p_store");
    ret_names.push_back("COI_mean_store");
    ret_names.push_back("p_accept");
    ret_names.push_back("disruptive_accept");
    ret_names.push_back("MC_accept");
    ret_names.push_back("runTime");
    //ret_names.push_back("MCMC_ptr");
    
    ret.names() = ret_names;
    return ret;
}

//------------------------------------------------
// TODO - text

// [[Rcpp::export]]
Rcpp::List continueMCMC_cpp(Rcpp::List args_params) {
    
    // start program
    print("Continuing MCMC");
    
    // start timer
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    
    // define MCMC object from MCMC_ptr
    //Rcpp::XPtr<MCMC_biallelic> mainMCMC = Rcpp::as<Rcpp::XPtr<MCMC_biallelic> > (args_params["MCMC_ptr"]);
    
    //mainMCMC->particleVec[0].x.push_back(1);
    
    //printVector(mainMCMC->particleVec[0].x);
    
    // end timer
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
    print("   program completed in", time_span.count(), "seconds");
    
    return Rcpp::List::create(Rcpp::Named("foo")=1);
}

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::NumericVector foobar_linear(Rcpp::NumericVector x) {
    
    // allocate the output vector
    Rcpp::NumericVector output(x.size());
    
    // do something
    for (int i=0; i<int(x.size()); i++) {
        for (int j=0; j<1e3; j++) {
            output[i] += rnorm1(0,1);
        }
    }
    
    // return the new matrix
    return output;
}

//------------------------------------------------
// TODO - text

// [[Rcpp::export]]
Rcpp::List fixLabels_cpp(Rcpp::List args_model) {
    
    // extract cost matrix
    vector< vector<double> > costMat = Rcpp_to_mat_double(args_model["costMat"]);
    
    // other objects needed for Hungarian algorithm
    int K = costMat.size();
    vector<int> bestPerm(K);
    vector<int> edgesLeft(K);
    vector<int> edgesRight(K);
    vector<int> blockedLeft(K);
    vector<int> blockedRight(K);
    
    // find best permutation of current labels using Hungarian algorithm
    bestPerm = hungarian(costMat, edgesLeft, edgesRight, blockedLeft, blockedRight);
    
    return Rcpp::List::create(Rcpp::Named("bestPerm")=bestPerm);
}


