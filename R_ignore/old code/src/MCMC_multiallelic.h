
#pragma once

#include "particle_multiallelic.h"

//------------------------------------------------
// class defining categorical MCMC
class MCMC_multiallelic {
    
public:
    
    // PUBLIC OBJECTS
    
    // data and basic data properties
    std::vector< std::vector<int> > obsData_raw;
    std::vector< std::vector< std::vector<int> > > obsData;
    int n;
    int L;
    std::vector<int> Jl;
    
    // model and MCMC parameters
    int burnin;
    int burnin_actual;
    bool autoConverge;
    int samples;
    int rungs;
    bool MetropolisCouple;
    bool disruptivePropose;
    int chains;
    bool parallel;
    double precision;
    int precisionSize;
    int reportIteration;
    bool flushConsole;
    bool solveLabelSwitching;
    bool estimateError;
    double lambda;
    
    int K;
    int COI_max;
    int COI_model;
    double GTI_pow;
    double p_propSD_raw;
    bool p_RobbinsMonro;
    
    // lookup table for lgamma function
    int lgamma_lookup_dim1;
    int lgamma_lookup_dim2;
    std::vector< std::vector<double> > lgamma_lookup;
    
    // objects over chains and rungs
    std::vector< std::vector<particle_multiallelic> > particleMat;
    std::vector< std::vector<double> > betaPowMat;
    std::vector< std::vector<int> > rungOrderMat;
    
    std::vector<particle_multiallelic> particleVec;
    std::vector<double> betaPowVec;
    std::vector<int> rungOrderVec;
    
    // running estimate of Qmatrix
    std::vector< std::vector<double> > logQmatrix_running;
    std::vector< std::vector<double> > Qmatrix_sum;
    
    // ordering of labels
    std::vector<int> labelOrder;
    std::vector<int> labelOrder_new;
    
    // objects for Hungarian algorithm
    std::vector< std::vector<double> > costMat;
    std::vector<int> bestPerm;
    std::vector<int> bestPermOrder;
    std::vector<int>edgesLeft;
    std::vector<int>edgesRight;
    std::vector<int>blockedLeft;
    std::vector<int>blockedRight;
    
    // objects for storing results
    std::vector< std::vector<double> > burnin_logLike;
    std::vector< std::vector<double> > logLike_store;
    std::vector< std::vector<int> > m_store;
    std::vector< std::vector< std::vector< std::vector<double> > > > p_store;
    std::vector< std::vector<double> > COI_mean_store;
    
    std::vector< std::vector<int> > p_accept;
    int disruptive_accept;
    int MC_accept;
    
    // PUBLIC FUNCTIONS
    
    // constructors
    MCMC_multiallelic(Rcpp::List args_data, Rcpp::List args_model);
    
    // other functions
    void burninMCMC(Rcpp::List args_functions);
    void runMCMC(Rcpp::List args_functions);
    void fixLabels(particle_multiallelic particle);
    void MetropolisCoupling(std::vector<particle_multiallelic> &particleVec, std::vector<int> &rungOrderVec);
    
};

