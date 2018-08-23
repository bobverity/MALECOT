
#pragma once

#include "particle_biallelic.h"

//------------------------------------------------
// class defining categorical MCMC
class MCMC_biallelic {
    
public:
    
    // PUBLIC OBJECTS
    
    // data and basic data properties
    std::vector< std::vector<int> > obsData;
    int n;
    int L;
    std::vector<bool> anyHet;
    
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
    
    int K;
    int COI_max;
    int COI_model;
    double GTI_pow;
    double e1;
    double e2;
    double e1_propSD_raw;
    double e2_propSD_raw;
    bool e1_RobbinsMonro;
    bool e2_RobbinsMonro;
    double p_propSD_raw;
    bool p_RobbinsMonro;
    
    // lookup tables
    std::vector< std::vector<double> > lookup_homo;
    std::vector< std::vector<double> > lookup_het;
    std::vector< std::vector<double> > lookup_log_homo;
    std::vector< std::vector<double> > lookup_log_het;
    
    // objects over chains and rungs
    std::vector< std::vector<particle_biallelic> > particleMat;
    std::vector< std::vector<double> > betaPowMat;
    std::vector< std::vector<int> > rungOrderMat;
    
    std::vector<particle_biallelic> particleVec;
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
    std::vector< std::vector< std::vector<double> > > p_store;
    std::vector<double> e1_store;
    std::vector<double> e2_store;
    std::vector< std::vector<double> > COI_mean_store;
    
    std::vector< std::vector<int> > p_accept;
    int e1_accept;
    int e2_accept;
    int disruptive_accept;
    int MC_accept;
    
    // PUBLIC FUNCTIONS
    
    // constructors
    MCMC_biallelic(Rcpp::List args_data, Rcpp::List args_model);
    
    // other functions
    void burninMCMC(Rcpp::List args_functions);
    void runMCMC(Rcpp::List args_functions);
    void fixLabels(particle_biallelic particle);
    void MetropolisCoupling(std::vector<particle_biallelic> &particleVec, std::vector<int> &rungOrderVec);
    
};

