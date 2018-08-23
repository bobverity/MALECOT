
#pragma once

//------------------------------------------------
// class defining single particle_multiallelic of categorical MCMC
class particle_multiallelic {
    
    public:
    
    // PUBLIC OBJECTS
    
    // pointers to observed data and basic data properties
    std::vector< std::vector< std::vector<int> > > * obsData_ptr;
    int n;
    int L;
    std::vector<int> Jl;
    
    // model and MCMC parameters
    bool parallel;
    double precision;
    int precisionSize;
    int M0;
    
    int K;
    int COI_max;
    int COI_model;
    double COI_dispersion;
    double betaPow;
    double lambda;
    
    // lookup table for lgamma function
    int lgamma_lookup_dim1;
    int lgamma_lookup_dim2;
    std::vector< std::vector<double> > * lgamma_lookup_ptr;
    
    // COI mean
    std::vector<double> COI_mean;
    std::vector<double> COI_mean_shape;
    std::vector<double> COI_mean_rate;
    std::vector<double> COI_mean_propSD;
    
    // grouping
    std::vector<int> group;
    
    // likelihood
    std::vector< std::vector<double> > logLike_new_group;
    double logLike_allData;
    
    // COI and allele frequencies
    std::vector<int> m;
    std::vector< std::vector< std::vector<double> > > p;
    
    // sample and population configuration
    std::vector< std::vector< std::vector<int> > > configuration;
    std::vector< std::vector< std::vector<int> > > pop_configuration;
    
    // probability vectors and matrices used when updating draws
    std::vector< std::vector<double> > logQmatrix_update;
    std::vector< std::vector<double> > Qmatrix_update;
    
    std::vector<double> sumLogLike_old_vec;
    std::vector<double> sumLogLike_new_vec;
    std::vector<double> p_prop;
    std::vector<double> COI_mean_prop;
    /*
    // objects for disruptive proposal
    std::vector<int> dp_group_key;
    std::vector<int> dp_index;
    std::vector<int> dp_group_old;
    std::vector<int> dp_group_new;
    std::vector< std::vector<double> > dp_p;
    
    std::vector< std::vector<double> > dp_shape1;
    std::vector< std::vector<double> > dp_shape2;
    std::vector<double> dp_logProbVec;
    std::vector<double> dp_probVec;
    */
    // store acceptance rates
    std::vector< std::vector<int> > p_accept;
    int disruptive_accept;
    
    // function pointers
    void (particle_multiallelic::*update_sample_ptr) ();
    void (particle_multiallelic::*update_p_ptr) ();
    void (particle_multiallelic::*update_group_ptr) ();
    
    
    // PUBLIC FUNCTIONS
    
    // constructors
    particle_multiallelic();
    
    // other functions
    void initialise(std::vector< std::vector< std::vector<int> > > &obsData, std::vector<int> Jl_, Rcpp::List args_model, std::vector< std::vector<double> > &lgamma_lookup, double betaPow_);
    
    void update_sample();
    void update_configuration();
    void update_p();
    void update_group();
    void update_COI_mean(bool RobbinsMonro, int iteration);
    
    double sum_over_configurations(std::vector<int> &S, int m, std::vector<int> &z);
    void draw_configuration(std::vector<int> &S, int m, std::vector<int> &z, std::vector<int> &c, double likelihood_sum);
    
    void update_sample_linear();
    void update_p_linear();
    void update_group_linear();
    
    void update_sample_parallel();
    void update_p_parallel();
    void update_group_parallel();
    
    void propose_disruptive();
    void calculateLogLike();
    
    double lgamma2(int x, int n, double alpha);
    
};

