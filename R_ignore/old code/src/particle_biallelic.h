
#pragma once

//------------------------------------------------
// class defining single particle_biallelic of categorical MCMC
class particle_biallelic {
    
    public:
    
    // PUBLIC OBJECTS
    
    // pointers to observed data and basic data properties
    std::vector< std::vector<int> > * obsData_ptr;
    int n;
    int L;
    std::vector<bool> * anyHet_ptr;
    
    // pointers to lookup tables
    std::vector< std::vector<double> > * lookup_homo_ptr;
    std::vector< std::vector<double> > * lookup_het_ptr;
    
    // model and MCMC parameters
    bool parallel;
    double precision;
    int precisionSize;
    bool estimateError;
    int M0;
    
    int K;
    int COI_model;
    int COI_max;
    double COI_dispersion;
    double betaPow;
    double e1;
    double e2;
    double e1_max;
    double e2_max;
    double e1_propSD;
    double e2_propSD;
    double p_propSD_raw;
    
    std::vector< std::vector<double> > p_propSD;
    
    std::vector<double> COI_mean;
    std::vector<double> COI_mean_shape;
    std::vector<double> COI_mean_rate;
    std::vector<double> COI_mean_propSD;
    
    // grouping
    std::vector<int> group;
    
    // likelihood
    std::vector< std::vector<double> > logLike_old;
    std::vector< std::vector<double> > logLike_new;
    std::vector< std::vector<double> > logLike_new_group;
    double logLike_allData;
    
    // COI and allele frequencies
    std::vector<int> m;
    std::vector< std::vector<double> > p;
    
    // probability vectors and matrices used when updating draws
    std::vector< std::vector<double> > logQmatrix_update;
    std::vector< std::vector<double> > Qmatrix_update;
    
    std::vector<double> sumLogLike_old_vec;
    std::vector<double> sumLogLike_new_vec;
    std::vector<double> p_prop;
    std::vector<double> COI_mean_prop;
    
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
    
    // store acceptance rates
    std::vector< std::vector<int> > p_accept;
    int e1_accept;
    int e2_accept;
    int disruptive_accept;
    
    // function pointers
    double (particle_biallelic::*logProb_genotype_ptr) (int, double, int, double, double);
    void (particle_biallelic::*update_e_ptr) ();
    void (particle_biallelic::*update_e1_ptr) (bool, int);
    void (particle_biallelic::*update_e2_ptr) (bool, int);
    void (particle_biallelic::*update_p_ptr) (bool, int);
    void (particle_biallelic::*update_m_ptr) ();
    void (particle_biallelic::*update_group_ptr) ();
    
    
    // PUBLIC FUNCTIONS
    
    // constructors
    particle_biallelic();
    
    // other functions
    void initialise(std::vector< std::vector<int> > &obsData, std::vector<bool> &anyHet, Rcpp::List args_model, std::vector< std::vector<double> > &lookup_homo, std::vector< std::vector<double> > &lookup_het, double betaPow_);
    
    double logProb_genotype(int S, double p, int m, double e1, double e2);
    void update_e();
    void update_e1(bool RobbinsMonro, int iteration);
    void update_e2(bool RobbinsMonro, int iteration);
    void update_p(bool RobbinsMonro, int iteration);
    void update_m();
    void update_COI_mean(bool RobbinsMonro, int iteration);
    void update_group();
    
    double logProb_genotype_exact(int S, double p, int m, double e1, double e2);
    double logProb_genotype_lookup(int S, double p, int m, double e1, double e2);
    double logProb_genotype_lookup_varE(int S, double p, int m, double e1, double e2);
    
    void update_e_linear();
    void update_e1_linear(bool RobbinsMonro, int iteration);
    void update_e2_linear(bool RobbinsMonro, int iteration);
    void update_p_linear(bool RobbinsMonro, int iteration);
    void update_m_linear();
    void update_group_linear();
    
    void update_e_parallel();
    void update_p_parallel(bool RobbinsMonro, int iteration);
    void update_m_parallel();
    void update_group_parallel();
    
    void propose_disruptive();
    
    void calculateLogLike();
    
};

