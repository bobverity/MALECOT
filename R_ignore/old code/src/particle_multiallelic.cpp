
#include <Rcpp.h>
#include "particle_multiallelic.h"
#include "misc.h"
#include "probability.h"

using namespace std;

//------------------------------------------------
// particle_multiallelic::
// default constructor for particle_multiallelic class
particle_multiallelic::particle_multiallelic() {}

//------------------------------------------------
// particle_multiallelic::
// text
void particle_multiallelic::initialise(vector< vector< vector<int> > > &obsData, vector<int> Jl_, Rcpp::List args_model, vector< vector<double> > &lgamma_lookup, double betaPow_) {
    
    // pointers to observed data and basic data properties
    obsData_ptr = &obsData;
    n = (*obsData_ptr).size();
    L = (*obsData_ptr)[0].size();
    Jl = Jl_;
    
    // model and MCMC parameters
    parallel = Rcpp_to_bool(args_model["parallel"]);
    M0 = Rcpp_to_int(args_model["M0"]);
    
    K = Rcpp_to_int(args_model["K"]);
    COI_max = Rcpp_to_int(args_model["COI_max"]);
    COI_model = Rcpp_to_int(args_model["COI_model_numeric"]);
    COI_dispersion = Rcpp_to_double(args_model["COI_dispersion"]);
    betaPow = betaPow_;
    lambda = 1.0;   // TODO - change lambda format to allow multiple values
    
    // lookup table for lgamma function
    lgamma_lookup_ptr = &lgamma_lookup;
    lgamma_lookup_dim1 = (*lgamma_lookup_ptr).size();
    lgamma_lookup_dim2 = (*lgamma_lookup_ptr)[0].size();
    
    // initialise COI_mean
    COI_mean = vector<double>(K,COI_max);
    COI_mean_shape = vector<double>(K);
    COI_mean_rate = vector<double>(K);
    COI_mean_propSD = vector<double>(K, 0.5);
    
    // initialise with random grouping
    group = vector<int>(n);
    for (int i=0; i<n; i++) {
        group[i] = sample2(1,K)-1;
    }
    
    // initialise likelihood
    logLike_new_group = vector< vector<double> >(K, vector<double>(L));
    
    // initialise COI
    m = vector<int>(n,M0);
    
    // initialise allele frequencies and population configuration
    p = vector< vector< vector<double> > >(K, vector< vector<double> >(L));
    pop_configuration = vector< vector< vector<int> > >(K, vector< vector<int> >(L));
    for (int k=0; k<K; k++) {
        for (int j=0; j<L; j++) {
            p[k][j] = vector<double>(Jl[j],1/double(Jl[j]));
            pop_configuration[k][j] = vector<int>(Jl[j]);
        }
    }
    
    // initialise sample configuration and update population configuration
    configuration = vector< vector< vector<int> > >(n, vector< vector<int> >(L));
    for (int i=0; i<n; i++) {
        for (int j=0; j<L; j++) {
            
            // initialise configuration vector
            configuration[i][j] = vector<int>((*obsData_ptr)[i][j].size());
            
            // if single observed allele then must have m[i] copies of this allele
            if ((*obsData_ptr)[i][j].size()==1) {
                if ((*obsData_ptr)[i][j][0] != 0) { // if not missing
                    configuration[i][j][0] = m[i];
                }
            } else {  // otherwise drop into categories at random
                for (int i2=0; i2<m[i]; i2++) {
                    int s = sample2(1,configuration[i][j].size())-1;
                    configuration[i][j][s] ++;
                }
            }
            
            // update population configuration
            for (int h=0; h<int(configuration[i][j].size()); h++) {
                int hap = (*obsData_ptr)[i][j][h];
                if (hap!=0) {
                    pop_configuration[group[i]][j][hap-1] += configuration[i][j][h];
                }
            }
            
        }
    }
    
    // probability vectors and matrices used when updating draws
    logQmatrix_update = vector< vector<double> >(n, vector<double>(K));
    Qmatrix_update = vector< vector<double> >(n, vector<double>(K));
    
    sumLogLike_old_vec = vector<double>(K);
    sumLogLike_new_vec = vector<double>(K);
    p_prop = vector<double>(K);
    COI_mean_prop = vector<double>(K);
    /*
    // objects for disruptive proposal
    dp_group_key = vector<int>(3);
    dp_index = vector<int>(n);
    dp_group_old = vector<int>(n);
    dp_group_new = vector<int>(n);
    dp_p = vector< vector<double> >(3, vector<double>(L));
    
    dp_shape1 = vector< vector<double> >(3, vector<double>(L));
    dp_shape2 = vector< vector<double> >(3, vector<double>(L));
    dp_logProbVec = vector<double>(3);
    dp_probVec = vector<double>(3);
    */
    // store acceptance rates
    p_accept = vector< vector<int> >(K, vector<int>(L));
    disruptive_accept = 0;
    
    // function pointers
    if (parallel) {
        update_p_ptr = &particle_multiallelic::update_p_parallel;
        update_sample_ptr = &particle_multiallelic::update_sample_parallel;
        update_group_ptr = &particle_multiallelic::update_group_parallel;
    } else {
        update_p_ptr = &particle_multiallelic::update_p_linear;
        update_sample_ptr = &particle_multiallelic::update_sample_linear;
        update_group_ptr = &particle_multiallelic::update_group_linear;
    }
    
    // calculate initial likelihood
    calculateLogLike();
}

//########################################################################################################
// switch between function types

//------------------------------------------------
// particle_multiallelic::
// switch between functions for updating m
void particle_multiallelic::update_sample() {
    (this->*update_sample_ptr)();
}

//------------------------------------------------
// particle_multiallelic::
// switch between functions for updating p
void particle_multiallelic::update_p() {
    (this->*update_p_ptr)();
}

//------------------------------------------------
// particle_multiallelic::
// switch between functions for updating group
void particle_multiallelic::update_group() {
    (this->*update_group_ptr)();
}

//########################################################################################################
// dealing with configurations

//------------------------------------------------
// particle_multiallelic::
// obtain marginal likelihood of observed data at a particular locus by summing over all possible configurations. Credit goes to to Inna Gerlovina for the neat way of looping through multinomial combinations.
double particle_multiallelic::sum_over_configurations(vector<int> &S, int m, vector<int> &z) {
    
    // skip over missing data
    if (S[0] == 0) {
        return 1;
    }
    
    // get sum over pop configuration
    int z_sum = sum(z);
    
    // special case if single allele observed at this locus - only one configuration possible
    int S_size = S.size();
    if (S_size == 1) {
        double logLike = lgamma2(z_sum, z.size(), lambda) - lgamma2(z_sum + m, z.size(), lambda) + lgamma2(z[S[0]-1] + m, 1 ,lambda) - lgamma2(z[S[0]-1], 1, lambda);
        return exp(logLike);
    }
    
    // create objects for looping through configurations
    int vs = S_size - 1;
    int vmax = m - vs;
    int no = vmax;
    vector<int> v(vs,1);
    double const1 = lgamma2(m+1, 0, 1) + lgamma2(z_sum, z.size(), lambda) - lgamma2(z_sum + m, z.size(), lambda);
    for (int k=0; k<=vs; k++) {
        const1 -= lgamma2(z[S[k]-1], 1, lambda);
    }
    
    // likelihood of first configuration
    double sp = const1;
    for (int k=0; k<vs; k++) {
        sp += lgamma2(z[S[k]-1] + v[k], 1, lambda) - lgamma2(v[k] + 1, 0, 1);
    }
    sp += lgamma2(z[S[vs]-1] + no, 1, lambda) - lgamma2(no + 1, 0, 1);
    double likelihood = exp(sp);
    
    // loop through all remaining configurations
    int i = 0;
    while (v[vs-1] < vmax) {
        
        // update configuration
        if (v[0] < vmax && no > 1) {
            v[0]++;
            no--;
            i = 0;
        } else {
            if (no <= 1) i++;
            while (v[i] == vmax && i < vs) {
                i++;
            }
            v[i]++;
            no = m - i;
            for (int k=i; k<vs; k++) {
                no -= v[k];
            }
            if (no < 1) continue;
            for (int k=0; k<i; k++) {
                v[k] = 1;
            }
        }
        
        // update likelihood
        sp = const1;
        for (int k=0; k<vs; k++) {
            sp += lgamma2(z[S[k]-1] + v[k], 1, lambda) - lgamma2(v[k] + 1, 0, 1);
        }
        sp += lgamma2(z[S[vs]-1] + no, 1, lambda) - lgamma2(no + 1, 0, 1);
        likelihood += exp(sp);
        
    }   // end while loop over configurations
    
    return(likelihood);
}

//------------------------------------------------
// particle_multiallelic::
// draw from all possible configurations in proportion to how much they contribute to the overall marginal likelihood of observed data at a particular locus.
void particle_multiallelic::draw_configuration(vector<int> &S, int m, vector<int> &z, vector<int> &c, double likelihood_sum) {
    
    // skip over missing data
    if (S[0] == 0) {
        return;
    }
    
    // special case if single allele observed at this locus - only one configuration possible
    int S_size = S.size();
    if (S_size == 1) {
        c = vector<int>(1,m);
        return;
    }
    
    // get sum over pop configuration
    int z_sum = sum(z);
    
    // create objects for looping through configurations
    int vs = S_size - 1;
    int vmax = m - vs;
    int no = vmax;
    vector<int> v(vs,1);
    double const1 = lgamma2(m+1, 0, 1) + lgamma2(z_sum, z.size(), lambda) - lgamma2(z_sum + m, z.size(), lambda);
    for (int k=0; k<=vs; k++) {
        const1 -= lgamma2(z[S[k]-1], 1, lambda);
    }
    
    // draw uniform value to dictate when to exit loop
    double rand1 = runif1(0,likelihood_sum);
    
    // likelihood of first configuration
    double sp = const1;
    for (int k=0; k<vs; k++) {
        sp += lgamma2(z[S[k]-1] + v[k], 1, lambda) - lgamma2(v[k]+1, 0, 1);
    }
    sp += lgamma2(z[S[vs]-1] + no, 1, lambda) - lgamma2(no+1, 0, 1);
    double likelihood = exp(sp);
    
    // loop through all remaining configurations
    int i = 0;
    while (v[vs-1] < vmax) {
        
        // update configuration
        if (v[0] < vmax && no > 1) {
            v[0]++;
            no--;
            i = 0;
        } else {
            if (no <= 1) i++;
            while (v[i] == vmax && i < vs) {
                i++;
            }
            v[i]++;
            no = m - i;
            for (int k=i; k<vs; k++) {
                no -= v[k];
            }
            if (no < 1) continue;
            for (int k=0; k<i; k++) {
                v[k] = 1;
            }
        }
        
        // update likelihood
        sp = const1;
        for (int k=0; k<vs; k++) {
            sp += lgamma2(z[S[k]-1] + v[k], 1, lambda) - lgamma2(v[k]+1, 0, 1);
        }
        sp += lgamma2(z[S[vs]-1] + no, 1, lambda) - lgamma2(no+1, 0, 1);
        likelihood += exp(sp);
        
        // exit if likelihood exceeds rand1
        if (likelihood > rand1) {
            break;
        }
        
    }   // end while loop over configurations
    
    // update configuration
    v.push_back(no);
    c = v;
}

//------------------------------------------------
// particle_multiallelic::
// linear update of configuration
void particle_multiallelic::update_configuration() {
    /*
     // loop through individuals
     for (int i=0; i<n; i++) {
     int thisGroup = group[i];
     
     // skip ahead if COI=1
     if (m[i]==1) {
     continue;
     }
     
     // loop through loci
     for (int j=0; j<L; j++) {
     
     // skip ahead if single observed allele at this locus (also skips over missing data)
     int c = configuration[i][j].size();
     if (c==1) {
     continue;
     }
     
     // draw two different alleles to change - will decrement s1 and increment s2
     int s1 = sample2(1,c) - 1;
     int s2 = sample2(1,c-1) - 1;
     s2 = (s2<s1) ? s2 : s2+1;
     
     // skip ahead if count will be dropped below 1
     if (configuration[i][j][s1]==1) {
     continue;
     }
     
     // calculate ratio of likelihoods
     int n1 = configuration[i][j][s1];
     int n2 = configuration[i][j][s2];
     int hap1 = (*obsData_ptr)[i][j][s1];
     int hap2 = (*obsData_ptr)[i][j][s2];
     double p1 = p[thisGroup][j][hap1-1];
     double p2 = p[thisGroup][j][hap2-1];
     double MH_ratio = n1/double(n2+1) * p2/p1;
     
     // Metropolis-Hastings step
     if (runif_0_1()<MH_ratio) {
     configuration[i][j][s1] --;
     configuration[i][j][s2] ++;
     
     alleleCount[thisGroup][j][hap1-1] --;
     alleleCount[thisGroup][j][hap2-1] ++;
     }
     
     }
     }
     */
}

//########################################################################################################
// linear update functions

//------------------------------------------------
// particle_multiallelic::
// linear update of m
void particle_multiallelic::update_sample_linear() {
    
    // loop through samples
    vector<double> like_new_allLoci(L);
    for (int i=0; i<n; i++) {
        int thisGroup = group[i];
        
        // propose new m
        int m_prop = rbernoulli1(0.5);
        m_prop = (m_prop==0) ? m[i]-1 : m[i]+1;
        
        int tmp1 = rbernoulli1(0.5);
        if (tmp1==0) {
            m_prop = m[i];
        }
        
        // skip if proposed m impossible based on COI limits
        if (m_prop==0 || m_prop>COI_max) {
            continue;
        }
        
        // skip if proposed m less than number of observed alleles at any locus
        bool skip = false;
        for (int j=0; j<L; j++) {
            if (m_prop < (*obsData_ptr)[i][j].size()) {
                skip = true;
                break;
            }
        }
        if (skip) continue;
        
        // subtract this sample configuration from population configuration
        for (int j=0; j<L; j++) {
            for (int h=0; h<int(configuration[i][j].size()); h++) {
                int hap = (*obsData_ptr)[i][j][h];
                if (hap!=0) {
                    pop_configuration[thisGroup][j][hap-1] -= configuration[i][j][h];
                }
            }
        }
        
        // calculate old and new likelihoods over all loci
        double logLike_old = 0;
        double logLike_new = 0;
        for (int j=0; j<L; j++) {
            
            // calculate old likelihood of this locus by summing over all possible configurations
            logLike_old += log( sum_over_configurations((*obsData_ptr)[i][j], m[i], pop_configuration[thisGroup][j]) );
            
            // same for new likelihood, but store raw likelihood for each locus
            like_new_allLoci[j] = sum_over_configurations((*obsData_ptr)[i][j], m_prop, pop_configuration[thisGroup][j]);
            logLike_new += log(like_new_allLoci[j]);
            
        }   // end loop over loci
        
        // catch impossible proposed values
        if (logLike_new <= -OVERFLO) {
            logLike_new = -OVERFLO;
        }
        
        // apply Poisson or negative Binomial prior
        if (COI_model==2) {
            logLike_old += dpois1(m[i]-1, COI_mean[thisGroup]-1);
            logLike_new += dpois1(m_prop-1, COI_mean[thisGroup]-1);
        } else if (COI_model==3) {
            logLike_old += dnbinom1(m[i]-1, COI_mean[thisGroup]-1, COI_dispersion);
            logLike_new += dnbinom1(m_prop-1, COI_mean[thisGroup]-1, COI_dispersion);
        }
        
        // Metropolis step
        if (log(runif_0_1()) < (logLike_new - logLike_old)) {
            m[i] = m_prop;
        }
        
        // draw new configuration at all loci
        for (int j=0; j<L; j++) {
            draw_configuration((*obsData_ptr)[i][j], m[i], pop_configuration[thisGroup][j], configuration[i][j], like_new_allLoci[j]);
        }
        
        // add new sample configuration to population configuration
        for (int j=0; j<L; j++) {
            for (int h=0; h<int(configuration[i][j].size()); h++) {
                int hap = (*obsData_ptr)[i][j][h];
                if (hap!=0) {
                    pop_configuration[thisGroup][j][hap-1] += configuration[i][j][h];
                }
            }
        }
        
    }   // end loop over samples
    
}

//------------------------------------------------
// particle_multiallelic::
// linear update of p
void particle_multiallelic::update_p_linear() {
    
    // loop through demes and loci
    for (int k=0; k<K; k++) {
        for (int j=0; j<L; j++) {
            
            // draw p from conditional posterior distribution
            rdirichlet2(p[k][j], pop_configuration[k][j], lambda);
        }
    }
}

//------------------------------------------------
// particle_multiallelic::
// linear update of group
void particle_multiallelic::update_group_linear() {
    
    // no change if K==1
    if (K==1) {
        return;
    }
    
    // loop through samples
    for (int i=0; i<n; i++) {
        int thisGroup = group[i];
        
        // reset logQmatrix_update for this individual
        fill(logQmatrix_update[i].begin(), logQmatrix_update[i].end(), 0);
        
        // calculate likelihood
        for (int k=0; k<K; k++) {
            for (int j=0; j<L; j++) {
                for (int h=0; h<int(configuration[i][j].size()); h++) {
                    int hap = (*obsData_ptr)[i][j][h];
                    if (hap!=0) {   // if not missing
                        logQmatrix_update[i][k] += betaPow*configuration[i][j][h]*log(p[k][j][hap-1]);
                    }
                }
            }
            
            // add prior information
            if (COI_model==2) {
                logQmatrix_update[i][k] += dpois1(m[i]-1, COI_mean[k]-1);
            } else if (COI_model==3) {
                logQmatrix_update[i][k] += dnbinom1(m[i]-1, COI_mean[k]-1, COI_dispersion);
            }
            
            // limit small (-inf) values
            if (logQmatrix_update[i][k] <= -OVERFLO) {
                logQmatrix_update[i][k] = -OVERFLO;
            }
        }
        
        // exponentiate without underflow
        double group_logProbVec_max = max(logQmatrix_update[i]);
        double group_probVec_sum = 0;
        for (int k=0; k<K; k++) {
            logQmatrix_update[i][k] -= group_logProbVec_max;
            Qmatrix_update[i][k] = exp(logQmatrix_update[i][k]);
            group_probVec_sum += Qmatrix_update[i][k];
        }
        for (int k=0; k<K; k++) {
            Qmatrix_update[i][k] /= group_probVec_sum;
        }
        
        // sample new group
        int newGroup = sample1(Qmatrix_update[i], 1.0) - 1;
        
        // if group has changed update allele counts
        if (newGroup!=thisGroup) {
            for (int j=0; j<L; j++) {
                for (int h=0; h<int(configuration[i][j].size()); h++) {
                    int hap = (*obsData_ptr)[i][j][h];
                    pop_configuration[thisGroup][j][hap-1] -= configuration[i][j][h];
                    pop_configuration[newGroup][j][hap-1] += configuration[i][j][h];
                }
            }
        }
        
        // update group
        group[i] = newGroup;
    }
    
}

//------------------------------------------------
// particle_multiallelic::
// calculate total log-likelihood
// TODO - function pointer?
void particle_multiallelic::update_COI_mean(bool RobbinsMonro, int iteration) {
    
    // split method between poisson and negative binomial
    if (COI_model==2) { // poisson model
        
        // reset shape and rate vectors
        fill(COI_mean_shape.begin(), COI_mean_shape.end(), 0.25);
        fill(COI_mean_rate.begin(), COI_mean_rate.end(), 0.25);
        
        // update shape and rate vectors to posterior values
        for (int i=0; i<n; i++) {
            int thisGroup = group[i];
            COI_mean_shape[thisGroup] += m[i]-1;
            COI_mean_rate[thisGroup] += 1;
        }
        
        // Gibbs sample new COI means from conditional posterior
        for (int k=0; k<K; k++) {
            COI_mean[k] = rgamma1(COI_mean_shape[k], COI_mean_rate[k])+1;
        }
        
    } else if (COI_model==3) { // negative binomial model
        
        // clear vectors storing sum over old and new likelihood
        fill(sumLogLike_old_vec.begin(), sumLogLike_old_vec.end(), 0);
        fill(sumLogLike_new_vec.begin(), sumLogLike_new_vec.end(), 0);
        
        // propose new COI_mean for all demes
        for (int k=0; k<K; k++) {
            COI_mean_prop[k] = rnorm1_interval(COI_mean[k], COI_mean_propSD[k], 1, COI_max);
        }
        
        // calculate likelihood
        for (int i=0; i<n; i++) {
            sumLogLike_old_vec[group[i]] += dnbinom1(m[i]-1, COI_mean[group[i]]-1, COI_dispersion);
            sumLogLike_new_vec[group[i]] += dnbinom1(m[i]-1, COI_mean_prop[group[i]]-1, COI_dispersion);
        }
        
        // Metropolis step for all demes
        for (int k=0; k<K; k++) {
            
            // Metropolis step
            if (log(runif_0_1())<(sumLogLike_new_vec[k] - sumLogLike_old_vec[k])) {
                COI_mean[k] = COI_mean_prop[k];
                
                // Robbins-Monro positive update
                if (RobbinsMonro) {
                    COI_mean_propSD[k] += (1-0.23)/sqrt(double(iteration));
                }
                
            } else {
                
                // Robbins-Monro negative update
                if (RobbinsMonro) {
                    COI_mean_propSD[k] -= 0.23/sqrt(double(iteration));
                    if (COI_mean_propSD[k] < UNDERFLO) {
                        COI_mean_propSD[k] = UNDERFLO;
                    }
                }
                
            }
        }  // end Metropolis step
        
    }
    
}

//########################################################################################################
// parallel functions

//------------------------------------------------
// particle_multiallelic::
// text
void particle_multiallelic::update_sample_parallel() {
    
    print("   TODO - update m parallel");
    
}

//------------------------------------------------
// particle_multiallelic::
// text
void particle_multiallelic::update_p_parallel() {
    
    print("   TODO - update p parallel");
    
}

//------------------------------------------------
// particle_multiallelic::
// text
void particle_multiallelic::update_group_parallel() {
    
    print("   TODO - update group parallel");
    
}

//########################################################################################################
// misc functions

//------------------------------------------------
// particle_multiallelic::
// calculate total log-likelihood
void particle_multiallelic::calculateLogLike() {
    
    logLike_allData = 0;
    for (int i=0; i<n; i++) {
        for (int j=0; j<L; j++) {
            logLike_allData += lgamma(m[i]+1);
            for (int h=0; h<int(configuration[i][j].size()); h++) {
                int hap = (*obsData_ptr)[i][j][h];
                if (hap!=0) {
                    logLike_allData += configuration[i][j][h]*log(p[group[i]][j][hap-1]) - lgamma(configuration[i][j][h]+1);
                }
            }
        }
    }
    
}

//------------------------------------------------
// particle_multiallelic::
// TODO
double particle_multiallelic::lgamma2(int x, int n, double alpha) {
    return (*lgamma_lookup_ptr)[n][x];
}
