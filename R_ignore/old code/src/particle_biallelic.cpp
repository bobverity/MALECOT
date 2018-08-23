
#include <Rcpp.h>
#include "particle_biallelic.h"
#include "misc.h"
#include "probability.h"

using namespace std;

//------------------------------------------------
// particle_biallelic::
// default constructor for particle_biallelic class
particle_biallelic::particle_biallelic() {}

//------------------------------------------------
// particle_biallelic::
// text
void particle_biallelic::initialise(vector< vector<int> > &obsData, vector<bool> &anyHet, Rcpp::List args_model, vector< vector<double> > &lookup_homo, vector< vector<double> > &lookup_het, double betaPow_) {
    
    // pointers to observed data and basic data properties
    obsData_ptr = &obsData;
    n = (*obsData_ptr).size();
    L = (*obsData_ptr)[0].size();
    anyHet_ptr = &anyHet;
    
    // pointers to lookup tables
    lookup_homo_ptr = &lookup_homo;
    lookup_het_ptr = &lookup_het;
    
    // model and MCMC parameters
    parallel = Rcpp_to_bool(args_model["parallel"]);
    precision = Rcpp_to_double(args_model["precision"]);
    precisionSize = 1/precision;
    estimateError = Rcpp_to_bool(args_model["estimateError"]);
    M0 = Rcpp_to_int(args_model["M0"]);
    
    K = Rcpp_to_int(args_model["K"]);
    COI_model = Rcpp_to_int(args_model["COI_model_numeric"]);
    COI_max = Rcpp_to_int(args_model["COI_max"]);
    COI_dispersion = Rcpp_to_double(args_model["COI_dispersion"]);
    betaPow = betaPow_;
    e1 = Rcpp_to_double(args_model["e1"]);
    e2 = Rcpp_to_double(args_model["e2"]);
    e1_max = Rcpp_to_double(args_model["e1_max"]);
    e2_max = Rcpp_to_double(args_model["e2_max"]);
    e1_propSD = Rcpp_to_double(args_model["e1_propSD"]);
    e2_propSD = Rcpp_to_double(args_model["e2_propSD"]);
    p_propSD_raw = Rcpp_to_double(args_model["p_propSD"]);
    
    // initialise e1_propSD and e2_propSD
    e1_propSD = (e1_propSD<0) ? 0.5 : e1_propSD;
    e2_propSD = (e2_propSD<0) ? 0.5 : e2_propSD;
    
    // initialise p_propSD
    p_propSD_raw = (p_propSD_raw<0) ? 0.5 : p_propSD_raw;
    p_propSD = vector< vector<double> >(K, vector<double>(L,p_propSD_raw));
    
    // initialise COI_mean
    COI_mean = vector<double>(K,COI_max);
    COI_mean_shape = vector<double>(K);
    COI_mean_rate = vector<double>(K);
    COI_mean_propSD = vector<double>(K, 0.5);
    
    // grouping
    group = vector<int>(n);
    
    // likelihood
    logLike_old = vector< vector<double> >(n, vector<double>(L));
    logLike_new = vector< vector<double> >(n, vector<double>(L));
    logLike_new_group = vector< vector<double> >(K, vector<double>(L));
    
    // COI and allele frequencies
    m = vector<int>(n,M0);
    p = vector< vector<double> >(K, vector<double>(L,0.5));
    
    // probability vectors and matrices used when updating draws
    logQmatrix_update = vector< vector<double> >(n, vector<double>(K));
    Qmatrix_update = vector< vector<double> >(n, vector<double>(K));
    
    sumLogLike_old_vec = vector<double>(K);
    sumLogLike_new_vec = vector<double>(K);
    p_prop = vector<double>(K);
    COI_mean_prop = vector<double>(K);
    
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
    
    // store acceptance rates
    p_accept = vector< vector<int> >(K, vector<int>(L));
    e1_accept = 0;
    e2_accept = 0;
    disruptive_accept = 0;
    
    // function pointers
    if (precision==0) {
        logProb_genotype_ptr = &particle_biallelic::logProb_genotype_exact;
    } else {
        if (estimateError) {
            logProb_genotype_ptr = &particle_biallelic::logProb_genotype_lookup_varE;
        } else {
            logProb_genotype_ptr = &particle_biallelic::logProb_genotype_lookup;
        }
    }
    if (parallel) {
        update_e_ptr = &particle_biallelic::update_e_parallel;
        update_p_ptr = &particle_biallelic::update_p_parallel;
        update_m_ptr = &particle_biallelic::update_m_parallel;
        update_group_ptr = &particle_biallelic::update_group_parallel;
    } else {
        update_e_ptr = &particle_biallelic::update_e_linear;
        update_e1_ptr = &particle_biallelic::update_e1_linear;
        update_e2_ptr = &particle_biallelic::update_e2_linear;
        update_p_ptr = &particle_biallelic::update_p_linear;
        update_m_ptr = &particle_biallelic::update_m_linear;
        update_group_ptr = &particle_biallelic::update_group_linear;
    }
    
    // initialise with random grouping
    for (int i=0; i<n; i++) {
        group[i] = sample2(1,K)-1;
    }
    
    // calculate initial likelihood
    logLike_allData = 0;
    for (int i=0; i<n; i++) {
        for (int j=0; j<L; j++) {
            logLike_old[i][j] = logProb_genotype((*obsData_ptr)[i][j], p[group[i]][j], m[i], e1, e2);
            logLike_allData += logLike_old[i][j];
        }
    }
    
}

//########################################################################################################
// switch between function types

//------------------------------------------------
// particle_biallelic::
// switch between functions for calculating log-probability of genotypes
double particle_biallelic::logProb_genotype(int S, double p, int m, double e1, double e2) {
    return (this->*logProb_genotype_ptr)(S, p, m, e1, e2);
}

//------------------------------------------------
// particle_biallelic::
// switch between functions for updating e
void particle_biallelic::update_e() {
    (this->*update_e_ptr)();
}

//------------------------------------------------
// particle_biallelic::
// switch between functions for updating e1
void particle_biallelic::update_e1(bool RobbinsMonro, int iteration) {
    (this->*update_e1_ptr)(RobbinsMonro, iteration);
}

//------------------------------------------------
// particle_biallelic::
// switch between functions for updating e2
void particle_biallelic::update_e2(bool RobbinsMonro, int iteration) {
    (this->*update_e2_ptr)(RobbinsMonro, iteration);
}

//------------------------------------------------
// particle_biallelic::
// switch between functions for updating p
void particle_biallelic::update_p(bool RobbinsMonro, int iteration) {
    (this->*update_p_ptr)(RobbinsMonro, iteration);
}

//------------------------------------------------
// particle_biallelic::
// switch between functions for updating m
void particle_biallelic::update_m() {
    (this->*update_m_ptr)();
}

//------------------------------------------------
// particle_biallelic::
// switch between functions for updating group
void particle_biallelic::update_group() {
    (this->*update_group_ptr)();
}

//########################################################################################################
// genotype probabilities

//------------------------------------------------
// particle_biallelic::
// exact log-probabilty of genotypes
double particle_biallelic::logProb_genotype_exact(int S, double p, int m, double e1, double e2) {
    double ret = 0; // return value for missing data
    if (S==1) {
        ret = log( (1.0-e1)*((double)pow(p,m)) + 0.5*e2*(1.0-(double)pow(p,m)-(double)pow(1.0-p,m)) );
    } else if (S==3) {
        ret = log( (1.0-e1)*((double)pow(1.0-p,m)) + 0.5*e2*(1.0-(double)pow(p,m)-(double)pow(1.0-p,m)) );
    } else if (S==2) {
        ret = log( e1*((double)pow(p,m)) + e1*((double)pow(1.0-p,m)) + (1.0-e2)*(1.0-pow(p,m)-pow(1.0-p,m)) );
    }
    ret = (ret < -OVERFLO) ? -OVERFLO : ret;    // catch -inf values
    return ret;
}

//------------------------------------------------
// particle_biallelic::
// log-probabilty of genotypes using lookup tables
double particle_biallelic::logProb_genotype_lookup(int S, double p, int m, double e1, double e2) {
    int p_index = round(p*precisionSize);
    double ret = 0; // return value for missing data
    if (S==1) {
        ret = (*lookup_homo_ptr)[p_index][m-1];
    } else if (S==3) {
        ret = (*lookup_homo_ptr)[precisionSize-p_index][m-1];
    } else if (S==2) {
        ret = (*lookup_het_ptr)[p_index][m-1];
    }
    ret = (ret < -OVERFLO) ? -OVERFLO : ret;    // catch -inf values
    return ret;
}

//------------------------------------------------
// particle_biallelic::
// log-probabilty of genotypes using lookup tables, with error terms specified
double particle_biallelic::logProb_genotype_lookup_varE(int S, double p, int m, double e1, double e2) {
    int p_index = round(p*precisionSize);
    double ret = 0; // return value for missing data
    if (S==1) {
        ret = log( (1.0-e1)*(*lookup_homo_ptr)[p_index][m-1] + 0.5*e2*(*lookup_het_ptr)[p_index][m-1]);
    } else if (S==3) {
        ret = log( (1.0-e1)*(*lookup_homo_ptr)[precisionSize-p_index][m-1] + 0.5*e2*(*lookup_het_ptr)[p_index][m-1]);
    } else if (S==2) {
        ret = log( e1*(*lookup_homo_ptr)[p_index][m-1] + e1*(*lookup_homo_ptr)[precisionSize-p_index][m-1] + (1.0-e2)*(*lookup_het_ptr)[p_index][m-1] );
    }
    ret = (ret < -OVERFLO) ? -OVERFLO : ret;    // catch -inf values
    return ret;
}

//########################################################################################################
// linear update functions

//------------------------------------------------
// particle_biallelic::
// linear update of both e1 and e2. By updating both error terms at the same time we avoid having to recalculate the likelihood for each parameter separately. However, the downside is that we cannot update the proposal standard deviations by Robbins-Monro, because both parameters are accepted and rejected together. This is why we also have update functions for e1 and e2 separately below.
void particle_biallelic::update_e_linear() {
    
    // define sum over old and new likelihood
    double sumLogLike_old = 0;
    double sumLogLike_new = 0;
    
    // propose new values for e1 and e2
    double e1_prop = rnorm1_interval(e1, e1_propSD, 0, e1_max);
    double e2_prop = rnorm1_interval(e2, e2_propSD, 0, e2_max);
    
    // calculate likelihood
    for (int i=0; i<n; i++) {
        for (int j=0; j<L; j++) {
            logLike_new[i][j] = logProb_genotype((*obsData_ptr)[i][j], p[group[i]][j], m[i], e1_prop, e2_prop);
            sumLogLike_old += logLike_old[i][j];
            sumLogLike_new += logLike_new[i][j];
        }
    }
    
    // catch impossible proposed values
    if (sumLogLike_new <= -OVERFLO) {
        return;
    }
    
    // Metropolis step
    if (log(runif_0_1())<betaPow*(sumLogLike_new - sumLogLike_old)) {
        
        // update e1 and e2
        e1 = e1_prop;
        e2 = e2_prop;
        
        // update logLike
        for (int i=0; i<n; i++) {
            for (int j=0; j<L; j++) {
                logLike_old[i][j] = logLike_new[i][j];
            }
        }
        
        // update acceptance rates
        e1_accept++;
        e2_accept++;
        
    }   // end Metropolis step
}

//------------------------------------------------
// particle_biallelic::
// linear update of e1. This function allows optional Robbins-Monro update step of e1 proposal standard deviation.
void particle_biallelic::update_e1_linear(bool RobbinsMonro, int iteration) {
    
    // define sum over old and new likelihood
    double sumLogLike_old = 0;
    double sumLogLike_new = 0;
    
    // propose new value for e1
    double e1_prop = rnorm1_interval(e1, e1_propSD, 0, e1_max);
    
    // calculate likelihood
    for (int i=0; i<n; i++) {
        for (int j=0; j<L; j++) {
            logLike_new[i][j] = logProb_genotype((*obsData_ptr)[i][j], p[group[i]][j], m[i], e1_prop, e2);
            sumLogLike_old += logLike_old[i][j];
            sumLogLike_new += logLike_new[i][j];
        }
    }
    
    // catch impossible proposed values
    if (sumLogLike_new <= -OVERFLO) {
        return;
    }
    
    // Metropolis step
    if (log(runif_0_1())<betaPow*(sumLogLike_new - sumLogLike_old)) {
        
        // update e1
        e1 = e1_prop;
        
        // update logLike
        for (int i=0; i<n; i++) {
            for (int j=0; j<L; j++) {
                logLike_old[i][j] = logLike_new[i][j];
            }
        }
        
        // Robbins-Monro positive update
        if (RobbinsMonro) {
            e1_propSD  += (1-0.23)/sqrt(double(iteration));
        }
        
        // update acceptance rates
        e1_accept++;
        
    } else {
        
        // Robbins-Monro negative update e1
        if (RobbinsMonro) {
            e1_propSD  -= 0.23/sqrt(double(iteration));
            if (e1_propSD < UNDERFLO) {
                e1_propSD = UNDERFLO;
            }
        }
        
    }    // end Metropolis step
    
}

//------------------------------------------------
// particle_biallelic::
// linear update of e2. This function allows optional Robbins-Monro update step of e2 proposal standard deviation.
void particle_biallelic::update_e2_linear(bool RobbinsMonro, int iteration) {
    
    // define sum over old and new likelihood
    double sumLogLike_old = 0;
    double sumLogLike_new = 0;
    
    // propose new value for e2
    double e2_prop = rnorm1_interval(e2, e2_propSD, 0, e2_max);
    
    // calculate likelihood
    for (int i=0; i<n; i++) {
        for (int j=0; j<L; j++) {
            logLike_new[i][j] = logProb_genotype((*obsData_ptr)[i][j], p[group[i]][j], m[i], e1, e2_prop);
            sumLogLike_old += logLike_old[i][j];
            sumLogLike_new += logLike_new[i][j];
        }
    }
    
    // catch impossible proposed values
    if (sumLogLike_new <= -OVERFLO) {
        return;
    }
    
    // Metropolis step
    if (log(runif_0_1())<betaPow*(sumLogLike_new - sumLogLike_old)) {
        
        // update e2
        e2 = e2_prop;
        
        // update logLike
        for (int i=0; i<n; i++) {
            for (int j=0; j<L; j++) {
                logLike_old[i][j] = logLike_new[i][j];
            }
        }
        
        // Robbins-Monro positive update
        if (RobbinsMonro) {
            e2_propSD  += (1-0.23)/sqrt(double(iteration));
        }
        
        // update acceptance rates
        e2_accept++;
        
    } else {
        
        // Robbins-Monro negative update e1
        if (RobbinsMonro) {
            e2_propSD  -= 0.23/sqrt(double(iteration));
            if (e2_propSD < UNDERFLO) {
                e2_propSD = UNDERFLO;
            }
        }
        
    }   // end Metropolis step
    
}

//------------------------------------------------
// particle_biallelic::
// linear update of p
void particle_biallelic::update_p_linear(bool RobbinsMonro, int iteration) {
    
    // loop through loci
    for (int j=0; j<L; j++) {
        
        // clear vectors storing sum over old and new likelihood
        fill(sumLogLike_old_vec.begin(), sumLogLike_old_vec.end(), 0);
        fill(sumLogLike_new_vec.begin(), sumLogLike_new_vec.end(), 0);
        
        // draw p for all demes
        for (int k=0; k<K; k++) {
            p_prop[k] = rnorm1_interval(p[k][j], p_propSD[k][j], 0, 1.0);
        }
        
        // calculate likelihood
        for (int i=0; i<n; i++) {
            logLike_new[i][j] = logProb_genotype((*obsData_ptr)[i][j], p_prop[group[i]], m[i], e1, e2);
            sumLogLike_old_vec[group[i]] += logLike_old[i][j];
            sumLogLike_new_vec[group[i]] += logLike_new[i][j];
        }
        
        // TODO - add lambda prior
        
        // Metropolis step for all K
        for (int k=0; k<K; k++) {
            
            // catch impossible proposed values
            if (sumLogLike_new_vec[k] <= -OVERFLO) {
                continue;
            }
            
            // Metropolis step
            if (log(runif_0_1())<betaPow*(sumLogLike_new_vec[k] - sumLogLike_old_vec[k])) {
                
                // update p
                p[k][j] = p_prop[k];
                
                // update logLike in all individuals from this group
                for (int i=0; i<n; i++) {
                    if (group[i]==k) {
                        logLike_old[i][j] = logLike_new[i][j];
                    }
                }
                
                // Robbins-Monro positive update
                if (RobbinsMonro) {
                    p_propSD[k][j]  += (1-0.23)/sqrt(double(iteration));
                }
                
                // update acceptance rates
                p_accept[k][j]++;
                
            } else {
                
                // Robbins-Monro negative update
                if (RobbinsMonro) {
                    p_propSD[k][j]  -= 0.23/sqrt(double(iteration));
                    if (p_propSD[k][j] < UNDERFLO) {
                        p_propSD[k][j] = UNDERFLO;
                    }
                }
                
            }
        }  // end Metropolis step
        
    }   // end loop through loci
    
}

//------------------------------------------------
// particle_biallelic::
// linear update of m
void particle_biallelic::update_m_linear() {
    
    // define sum over old and new likelihood
    double sumLogLike_old = 0;
    double sumLogLike_new = 0;
    
    // loop through individuals
    for (int i=0; i<n; i++) {
        int thisGroup = group[i];
        
        // propose new m
        int m_prop = rbernoulli1(0.5);
        m_prop = (m_prop==0) ? m[i]-1 : m[i]+1;
        
        // TODO - skip if proposed m impossible
        //if (m_prop==1 && (*anyHet_ptr)[i]) {
        //    continue;
        //}
        
        // calculate likelihood
        if (m_prop>0 && m_prop<=COI_max) {
            
            // reset sum over old and new likelihood
            sumLogLike_old = 0;
            sumLogLike_new = 0;
            
            // calculate likelihood
            for (int j=0; j<L; j++) {
                logLike_new[i][j] = logProb_genotype((*obsData_ptr)[i][j], p[thisGroup][j], m_prop, e1, e2);
                sumLogLike_old += betaPow*logLike_old[i][j];
                sumLogLike_new += betaPow*logLike_new[i][j];
            }
            
            // catch impossible proposed values
            if (sumLogLike_new <= -OVERFLO) {
                continue;
            }
            
            // apply Poisson or negative Binomial prior
            if (COI_model==2) {
                sumLogLike_old += dpois1(m[i]-1, COI_mean[thisGroup]-1);
                sumLogLike_new += dpois1(m_prop-1, COI_mean[thisGroup]-1);
            } else if (COI_model==3) {
                sumLogLike_old += dnbinom1(m[i]-1, COI_mean[thisGroup]-1, COI_dispersion);
                sumLogLike_new += dnbinom1(m_prop-1, COI_mean[thisGroup]-1, COI_dispersion);
            }
            
            // Metropolis step
            if (log(runif_0_1()) < (sumLogLike_new - sumLogLike_old)) {
                m[i] = m_prop;
                for (int j=0; j<L; j++) {
                    logLike_old[i][j] = logLike_new[i][j];
                }
            }
        }
        
    }   // end loop through individuals
    
}

//------------------------------------------------
// particle_biallelic::
// calculate total log-likelihood
// TODO - function pointer?
void particle_biallelic::update_COI_mean(bool RobbinsMonro, int iteration) {
    
    // split method between poisson and negative binomial
    // poisson model
    if (COI_model==2) {
        
        // reset shape and rate vectors
        fill(COI_mean_shape.begin(), COI_mean_shape.end(), 0.25);
        fill(COI_mean_rate.begin(), COI_mean_rate.end(), 0.25);
        
        // update shape and rate vectors to posterior values
        for (int i=0; i<n; i++) {
            int thisGroup = group[i];
            COI_mean_shape[thisGroup] += m[i]-1;
            COI_mean_rate[thisGroup] += 1;
        }
        
        // draw new COI means from conditional posterior
        for (int k=0; k<K; k++) {
            COI_mean[k] = rgamma1(COI_mean_shape[k], COI_mean_rate[k])+1;
        }
        
    }
    // poisson model
    else if (COI_model==3) {
        
        // clear vectors storing sum over old and new likelihood
        fill(sumLogLike_old_vec.begin(), sumLogLike_old_vec.end(), 0);
        fill(sumLogLike_new_vec.begin(), sumLogLike_new_vec.end(), 0);
        
        // draw COI_mean for all demes
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

//------------------------------------------------
// particle_biallelic::
// text
void particle_biallelic::update_group_linear() {
    
    // no change if K==1
    if (K==1) {
        return;
    }
    
    // loop through individuals
    for (int i=0; i<n; i++) {
        int thisGroup = group[i];
        
        // reset logQmatrix_update for this individual
        fill(logQmatrix_update[i].begin(), logQmatrix_update[i].end(), 0);
        
        // calculate likelihood
        for (int k=0; k<K; k++) {
            fill(logLike_new_group[k].begin(), logLike_new_group[k].end(), 0);
            if (k==thisGroup) {  // likelihood already calculated for current group
                for (int j=0; j<L; j++) {
                    logLike_new_group[k][j] = logLike_old[i][j];
                    logQmatrix_update[i][k] += betaPow*logLike_new_group[k][j];
                }
            } else {    // calculate likelihood for group k
                for (int j=0; j<L; j++) {
                    logLike_new_group[k][j] = logProb_genotype((*obsData_ptr)[i][j], p[k][j], m[i], e1, e2);
                    logQmatrix_update[i][k] += betaPow*logLike_new_group[k][j];
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
        int newGroup = sample1(Qmatrix_update[i], 1.0)-1;
        
        // if group has changed update likelihood
        if (newGroup!=thisGroup) {
            for (int j=0; j<L; j++) {
                logLike_old[i][j] = logLike_new_group[newGroup][j];
            }
        }
        
        // update group
        group[i] = newGroup;
    }
    
}

//########################################################################################################
// methods for improving mixing

//------------------------------------------------
// particle_biallelic::
// Propose new sample state by drawing from fairly naive proposal. Start by choosing 3 demes at random (or 2 demes if K=2). We will update alleles and grouping for these demes only. Start by assuming that the conditional posterior allele frequencies can be adequately captured by a Beta distribution that approximates the true distribution. This beta approximation uses shape parameters that match the true distribution in both the mean and variance. Turning to the probability of a sample belonging in each deme, we can integrate over this Beta distribution analytically to give the probabiltiy unconditional on allele frequencies, giving a sort of approximate Gibbs sampler. We remove all samples from the 3 target demes and drop samples into groups one at a time according to these probabilities. Note that samples are allocated into groups only conditional on the samples that have been dropped in already, and so this is not a complete Gibbs sampling step, but rather a sort of naive proposal that hopefully tends to drop samples into decent groupings. Finally, we calculate the likelihood of the new grouping and combine this with forwards and backwards proposal probabilities to arrive at a correct Metropolis-Hastings acceptance ratio for the proposed move.
// Some description of object formats:
// Elements in dp_group are always {0,1,2} (or just {0,1} if K==2), and so do not match up exactly with values in the true MCMC group. To get to final MCMC group values you must use the dp_group_key.
// The object dp_index is a vector that stores all i for which the MCMC group is within dp_group_key. Although this vector is of length n, it only contains meaningful values up to length dp_index_length, which will vary from iteration to iteration.

// NOTE - some improvements are needed for this step to really work well.
//  1) Samples should be dropped into groups in a random order each time, rather than sequentially.
//  2) Efficiency should be improved through lookup tables of the lgamma function.
//  3) Potentially a second (or even several) Gibbs sampling steps should preceed the initial assignment to improve the likelihood of the proposed group, although this will come at the cost of a reduced proposal probability.
//  4) Perform tests to determine whether all chains should be updated, or just the cold chain.

void particle_biallelic::propose_disruptive() {
    
    // limit to K>1
    if (K==1) {
        return;
    }
    
    // limit to cold chain only?
    if (betaPow<1) {
        //return;
    }
    
    // get number of groups that will be resampled
    int dp_K = (K==2) ? 2 : 3;
    
    // sample dp_K groups from 0:(K-1) without replacement
    vector<int> all_K = seq_int(0,K-1);
    dp_group_key[2] = -1; // essentially masks third element of this vector for cases when dp_K==2
    for (int i=0; i<dp_K; i++) {
        int samp1 = sample2(0,all_K.size()-1);
        dp_group_key[i] = all_K[samp1];
        all_K.erase(all_K.begin()+samp1);
    }
    
    // get index of all elements assigned to any of these groups
    int dp_index_length = 0;
    for (int i=0; i<n; i++) {
        if (group[i]==dp_group_key[0]) {
            dp_group_old[i] = 0;
            dp_index[dp_index_length] = i;
            dp_index_length++;
        } else if (group[i]==dp_group_key[1]) {
            dp_group_old[i] = 1;
            dp_index[dp_index_length] = i;
            dp_index_length++;
        } else if (group[i]==dp_group_key[2]) {
            dp_group_old[i] = 2;
            dp_index[dp_index_length] = i;
            dp_index_length++;
        }
    }
    
    // calculate current log-likelihood
    double logLike_old = 0;
    for (int i=0; i<dp_index_length; i++) {
        int i2 = dp_index[i];
        for (int j=0; j<L; j++) {
            logLike_old += betaPow * logProb_genotype((*obsData_ptr)[i2][j], p[dp_group_key[dp_group_old[i2]]][j], m[i2], e1, e2);
        }
    }
    
    // reset beta distribution shape parameters
    for (int k=0; k<dp_K; k++) {
        fill(dp_shape1[k].begin(), dp_shape1[k].end(), 1);
        fill(dp_shape2[k].begin(), dp_shape2[k].end(), 1);
    }
    
    // draw new grouping
    double propProb_forward = 0;
    for (int i=0; i<dp_index_length; i++) {
        int i2 = dp_index[i];
        
        // calculate logProb of data from each group, integrated over p
        fill(dp_logProbVec.begin(), dp_logProbVec.end(), 0);
        for (int j=0; j<L; j++) {
            if ((*obsData_ptr)[i2][j]==1) {
                for (int k=0; k<dp_K; k++) {
                    dp_logProbVec[k] += lgamma(dp_shape1[k][j]+dp_shape2[k][j]) - lgamma(dp_shape1[k][j]+dp_shape2[k][j]+m[i2]) + lgamma(dp_shape1[k][j]+m[i2]) - lgamma(dp_shape1[k][j]);
                }
            } else if ((*obsData_ptr)[i2][j]==3) {
                for (int k=0; k<dp_K; k++) {
                    dp_logProbVec[k] += lgamma(dp_shape1[k][j]+dp_shape2[k][j]) - lgamma(dp_shape1[k][j]+dp_shape2[k][j]+m[i2]) + lgamma(dp_shape2[k][j]+m[i2]) - lgamma(dp_shape2[k][j]);
                }
            } else if ((*obsData_ptr)[i2][j]==2) {
                for (int k=0; k<dp_K; k++) {
                    double const1 = lgamma(dp_shape1[k][j]+dp_shape2[k][j]) - lgamma(dp_shape1[k][j]+dp_shape2[k][j]+m[i2]) + lgamma(dp_shape1[k][j]+m[i2]) - lgamma(dp_shape1[k][j]);
                    double const2 = lgamma(dp_shape1[k][j]+dp_shape2[k][j]) - lgamma(dp_shape1[k][j]+dp_shape2[k][j]+m[i2]) + lgamma(dp_shape2[k][j]+m[i2]) - lgamma(dp_shape2[k][j]);
                    double const3 = 1 - exp(const1) - exp(const2);
                    if (const3 < UNDERFLO) { const3 = UNDERFLO; }
                    dp_logProbVec[k] += log(const3);
                }
            }
        }   // end loop over j
        
        // exponentiate without underflow
        fill(dp_probVec.begin(), dp_probVec.end(), 0);
        double logProbVec_max = max(dp_logProbVec);
        double probVec_sum = 0;
        for (int k=0; k<dp_K; k++) {
            dp_logProbVec[k] -= logProbVec_max;
            dp_probVec[k] = exp(dp_logProbVec[k]);
            probVec_sum += dp_probVec[k];
        }
        double logProbVec_sum = log(probVec_sum);
        for (int k=0; k<dp_K; k++) {
            dp_probVec[k] /= probVec_sum;
            dp_logProbVec[k] -= logProbVec_sum;
        }
        
        // sample new group
        dp_group_new[i2] = sample1(dp_probVec, 1.0)-1;
        
        // calculate probability of new grouping
        propProb_forward += dp_logProbVec[dp_group_new[i2]];
        
        // update shape parameters
        for (int j=0; j<L; j++) {
            if ((*obsData_ptr)[i2][j]==1) {
                dp_shape1[dp_group_new[i2]][j] += betaPow*m[i2];
            } else if ((*obsData_ptr)[i2][j]==3) {
                dp_shape2[dp_group_new[i2]][j] += betaPow*m[i2];
            } else if ((*obsData_ptr)[i2][j]==2) {
                double shape_modifier = 1.5*(m[i2]+2)*(m[i2]+3)/double((m[i2]+2)*(m[i2]+3)-4*m[i2]) - 1.5;
                dp_shape1[dp_group_new[i2]][j] += betaPow*shape_modifier;
                dp_shape2[dp_group_new[i2]][j] += betaPow*shape_modifier;
            }
        }   // end loop over j
        
    }   // end loop over i
    
    // draw new allele frequencies from approximate beta distribution and update proposal probability
    for (int j=0; j<L; j++) {
        for (int k=0; k<dp_K; k++) {
            dp_p[k][j] = rbeta1(dp_shape1[k][j], dp_shape2[k][j]);
            dp_p[k][j] = round(dp_p[k][j]*precisionSize)/double(precisionSize);
            if (dp_p[k][j]==0) { dp_p[k][j] = 1/double(precisionSize); }
            if (dp_p[k][j]==1) { dp_p[k][j] = 1.0 - 1/double(precisionSize); }
            propProb_forward += dbeta1(dp_p[k][j], dp_shape1[k][j], dp_shape2[k][j]);
        }
    }
    
    // calculate new likelihood
    double logLike_new = 0;
    for (int i=0; i<dp_index_length; i++) {
        int i2 = dp_index[i];
        for (int j=0; j<L; j++) {
            logLike_new += betaPow * logProb_genotype((*obsData_ptr)[i2][j], dp_p[dp_group_new[i2]][j], m[i2], e1, e2);
        }
    }
    
    // reset beta distribution shape parameters
    for (int k=0; k<dp_K; k++) {
        fill(dp_shape1[k].begin(), dp_shape1[k].end(), 1);
        fill(dp_shape2[k].begin(), dp_shape2[k].end(), 1);
    }
    
    // calculate backwards probability of old grouping
    double propProb_backward = 0;
    for (int i=0; i<dp_index_length; i++) {
        int i2 = dp_index[i];
        
        // calculate logProb of data from each group, integrated over p
        fill(dp_logProbVec.begin(), dp_logProbVec.end(), 0);
        for (int j=0; j<L; j++) {
            if ((*obsData_ptr)[i2][j]==1) {
                for (int k=0; k<dp_K; k++) {
                    dp_logProbVec[k] += lgamma(dp_shape1[k][j]+dp_shape2[k][j]) - lgamma(dp_shape1[k][j]+dp_shape2[k][j]+m[i2]) + lgamma(dp_shape1[k][j]+m[i2]) - lgamma(dp_shape1[k][j]);
                }
            } else if ((*obsData_ptr)[i2][j]==3) {
                for (int k=0; k<dp_K; k++) {
                    dp_logProbVec[k] += lgamma(dp_shape1[k][j]+dp_shape2[k][j]) - lgamma(dp_shape1[k][j]+dp_shape2[k][j]+m[i2]) + lgamma(dp_shape2[k][j]+m[i2]) - lgamma(dp_shape2[k][j]);
                }
            } else if ((*obsData_ptr)[i2][j]==2) {
                for (int k=0; k<dp_K; k++) {
                    double const1 = lgamma(dp_shape1[k][j]+dp_shape2[k][j]) - lgamma(dp_shape1[k][j]+dp_shape2[k][j]+m[i2]) + lgamma(dp_shape1[k][j]+m[i2]) - lgamma(dp_shape1[k][j]);
                    double const2 = lgamma(dp_shape1[k][j]+dp_shape2[k][j]) - lgamma(dp_shape1[k][j]+dp_shape2[k][j]+m[i2]) + lgamma(dp_shape2[k][j]+m[i2]) - lgamma(dp_shape2[k][j]);
                    double const3 = 1 - exp(const1) - exp(const2);
                    if (const3 < UNDERFLO) { const3 = UNDERFLO; }
                    dp_logProbVec[k] += log(const3);
                }
            }
        }   // end loop over j
        
        // exponentiate without underflow
        fill(dp_probVec.begin(), dp_probVec.end(), 0);
        double logProbVec_max = max(dp_logProbVec);
        double probVec_sum = 0;
        for (int k=0; k<dp_K; k++) {
            dp_logProbVec[k] -= logProbVec_max;
            dp_probVec[k] = exp(dp_logProbVec[k]);
            probVec_sum += dp_probVec[k];
        }
        double logProbVec_sum = log(probVec_sum);
        for (int k=0; k<dp_K; k++) {
            dp_probVec[k] /= probVec_sum;
            dp_logProbVec[k] -= logProbVec_sum;
        }
        
        // calculate bakwards probability of old grouping
        propProb_backward += dp_logProbVec[dp_group_old[i2]];
        
        // update shape parameters
        for (int j=0; j<L; j++) {
            if ((*obsData_ptr)[i2][j]==1) {
                dp_shape1[dp_group_old[i2]][j] += betaPow*m[i2];
            } else if ((*obsData_ptr)[i2][j]==3) {
                dp_shape2[dp_group_old[i2]][j] += betaPow*m[i2];
            } else if ((*obsData_ptr)[i2][j]==2) {
                double shape_modifier = 1.5*(m[i2]+2)*(m[i2]+3)/double((m[i2]+2)*(m[i2]+3)-4*m[i2]) - 1.5;
                dp_shape1[dp_group_old[i2]][j] += betaPow*shape_modifier;
                dp_shape2[dp_group_old[i2]][j] += betaPow*shape_modifier;
            }
        }   // end loop over j
        
    }   // end loop over i
    
    // calculate backward proposal probability of old allele frequences
    for (int j=0; j<L; j++) {
        for (int k=0; k<dp_K; k++) {
            propProb_backward += dbeta1(p[dp_group_key[k]][j], dp_shape1[k][j], dp_shape2[k][j]);
        }
    }
    
    // Metropolis-Hastings step
    double MH_ratio = (propProb_backward + logLike_new) - (propProb_forward + logLike_old);
    if ( log(runif_0_1()) < MH_ratio ) {
        
        // update grouping
        for (int i=0; i<dp_index_length; i++) {
            int i2 = dp_index[i];
            group[i2] = dp_group_key[dp_group_new[i2]];
        }
        
        // update acceptance rates
        disruptive_accept ++;
    }
    
}

//########################################################################################################
// parallel functions

//------------------------------------------------
// particle_biallelic::
// text
void particle_biallelic::update_e_parallel() {
    
    print("   TODO - update e parallel");
    
}

//------------------------------------------------
// particle_biallelic::
// text
void particle_biallelic::update_p_parallel(bool RobbinsMonro, int iteration) {
    
    print("   TODO - update p parallel");
    
}

//------------------------------------------------
// particle_biallelic::
// text
void particle_biallelic::update_m_parallel() {
    
    print("   TODO - update m parallel");
    
}

//------------------------------------------------
// particle_biallelic::
// text
void particle_biallelic::update_group_parallel() {
    
    print("   TODO - update group parallel");
    
}

//########################################################################################################
// misc functions

//------------------------------------------------
// particle_biallelic::
// calculate total log-likelihood
void particle_biallelic::calculateLogLike() {
    logLike_allData = 0;
    for (int i=0; i<n; i++) {
        for (int j=0; j<L; j++) {
            logLike_allData += logLike_old[i][j];
        }
    }
}

