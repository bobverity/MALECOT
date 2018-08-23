
#include <Rcpp.h>
#include <math.h>
#include "probability.h"
#include "misc.h"

using namespace std;

// comment this line out to use R default random functions
#define USE_MY_RANDOM

//-- set random seed --
random_device rd;
default_random_engine generator(rd());

//------------------------------------------------
// draw from continuous uniform distribution on interval [0,1)
// NB: execution time found to be ~69% of R::runif() method
double runif_0_1() {
    uniform_real_distribution<double> uniform_0_1(0.0,1.0);
    return uniform_0_1(generator);
}

//------------------------------------------------
// draw from continuous uniform distribution on interval [a,b)
// NB: execution time found to be ~69% of R::runif() method
double runif1(double a, double b) {
    uniform_real_distribution<double> uniform_a_b(a,b);
    return uniform_a_b(generator);
}

//------------------------------------------------
// draw from Bernoulli(p) distribution
// NB: execution time found to be ~34% of R::rbinom() method
bool rbernoulli1(double p) {
    bernoulli_distribution dist_bernoulli(p);
    return dist_bernoulli(generator);
}

//------------------------------------------------
// draw from univariate normal distribution
// NB: execution time found to be ~78% of R::rnorm() method
double rnorm1(double mean, double sd) {
    normal_distribution<double> dist_norm(mean,sd);
    return dist_norm(generator);
}

//------------------------------------------------
// draw from univariate normal distribution and reflect to interval (a,b)
double rnorm1_interval(double mean, double sd, double a, double b) {
    
    // draw raw value relative to a
    double ret = rnorm1(mean, sd) - a;
    
    // reflect off boundries at 0 and (b-a)
    if (ret<0 || ret>(b-a)) {
        // use multiple reflections to bring into range [-(b-a), 2(b-a)]
        while (ret < -(b-a)) {
            ret += 2*(b-a);
        }
        while (ret > 2*(b-a)) {
            ret -= 2*(b-a);
        }
        
        // use one more reflection to bring into range [0, (b-a)]
        if (ret < 0) {
            ret = -ret;
        }
        if (ret > (b-a)) {
            ret = 2*(b-a) - ret;
        }
    }
    
    // no longer relative to a
    ret += a;
    
    // don't let ret equal exactly a or b
    if (ret==a) {
        ret += UNDERFLO;
    } else if (ret==b) {
        ret -= UNDERFLO;
    }
    
    return ret;
}

//------------------------------------------------
// sample single value x that lies between a and b (inclusive) with equal probability. Works on positive or negative values of a or b, and works irrespective of which of a or b is larger.
int sample2(int a, int b) {
    if (a<b) {
        return floor(runif1(a, b+1));
    } else {
        return floor(runif1(b, a+1));
    }
}

//------------------------------------------------
// sample single value from given probability vector (that sums to pSum)
int sample1(vector<double> &p, double pSum) {
    double rand = pSum*runif_0_1();
    double z = 0;
    for (int i=0; i<int(p.size()); i++) {
        z += p[i];
        if (rand<z) {
            return i+1;
        }
    }
    return 0;
}

//------------------------------------------------
// draw from gamma(shape,rate) distribution
// NB: execution time found to be ~22% of R::rgamma() method
double rgamma1(double shape, double rate) {
    gamma_distribution<double> rgamma(shape,1.0/rate);
    double x = rgamma(generator);
    
    // check for zero or infinite values (catches bug present in Visual Studio 2010)
    if (x<UNDERFLO) {
        x = UNDERFLO;
    }
    if (x>OVERFLO) {
        x = OVERFLO;
    }
    
    return x;
}

//------------------------------------------------
// draw from beta(alpha,beta) distribution
// NB: execution time found to be ~50% of cpp ratio of gamma draws method
double rbeta1(double shape1, double shape2) {
    if (shape1==1 && shape2==1) {
        return runif_0_1();
    }
    return R::rbeta(shape1,shape2);
}

//------------------------------------------------
// probability density of beta(shape1,shape2) distribution
double dbeta1(double x, double shape1, double shape2, bool returnLog) {
    return R::dbeta(x, shape1, shape2, returnLog);
}

//------------------------------------------------
// draw from dirichlet distribution using vector of shape parameters. Return vector of values.
vector<double> rdirichlet1(vector<double> &shapeVec) {
    
    int n = shapeVec.size();
    vector<double> ret(n);
    double retSum = 0;
    for (int i=0; i<n; i++) {
        ret[i] = rgamma1(shapeVec[i], 1.0);
        retSum += ret[i];
    }
    double retSumInv = 1.0/retSum;
    for (int i=0; i<n; i++) {
        ret[i] *= retSumInv;
    }
    return(ret);
}

//------------------------------------------------
// draw from dirichlet distribution using bespoke inputs. Outputs are given in x (passed by reference for speed). Shape parameters are equal to alpha+beta, where alpha is an integer vector, and beta is a single double.
void rdirichlet2(std::vector<double> &x, std::vector<int> &alpha, double beta) {
    
    int n = x.size();
    double xSum = 0;
    for (int i=0; i<n; i++) {
        x[i] = rgamma1(alpha[i]+beta, 1.0);
        xSum += x[i];
    }
    double xSumInv = 1.0/xSum;
    for (int i=0; i<n; i++) {
        x[i] *= xSumInv;
        if (x[i] <= UNDERFLO) {
            x[i] = UNDERFLO;
        }
    }
}

//------------------------------------------------
// draw from Poisson(rate) distribution
// NB: execution time found to be generally equal or faster than c++ poisson_distribution method (e.g. faster when rate~=5)
int rpois1(double rate) {
    return R::rpois(rate);
}

//------------------------------------------------
// probability mass of Poisson(rate) distribution
double dpois1(int n, double rate, bool returnLog) {
    return R::dpois(n,rate,returnLog);
}

//------------------------------------------------
// draw from negative binomial distribution with mean lambda and variance gamma*lambda (gamma must be >1)
int rnbinom1(double lambda, double gamma) {
    return R::rnbinom(lambda/(gamma-1), 1/gamma);
}

//------------------------------------------------
// probability mass of negative binomial distribution with mean lambda and variance gamma*lambda (gamma must be >1)
double dnbinom1(int n, double lambda, double gamma, bool returnLog) {
    return R::dnbinom(n, lambda/(gamma-1), 1/gamma, returnLog);
}

