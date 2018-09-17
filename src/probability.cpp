
#include <Rcpp.h>
#include <math.h>

#include "probability.h"
#include "misc_v1.h"

using namespace std;

//------------------------------------------------
// draw from continuous uniform distribution on interval [0,1)
double runif_0_1() {
  return R::runif(0,1);
}

//------------------------------------------------
// draw from continuous uniform distribution on interval [a,b)
double runif1(double a, double b) {
  return R::runif(a,b);
}

//------------------------------------------------
// draw from Bernoulli(p) distribution
bool rbernoulli1(double p) {
  return R::rbinom(1, p);
}

//------------------------------------------------
// draw from geometric(p) distribution with mean (1-p)/p
int rgeom1(double p) {
  return R::rgeom(p);
}

//------------------------------------------------
// draw from univariate normal distribution
double rnorm1(double mean, double sd) {
  return R::rnorm(mean, sd);
}

//------------------------------------------------
// density of univariate normal distribution
double dnorm1(double x, double mean, double sd, bool return_log) {
  return R::dnorm(x, mean, sd, return_log);
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
// draw from multivariate logit-normal distribution, given mean and standard 
// deviation on the logit scale. Produces one more draw than the number of
// meanlog values. Assumes same standard deviation in all dimensions and zero
// correlation
vector<double> rmlogitnorm1(const vector<double> &meanlog, double sdlog) {
  
  int n = int(meanlog.size());
  vector<double> ret(n+1);
  double tmp1 = 0;
  for (int i=0; i<n; ++i) {
    ret[i] = exp(rnorm1(meanlog[i], sdlog));
    tmp1 += ret[i];
  }
  double tmp2 = 1.0/(1.0 + tmp1);
  for (int i=0; i<n; ++i) {
    ret[i] *= tmp2;
  }
  ret[n] = tmp2;
  
  return ret;
}

//------------------------------------------------
// probability density of rmlogitnorm1 distribution
double dmlogitnorm1(const std::vector<double> &x, const vector<double> &meanlog, double sdlog, bool return_log) {
  
  int n = int(meanlog.size());
  double ret = -log(x[n]);
  for (int i=0; i<n; ++i) {
    ret += dnorm1(log(x[i]/x[n]), meanlog[i], sdlog, true) - log(x[i]);
  }
  if (!return_log) {
    return exp(ret);
  }
  return ret;
}

//------------------------------------------------
// equivalent to rmlogitnorm2, except proportions p are passed in and 
// transformed to obtain the meanlog values. Note that p should sum to 1, and 
// therefore should be one element longer than the meanlog argument in 
// rmlogitnorm1. If sdlog is zero then draws will equal p, although generally p 
// does not equal the mean of the distribution.
std::vector<double> rmlogitnorm2(const std::vector<double> &p, double sdlog) {
  
  int n = int(p.size()) - 1;
  vector<double> ret(n+1);
  double tmp1 = 0;
  for (int i=0; i<n; ++i) {
    ret[i] = exp(rnorm1(log(p[i]/p[n]), sdlog));
    tmp1 += ret[i];
  }
  double tmp2 = 1.0/(1.0 + tmp1);
  for (int i=0; i<n; ++i) {
    ret[i] *= tmp2;
  }
  ret[n] = tmp2;
  
  return ret;
}

//------------------------------------------------
// probability density of rmlogitnorm1 distribution
double dmlogitnorm2(const std::vector<double> &x, const vector<double> &p, double sdlog, bool return_log) {
  
  int n = int(p.size()) - 1;
  double ret = -log(x[n]);
  for (int i=0; i<n; ++i) {
    ret += dnorm1(log(x[i]/x[n]), log(p[i]/p[n]), sdlog, true) - log(x[i]);
  }
  if (!return_log) {
    return exp(ret);
  }
  return ret;
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
// sample single value x that lies between a and b (inclusive) with equal 
// probability. Works on positive or negative values of a or b, and works 
// irrespective of which of a or b is larger.
int sample2(int a, int b) {
  if (a<b) {
    return floor(runif1(a, b+1));
  } else {
    return floor(runif1(b, a+1));
  }
}

//------------------------------------------------
// sample a given number of values from a vector without replacement (templated
// for different data types). Note, this function re-arranges the original
// vector (passed in by reference), and the result is stored in the first n
// elements.
// sample3
// DEFINED IN HEADER

//------------------------------------------------
// draw from gamma(shape,rate) distribution
double rgamma1(double shape, double rate) {
  double x = R::rgamma(shape, 1/rate);
  
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
double rbeta1(double shape1, double shape2) {
  if (shape1==1 && shape2==1) {
    return runif_0_1();
  }
  return R::rbeta(shape1, shape2);
}

//------------------------------------------------
// probability density of beta(shape1,shape2) distribution
double dbeta1(double x, double shape1, double shape2, bool return_log) {
  return R::dbeta(x, shape1, shape2, return_log);
}

//------------------------------------------------
// draw from symmetric dirichlet distribution using single shape parameter
// repeated d times
vector<double> rsym_dirichlet1(double shape, int d) {
  // draw a series of gamma random variables
  vector<double> ret(d);
  double ret_sum = 0;
  for (int i=0; i<d; i++) {
    ret[i] = rgamma1(shape, 1.0);
    ret_sum += ret[i];
  }
  // divide all by the sum
  double ret_inv_sum = 1.0/ret_sum;
  for (int i=0; i<d; i++) {
    ret[i] *= ret_inv_sum;
  }
  return(ret);
}

//------------------------------------------------
// probability density of symmetric Dirichlet distribution
double dsym_dirichlet1(const std::vector<double> &p, double shape, bool return_log) {
  int d = p.size();
  double ret = lgamma(d*shape);
  double lgs = lgamma(shape);
  for (int i=0; i<d; ++i) {
    ret += (shape - 1)*log(p[i]) - lgs;
  }
  if (!return_log) {
    return exp(ret);
  }
  return ret;
}

//------------------------------------------------
// draw from dirichlet distribution using vector of shape parameters
vector<double> rdirichlet1(const vector<double> &shape_vec) {
  // draw a series of gamma random variables
  int n = shape_vec.size();
  vector<double> ret(n);
  double retSum = 0;
  for (int i=0; i<n; i++) {
    ret[i] = rgamma1(shape_vec[i], 1.0);
    retSum += ret[i];
  }
  // divide all by the sum
  double retSumInv = 1.0/retSum;
  for (int i=0; i<n; i++) {
    ret[i] *= retSumInv;
  }
  return(ret);
}

//------------------------------------------------
// draw from dirichlet distribution using bespoke inputs. Outputs are stored in 
// x, passed by reference for speed. Shape parameters are equal to alpha+beta, 
// where alpha is an integer vector, and beta is a single double.
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
int rpois1(double rate) {
  return R::rpois(rate);
}

//------------------------------------------------
// probability mass of Poisson(rate) distribution
double dpois1(int n, double rate, bool return_log) {
  return R::dpois(n,rate,return_log);
}

//------------------------------------------------
// draw from negative binomial distribution with mean lambda and variance gamma*lambda (gamma must be >1)
int rnbinom1(double lambda, double gamma) {
  return R::rnbinom(lambda/(gamma-1), 1/gamma);
}

//------------------------------------------------
// probability mass of negative binomial distribution with mean lambda and
// variance gamma*lambda (gamma must be >1)
double dnbinom1(int n, double lambda, double gamma, bool return_log) {
  return R::dnbinom(n, lambda/(gamma-1), 1/gamma, return_log);
}

