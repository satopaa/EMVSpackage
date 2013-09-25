//log_prior.cpp
#include "log_prior.h"
#include "beta.h"
#include <math.h> 
#include <RcppArmadillo.h>

double log_prior(uvec& gamma, const string &type, double a, double b, int n){
  double res;
  if(type.compare("fixed") == 0){
    // NOT IMPLEMENTED YET
  } else if (type.compare("betabinomial") == 0){
    double x = gamma.n_elem + a;
    double y = n - gamma.n_elem + b;
    res = log(beta(x,y)) - log(beta(a,b));
    
    // Stirling approximation:
    if(!is_finite(res)){
      res = 0.5*log(2*PI)+(x-0.5)*log(x)+(y-0.5)*log(y)-(x+y-0.5)*log(x+y);
    }
    
  } else if (type.compare("MRF") == 0){
    // NOT IMPLEMENTED YET
  }
  return res;
}
