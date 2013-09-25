//log_prior.h
#ifndef _EMVSpackage_LOG_PRIOR_H
#define _EMVSpackage_LOG_PRIOR_H
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

double log_prior(uvec& gamma, const string &type, double a, double b, int n);
#endif
