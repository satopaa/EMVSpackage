//density_norm.h
#ifndef _EMVSpackage_DNORM_H
#define _EMVSpackage_DNORM_H
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
vec density_norm(vec &x, double mu, double sigma);
#endif
