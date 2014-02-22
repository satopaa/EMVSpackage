//E_logistic.h
#ifndef _EMVSpackage_E_LOGISTIC_H
#define _EMVSpackage_E_LOGISTIC_H
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
mat E_logistic(vec &beta_k, double sigma_k, double v0, double v1, vec &p, double t);
#endif
