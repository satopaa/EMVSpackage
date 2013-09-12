//E_beta_binom.h
#ifndef _EMVSpackage_E_BETA_BINOM_H
#define _EMVSpackage_E_BETA_BINOM_H
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
mat E_beta_binom(vec &beta_k, double sigma_k, double v0, double v1, double p, double t);
#endif
