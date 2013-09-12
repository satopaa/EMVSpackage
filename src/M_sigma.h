//M_sigma.h
#ifndef _EMVSpackage_M_SIGMA_H
#define _EMVSpackage_M_SIGMA_H
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
double M_sigma(mat &Y, mat &X, vec &beta_k, vec &inv_var, double eta, double lambda);
#endif
