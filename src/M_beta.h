//M_beta.h
#ifndef _EMVSpackage_M_BETA_H
#define _EMVSpackage_M_BETA_H
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
vec M_beta(mat &XtY, mat &X, mat &XtX, vec &inv_var);
#endif
