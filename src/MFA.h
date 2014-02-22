//MFA.h
#ifndef _EMVSpackage_MFA_H
#define _EMVSpackage_M_MFA_H
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
vec MFA(vec &start, vec &mu, mat &Sigma);
#endif
