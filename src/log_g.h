//log_g.h
#ifndef _EMVSpackage_LOG_G_H
#define _EMVSpackage_LOG_G_H
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
double log_g(uvec &gamma, mat *X, mat *Y, double eta, double lambda, double v0, double v1, string type, );
#endif
