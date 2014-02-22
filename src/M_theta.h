//M_theta.h
#ifndef _EMVSpackage_M_THETA_H
#define _EMVSpackage_M_THETA_H
#include <RcppArmadillo.h>
#include <nlopt.hpp>
using namespace Rcpp;
using namespace arma;
using namespace nlopt;
vec M_theta(vec &b_k, mat &Z, vec &post, double a, double b);
#endif
