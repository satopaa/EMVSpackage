//Q_logistic.h
#ifndef _EMVSpackage_Q_LOGISTIC_H
#define _EMVSpackage_Q_LOGISTIC_H
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
double Q_logistic(const std::vector<double> &x, std::vector<double> &grad, void* f_data);
#endif
