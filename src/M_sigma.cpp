//M_sigma.cpp
#include "M_sigma.h"
#include <math.h> 

double M_sigma(mat &Y, mat &X, vec &beta_k, vec &inv_var, double eta, double lambda){
  vec e = Y - X*beta_k;
  vec e2 = diagmat(sqrt(inv_var)) * beta_k;
  
  double res = as_scalar(e.t() * e);
  res += as_scalar(e2.t() * e2);
  res += eta*lambda;
  res /= X.n_rows +  X.n_cols + eta;
  res = sqrt(res);
  return res;
}
