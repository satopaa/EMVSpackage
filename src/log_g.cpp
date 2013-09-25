//log_g.cpp
#include "log_g.h"
#include "log_prior.h"
#include <math.h> 

double log_g(uvec &gamma, mat &X, mat &Y, double eta, double lambda, double v0, double v1, const string &type, double a, double b){
  int q = gamma.n_elem;
  double Ssq, log_val;
  int n = Y.n_elem;
  int N = X.n_cols;
  if(q > 0){
    mat X_gamma = X.cols(gamma);
    vec diag = ones<vec>(q)*v1;
    mat X_tilde = join_cols(X_gamma, diagmat(sqrt(1/diag)));
    Ssq = as_scalar(Y.t()*Y - Y.t()*X_gamma*inv(X_gamma.t() *X_gamma + diagmat(1/diag)) * X_gamma.t() * Y);
    log_val = -0.5*log(det(X_tilde.t()*X_tilde));
    log_val -= q/2*log(v1);
    log_val -= (n+eta)/2*log(eta*lambda+Ssq);
    log_val += log_prior(gamma,type,a,b,N);
  } else {
    Ssq = as_scalar(Y.t()*Y);
    log_val =  -(n+eta)/2*log(eta*lambda+Ssq);
    log_val += log_prior(gamma,type,a,b,N);
  }
  return log_val;
}
