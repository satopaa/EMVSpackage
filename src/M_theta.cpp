//M_theta.cpp
#include "M_theta.h"
#include "Q_logistic.h"
#include <math.h> 

vec M_theta(vec &b_k, mat &Z, vec &post, double a, double b){
  int p = b_k.n_elem;

  opt opt(nlopt::LN_NELDERMEAD, p); // This is used by R's optim(...)
  mat auxiliary[4];
  auxiliary[0] = Z;
  auxiliary[1] = post;
  auxiliary[2] = a;
  auxiliary[3] = b;
  
  //std::vector<double> lb(p, 0);
  //opt.set_lower_bounds(lb);
  opt.set_min_objective(Q_logistic, &auxiliary);  
  opt.set_xtol_rel(1e-8);
  
  std::vector<double> init = arma::conv_to< std::vector< double > >::from(b_k);
  //std::vector<double> init = <vector> b_k;
  double minf;
  nlopt::result result = opt.optimize(init, minf);

  vec res(init);  
  return res;
}
