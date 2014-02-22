//E_logistic.cpp
#include "E_logistic.h"
#include "density_norm.h"
#include <math.h> 

mat E_logistic(vec &beta_k, double sigma_k, double v0, double v1, vec &p, double t){
  mat stars(beta_k.n_elem, 2);

  // Compute p_star
  vec dens1 = density_norm(beta_k, 0.0, sigma_k * sqrt(v1));
  vec dens0 = density_norm(beta_k, 0.0, sigma_k * sqrt(v0));
  stars.col(1) = pow(p%dens1, t);
  stars.col(1) /= (pow(p%dens1,t) + pow((1-p)%dens0, t));

  // Compute d_star
  stars.col(0) =  stars.col(1)/v1 + (1-stars.col(1))/v0;  
  return stars;
}
