//density_norm.cpp
#include "density_norm.h"
#include <math.h> 

vec density_norm(vec &x, double mu, double sigma){
  vec dens = 1/(sigma*sqrt(2*PI)) * ones<vec>(x.n_elem);    
  dens %= exp(-square(x-mu)/(2 * pow(sigma, 2))); 
  return dens;
}
