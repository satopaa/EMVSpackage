//density_norm.cpp
#include "density_norm.h"
#include <math.h> 
#define PI 3.1415926535897

vec density_norm(vec &x, double mu, double sigma){
  vec temp = 1/(sigma*sqrt(2*PI)) * ones<vec>(x.n_elem);    
  temp %= exp(-square(x-mu)/(2 * pow(sigma, 2))); 
  return temp;
}
