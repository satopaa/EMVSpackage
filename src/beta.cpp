//beta.cpp
#include "beta.h"
#include <math.h> 

double beta(double x, double y){
  double temp = tgamma(x);
  temp *= tgamma(y);
  temp /= tgamma(x + y);
  return temp;
}
