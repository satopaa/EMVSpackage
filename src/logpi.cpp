//logpi.cpp
#include "logpi.h"
#include <math.h> 

double logpi(vec &x, double a, double b){
  double res = sum(a*x - (a+b)*log(1+exp(x)));
  return res;
}
