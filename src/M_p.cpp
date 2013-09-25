//M_p.cpp
#include "M_p.h"
#include <math.h> 

double M_p(vec &post, double a, double b){
  double res = accu(post) + a -1;
  res /= (b+a+post.n_elem-2);
  return res;
}
