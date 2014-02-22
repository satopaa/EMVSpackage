//MFA.cpp
#include "MFA.h"
#include <math.h> 

vec MFA(vec &start, vec &mu, mat &Sigma){
  vec E = start;
  double Eps = 1;
  vec E_new = E;
  vec f;
  while(Eps > 0.001){
    f = mu + Sigma*E;
    E = 1/(1+exp(-f));
    Eps = sum(abs(E_new-E));
    E_new = E;
  }
  return E_new;
}
