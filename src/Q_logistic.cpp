//Q_logistic.cpp
#include "Q_logistic.h"
#include "logpi.h"
#include <math.h> 

double Q_logistic(const std::vector<double> &x, std::vector<double> &grad, void* f_data){
  vec theta(x);
  mat* pAuxiliary = (mat*)(f_data);
  mat Z = pAuxiliary[0];
  //Z.print();
  
  vec probs = pAuxiliary[1];
  //probs.print();

  double a = as_scalar(pAuxiliary[2]);
  double b = as_scalar(pAuxiliary[3]);
  //Rcout << a << b << endl;

  double res = sum(probs%(Z*theta));
  res -= sum(log(1+exp(Z*theta)));
  res += sum( logpi(theta,a,b) );
  //temp += logpi(theta,a,b);
  //temp.print();
  //theta.print();
  return (-res);
}
