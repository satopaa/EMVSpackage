//M_beta.cpp
#include "M_beta.h"
#include <math.h> 

vec M_beta(mat &XtY, mat &X, mat &XtX, vec &inv_var){
  
  if(X.n_cols > X.n_rows){
    return zeros<vec>(X.n_cols); /// NOT IMPLEMENTED YET!!

  } else {
    mat Psi = diagmat(inv_var);
    Psi += XtX;
    Psi = symmatu(Psi);
    Psi = solve(Psi, XtY);
    return Psi;
  }
}
