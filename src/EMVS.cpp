#include "EMVS.h"
#include "E_beta_binom.h"
#include "density_norm.h"
#include "delogit.h"
#include "beta.h"
#include "log_g.h"
#include "M_beta.h"
#include "M_sigma.h"
#include "M_p.h"

#include <math.h> 

using namespace std;
using namespace Rcpp;
using namespace arma;

/*
# Y          ... n x 1 vector of centered responses
# X          ... n x p matrix of standardized predictors
# v0s        ... sequence of spike variance parameters
# v1         ... slab variance parameter; if missing then imposed a prior and estimated
# beta_init  ... starting values for the beta vector; if missing then the limiting
#		     case of deterministic annealing used
# sigma_init ... a starting value for sigma
# epsilon    ... convergence margin parameter
# type 	 ... type of the prior on the model space
#                type="betabinomial"  for the betabinomial prior with hyperparameters a and b
#	           type="MRF"           for the MRF prior  with hyperparameters theta and Sigma
#                type="fixed"         for the fixed inclusion probability p
#		     type="logistic"	  for the logistic regression prior with a grouping matrix Z
# temperature... inverse temperature parameter for annealing; if missing then
#		     t=1 without annealing
# mu         ... sparsity MRF hyper-parameter, where mu(1,...1)' is the mean vector of the MRF prior
# Sigma      ... smoothness matrix in MRF prior
# p          ... fixed prior probability of inclusion, if type="fixed"			
# a,b        ... hyperparameters for the betabinomial prior, if type="betabinomial"
# a_v1,b_v1  ... parameters of the prior for v1, if v1 missing
# v1_g       ... v1 value for the model evaluation by the g-function
# Z          ... group identification matrix for the logistic prior
# MRF implemented here with fixed sparsity parameter mu
*/

SEXP EMVS(SEXP Y_R,
	  SEXP X_R,
	  SEXP v0s_R,
	  SEXP v1_R,
	  SEXP type_R,
	  SEXP beta_init_R,
	  SEXP sigma_init_R,
	  SEXP epsilon_R,
	  SEXP temperature_R,
	  SEXP Z_R,
	  SEXP mu_R,
	  SEXP Sigma_R,
	  SEXP p_R,
	  SEXP a_R,
	  SEXP b_R,
	  SEXP a_v1_R,
	  SEXP b_v1_R,
	  SEXP v1_g_R){

  // First recast all R-types to C++-types:
  vec Y = as<vec>(Y_R);  
  mat X = as<mat>(X_R);  
  mat Xt = X.t();  

  vec v0s = as<vec>(v0s_R);
  double v0, v1, v_k, p_k, p, eps;
  int niter, i;
  double sigma_init = as<double>(sigma_init_R);
  string type = as<string>(type_R);
  double epsilon = as<double>(epsilon_R);
    
  int L = v0s.n_elem;
  int dim = X.n_cols;
  int n = X.n_rows;

  bool v1_missing = Rf_isNull( v1_R );
  vec v1s;
  if(v1_missing){
    v1 = 1000; // NOTE NOTE NOTE NOTE
    v1s = zeros<vec>(L);
  } else {
    v1 = as<double>(v1_R);
    v1s = ones<vec>(L)*v1;    
  }

  if(!Rf_isNull(p_R)){
    p = as<double>(p_R);
  }
  double a, b, v1_g;
   if(!Rf_isNull(a_R)){
    a = as<double>(a_R);
  }
 if(!Rf_isNull(b_R)){
    b = as<double>(b_R);
  }
 if(!Rf_isNull(v1_g_R)){
   v1_g = as<double>(v1_g_R);
  }

  bool temperature_missing = Rf_isNull( temperature_R );
  double temperature;
  if(temperature_missing){
    temperature = 1;
  } else {
    temperature = as<double>(temperature_R);
  }

  bool beta_init_missing = Rf_isNull( beta_init_R );
  vec beta_init;
  if(!beta_init_missing){
    beta_init = as<vec>(beta_init_R);
  } 

  vec intersects = zeros<vec>(L);
  vec sigmas = zeros<vec>(L);
  mat betas = zeros<mat>(L,dim);
  mat posts = zeros<mat>(L,dim);
  
  vec log_post = zeros<vec>(L);
  vec niters = zeros<vec>(L);
  mat E_step;

  mat XtY = Xt*Y;
  mat XtX = Xt*X;
  uvec index;
  vec beta_k, beta_new;
  double sigma_k, c, w;
  vec inv_var;
  vec post;

  Rcout << " Iteration begins: " << endl;
  for(i = 0; i < L; ++i){
    v0 = v0s[i];
    
    if(beta_init_missing){
      inv_var = ones<vec>(dim)*( (v0s[i]+v1)/(2*v0s[i]*v1) );
      beta_k = M_beta(XtY, X, XtX, inv_var); // NEEDS TO IMPLEMENT WOODBURRY FORMULA
    } else {
      beta_k = beta_init;
    }
    
    beta_new = beta_k;
    sigma_k = sigma_init;
    
    if(v1_missing){
      v_k = 1000;
    } else {
      v_k = v1;
    }
    
    if(type.compare("betabinomial")){
      p_k = 0.5;
    } else if (type.compare("fixed")){
      p_k = p;
    } else if(type.compare("logistic")){
      // NOT IMPLEMENTED YET
    }
    eps = epsilon + 1;
    niter = 1;
    
    while(eps > epsilon){
      Rcout << "v0 = " << v0 << "; iter = " << niter++ << endl;
      
      if(type.compare("betabinomial") || type.compare("fixed")){
	E_step = E_beta_binom(beta_k, sigma_k, v0, v_k, p_k, temperature);
      }  else if(type.compare("MRF")){
	// NOT IMPLEMENTED YET
      }  else if(type.compare("logistic")){
	// NOT IMPLEMENTED YET
      }
      inv_var = E_step.col(0);
      post = E_step.col(0);
      
      beta_k = M_beta(XtY, X, XtX, inv_var);
      sigma_k = M_sigma(Y,X,beta_k, inv_var, 1, 1);
      
      if(type.compare("betabinomial")){
	p_k = M_p(post, a, b);
      }  else if(type.compare("logistic")){
	// NOT IMPLEMENTED YET
      }
      
      if(!Rf_isNull(v1_R)){
	// NOT IMPLEMENTED YET
      }
      
      eps = max(abs(beta_new - beta_k));
      beta_new = beta_k;
      Rcout << eps << endl;
    }
    
    // Store values:
    posts.row(i) = post.t();
    betas.row(i) = beta_new.t();
    sigmas[i] = sigma_k;
    v1s[i] = v_k;

    c = sqrt(v1s[i]/v0s[i]);
    if(type.compare("betabinomial") || type.compare("fixed")){
      w = (1-p_k)/p_k;
      intersects[i] = sigmas[i];
      intersects[i] *= sqrt(v0s[i]);
      intersects[i] *= sqrt(2*log(w*c)*c*c/(c*c-1));
    } 
    index = find(post > 0.5);

    log_post[i] = log_g(index,X,Y,1,1,0,v1_g,type,a,b);
    //log_post[i] = 1;
    niters[i] = niter;
  }
  

  // Wrap the results into a list and return.
  List list;
  list["betas"] = betas;
  list["log_post"] = log_post;
  list["intersects"] = intersects;
  list["sigmas"] = sigmas;
  list["v1s"] = v1s;
  list["v0s"] = v0s;
  list["niters"] = niters;
  list["posts"] = posts;
  list["inv_var"] = inv_var;
  list["type"] = type;
  
  if(type.compare("betabinomial") || (type.compare("fixed"))){
    list["p_k"] = p_k;
  } else if (type.compare("logistic")){
      //list["theta"] = "theta_k";
  } 
  return list;
}

