#include "EMVS.h"
#include "E_beta_binom.h"
#include "E_logistic.h"
#include "density_norm.h"
#include "delogit.h"
#include "beta.h"
#include "log_g.h"
#include "M_beta.h"
#include "M_theta.h"
#include "M_sigma.h"
#include "M_p.h"
#include "Q_logistic.h"

#include <math.h> 
#include <ctime>

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

  clock_t tstart, tend; // Debugging
  bool debug = true;

  // First recast all R-types to C++-types: MUST BE DONE FASTER!
  if(debug)tstart = clock(); // START
  vec Y = as<vec>(Y_R);  
  mat X = as<mat>(X_R);  
  //mat Xt = X.t();  
  if(debug)Rcout << "Setup (Y, X, Xt) took "<< difftime(clock(), tstart)/CLOCKS_PER_SEC <<" second(s)."<< endl;


  if(debug)tstart = clock(); // START
  vec v0s = as<vec>(v0s_R);
  double v0, v_k, p_k, p, eps;
  int niter, i;
  double sigma_init = as<double>(sigma_init_R);
  double epsilon = as<double>(epsilon_R);

  const string type = as<string>(type_R);
  bool betabinomial = (type.compare("betabinomial") == 0);
  bool fixed = (type.compare("fixed") == 0);
  bool logistic = (type.compare("logistic") == 0);
  bool MRF = (type.compare("MRF") == 0);


    
  const int L = v0s.n_elem;
  const int dim = X.n_cols;
  const int n = X.n_rows;

  bool v1_missing = Rf_isNull(v1_R);
  double v1;
  vec v1s;
  if(v1_missing){
    //v1 = 1000; // NOTE NOTE NOTE NOTE
    v1s = zeros<vec>(L); 
  } else {
    v1 = as<double>(v1_R);
    v1s = ones<vec>(L)*v1;    
  }

  if(!Rf_isNull(p_R)){
    p = as<double>(p_R);
  }
  
  mat Z;
  if(logistic){
    if(!Rf_isNull(Z_R)){
      Z = as<mat>(Z_R);
    } else {
      Rcout << "For the logistic prior Z must be specified" << endl;
    }
  }
  
  if(debug)Rcout << "Section 2 took "<< difftime(clock(), tstart)/CLOCKS_PER_SEC <<" second(s)."<< endl;

  if(debug)tstart = clock(); // START

  double a,b;
  if(!Rf_isNull(a_R) && !Rf_isNull(b_R)){
     a = as<double>(a_R);
     b = as<double>(b_R);
  }
  
  double v1_g;
  if(!Rf_isNull(v1_g_R)){
    v1_g = as<double>(v1_g_R);
  }
  
  bool temperature_missing = Rf_isNull(temperature_R);
  double temperature;
  if(temperature_missing){
    temperature = 1;
  } else {
    temperature = as<double>(temperature_R);
  }

  bool beta_init_missing = Rf_isNull(beta_init_R);
  vec beta_init;
  if(!beta_init_missing){
    beta_init = as<vec>(beta_init_R);
  } 
  if(debug)Rcout << "Section 3 took "<< difftime(clock(), tstart)/CLOCKS_PER_SEC <<" second(s)."<< endl;

  vec intersects(L); 
  vec sigmas(L);
  mat betas(L,dim);
  mat posts(L,dim);  
  vec log_post(L);
  vec niters(L);

  vec theta_k(dim);
  vec linpred(dim);
  mat E_step;

  if(debug)tstart = clock(); // START
  mat XtY = X.t()*Y; 
  if(debug)Rcout << "Section 4 took "<< difftime(clock(), tstart)/CLOCKS_PER_SEC <<" second(s)."<< endl;
  if(debug)tstart = clock(); // START
  mat XtX = symmatu(X.t()*X); // Transpose with itself is optimized!
  if(debug) Rcout << "Section 5 took "<< difftime(clock(), tstart)/CLOCKS_PER_SEC <<" second(s)."<< endl;
  uvec index;
  vec beta_k, beta_new;
  double sigma_k, c, w;
  vec inv_var, inv_var_temp;
  vec post;


  //  Rcout << " Iteration begins: " << endl;
  for(i = 0; i < L; ++i){
    v0 = v0s[i];
    
    if(beta_init_missing){
      inv_var_temp = ones<vec>(dim)*( (v0s[i]+v1)/(2*v0s[i]*v1) );
      beta_k = M_beta(XtY, X, XtX, inv_var_temp); // NEEDS TO IMPLEMENT WOODBURRY FORMULA
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
    
    if(betabinomial){
      p_k = 0.5;
    } else if (fixed){
      p_k = p;
    } else if(logistic){
      theta_k = zeros<vec>(Z.n_cols);
    }
    eps = epsilon + 1;
    niter = 1;
    
    /////////////////////////// TESTING ////////////////////////////
    /*     
	 vec beta_k_test;
    beta_k_test << 1 << 2 << 3;
    double sigma_test = 2;
    double   v0_test = 0.2;
    double   v1_test = 0.8;
    double   t_test = 0.8;
    double   p_test = 0.5;

    mat test_mat = E_beta_binom(beta_k_test, sigma_test, v0_test, v1_test, p_test, t_test);
    test_mat.print("testing:");


    vec inv_var_test = ones<vec>(500);
    vec test_mat2 = M_beta(XtY, X, XtX, inv_var_test);
    test_mat2.print("testing2:");


    beta_k_test = ones<vec>(500)*2.1;
    double test_double = M_sigma(Y, X, beta_k_test, inv_var_test, 1, 1);
    Rcout << "Testing: " << test_double << endl;
    
    beta_k_test << 1 << 2 << 3;
    test_double = M_p(beta_k_test, 1.2, 3.5);
    Rcout << "Testing: " << test_double << endl;
    */
       
    vec beta_k_test;
    beta_k_test  << 1 << 2 << 3;
    double sigma_test = 2;
    double v0_test = 0.2;
    double v_k_test = 0.4;
    vec exp_linpred_test;
    exp_linpred_test << 0.3 << 0.5 << 0.8;
    double temperature_test = 0.8;
    mat test_result = E_logistic(beta_k_test, sigma_test, v0_test, v_k_test, exp_linpred_test, temperature_test);
    //test_result.print();

    double a_test = 1.3;
    double b_test = 1.5;
    mat Z_test = diagmat(exp_linpred_test);
    vec theta_k_test;
    theta_k_test << 1.3 << 0.5 << 0.7;
    vec post_test;
    post_test << 0.2 << 0.3 << 0.5;
    vec theta_test = M_theta(theta_k_test, Z_test, post_test, a_test, b_test);
    theta_test.print();
    
    //mat auxiliary[4];
    //auxiliary[0] = Z_test;
    //auxiliary[1] = post_test;
    //auxiliary[2] = a_test;
    //auxiliary[3] = b_test;
    //std::vector<double> grad_test (3,0.6);
    //std::vector<double> x_test (3,0.8);
    //double res_test = Q_logistic(x_test, grad_test, auxiliary);
    //Rcout << res_test << endl;
      

    /////////////////////////// TESTING ////////////////////////////

    // test_a = 2.0;
    // test_b = 4.0;
    

    //sleep(10);

    //while(FALSE){

    while(eps > epsilon){
      //Rcout << "v0 = " << v0 << "; iter = " << niter++ << endl;
      
      if(betabinomial || fixed){
	if(debug)tstart = clock(); // START
	E_step = E_beta_binom(beta_k, sigma_k, v0, v_k, p_k, temperature);
	if(debug)Rcout << "E-step took "<< difftime(clock(), tstart)/CLOCKS_PER_SEC <<" second(s)."<< endl;
	//Rcout << "EXPECT STEP!!" << endl;
      }  else if(MRF){
	// NOT IMPLEMENTED YET
      }  else if(logistic){
	linpred = Z * theta_k;
	vec exp_linpred =  delogit(linpred);
	E_step = E_logistic(beta_k, sigma_k, v0, v_k, exp_linpred, temperature); // This works!	
      }
      inv_var = E_step.col(0);
      post = E_step.col(1);
      
      if(debug) tstart = clock(); // START
      beta_k = M_beta(XtY, X, XtX, inv_var);
      if(debug)Rcout << "M-step beta took "<< difftime(clock(), tstart)/CLOCKS_PER_SEC <<" second(s)."<< endl;

      if(debug)tstart = clock(); // START
      sigma_k = M_sigma(Y,X,beta_k, inv_var, 1, 1);
      if(debug)Rcout << "M-step sigma_k took "<< difftime(clock(), tstart)/CLOCKS_PER_SEC <<" second(s)."<< endl;
      
      if(betabinomial){
	if(debug)tstart = clock(); // START
	p_k = M_p(post, a, b);
	if(debug)Rcout << "M-step p_k took "<< difftime(clock(), tstart)/CLOCKS_PER_SEC <<" second(s)."<< endl;
      }  else if(logistic){
	if(debug)tstart = clock(); // START
	theta_k = M_theta(theta_k, Z, post, a, b); // This works!
	if(debug)Rcout << "M-step theta took "<< difftime(clock(), tstart)/CLOCKS_PER_SEC <<" second(s)."<< endl;
      }
      
      if(v1_missing){
	// NOT IMPLEMENTED YET
      }
      
      eps = max(abs(beta_new - beta_k));
      beta_new = beta_k;
      //Rcout << eps << endl;
    }
    
    // Store values:
    posts.row(i) = post.t();
    betas.row(i) = beta_new.t();
    sigmas[i] = sigma_k;
    v1s[i] = v_k;

    c = sqrt(v1s[i]/v0s[i]);
    if(betabinomial || fixed){
      w = (1-p_k)/p_k;
      intersects[i] = sigmas[i];
      intersects[i] *= sqrt(v0s[i]);
      intersects[i] *= sqrt(2*log(w*c)*c*c/(c*c-1));
    } 
    index = find(post > 0.5);

    log_post[i] = log_g(index,X,Y,1,1,0,v1_g,"betabinomial",a,b);
    niters[i] = niter;
  }
 

 
  /////////////////////////////////////////////////////////////////////////////
  // Wrap the results into a list and return.
  /////////////////////////////////////////////////////////////////////////////
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
  
  if(betabinomial || fixed){
    list["p_k"] = p_k;
  } else if (logistic){
    list["theta"] = theta_k;
  } 
  return list;
}

