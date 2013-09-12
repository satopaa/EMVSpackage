#ifndef _EMVSpackage_EMVS_H
#define _EMVSpackage_EMVS_H

#include <RcppArmadillo.h>

RcppExport SEXP EMVS(SEXP Y_R,
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
		     SEXP v1_g_R);
#endif
