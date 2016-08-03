
#include <glib.h>
#include <math.h>

#include "util.h"

#include "bkgd_evo_mdl_param.h"

/**
 * Sets all attributes of provided parameter
 * structure to 0.
 */
void param_set_zero(BkgdEvoMdlParam *param) {
  param->u_ex_scale = 0.0;
  param->u_nex_scale = 0.0;

  param->T_hc = 0.0;
  param->T_hcg = 0.0;
  param->T_hcgo = 0.0;
  param->T_hcgom = 0.0;
  param->T_hcm = 0.0;
  param->T_hco = 0.0;
  param->T_hcom = 0.0;

  param->N_hc = 0.0;
  param->N_hcg = 0.0;
  param->N_hcgo = 0.0;
  param->N_hcgom = 0.0;
  param->N_hcm = 0.0;
  param->N_hco = 0.0;
  param->N_hcom = 0.0;

  param->mu = 0.0;
  param->mu_i = 0.0;
  param->mu_v = 0.0;
  param->mu_a = 0.0;
  param->mu_b = 0.0;
  param->mu_c = 0.0;

  param->B = 0.0;
  param->k_hcg = 0.0;

  param->lambda = 0.0;
  param->lambda_i = 0.0;
  param->lambda_v = 0.0;
}


/**
 * adds from attributes (multiplied by coef) to to attributes 
 */
void param_add(const BkgdEvoMdlParam *from, const double coef,
		      BkgdEvoMdlParam *to) {

  if(isnan(coef)) {
    g_error("param_add: coef is nan");
  }
  if(isinf(coef)) {
    g_error("param_add: coef is inf");
  }

/*   if(((coef > 0.0 )&& ((coef < 1e-50) || (coef > 1e50))) || */
/*      ((coef < 0.0) && ((coef < -1e50) || (coef > -1e-50)))) { */
/*     g_warning("param_add: large or small coef: %g", coef); */
/*   } */
  
  to->T_hc    += coef * from->T_hc;
  to->T_hcg   += coef * from->T_hcg;
  to->T_hcgo  += coef * from->T_hcgo;
  to->T_hcgom += coef * from->T_hcgom;
  to->T_hcm   += coef * from->T_hcm;
  to->T_hco   += coef * from->T_hco;
  to->T_hcom  += coef * from->T_hcom;

  to->N_hc    += coef * from->N_hc;
  to->N_hcg   += coef * from->N_hcg;
  to->N_hcgo  += coef * from->N_hcgo;
  to->N_hcgom += coef * from->N_hcgom;
  to->N_hcm   += coef * from->N_hcm;
  to->N_hco   += coef * from->N_hco;
  to->N_hcom   += coef * from->N_hcom;

  to->u_ex_scale  += coef * from->u_ex_scale;
  to->u_nex_scale += coef * from->u_nex_scale;

  to->mu   += coef * from->mu;
  to->mu_i += coef * from->mu_i;
  to->mu_v += coef * from->mu_v;
  to->mu_a += coef * from->mu_a;
  to->mu_b += coef * from->mu_b;
  to->mu_c += coef * from->mu_b;

  to->lambda += coef * from->lambda;
  to->lambda_i += coef * from->lambda_i;
  to->lambda_v += coef * from->lambda_v;
}


/**
 * Creates new BkgdEvoMdlParam structure with
 * attributes initialized to 0
 */
BkgdEvoMdlParam *bkgd_evo_param_new() {
  BkgdEvoMdlParam *param;

  param = g_new(BkgdEvoMdlParam, 1);
  param_set_zero(param);
  
  return param;
}
