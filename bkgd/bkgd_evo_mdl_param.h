#ifndef __BKGD_EVO_MDL_PARAM_H__
#define __BKGD_EVO_MDL_PARAM_H__

#include <math.h>
typedef struct {
  double T_hc;    /* time to H/C speciation in generations */
  double T_hcg;   /* time b/w HC and HC/G speciation */
  double T_hcgo;  /* time b/w HCG and HCG/O speciation */
  double T_hcgom; /* time b/w HCGO and HCGO/M speciation */
    
  double N_hc;    /* effective size of ancestral HC population */
  double N_hcg;   /* effective size of ancestral HCG population */
  double N_hcgo;  /* effective size of ancestral HCGO population */
  double N_hcgom; /* effective size of ancestral HCGOM population */

  /* HCM model only */
  double T_hcm;   /* time b/w HC and HC/M speciation (used for HCM model) */
  double N_hcm;

  /* HCOM model only */
  double T_hco;
  double T_hcom;
  double N_hco;
  double N_hcom;

  /* mutation rate coefficients */
  double mu_a;
  double mu_b;
  double mu_c;

  double mu; /* mutation rate */

  /* transition/transversion rates of mutation */
  double mu_i;
  double mu_v;  

  double u_ex_scale;
  double u_nex_scale;

  double B;

  double k_hcg; /* prob HC coalescent predates HC/G speciation */

  /* scaling factors for double substitution probabilities */
  double lambda; 
  double lambda_i;
  double lambda_v;

} BkgdEvoMdlParam;


/**
 * Sets all attributes of provided parameter
 * structure to 0.
 */
void param_set_zero(BkgdEvoMdlParam *param);

void param_add(const BkgdEvoMdlParam *from, const double coef,
	       BkgdEvoMdlParam *to);

BkgdEvoMdlParam *bkgd_evo_param_new();


#endif
