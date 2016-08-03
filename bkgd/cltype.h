#ifndef __CLTYPE_H__
#define __CLTYPE_H__

#include <math.h>
#include "branch.h"
#include "bkgd_evo_mdl_param.h"

ColType *cltype_new(const char *name,  int subst_type,
		    void (*set_n)(ColType *, const BkgdBin *));

void cltype_set_br(ColType *cltype, Branch *br);

void cltype_add_dbl_br(ColType *cltype, Branch *br1, Branch *br2);


void cltype_set_prob(ColType *cltype, const BkgdEvoMdlParam *param,
		     const BkgdEvoMdlConfig *conf);

void cltype_set_dprob(ColType *cltype, const BkgdEvoMdlParam *param, 
		      const BkgdEvoMdlConfig *conf);

void cltype_cons_set_prob(ColType *cltype, const BkgdEvoMdlParam *param,
			  const BkgdEvoMdlConfig *conf);

void cltype_cons_set_dprob(ColType *cltype, const BkgdEvoMdlParam *param, 
			   const BkgdEvoMdlConfig *conf);


void cltype_free(ColType *cltype);


/**********************************
 * column type site counts for HC model
 **********************************/

/* sets number of HC column types in HC model */
void cltype_set_n__hc_h_c(ColType *cltype, const BkgdBin *bin);

/* sets number of CONS column types in HC model */
void cltype_set_n__hc_cons(ColType *cltype, const BkgdBin *bin);



/*********************************
 * column type site counts for HCM model
 *********************************/

void cltype_set_n__hcm_h_c(ColType *cltype, const BkgdBin *bin);

void cltype_set_n__hcm_m(ColType *cltype, const BkgdBin *bin);

void cltype_set_n__hcm_cons(ColType *cltype, const BkgdBin *bin);


/*********************************
 * column type site counts for HCOM model
 *********************************/

void cltype_set_n__hcom_h_c(ColType *cltype, const BkgdBin *bin);

void cltype_set_n__hcom_hc(ColType *cltype, const BkgdBin *bin);

void cltype_set_n__hcom_o(ColType *cltype, const BkgdBin *bin);

void cltype_set_n__hcom_m(ColType *cltype, const BkgdBin *bin);

void cltype_set_n__hcom_cons(ColType *cltype, const BkgdBin *bin);


/************************************
 * column type site counts for HCGOM model
 ************************************/

void cltype_set_n__hcgom_h(ColType *cltype, const BkgdBin *bin);

void cltype_set_n__hcgom_c(ColType *cltype, const BkgdBin *bin);

void cltype_set_n__hcgom_g(ColType *cltype, const BkgdBin *bin);

void cltype_set_n__hcgom_o(ColType *cltype, const BkgdBin *bin);

void cltype_set_n__hcgom_m(ColType *cltype, const BkgdBin *bin);

void cltype_set_n__hcgom_hc(ColType *cltype, const BkgdBin *bin);

void cltype_set_n__hcgom_hcg(ColType *cltype, const BkgdBin *bin);

void cltype_set_n__hcgom_hg_cg(ColType *cltype, const BkgdBin *bin);


void cltype_set_n__hcgom_hg(ColType *cltype, const BkgdBin *bin);

void cltype_set_n__hcgom_cg(ColType *cltype, const BkgdBin *bin);

void cltype_set_n__hcgom_ho_co(ColType *cltype, const BkgdBin *bin);

void cltype_set_n__hcgom_cons(ColType *cltype, const BkgdBin *bin);

void cltype_set_n__hcgo_cons(ColType *cltype, const BkgdBin *bin);



#endif
