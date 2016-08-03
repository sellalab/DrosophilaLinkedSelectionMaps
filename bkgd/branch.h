
#ifndef __BRANCH_H__
#define __BRANCH_H__

#include <math.h>
#define BRANCH_LOW_PROB 1e-20
#define BRANCH_LEN_SMALL 1e-12


#include "bkgd_evo_mdl_param.h"
#include "bkgd_evo_mdl.h"


void branch_set_prob(Branch *, const BkgdEvoMdlParam *,
		     const BkgdEvoMdlConfig *conf);

void branch_set_dprob(Branch *br, const BkgdBin *bin, 
		      const BkgdEvoMdlParam *param, 
		      const BkgdEvoMdlConfig *conf);



/********************
 * HC model functions
 ********************/

/* branch length & probabilities */
void branch_set_len__hc_h_c(Branch *, const BkgdEvoMdlParam *);

/* branch length partial derivatives */
void branch_set_dlen__hc_h_c(Branch *, const BkgdBin *, 
			     const BkgdEvoMdlParam *);


/********************
 * HCM model functions
 ********************/

/* branch length & probabilities */
void branch_set_len__hcm_h_c(Branch *, const BkgdEvoMdlParam *);

void branch_set_len__hcm_m(Branch *, const BkgdEvoMdlParam *);

/* branch length partial derivatives */
void branch_set_dlen__hcm_h_c(Branch *, const BkgdBin *, 
			      const BkgdEvoMdlParam *);

void branch_set_dlen__hcm_m(Branch *, const BkgdBin *, 
			    const BkgdEvoMdlParam *);



/********************
 * HCOM model functions
 ********************/

/* branch length & probabilities */
void branch_set_len__hcom_h_c(Branch *, const BkgdEvoMdlParam *);

void branch_set_len__hcom_hc(Branch *, const BkgdEvoMdlParam *);

void branch_set_len__hcom_o(Branch *, const BkgdEvoMdlParam *);

void branch_set_len__hcom_m(Branch *, const BkgdEvoMdlParam *);

/* branch length partial derivatives */
void branch_set_dlen__hcom_h_c(Branch *, const BkgdBin *, 
			       const BkgdEvoMdlParam *);

void branch_set_dlen__hcom_hc(Branch *, const BkgdBin *, 
			      const BkgdEvoMdlParam *);

void branch_set_dlen__hcom_o(Branch *, const BkgdBin *, 
			     const BkgdEvoMdlParam *);

void branch_set_dlen__hcom_m(Branch *, const BkgdBin *, 
			     const BkgdEvoMdlParam *);



/********************************************************
 * HCGOM model functions
 ********************************************************/

/* branch length & probabilities */
void branch_set_len__hcgom_h(Branch *br, const BkgdEvoMdlParam *param);

void branch_set_len__hcgom_c(Branch *br, const BkgdEvoMdlParam *param);

void branch_set_len__hcgom_g(Branch *br, const BkgdEvoMdlParam *param);

void branch_set_len__hcgom_hc(Branch *br, const BkgdEvoMdlParam *param);

void branch_set_len__hcgom_hg_cg(Branch *br, const BkgdEvoMdlParam *param);


void branch_set_len__hcgom_hg(Branch *br, const BkgdEvoMdlParam *param);

void branch_set_len__hcgom_cg(Branch *br, const BkgdEvoMdlParam *param);


void branch_set_len__hcgom_hcg(Branch *br, const BkgdEvoMdlParam *param);

void branch_set_len__hcgom_o(Branch *br, const BkgdEvoMdlParam *param);

void branch_set_len__hcgom_m(Branch *br, const BkgdEvoMdlParam *param);


/* partial derivatives */
void branch_set_dlen__hcgom_h(Branch *br, const BkgdBin *bin,
			      const BkgdEvoMdlParam *param);

void branch_set_dlen__hcgom_c(Branch *br, const BkgdBin *bin,
			      const BkgdEvoMdlParam *param);

void branch_set_dlen__hcgom_g(Branch *br, const BkgdBin *bin, 
			      const BkgdEvoMdlParam *param);

void branch_set_dlen__hcgom_hc(Branch *, const BkgdBin *bin, 
			       const BkgdEvoMdlParam *param);

void branch_set_dlen__hcgom_hg_cg(Branch *br, const BkgdBin *bin, 
				  const BkgdEvoMdlParam *param);

void branch_set_dlen__hcgom_hg(Branch *br, const BkgdBin *bin, 
			       const BkgdEvoMdlParam *param);

void branch_set_dlen__hcgom_cg(Branch *br, const BkgdBin *bin, 
			       const BkgdEvoMdlParam *param);

void branch_set_dlen__hcgom_hcg(Branch *br, const BkgdBin *bin, 
				const BkgdEvoMdlParam *param);

void branch_set_dlen__hcgom_o(Branch *br, const BkgdBin *bin, 
			      const BkgdEvoMdlParam *param);

void branch_set_dlen__hcgom_m(Branch *br, const BkgdBin *bin, 
			      const BkgdEvoMdlParam *param);


#endif

