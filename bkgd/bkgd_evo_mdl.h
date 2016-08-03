#ifndef __BKGD_EVO_MDL__
#define __BKGD_EVO_MDL__

#include <math.h>
#include "config.h"
#include "model.h"

#include "bkgd_evo_mdl_param.h"

#define BKGD_BIN_SCALE 1000

#define CLTYPE_MAX_DBL_BRANCH 10


/* use human/chimp data */
#define MDL_TYPE_HC 1
/* use human/chimp/macaque data */    
#define MDL_TYPE_HCM 2  
/* human/chimp/gorilla/orang/macaque */ 
#define MDL_TYPE_HCGOM 3
/* human/chimp/orang/macaque */
#define MDL_TYPE_HCOM 4
/* human/chimp/gorilla/orang */
#define MDL_TYPE_HCGO 5


/* mu is defined by single param for all category values */
#define MU_TYPE_SINGLE 1 
/* mu is a quadratic function of the category value */
#define MU_TYPE_CAT_QUAD 2
/* mu is a linear function of the category value */
#define MU_TYPE_CAT_LIN 3


#define MAX_LINE 32768

#define MAX_BRANCH 100
#define MAX_COL_TYPE 100

#define SUBST_MDL_JUKES_CANTOR 1
#define SUBST_MDL_KIMURA 2


/* basic substitution used for jukes-cantor model */
#define SUBST_TYPE_SUBST 0
/* substitution types used for kimura 2-param model */
#define SUBST_TYPE_TRANSITION 1
#define SUBST_TYPE_TRANSVERSION 2
/* substitution type for conserved (invarient) column */
#define SUBST_TYPE_CONSERVED 3


typedef struct Branch_t Branch;

typedef struct BkgdEvoMdlConfig_t BkgdEvoMdlConfig;


typedef struct {
  double B_ex;   /* exonic B value for bin */
  double B_nex;  /* non-exonic B value for bin */

  double lB_ex;  /* log exonic B value */
  double lB_nex; /* log non-exonic B value */

  long h;    /* count of sites where human differs from other species  */
  long c;    /* count of sites where chimp differs from other species */
  long g;    /* count of sites where gorilla differs from other species */
  long o;    /* count of sites where orang differs from other species */
  long m;    /* count of sites where macaque differs from other species */
  long hc;   /* count of sites where human+chimp differ from other species */
  long hg;   /* count of sites where human+gorilla differ */
  long cg;   /* count of sites where chimp+gorilla differ */
  long ho;   /* count of sites where human+orang differ */
  long co;   /* count of sites where chimp+orang differ */
  long hcg;  /* count of sites where human+chimp+orang differ */  

  /* counts of TRANSITION differences */
  long h_i;
  long c_i;
  long g_i;
  long o_i;
  long m_i;
  long hc_i;
  long hg_i;
  long cg_i;
  long ho_i;
  long co_i;
  long hcg_i;

  /* counts of TRANSVERSION differences */
  long h_v;
  long c_v;
  long g_v;
  long o_v;
  long m_v;
  long hc_v;
  long hg_v;
  long cg_v;
  long ho_v;
  long co_v;
  long hcg_v;

  long cons;  /* number of sites where all species have same base */

  double cat; /* category value (e.g. GC content or HM divergence) */
} BkgdBin;


struct Branch_t {
  int id;           /* identifier of branch, also index in conf branch array */
  char *name;       /* name of branch */
  double len;       /* length of branch in generations */

  double jc;        /* jukes-cantor substitution constant for branch */
  double k_v;       /* kimura transition constant for branch */
  double k_i;       /* kimura transversion constant for branch */

  double prob;     /* probability of observing subst on branch */
  double prob_i;   /* probability of observing transition */
  double prob_v;   /* probability of observing transversion */

  /* partial derivatives of branch len wrt each param */
  BkgdEvoMdlParam dlen; 

  /* partial derivatives of branch subst prob wrt to each param */
  BkgdEvoMdlParam dprob;
  BkgdEvoMdlParam dprob_i;
  BkgdEvoMdlParam dprob_v;
  
  /* ptr to function used to set branch length and probabilities 
   * from parameters and data
   */
  void (*set_len)(Branch *, const BkgdEvoMdlParam *);

  /* ptr to function used to set branch length partial derivatives */
  void (*set_dlen)(Branch *, const BkgdBin *, const BkgdEvoMdlParam *); 
};



typedef struct ColType_t ColType;


/**
 * Structure representing an alignment column type
 * used in the evolutionary bkgd likelihood model
 */
struct ColType_t {
  long n; /* number of sites with this col type */
  char *name; /* name of column type */

  /* ptr to branches that can give this column type 
   * when mutation occurs on them
   */
  Branch *br;

  /* ptr to branch pairs give col type with double mutations */
  int n_dbl_br;
  Branch **dbl_br[CLTYPE_MAX_DBL_BRANCH];

  /* type of substitution (transition, transversion, etc.)*/
  int subst_type;

  /* total prob of observing this col type */
  double prob; 

  /* prob of observing this column due to single-branch subst */
  double prob_single;  

  /* prob of observing this column due to double-branch substs
   * (one for each pair of possible double substs)
   */
  double prob_double[CLTYPE_MAX_DBL_BRANCH];

  /*
   * total prob that this column is due to double-branch substs
   * (summed over all possible double-branch types)
   */
  double prob_ttl_double;

  /* partial derivatives of col-type prob w.r.t. each model param */
  BkgdEvoMdlParam dprob;

  /* ptr to function that sets n attribute given data bin */
  void (*set_n)(ColType *cltype, const BkgdBin *bin);  
};



typedef struct {
  int id;
  char *name;

  int n_subst_type;
  int *subst_type_ids;
  char **subst_type_names;
} SubstModel;


struct BkgdEvoMdlConfig_t {
  int mdl_type;  /* species used in model */

  SubstModel *subst_mdl; /* DNA substitution model */
  int use_double_substs; /* flag indicating whether to model double susbts */

  int mu_type;
  double min_cat; /* minimum category value used for analysis */
  double max_cat; /* maximum category value used for analysis */

  /* gene tree branches associated with model: */
  Branch branches[MAX_BRANCH];

  /* alignment column types associated with model: */
  ColType *cltypes[MAX_COL_TYPE];

  /* points to conserved column type */
  ColType *cons_cltype;

  /* number of branches in model */
  int n_branch;

  /* number of column types in model */
  int n_cltype;
};


void bkgd_evo_mdl_read_data(Model *mdl, Config *config);
void bkgd_evo_mdl_free_data(Model *mdl);
void bkgd_evo_mdl_bin_data(Model *mdl, Config *config);

void set_bkgd_evo_param_from_model(const Model *mdl, 
				   const BkgdEvoMdlConfig *conf,
				   BkgdEvoMdlParam *param);

void bkgd_evo_mdl_set_mdl_deriv(Model *mdl, BkgdEvoMdlParam *deriv);

Model *bkgd_evo_mdl_new(Config *config);
void bkgd_evo_mdl_free(Model *mdl);
double bkgd_evo_mdl_calc_ll(Model *mdl, int calc_ll, int calc_grad);

void bkgd_evo_mdl_calc_mu(BkgdEvoMdlParam *param, Model *mdl, BkgdBin *bin);

#endif
