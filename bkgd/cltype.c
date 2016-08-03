
#include <glib.h>

#include "cltype.h"
#include "bkgd_evo_mdl.h"



/**
 * Allocates mem for new ColType structure.
 */
ColType *cltype_new(const char *name, int subst_type,
		    void (*set_n)(ColType *, const BkgdBin *)) {
  int i;
  ColType *cltype;

  cltype = g_new(ColType, 1);

  cltype->name = g_strdup(name);

  cltype->br = NULL;
  cltype->n_dbl_br = 0;
  cltype->n = 0;
  cltype->prob = 0.0;

  cltype->prob_single = 0.0;
  for(i = 0; i < CLTYPE_MAX_DBL_BRANCH; i++) {
    cltype->prob_double[i] = 0.0;
  }
  cltype->prob_ttl_double = 0.0;

  cltype->set_n = set_n;
  cltype->subst_type = subst_type;
  
  return cltype;
}


/**
 * Frees mem allocated for ColType structure
 */
void cltype_free(ColType *cltype) {
  int i;

  g_free(cltype->name);

  for(i = 0; i < cltype->n_dbl_br; i++) {
    g_free(cltype->dbl_br[i]);
  }

  g_free(cltype);
}



/**
 * Sets the branch on which single subst results in col type.
 * Can be NULL
 */
void cltype_set_br(ColType *cltype, Branch *br) {
  cltype->br = br;
}



/**
 * Adds pair of branches that give col type when mutations
 * occurs on *each* of them.
 */
void cltype_add_dbl_br(ColType *cltype, Branch *br1, Branch *br2) {
  int i;
    
  if(cltype->n_dbl_br >= CLTYPE_MAX_DBL_BRANCH) {
    g_error("cltype_add_dbl_br: maximum number of double branches (%d)"
	    " exceeded", CLTYPE_MAX_DBL_BRANCH);
  }

  cltype->n_dbl_br++;
  i = cltype->n_dbl_br-1;
  
  cltype->dbl_br[i] = g_new(Branch *, 2);
  cltype->dbl_br[i][0] = br1;
  cltype->dbl_br[i][1] = br2;
}


/**
 * Sets the prob of a conserved column type. Probabilities of
 * all non-conserved columns should be set first using the
 * cltype_set_prob function.
 *
 */
void cltype_cons_set_prob(ColType *cltype, const BkgdEvoMdlParam *param,
			  const BkgdEvoMdlConfig *conf) {
  double p;
  int i;
  
  if(cltype->subst_type != SUBST_TYPE_CONSERVED) {
    g_error("cltype_cons_set_prob: expected susbt type "
	    "to be SUBST_TYPE_CONSERVED");
  }

  /* prob of conserved col is 1 - (sum of all other column probs) */
  p = 1.0;
  for(i = 0; i < conf->n_cltype; i++) {
    if(conf->cltypes[i]->subst_type != SUBST_TYPE_CONSERVED) {
      p -= conf->cltypes[i]->prob;
    }
  }

  if(p < 0.0) {
    g_error("cltype_cons_set_prob: probability of "
	      "conserved column is < 0.0 (%g)", p);
    p = 0.0;
  }

  cltype->prob = p;

  return;
}


/**
 * Sets the probability of this column type. Probabilities of branches
 * provided branches must be set before this is called.
 */
void cltype_set_prob(ColType *cltype, const BkgdEvoMdlParam *param,
		     const BkgdEvoMdlConfig *conf) {
  int i, j;
  double p, ttl_p;
  const Branch *br = conf->branches;

  
  /* fprintf(stderr, "col type: %s\n", cltype->name); */
  if(((cltype->br == NULL) && (cltype->n_dbl_br == 0)) ||
     (cltype->subst_type == SUBST_TYPE_CONSERVED)) {
    g_error("cltype_set_prob: column type should be associated "
	    "single or double branch substitutions; cltype_cons_set_prob "
	    "should be used for conserved columns");
  }

  ttl_p = 0.0;

  if(cltype->br) {
    p = 1.0;

    for(i = 0; i < conf->n_branch; i++) {
      if(br[i].id == cltype->br->id) {

	if(conf->subst_mdl->id == SUBST_MDL_JUKES_CANTOR) {
	  p *= br[i].prob / 3.0;
	}
	else if(conf->subst_mdl->id == SUBST_MDL_KIMURA) {
	  if(cltype->subst_type == SUBST_TYPE_TRANSITION) {
	    p *= br[i].prob_i;
	  } 
	  else if(cltype->subst_type == SUBST_TYPE_TRANSVERSION) {
	    /* if we want JUKES and KIMURA to be nested models:
	     *   p *= 0.5 * br[i].prob_v; 
	     */
	    p *= br[i].prob_v;
	  } else {
	    g_error("cltype_set_prob: unknown subst type for column");
	  }
	}
	else {
	  g_error("cltype_set_prob: unknown substitution model");
	}
      } else {
	p *= (1.0 - br[i].prob);
      }
    }

    ttl_p += p;

    /* set prob of observing this col due to single subst */
    cltype->prob_single = p;
  } else {
    cltype->prob_single = 0.0;
  }

  /* calculate prob of observing this col due to double subst */ 
  cltype->prob_ttl_double = 0.0;

  if(cltype->n_dbl_br > 0) {
    for(i = 0; i < cltype->n_dbl_br; i++) {

      p = 0.0;
      if(conf->subst_mdl->id == SUBST_MDL_JUKES_CANTOR) {
	/* Under jukes-cantor model, prob of double subst is product of
	 * two branch probabilities times 1/9
	 */
	p = param->lambda * 
	  (cltype->dbl_br[i][0]->prob * cltype->dbl_br[i][1]->prob)/9.0;
      } else if(conf->subst_mdl->id == SUBST_MDL_KIMURA) {
	switch(cltype->subst_type) {
	case(SUBST_TYPE_TRANSITION):
	  p = param->lambda_i * 
	    (cltype->dbl_br[i][0]->prob_i * cltype->dbl_br[i][1]->prob_i);
	  break;

	case(SUBST_TYPE_TRANSVERSION):
	  /* if we want JUKES and KIMURA to be nested models: */
	  /* p = 0.25 * (cltype->dbl_br[i][0]->prob_v * 
           *              cltype->dbl_br[i][1]->prob_v);
	   */
	  p = param->lambda_v * 0.5 * (cltype->dbl_br[i][0]->prob_v * 
				       cltype->dbl_br[i][1]->prob_v);
	  break;
	default:
	  g_error("cltype_set_prob: unknown substitution type");
	}
      } else {
	g_error("cltype_set_prob: unknown substitution model");
      }

      /* multiply in prob of non-substitution on other branches */
      for(j = 0; j < conf->n_branch; j++) {
	if((br[j].id != cltype->dbl_br[i][0]->id) && 
	   (br[j].id != cltype->dbl_br[i][1]->id)) {
	  p *= (1.0 - br[j].prob);
	}
      }

      cltype->prob_double[i] = p;
      cltype->prob_ttl_double += p;
      ttl_p += p;
    }
  } else {
    cltype->prob_ttl_double = 0.0;
  }

  cltype->prob = ttl_p;
}





static void cltype_set_dprob_kimura(ColType *cltype, 
				    const BkgdEvoMdlParam *param, 
				    const BkgdEvoMdlConfig *conf) {
  double coef;
  int i,j;
  Branch *br1, *br2;
  const Branch *br = conf->branches;

  /* now calculate partial derivs for each column type prob (pi)
   * using branch subst prob partials
   */
  param_set_zero(&cltype->dprob);

  if((cltype->subst_type != SUBST_TYPE_TRANSITION) &&
     (cltype->subst_type != SUBST_TYPE_TRANSVERSION)) {
    g_error("cltype_set_dprob_kimura: unknown subst type for column %s",
	    cltype->name);
  }

  if(cltype->br) {
    /* add in contribution from single-branch substitutions */
    if(cltype->subst_type == SUBST_TYPE_TRANSITION) {
      for(i = 0; i < conf->n_branch; i++) {
	if(br[i].id == cltype->br->id) {
	  coef = cltype->prob_single / (br[i].prob_i);
	  param_add(&br[i].dprob_i, coef, &cltype->dprob);
	} else {
	  coef = cltype->prob_single / (br[i].prob - 1.0);
	  param_add(&br[i].dprob, coef, &cltype->dprob);
	}
      }
    } else {
      /* Transversion */
      for(i = 0; i < conf->n_branch; i++) {
	if(br[i].id == cltype->br->id) {
	  coef = cltype->prob_single / br[i].prob_v;
	  param_add(&br[i].dprob_v, coef, &cltype->dprob);
	} else {
	  coef = cltype->prob_single / (br[i].prob - 1.0);
	  param_add(&br[i].dprob, coef, &cltype->dprob);
	}
      }
    }
  }

  /* add in contributions from double-branch substitutions */
  for(i = 0; i < cltype->n_dbl_br; i++) {
    br1 = cltype->dbl_br[i][0];
    br2 = cltype->dbl_br[i][1];

    if(cltype->subst_type == SUBST_TYPE_TRANSITION) {
      /* transition */
      for(j = 0; j < conf->n_branch; j++) {
	if((br[j].id == br1->id) || (br[j].id == br2->id)) {
	  coef = cltype->prob_double[i] / br[j].prob_i;
	  param_add(&br[j].dprob_i, coef, &cltype->dprob);
	} else {
	  coef = cltype->prob_double[i] / (br[j].prob-1.0);
	  param_add(&br[j].dprob, coef, &cltype->dprob);
	}
      }
    } else {
      /* transversion */
      for(j = 0; j < conf->n_branch; j++) {
	if((br[j].id == br1->id) || (br[j].id == br2->id)) {
	  coef = cltype->prob_double[i] / br[j].prob_v;
	  param_add(&br[j].dprob_v, coef, &cltype->dprob);
	} else {
	  coef = cltype->prob_double[i] / (br[j].prob-1.0);
	  param_add(&br[j].dprob, coef, &cltype->dprob);
	}
      }
    }
  }

  /* partial deriv wrt lambda is just double subst prob term 
   * but dividing out the lambda
   */
  if(cltype->subst_type == SUBST_TYPE_TRANSITION) {
    cltype->dprob.lambda_i = cltype->prob_ttl_double / param->lambda_i;
    cltype->dprob.lambda_v = 0.0;
  } else {
    cltype->dprob.lambda_i = 0.0;
    cltype->dprob.lambda_v = cltype->prob_ttl_double / param->lambda_v;
  }
}



static void cltype_set_dprob_jukes_cantor(ColType *cltype, 
					  const BkgdEvoMdlParam *param, 
					  const BkgdEvoMdlConfig *conf) {
				   
  double coef;
  int i, j;
  Branch *br1, *br2;
  const Branch *br = conf->branches;

  /* now calculate partial derivs for each column type prob (pi)
   * using branch subst prob partials
   */
  param_set_zero(&cltype->dprob);

  if(cltype->br) {
    /* add in contribution from single-branch substitutions */

    for(i = 0; i < conf->n_branch; i++) {
      if(br[i].id == cltype->br->id) {
	coef = cltype->prob_single/ (br[i].prob);
      } else {
	coef = cltype->prob_single / (br[i].prob - 1.0);
      }

      param_add(&br[i].dprob, coef, &cltype->dprob);
    }
  }

  /* add in contributions from double-branch mutations */
  for(i = 0; i < cltype->n_dbl_br; i++) {
    br1 = cltype->dbl_br[i][0];
    br2 = cltype->dbl_br[i][1];

    /* contribution from (1-p) partial terms */
    for(j = 0; j < conf->n_branch; j++) {
      if((br[j].id == br1->id) || (br[j].id == br2->id)) {
	coef = cltype->prob_double[i] / br[j].prob;
      } else {
	coef = cltype->prob_double[i] / (br[j].prob-1.0);
      }
	
      param_add(&br[j].dprob, coef, &cltype->dprob);
    }
  }
  cltype->dprob.lambda = cltype->prob_ttl_double / param->lambda;
}



/**
 * Sets partial derivs of conserved column probability wrt model
 * parameters. Partial derivs for other column types should
 * be set first using cltype_set_dprob function.
 */
void cltype_cons_set_dprob(ColType *cltype, const BkgdEvoMdlParam *param,
			   const BkgdEvoMdlConfig *conf) {
  double coef;
  int i;

  param_set_zero(&cltype->dprob);
  
  if(cltype->subst_type != SUBST_TYPE_CONSERVED) {
    g_error("cltype_cons_set_dprob: expected a conserved col type");
  }

  coef = -1.0;
  for(i = 0; i < conf->n_cltype; i++) {
    if(conf->cltypes[i]->subst_type != SUBST_TYPE_CONSERVED) {
      param_add(&conf->cltypes[i]->dprob, coef, &cltype->dprob);
    }
  }
}


/**
 * Sets first partial derivatives of probabilities for column types
 * with substitutions
 */
void cltype_set_dprob(ColType *cltype, const BkgdEvoMdlParam *param, 
		      const BkgdEvoMdlConfig *conf) {


  if(((cltype->br == NULL) && (cltype->n_dbl_br == 0)) ||
     (cltype->subst_type == SUBST_TYPE_CONSERVED)) {
    /* conserved column type */
    g_error("cltype_set_dprob: conserved column types should be "
	    "passed to cltype_cons_set_dprob instead");
  }

  if(conf->subst_mdl->id == SUBST_MDL_JUKES_CANTOR) {
    cltype_set_dprob_jukes_cantor(cltype, param, conf);
  }
  else if(conf->subst_mdl->id == SUBST_MDL_KIMURA) {
    cltype_set_dprob_kimura(cltype, param, conf);
  } 
  else {
    g_error("cltype_set_dprob: unknown substutiton model");
  }
}



/**********************************
 * column type site counts for HC model
 **********************************/

/* sets number of HC column types in HC model */
void cltype_set_n__hc_h_c(ColType *cltype, const BkgdBin *bin) {

  if(cltype->subst_type == SUBST_TYPE_TRANSITION) {
    cltype->n = bin->h_i + bin->c_i + bin->hg_i + bin->cg_i +
      bin->ho_i + bin->co_i;
  }
  else if(cltype->subst_type == SUBST_TYPE_TRANSVERSION) {
    cltype->n = bin->h_v + bin->c_v + bin->hg_v + bin->cg_v +
      + bin->ho_v + bin->co_v;
  }
  else {
    cltype->n = bin->h + bin->c + bin->hg + bin->cg + bin->ho + bin->co;
  } 
}

/* sets number of CONS column types in HC model */
void cltype_set_n__hc_cons(ColType *cltype, const BkgdBin *bin) {
  cltype->n = bin->g + bin->o + bin->m + bin->hc + bin->m + bin->cons;
}



/*********************************
 * column type site counts for HCM model
 *********************************/

void cltype_set_n__hcm_h_c(ColType *cltype, const BkgdBin *bin) {

  if(cltype->subst_type == SUBST_TYPE_TRANSITION) {
    cltype->n = bin->h_i + bin->c_i + bin->hg_i + bin->cg_i +
      bin->ho_i + bin->co_i;
  }
  else if(cltype->subst_type == SUBST_TYPE_TRANSVERSION) {
    cltype->n = bin->h_v + bin->c_v + bin->hg_v + bin->cg_v +
      bin->ho_v + bin->co_v;
  } else {
    cltype->n = bin->h + bin->c + bin->hg + bin->cg + bin->ho + bin->co;
  }
}

void cltype_set_n__hcm_m(ColType *cltype, const BkgdBin *bin) {
  /* count sites where H+C same, but M different */
  if(cltype->subst_type == SUBST_TYPE_TRANSITION) {
    cltype->n = bin->hc_i + bin->hcg_i + bin->m_i;
  }
  else if(cltype->subst_type == SUBST_TYPE_TRANSVERSION) {
    cltype->n = bin->hc_v + bin->hcg_v + bin->m_v;
  } else {
    cltype->n = bin->hc + bin->hcg + bin->m;
  }
}

void cltype_set_n__hcm_cons(ColType *cltype, const BkgdBin *bin) {
  cltype->n = bin->g + bin->o + bin->cons;
}


/*********************************
 * column type site counts for HCOM model
 *********************************/

void cltype_set_n__hcom_h_c(ColType *cltype, const BkgdBin *bin) {
  switch(cltype->subst_type) {
  case(SUBST_TYPE_TRANSITION):
    cltype->n = bin->h_i + bin->c_i + bin->hg_i + bin->cg_i;
    break;

  case(SUBST_TYPE_TRANSVERSION):
    cltype->n = bin->h_v + bin->c_v + bin->hg_v + bin->cg_v;
    break;

  default:
    cltype->n = bin->h + bin->c + bin->hg + bin->cg;
  }
}


void cltype_set_n__hcom_hc(ColType *cltype, const BkgdBin *bin) {
  switch(cltype->subst_type) {
  case(SUBST_TYPE_TRANSITION):
    cltype->n = bin->hc_i + bin->hcg_i;
    break;

  case(SUBST_TYPE_TRANSVERSION):
    cltype->n = bin->hc_v + bin->hcg_v;
    break;

  default:
    cltype->n = bin->hc + bin->hcg;
  }
}

void cltype_set_n__hcom_o(ColType *cltype, const BkgdBin *bin) {
  switch(cltype->subst_type) {
  case(SUBST_TYPE_TRANSITION):
    cltype->n = bin->o_i;
    break;

  case(SUBST_TYPE_TRANSVERSION):
    cltype->n = bin->o_v;    
    break;

  default:
    cltype->n = bin->o;
  }
}

void cltype_set_n__hcom_m(ColType *cltype, const BkgdBin *bin) {
  switch(cltype->subst_type) {
  case(SUBST_TYPE_TRANSITION):
    cltype->n = bin->m_i;
    break;

  case(SUBST_TYPE_TRANSVERSION):
    cltype->n = bin->m_v;
    break;

  default:
    cltype->n = bin->m;
  }
}

void cltype_set_n__hcom_cons(ColType *cltype, const BkgdBin *bin) {
    cltype->n = bin->g + bin->cons;
}



/************************************
 * column type site counts for HCGOM model
 ************************************/

void cltype_set_n__hcgom_h(ColType *cltype, const BkgdBin *bin) {
  switch(cltype->subst_type) {
  case(SUBST_TYPE_TRANSITION):
    cltype->n = bin->h_i;
    break;

  case(SUBST_TYPE_TRANSVERSION):
    cltype->n = bin->h_v;
    break;

  default:
    cltype->n = bin->h;
  }
}

void cltype_set_n__hcgom_c(ColType *cltype, const BkgdBin *bin) {
  switch(cltype->subst_type) {
  case(SUBST_TYPE_TRANSITION):
    cltype->n = bin->c_i;
    break;
  case(SUBST_TYPE_TRANSVERSION):
    cltype->n = bin->c_v;
    break;

  default:
    cltype->n = bin->c;
  }
}

void cltype_set_n__hcgom_g(ColType *cltype, const BkgdBin *bin) {
  switch(cltype->subst_type) {
  case(SUBST_TYPE_TRANSITION):
    cltype->n = bin->g_i;
    break;

  case(SUBST_TYPE_TRANSVERSION):
    cltype->n = bin->g_v;
    break;

  default:
    cltype->n = bin->g;
  }
}


void cltype_set_n__hcgom_o(ColType *cltype, const BkgdBin *bin) {
  switch(cltype->subst_type) {
  case(SUBST_TYPE_TRANSITION):
    cltype->n = bin->o_i;
    break;

  case(SUBST_TYPE_TRANSVERSION):
    cltype->n = bin->o_v;
    break;

  default:
    cltype->n = bin->o;
  }
}

void cltype_set_n__hcgom_m(ColType *cltype, const BkgdBin *bin) {
  switch(cltype->subst_type) {
  case(SUBST_TYPE_TRANSITION):
    cltype->n = bin->m_i;
    break;

  case(SUBST_TYPE_TRANSVERSION):
    cltype->n = bin->m_v;
    break;

  default:
    cltype->n = bin->m;
  }
}

void cltype_set_n__hcgom_hc(ColType *cltype, const BkgdBin *bin) {
  switch(cltype->subst_type) {
  case(SUBST_TYPE_TRANSITION):
    cltype->n = bin->hc_i;
    break;

  case(SUBST_TYPE_TRANSVERSION):
    cltype->n = bin->hc_v;
    break;

  default:
    cltype->n = bin->hc;
  }
}

void cltype_set_n__hcgom_hcg(ColType *cltype, const BkgdBin *bin) {
  switch(cltype->subst_type) {
  case(SUBST_TYPE_TRANSITION):
    cltype->n = bin->hcg_i;
    break;

  case(SUBST_TYPE_TRANSVERSION):
    cltype->n = bin->hcg_v;
    break;

  default:
    cltype->n = bin->hcg;
  }
}


void cltype_set_n__hcgom_hg_cg(ColType *cltype, const BkgdBin *bin) {
  switch(cltype->subst_type) {
  case(SUBST_TYPE_TRANSITION):
    cltype->n = bin->hg_i + bin->cg_i;
    break;

  case(SUBST_TYPE_TRANSVERSION):
    cltype->n = bin->hg_v + bin->cg_v;
    break;
    
  default:
    cltype->n = bin->hg + bin->cg;

  }
}



void cltype_set_n__hcgom_hg(ColType *cltype, const BkgdBin *bin) {
  switch(cltype->subst_type) {
  case(SUBST_TYPE_TRANSITION):
    cltype->n = bin->hg_i;
    break;

  case(SUBST_TYPE_TRANSVERSION):
    cltype->n = bin->hg_v;
    break;
    
  default:
    cltype->n = bin->hg;
  }
}


void cltype_set_n__hcgom_cg(ColType *cltype, const BkgdBin *bin) {
  switch(cltype->subst_type) {
  case(SUBST_TYPE_TRANSITION):
    cltype->n = bin->cg_i;
    break;

  case(SUBST_TYPE_TRANSVERSION):
    cltype->n = bin->cg_v;
    break;
    
  default:
    cltype->n = bin->cg;
  }
}


void cltype_set_n__hcgom_ho_co(ColType *cltype, const BkgdBin *bin) {
  switch(cltype->subst_type) {
  case(SUBST_TYPE_TRANSITION):
    cltype->n = bin->ho_i + bin->co_i;
    break;

  case(SUBST_TYPE_TRANSVERSION):
    cltype->n = bin->ho_v + bin->co_v;
    break;

  default:
    cltype->n = bin->ho + bin->co;
  }
}

void cltype_set_n__hcgom_cons(ColType *cltype, const BkgdBin *bin) {
  cltype->n = bin->cons;
}

void cltype_set_n__hcgo_cons(ColType *cltype, const BkgdBin *bin) {
  cltype->n = bin->cons + bin->m;
}


