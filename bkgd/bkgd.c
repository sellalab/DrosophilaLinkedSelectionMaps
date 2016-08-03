#include <stdio.h>
#include <math.h>

#include <glib.h>

#include <gsl/gsl_integration.h>

#include "bkgd.h"
#include "bkgd_param.h"
#include "interp_tab.h"


/* amount to scale b values to before rounding to integer */
double BKGD_SCALE = 1000.0;
//#define BKGD_SCALE 1000.0


float gammln(float xx)
{
    /* (C) Copr. 1986-92 Numerical Recipes Software . */
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
    
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}


double bkgd_t_dist_exp(double t, void *v) {
  BkgdParam *p = v;
  return p->exp_c * exp(-p->exp_lambda * t);
}


double bkgd_t_dist_gamma(double t, void *v) {
  BkgdParam *p = v;

  double k = p->gamma_shape;
  double m = p->gamma_scale;
  double tgamma;
  
  tgamma = exp(gammln(k));

  return (p->gamma_c * pow(t, k-1.0) * exp(-t/m)) / (tgamma * pow(m,k));
}


/**
 * Integrand for integration w.r.t. selection coef t, used to
 * calculate contribution of conserved BLOCK to B exp sum. This
 * function uses (analytic) integral approximation to sum over bases
 * in conserved block.
 *
 * This function does not multiply integrand by u/r. This must be done
 * to integral after.
 */
double bkgd_blk_integrand(double t, void *v) {
  BkgdParam *p = v;
  double rho, f;

  rho = (1.0-t)/t;

  if(p->t_dist == NULL) {
    /* point distribution */
    f = (t == p->t) ? 1.0 : 0.0;
  } else {
    f = p->t_dist(t, v);
  }

  if(t == 1.0) {
    /* function value is undefined at t=1.0 */
    return 0.0;
  }
  
  return (f/(1.0-t)) * (1.0/(1.0 + rho*p->r_near) - 1.0/(1.0 + rho*p->r_far));
}


/**
 * Integrand for integration w.r.t. selection coef t. Used to
 * calculate contribution of conserved SITE to B exp sum.
 *
 * This function does not multiply integrand by u. This must be done
 * to integral after.
 */
double bkgd_site_integrand(double t, void *v) {
  BkgdParam *p = v;
  double rho, f;
  double d;

  rho = (1.0-t)/t;

  if(p->t_dist == NULL) {
    /* point distribution */
    f = (t == p->t) ? 1.0 : 0.0;
  } else {
    f = p->t_dist(t, v);
  }

  d = (1.0 + rho * p->r_near);

  if(d == 0.0) {
    return 0.0;
  }

  return f/(t*d*d);
}


/**
 * Function representing the contribution of a conserved block to the
 * B-value 1st derivative sum. This function is integrated over t
 * numerically to calculate the contribution.
 *
 * The expression is not multiplied by (u*B*(delta))/r; this must
 * done outside of the integration.
 */
double bkgd_drv1_blk_integrand(double t, void *parm) {
  BkgdParam *p = parm;
  double rho, f, d1, d2;
  
  rho = (1.0-t)/t;

  if(p->t_dist == NULL) {
    /* point distribution */
    f = (t == p->t) ? 1.0 : 0.0;
  } else {
    f = p->t_dist(t, parm);
  }
  
  d1 = 1.0 + rho * p->r_near;
  d2 = 1.0 + rho * p->r_far;

  if(d1 == 0.0 || d2 == 0.0) {
    return 0.0;
  }

  return f * (1.0/(t*d1*d1) - 1.0/(t*d2*d2));
}


/**
 * Function representing the contribution of a single site to the
 * B-value 1st derivative sum.
 *
 * The expression is not multiplied by (u*B*(delta)) this must done
 * outside of the integration.
 */
double bkgd_drv1_site_integrand(double t, void *parm) {
  BkgdParam *p = parm;
  double d, rho, f;
  
  rho = (1.0-t)/t;

  if(p->t_dist == NULL) {
    /* point distribution */
    f = (t == p->t) ? 1.0 : 0.0;
  } else {
    f = p->t_dist(t, parm);
  }

  d = (1.0 + rho*p->r_near);

  if(d == 0.0) {
    return 0.0;
  }

  return rho*f*(2.0/(t*d*d*d));
}



/**
 * Function representing the contribution of a conserved block to the
 * B-value 2nd derivative sum. This function is integrated over t
 * numerically to calculate the contribution.
 *
 * The expression is not multiplied by (u*B*(delta))/r; this must
 * done outside of the integration.
 */
double bkgd_drv2_blk_integrand(double t, void *parm) {
  BkgdParam *p = parm;
  double rho, f, d1, d2;
  double y;
  
  rho = (1.0-t)/t;

  if(p->t_dist == NULL) {
    /* point distribution */
    f = (t == p->t) ? 1.0 : 0.0;
  } else {
    f = p->t_dist(t, parm);
  }
  
  d1 = 1.0 + rho * p->r_near;
  d2 = 1.0 + rho * p->r_far;

  if(d1 == 0.0 || d2 == 0.0) {
    return 0.0;
  }

  y = f*rho*(2.0/(t*d1*d1*d1) - 2.0/(t*d2*d2*d2));

  
/* fprintf(stderr, "d1=%g, d2=%g, f=%g, rho=%g, t=%g\n", d1, d2, f, rho, t); */
/* fprintf(stderr, "integrand=%g\n", y); */

  return y;
}



/**
 * Function representing the contribution of a single site to the
 * B-value 2nd derivative sum.
 *
 * The expression is not multiplied by (u*B*(delta)) this must done
 * outside of the integration.
 */
double bkgd_drv2_site_integrand(double t, void *parm) {
  BkgdParam *p = parm;
  double d, rho, f;
  
  rho = (1.0-t)/t;

  if(p->t_dist == NULL) {
    /* point distribution */
    f = (t == p->t) ? 1.0 : 0.0;
  } else {
    f = p->t_dist(t, parm);
  }

  d = (1.0 + rho*p->r_near);

  if(d == 0.0) {
    return 0.0;
  }


  return (rho*rho*f*6.0)/(t*d*d*d*d);
}





/**
 * Calculates the contribution of the provided conserved block to the
 * the b-value summation for the provided position.
 */
static inline double get_blk_term(ConsBlock *cblk, BkgdPoint *bpoint,
				  double *drv1_term, double *drv2_term, 
				  BkgdParam *p) {
  double term, r_dist, r_len;
  long len;

  /* add in calculations for 1st and 2nd derivative terms */

  if(cblk->end < bpoint->pos) {
    /* block is to left of current position */
    r_dist = bpoint->r_pos - cblk->r_end;
    r_len  = cblk->r_end - cblk->r_start;

    if(r_len < 0.0) {
      g_error("get_blk_term: left block len (%g) should not be < 0.0", 
	      r_len);
    }
    if(r_dist < 0.0) {
      g_error("get_blk_term: left block dist (%g) should not be < 0.0", 
	      r_dist);
    }

    if(r_len == 0.0) {
      len = cblk->end - cblk->start + 1;
      term = (double)len * interp_tab_lookup_1d(p->intg_tab_site, r_dist);
      if(drv1_term) {
	*drv1_term = len * interp_tab_lookup_1d(p->intg_tab_drv1_site, r_dist);
      }
      if(drv2_term) {
	*drv2_term = len * interp_tab_lookup_1d(p->intg_tab_drv2_site, r_dist);
      }
    } else {
      term = interp_tab_lookup(p->intg_tab_blk, r_dist, r_len) / cblk->r;
      if(drv1_term) {
	*drv1_term = 
	  interp_tab_lookup(p->intg_tab_drv1_blk, r_dist, r_len) / cblk->r;
      }	
      if(drv2_term) {
	*drv2_term =
	    interp_tab_lookup(p->intg_tab_drv2_blk, r_dist, r_len) / cblk->r;
      }
    }
  }
  else if(cblk->start < bpoint->pos) {
    /* current pos is in middle of block, split block into two parts */

    r_dist = 0.0;

    if(cblk->r > 0.0) {
      /* compute contribution from left portion of block */
      r_len  = bpoint->r_pos - cblk->r_start;

      if(r_len < 0.0) {
	g_error("get_blk_term: overlap block left len (%g) "
		"should not be < 0.0", r_len);
      }

      term = interp_tab_lookup(p->intg_tab_blk, r_dist, r_len);
      if(drv1_term) {
	*drv1_term = interp_tab_lookup(p->intg_tab_drv1_blk, r_dist, r_len);
      }
      if(drv2_term) {
	*drv2_term = interp_tab_lookup(p->intg_tab_drv2_blk, r_dist, r_len);
      }
 
      if(cblk->end > bpoint->pos) {
	/* compute contribution from right portion of block */
	r_len = cblk->r_end - bpoint->r_pos;

	if(r_len < 0.0) {
	  g_error("get_blk_term: overlap block right len (%g) "
		  "should not be < 0.0", r_len);
	}

	term += interp_tab_lookup(p->intg_tab_blk, r_dist, r_len);

	if(drv1_term) {
	  /* 1st deriv terms are negative when block > r_pos */
	  *drv1_term -= interp_tab_lookup(p->intg_tab_drv1_blk, r_dist, r_len);
	}
	if(drv2_term) {
	  *drv2_term -= interp_tab_lookup(p->intg_tab_drv2_blk, r_dist, r_len);
	}
      }

      term /= cblk->r;

      if(drv1_term) {
	*drv1_term  = *drv1_term/cblk->r;
      }
      if(drv2_term) {
	*drv2_term = *drv2_term/cblk->r;
      }
    } else {
      /*  no recombination within block: no need to split into two parts */ 
      len = cblk->end - cblk->start + 1;
      term = len * interp_tab_lookup_1d(p->intg_tab_site, r_dist);

      /* no contribution to 1st deriv */
      if(drv1_term) {
	*drv1_term = 0.0;
      }
      if(drv2_term) {
	*drv2_term = 0.0;
      }
    }
  }
  else {    
    /* block is to right of current position */
    r_dist = cblk->r_start - bpoint->r_pos;
    r_len  = cblk->r_end - cblk->r_start;

    if(r_len < 0.0) {
      g_error("get_blk_term: right block len (%g) should not be < 0.0", 
	      r_len);
    }
    if(r_dist < 0.0) {
      fprintf(stderr, "cblk->start=%ld, cblk->end=%ld, "
	      "cblk->r_start=%.10f, cblk->r_end=%.10f\n",
	      cblk->start, cblk->end, cblk->r_start, cblk->r_end);
      fprintf(stderr, "bpoint->pos=%ld, bpoint->r_pos=%.10f\n",
	      bpoint->pos, bpoint->r_pos);

      g_error("get_blk_term: right block dist (%g) should not be < 0.0",
	      r_dist);
    }

    if(r_len == 0.0) {
      len = cblk->end - cblk->start + 1;
      term = len * interp_tab_lookup_1d(p->intg_tab_site, r_dist);
      if(drv1_term) {
	*drv1_term = -len*interp_tab_lookup_1d(p->intg_tab_drv1_site, r_dist);
      }
      if(drv2_term) {
	*drv2_term = -len*interp_tab_lookup_1d(p->intg_tab_drv2_site, r_dist);
      }
    } else {
      term = interp_tab_lookup(p->intg_tab_blk, r_dist, r_len) / cblk->r;
      if(drv1_term) {
	*drv1_term = 
	  -interp_tab_lookup(p->intg_tab_drv1_blk, r_dist, r_len)/cblk->r;
      }
      if(drv2_term) {
	*drv2_term =
	  -interp_tab_lookup(p->intg_tab_drv2_blk, r_dist, r_len)/cblk->r;
      }
    }
  }


  return term;
}



/**
 * Helper function. Calculates the upper bound on remaining
 * contribution to the B value sum
 */
static inline double calc_sum_up_bound(BkgdPoint *bpoint, GList *cur_cons, 
				       BkgdParam *parm) {
  GList *next;
  ConsBlock *cblk;
  long cons_ttl;
  double r_dist;
  double intg;

  if(cur_cons == NULL) {
    return 0.0;
  }

  /* get count of remaining conserved bases */
  cblk = cur_cons->data;

  /* get position of next conserved base */
  if(bpoint->pos > cblk->end) {
    /* we are moving left on chr */
    cons_ttl = cblk->left_ttl;
    next = g_list_previous(cur_cons);
    if(next == NULL) {
      if(cons_ttl != 0) {
	g_error("calc_sum_up_bound: expected 0 conserved bases not %ld",
		cons_ttl);
      }
      return 0.0;
    }
    cblk = next->data;
    r_dist = bpoint->r_pos - cblk->r_end;
  } else {
    /* we are moving right on chr */
    cons_ttl = cblk->right_ttl;
    next = g_list_next(cur_cons);
    if(next == NULL) {
      if(cons_ttl != 0) {
	g_error("calc_sum_up_bound: expected 0 conserved bases not %ld",
		cons_ttl);
      }
      return 0.0;
    }
    cblk = next->data;
    r_dist = cblk->r_start - bpoint->r_pos;

    if(r_dist < 0.0) {
      fprintf(stderr, "bpoint->pos=%ld, bpoint->r_pos=%.10f\n", bpoint->pos,
	      bpoint->r_pos);
      fprintf(stderr, "cblk->start=%ld, cblk->end=%ld, "
	      "cblk->r_start=%.10f, cblk->r_end=%.10f\n", 
	      cblk->start, cblk->end, cblk->r_start, cblk->r_end);

      g_error("calc_sum_up_bound: r_dist %g should not be < 0.0", 
	      r_dist);
    }

  }
  
  /* to calc upper bound on remaining sum assume the next conserved
   * element contains all of the conserved bases and has rec rate 0.
   * I.e. treat it as a single site with a delterious rate c times
   * higher (where c is the number of remaing sites).
   */

  intg = interp_tab_lookup_1d(parm->intg_tab_site, r_dist); 
  
  return (double)cons_ttl * intg;
}



static inline double calc_sum_remaining(BkgdPoint *bpoint, GList *cur_cons,
					double r_chr_len, BkgdParam *parm) {
  ConsBlock *cblk;
  GList *next;
  long c;
  double r_dist, r_len;
  double m, est;

  /* If we assume that conserved sites are evenly distributed along
   * remainder of conserved distance, then we can treat the whole
   * remaining chromosome as a conserved block with a delterious mutation
   * rate that is scaled by c/m (where c is the number of conserved bases
   * and m is the remaining rec dist). 
   *
   * r_near for the block is simply the distance to the next cons
   * block, and r_far is dist to end of chr.
   */

  /* get count of remaining conserved bases */
  cblk = cur_cons->data;

  /* get position of next conserved base */
  if(bpoint->pos > cblk->end) {
    /* we are moving left on chr */
    c = cblk->left_ttl;
    next = g_list_previous(cur_cons);
    if(next == NULL) {
      if(c != 0) {
	g_error("calc_sum_up_bound: expected 0 conserved bases not %ld",c);
      }
      return 0.0;
    }
    cblk = next->data;
    m = cblk->r_end;
    r_dist = bpoint->r_pos - cblk->r_end;
    r_len  = cblk->r_end; /* len is entire chr up to end of blk */
  } else {
    /* we are moving right on chr */
    c = cblk->right_ttl;
    next = g_list_next(cur_cons);
    if(next == NULL) {
      if(c != 0) {
	g_error("calc_sum_up_bound: expected 0 conserved bases not %ld",c);
      }
      return 0.0;
    }
    cblk = next->data;
    m = r_chr_len - cblk->r_start;
    r_dist = cblk->r_start - bpoint->r_pos;
    r_len = r_chr_len - cblk->r_start; /* len is from blk to end of chr */
  }

  if(r_dist < 0.0) {
    g_error("calc_sum_remaining: r_dist (%g) should not be < 0.0", r_dist);
  }
  if(r_len < 0.0) {
    g_error("calc_sum_remaining: r_len (%g) should not be < 0.0", r_len);
  }

  est = ((double)c/m) * interp_tab_lookup(parm->intg_tab_blk, r_dist, r_len);

  /*fprintf(stderr, "remaining sum estimate: %g\n", est);*/

  return est;
}


/**
 * Calculates background selection B, at the position represented
 * by the provided BkgdPoint structure. The B-value and first and
 * second derivatives of B are all set on the provided point.
 *
 * cons_list must point to a list of conserved blocks, and next_cons
 * should point to the element in the list which is the next conserved
 * block (or currently overlapping block) on the chromosome.
 */
void bkgd_calc_b(BkgdPoint *bpoint, GList *cons_list, GList *next_cons, 
		 BkgdParam *parm, double r_chr_len) {
  double sum, b, max_sum;
  double drv1_sum, drv1_term;
  double drv2_sum, drv2_term;
  GList *cur;
  ConsBlock *cblk;
  
  sum = drv1_sum = drv2_sum = 0.0;


  /* first consider conserved blocks to left of current site */
  if(next_cons == NULL) {
    cur = g_list_last(cons_list);
  } else {
    cur = g_list_previous(next_cons);
  }
  while(cur != NULL) {
    cblk = cur->data;

    if(cblk->end > bpoint->pos) {
      g_error("expected pos (%ld) to be >= cblk (%ld-%ld) when moving"
	      " leftwards on chr", bpoint->pos, cblk->start, cblk->end);
    }
    
    sum += get_blk_term(cblk, bpoint, &drv1_term, &drv2_term, parm);

    drv1_sum += drv1_term;
    drv2_sum += drv2_term;

    if(parm->apprx_sum) {
      /* calculate the maximum sum remaining in this direction */
      max_sum = calc_sum_up_bound(bpoint, cur, parm);
      
      if(max_sum < parm->max_sum_thresh) {
	/* contrib from remaining terms is very small.
	 * estimate contribution and terminate sum
	 */
	sum += calc_sum_remaining(bpoint, cur, r_chr_len, parm);
	/* TODO: should add remainder for  1st/2nd derivatives */
	break;
      }
    }

    cur = g_list_previous(cur);
  }

  /* now consider conserved blocks to right of current site */
  cur = next_cons;
  while(cur != NULL) {
    cblk = cur->data;

    if(cblk->end < bpoint->pos) {
      g_error("expected pos (%ld) to be < cblk->end (%ld) when moving"
	      " rigthwards on chr", bpoint->pos, cblk->end);
    }

    sum += get_blk_term(cblk, bpoint, &drv1_term, &drv2_term, parm);
    drv1_sum += drv1_term;
    drv2_sum += drv2_term;


    if(parm->apprx_sum) {
      /* calculate the maximum sum remaining in this direction */
      max_sum = calc_sum_up_bound(bpoint, cur, parm);
      
      if(max_sum < parm->max_sum_thresh) {
	/* contrib from remaining terms is very small.
	 * estimate contribution and terminate sum
	 */
	sum += calc_sum_remaining(bpoint, cur, r_chr_len, parm);
	/* TODO: should add remainder for  1st/2nd derivatives */
	break;
      }
    }

    cur = g_list_next(cur);
  }

  /* calculate b */
  b = exp(-parm->u * sum);
  bpoint->b = b;

  /* calculate first derivative of b */
  bpoint->b_drv1 = 2.0 * parm->u * b * drv1_sum;

  /* calculate second derivative */
  bpoint->b_drv2 = (bpoint->b_drv1 * bpoint->b_drv1)/b - parm->u * b*drv2_sum;

  return;
}



