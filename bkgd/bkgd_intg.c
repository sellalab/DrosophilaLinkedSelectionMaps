#include <gsl/gsl_integration.h>
#include <stdio.h>

#include "bkgd_intg.h"
#include "bkgd_param.h"

#define INTG_SCAN_SZ 100
#define INTG_DIFF_THRESH 1e-6

#define INTG_MIN 1e-50
#define INTG_LOW 1e-6
#define INTG_HIGH 1e6

#define INTG_RANGE_HIGH 1
#define INTG_RANGE_MID 2
#define INTG_RANGE_LOW 3
#define INTG_RANGE_ZERO 4

typedef struct {
  double a;
  double b;
  double a_y;
  double b_y;

  int type; /* HIGH, MID, LOW or ZERO */

} BkgdIntgRange;


/**
 * frees provided list of integration ranges
 */
static void free_intg_ranges(GSList *ranges) {
  GSList *cur;

  cur = ranges;
  while(cur != NULL) {
    g_free(cur->data);
    cur = g_slist_next(cur);
  }

  if(ranges != NULL) {
    g_slist_free(ranges);
  }
}




int intg_range_type(double val) {
  if(val < INTG_MIN || isnan(val)) {
    return INTG_RANGE_ZERO;
  }
  if(val < INTG_LOW) {
    return INTG_RANGE_LOW;
  }
  if(val < INTG_HIGH) {
    return INTG_RANGE_MID;
  }

  return INTG_RANGE_HIGH;
}



/**
 * Breaks the range of integration down into sub-regions with high-,
 * mid- and low- and "near zero" function values. This is so that
 * numerical integration can be performed separately on each of the
 * regions reducing the possibilty of round-off error
 * occuring. 
 * 
 * Assumes that function is positive over entire range of
 * integration.
 *
 * May not work well for especially tricky functions where there are
 * oscillations between low to high function values in a short range.
 */
static GSList *find_intg_ranges(gsl_function *f, double a, double b, 
				void *param) {

  BkgdIntgRange *r;
  GSList *ranges;
  double delta, x, y, mid, start, end, mid_y, start_y;
  int new_type, mid_type;

  ranges = NULL;
  r = g_new(BkgdIntgRange,1);
  r->a = a;
  r->b = a;

  /* start a range of integration at the start */
  r->a_y = f->function(r->a, param);

  r->type = intg_range_type(r->a_y);
          
  delta = (b - a) / (double)INTG_SCAN_SZ;

  /* look for areas where integration will run into problems
   * and create separate regions for these
   */
  x = a;

  while(x < b) {
    x = x + delta;
    if(x > b) {
      x = b;
    }

    y = f->function(x, param);

    new_type = intg_range_type(y);
    if(r->type == new_type) {
      /* next test point is still of same type so extend current range */
      r->b = x;
      r->b_y = y;
    } else {
      /* search for transition point between ranges */
      start = r->b;
      start_y = r->b_y;
      end = x;
	
      mid = x;
      mid_y = y;
      while((end - start) > INTG_DIFF_THRESH) {
	mid = (start+end)/2.0;
	mid_y = f->function(mid, param);
	mid_type = intg_range_type(mid_y);

	if(mid_type == r->type) {
	  /* mid point is in range of current type */
	  start = mid;
	  start_y = mid_y;
	} else {
	  /* mid point is still too far away */
	  new_type = mid_type;
	  end = mid;
	}
      }

      /* end old range  */
      r->b = start;
      r->b_y = start_y;
      ranges = g_slist_append(ranges, r);

      /* start new range */
      r = g_new(BkgdIntgRange, 1);
      r->a = start;
      r->b = start;
      r->a_y = start_y;
      r->b_y = start_y;
      r->type = new_type;

      x = start;
    }
  }

  /* add last valid range to list of ranges */
  if(r != NULL) {
    if(r->a == r->b) {
      g_free(r);
    } else {
      r->b_y = f->function(r->b, param);

      ranges = g_slist_append(ranges, r);
    }
  }

  return ranges;
}






void bkgd_gsl_integration_wrapper(gsl_function *f, double a, double b,
				  gsl_integration_workspace *w,
				  void *param, 
				  double *intg_ptr, double *err_ptr) {

  GSList *ranges, *cur;
  BkgdIntgRange *r;
  double intg, intg_ttl, err_ttl, err = 0;

  /* find set of ranges over which numerical integration can be performed */
  ranges = find_intg_ranges(f, a, b, param); 
    
  cur = ranges;

  intg_ttl = 0.0;
  err_ttl = 0.0;

  while(cur != NULL) {
    r = cur->data;

    /* integrate numerically over range t=[a..b]*/

/*     switch(r->type) { */
/*     case(INTG_RANGE_ZERO): */
/*       fprintf(stderr, "  zero-valued range:"); */
/*       break; */
/*     case(INTG_RANGE_LOW): */
/*       fprintf(stderr, "  low-valued range:"); */
/*       break; */
/*     case(INTG_RANGE_MID): */
/*       fprintf(stderr, "  mid-valued range:"); */
/*       break; */
/*     case(INTG_RANGE_HIGH): */
/*       fprintf(stderr, "  high-valued range:"); */
/*       break; */
/*     default: */
/*       g_error("unknown range type"); */
/*     } */
/*     fprintf(stderr, "[%g,%g], edge vals:[%g,%g]\n", */
/* 	    r->a, r->b, r->a_y, r->b_y); */

    if(r->type == INTG_RANGE_ZERO) {
      /* skip ranges that are very near zero */
      intg = 0.0;
    } else {
      gsl_integration_qags(f, r->a, r->b, BKGD_INTG_EPS_ABS,
			   BKGD_INTG_EPS_REL, BKGD_INTG_LIMIT,
			   w, &intg, &err);
    }
    /* fprintf(stderr,"    intg area=%g\n", intg);*/

    intg_ttl += intg;
    err_ttl += err;
    cur = g_slist_next(cur);
  }

  free_intg_ranges(ranges);  

  /* fprintf(stderr, "done intg=%g\n", intg_ttl); */

  *intg_ptr = intg_ttl;
  *err_ptr = err_ttl;
}




/**
 * Numerically evaluates an integral for a cons BLOCK of length and
 * distance r_len and r_dist (in morgans). A BkgdIntg argument must be
 * provided that specifies the function to be integrated and the
 * relevant parameters.
 *
 * This function is essentially a wrapper that is used to populate an
 * interpolation table for rapid lookup of integrals.
 *
 */
double bkgd_calc_blk_integral(double r_dist, double r_len, void *parm) {
  BkgdIntg *bi = parm;
  double intg, err;

  bi->p->r_near = r_dist;
  bi->p->r_far = r_dist + r_len;

  /*  fprintf(stderr, "r_dist=%.20f, r_len=%g\n", r_dist, r_len);*/

  if(bi->p->t_dist) {
    bkgd_gsl_integration_wrapper(&bi->f, bi->p->a, bi->p->b, 
				 bi->w, bi->p, &intg, &err);

  } else {
    /* No integration needed since we are using point distribution,
     * Just call integrand with t value.
     */
    intg = bi->f.function(bi->p->t, bi->p);
  }
  

  return intg;
}




/**
 * Numerically evaluates an integral for a cons SITE at distance
 * r_dist (in morgans). A BkgdIntg argument must be provided that
 * specifies the function to be integrated and the relevant
 * parameters.
 *
 * This function is essentially a wrapper that is used to populate an
 * interpolation table for rapid lookup of integrals.
 *
 */
double bkgd_calc_site_integral(double r_dist, void *parm) {
  BkgdIntg *bi = parm;
  double intg, err;

  bi->p->r_near = bi->p->r_far = r_dist;

  if(bi->p->t_dist) {
    bkgd_gsl_integration_wrapper(&bi->f, bi->p->a, bi->p->b, 
				 bi->w, bi->p, &intg, &err);

  } else {
    /* No integration needed since we are using point distribution,
     * Just call integrand with t value.
     */
    intg = bi->f.function(bi->p->t, bi->p);
  }
  
  return intg;
}
