#include <glib.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "config.h"
#include "numer.h"
#include "bkgd_interp.h"
#include "bkgd_param.h"
#include "bkgd_point.h"
#include "bkgd.h"
#include "bkgd_intg.h"


/* maximum distance between interpolation points in morgans */
double BKGD_INTERP_MAX_DIST = 0.0001;
//#define BKGD_INTERP_MAX_DIST 1e-4

/* if the bkgd parameter changes by this much between chosen
 * interpolation points, find closer points to use
 */
double BKGD_CHANGE_THRESH = 0.002;
//#define BKGD_CHANGE_THRESH 0.002

/* when doing quadratic extrapolation aim for this much change in B */
double BKGD_OFFSET = 0.001;
//#define BKGD_OFFSET 0.001





/**
 * Returns the location of the next sample that we want to take given
 * the current location by fitting a quadratic at the current
 * position.
 */
static double get_r_next(double x, double bkgd, double b_drv1, double b_drv2) {
  double a, b, c, c_new, discr, sq, x1,x2,x3,x4, closest;

  /* Fit a curve to B at x by finding a with the same value, 1st and
   * 2nd derivates at this point.
   */

  /* find coefficients of quadtratic */
  a = 0.5 * b_drv2;
  b = b_drv1 - b_drv2 * x;
  c = a*x*x - b_drv1*x + bkgd;


  /* now we want to find point where we estimate B will have changed
   * by BKGD_OFFSET offset using the quadratic equation
   */

  /* find next point that is offset GREATER */
  c_new = c - bkgd-BKGD_OFFSET;
  discr = b*b - 4*a*c_new;

  if(discr < 0.0) {
    x1 = x2 = x + BKGD_INTERP_MAX_DIST;
  } else {
    sq = sqrt(discr);
    x1 = (-b + sq)/(2.0*a);
    x2 = (-b - sq)/(2.0*a);
  }
  
  /* find next point that is offset LOWER */
  c_new = c - bkgd+BKGD_OFFSET;

  discr = b*b - 4*a*c_new;

  if(discr < 0.0) {
    x3 = x4 = x + BKGD_INTERP_MAX_DIST;
  } else {
    sq = sqrt(discr);
    x3 = (-b + sq)/(2.0*a);
    x4 = (-b + sq)/(2.0*a);
  }

  /* take closest point which is greater than x */
  closest = x+BKGD_INTERP_MAX_DIST;
  if(x1 > x && x1 < closest) {
    closest = x1;
  }
  if(x2 > x && x2 < closest) {
    closest = x2;
  }
  if(x3 > x && x3 < closest) {
    closest = x3;
  }
  if(x4 > x && x4 < closest) {
    closest = x4;
  }

  
/*   fprintf(stderr, "x1=%g, x2=%g, x3=%g, x4=%g closest=%g\n",  */
/* 	  x1, x2, x3, x4, closest); */

  return closest;
}






/**
 * Decides on what position should be evaluated next using quadratic
 * extrapolation from last_p, and sets the attributes of p
 * appropriately.
 */
static void get_next_point(BkgdInterp *bgi, BkgdPoint *last_p, BkgdPoint *p) {
  double next_r, r_delta, b_delta;
  long pos_delta;
  int keep_point;
  ConsBlock *cblk;
  BkgdPoint *saved_p;

  /* predict a good next r position */
  next_r = get_r_next(last_p->r_pos,  last_p->b, 
		      last_p->b_drv1, last_p->b_drv2);
  r_delta = 0.0;

  /* check the queue of saved positions before evaluating B at new positions */
  if(bgi->p_queue->length > 0) {
    saved_p = g_queue_peek_head(bgi->p_queue);

    if(saved_p->r_pos <= next_r) {
      /* saved position is closer than predicted position, so try
       * it instead
       */

      /* calculate differences in saved position and last position */
      b_delta   = fabs(saved_p->b - last_p->b);
      pos_delta = saved_p->pos - last_p->pos;
      r_delta   = (saved_p->r_pos - last_p->r_pos);
      
      if(pos_delta > 1 && b_delta > BKGD_CHANGE_THRESH && r_delta > 0.0) {
	/* Don't use saved position because still too far away */
	/* instead try a position that is 1/4 as far away */
	next_r = last_p->r_pos + r_delta * 0.25;
	keep_point = FALSE;
      } else {
	/* saved position is good */
	g_queue_pop_head(bgi->p_queue);

/*  	fprintf(stderr, "using saved point at %g, queue len=%d\n", */
/*  		saved_p->r_pos, bgi->p_queue->length); */
/* 	fprintf(stderr, "last_p->b=%g, saved_p->b: %g, " */
/* 		"b_delta=%g, pos_delta=%ld, r_delta=%g\n", */
/* 		last_p->b, saved_p->b, b_delta, pos_delta, r_delta); */

	bkgd_point_copy(saved_p, p);
	g_free(saved_p);
	keep_point = TRUE;
      }
    } else {
      /* predicted position is closer, try it instead */
      keep_point = FALSE;
    }
  } else {
    /* there are no saved positions so try predicted position */
    keep_point = FALSE;
  }

  /* keep trying closer points until there are no closer points,
   * or the change threshold is within our tolerance
   */
  while(!keep_point) {
    p->pos = last_p->pos+1;


    /* find first base that is greater than next r position */
    while((p->pos < bgi->chr_len) && 
	  (rectab_rpos(bgi->rtab,p->pos) < next_r)) {
      p->pos += 1;
    }

    if((p->pos > last_p->pos+1) && 
       (rectab_rpos(bgi->rtab,p->pos) > next_r)) {
      /* backup one base so we don't exceed r pos we were aiming for */
      p->pos -= 1;
    }
  
    p->r_pos = rectab_rpos(bgi->rtab,p->pos);

    /* Update next_cons so that it points to next conserved block.
     * Remember that p can move backwards if we chose a site that
     * was too far ahead. When this happens we may have to backup in
     * the conserved block list instead of moving forwards.
     */
    /* first backup */
    while(bgi->next_cons != NULL) {
      cblk = bgi->next_cons->data;
      if(cblk->start <= p->pos) {
	break;
      }
      bgi->next_cons = g_list_previous(bgi->next_cons);
    }
    if(bgi->next_cons == NULL) {
      bgi->next_cons = bgi->cons_list;
    }
    /* now advance */
    while(bgi->next_cons != NULL) {
      cblk = bgi->next_cons->data;
      
      if(cblk->end >= p->pos) {
	break;
      }
      bgi->next_cons = g_list_next(bgi->next_cons);
    }
  
    /* calculate B value and 1st/2nd derivatives at new position */
    bkgd_calc_b(p, bgi->cons_list, bgi->next_cons, 
		bgi->parm, bgi->rtab->chr_r_len);

    /* calculate differences in position and b between points */
    b_delta = fabs(p->b - last_p->b);
    pos_delta = p->pos - last_p->pos;
    r_delta = (p->r_pos - last_p->r_pos);
   
/*     fprintf(stderr, "B EVAL:p->pos=%ld, p->r_pos=%g,  p->b=%g, " */
/* 	    "last_p->b=%g, pos_delta=%ld,r_delta=%g, b_delta=%g\n", */
/*  	    p->pos, p->r_pos, p->b, last_p->b, pos_delta, r_delta, b_delta); */

    if(pos_delta > 1 && b_delta > BKGD_CHANGE_THRESH && r_delta > 0.0) {
      /* We don't want use this point because it is still too far
       * away.  Push it onto queue for later use.
       */
      g_queue_push_head(bgi->p_queue, bkgd_point_dup(p));
      
      /* Try a position that is 1/4 as far away as last position we
       * aimed for. Do not use r_delta here (which is always less than
       * or equal to the difference we were aiming for) because if
       * there is a jump in rec-dist at a particular site, we can end
       * up trying the same position over and over again.
       */
      next_r = last_p->r_pos + (next_r - last_p->r_pos)*0.25;
    } else {
      /* keep this position */
      keep_point = TRUE;
    }
  }

  /* update slope */
  if(r_delta == 0.0) {
    /* This occurs when two points are immediately adjacent (at least
     * in rec dist). In this case there is no change in the B values
     * between the points, so just call slope 0.
     */    
    bgi->slope = 0.0;
  } else {
    bgi->slope = (p->b - last_p->b) / r_delta;
  }

  // DBG : if(last_p->pos >= 11000000) fprintf(stderr, "last-p-pos = %ld, p-pos = %ld\n", last_p->pos, p->pos); 
}




/**
 * Creates a new BkgdInterp structure, that is intended to
 * transparently decide on appropriate locations to evaluate the B
 * function and to perform interpolation between the locations where B
 * was actually evaluated.
 */
BkgdInterp *bkgd_interp_new(RecRateTable *rtab, long chr_len, 
			    GList *cons_list, BkgdParam *parm) {
  BkgdInterp *bgi;

  bgi = g_new(BkgdInterp, 1);

  bgi->rtab = rtab;
  bgi->chr_len = chr_len;
  bgi->cons_list = cons_list;
  bgi->next_cons = cons_list;
  bgi->parm = parm;

  bgi->p_queue = g_queue_new();

  /* use first base on chr as first point */
  bgi->p1 = g_new(BkgdPoint, 1);
  bgi->p1->pos = 1;
  bgi->p1->r_pos = rectab_rpos(rtab,1);

  bkgd_calc_b(bgi->p1, bgi->cons_list, bgi->next_cons, parm, rtab->chr_r_len);

  bgi->p2 = g_new(BkgdPoint, 1);

  /* compute position of second point */
  get_next_point(bgi, bgi->p1, bgi->p2);

  return bgi;
}



/**
 * Frees memory allocated for BkgdInterp. 
 */
void bkgd_interp_free(BkgdInterp *bgi) {
  BkgdPoint *p;
  g_free(bgi->p1);
  g_free(bgi->p2);

  while(bgi->p_queue->length > 0) {
    p = g_queue_pop_head(bgi->p_queue);
    g_free(p);
  }
  g_queue_free(bgi->p_queue);

  g_free(bgi);
}


/**
 * Calculates an interpolated value at the desired position using the
 * provided initialised interpolator.
 */
double bkgd_interp_eval(BkgdInterp *bgi, long pos) {
  double r_delta;
  BkgdPoint *p;

  if(bgi->p2->pos < pos) {
    /* need to find next point and update slope */
    p = bgi->p1;
    bgi->p1 = bgi->p2; /* second point becomes first point */
    get_next_point(bgi, bgi->p1, p);
    bgi->p2 = p; /* new point becomes second point */
  }

  /* perform linear interpolation (in recombination space) */
  r_delta = rectab_rpos(bgi->rtab, pos) - bgi->p1->r_pos;

/*   fprintf(stderr, "prev=%g, next=%g, slope=%g, r_delta=%g\n", */
/* 	  bgi->p1->b, bgi->p2->b, bgi->slope, r_delta); */

  return bgi->p1->b + r_delta*bgi->slope;
}



