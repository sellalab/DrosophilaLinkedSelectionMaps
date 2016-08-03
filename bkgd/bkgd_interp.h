#ifndef __BKGD_INTERP_H__
#define __BKGD_INTERP_H__

#include <math.h>
#include <glib.h>
#include <gsl/gsl_spline.h>
#include "config.h"
#include "bkgd.h"
#include "bkgd_point.h"
#include "rectab.h"


/* maximum distance between interpolation points in morgans */
extern double BKGD_INTERP_MAX_DIST;
//#define BKGD_INTERP_MAX_DIST 1e-4

/* if the bkgd parameter changes by this much between chosen
 * interpolation points, find closer points to use
 */
extern double BKGD_CHANGE_THRESH;
//#define BKGD_CHANGE_THRESH 0.002

/* when doing quadratic extrapolation aim for this much change in B */
extern double BKGD_OFFSET;
//#define BKGD_OFFSET 0.001

typedef struct {
  long chr_len;    /* number of bases on chr */
  RecRateTable *rtab;  /* recombination map */

  BkgdPoint *p1; /* first point in interpolation region */
  BkgdPoint *p2; /* last point in interpolation region */

  /* stores points that we evaluated that turned out to be too distant */
  GQueue *p_queue; 

  double slope; /* slope of line between two points (in B/rec dist) */

  GList *cons_list;
  GList *next_cons;
  BkgdParam *parm;
} BkgdInterp;


BkgdInterp *bkgd_interp_new(RecRateTable *rtab, long chr_len, 
			    GList *cons_list, BkgdParam *parm);

void bkgd_interp_free(BkgdInterp *bgi);

double bkgd_interp_eval(BkgdInterp *bgi, long pos);


#endif 
