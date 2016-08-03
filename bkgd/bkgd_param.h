#ifndef __BKGD_PARAM_H__
#define __BKGD_PARAM_H__

#include <math.h>
#include "config.h"
#include "util.h"

#include "interp_tab.h"



#define BKGD_PARAM_EPS_REL 1e-6

/* The range of integration is constrained so that the minimum value
 * that we allow to be passed to exp (when using the exponential t
 * distribution) is defined by this constant. When values get too low
 * (below about 1e-15) we start running into problems with the
 * numerical integration.
 */
#define BKGD_PARAM_EXP_LOW_BOUND -12

#define BKGD_PARAM_DIST_TYPE_POINT 1
#define BKGD_PARAM_DIST_TYPE_EXP 2
#define BKGD_PARAM_DIST_TYPE_GAMMA 3

extern double BKGD_PARAM_MAX_SUM_THRESH;
//#define BKGD_PARAM_MAX_SUM_THRESH 0.001

typedef struct {
  double t;   /* mean value of selection coef distribution */
  double u;   /* deleterious mutation rate */
  
  /* pointer to selcetion coefficient distribution function */
  double (*t_dist)(double t, void *parm);

  /* a and b define truncation of selection coef distribution f(t):
   *   f(t) = C h(t)  for a <= t <= b
   *   f(t) = 0       otherwise
   */
  double a;
  double b;  

  /* parameters used by exp distribution */
  double exp_c; 
  double exp_lambda;

  /* parameters used by gamma distribution */
  double gamma_c;      /* normalizing constant (b/c PDF is truncated) */
  double gamma_scale;  /* scale parameter (theta) */
  double gamma_shape;  /* shape parameter (k) */

  double r_near;
  double r_far;

  /* approximate B sums by terminating summation when remaining
   * contribution is small?
   */
  int apprx_sum;
  double max_sum_thresh;

  /* fast lookup tables for integrals */
  InterpTab *intg_tab_blk;
  InterpTab *intg_tab_site;
  InterpTab *intg_tab_drv1_blk;
  InterpTab *intg_tab_drv1_site;
  InterpTab *intg_tab_drv2_blk;
  InterpTab *intg_tab_drv2_site;
} BkgdParam;



BkgdParam *bkgd_param_new(Config *config);
void bkgd_param_free(BkgdParam *param);

#endif
