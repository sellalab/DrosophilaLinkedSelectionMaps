#ifndef __BKGD_H__
#define __BKGD_H__

//#include <glib.h>
#include <math.h>
#include "bkgd_param.h"
#include "bkgd_point.h"




/* amount to scale b values to before rounding to integer */
extern double BKGD_SCALE;
//#define BKGD_SCALE 1000.0



typedef struct {
  long start;
  long end;

  long left_ttl; /* how many conserved sites are to left of this block */
  long right_ttl; /* how many conserved sites are to right of this block */

  double r; /* rec-rate for block */
  double r_start;
  double r_end;
} ConsBlock;


double bkgd_t_dist_exp(double t, void *v);
double bkgd_t_dist_gamma(double t, void *v);

double bkgd_blk_integrand(double t, void *v);
double bkgd_site_integrand(double t, void *v);
double bkgd_drv1_blk_integrand(double t, void *parm);
double bkgd_drv1_site_integrand(double t, void *parm);
double bkgd_drv2_blk_integrand(double t, void *parm);
double bkgd_drv2_site_integrand(double t, void *parm);


void bkgd_calc_b(BkgdPoint *bpoint, GList *cons_list, 
		 GList *next_cons, BkgdParam *parm, 
		 double r_chr_len);



#endif
