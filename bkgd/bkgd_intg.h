#ifndef __BKGD_INTG_H__
#define __BKGD_INTG_H__

#include <math.h>
#include <gsl/gsl_integration.h>

#include "bkgd_intg.h"
#include "bkgd_param.h"


/* constants for numerical integration */
#define BKGD_INTG_EPS_REL 1e-3
#define BKGD_INTG_EPS_ABS 0
#define BKGD_INTG_LIMIT 1000


typedef struct {
  BkgdParam *p;
  gsl_function f;
  gsl_integration_workspace *w;
} BkgdIntg;


double bkgd_calc_blk_integral(double r_dist, double r_len, void *parm);
double bkgd_calc_site_integral(double r_dist, void *parm);


void bkgd_gsl_integration_wrapper(gsl_function *f, double a, double b,
				  gsl_integration_workspace *w, void *param, 
				  double *intg_ptr, double *err_ptr);


#endif
