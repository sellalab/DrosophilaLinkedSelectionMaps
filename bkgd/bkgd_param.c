#include <glib.h>
#include <stdio.h>
#include <string.h>

#include "config.h"
#include "bkgd_param.h"
#include "bkgd.h"
#include "interp_tab.h"
#include "bkgd_intg.h"
#include "bkgd_interp.h"

#define BKGD_PARAM_MIN_T_DIST 1e-6
#define BKGD_PARAM_B_STEP 0.001
#define BKGD_PARAM_A_STEP 0.001


double BKGD_PARAM_MAX_SUM_THRESH = 0.001;
//#define BKGD_PARAM_MAX_SUM_THRESH 0.001



void adjust_t_upbound(BkgdParam *p) {
  double f;

  f = p->t_dist(p->b, p);
  while(f < BKGD_PARAM_MIN_T_DIST && p->b > p->a) {
    p->b -= BKGD_PARAM_B_STEP;
    f = p->t_dist(p->b, p);
  }
  if(p->b < p->a) {
    g_error("adjust_t_upbound: could not find reasonable upper bound on t");
  }
}


void adjust_t_lowbound(BkgdParam *p) {
  double f;

  f = p->t_dist(p->a, p);
  while(f < BKGD_PARAM_MIN_T_DIST && p->a < p->b) {
    p->a += BKGD_PARAM_B_STEP;
    f = p->t_dist(p->a, p);
  }
  if(p->a > p->b) {
    g_error("adjust_t_upbound: could not find reasonable lower bound on t");
  }
}



/**
 * The exponential distribution we use is truncated.
 * It is defined as:
 *
 *   f(t) = C e^(-lambda*t)  for a <= t <= b
 *   f(t) = 0                otherwise
 *
 * We need to find a lambda that gives mean t for this truncated
 * distribution, and C that normalizes the distribution so the
 * integral over its range is 1.
 *
 * Use this definition of mean to find value:
 * 
 */
static void calc_exp_dist_param(BkgdParam *parm) {
  double upper, lower, lambda, mean, c, a, b, err;
  double upper_bound;

  err = -1.0;
  a = parm->a;
  b = parm->b;

  fprintf(stderr, "estimating C and lambda for truncated exponential "
	  "distribution:\n  a=%g <= t <= b=%g, mean t=%g\n", a, b, parm->t);

  /* start with 1/t as upper bound, and move it upwards
   * by doubling as needed
   */
  lower = 0.0;
  upper = 1.5/parm->t;
  mean = (1.0 + a * exp(-upper*a) - b*exp(-upper*b))/upper;
  lambda = upper;
  while(mean > parm->t) {
    lower = upper;
    upper *= 2.0;
    mean = (1.0 + a * exp(-upper*a) - b*exp(-upper*b))/upper;
  }

  /* do binary search for value that gives desired mean, taking into
   * account that function is monotonically decreasing with increasing
   * lambda
   */
  while(err < 0.0 || err > mean*BKGD_PARAM_EPS_REL) {
    /* next point is midpoint between bounds */
    lambda = (lower + upper) / 2;
    
    mean = (1.0 + a * exp(-lambda*a) - b*exp(-lambda*b))/lambda;

    err = fabs(mean - parm->t);

    if(mean > parm->t) {
      /* cur value is too low, so need higher lambda, 
       * set midpoint as new lower bound 
       */
      lower = lambda;
    } else {
      upper = lambda;
    }
  }

  c = lambda / (exp(-lambda*a) - exp(-lambda*b));
  parm->exp_lambda = lambda;
  parm->exp_c = c;
  
  fprintf(stderr, "  final estimates:\n  lambda=%g, c=%g, err=%g\n", 
	  lambda, c, err);


  /* integration gets into trouble once exp(-x) drops below minimum
   * representable value on machine, define upper bound on t to
   * account for this
   */
  upper_bound = BKGD_PARAM_EXP_LOW_BOUND / -lambda;
  if(upper_bound < parm->a) {
    g_error("calc_exp_dist_param: cannot integrate exponential\n"
	    "with mean %g and specified range [%g..%g] because\n"
	    "too close to 0. upper bound=%g\n",
	    parm->t, parm->a, parm->b, upper_bound);
  }

  if(upper_bound < parm->b) {
    fprintf(stderr, "  using b=%g as upper bound on range of integration "
	    "instead of %g\n", upper_bound, parm->b);
    parm->b = upper_bound;
  }
}





/*
 * integrand used to calculate expectation of truncated gamma distr
 */
double gamma_expect_integrand(double x, void *v) {
  double f;

  f = bkgd_t_dist_gamma(x, v);
  return x * f;
}


/**
 * Used to find scale parameter for truncated gamma distribution.
 * This function will return 0, when scale parameter gives desired
 * mean.
 */
double gamma_root_func(double scale, void *v) {
  BkgdIntg *bi = v;
  double expect, err;

  bi->p->gamma_scale = scale;
  
  /* calculate mean */


  bkgd_gsl_integration_wrapper(&bi->f, bi->p->a, bi->p->b, bi->w,
			       bi->p, &expect, &err);


  /* will be 0 when mean is t and pdf integral is 1 */
  return expect - bi->p->t;
}




static void calc_gamma_dist_param(BkgdParam *parm) {
  BkgdIntg b_intg;
  double low_brk, up_brk, mid, y, err, intgl;

  parm->gamma_c = 1.0;

  b_intg.w = gsl_integration_workspace_alloc(BKGD_INTG_LIMIT);
  b_intg.f.params = parm;
  b_intg.p = parm;
  b_intg.f.function = &gamma_expect_integrand;

  fprintf(stderr, "estimating scale for truncated gamma "
	  "distribution:\n  shape=%g, a=%g <= t <= b=%g, mean t=%g\n",
	  parm->gamma_shape, parm->a, parm->b, parm->t);

  /* fit scale parameter by searching for root of function */

  up_brk = 1.0;
  low_brk = 1e-4;
  
  /* find lower bracket that is on left side of root */

  /* fprintf(stderr, "finding low bracket\nlow_brk=%g\n", low_brk); */

  y = gamma_root_func(low_brk, &b_intg);
  while(y > 0.0) {
    low_brk *= 0.5;
    y = gamma_root_func(low_brk, &b_intg);
    /* fprintf(stderr, "low_brk=%g, y=%g\n", low_brk, y); */
  }

  /* find upper bracket that is on right side of root */
  /* fprintf(stderr, "finding upper bracket\n"); */
  y = gamma_root_func(up_brk, &b_intg);
  while(y < 0.0) {
    up_brk *= 2.0;
    y = gamma_root_func(up_brk, &b_intg);
    /* fprintf(stderr, "up_brk=%g, y=%g\n", up_brk, y); */
  }

  /* fprintf(stderr, "low_brk=%g, up_brk=%g\n", low_brk, up_brk); */

  /* perform binary search for root, halving size of bracket each time */
  err = up_brk - low_brk;
  mid = (low_brk + up_brk) * 0.5;
  while(err > mid*BKGD_PARAM_EPS_REL) {
    y = gamma_root_func(mid, &b_intg);
        
    if(y < 0.0) {
      low_brk = mid;
    } else {
      up_brk = mid;
    }

    mid = (low_brk + up_brk) * 0.5;
    err = up_brk - low_brk;
    /* fprintf(stderr, "low_brk=%g, up_brk=%g, err=%g\n",
     * low_brk, up_brk, err); */
  }

  fprintf(stderr, "  gamma_scale=%g, err=%g\n", parm->gamma_scale, err);

  /* now find normalizing constant to ensure PDF integral = 1 despite
   * truncation
   */
  b_intg.f.function = &bkgd_t_dist_gamma;  

  bkgd_gsl_integration_wrapper(&b_intg.f, parm->a, parm->b, b_intg.w,
			       parm, &intgl, &err);

/*   gsl_integration_qags(&b_intg.f, parm->a, parm->b, */
/* 		       BKGD_INTG_EPS_ABS, BKGD_INTG_EPS_REL, BKGD_INTG_LIMIT, */
/* 		       b_intg.w, &intgl, &err); */
  
  parm->gamma_c = 1.0 / intgl;

  fprintf(stderr, "  gamma_c=%g\n", parm->gamma_c);
  
  gsl_integration_workspace_free(b_intg.w);


  /* integration gets into trouble when f(t) drops below min
   * representable value on machine, define upper bound on t to
   * account for this
   *
   * This is commented out because it is now handled by 
   * the integration routine in a better way.
   */
  /*   adjust_t_upbound(parm); */
  /*   adjust_t_lowbound(parm); */

  /*   fprintf(stderr, "range of integration is now [a=%g..b=%g]\n",  */
  /* 	  parm->a, parm->b); */

}







static void create_intg_tabs(BkgdParam *parm) {
  BkgdIntg b_intg;

  b_intg.w = gsl_integration_workspace_alloc(BKGD_INTG_LIMIT);
  b_intg.f.params = parm;
  b_intg.p = parm;

  fprintf(stderr, "creating interp tab for blk 2nd deriv integrals\n");
  b_intg.f.function = &bkgd_drv2_blk_integrand;
  parm->intg_tab_drv2_blk = 
    interp_tab_create(&bkgd_calc_blk_integral, &b_intg);

  fprintf(stderr, "creating interp tab for site 2nd deriv integrals\n");
  b_intg.f.function = &bkgd_drv2_site_integrand;
  parm->intg_tab_drv2_site = 
    interp_tab_create_1d(&bkgd_calc_site_integral, &b_intg);


  fprintf(stderr, "creating interp tab for blk integrals\n");
  b_intg.f.function = &bkgd_blk_integrand;
  parm->intg_tab_blk = interp_tab_create(&bkgd_calc_blk_integral, &b_intg);

  fprintf(stderr, "creating interp tab for site integrals\n");
  b_intg.f.function = &bkgd_site_integrand;
  parm->intg_tab_site = interp_tab_create_1d(&bkgd_calc_site_integral, 
					     &b_intg);

  fprintf(stderr, "creating interp tab for blk 1st deriv integrals\n");
  b_intg.f.function = &bkgd_drv1_blk_integrand;
  parm->intg_tab_drv1_blk = 
    interp_tab_create(&bkgd_calc_blk_integral, &b_intg);

  fprintf(stderr, "creating interp tab for site 1st deriv integrals\n");
  b_intg.f.function = &bkgd_drv1_site_integrand;
  parm->intg_tab_drv1_site = 
    interp_tab_create_1d(&bkgd_calc_site_integral, &b_intg);

  gsl_integration_workspace_free(b_intg.w);
}



/**
 * Retrieves parameters that are used for calculating B values from
 * the config and creates a new BkgdParam structure.
 */
BkgdParam *bkgd_param_new(Config *config) {
  BkgdParam *parm;
  char *dist_type;
  double local_BKGD_SCALE, local_BKGD_INTERP_MAX_DIST, local_BKGD_CHANGE_THRESH, local_BKGD_OFFSET, local_BKGD_PARAM_MAX_SUM_THRESH;

  parm = g_new(BkgdParam,1);



  /* load approximation parameters. if they are not present in the file, default values are used and a warning is issued */
  local_BKGD_SCALE = config_get_double(config, "BKGD_SCALE");
  if(local_BKGD_SCALE != 0.0)
	  BKGD_SCALE = local_BKGD_SCALE;

  local_BKGD_INTERP_MAX_DIST = config_get_double(config, "BKGD_INTERP_MAX_DIST");
  if(local_BKGD_INTERP_MAX_DIST != 0.0)
	  BKGD_INTERP_MAX_DIST = local_BKGD_INTERP_MAX_DIST;

  local_BKGD_CHANGE_THRESH = config_get_double(config, "BKGD_CHANGE_THRESH");
  if(local_BKGD_CHANGE_THRESH != 0.0)
	  BKGD_CHANGE_THRESH = local_BKGD_CHANGE_THRESH;

  local_BKGD_OFFSET = config_get_double(config, "BKGD_OFFSET");
  if(local_BKGD_OFFSET != 0.0)
	  BKGD_OFFSET = local_BKGD_OFFSET;

  local_BKGD_PARAM_MAX_SUM_THRESH = config_get_double(config, "BKGD_PARAM_MAX_SUM_THRESH");
  if(local_BKGD_PARAM_MAX_SUM_THRESH != 0.0)
	  BKGD_PARAM_MAX_SUM_THRESH = local_BKGD_PARAM_MAX_SUM_THRESH;
  
  fprintf(stderr, "BKGD_SCALE=%g\nBKGD_INTERP_MAX_DIST=%g\nBKGD_CHANGE_THRESH=%g\nBKGD_OFFSET=%g\n", BKGD_SCALE, BKGD_INTERP_MAX_DIST, BKGD_CHANGE_THRESH, BKGD_OFFSET);


  /* t=sh is selection coefficient times heterozygosity coef */
  parm->t = config_get_double(config, "PARAM_T");
  fprintf(stderr, "t=%g\n", parm->t);

  /* u is deleterious mutation rate */
  parm->u = config_get_double(config, "PARAM_U");

  fprintf(stderr, "u=%g\n", parm->u);


  if(config_get_boolean(config, "USE_SUM_APPROXIMATION")) {
    parm->apprx_sum = TRUE;
    /* convert B contrib thresh to sum thresh */
    parm->max_sum_thresh = log(1.0-BKGD_PARAM_MAX_SUM_THRESH) / -parm->u;
    fprintf(stderr, "using approx sum upper bound threshold=%g\n",
	    parm->max_sum_thresh);
  } else {
    parm->apprx_sum = FALSE;
    parm->max_sum_thresh = 0;
  }



  dist_type = config_get_str(config, "PARAM_T_DIST_TYPE");


  if(strcmp(dist_type, "POINT")==0) {
    parm->a = 0;
    parm->b = 0;
    parm->t_dist = NULL;
  }
  else if(strcmp(dist_type, "EXPONENTIAL")==0) {
    parm->t_dist = &bkgd_t_dist_exp;
    parm->a = config_get_double(config, "PARAM_T_DIST_TRUNC_LOW");
    parm->b = config_get_double(config, "PARAM_T_DIST_TRUNC_HIGH");
    calc_exp_dist_param(parm);
  }
  else if(strcmp(dist_type, "GAMMA")==0) {
    parm->t_dist = &bkgd_t_dist_gamma;
    parm->a = config_get_double(config, "PARAM_T_DIST_TRUNC_LOW");
    parm->b = config_get_double(config, "PARAM_T_DIST_TRUNC_HIGH");

    parm->gamma_shape = config_get_double(config, "PARAM_GAMMA_SHAPE");

    calc_gamma_dist_param(parm);
  }
  else {
    g_error("get_bkgd_param: PARAM_T_DIST_TYPE must be one of "
	    "POINT, EXPONENTIAL, GAMMA");
  }
  

  create_intg_tabs(parm);

  return parm;
}


/**
 * Frees memory associated with provided structure
 */
void bkgd_param_free(BkgdParam *param) {
  /* free interpolation tables if they are defined */
  if(param->intg_tab_blk) {
    interp_tab_free(param->intg_tab_blk);
  }
  if(param->intg_tab_site) {
    interp_tab_free(param->intg_tab_site);
  }
  if(param->intg_tab_drv1_blk) {
    interp_tab_free(param->intg_tab_drv1_blk);
  }
  if(param->intg_tab_drv1_site) {
    interp_tab_free(param->intg_tab_drv1_site);
  }
  if(param->intg_tab_drv2_blk) {
    interp_tab_free(param->intg_tab_drv2_blk);
  }
  if(param->intg_tab_drv2_site) {
    interp_tab_free(param->intg_tab_drv2_site);
  }

  g_free(param);
}
