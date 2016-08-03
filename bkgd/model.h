#ifndef __MODEL_H__
#define __MODEL_H__

#include <math.h>
#include <glib.h>
#include <stdio.h>

#define MODEL_INIT_SZ 20


typedef struct Model_t Model;


typedef struct {
  /* name of parameter */
  char *name;

  /* current value of paramter */
  double val;
  
  /* partial derivative of ll function wrt this param  */
  double deriv;

  /* scaling coefficient that parameter value and deriv are multiplied
   * by before being passed to minimizer in order to put parameters in
   * reasonably similar range
   */
  double scale;

  /* TRUE if parameter has upper bound, FALSE otherwise */
  int up_bound;

  /* upper bound on parameter value */
  double max_val;

  /* TRUE if parameter is allowed to be negative */
  int allow_neg;

  /* TRUE if parameter is locked, FALSE otherwise */
  int locked;
  
} ModelParam;


struct Model_t {
  /* ptr to function that calculates log-likelihood */
  double (*llhood) (Model *m, int, int);

  /* array of parameters, n_param long */
  ModelParam *param;

  /* current size of parameter array (is kept >= number of param) */
  int param_sz;

  /* number of parameters */
  int n_param;
  
  /* number of parameters that are currently unlocked:
   * value between 0 and n_param
   */
  int n_unlocked;

  /* ptr to data for model */
  void *data;

  /* ptr to model configuration information, can be NULL */
  void *config;

  /* table that maps parameter names to indices */
  GHashTable *name2idx;
};



Model *model_new();
void model_free(Model *mdl);

ModelParam *model_get_param(const Model *mdl, const char *name);

double model_get_param_val(const Model *mdl, const char *name);
double model_get_param_deriv(const Model *mdl, const char *name);

void model_set_param_val(Model *mdl, const char *name, const double val);
void model_set_param_deriv(Model *mdl, const char *name, const double deriv);

void model_write_param(const Model *mdl, FILE *fh);
void model_write_param_ln(const Model *mdl, FILE *fh);

void model_write_grad(const Model *mdl, FILE *fh);
void model_write_grad_ln(const Model *mdl, FILE *fh);

ModelParam *model_add_param(Model *mdl, const char *name,
			    const int is_locked, const int allow_neg,
			    const double scale);

void model_lock_param(Model *mdl, const char *param_name);

GList *model_get_rand_unlocked_param(const Model *m, int n);

#endif
