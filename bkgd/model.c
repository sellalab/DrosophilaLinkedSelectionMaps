#include <math.h>
#include <glib.h>
#include <stdio.h>
#include <stdlib.h>

#include "util.h"

#include "model.h"


/**
 * Creates a new model structure.  llhood, grad, data and config
 * attributes are all initialized to NULL.
 */
Model *model_new() {
  Model *mdl;
  
  mdl = g_new(Model, 1);

  mdl->llhood = NULL;

  mdl->param_sz = MODEL_INIT_SZ;
  mdl->param = g_new(ModelParam, mdl->param_sz);

  /* initialize partial derivative functions to NULL */

  mdl->n_param = 0;
  mdl->n_unlocked = 0;

  mdl->data = NULL;
  mdl->config = NULL;

  mdl->name2idx = g_hash_table_new(g_str_hash, g_str_equal);
  
  return mdl;
}


/**
 * Frees memory for provided model and associated parameters. Does not
 * free mem associated with data or config attributes.
 */
void model_free(Model *mdl) {
  int i;
  for(i = 0; i < mdl->n_param; i++) {
    g_free(mdl->param[i].name);
  }
  g_free(mdl->param);
  util_hash_table_free(mdl->name2idx);
  g_free(mdl);
}



/**
 * Retrieves a pointer to the ModelParam object with the provided
 * name.  Raises error if the paramter does not exist
 */
ModelParam *model_get_param(const Model *mdl, const char *name) {
  int *param_idx;

  param_idx = g_hash_table_lookup(mdl->name2idx, name);

  if(param_idx == NULL) {
    g_error("model_get_param: parameter with name '%s' not found", name);
    return NULL;
  }

  return &mdl->param[*param_idx];
}


/**
 * Adds a new parameter with the given attributes to the model.  By
 * default the parameter is assigned a value of 0.0. To set a
 * different value, model_set_param_value should be used. Returns a
 * ptr to the param structure.
 */
ModelParam *model_add_param(Model *mdl, const char *name,
			    const int is_locked,
			    const int allow_neg,
			    const double scale) {
  ModelParam *param;
  int *idx;

  if(mdl->n_param >= mdl->param_sz) {
    /* we are out of mem for parameters, allocate more  */
    mdl->param_sz *= 2; /* double size of array */
    mdl->param = g_renew(ModelParam, mdl->param, mdl->param_sz);
  }


  /* store index of new paramter in hash table, checking to make sure
   * it doesn't already exist
   */  
  idx = g_hash_table_lookup(mdl->name2idx, name);
  if(idx != NULL) {
    g_error("model_add_param: parameter with name '%s' already exists",
	    name);
  }
  idx = g_new(int, 1);
  *idx = mdl->n_param;
  g_hash_table_insert(mdl->name2idx, g_strdup(name), idx);

  /* initialize parameter */
  param = &mdl->param[*idx];
  param->name = g_strdup(name);
  param->locked = is_locked;
  param->up_bound = FALSE;
  param->max_val = 0.0;
  param->deriv = 0.0;
  param->allow_neg = allow_neg;

  param->scale = scale;

  /* update related model attributes */
  mdl->n_param += 1;  
  if(!is_locked) {
    mdl->n_unlocked++;
  }

  return &mdl->param[*idx];
}



/**
 * Retrieves the current parameter value for the parameter with the
 * given name.
 */
double model_get_param_val(const Model *mdl, const char *name) {
  ModelParam *param;
  param = model_get_param(mdl, name);
  return param->val;
}


/**
 * sets the parameter value for the parameter with the given name.
 */
void model_set_param_val(Model *mdl, const char *name, const double val) {
  ModelParam *param;
  param = model_get_param(mdl, name);

  if(isnan(val) || isinf(val)) {
    g_error("model_set_param_val: invalid param value: %s=%g",
	    name, val);
  }

  param->val = val;
}



/**
 * Retrieves the derivative for the parameter with the given name.
 */
double model_get_param_deriv(const Model *mdl, const char *name) {
  ModelParam *param;
  param = model_get_param(mdl, name);
  return param->deriv;
}



/**
 * sets the value of the partial derivative of the LL function 
 * wrt the named paramter 
 */
void model_set_param_deriv(Model *mdl, const char *name, const double deriv) {
  ModelParam *param;
  param = model_get_param(mdl, name);


  if(isnan(deriv) || isinf(deriv)) {
    g_error("model_set_param_deriv: invalid param deriv: d%s=%g",
	    name, deriv);
  }


  param->deriv = deriv;
}


/**
 * Writes string of paramters and values to provided filehandle. Does
 * not write trailing newline.
 */
void model_write_param(const Model *mdl, FILE *fh) {
  int i;
  char sep[2];

  sep[0] = '\0';
  sep[1] = '\0';

  for(i = 0; i < mdl->n_param; i++) {
    if(mdl->param[i].locked) {
      fprintf(fh, "%s*", sep);
    } else {
      fprintf(fh, sep);
    }

    fprintf(fh, "%s=%g", mdl->param[i].name, mdl->param[i].val);
    sep[0] = ' ';
  }
}



/**
 * Locks the specified parameter
 */
void model_lock_param(Model *mdl, const char *param_name) {
  ModelParam *param;

  param = model_get_param(mdl, param_name);

  if(param->locked) {
    g_warning("param '%s' is already locked", param_name);
  } else {
    param->locked = TRUE;
    mdl->n_unlocked--;
  }
}



/**
 * Writes a representation of the gradient to the provided filehandle.
 */
void model_write_grad(const Model *mdl, FILE *fh) {
  int i;
  char *sep;

  sep = "";

  fprintf(fh, "<");
  for(i = 0; i < mdl->n_param; i++) {
    if(!mdl->param[i].locked) {
      fprintf(fh, "%s%g",sep, mdl->param[i].deriv);
      sep = ", ";
    }
  }
  fprintf(fh, ">");
}


/**
 * Writes a representation of the gradient to the provided filehandle.
 */
void model_write_grad_ln(const Model *mdl, FILE *fh) {
  model_write_grad(mdl, fh);
  fprintf(fh, "\n");
}



/**
 * Writes string of paramters and values with trailing newline to
 * provided filehandle.
 */
void model_write_param_ln(const Model *mdl, FILE *fh) {

  model_write_param(mdl,fh);
  fprintf(fh, "\n");
}




/**
 * Returns a list of randomly selected unlocked ModelParam structures.
 * This can be used to lock a subset of parameters when attempting to
 * escape a ridge on the likelihood surface.
 */
GList *model_get_rand_unlocked_param(const Model *m, int n) {
  GList *all_unlocked;
  GList *selected, *elem;
  long i, n_remaining, idx;

  if(n > m->n_unlocked) {
    g_error("model_get_rand_unlocked_param: number of requested "
	    "parameters (%d) exceeds number of unlocked parameters(%d)",
	    m->n_unlocked, n);
  }

  if(m->n_unlocked == 0) {
    return NULL;
  }

  /* build list of unlocked parameters */
  all_unlocked = NULL;
  selected = NULL;
  n_remaining = 0;
  for(i = 0; i < m->n_param; i++) {
    if(!m->param[i].locked) {
      all_unlocked = g_list_append(all_unlocked, &m->param[i]);
      n_remaining++;
    }
  }

  if(n_remaining != m->n_unlocked) {
    g_error("model_get_rand_unlocked_param: expected %d unlocked param "
	    "but only found %ld", m->n_unlocked, n_remaining);
  }

  /* now remove n random parameters from the list of all unlocked param
   * and add them to the sampled list
   */
  for(i = 0; i < n; i++) {
    idx = rand() % n_remaining;
	//    idx = random() % n_remaining;
    
    elem = g_list_nth(all_unlocked,idx);

    all_unlocked = g_list_remove_link(all_unlocked, elem);
    selected = g_list_append(selected, elem->data);
    g_list_free(elem);
    n_remaining--;
  }
  
  return selected;
}
