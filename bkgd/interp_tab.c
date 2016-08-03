#include <glib.h>
#include <stdio.h>
#include <math.h>

#include "interp_tab.h"


/**
 * Converts the provided interpolation table index into a double
 */
static double idx2dbl(int idx) {
  int fr = (idx-1) % (int)INTERP_TAB_FR_MULT;
  int x  = (idx-1) / INTERP_TAB_FR_MULT;

  if(idx < 0) {
    g_error("idx2dbl: provided idx(%d) must be positive", idx);
  }

  if(idx == 0) {
    return 0.0;
  }

  return (((INTERP_TAB_FR_OFFSET*fr)/INTERP_TAB_FR_MULT)+INTERP_TAB_FR_OFFSET) * pow(2.0, x+INTERP_TAB_MIN_EXP);
}


/**
 * Converts the provided double to an interpolation table index.
 */
static inline int dbl2idx(double d) {
  double fr;
  int x,idx;

  if(d < 0.0) {
    g_error("dbl2idx: provided double (%g) must be positive", d);
  }

  if(d == 0.0) {
    return 0;
  }

  fr = frexp(d, &x);

  if(x < INTERP_TAB_MIN_EXP) {
    return 0;
  }

  if(x > INTERP_TAB_MAX_EXP) {
    g_warning("dbl2idx: value %g exceeds interpolation table maximum\n",d);
    x = INTERP_TAB_MAX_EXP;
  }
    
  idx = (((fr-INTERP_TAB_FR_OFFSET)*INTERP_TAB_FR_OFFSET_INV) + (x-INTERP_TAB_MIN_EXP))*INTERP_TAB_FR_MULT + 1;
  return idx;
}



/**
 * Looks up entries in table that flank requested data point and
 * performs interpolation to infer value at requested point.
 * This is currently implemented as linear interpolation, but
 * something more sophisticated could be introduced.
 * 
 * This function is for interpolation of a function with 1 variable,
 * for 2d interpolation use interp_tab_lookup().
 */
double interp_tab_lookup_1d(InterpTab *tab, double x) {
  double x1, x2, y1, y2, y;
  int x_idx1, x_idx2;

  if(tab->z != NULL) {
    g_error("interp_tab_lookup_1d: table initialized for 2d interpolation. "
	    "interp_tab_lookup function should be used instead.");
  }

  /* get indexes that correspond to table entries on either side
   * of requested value
   */
  x_idx2 = dbl2idx(x);
  if(x_idx2 == 0) {
    x_idx1 = 0;
    x_idx2 = 1;
  } else {
    if(tab->x[x_idx2] >= x || x_idx2 == INTERP_TAB_N_ENTRIES-1) {
      x_idx1 = x_idx2 - 1;
    } else {
      x_idx1 = x_idx2;
      x_idx2 = x_idx1+1;
    }
  }

  x1 = tab->x[x_idx1];
  x2 = tab->x[x_idx2];
  y1 = tab->y[x_idx1];
  y2 = tab->y[x_idx2];

  /* perform linear interpolation */
  y = y1 + ((y2-y1)/(x2-x1)) * (x-x1);

  /*
   * fprintf(stderr, "y1=%g, y2=%g, x1=%g, x2=%g, x=%g, y=%g\n", 
   * y1,y2,x1,x2,x,y);
   */

  return y;
}



/**
 * Looks up nearest entries in table to the requested x and y values,
 * and performs bilinear interpolation to estimate z value at
 * requested point.
 */
/*inline */double interp_tab_lookup(InterpTab *tab, double x, double y) {
  double x1, x2, y1, y2, z;
  int x_idx1, x_idx2, y_idx1, y_idx2;

  if(tab->z == NULL) {
    g_error("interp_tab_lookup: table initialized for 1d interpolation. "
	    "interp_tab_lookup_1d function should be used instead.");
  }

  /* get indexes that correspond to table entries on either side
   * of requested value
   */
  x_idx2 = dbl2idx(x);
  if(x_idx2 == 0) {
    x_idx1 = 0;
    x_idx2 = 1;
  } else {
    if(tab->x[x_idx2] >= x || x_idx2 == INTERP_TAB_N_ENTRIES-1) {
      x_idx1 = x_idx2 - 1;
    } else {
      x_idx1 = x_idx2;
      x_idx2 = x_idx1+1;
    }
  }

  y_idx2 = dbl2idx(y);
  if(y_idx2 == 0) {
    y_idx1 = 0;
    y_idx2 = 1;
  } else {
    if(tab->y[x_idx2] >= y || y_idx2 == INTERP_TAB_N_ENTRIES-1) {
      y_idx1 = y_idx2 - 1;
    } else {
      y_idx1 = y_idx2;
      y_idx2 = y_idx1+1;
    }
  }

  /* perform bilinear interpolation */
  x1 = tab->x[x_idx1];
  x2 = tab->x[x_idx2];
  y1 = tab->x[y_idx1];
  y2 = tab->x[y_idx2];

  /* d_inv is 1/ ((x2-x1)*(y2-y1)); */

  z = (tab->z[x_idx1][y_idx1] * (x2-x)*(y2-y) + 
       tab->z[x_idx2][y_idx1] * (x-x1)*(y2-y) +
       tab->z[x_idx1][y_idx2] * (x2-x)*(y-y1) +
       tab->z[x_idx2][y_idx2] * (x-x1)*(y-y1)) * tab->d_inv[x_idx1][y_idx1];

  return z;
}



/**
 * Creates an interpolation table for interpolation over 1 dimension
 */
InterpTab *interp_tab_create_1d(double (*f)(double, void *), void *param) {
  InterpTab *tab;
  int idx_x;
  double x;

  tab = g_new(InterpTab, 1);

  tab->x  = g_new(double, INTERP_TAB_N_ENTRIES);
  tab->y  = g_new(double, INTERP_TAB_N_ENTRIES);
  tab->z = NULL;
  tab->d_inv = NULL;
  
  for(idx_x = 0; idx_x < INTERP_TAB_N_ENTRIES; idx_x++) {
    x = idx2dbl(idx_x);

    tab->x[idx_x] = x;
    tab->y[idx_x] = f(x, param);

    /*
     * fprintf(stderr, "interp_tab_create_1d: x=%g, y=%g\n", 
     * tab->x[idx_x], tab->y[idx_x]);
     */
  }
  
  return tab;
}



/**
 * Creates an interpolation table for the provided (2d) function
 */
InterpTab *interp_tab_create(double (*f)(double,double, void *), void *param) {
  InterpTab *tab;
  int idx_x, idx_y;
  double x,y,z, x_diff, y_diff;

  tab = g_new(InterpTab, 1);

  tab->x  = g_new(double, INTERP_TAB_N_ENTRIES);
  tab->y  = g_new(double, INTERP_TAB_N_ENTRIES);
  tab->z = g_new(double *, INTERP_TAB_N_ENTRIES);

  tab->d_inv = g_new(double *, INTERP_TAB_N_ENTRIES);
  
  for(idx_x = 0; idx_x < INTERP_TAB_N_ENTRIES; idx_x++) {
    x = idx2dbl(idx_x);

    tab->x[idx_x] = x;

    tab->z[idx_x] = g_new(double, INTERP_TAB_N_ENTRIES);

    for(idx_y = 0; idx_y < INTERP_TAB_N_ENTRIES; idx_y++) {
      y = idx2dbl(idx_y);

      tab->y[idx_y] = y;

      /* calculate function value at this point */
      z = f(x,y, param);

      tab->z[idx_x][idx_y] = z;
    }
  }

  tab->d_inv = g_new(double *, INTERP_TAB_N_ENTRIES);
  for(idx_x = 0; idx_x < INTERP_TAB_N_ENTRIES; idx_x++) {
    tab->d_inv[idx_x] = g_new(double, INTERP_TAB_N_ENTRIES);

    if(idx_x < INTERP_TAB_N_ENTRIES-1) {
      x_diff = tab->x[idx_x+1] - tab->x[idx_x];
    } else {
      /* this becomes extrapolation past end of table */
      x_diff = tab->x[idx_x] - tab->x[idx_x-1];
    }

    for(idx_y = 0; idx_y < INTERP_TAB_N_ENTRIES; idx_y++) {
      if(idx_y < INTERP_TAB_N_ENTRIES-1) {
	y_diff = tab->y[idx_y+1] - tab->y[idx_y];
      } else {
	/* extrapolation rather than interpolation */
	y_diff = tab->y[idx_y] - tab->y[idx_y-1];
      }
      tab->d_inv[idx_x][idx_y] = 1.0 / (x_diff*y_diff);
    }
  }
  
  return tab;
}



/**
 * Frees the memory allocated for the provided interpolation table.
 */
void interp_tab_free(InterpTab *tab) {
  int i;

  g_free(tab->x);
  g_free(tab->y);

  if(tab->z) {
    for(i = 0; i < INTERP_TAB_N_ENTRIES; i++) {
      g_free(tab->z[i]);
    }
    g_free(tab->z);
  }
  if(tab->d_inv) {
    for(i = 0; i < INTERP_TAB_N_ENTRIES; i++) {
      g_free(tab->d_inv[i]);
    }
    g_free(tab->d_inv);
  }
  
  g_free(tab);
}
