#ifndef __INTERP_TAB_H__
#define __INTERP_TAB_H__

#include <math.h>
/* multiplier for fractional portion of decomposed double */
#define INTERP_TAB_FR_MULT 10.0
/* offset of fractional portion of decomposed double */
#define INTERP_TAB_FR_OFFSET 0.5
/* inverse of FR_OFFSET */
#define INTERP_TAB_FR_OFFSET_INV 2.0
/* minimum exponent stored in interpolation table */
#define INTERP_TAB_MIN_EXP -20
/* maximum exponent stored in interpolation table */
#define INTERP_TAB_MAX_EXP   5
/* number of entries in the interpolation table */
#define INTERP_TAB_N_ENTRIES ((INTERP_TAB_MAX_EXP-(INTERP_TAB_MIN_EXP))*INTERP_TAB_FR_MULT + 1)



typedef struct {
  double *x;  /* distance */
  double *y;  /* len */
  double **z; /* blck sum */
  double **d_inv; /* 1/((x2-x1)(y2-y1), distance between points inverted */
} InterpTab;


double interp_tab_lookup(InterpTab *tab, double x, double y);
double interp_tab_lookup_1d(InterpTab *tab, double x);

InterpTab *interp_tab_create(double (*f)(double,double,void *), void *);
InterpTab *interp_tab_create_1d(double (*f)(double,void *), void *);

void interp_tab_free(InterpTab *tab);


#endif
