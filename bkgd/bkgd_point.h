
#ifndef __BKGD_POINT_H__
#define __BKGD_POINT_H__

#include <math.h>
typedef struct {
  long pos;       /* chr position in bp */
  double r_pos;   /* chr position in M */
  double b;       /* B value */
  double b_drv1;  /* first deriv of B */
  double b_drv2;  /* second deriv of B */
} BkgdPoint;


BkgdPoint *bkgd_point_dup(BkgdPoint *p);
void bkgd_point_copy(BkgdPoint *p_from, BkgdPoint *p_to);

#endif
