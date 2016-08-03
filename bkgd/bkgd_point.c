
#include <glib.h>
#include "bkgd_point.h"


/**
 * Creates a new copy of a point with the same attributes
 */
BkgdPoint *bkgd_point_dup(BkgdPoint *p) {
  BkgdPoint *new_p;

  new_p = g_new(BkgdPoint,1);
  new_p->pos = p->pos;
  new_p->r_pos = p->r_pos;
  new_p->b = p->b;
  new_p->b_drv1 = p->b_drv1;
  new_p->b_drv2 = p->b_drv2;

  return new_p;
}



/**
 * Copies attributes of first provided point into second provided
 * point.
 */
void bkgd_point_copy(BkgdPoint *p_from, BkgdPoint *p_to) {
  p_to->pos    = p_from->pos;
  p_to->r_pos  = p_from->r_pos;
  p_to->b      = p_from->b;
  p_to->b_drv1 = p_from->b_drv1;
  p_to->b_drv2 = p_from->b_drv2;
}
