
#ifndef __DIST_H__
#define __DIST_H__

#include <math.h>
#include <glib.h>
#include "seqmask.h"


#define DIST_NA G_MAXINT
#define DIST_NEAR 1
#define DIST_FAR  2


int *dist_calc_dists(SeqMask *mask, unsigned char mask_id);
int *dist_calc_dists_sign(SeqMask *mask, unsigned char mask_id, 
			  int use_neg_dists);

float *dist_calc_dists_float(SeqMask *mask, unsigned char mask_id);

float *dist_calc_genet_dists(SeqMask *mask, unsigned char mask_id,
			      float *recomb_rates);

float *dist_read_hk(char *hk_file, long seq_len);

float *dist_read_recomb_rate_feats(char *filename, long seq_len);

void dist_recomb_rates_to_dists(double *rates, long seq_len);

float *dist_read_recomb_rates(char *filename, long seq_len, double scale);

double *dist_recomb_rates_from_feats(SeqFeature *sf, long n_sf, long seq_len,
				     double scale, int interpolate);

float *dist_combine_hk(float *hk1, float u_scale1, 
			float *hk2, float u_scale2, 
			long seq_len, float max_hk);

#endif
