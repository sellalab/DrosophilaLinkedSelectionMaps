#ifndef __BKGD_DATA_H__
#define __BKGD_DATA_H__

#include <math.h>
#include "config.h"
#include "bkgd_evo_mdl.h"


typedef struct {
  int n_bin;
  BkgdBin *bin;
} BkgdEvoMdlData;


void bkgd_data_combine_bin(BkgdEvoMdlData *data, size_t n_bin);

void bkgd_data_combine_cat_bin(BkgdEvoMdlData *data, double min_cat,
			       double max_cat);

void bkgd_data_collapse_nex_bin(BkgdEvoMdlData *data);

double bkgd_data_mean_cat_val(BkgdEvoMdlData *data);

long bkgd_data_n_bin_sites(const BkgdBin *bin);

long bkgd_data_n_sites(BkgdEvoMdlData *data);

BkgdEvoMdlData *bkgd_data_read_data(Config *config);

BkgdEvoMdlData *bkgd_data_read_i_v_data(Config *config);

void bkgd_data_write_site_counts(FILE *fh, BkgdEvoMdlData *data);

void bkgd_data_jukes_cantor_correct_cats(BkgdEvoMdlData *data);

void bkgd_data_set_mdiv_cats(BkgdEvoMdlData *data);

#endif
