#ifndef __BKGD_FILES_H__
#define __BKGD_FILES_H__

#include <math.h>
// #include "config.h"
#include "seqfeat.h"


SeqFeature* load_conserved_from_file(char *filename, char *chr_name, long *pn_sf);
SeqFeature* load_genetic_map_from_file(char *filename, Chromosome *chr, long *pn_sf);
// GUY : long	get_chr_features(Chromosome* chr, char* chr_token, char* chr_features_file);
long	print_params(FILE *f_out, Chromosome *chr, char *params_file);

#endif
