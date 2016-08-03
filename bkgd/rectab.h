#ifndef __REC_TAB_H__
#define __REC_TAB_H__

#include <math.h>
#include "seqfeat.h"

typedef struct {
  long start; /* start of block in bp */
  long end;   /* end of block in bp */
  double r_start;  /* start of block in morgans */
  double r_end;    /* end of block in morgans   */
  double rate;     /* recombination rate in block */
} RecRateBlock;


typedef struct {
  double chr_r_len;    /* total recombination distance */
  long n;            /* number of blocks of different rates */
  long cur;           /* index of current block */
  RecRateBlock *blk; /* array of blocks of different rates */
} RecRateTable;


RecRateTable *rectab_from_feats(SeqFeature *sf, long n_sf, long chr_len, 
				double scale);

double rectab_rpos(RecRateTable *rt, long pos);
double rectab_rate(RecRateTable *rt, long pos);
void rectab_free(RecRateTable *rt);

#endif
