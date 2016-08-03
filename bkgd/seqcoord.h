#ifndef __SEQ_COORD_H__
#define __SEQ_COORD_H__

#include <math.h>
#include <glib.h>
#include <stdlib.h>
#include <stdio.h>

#include "chr.h"

#define STRAND_FWD 1
#define STRAND_NONE 0
#define STRAND_REV -1

#define strand_to_char(s) \
  ((s == STRAND_FWD) ? '+' : \
  ((s == STRAND_REV) ? '-' : \
                           '.'))

#define char_to_strand(c) \
  ((c == '+') ? STRAND_FWD : \
  ((c == '-') ? STRAND_REV : \
                STRAND_NONE  ))


#define seq_coord_len(c_ptr) ((c_ptr)->end - (c_ptr)->start + 1)


typedef struct {
  Chromosome *chr;
  long start;
  long end;
  short strand;
  char *seqname;
} SeqCoord;


void seq_coord_copy(const SeqCoord *src, SeqCoord *dst);

int seq_coord_cmp(const void *p1, const void *p2);
int seq_coord_cmp_end(const void *p1, const void *p2);
int seq_coord_cmp_nostrand(const void *p1, const void *p2);

long seq_coord_array_len(SeqCoord *c, long num_coords);

int seq_coord_ovlp(SeqCoord *sc1, SeqCoord *sc2, int cmp_strand);

void seq_coord_array_free(SeqCoord *scs, long num);

char *seq_coord_str(SeqCoord *sc);

#endif
