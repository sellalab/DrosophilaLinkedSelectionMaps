#ifndef __CHR_H__
#define __CHR_H__

#include <math.h>
#include <glib.h>


typedef struct {
  int id;
  char *name;
  long len;
  char *assembly;
} Chromosome;


void chr_array_free(Chromosome *chr, int n_chr);
Chromosome *chr_copy(const Chromosome *chr);
void chr_free(Chromosome *chr);


#endif
