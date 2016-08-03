
#include <glib.h>

#include "chr.h"



/**
 * Returns a (deep) copy of the provided chromosome
 */
Chromosome *chr_copy(const Chromosome *chr) {
  Chromosome *new_chr;

  new_chr = g_new(Chromosome, 1);
  new_chr->id = chr->id;

  if(chr->name) {
    new_chr->name = g_strdup(chr->name);
  } else {
    new_chr->name = NULL;
  }

  if(chr->assembly) {
    new_chr->assembly = g_strdup(chr->assembly);
  } else {
    new_chr->assembly = NULL;
  }

  new_chr->len = chr->len;
  
  return new_chr;
}


/**
 * Frees memory allocated for an array of chromosomes
 */
void chr_array_free(Chromosome *chrs, int n_chr) {
  int i;

  for(i = 0; i < n_chr; i++) {
    if(chrs[i].name) {
      g_free(chrs[i].name);
    }
    if(chrs[i].assembly) {
      g_free(chrs[i].assembly);
    }
  }

  g_free(chrs);
}


/**
 * Frees memory that was allocated for a single chromosome
 */
void chr_free(Chromosome *chr) {
  if(chr->name != NULL) {
    g_free(chr->name);
  }

  if(chr->assembly != NULL) {
    g_free(chr->assembly);
  }

  g_free(chr);
}

