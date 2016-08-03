
#ifndef __FASTA_H__
#define __FASTA_H__

#include <math.h>
#include <glib.h>


#define FASTA_BUF_SZ 1024

typedef struct {
  char *path;
  char *header;
  char *seqstr;
  long seqlen;
} FASTA;


FASTA *fasta_read_file_array(char *filename, long *num_read);
FASTA *fasta_read_file(char *filename);
void fasta_free(FASTA *fasta);

#endif
