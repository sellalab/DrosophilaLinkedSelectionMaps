#ifndef __GENE_H__
#define __GENE_H__

#include <math.h>
#include <glib.h>
#include "transcript.h"
#include "seqcoord.h"
#include "method.h"

#define GENE_LINE_PREFIX "GENE"


typedef struct {
  long id;             /* unique identifier */
  SeqCoord c;          /* gene coords */
  Method *method;      /* method used to predict gene, can be NULL */
  GSList *transcripts; /* transcript list, one for each splice form */
  char *name;
} Gene;


Gene *gene_group_transcripts(Transcript *trs, long num_trs, long *num_genes);

void gene_array_free(Gene *genes, long num_genes, int free_transcripts);

void gene_count_coding_substs(Gene *gene, Seq *q_seq, Seq *t_seq);

Transcript *gene_longest_transcript(Gene *gene);


inline void gene_write_flatfile_line(FILE *fh, Gene *gene);

void gene_write_flatfile(FILE *fh, Gene *gene_array, long num_genes);

Gene *gene_read_flatfile(char *filename, Transcript *trs, long num_trs,
			 long *num_read);


int gene_cmp(const void *p1,const void *p2);
int gene_cmp_nostrand(const void *p1,const void *p2);

#endif
