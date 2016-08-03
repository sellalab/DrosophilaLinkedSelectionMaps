#ifndef __TRANSCRIPT_H__
#define __TRANSCRIPT_H__

#include <math.h>
#include <stdio.h>
#include <glib.h>

#include "seq.h"
#include "seqcoord.h"
#include "method.h"

#define TRANSCRIPT_LINE_PREFIX "TRANSCRIPT"

typedef struct {
  long id;     /* unique identifier */
  SeqCoord c;  /* transcript coords */

  /* genomic coordinate representing CDS region */
  SeqCoord cds;
  int num_exons; 

  Method *method;

  /* array of exon coords, ordered from 5' to 3': on fwd strand
   * first exon has lowest coordinate, on rev strand last exon has
   * lowest coordinate.
   */
  SeqCoord *exons;

  char *name;
} Transcript;


void transcript_free(Transcript *tr);
void transcript_array_free(Transcript *trs, long num);

void transcript_parse_exons(Transcript *tr, char *ex_starts, char *ex_ends,
			    int zero_based);

Transcript *transcript_read_gene_file(char *table_file, long *num_read);
inline void transcript_write_bed(FILE *fh, Transcript *tr);
void transcript_write_flatfile(FILE *fh, Transcript *tr_array, long num_trs);
inline void transcript_write_flatfile_line(FILE *fh, Transcript *tr);
Transcript *transcript_read_flatfile(char *filename, long *num_read);

int transcript_has_cds(Transcript *tr);
SeqCoord *transcript_cds(Transcript *tr, int *n_coords);
long transcript_cds_len(Transcript *tr);
Seq *transcript_cds_seq(Transcript *tr, Seq *seq);
long *transcript_cds2genome_map(Transcript *tr);

SeqCoord *transcript_utr3(Transcript *tr, int *num_coords);
long transcript_utr3_len(Transcript *tr);
Seq *transcript_utr3_seq(Transcript *tr, Seq *seq);

SeqCoord *transcript_utr5(Transcript *tr, int *num_coords);
long transcript_utr5_len(Transcript *tr);
Seq *transcript_utr5_seq(Transcript *tr, Seq *seq);

Seq *transcript_cdna_seq(Transcript *tr, Seq *seq);

SeqCoord *transcript_upstream(Transcript *tr, Seq *seq, long num_bp);
Seq *transcript_upstream_seq(Transcript *tr, Seq *seq, long num_bp);

SeqCoord *transcript_tss(Transcript *tr, Seq *seq, 
			 long upstream_bp, long dnstream_bp);
Seq *transcript_tss_seq(Transcript *tr, Seq *seq, long upstream_bp,
			long dnstream_bp);

int transcript_cmp(const void *p1, const void *p2);
int transcript_cmp_end(const void *p1, const void *p2);
int transcript_cmp_nostrand(const void *p1, const void *p2);

SeqCoord *transcript_introns(Transcript *tr, int *num_introns);


#endif
