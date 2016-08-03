#ifndef __SEQ_H__
#define __SEQ_H__

#include <math.h>
#include <glib.h>


#include "seqcoord.h"
#include "fasta.h"

/* 
 * A DNA, RNA, or amino acid sequence 
 */
typedef struct {
  SeqCoord c;          /* coordinates of sequence */
  char *name;          /* name of the sequence    */
  unsigned char *sym;  /* array of nucleotide ids that make up the sequence */
  long len;   /* length of the sequence */
} Seq;


Seq *seq_read_fasta_file(char *filename);
Seq *seq_from_fasta(FASTA *fasta);

void seq_rev(Seq *fwd_seq);
void seq_comp(Seq *seq);
void seq_revcomp(Seq *fwd_seq);
void seq_nucs_revcomp(guchar *nuc_ids, long len);
Seq *seq_subseq(const Seq *seq, const SeqCoord *coord);

Seq *seq_subseq_coords(const Seq *seq, const SeqCoord *coords, 
		       const long num_coords);

char *seq_get_seqstr_buf(Seq *seq, char *buf);
gchar *seq_get_seqstr(Seq *seq);



Seq *seq_dup(Seq *seq);
void seq_free(Seq *seq);
void seq_array_free(Seq *seqs, long num);


#endif
