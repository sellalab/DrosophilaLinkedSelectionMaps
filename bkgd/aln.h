#ifndef __ALN_H__
#define __ALN_H__

#include <math.h>
#include <glib.h>
#include <stdio.h>
#include "seq.h"
#include "method.h"

#define ALN_SEQNAME_SZ 51
#define ALN_AXT_LINE_BUF_SZ 1024
#define ALN_AXT_NUM_TOKS 9

#define ALN_QUERY  0
#define ALN_TARGET 1

typedef struct {
  long id;

  /* number of sequences in this alignment block */
  int n_seq; 

  /* array of SeqCoord structs n_seq long. 
   * The coords in each sequence spanned by this alignment block 
   */
  SeqCoord *coords;   

  /* Length of this alignment block. May not be same length as regions
   * spanned by coords, because of gaps and inserts.
   */
  long len;

  /* Score of this alignment block, if applicable */
  float score;

  /* A CIGAR string that represents the deletions, insertions and
   * matches in this alignment block, from the perspective of the
   * query sequence. The string contains a type character followed by
   * a length. Type characters can be M for match, I for insert or D
   * for deletion. If the size of the match/insert/deletion is 1 the
   * following number can be omitted.
   *
   * This is useful because it means that an alignment block can be
   * represented using coordinates and this relatively compressed
   * format string, rather than having all the symbols of the
   * alignment specified.
   *
   * This is the same as the CIGAR format used by the Ensembl database
   * and similar to the CIGAR format from exonerate.
   *
   * The following is an example alignment and the corresponding
   * CIGAR string:
   *
   * Alignment:
   *   ACTGA--CCCT-A
   *   ACTCAAACC-TTA
   *
   * CIGAR: 
   *   M5D2M2IMDM
   * 
   */
  char *cigar;

  /* NUC symbols (including gaps) as defined in nuc.h */
  unsigned char **sym;

  /* Lookup-table for genomic coords from alignment coords.
   * Coordinates are 1-based like those in SeqCoord objects, 
   * and are set to 0 when there is a gap (-) in the target sequence
   */
  long *aln2seq;

  /* The method that was used to generate this alignment 
   * can be set to NULL
   */
  Method *method;

} Alignment;




GSList *aln_read_axt(char *filename, GHashTable *trg_chrlen_tab);
void aln_write_axt(FILE *file, Alignment *aln, GHashTable *trg_chrlen_tab);
Seq *aln_read_axt_target_seq(char *filename, const long qseq_len);
Seq *aln_list_to_target_seq(GSList *alns, long qseq_len);

Seq **aln_read_tba_seqs(char *filename, long *n_seq);

Alignment *aln_read_maf(char *filename, char **seq_labels,
			int n_seq, long *num_alns);

void aln_free(Alignment *aln);
void aln_slist_free(GSList *alns);
void aln_array_free(Alignment *aln, long num_alns);

gchar *aln_seq2cigar(unsigned char *q_sym, unsigned char *h_sym, int aln_len);

#endif
