#include <stdio.h>
#include <stdlib.h>
#include <glib.h>
#include <string.h>
#include <ctype.h>

#include "nuc.h"
#include "seq.h"
#include "seqcoord.h"
#include "fasta.h"

#include "util.h"


/*
 * Converts a FASTA data structure into a Seq data structure.  The
 * FASTA data structure cannot be used after the conversion as the
 * same memory is reused for the sequence (and the pointer to the
 * sequence in the FASTA data structure is set to NULL).  This is done
 * to save memory because sequences can be very large.
 */
Seq *seq_from_fasta(FASTA *fasta) {
  gulong i;
  Seq *seq;

  seq = g_new(Seq, 1);

  seq->len = fasta->seqlen;
  seq->sym = (unsigned char *)fasta->seqstr;

  seq->c.start = 1;
  seq->c.end = seq->len;
  seq->c.strand = STRAND_FWD;
  seq->c.chr = NULL;

  /* convert characters in sequence string to nucleotide symbols */
  for(i = 0; i < seq->len; i++) {
    seq->sym[i] = nuc_char_to_id(seq->sym[i]);
  }

  if(fasta->header == NULL) {
    seq->c.seqname = NULL;
  } else {
    seq->name = g_strdup(fasta->header);

    /* strip leading '>' */
    if(seq->name[0] == '>') {
      i = 1;
      while(seq->name[i] != '\0') {
	seq->name[i-1] = seq->name[i];
	i++;
      }
      seq->name[i-1] = '\0';
    }
  }

  seq->c.seqname = g_strdup(seq->name);
  fasta->seqstr = NULL;
  
  return seq;
}



/*
 * Convenience function. Allows a sequence to be created using a FASTA
 * filename without having to create a FASTA object first.
 */
Seq *seq_read_fasta_file(gchar *filename) {
  FASTA *fasta;
  Seq *seq;

  fasta = fasta_read_file(filename);
  seq = seq_from_fasta(fasta);
  fasta_free(fasta);
  
  return seq;
}



/*
 * Returns a string representation of the nucleotides of this
 * sequence. The returned string should be freed when it is no longer
 * needed.
 */
gchar *seq_get_seqstr(Seq *seq) {
  gulong i;

  gchar *seqstr;
  seqstr = g_new(gchar, seq->len+1);

  for(i = 0; i < seq->len; i++) {
    seqstr[i] = nuc_id_to_char(seq->sym[i]);
  }
  seqstr[seq->len] = '\0';

  return seqstr;
}


/*
 * Fills the provided buffer with a null-terminated string
 * representation of the nucleotides in this sequence. The buffer must
 * be of length seq->len + 1 or greater.
 */
gchar *seq_get_seqstr_buf(Seq *seq, gchar *buf) {
  gulong i;

  for(i = 0; i < seq->len; i++) {
    buf[i] = nuc_id_to_char(seq->sym[i]);
  }
  buf[seq->len] = '\0';

  return buf;
}




/*
 * Returns a new Seq object that represents the subsequence defined by
 * the provided coordinates. If the coordinates are on the opposite
 * strand of the original sequence, the returned sequence is the
 * reverse complement of the region (i.e. the sequence is concatenated
 * in order of the provided coordinates and then
 * reverse-complemented).
 *
 * The provided coordinates should be entirely within the region
 * specified by the sequence, otherwise an error is raised.
 *
 * The returned sequence will have coords 1-len where len is the
 * length of the returned sequence. This is because there is no simple
 * way to represent the genomic coordinates of a sequence that may be
 * composed of several disjoint regions. If the genomic coordinates
 * are required the seq_subseq function, which takes only a single
 * coordinate, should be used instead.
 *
 * The returned Seq should be freed when it is no longer needed by
 * using the seq_free function.
 */
Seq *seq_subseq_coords(const Seq *seq, const SeqCoord *coords, 
		       const long num_coords) {
  long len;
  Seq *new_seq;
  long i, j, array_start;
  short int strand;

  if(num_coords < 1) {
    g_error("seq_subseq_coords: must provide at least 1 coordinate");
  }

  strand = coords[0].strand;

  /* figure out length of subseq and check validity of coords */
  len = 0;
  for(i = 0; i < num_coords; i++) {
    if(coords[i].start > coords[i].end || 
       coords[i].start < seq->c.start || coords[i].end > seq->c.end) {
      g_error("seq_subseq_coords: Request for bad coordinates %lu-%lu "
	      "from seq with coordinates %lu-%lu", 
	      coords[i].start, coords[i].end, seq->c.start, seq->c.end);
    }
    if(coords[i].strand != strand) {
      g_error("seq_subseq_coords: retrieval of subseq from multiple "
	      "strands is not implemented.");
    }

    len += coords[i].end - coords[i].start + 1;
  }
  
  new_seq = g_new(Seq, 1);
  new_seq->c.seqname = NULL;
  new_seq->c.start = 1;
  new_seq->c.end = len;
  new_seq->c.strand = strand;
  new_seq->len = len;
  new_seq->sym = g_new(unsigned char, len);
  new_seq->name = NULL;



  /* populate new seq from old seq */
  j = 0;
  for(i = 0; i < num_coords; i++) {
    len = coords[i].end - coords[i].start + 1;

    /* figure out start position in array */
    if(seq->c.strand == STRAND_REV) {
      /* Sequence is reversed, so coords are relative to end.
       * Array starts at 0, so no need for +1.
       */
      array_start = seq->c.end - coords[i].end;
    } else {
      /* Sequence is not reversed, coords relative to start.
       * Array starts at 0, so no need for +1.
       */
      array_start = coords[i].start - seq->c.start;
    }

    /* copy region of sequence */
    memcpy(&new_seq->sym[j], &seq->sym[array_start], len);

    if((seq->c.strand == STRAND_FWD && strand == STRAND_REV) ||
       (seq->c.strand == STRAND_REV && strand == STRAND_FWD)) {
      /* reverse complement each copied portion if sequences are on
       * opposite strands
       */
      seq_nucs_revcomp(&new_seq->sym[j], len);
    }

    j += len;
  }
  
  return new_seq;
}


/*
 * Returns an new Seq object that represents the subsequence defined
 * by the provided SeqCoord struct. If the coordinate is on the
 * reverse strand the returned sequence is the reverse complement of
 * the region. 
 *
 * The coordinates of the returned sequence will be the same as those
 * requested (unlike the seq_subseq_coords function).
 *
 */
Seq *seq_subseq(const Seq *seq, const SeqCoord *coord) {
  Seq *new_seq;

  new_seq = seq_subseq_coords(seq, coord, 1);

  /* set coordinates */
  new_seq->c.start  = coord->start;
  new_seq->c.end    = coord->end;
  new_seq->c.strand = coord->strand;
  new_seq->c.seqname = NULL;
  new_seq->c.chr = coord->chr;

  return new_seq;
}


/*
 * Does an in-place complement (but not reverse) of the provided sequence.
 * Does not alter the coordinates of the sequence.
 */
void seq_comp(Seq *seq) {
  long i;
  for(i = 0; i < seq->len; i++) {
    seq->sym[i] = nuc_comp(seq->sym[i]);
  }
}


/**
 * Performs an in-place of reversal (but not complement!) of the
 * provided sequence. Does not alter the coordinates of the sequence.
 */
void seq_rev(Seq *fwd) {
  util_breverse(fwd->sym, fwd->len);
}


/**
 * Does an in-place reverse complement of a provided string of
 * nucleotide ids.
 */
void seq_nucs_revcomp(unsigned char *nuc_ids, long len) {
  long i;

  util_breverse(nuc_ids, len);
  for(i = 0; i < len; i++) {
    nuc_ids[i] = nuc_comp(nuc_ids[i]);
  }
}


/**
 * Performs an in-place reverse-compliment of the provided sequence.
 * and flips the strand of the sequence coordinate.
 */
void seq_revcomp(Seq *fwd_seq) {
  seq_rev(fwd_seq);
  seq_comp(fwd_seq);
  
  if(fwd_seq->c.strand == STRAND_FWD) {
    fwd_seq->c.strand = STRAND_REV;
  }
  else if(fwd_seq->c.strand == STRAND_REV) {
    fwd_seq->c.strand = STRAND_FWD;
  }
}


/**
 * Creates a new copy of a sequence and returns it. The sequence
 * should be freed when it is no longer needed.
 */
Seq *seq_dup(Seq *seq) {
  Seq *new_seq;

  new_seq = g_new(Seq, 1);
  new_seq->len = seq->len;

  if(seq->name == NULL) {
    new_seq->name = NULL;
  } else {
    new_seq->name = g_strdup(seq->name);
  }

  new_seq->c.start = seq->c.start;
  new_seq->c.end   = seq->c.end;
  new_seq->c.strand = seq->c.strand;
  new_seq->c.chr = seq->c.chr;

  if(seq->c.seqname == NULL) {
    new_seq->c.seqname = NULL;
  } else {
    new_seq->c.seqname = g_strdup(seq->c.seqname);
  }
  new_seq->sym = g_new(unsigned char, seq->len);
  new_seq->sym = memcpy(new_seq->sym, seq->sym, seq->len);

  return new_seq;
}


/**
 * Frees an allocated sequence structure. 
 */  
void seq_free(Seq *seq) {
  g_free(seq->sym);

  if(seq->name != NULL) {
    g_free(seq->name);
  }

  if(seq->c.seqname != NULL) {
    g_free(seq->c.seqname);
  }

  g_free(seq);
}


/**
 * Frees an array of sequence structures. 
 */  
void seq_array_free(Seq *seqs, long num) {
  gulong i;

  for(i = 0; i < num; i++) {
    g_free(seqs[i].sym);
    if(seqs[i].name != NULL) {
      g_free(seqs[i].name);
    }
  }
  g_free(seqs);
}

