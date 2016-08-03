#include <glib.h>
#include <ctype.h>
#include <string.h>

#include "aa.h"
#include "aln.h"
#include "transcript.h"
#include "seqmask.h"
#include "nuc.h"
#include "dist.h"

#include "util.h"

/*
 * Returns a new sequence mask for the region represented by the
 * provided SeqCoord. A shallow copy of the provided coord is made (
 * just the chr ptr is copied, not the chr itself).
 */
SeqMask *seqmask_new(SeqCoord *coord) {
  SeqMask *mask;
  unsigned char i;
  long j;

  if(coord->start > coord->end) {
    g_error("seqmask_new: coord start must be <= end");
  }

  mask = g_new(SeqMask, 1);

  mask->len = coord->end - coord->start + 1;
  mask->mask = g_new(unsigned char, mask->len);

  /* initialize coordinate */
  mask->c.start   = coord->start;
  mask->c.end     = coord->end;
  mask->c.strand  = coord->strand;

  if(coord->seqname == NULL) {
    mask->c.seqname = NULL;
  } else {
    mask->c.seqname = g_strdup(coord->seqname);
  }

  mask->c.chr = coord->chr;

  /* clear sequence mask */
  for(j = 0; j < mask->len; j++) {
    mask->mask[j] = SEQMASK_CLEAR;
  }

  /* initialize mask ids */
  for(i = 0; i < SEQMASK_NUM_IDS; i++) {
    mask->used_ids[i] = FALSE;
  }
  
  return mask;
}



/*
 * Shifts the sequence mask so that it is associated with a new
 * region. The mask is cleared, as are the mask identifiers, but
 * memory is only re-allocated if the size of the region has changed.
 */
void seqmask_move(SeqMask *mask, SeqCoord *coord) {  
  long  i;
  long new_len;

  if(coord->start > coord->end) {
    g_error("seqmask_new: coord start must be <= end");
  }

  new_len = coord->end - coord->start + 1;
  
  if(new_len != mask->len) {
    mask->mask = g_realloc(mask->mask, sizeof(unsigned char) * new_len);
    mask->len = new_len;
  }

  mask->c.start   = coord->start;
  mask->c.end     = coord->end;
  mask->c.strand  = coord->strand;


  if(coord->seqname == NULL) {
    mask->c.seqname = NULL;
  } else {
    if(mask->c.seqname == NULL) {
      mask->c.seqname = g_strdup(coord->seqname);
    } else {
      /* re-copy sequence name, if it has changed */
      if(strcmp(mask->c.seqname, coord->seqname) != 0) {
	g_free(mask->c.seqname);
	mask->c.seqname = g_strdup(coord->seqname);
      }
    }
  }

  mask->c.chr = coord->chr;

  /* clear sequence mask */
  for(i = 0; i < mask->len; i++) {
    mask->mask[i] = SEQMASK_CLEAR;
  }

  /* initialize mask ids */
  for(i = 0; i < SEQMASK_NUM_IDS; i++) {
    mask->used_ids[i] = FALSE;
  }
}



/**
 * Frees the memory allocated for a sequence mask.
 */
void seqmask_free(SeqMask *mask) {
  g_free(mask->mask);

  if(mask->c.seqname != NULL) {
    g_free(mask->c.seqname);
  }

  g_free(mask);
}



/*
 * Returns the value of a free identifier in the provided mask that
 * can be used for masking and sets it as used.  If this mask has no
 * more free identifiers an error occurs and the program is exited.
 */
unsigned char seqmask_get_id(SeqMask *mask) {
  unsigned char i;

  for(i = 0; i < SEQMASK_NUM_IDS; i++) {
    if(mask->used_ids[i] == FALSE) {
      mask->used_ids[i] = TRUE;
      return (0x01 << i);
    }
  }

  g_error("seqmask_get_id(): All seqmask ids are already taken.");

  return 0;
}



/*
 * Returns a new seqmask which consists of the regions of the original
 * seqmask specified by the provided coordinate. The flagged regions
 * and ids from the original mask are copied as well. If the requested
 * coordinate is on the opposite strand of the original mask, the
 * submask is reversed (for consistancy with sequences from the
 * opposite strand).
 *
 * The new mask should be freed when it is no longer needed.
 */
SeqMask *seqmask_submask(SeqMask *mask, SeqCoord *coord) {
  SeqMask *new_mask;
  long offset;
  long i;
  unsigned char j;

  if(coord->start < mask->c.start || coord->end > mask->c.end) {
    g_error("seqmask_submask: requested region (%lu-%lu) must be entirely "
	    "contained by seqmask (%lu-%lu)\n", coord->start, coord->end,
	    mask->c.start, mask->c.end);
  }

  if(coord->start > coord->end) {
    g_error("seqmask_submask: invalid coordinates");
  }

  new_mask = seqmask_new(coord);

  offset = coord->start - mask->c.start;

  /* copy old sequence mask */
  for(i = 0; i < new_mask->len; i++) {
    new_mask->mask[i] = mask->mask[i + offset];
  }

  /* copy old mask ids */
  for(j = 0; j < SEQMASK_NUM_IDS; j++) {
    new_mask->used_ids[j] = mask->used_ids[j];
  }
  
  if((mask->c.strand == STRAND_FWD && new_mask->c.strand == STRAND_REV) ||
     (mask->c.strand == STRAND_REV && new_mask->c.strand == STRAND_FWD)) {
    util_breverse(new_mask->mask, new_mask->len);    
  }

  return new_mask;
}


/*
 * "Inverts" a sequence mask, by flagging bases using new_id, which
 * do *not* have the old_id bit-flag set.
 */
void seqmask_invert_mask(SeqMask *mask, unsigned char old_id, 
			 unsigned char new_id) {
  long i;

  for(i = 0; i < mask->len; i++) {
    if((mask->mask[i] & old_id) == 0) {
      mask->mask[i] = mask->mask[i] | new_id;
    }
  }
}



/*
 * Given a mask and an id, sets the bits corresponding to the given id
 * to 1 in positions of the mask's array that overlap with softmasked
 * repeats.  Softmasked repeats are identified as lower-case
 * nucleotides in the provided sequence.
 */
void seqmask_soft_repeats(SeqMask *mask, unsigned char mask_id, 
			  char *seq_str) {
  long i;
  
  for(i = 0; i < mask->len; i++) {
    if(islower(seq_str[i])) {
      mask->mask[i] = mask->mask[i] | mask_id;
    }
  }
}



/*
 * Given a mask and an id, sets the bits corresponding to the given id
 * to 1 in positions of the mask's array that DON'T overlap with
 * softmasked repeats.  I.e. all uppercase bases will be
 * masked.
 */
void seqmask_non_soft_repeats(SeqMask *mask, unsigned char mask_id, 
			      char *seq_str) {			      
  long i;
  
  for(i = 0; i < mask->len; i++) {
    if(isupper(seq_str[i])) {
      mask->mask[i] = mask->mask[i] | mask_id;
    }
  }
}




/**
 * A Function to mask bases that are flanked by more or less than a specified
 * number of flanking G+C bases. Bases that have any flanking Ns, are
 * also omitted.
 */
void seqmask_gc_flanking(SeqMask *mask, unsigned char mask_id, Seq *seq,
			 int num_flanking, int min_gc, int max_gc) {
  
  long i, pos;
  int gc_count, n_count, win_sz;
  unsigned char base;
  

  if(seq->len != mask->len) {
    g_error("seqmask_gc_flanking: seq len (%lu) should equal mask len (%lu)",
	    seq->len, mask->len);
  }

  win_sz = num_flanking*2 + 1;

  if((seq->len < win_sz) || (num_flanking < 1)) {
    return;
  }

  /* get counts for first window (exclude base in middle of window) */
  gc_count = 0;
  n_count = 0;
  pos = num_flanking;

  /* upstream flanking */
  for(i = 0; i < num_flanking; i++) {
    if(seq->sym[i] == NUC_G || seq->sym[i] == NUC_C) {
      gc_count++;
    } 
    else if(seq->sym[i] == NUC_N) {
      n_count++;
    }
  }
  /* downstream flanking */
  for(i = num_flanking+1; i < win_sz; i++) {
    if(seq->sym[i] == NUC_G || seq->sym[i] == NUC_C) {
      gc_count++;
    } 
    else if(seq->sym[i] == NUC_N) {
      n_count++;
    }
  }

  if((n_count > 0) || (gc_count > max_gc) || (gc_count < min_gc)) {
    mask->mask[pos] = mask->mask[pos] | mask_id;
  }
  
  /* now move along rest of sequence, updating counts of flanking
   * G+Cs and Ns as we go
   */
  for(pos = num_flanking+1; pos < seq->len - num_flanking; pos++) {
    /* peek at upstream flanking base that is exiting window */
    base = seq->sym[pos - num_flanking - 1];
    if(base == NUC_G || base == NUC_C) {
      gc_count--;
    }
    else if(base == NUC_N) {
      n_count--;
    }

    /* peek at base that was center base last iter (now upstream flanking) */
    base = seq->sym[pos-1];
    if(base == NUC_G || base == NUC_C) {
      gc_count++;
    }
    else if(base == NUC_N) {
      n_count++;
    }

    /* current center base was a downstream flanking base last iter */
    base = seq->sym[pos];
    if(base == NUC_G || base == NUC_C) {
      gc_count--;
    }
    if(base == NUC_N) {
      n_count--;
    }


    /* peek at downstream flanking base that is entering window */
    base = seq->sym[pos + num_flanking];
    if(base == NUC_G || base == NUC_C) {
      gc_count++;
    }
    if(base == NUC_N) {
      n_count++;
    }
    
    /* mask center base if there are too many or too few flanking G+Cs
     * or if there are any flanking Ns
     */
    if((n_count > 0) || (gc_count > max_gc) || (gc_count < min_gc)) {
      mask->mask[pos] = mask->mask[pos] | mask_id;
    }
  }

  return;
}

  



/**
 * Masks regions where a CpG occurs in *any* of the provided sequences.
 */
void seqmask_cpg(SeqMask *mask, unsigned char mask_id, Seq **seqs, 
		 int num_seqs) {
  long i, seq_len;
  int j;
  int is_cpg;

  if(num_seqs < 1) {
    return;
  }

  seq_len = seqs[0]->len;

  if(seq_len != mask->len) {
    g_error("seqmask_cpg: seq len (%lu) should equal mask len (%lu)",
	    seq_len, mask->len);
  }

  for(i = 1; i < seq_len; i++) {
    is_cpg = FALSE;
    
    /* look for presence of CpG in sequences */
    for(j = 0; j < num_seqs; j++) {
      if(seqs[j]->sym[i-1] == NUC_C && seqs[j]->sym[i] == NUC_G) {
	is_cpg = TRUE;
	break;
      }
    }

    if(is_cpg) {
      /* Mask both bases that make up CpG */
      mask->mask[i-1] = mask->mask[i-1] | mask_id;
      mask->mask[i] = mask->mask[i] | mask_id;
    }
  }
}





/**
 * Masks all bases that are in a CpG context. See
 * seqmask_cpg_context_region for details.
 */
void seqmask_cpg_context(SeqMask *mask, unsigned char mask_id, Seq *seq) {
  long i;

  if(seq->len != mask->len) {
    g_error("seqmask_cpg_context: seq len (%lu) should equal mask len (%lu)",
	    seq->len, mask->len);
  }

  for(i = 1; i < mask->len; i++) {
    if(seq->sym[i-1] == NUC_C) {
      mask->mask[i] = mask->mask[i] | mask_id;
    } 
    else if((i < seq->len - 1) && (seq->sym[i+1] == NUC_G)) {
      mask->mask[i] = mask->mask[i] | mask_id;
    } 
  }
}




/*
 * Given a set of aligned sequences (all of the same length), masks
 * dinucleotides that may have undergone methylated cytosine CpG
 * deamination.  These are dinucleotides which are CpG, TpG, or CpA in
 * *all* of the species in the alignment. Additionally, if all species
 * have the same dinucleotide it is not considered to be a CpG
 * deamination.  For example, CG/CG/TG would be masked, but CG/CG/CG
 * would not be.
 */
void seqmask_deam(SeqMask *mask, unsigned char mask_id, 
		  Seq **seqs, int num_seqs) {
  long i, seq_len;
  int is_cpg, all_same;
  int j;

  if(num_seqs < 1) {
    return;
  }

  seq_len = seqs[0]->len;

  if(seq_len != mask->len) {
    g_error("seqmask_deam: seq len (%lu) should equal mask len (%lu)",
	    seq_len, mask->len);
  }

  for(i = 1; i < seq_len; i++) {
    is_cpg = TRUE;
    all_same = TRUE;
    
    /* look for presence of CpG/TpG/CpA in sequences */
    for(j = 0; j < num_seqs; j++) {
      if((seqs[j]->sym[i-1] == NUC_C && seqs[j]->sym[i] == NUC_G) ||
	 (seqs[j]->sym[i-1] == NUC_T && seqs[j]->sym[i] == NUC_G) ||
	 (seqs[j]->sym[i-1] == NUC_C && seqs[j]->sym[i] == NUC_A)) {
	/* check if same dinucleotide as previous seq in alignment */
	if((j > 0) && ((seqs[j]->sym[i-1] != seqs[j-1]->sym[i-1]) ||
		       (seqs[j]->sym[i] != seqs[j-1]->sym[i]))) {
	  all_same = FALSE;
	}
      } else {
	is_cpg = FALSE;
	break;
      }
    }

    if(is_cpg && !all_same) {
      /* mask both bases that are CpG/TpG/CpA */
      mask->mask[i-1] = mask->mask[i-1] | mask_id;
      mask->mask[i] = mask->mask[i] | mask_id;
    }
  }
}



/*
 * Masks all bases that are not part of a CpG dinucleotide in at least
 * one of the provided sequences.
 */
void seqmask_non_cpg(SeqMask *mask, unsigned char mask_id, Seq **seqs, 
		     int num_seqs) {

  long i, seq_len;
  int j;
  int is_cpg;

  seq_len = seqs[0]->len;

  if(seq_len != mask->len) {
    g_error("seqmask_non_cpg: seq len (%lu) should equal mask len (%lu)",
	    seq_len, mask->len);
  }
  
  for(i = 0; i < seq_len; i++) {
    is_cpg = FALSE;
    
    /* look for presence of CpG in sequences */
    for(j = 0; j < num_seqs; j++) {
      if(i < seqs[0]->len-1 && 
	 (seqs[j]->sym[i] == NUC_C && seqs[j]->sym[i+1] == NUC_G)) {
	/* this nuc is C and next is G */
	is_cpg = TRUE;
	break;
      }
      if(i > 0 && 
	 (seqs[j]->sym[i-1] == NUC_C && seqs[j]->sym[i] == NUC_G)) {
	/* this nuc is G and prev was C */
	is_cpg = TRUE;
	break;
      }
    }

    if(!is_cpg) {
      mask->mask[i] = mask->mask[i] | mask_id;
    }
  }

  return;
}


/*
 * Masks bases which are Ns in any of the provided sequences
 */
void seqmask_ns(SeqMask *mask, unsigned char mask_id, Seq **seqs,
		int n_seqs) {
  long i;
  int j;

  if(seqs[0]->len != mask->len) {
    g_error("seqmask_ns: seq len (%lu) should equal mask len (%lu)",
	    seqs[0]->len, mask->len);
  }
  
  for(i = 0; i < mask->len; i++) {
    for(j = 0; j < n_seqs; j++) {
      if(seqs[j]->sym[i] == NUC_N) {
	mask->mask[i] = mask->mask[i] | mask_id;
	break;
      }
    }
  }

  return;
}


/*
 * Masks bases which are GAPs in any of the provided sequences
 */
void seqmask_gaps(SeqMask *mask, unsigned char mask_id, Seq **seqs,
		int n_seqs) {
  long i;
  int j;

  if(seqs[0]->len != mask->len) {
    g_error("seqmask_ns: seq len (%lu) should equal mask len (%lu)",
	    seqs[0]->len, mask->len);
  }
  
  for(i = 0; i < mask->len; i++) {
    for(j = 0; j < n_seqs; j++) {
      if(seqs[j]->sym[i] == NUC_GAP) {
	mask->mask[i] = mask->mask[i] | mask_id;
	break;
      }
    }
  }

  return;
}



/*
 * Given an array of aligned sequences masks all bases which are
 * adjacent to a mismatch, GAP or N, not including the mismatch itself
 * (unless, of course, it is also adjacent to a mismatch or gap).
 *
 * This allows only mismatches flanked by conserved bases to be
 * counted and is used to avoid questionable portions of an alignment.
 */
void seqmask_mismatch_adjacent(SeqMask *mask, unsigned char mask_id, 
			       Seq **seqs, int n_seqs,  int mask_n_adj, 
			       int mask_gap_adj, 
			       int mask_mismatch_adj) {

  seqmask_mismatch_near(mask, mask_id, seqs, n_seqs,  
			mask_n_adj, mask_gap_adj, mask_mismatch_adj, 1);
}




/**
 * Masks around alignment block ends
 */
void seqmask_aln_blk_ends(SeqMask *mask, unsigned char mask_id, 
			  GSList *aln_lst, int n_flnk) {

  long i, start, end;
  GSList *cur;
  Alignment *aln;
  
  cur = aln_lst;
  while(cur != NULL) {
    aln = cur->data;

    /* mask portion of alignment by block start */
    start = aln->coords[ALN_QUERY].start;
    end = start + n_flnk - 1;
    if(end > mask->len) {
      end = mask->len;
    }
    for(i = start-1; i < end; i++) {
      mask->mask[i] = mask->mask[i] | mask_id;
    }

    /* mask portion of sequence near block end  */
    end = aln->coords[ALN_QUERY].end;
    start = end - n_flnk + 1;
    if(start < 1) {
      start = 1;
    }
    for(i = start-1; i < end; i++) {
      mask->mask[i] = mask->mask[i] | mask_id;
    }

    cur = cur->next;
  }
}





/**
 * Masks around alignment block ends and internal gaps that are present in
 * either query or target sequence.
 */
void seqmask_aln_gaps(SeqMask *mask, unsigned char mask_id, 
		      GSList *aln_lst, int n_flnk) {

  long i,j, qpos, start, end;
  GSList *cur;
  Alignment *aln;
  unsigned char qnuc,tnuc;
  
  cur = aln_lst;
  while(cur != NULL) {
    aln = cur->data;

    /* mask portion of alignment by block start */
    start = aln->coords[ALN_QUERY].start;
    end = start + n_flnk - 1;
    if(end > mask->len) {
      end = mask->len;
    }
    for(i = start-1; i < end; i++) {
      mask->mask[i] = mask->mask[i] | mask_id;
    }

    /* mask portion of sequence near block end  */
    end = aln->coords[ALN_QUERY].end;
    start = end - n_flnk + 1;
    if(start < 1) {
      start = 1;
    }
    for(i = start-1; i < end; i++) {
      mask->mask[i] = mask->mask[i] | mask_id;
    }


    /* mask around gaps in the alignment */
    qpos = aln->coords[ALN_QUERY].start;
    for(i = 0; i < aln->len; i++) {
      qnuc = aln->sym[ALN_QUERY][i];
      tnuc = aln->sym[ALN_TARGET][i];

      if((tnuc == NUC_GAP) || (qnuc == NUC_GAP)) {
	/* mask around gap in query or target sequence */
	start = qpos-n_flnk;
	end   = qpos+n_flnk;

	if(start < 1) {
	  start = 1;
	}
	if(end > mask->len) {
	  end = mask->len;
	}
	for(j = start-1; j < end; j++) {
	  mask->mask[j] = mask->mask[j] | mask_id;
	}
      }

      if(qnuc != NUC_GAP) {
	/* update query sequence position */
	qpos++;
      }
    }

    cur = cur->next;
  }
}



/**
 * Masks around alignment block ends and gaps that are present in
 * either query or target sequence.
 */
void seqmask_aln_mismatch(SeqMask *mask, unsigned char mask_id, 
			  GSList *aln_lst, int n_flnk) {

  long i,j, qpos, start, end;
  GSList *cur;
  Alignment *aln;
  unsigned char qnuc,tnuc;
  
  cur = aln_lst;
  while(cur != NULL) {
    aln = cur->data;

    /* mask around mismatches in the alignment */
    qpos = aln->coords[ALN_QUERY].start;
    for(i = 0; i < aln->len; i++) {
      qnuc = aln->sym[ALN_QUERY][i];
      tnuc = aln->sym[ALN_TARGET][i];

      if((qnuc != NUC_GAP) && (tnuc != NUC_GAP)) {
	if(qnuc != tnuc) {
	  /* mask around mismatch */
	  start = qpos-n_flnk;
	  end   = qpos+n_flnk;

	  if(start < 1) {
	    start = 1;
	  }
	  if(end > mask->len) {
	    end = mask->len;
	  }

	  /* mask sites to left */
	  for(j = start-1; j < qpos-1; j++) {
	    mask->mask[j] = mask->mask[j] | mask_id;
	  }
	  /* mask sites to right */
	  for(j = qpos; j < end; j++) {
	    mask->mask[j] = mask->mask[j] | mask_id;
	  }
	}
      }

      if(qnuc != NUC_GAP) {
	/* update query sequence position */
	qpos++;
      }
    }

    cur = cur->next;
  }
}



/**
 * function used for debugging quality score masking
 */
static void write_qscore_region(unsigned char *qvals, long pos, long flnk) {
  long i;

  for(i = pos-flnk; i < pos; i++) {
    fprintf(stderr, "%2d ", qvals[i-1]);
  }
  fprintf(stderr, " %2d  ", qvals[pos-1]);

  for(i = pos+1; i < pos+flnk+1; i++) {
    fprintf(stderr, "%2d ", qvals[i-1]);
  }

  fprintf(stderr, "\n");
}


/**
 * Masks sites that have low sequence quality scores or are flanked by
 * low quality scores. 
 * 
 * qvals - array of quality scores, same length as provided mask
 * n_flnk - number of flanking sites checked for low scores.
 * min_flnk - min quality score for flanking sites
 * min_site - min quality score for site of interest
 *
 */
void seqmask_low_qscores(SeqMask *mask, unsigned char mask_id,
		    unsigned char *qvals, int n_flnk, 
		    unsigned char min_flnk, unsigned char min_site) {

  long i;
  int n_low_left;   /* number of low-scoring left-flanking  */
  int n_low_right;  /* number of low-scoring right-flanking */

  /* count number of low-scoring sites in left-flanking region 
   * and mask these first sites automatically
   */
  n_low_left = 0;
  for(i = 0; i < n_flnk && i < mask->len; i++) {
    if(qvals[i] < min_flnk) {
      n_low_left++;
    }
    mask->mask[i] = mask->mask[i] | mask_id;
  }
  /* count number of low-scoring sites in right-flanking bases */
  n_low_right = 0;
  for(i = n_flnk+1; (i < n_flnk+n_flnk+1) && (i < mask->len); i++) {
    if(qvals[i] < min_flnk) {
      n_low_right++;
    }
  }

  /* move along entire mask, masking sites that have any low-scoring
   * flanking sites, or that are below a certain threshold themselves
   */
  for(i = n_flnk; i < mask->len - n_flnk; i++) {
    if(n_low_right || n_low_left || qvals[i] < min_site) {
      mask->mask[i] = mask->mask[i] | mask_id;
    } 

    if(i < mask->len - n_flnk - 1) {
      /* update count of low-scoring flanking sites */
      if(qvals[i-n_flnk] < min_flnk) {
	/* old left-most flanking site was low scoring */
	n_low_left--;
      }
      if(qvals[i] < min_flnk) {
	/* new left flanking site is low scoring */
	n_low_left++;
      }
      if(qvals[i+1] < min_flnk) {
	/* old right flanking site (new center site) was low scoring */
	n_low_right--;
      }
      if(qvals[i+n_flnk+1] < min_flnk) {
	/* new right-most flanking site is low scoring */
	n_low_right++;
      }
    }
  }

  /* mask last sites in region automatically */
  for(i = mask->len-n_flnk; i < mask->len; i++) {
    mask->mask[i] = mask->mask[i] | mask_id;
  }

}


/*
 * Given an array of aligned sequences masks all bases which are within
 * dist bases from a mismatch, GAP or N, not including the mismatch itself
 * (unless, of course, it is near another mismatch or gap).
 *
 * This allows only mismatches flanked by conserved bases to be
 * counted and is used to avoid questionable portions of an alignment.
 */
void seqmask_mismatch_near(SeqMask *mask, unsigned char mask_id, 
			   Seq **seqs, int n_seqs,  int mask_n_adj,
			   int mask_gap_adj, 
			   int mask_mismatch_adj, int dist) {
  long i, start, end;
  int j;
  int is_mismatch, is_gap, is_n;
  
  if(n_seqs == 0) {
    return;
  }

  if(mask->len != seqs[0]->len) {
    g_error("seqmask_mismatch_near: seq len (%lu) should equal mask len (%lu)",
	    seqs[0]->len, mask->len);
  }

  for(j = 0; j < n_seqs; j++) {
    if(seqs[j]->len != mask->len) {
      g_error("seqmask_mismatch_adjacent: expected all sequences to be same "
	      "length as mask (mask=%lu, seq=%lu)", mask->len, seqs[j]->len);
    }
  }

  for(i = 1; i < mask->len-1; i++) {
    is_mismatch = FALSE;
    is_gap = FALSE;
    is_n = FALSE;
      
    for(j = 0; j < n_seqs; j++) {
      /* check if is GAP or N in any sequence, or if mismatch */
      if(seqs[j]->sym[i] == NUC_GAP) { 
	is_gap = TRUE;
      } 
      if(seqs[j]->sym[i] == NUC_N) {
	is_n = TRUE;
      } 
      if(j > 0 && (seqs[j]->sym[i] != seqs[j-1]->sym[i])) {
	is_mismatch = TRUE;
      }
    }

    if((mask_mismatch_adj && is_mismatch) ||
       (mask_n_adj && is_n) || (mask_gap_adj && is_gap)) {
      
      /* mask bases to left of mismatch/gap/N */
      start = (dist > i) ? 0 : i - dist;
      for(j = start; j < i; j++) {
	mask->mask[j] = mask->mask[j] | mask_id;
      }

      /* mask bases to right of mismatch/gap/N */
      end = i + dist;
      if(end >= mask->len) {
	end = mask->len-1;
      }
      for(j = i+1; j <= end; j++) {
	mask->mask[j] = mask->mask[j] | mask_id;
      }
    }
  }
}





/*
 * Given a transcript, flags positions of exons in the provided mask.
 */
void seqmask_exons(SeqMask *mask, unsigned char mask_id, Transcript *tr) {
  int i;

  for(i = 0; i < tr->num_exons; i++) {
    seqmask_region(mask, mask_id, &tr->exons[i]);
  }
}


/*
 * Given a transcript, flags positions of CDS in the provided mask.
 */
void seqmask_cds(SeqMask *mask, unsigned char mask_id, Transcript *tr) {
  int i, n;
  SeqCoord *cds;

  cds = transcript_cds(tr, &n);

  if(n > 0) {
    for(i = 0; i < n; i++) {
      seqmask_region(mask, mask_id, &cds[i]);
    }

    seq_coord_array_free(cds, n);
  }
}


/*
 * Given a transcript, flags positions of the 5' UTR in the provided
 * mask.
 */
void seqmask_utr5(SeqMask *mask, unsigned char mask_id, Transcript *tr) {
  int i,n;
  SeqCoord *utr5;

  utr5 = transcript_utr5(tr, &n);

  if(n > 0) {
    for(i = 0; i < n; i++) {
      seqmask_region(mask, mask_id, &utr5[i]);
    }
    seq_coord_array_free(utr5, n);
  }
}


/*
 * Given a transcript, flags positions of the 3' UTR in the provided
 * mask.
 */
void seqmask_utr3(SeqMask *mask, unsigned char mask_id, Transcript *tr) {
  int i,n;
  SeqCoord *utr3;

  utr3 = transcript_utr3(tr, &n);

  if(n > 0) {
    for(i = 0; i < n; i++) {
      seqmask_region(mask, mask_id, &utr3[i]);
    }
    seq_coord_array_free(utr3, n);
  }
}



/**
 * Given a transcript, flags the positions of the introns.
 */
void seqmask_introns(SeqMask *mask, unsigned char mask_id, Transcript *tr) {
  int i, n_introns;
  SeqCoord *introns;
  
  introns = transcript_introns(tr, &n_introns);

  for(i = 0; i < n_introns; i++) {
    seqmask_region(mask, mask_id, &introns[i]);
  }

  seq_coord_array_free(introns, n_introns);
}



/**
 * Given a transcript, masks the first intron.
 */
void seqmask_first_intron(SeqMask *mask, unsigned char mask_id, 
			  Transcript *tr) {
  SeqCoord *introns;
  int n_introns;

  introns = transcript_introns(tr, &n_introns);

  if(n_introns > 0) {
    seqmask_region(mask, mask_id, &introns[0]);
  }

  seq_coord_array_free(introns, n_introns);

  return;
}



/**
 * Given a transcript, masks all introns except first introns.
 */
void seqmask_non_first_introns(SeqMask *mask, unsigned char mask_id, 
			       Transcript *tr) {
  int i, n_introns;
  SeqCoord *introns;
  
  introns = transcript_introns(tr, &n_introns);

  for(i = 1; i < n_introns; i++) {
    seqmask_region(mask, mask_id, &introns[i]);
  }

  seq_coord_array_free(introns, n_introns);
}





/*
 * Masks a region specified by a provided SeqCoord.  Portions of the
 * provided region which are outside of the coordinate of the sequence
 * mask are ignored. Chromosome or seqname information is not taken
 * into consideration. Currently strand information is also ignored.
 */
void seqmask_region(SeqMask *mask, const unsigned char mask_id, 
		    const SeqCoord *region) {
  long i;
  long start, end;

  
  if(region->end < region->start) {
    g_error("seqmask_region: region end (%ld) must be >= region start (%ld)",
	    region->end, region->start);
  }

  if(region->end < mask->c.start) {
    /* entire region is before seqmask */
    return;
  }

  if(region->start > mask->c.end) {
    /* entire region is after seqmask */
    return;
  }

  if(region->start < mask->c.start) {
    start = 1;
  } else {
    start = region->start - mask->c.start + 1;
  }

  if(region->end > mask->c.end) {
    end = mask->len;
  } else {
    end = region->end - mask->c.start + 1;
  }

  if(start < 1 || end > mask->len) {
    g_error("seqmask_region: coords %ld-%ld out of array bounds 1-%ld",
	    start, end, mask->len);
  }

  for(i = start-1; i < end; i++) {
    mask->mask[i] = mask->mask[i] | mask_id;
  }
}


/* Convenience function to mask the regions represented by the
 * coordinates of the provided SeqFeatures
 */
void seqmask_features(SeqMask *mask, const unsigned char mask_id, 
		      const SeqFeature *sfeats, const long n_sfeat) {
  long i;
  for(i = 0; i < n_sfeat; i++) {
    seqmask_region(mask, mask_id, &sfeats[i].c);
  }
}



/*
 * Masks coordinates of features from the provided array which have
 * an attribute with the specified name and value.
 */
void seqmask_features_attrib(SeqMask *mask, const unsigned char mask_id,
			     const SeqFeature *sf, const long n_sf, 
			     const char *attr_name, const char *attr_val) {
  long i;
  char *val;
  int is_name;

  is_name = (strcmp("name", attr_name) == 0);

  for(i = 0; i < n_sf; i++) {

    val = seqfeat_get_attrib_str(&sf[i], attr_name);
    
    if(val == NULL && is_name) {
      /* if the attrib name is "name" check if name is set 
       * instead of an actual attribute
       */
      if(sf[i].name && (strcmp(sf[i].name, attr_val)==0)) {
	seqmask_region(mask, mask_id, &sf[i].c);
      }
    } else {
      if((strcmp(val, attr_val) == 0)) {
	seqmask_region(mask, mask_id, &sf[i].c);
      }
    }
  }  
}




/*
 * Masks segments which are duplicated elsewhere in the genome
 * by reading their coordinates and lengths from a file
 */
void seqmask_dupsegs(SeqMask *mask, unsigned char mask_id, gchar *file) {
  FILE *fh;
  gchar *line, **tok;
  long start, end, len, i;

  fh = fopen(file, "r");
  if(fh == NULL) {
    g_error("seqmask_dupsegs: could not open file '%s'\n", file);
  }

  /* skip header */
  line = util_fgets_line(fh);
  g_free(line);

  while((line = util_fgets_line(fh)) != NULL) {
    tok = g_strsplit(line, " ", 2);

    if(tok[0] == NULL || tok[1] == NULL) {
      g_error("seqmask_dupsegs: expected two tokens per line");
    }
        
    start = strtol(tok[0], NULL, 10);
    len   = strtol(tok[1], NULL, 10);

    end = start + len - 1;

    for(i = start-1; i < end; i++) {
      mask->mask[i] = mask->mask[i] | mask_id;
    }

    g_strfreev(tok);
    g_free(line);
  }

  fclose(fh);

}




/**
 * If near_or_far is set to DIST_NEAR, sites that are within
 * dist_thresh bases of sites flagged with site_mask_id are masked.
 * If near_or_far is set to DIST_FAR, then sites that are farther away
 * than dist_thresh bases are masked.
 */
void seqmask_sites_at_dist(SeqMask *mask, unsigned char mask_id, 
			   unsigned char site_mask_id, 
			   long dist_thresh,
			   int near_or_far) {

  long i;
  int *dists;

  /* calculate distance to masked sites */
  fprintf(stderr, "calcing dists\n");
  dists = dist_calc_dists(mask, site_mask_id);

  if(near_or_far == DIST_NEAR) {
    fprintf(stderr, "masking sites within %ld bp of sites with id %d\n",
	    dist_thresh, site_mask_id);
    for(i = 0; i < mask->len; i++) {
      if(dists[i] <= dist_thresh) {
	mask->mask[i] = mask->mask[i] | mask_id;
      }
    }
  } 
  else if(near_or_far == DIST_FAR) {
    for(i = 0; i < mask->len; i++) {
      if(dists[i] > dist_thresh) {
	mask->mask[i] = mask->mask[i] | mask_id;
      }
    }
  } else {
    g_error("seqmask_sites_at_dist: near_or_far must be "
	    "DIST_FAR or DIST_NEAR");
  }

  g_free(dists);

  return;
}









static void mask_ffd_codon_sites(Seq **seqs, const int n_seq, 
				 const short strand, const long cd2gn_pos[3],
				 SeqMask *mask, const unsigned char mask_id) {
  unsigned char **codons, nuc;
  int is_ffd, degen, i, j;
  long pos, idx;
  long seq_start;

  if(n_seq < 1) {
    return;
  }


  seq_start = seqs[0]->c.start;
  /* create array of codons, one for each sequence */
  codons = g_new(unsigned char *, n_seq);
  for(i = 0; i < n_seq; i++) {
    codons[i] = g_new(unsigned char, 3);

    for(j = 0; j < 3; j++) {
      /* get genomic position and set codon nucleotide */
      pos = cd2gn_pos[j];
      idx = pos - seq_start;
      if((idx < 0) || (idx > seqs[i]->len)) {
        g_error("%s:%d: invalid sequence index %ld", __FILE__, __LINE__, idx);
      }
      if(strand == STRAND_FWD) {
	nuc =  seqs[i]->sym[idx];
      } else { 
	nuc = nuc_comp(seqs[i]->sym[idx]);
      }

      /* skip if codon from any sequence contains N or GAP */
      if((nuc == NUC_N) || (nuc == NUC_GAP)) {
	for(j = 0; j <= i; j++) {
	  g_free(codons[j]);
	}
	g_free(codons);
	return;
      }

      codons[i][j] = nuc;
    }
  }



  /* now check degeneracy of each position
   * only calling a site four-fold-degenerate if
   * it has degeneracy=3 in all sequences
   */
  
  for(i = 0; i < 3; i++) {
    is_ffd = TRUE;

    for(j = 0; j < n_seq; j++) {
      degen = aa_codon_degeneracy(codons[j], i);
      if(degen < 3) {
	is_ffd = FALSE;
	break;
      }
    }

    if(is_ffd) {
      /* this is a four-fold-degenerate site in all sequences */
      pos = cd2gn_pos[i];
      idx = pos - seq_start;
      mask->mask[idx] = mask->mask[idx] | mask_id;

      /* fprintf(stderr, "FFD site: codon_pos: %d, genome_pos: %ld, idx=%ld\n",
       *	      i+1, pos, idx);
       */
    }
  }


  /* free codon array */
  for(i = 0; i < n_seq; i++) {
    g_free(codons[i]);
  }
  g_free(codons);
}







/**
 * Masks four-fold-degenerate sites for provided transcript. Only sites
 * which are FFD in ALL provided sequences are masked.
 */
void seqmask_ffd_sites(SeqMask *mask, const unsigned char mask_id,
		       Seq **seqs, const int n_seq, Transcript *tr) {
  SeqCoord *cds;
  int cd_idx, gn_pos, n_cds_coord, cur_coord, strand;
  long cd2gn_pos[3];
 
  if(n_seq < 1) {
    return;
  }

  /* get array of genomic coordinates representing CDS */
  cds = transcript_cds(tr, &n_cds_coord);
  
  if(!cds || n_cds_coord < 1) {
    return;
  }

  strand = tr->c.strand;
  cur_coord = 0;
  cd_idx = 0;

  if(strand == STRAND_FWD) {
    gn_pos = cds[cur_coord].start;
  } else {
    gn_pos = cds[cur_coord].end;
  }

  while(TRUE) {
    /* maintain map of codon -> genome coordinates */
    cd2gn_pos[cd_idx] = gn_pos;

    /* update genomic position */
    if(strand == STRAND_FWD) {
      gn_pos++;
      if(gn_pos > cds[cur_coord].end) {
	/* move ahead to next exon */
	cur_coord++;
	if(cur_coord >= n_cds_coord) {
	  /* terminate, at end of CDS */
	  break;
	}
	gn_pos = cds[cur_coord].start;
      }
    } else {
      gn_pos--;
      if(gn_pos < cds[cur_coord].start) {
	/* advance to next exon */
	cur_coord++;
	if(cur_coord >= n_cds_coord) {
	  /* already at last exon in CDS */
	  break;
	}
	gn_pos = cds[cur_coord].end;
      }
    }

    /* update codon position */
    cd_idx++;
      
    if(cd_idx > 2) {
      /* mask any four-fold degenerate sites in this codon */
      mask_ffd_codon_sites(seqs, n_seq, strand, cd2gn_pos, mask, mask_id);

      /* start a new codon */
      cd_idx = 0;      
    }
  }

  if(cd_idx != 2) {
    /* sanity check: CDS should end on third position of codon */
    g_error("%s:%d: expected fwd strand CDS to end on third codon position", 
	    __FILE__, __LINE__);
  }

  seq_coord_array_free(cds, n_cds_coord);
}
