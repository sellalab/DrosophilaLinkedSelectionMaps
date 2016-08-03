#ifndef __SEQMASK_H__
#define __SEQMASK_H__

#include <math.h>
#include <glib.h>

#include "transcript.h"
#include "seq.h"
#include "seqcoord.h"
#include "seqfeat.h"

#define SEQMASK_NUM_IDS 8
#define SEQMASK_CLEAR 0x00


typedef struct {
  int used_ids[SEQMASK_NUM_IDS];
  unsigned char *mask;
  long len;
  SeqCoord c;
} SeqMask;


/*
 * Returns a new sequence mask for a given length of sequence. The
 * mask must be freed when it is no longer needed.
 */
SeqMask *seqmask_new(SeqCoord *coord);
void seqmask_move(SeqMask *mask, SeqCoord *coord);
void seqmask_free(SeqMask *mask);
unsigned char seqmask_get_id(SeqMask *mask);

SeqMask *seqmask_submask(SeqMask *mask, SeqCoord *coord);

void seqmask_invert_mask(SeqMask *mask, unsigned char old_id, 
			 unsigned char new_id);

void seqmask_soft_repeats(SeqMask *mask, unsigned char mask_id, 
			  char *seqstr);
void seqmask_non_soft_repeats(SeqMask *mask, unsigned char mask_id, 
			      char *seqstr);

void seqmask_gc_flanking(SeqMask *mask, unsigned char mask_id, Seq *seq, 
			int num_flanking, int min_gc, int max_gc);

void seqmask_exons(SeqMask *mask, unsigned char mask_id, Transcript *tr);
void seqmask_cds(SeqMask *mask, unsigned char mask_id, Transcript *tr);
void seqmask_utr5(SeqMask *mask, unsigned char mask_id, Transcript *tr);
void seqmask_utr3(SeqMask *mask, unsigned char mask_id, Transcript *tr);
void seqmask_introns(SeqMask *mask, unsigned char mask_id, Transcript *tr);
void seqmask_first_intron(SeqMask *mask, unsigned char mask_id, 
			  Transcript *tr);
void seqmask_non_first_introns(SeqMask *mask, unsigned char mask_id, 
			       Transcript *tr);



void seqmask_ns(SeqMask *mask, unsigned char mask_id, Seq **seqs,
		int num_seqs);

void seqmask_gaps(SeqMask *mask, unsigned char mask_id, Seq **seqs,
		  int num_seqs);

void seqmask_cpg(SeqMask *mask, unsigned char mask_id, Seq **seqs, 
		 int num_seqs);

void seqmask_cpg_context(SeqMask *mask, unsigned char mask_id, Seq *seq);

void seqmask_non_cpg(SeqMask *mask,unsigned char mask_id, Seq **seqs,
		     int num_seqs);

void seqmask_deam(SeqMask *mask, unsigned char mask_id, Seq **seqs, 
		  int num_seqs);



void seqmask_aln_blk_ends(SeqMask *mask, unsigned char mask_id,
			  GSList *aln_lst, int n_flnk);


void seqmask_aln_gaps(SeqMask *mask, unsigned char mask_id,
		      GSList *aln_lst, int n_flnk);


void seqmask_aln_mismatch(SeqMask *mask, unsigned char mask_id,
		      GSList *aln_lst, int n_flnk);


void seqmask_mismatch_adjacent(SeqMask *mask, unsigned char mask_id, 
			       Seq **seqs, int num_seqs,  
			       int mask_n_adj, 
			       int mask_gap_adj, 
			       int mask_mismatch_adj);

void seqmask_mismatch_near(SeqMask *mask, unsigned char mask_id, 
			   Seq **seqs, int num_seqs,  int mask_n_adj, 
			   int mask_gap_adj, 
			   int mask_mismatch_adj, int dist);

void seqmask_low_qscores(SeqMask *mask, unsigned char mask_id,
		       unsigned char *qvals, int n_flnk, 
		       unsigned char min_flnk, unsigned char min_site);

void seqmask_region(SeqMask *mask, const unsigned char mask_id, 
		    const SeqCoord *coord);


void seqmask_features(SeqMask *mask, const unsigned char mask_id, 
		      const SeqFeature *sfeats, const long n_sfeat);


void seqmask_features_attrib(SeqMask *mask, const unsigned char mask_id, 
			    const SeqFeature *sfeats, const long n_sfeat,
			    const char *attrib_name, const char *attrib_val);



void seqmask_dupsegs(SeqMask *mask, unsigned char mask_id, char *file);

void seqmask_sites_at_dist(SeqMask *mask, unsigned char site_mask_id,
			   unsigned char mask_id, long dist_thresh,
			   int near_or_far);

void seqmask_ffd_sites(SeqMask *mask, const unsigned char mask_id,
		       Seq **seqs, const int n_seq, Transcript *tr);


#endif
