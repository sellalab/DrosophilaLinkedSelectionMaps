
#include "rectab.h"

#include "seqfeat.h"



/**
 * Creates a table from a set of features that can be used to obtain
 * recombination rates or distances at any position on a chromomsome.
 * If the features are reasonably large, this uses much less memory
 * than creating an array of recombination distances or rates.
 */
RecRateTable *rectab_from_feats(SeqFeature *orig_sf, long n_orig_sf, 
				long chr_len, double scale) {
  RecRateTable *rt;
  RecRateBlock *b, *prev_b;
  SeqFeature *sf;
  long i, len, sf1_len, sf2_len, last_end, n_sf;

  if(n_orig_sf < 1) {
    return NULL;
  }

  /* sort features */
  qsort(orig_sf, n_orig_sf, sizeof(SeqFeature), seqfeat_cmp_nostrand);

  sf = g_new(SeqFeature, n_orig_sf);

  /* Overlapping blocks create problems so eliminate them. When there
   * are overlapping blocks arbitrarily take longer block.  Presumably
   * these exist because of the mapping process from NCBI35->36, or
   * because of some bug at UCSC.
   */
  last_end = 0;
  n_sf = 0;
  for(i = 0; i < n_orig_sf; i++) {
    if(orig_sf[i].c.start < last_end) {
      g_warning("rectab_from_feats: resolving overlapping rec features\n");

      if(n_sf < 1) {
	g_error("rectab_from_feats: first feature has negative coord?");
      }

      sf1_len = sf[n_sf-1].c.end - sf[n_sf-1].c.start + 1;
      sf2_len = orig_sf[i].c.end - orig_sf[i].c.start + 1;

      if(sf2_len > sf1_len) {
	/* replace previous feature with this one */
	sf[n_sf-1].c.start = orig_sf[i].c.start;
	sf[n_sf-1].c.end = orig_sf[i].c.end;
	sf[n_sf-1].score = orig_sf[i].score;

	/* this could potentially overlap with even earlier
	 * feature now. In this case just shorten this feature
	 */
	if(n_sf > 1) {
	  g_warning("rectab_from_feats: resolving double overlap\n");
	  if(sf[n_sf-2].c.end >= sf[n_sf-1].c.start) {
	    sf[n_sf-1].c.start = sf[n_sf-2].c.end+1;
	  }
	}
	last_end = sf[n_sf-1].c.end;
      } else {
	/* keep previous feature, nothing needs doing */
      }
    } else {
      /* no overlap with previous, keep this feature (for now at least) */
      sf[n_sf].c.start = orig_sf[i].c.start;
      sf[n_sf].c.end   = orig_sf[i].c.end;
      sf[n_sf].score = orig_sf[i].score;
      last_end = sf[n_sf].c.end;
      n_sf += 1;
    }
  }


  rt = g_new(RecRateTable, 1);
  rt->n = n_sf;

  /* count number of extra blocks that will be needed */
  last_end = 0;
  for(i = 0; i < n_sf; i++) {
    if(sf[i].c.start > last_end+1) {
      /* need an extra block to fill gap */
      rt->n += 1;
    }
    last_end = sf[i].c.end;
  }

  if(last_end < chr_len) {
    /* need extra block for gap at end of chr */
    rt->n += 1;
  }

  /* allocate space for blocks */
  rt->blk = g_new(RecRateBlock, rt->n);
  rt->cur = 0;

  if(sf[0].c.start > 1) {
    /* there is a gap at the beginning of chr, assume same rec-rate
     * as first block
     */
    b = rt->blk;
	
    b->start = 1;
    b->end  = sf[0].c.start-1;
    len = b->end - b->start + 1;
    b->rate = sf[0].score * scale;
    b->r_start = 0.0;
    b->r_end = (double)(len-1) * b->rate;
    rt->cur = 1;
  }  else {
    rt->cur = 0;
  }


  for(i = 0; i < n_sf; i++) {
    if(rt->cur > 0 && sf[i].c.start > rt->blk[rt->cur-1].end+1) {
      /* need to add a block to fill gap */
      b = &rt->blk[rt->cur];
      prev_b = &rt->blk[rt->cur-1];

      if(i < n_sf-1) {
	/* take avg of prev blk and next */
	rt->blk[rt->cur].rate = (prev_b->rate + (sf[i+1].score*scale)) * 0.5;
      } else {
	/* this is gap at end of chr, just use prev rate */
	b->rate = prev_b->rate;
      }

      b->start = prev_b->end + 1;
      b->end = sf[i].c.start-1;
      len = b->end - b->start + 1;

      b->r_start = prev_b->r_end + prev_b->rate;
      b->r_end = b->r_start + (b->rate * (double)(len-1));

      rt->cur += 1;      
    }

    b = &rt->blk[rt->cur];    
    b->start = sf[i].c.start;
    b->end = sf[i].c.end;
    b->rate = sf[i].score * scale;
    len = b->end - b->start + 1;

    if(rt->cur > 0) {
      prev_b = &rt->blk[rt->cur-1];
      b->r_start = prev_b->r_end + prev_b->rate;
    } else {
      b->r_start = 0.0;
    }

    b->r_end = b->r_start + (b->rate * (double)(len-1));
    rt->cur += 1;
  }


  if(rt->blk[rt->cur-1].end < chr_len) {
    /* create final block to fill up to end of chr */
    b = &rt->blk[rt->cur];
    prev_b = &rt->blk[rt->cur-1];
    b->rate = prev_b->rate;
    b->start = prev_b->end + 1;
    b->end = chr_len;
    len = b->end - b->start + 1;

    b->r_start = prev_b->r_end + prev_b->rate;
    b->r_end = b->r_start + (b->rate * (double)(len-1));

    rt->cur += 1;
  }


  if(rt->cur != rt->n) {
    g_error("rectab_from_feats: expected %ld features, but got %ld",
	    rt->n, rt->cur);
  }

  rt->chr_r_len = rt->blk[rt->n-1].r_end;

  /* set cur to beginning of chr */
  rt->cur = 0;


  g_free(sf);

  return rt;
}


/* 
 * Helper function: moves cur blk index to the block
 * that overlaps with the provided position.
 */
static void find_cur_rblk(RecRateTable *rt, long pos) {

  /* TODO: If start in the wrong block, could do a quick estimation
   * of what block we should be near based on the number of blocks and
   * the size of the chr.
   */

  /* move back along blocks */
  while(rt->cur > 0 && rt->blk[rt->cur].start > pos) {
    rt->cur -= 1;
  }

  /* move fwd along blocks */
  while((rt->cur < (rt->n - 1)) && (rt->blk[rt->cur].end < pos)) {
    rt->cur += 1;
  }

  if(rt->blk[rt->cur].start > pos || rt->blk[rt->cur].end < pos) {
    g_error("rectab_r_pos: could not find recombination block overlapping "
	    "chromosome position %ld", pos);
  }

}

/*
 * Returns the recombination distance in morgans for a given physical
 * distance along the chromosome.
 */
double rectab_rpos(RecRateTable *rt, long pos) {
  RecRateBlock *b;

  find_cur_rblk(rt, pos);

  b = &rt->blk[rt->cur];
  return (double)(pos - b->start)*b->rate + b->r_start;
}


/**
 * Returns the recombination rate at the given physical position
 */
double rectab_rate(RecRateTable *rt, long pos) {
  find_cur_rblk(rt, pos);
  return rt->blk[rt->cur].rate;
}



/**
 * Frees memory allocated for a recombination rate table
 */
void rectab_free(RecRateTable *rt) {
  g_free(rt->blk);
  g_free(rt);
}
