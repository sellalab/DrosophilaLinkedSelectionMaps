
#include <glib.h>
#include <math.h>

#include "dist.h"
#include "seqfeat.h"

#include "util.h"
#include <zlib.h>

#define MAX_LINE 1024

/**
 * Returns an array of floats the length of a sequence. Each number
 * represents the physical distance to the nearest site that is masked
 * in the provided seqmask.  These are returned as floats rather than
 * ints, so that the return value is consistant with the
 * dist_calc_genet_dists function.
 */
float *dist_calc_dists_float(SeqMask *mask, unsigned char mask_id) {
  float *dist_array;
  int dist;
  int seen_first;
  long i;

  dist_array = g_new(float, mask->len);

  /* go through the scores fwd direction setting distances to the
   * distance to the previous conserved base.
   */
  seen_first = FALSE;
  dist = G_MAXINT; /* don't know how far away until we see first conserved */
  for(i = 0; i < mask->len; i++) {
    if(mask->mask[i] & mask_id) {
      /* we hit a conserved base */
      dist = 0;
      seen_first = TRUE;
    } else if(seen_first) {
      dist++;
    }

    dist_array[i] = (float)dist;
  }

  /* go through the scores in rev direction, updating those that are
   * closer to subsequent conserved bases
   */
  seen_first = FALSE;
  for(i = mask->len-1; i >= 0; i--) {
    if(mask->mask[i] & mask_id) {
      /* this is a masked site */
      dist = 0;
      seen_first = TRUE;
    } else if(seen_first) {
      dist++;
    }

    if(seen_first && dist_array[i] > dist) {
      dist_array[i] = (float)dist;
    }
  }

  return dist_array;
}



/**
 * Returns an array of integers that give, for each site, the physical
 * distance to the nearest masked site in bp. If use_neg_dists is
 * TRUE, positive and negative distance integers are
 * returned. Positive integers indicate that the nearest masked site
 * is upstream, negative integers indicate that the nearest masked
 * site is downstream. If use_neg_dists is FALSE, only positive
 * (absolute) distances are returned.
 */
int *dist_calc_dists_sign(SeqMask *mask, unsigned char mask_id, 
			  int use_neg_dists) {

  int *dist_array;
  int dist;
  int seen_first;
  long i;

  dist_array = g_new(int, mask->len);

  /* go through the scores fwd direction setting distances to the
   * distance to the previous conserved base.
   */
  seen_first = FALSE;
  dist = DIST_NA; /* don't know how far away until we see first conserved */
  for(i = 0; i < mask->len; i++) {
    if(mask->mask[i] & mask_id) {
      /* we hit a conserved base */
      dist = 0;
      seen_first = TRUE;
    } else if(seen_first) {
      dist++;
    }

    dist_array[i] = dist;
  }

  /* go through the scores in rev direction, updating those that are
   * closer to subsequent conserved bases
   */
  seen_first = FALSE;
  for(i = mask->len-1; i >= 0; i--) {
    if(mask->mask[i] & mask_id) {
      /* this is a masked site */
      dist = 0;
      seen_first = TRUE;
    } else if(seen_first) {
      dist++;
    }

    if(seen_first && dist_array[i] > dist) {
      if(use_neg_dists) {
	dist_array[i] = -dist;
      } else {
	dist_array[i] = dist;
      }
    }
  }

  return dist_array;
  
}


/**
 * Returns an array of ints the length of a sequence. Each number
 * represents the physical distance to the nearest site that is masked
 * in the provided seqmask.
 */
int *dist_calc_dists(SeqMask *mask, unsigned char mask_id) {
  return dist_calc_dists_sign(mask, mask_id, FALSE);
}



/**
 * Returns an array of floats representing the genetic distance of
 * each base to the nearest masked base in 10^-8 M (10 nM).
 */
float *dist_calc_genet_dists(SeqMask *mask, unsigned char mask_id,
			      float *recomb_rates) {
  float *dist_array, dist, last_dist;
  int seen_first;
  long i;

  dist_array = g_new(float, mask->len);

  /* go through the scores fwd direction setting distances to the
   * distance to the previous conserved base.
   */
  seen_first = FALSE;
  dist = G_MAXFLOAT; /* don't know how far away until we see first conserved */
  for(i = 0; i < mask->len; i++) {
    if(mask->mask[i] & mask_id) {
      /* we hit a conserved base */
      dist = 0;
      seen_first = TRUE;
    } else if(seen_first) {
      if(recomb_rates[i] == DIST_NA) {
	dist = DIST_NA;
      } else {
	dist += recomb_rates[i];
      }
    }

    dist_array[i] = dist;
  }

  /* go through the scores in rev direction, updating those that are
   * closer to subsequent conserved bases
   */
  seen_first = FALSE;
  last_dist = 0;
  for(i = mask->len-1; i >= 0; i--) {
    if(mask->mask[i] & mask_id) {
      /* this is a conserved base */
      dist = 0;
      seen_first = TRUE;
    } else if(seen_first) {
      if(recomb_rates[i] == DIST_NA) {
	if(dist != DIST_NA) { 
	  /* we are entering a region with undefined recombination rates
	   * record the last defined distance from this direction so we
	   * have a lower bound on the distance
	   */
	  last_dist = dist;
	}
	dist = DIST_NA;
      } else {
	dist += recomb_rates[i];
      }
    }

    if(seen_first) {
      if(dist == DIST_NA) {
	/* all we know about distance from this direction is that
	 * it must be > last_dist
	 */
	if(dist_array[i] > last_dist) {
	  /* we don't know which direction had closer conserved base
	   * so set dist to undefined
	   */
	  dist_array[i] = DIST_NA;
	}
      } else {
	if(dist_array[i] > dist) {
	  /* conserved base is closer from this direction  */
	  dist_array[i] = dist;
	}
      }
    }
  }
  
  return dist_array;
}




/**
 * Reads HK (B) values from the specified file. Even though these are
 * stored as integers, they are returned as floats for consistancy
 * with the recombination-rate distance-based methods that this
 * program also uses.
 */
float *dist_read_hk(gchar *hk_file, long seq_len) {
  gzFile *gzf;
  float *dists, hk;
  char line[MAX_LINE];
  char *l_ptr;
  long start, end, len, i, line_num;

  gzf = gzopen(hk_file, "rb");
  if(gzf == NULL) {
    g_error("read_hk_dists: could not open file '%s'", hk_file);
  }

  dists = g_new(float, seq_len);

  start = 1;
  end = 0;
  line_num = 0;
  while((gzgets(gzf, line, MAX_LINE) != NULL)) {
    line_num++;

    /* skip comment and header lines */
    if(line[0] == '#' || line[0] == 'f') {
      continue;
    }


    /* read b-value and len from line */
    hk = (float)strtod(line, &l_ptr);
    len = strtol(l_ptr, NULL, 10);

    end = start + len - 1;

    if(end > seq_len) {
      g_error("read_hk_dists: end coordinate (%lu) exceeds number of"
	      " bases in chromosome (%lu), line num:%lu", end, seq_len,
	      line_num);
    }

    for(i = start-1; i < end; i++) {
      dists[i] = hk;
    }
    start = end+1;
  }

  if(end != seq_len) {
    g_error("read_hk_dists: end coordinate (%ld) does not equal "
	    "length of sequence (%ld)", end, seq_len);
  }

  gzclose(gzf);

  return dists;
}


/**
 * Calculates combined HK values by multiplying together two other
 * sets of HK values. To save memory, one of the arrays is freed and
 * the values in the other one are replaced (and it is returned)
 */
float *dist_combine_hk(float *hk1, float u_scale1, 
		       float *hk2, float u_scale2, 
		       long seq_len, float max_hk) {
  long i;
  float b1, b2;
  float rescale;

  rescale = max_hk / 1000000.0;

  for(i = 0; i < seq_len; i++) {
    if(u_scale1 != 1.0 || u_scale2 != 1.0) {
      /* convert from range 0-999 back to 0-1.0 */
      b1 = (hk1[i]+0.5)/1000.0; 
      b2 = (hk2[i]+0.5)/1000.0;
      
      /* combine two B values, scaling them by their deleterious
       * mutation rate factors, and rescale again so between 0-max_hk
       */
      hk1[i] = exp(log(b1)*u_scale1 + log(b2)*u_scale2) * max_hk;
    } else {
      hk1[i] = (hk1[i]+0.5)*(hk2[i]+0.5) * rescale;
    }
  }
  g_free(hk2);

  return hk1;
}



/**
 * Reads recombination rates from a feature-style file where the
 * feature score is the rate in 10-8M and the coordinates the region.
 *
 * The returned rate array should be freed when it is no longer
 * needed.
 */
float *dist_read_recomb_rate_feats(char *filename, long seq_len) {
  long i,j, cur_pos, n_regions;
  SeqFeature *rec_regions;
  float *rates;

  rates = g_new(float, seq_len);
  
  /* read recombination rates from file */
  rec_regions = seqfeat_read_flatfile(filename, &n_regions);

  /** assume that regions in file are ordered */
  cur_pos = 1;
  
  for(i = 0; i < n_regions; i++) {
    /* fill in gaps between regions with undefined recombination rates */
    for(j = cur_pos; j < rec_regions[i].c.start; j++) {
      rates[j-1] = DIST_NA;
    }

    /* fill in recombination rate of region */
    for(j = rec_regions[i].c.start; j < rec_regions[i].c.end; j++) {
      rates[j-1] = rec_regions[i].score;
    }

    cur_pos = rec_regions[i].c.end+1;
  }

  /* set remaining region at chr end to undefined rec rates */
  for(j = cur_pos; j < seq_len; j++) {
    rates[j-1] = DIST_NA;
  }

  seqfeat_array_free(rec_regions, n_regions);

  return rates;
}




/**
 * Converts recombination rates to recombination distance along the
 * chromosome (i.e. cumulative recombination rate). Regions with
 * undefined rates are considered to have a rate of 0.  To save
 * memory, the same array is used to store the distances (i.e.  the
 * values in the passed rate array are replaced)
 */
void dist_recomb_rates_to_dists(double *rates, long seq_len) {
  long i;
  double ttl_r, cur_r;

  ttl_r = 0;
  
  for(i = 0; i < seq_len; i++) {
    cur_r = rates[i];
    if(cur_r == DIST_NA) {
      cur_r = 0.0;
    }
    ttl_r += cur_r;

    /* set value in array to the total chromosome recombination distance */
    rates[i] = ttl_r;
  }
}



/**
 * Converts a provided array of features to an array of recombination
 * rates the length of a chromosome where each position represents a
 * physical position (bp) along the chromosome.
 *
 * The per-nucleotide rate is taken as the feature score muiltiplied
 * by the scale argument. If the interpolate flag is TRUE, regions
 * where there are no features are assumed to have a recombination
 * rate that is the average of the two flanking features (or single
 * flanking feature on the edge of the chromosome).
 *
 */
double *dist_recomb_rates_from_feats(SeqFeature *sf, long n_sf, long seq_len,
				     double scale, int interpolate) {
  long i, pos, start, end;
  double *rec_rates, rate;

  rec_rates = g_new(double, seq_len);

  if(n_sf == 0) {
    g_warning("dist_recomb_rates_from_feats: no recombination "
	      "features provided, assuming rec rate = 0.0");

    /* fill in region, assuming no recombination */
    for(pos = 1; pos <= seq_len; pos++) {
      rec_rates[pos-1] = 0.0;
    }

    return rec_rates;
  }

  /* sort features */
  qsort(sf, n_sf, sizeof(SeqFeature), seqfeat_cmp_nostrand);
  
  pos = 1; /* current position along chromosome */

  for(i = 0; i < n_sf; i++) {
    if(sf[i].c.start > pos) {
      /* recombination rate is not defined in this region */

      start = pos;
      end = sf[i].c.start-1;

      if(end >= seq_len) {
	g_warning("dist_recomb_rates_from_feats: feature start (%ld) is "
		  "> seq_len (%ld)", sf[i].c.start, seq_len);
	end = seq_len;
      }

      if(interpolate) {
	if(i == 0) {
	  /* use the rate of the first defined feature */
	  rate = sf[i].score * scale;
	} 
	else {
	  /* use the average of the previous and next defined feature */
	  rate = ((sf[i-1].score + sf[i].score) / 2.0) * scale;
	}

	if(rate < 0.0) {
	  g_error("dist_recomb_rates_from_feats: recombination rate (%g)"
		  "in region %ld-%ld is less than 0.0", rate, start, end);
	}
      } else {
	/* flag rates as unknown */
	rate = DIST_NA;
      }

      /* fill in rec rates up to current feature */
      for(pos = start; pos <= end; pos++) {
	rec_rates[pos-1] = rate;
      }
      
    }

    /* fill in rec dists for region defined by feature */
    rate = sf[i].score * scale;

    start = sf[i].c.start;
    end = sf[i].c.end;
    if(start < 1) {
      g_warning("dist_recomb_rates_from_feats: feature start (%ld) is < 1",
		sf[i].c.start);
      start = 1;
    }
    if(end > seq_len) {
      g_warning("dist_recomb_rates_from_feats: feature start (%ld) is "
		"> seq_len (%ld)", sf[i].c.start, seq_len);
      end = seq_len;
    }

    if(rate < 0.0) {
      g_error("dist_recomb_rates_from_feats: recombination rate (%g) "
	      "in region %ld-%ld is less than 0.0", rate, start, end);
    }

    for(pos = start; pos <= end; pos++) {
      rec_rates[pos-1] = rate;
    }
  }

  /* fill in rates to end of chr if not defined */
  if(pos < seq_len) {
    if(interpolate) {
      /* use the rate of the last defined feature */
      rate = sf[n_sf-1].score * scale;

      if(rate < 0.0) {
	g_error("dist_recomb_rates_from_feats: recombination rate (%g) "
		"in region %ld-%ld is less than 0.0", 
		rate, pos, seq_len);
      }
    } else {
      rate = DIST_NA;
    }

    while(pos < seq_len) {
      rec_rates[pos-1] = rate;
      pos++;
    }
  }

  return rec_rates;
}




/**
 * Reads recombination rates from a file.
 */
float *dist_read_recomb_rates(char *filename,long seq_len, double scale) {
  long i, prev, start, end;
  FILE *fh;
  int n_tok;
  char **tok, *line;
  float *rates, rate;

  rates = g_new(float, seq_len);

  fh = fopen(filename, "r");
  if(fh == NULL) {
    g_error("read_recomb_rates: could not read file '%s'", filename);
  }

  
  prev = 1;
  end  = 1;
  while((line = util_fgets_line(fh)) != NULL) {
    tok = g_strsplit(line, "\t", 5);

    n_tok = 0;
    while(tok[n_tok] != NULL) { n_tok++; }
    if(n_tok != 5) {
      g_error("read_recomb_rates: expected 5 tokens per line, got %d", n_tok);
    }

    start = strtol(tok[2], NULL, 10) + 1; /* 0-based UCSC coords */
    end   = strtol(tok[3], NULL, 10);
    rate  = strtod(tok[4], NULL);

    if(prev < start) {
      /* undefined region before this, pad with flags */
      for(i = prev-1; i < start; i++) {
	rates[i] = DIST_NA;
      }
    }
    for(i = start-1; i < end; i++) {
      rates[i] = rate * scale;
    }

    g_strfreev(tok);
    g_free(line);
    
    prev = end+1;
  }

  /* pad remainder of chromomsome with flag */
  for(i = end; i < seq_len; i++) {
    rates[i] = DIST_NA;
  }
  
  fclose(fh);

  return rates;
}
