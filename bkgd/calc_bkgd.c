#include <math.h>
#include <string.h>

#include <glib.h>

#include "config.h"
#include "util.h"
#include "dist.h"
#include "rectab.h"

//#include <genome_db.h>

//#include <analysis.h>

#include "bkgd.h"
#include "bkgd_param.h"
#include "bkgd_interp.h"
#include "bkgd_files.h"



/**
 * Retrieves a table of recombination rates for current chromosome
 * from database
 */
static RecRateTable *get_rectab(Config *config, Chromosome *chr) {
  long n_sf;
//  SeqCoord *region;
//  SeqFeatureDBA *sf_dba;
  SeqFeature *sf;
  char *rec_tab_name;
  double rec_scale;
  RecRateTable *rtab;
  
//  region = analysis_get_region(an);

//  sf_dba = seqfeat_dba_create(an->dbc);

  rec_tab_name = config_get_str(config, "RECOMB_RATE_TABLE");
  rec_scale = config_get_double(config, "RECOMB_RATE_SCALE");

  /* get recombination rates from database */
  sf = load_genetic_map_from_file(rec_tab_name, chr, &n_sf);
//  sf = seqfeat_dba_fetch_by_region(sf_dba, region, rec_tab_name, &n_sf);
  rtab = rectab_from_feats(sf, n_sf, chr->len, rec_scale);
  seqfeat_array_free(sf, n_sf);

  return rtab;
}





/**
 * Converts provided recombination rates to recombination distances
 * along chromosome and retrieves a list of rec-dist coords
 * representing conserved blocks. Conserved blocks are split when
 * recombination rates change in order to make integral approximation
 * of background selection sum more efficient.
 */
static GList *get_cons_rec_dists(Config *config, RecRateTable *rtab) {
  long i, j, n_sf, ttl_left, ttl_right, ttl_cons, last_end;
//  long len; 
//  SeqCoord *region;
//  SeqFeatureDBA *sf_dba;
  SeqFeature *sf;
  ConsBlock *cblk;
  GList *cons_list, *cur;
  char *cons_tab;
  
//  region = analysis_get_region(an);
//  len = region->end - region->start + 1;

//  sf_dba = seqfeat_dba_create(an->dbc);

  cons_tab = config_get_str(config, "CONS_TABLE");

  /* get conserved elements from database */
  sf = load_conserved_from_file(cons_tab, config_get_str(config, "CHROMOSOME_NAME"), &n_sf);
//  sf = seqfeat_dba_fetch_by_region(sf_dba, region, cons_tab, &n_sf);

//  /* order conserved elements */
//  qsort(sf, n_sf, sizeof(SeqFeature), seqfeat_cmp_nostrand);

  /* Create list of conserved block coords splitting elements when
   * they span a change in recombination rate. The coordinates are
   * initially set to physical positions, and then updated to be
   * recombination distances.
   */
  cons_list = NULL;
  ttl_cons = 0;


  last_end = -1;
  for(i = 0; i < n_sf; i++) {
    cblk = g_new(ConsBlock, 1);
    cblk->start = sf[i].c.start;

    /* sanity check: we don't want overlapping blocks */
    if(cblk->start <= last_end) {
      g_error("get_cons_rec_dists: conserved blocks should not overlap: "
	      "cur_start=%ld, last_end=%ld", cblk->start, last_end);
    }
    
    ttl_cons += sf[i].c.end - sf[i].c.start + 1;

    cblk->r = rectab_rate(rtab, sf[i].c.start);
    for(j = sf[i].c.start; j < sf[i].c.end; j++) {
      if(rectab_rate(rtab,j+1) != cblk->r) {
	/* change in rate, end previous block */
	cblk->end = j;
	last_end = cblk->end;

	cons_list = g_list_append(cons_list, cblk);
	
/* 	fprintf(stderr, "splitting cons segment %ld-%ld into " */
/* 		"%ld-%ld and %ld-%ld (rec-rate change from %g to %g)\n", */
/* 		cblk->start, sf[i].c.end, */
/* 		cblk->start, cblk->end, */
/* 		j+1, sf[i].c.end, cblk->r, rectab_rate(rtab,j+1)); */

	/* start new block */
	cblk = g_new(ConsBlock, 1);
	cblk->start = j+1;
	cblk->r = rectab_rate(rtab, j+1);
      }
    }
    /* end current block */
    cblk->end = sf[i].c.end;
    last_end = cblk->end;

    cons_list = g_list_append(cons_list, cblk);
  }

  seqfeat_array_free(sf, n_sf);
  
  /* convert cons block phys coords to rec-dist coords */
  ttl_left = 0;
  ttl_right = ttl_cons;

  cur = cons_list;
  while(cur != NULL) {
    cblk = cur->data;
    cblk->r_start = rectab_rpos(rtab, cblk->start);
    cblk->r_end   = rectab_rpos(rtab, cblk->end);

    /* keep track of number of cons sites to left and right of each
     * blk (not including the block itself)
     */
    cblk->left_ttl = ttl_left;
    ttl_left  += cblk->end - cblk->start + 1;

    ttl_right -= cblk->end - cblk->start + 1;
    cblk->right_ttl = ttl_right;

    cur = g_list_next(cur);
  }

  return cons_list;
}



/* 
 * Calculates background selection strength at each position along a
 * chromosome using provided parameters, recombination map and list of
 * conserved elements.
 * 
 */
void calc_bkgd_chr(Config *config, Chromosome *chr, RecRateTable *rtab, GList *cons_list,
		   BkgdParam *parm,  FILE *out_fh) {
  double b;
  long pos;/*, chr_len;*/
  GList *next_cons;
  int b_int, prev_b_int, b_len;
//  SeqCoord *region;
  BkgdInterp *bgi;

//  region = analysis_get_region(an);
//  chr_len = rtab->blk[rtab->n-1].end ;//region->end - region->start + 1;

  next_cons = cons_list;

  /* b is background selection strength */
  prev_b_int = b_int = -1;

  b_len = 0;

  /* Create interpolator to estimate B values at positions along chr */
  bgi = bkgd_interp_new(rtab, chr->len, cons_list, parm);

  pos = 1;
  while(pos <= chr->len) {
    b = bkgd_interp_eval(bgi, pos);
    
    if((pos % 1000000)==0) {
      fprintf(stderr, ".");
    }
    /* fprintf(stderr, "pos=%ld, b=%g\n", pos, b); */

    /* truncate to 3 digit integer */
    b_int = (int)floor(b* BKGD_SCALE + 0.5);

    /* only print out value if rounded value is different from prev one */
    if(prev_b_int != b_int) {
      if(prev_b_int >= 0) {
	fprintf(out_fh, "%d %d\n", prev_b_int, b_len);
	/* fprintf(stderr, "%ld %d %d\n", pos, prev_b_int, b_len); */
      }
      prev_b_int = b_int;
      b_len = 0;
    }

    pos++;
    b_len++;
  }

  /* print out final value */
  if(b_len > 0) {
    fprintf(out_fh, "%d %d\n", b_int, b_len);
    /* fprintf(stderr, "%ld %d %d\n", pos, prev_b_int, b_len); */
  }

  fprintf(stderr, "\n");

  bkgd_interp_free(bgi);

}






int main(int argc, char **argv) {
  Config *config;
  Chromosome *chr;
//  Analysis *analysis;
  GList *cons_list, *cur;
  char *out_dir, *out_path, *out_token;
//  SeqCoord *region;
  RecRateTable *rtab;
  FILE *out_fh;
  BkgdParam *parm;
  
  if(argc != 3) {
    fprintf(stderr, "usage: %s <config_file> <output_file_name>\n",
	    argv[0]);
    exit(2);
  }

  fprintf(stderr, "Reading config\n");

  config = config_read_file(argv[1], CONFIG_MISSING_KEY_ERROR);
//  analysis = analysis_new(config);
  chr = g_new(Chromosome, 1);
/* GUY: Changed the length to be a direct configuration input
 *  chr_features_file = config_get_str(config, "CHROMOSOME_FEATURES");
 *  get_chr_features(chr, argv[2], chr_features_file);
 */
  chr->len = atol(config_get_str(config, "CHROMOSOME_LENGTH"));
  chr->id = -1;
  chr->assembly = NULL;
  chr->name = g_new(char, 256);
  sscanf(config_get_str(config, "CHROMOSOME_NAME"), "%s", chr->name);

  out_dir = config_get_str(config, "OUTPUT_DIR");
  // GUY: changed filneame to be a commandline argument
  // out_token = config_get_str(config, "OUTPUT_TOKEN");
  out_token = argv[2];


  parm = bkgd_param_new(config);

//  while(analysis_next_region(analysis)) {
//    region = analysis_get_region(analysis);

//    fprintf(stderr, "\nanalysis region: %s\n", 
//	    analysis_get_region_name(analysis));

//    out_path = g_strconcat(out_dir, out_token, ".bkgd", NULL);
	out_path = g_strconcat(out_dir, out_token, ".bkgd", NULL);
    out_fh = fopen(out_path, "w");
    if(out_fh == NULL) {
      g_error("main: could not open output file '%s'", out_path); 
    }
    g_free(out_path);

	print_params(out_fh, chr, argv[1]);

    fprintf(stderr, "retrieving recombination rates\n");
    rtab = get_rectab(config, chr);
//	rtab = get_rectab(analysis);
    
    /* get recombination distances of conserved sites,
     * and convert rec rates to dists
     */
    fprintf(stderr, "calculating conserved rec dists\n");
    cons_list = get_cons_rec_dists(config, rtab);
//    cons_list = get_cons_rec_dists(analysis, rtab);

    fprintf(stderr, "  total recomb dist for %s: %gM\n", 
	    chr->name, rtab->chr_r_len);

    /* calculate strength of background selection at each position on chr */
    fprintf(stderr, "calculating b vals\n");
    calc_bkgd_chr(config, chr, rtab, cons_list, parm, out_fh);
    
    fprintf(stderr, "freeing recombination distances\n");
    rectab_free(rtab);
	    
    /* free conserved blocks */
    fprintf(stderr, "freeing conserved blocks\n");
    cur = cons_list;
    while(cur != NULL) {
      g_free(cur->data);
      cur = g_list_next(cur);
    }
    g_list_free(cons_list);

    fclose(out_fh);
//  }

  fprintf(stderr, "freeing parameters and interpolation tables\n");
  bkgd_param_free(parm);

  fprintf(stderr, "freeing analysis and config\n");
//  analysis_free(analysis);
  config_free(config);

  return 0;  
}

