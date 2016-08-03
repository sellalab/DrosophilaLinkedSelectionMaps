#include <math.h>
#include <glib.h>
#include <zlib.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "util.h"
#include "elog.h"


#include "bkgd_data.h"


#define MAX_LINE_TOK 100
#define BKGD_COUNT_FILE_TYPE_HCM 1
#define BKGD_COUNT_FILE_TYPE_HCGOM 2


/*
 * Checks if two doubles are "close" to each other.
 * Returns: -1 if d2 > d1, 
 *           1 if d1 > d2
 *           0 if d1 close to d2
 */
static int dbl_cmp(double d1, double d2) {
  const double delta = 1e-6;

  if(d1 + delta < d2) {
    return -1;
  }
  if(d1 - delta > d2) {
    return 1;
  }

  return 0;
}




/* bin comparison function, used to sort bkgd bins 
 * by B_ex, B_nex, cat
 */
static int cmp_bin(const void *v1, const void *v2) {
  const BkgdBin *b1 = v1;
  const BkgdBin *b2 = v2;
  int cmp;

  b1 = v1;
  b2 = v2;

  cmp = dbl_cmp(b1->B_ex, b2->B_ex);
  if(cmp != 0) {
    return cmp;
  }

  cmp = dbl_cmp(b1->B_nex, b2->B_nex);
  if(cmp != 0) {
    return cmp;
  }

  cmp = dbl_cmp(b1->cat, b2->cat);
  return cmp;
}



/* bin comparison function, used to sort bkgd bins 
 * by B_ex, cat (but not B_nex)
 */
static int cmp_bin_no_B_nex(const void *v1, const void *v2) {
  const BkgdBin *b1 = v1;
  const BkgdBin *b2 = v2;
  int cmp;

  b1 = v1;
  b2 = v2;

  cmp = dbl_cmp(b1->B_ex, b2->B_ex);
  if(cmp != 0) {
    return cmp;
  }

  cmp = dbl_cmp(b1->cat, b2->cat);
  return cmp;
}




/**
 * Copies values from one bin into another
 */
static void bin_copy(BkgdBin *to, const BkgdBin *from) {
  to->B_ex = from->B_ex;
  to->lB_ex = from->lB_ex;
  to->B_nex = from->B_nex;
  to->lB_nex = from->lB_nex;
  to->cat = from->cat;

  to->h = from->h;
  to->c = from->c;
  to->g = from->g;
  to->o = from->o;
  to->m = from->m;
  to->hc = from->hc;
  to->hg = from->hg;
  to->cg = from->cg;
  to->ho = from->ho;
  to->co = from->co;
  to->hcg = from->hcg;

  to->h_i = from->h_i;
  to->c_i = from->c_i;
  to->g_i = from->g_i;
  to->o_i = from->o_i;
  to->m_i = from->m_i;
  to->hc_i = from->hc_i;
  to->hg_i = from->hg_i;
  to->cg_i = from->cg_i;
  to->ho_i = from->ho_i;
  to->co_i = from->co_i;
  to->hcg_i = from->hcg_i;

  to->h_v = from->h_v;
  to->c_v = from->c_v;
  to->g_v = from->g_v;
  to->o_v = from->o_v;
  to->m_v = from->m_v;
  to->hc_v = from->hc_v;
  to->hg_v = from->hg_v;
  to->cg_v = from->cg_v;
  to->ho_v = from->ho_v;
  to->co_v = from->co_v;
  to->hcg_v = from->hcg_v;

  to->cons = from->cons;
}


static void bin_add_counts(BkgdBin *to, const BkgdBin *from) {
  to->h += from->h;
  to->c += from->c;
  to->g += from->g;
  to->o += from->o;
  to->m += from->m;
  to->hc += from->hc;
  to->hg += from->hg;
  to->cg += from->cg;
  to->ho += from->ho;
  to->co += from->co;
  to->hcg += from->hcg;

  to->h_i += from->h_i;
  to->c_i += from->c_i;
  to->g_i += from->g_i;
  to->o_i += from->o_i;
  to->m_i += from->m_i;
  to->hc_i += from->hc_i;
  to->hg_i += from->hg_i;
  to->cg_i += from->cg_i;
  to->ho_i += from->ho_i;
  to->co_i += from->co_i;
  to->hcg_i += from->hcg_i;

  to->h_v += from->h_v;
  to->c_v += from->c_v;
  to->g_v += from->g_v;
  to->o_v += from->o_v;
  to->m_v += from->m_v;
  to->hc_v += from->hc_v;
  to->hg_v += from->hg_v;
  to->cg_v += from->cg_v;
  to->ho_v += from->ho_v;
  to->co_v += from->co_v;
  to->hcg_v += from->hcg_v;

  to->cons += from->cons;      
}


static void bin_set_counts_zero(BkgdBin *bin) {
  bin->h = 0;
  bin->c = 0;
  bin->g = 0;
  bin->o = 0;
  bin->m = 0;
  bin->hc = 0;
  bin->hg = 0;
  bin->cg = 0;
  bin->hcg = 0;
  bin->ho = 0;
  bin->co = 0;

  bin->h_i = 0;
  bin->c_i = 0;
  bin->g_i = 0;
  bin->o_i = 0;
  bin->m_i = 0;
  bin->hc_i = 0;
  bin->hg_i = 0;
  bin->cg_i = 0;
  bin->hcg_i = 0;
  bin->ho_i = 0;
  bin->co_i = 0;

  bin->h_v = 0;
  bin->c_v = 0;
  bin->g_v = 0;
  bin->o_v = 0;
  bin->m_v = 0;
  bin->hc_v = 0;
  bin->hg_v = 0;
  bin->cg_v = 0;
  bin->hcg_v = 0;
  bin->ho_v = 0;
  bin->co_v = 0;

  bin->cons = 0;
}



/**
 * Returns the number of sites in the provided data bin
 */
long bkgd_data_n_bin_sites(const BkgdBin *bin) {
  long n_sites;
  
  /* first pass: count number of sites */
  n_sites = 0;
  n_sites += bin->h;
  n_sites += bin->c;
  n_sites += bin->g;
  n_sites += bin->o;
  n_sites += bin->m;
  n_sites += bin->hc;
  n_sites += bin->hg;
  n_sites += bin->cg;
  n_sites += bin->hcg;
  n_sites += bin->ho;
  n_sites += bin->co;
  n_sites += bin->cons;
  
  return n_sites;  
}



/**
 * Returns the total number of sites in the provided data set
 */
long bkgd_data_n_sites(BkgdEvoMdlData *data) {
  long n_sites;
  int i;

  n_sites = 0;
  for(i = 0; i < data->n_bin; i++) {
    n_sites += bkgd_data_n_bin_sites(&data->bin[i]);
  }
  return n_sites;
}





/**
 * Calculates the mean category value (weighted by number sites with
 * each cat val) for the provided data.
 */
double bkgd_data_mean_cat_val(BkgdEvoMdlData *data) {
  int i;
  long n;
  double cat_ttl, ttl;

  ttl = 0;
  cat_ttl = 0.0;
  for(i = 0; i < data->n_bin; i++) {
    n = bkgd_data_n_bin_sites(&data->bin[i]);
    cat_ttl += data->bin[i].cat * (double)n;
    ttl += n;
  }

  return cat_ttl / (double)ttl;
}



/**
 * Merges data bins into specified number of bins that have
 * approximately same numbers of data points. The actual number of
 * bins created may be less than n_bin depending how sites are
 * distributed.
 */
void bkgd_data_combine_bin(BkgdEvoMdlData *data, size_t n_bin) {
  double b_ex_sum, b_nex_sum, cat_sum;
  BkgdBin *new_bin, *old_bin;
  long i, j, n_sites, sites_per_bin, sites_in_bin, n;
  size_t n_old_bin;

  old_bin = data->bin;
  n_old_bin = data->n_bin;

  /* count number of sites */
  n_sites = bkgd_data_n_sites(data);

  if(n_sites <= 0) {
    g_error("bkgd_evo_mdl_combine_bin: bins do not contain any sites");
  }

  sites_per_bin = n_sites / n_bin;

  /* second pass: allocate new bins and copy counts from old bins */
  fprintf(stderr, "creating %zd new bins with %ld sites each\n", n_bin,
                  sites_per_bin);

  new_bin = g_new(BkgdBin, n_bin);
  for(i = 0; i < n_bin; i++) {
    bin_set_counts_zero(&new_bin[i]);
  }

  j = 0;
  sites_in_bin = 0;
  b_ex_sum = 0.0;
  b_nex_sum = 0.0;
  cat_sum = 0.0;

  for(i = 0; i < n_old_bin; i++) {
    /* start a new bin if we exceed the number of sites per bin */

    if((sites_in_bin >= sites_per_bin) && (j < n_bin-1)) {
      /* update B values and category value for bin */
      new_bin[j].B_ex = b_ex_sum / (double)sites_in_bin;
      new_bin[j].B_nex = b_nex_sum / (double)sites_in_bin;
      new_bin[j].lB_ex = log(new_bin[j].B_ex);
      new_bin[j].lB_nex = log(new_bin[j].B_nex);
      new_bin[j].cat = cat_sum / (double)sites_in_bin;

      fprintf(stderr, "bin %ld has %ld sites\n", j, sites_in_bin);

      fprintf(stderr, "B_ex=%g, B_nex=%g, cat=%g,"
	      "HG=%ld, CG=%ld, HO=%ld, CO=%ld\n", new_bin[j].B_ex,
	      new_bin[j].B_nex, new_bin[j].cat, new_bin[j].hg,
	      new_bin[j].cg, new_bin[j].ho, new_bin[j].co);

      /* move on to next bin */
      b_ex_sum = 0.0;
      b_nex_sum = 0.0;
      cat_sum = 0.0;
      sites_in_bin = 0;
      j++;
    }
      
    /* add counts to current bin */
    bin_add_counts(&new_bin[j], &old_bin[i]);

    n = bkgd_data_n_bin_sites(&old_bin[i]);

    b_ex_sum += old_bin[i].B_ex * (double)n;
    b_nex_sum += old_bin[i].B_nex * (double)n;
    cat_sum += old_bin[i].cat * (double)n;

    sites_in_bin += n;
  }

  if(sites_in_bin > 0) {
    /* update final bin */
    new_bin[j].B_ex = b_ex_sum / (double)sites_in_bin;
    new_bin[j].B_nex = b_nex_sum / (double)sites_in_bin;
    new_bin[j].cat = cat_sum / (double)sites_in_bin;
    new_bin[j].lB_ex = log(new_bin[j].B_ex);
    new_bin[j].lB_nex = log(new_bin[j].B_nex);

    fprintf(stderr, "bin %ld has %ld sites\n", j, sites_in_bin);
    n_bin = j+1;
  } else {
    n_bin = j;
  }

  data->n_bin = n_bin;
  data->bin = new_bin;

  g_free(old_bin);
}





/**
 * Discards information about B_nex values and collapses bins into their
 * B_ex values. B_nex is set to 1.0.
 */
void bkgd_data_collapse_nex_bin(BkgdEvoMdlData *data) {
  double cur_b_ex, cur_cat;
  BkgdBin *new_bin, *old_bin;
  long i, j;
  size_t n_old_bin, n_new_bin;

  old_bin = data->bin;
  n_old_bin = data->n_bin;

  /* make sure old bins are sorted */
  fprintf(stderr, "sorting old bins\n");
  qsort(old_bin, n_old_bin, sizeof(BkgdBin), &cmp_bin_no_B_nex);
 
  cur_b_ex = cur_cat = -1.0;

  n_new_bin = 0;

  /* first pass: count number of new bins that will be required */
  fprintf(stderr, "counting number of bins required\n");
  for(i = 0; i < n_old_bin; i++) {

    if((dbl_cmp(old_bin[i].B_ex, cur_b_ex) != 0) ||
       (dbl_cmp(old_bin[i].cat, cur_cat) != 0)) {
      /* last bin did not have same B_ex and category, will need new bin */
      n_new_bin += 1;
      cur_b_ex = old_bin[i].B_ex;
      cur_cat = old_bin[i].cat;
    }
  }

  /* second pass: allocate new bins and copy counts from old bins */
  fprintf(stderr, "creating %zd new bins\n", n_new_bin);
  new_bin = g_new(BkgdBin, n_new_bin);
  j = -1;

  cur_b_ex = cur_cat = -1.0;
  for(i = 0; i < n_old_bin; i++) {
    if((dbl_cmp(old_bin[i].B_ex, cur_b_ex) != 0) ||
       (dbl_cmp(old_bin[i].cat, cur_cat) != 0)) {
	/* different B_ex or cat than cur bin, start new bin */
	j++;

	if(j >= n_new_bin) {
	  g_error("bkgd_evo_mdl_combine_cat_bin: expected %zd new bins, "
		  "but there are at least %ld",
		  n_new_bin, j+1);
	}

	bin_copy(&new_bin[j], &old_bin[i]);
	new_bin[j].B_nex = 1.0;
	new_bin[j].lB_nex = 0.0;

	cur_b_ex = new_bin[j].B_ex;
	cur_cat = new_bin[j].cat;
    } else {
      /* add counts to current bin */
      bin_add_counts(&new_bin[j], &old_bin[i]);
    }
  }

  /* sanity check */
  if(j != n_new_bin-1) {
    g_error("%s:%d: expected %zd new bins, but there were only %ld", 
	    __FILE__, __LINE__, n_new_bin, j+1);
  }

  fprintf(stderr, "created %zd bins from %zd old bins\n", n_new_bin, 
	  n_old_bin);

  data->n_bin = n_new_bin;
  data->bin = new_bin;

  g_free(old_bin);
}


/**
 * Collapses separate category bins, and discards cat bins that are outside
 * specified range of allowed category values. The old bins pointed to
 * by data->bin are freed and replaced with newly created bins.
 */
void bkgd_data_combine_cat_bin(BkgdEvoMdlData *data,
			       double min_cat, double max_cat) {
  double cur_b_ex, cur_b_nex;
  BkgdBin *new_bin, *old_bin;
  long i, j;
  size_t n_old_bin, n_new_bin;

  old_bin = data->bin;
  n_old_bin = data->n_bin;

  /* make sure old bins are sorted */
  fprintf(stderr, "sorting old bins\n");
  qsort(old_bin, n_old_bin, sizeof(BkgdBin), &cmp_bin);
 
  cur_b_ex = cur_b_nex = -1.0;

  n_new_bin = 0;

  /* first pass: count number of new bins that will be required */
  fprintf(stderr, "counting number of bins required\n");
  for(i = 0; i < n_old_bin; i++) {
    if((old_bin[i].cat < min_cat) || (old_bin[i].cat > max_cat)) {
      /* skip this bin, outside of specified category range */
      continue;
    }

    if((dbl_cmp(old_bin[i].B_ex, cur_b_ex) != 0) ||
       (dbl_cmp(old_bin[i].B_nex, cur_b_nex) != 0)) {
      /* last bin did not have same B_ex and B_nex, will need new bin */
      n_new_bin += 1;
      cur_b_ex = old_bin[i].B_ex;
      cur_b_nex = old_bin[i].B_nex;
    }
  }

  /* second pass: allocate new bins and copy counts from old bins */
  fprintf(stderr, "creating %zd new bins\n", n_new_bin);
  new_bin = g_new(BkgdBin, n_new_bin);
  j = -1;

  cur_b_ex = cur_b_nex = -1.0;
  for(i = 0; i < n_old_bin; i++) {
    if((old_bin[i].cat < min_cat) || (old_bin[i].cat > max_cat)) {
      /* category value is not in range we are using, discard bin */
      continue;
    }

    if((dbl_cmp(old_bin[i].B_ex, cur_b_ex) != 0) ||
       (dbl_cmp(old_bin[i].B_nex, cur_b_nex) != 0)) {
	/* different B_ex or B_nex than cur bin, start new bin */
	j++;

	if(j >= n_new_bin) {
	  g_error("bkgd_evo_mdl_combine_cat_bin: expected %zd new bins, "
		  "but there are at least %ld",
		  n_new_bin, j+1);
	}

	bin_copy(&new_bin[j], &old_bin[i]);
	new_bin[j].cat = 0.0;

	cur_b_ex = new_bin[j].B_ex;
	cur_b_nex = new_bin[j].B_nex;
    } else {
      /* add counts to current bin */
      bin_add_counts(&new_bin[j], &old_bin[i]);
    }
  }

  /* sanity check */
  if(j != n_new_bin-1) {
    g_error("bkgd_evo_mdl_combine_cat_bin: expected %zd new bins, "
	    "but there were only %ld", n_new_bin, j+1);
  }

  fprintf(stderr, "created %zd bins from %zd old bins\n", n_new_bin, 
	  n_old_bin);

  data->n_bin = n_new_bin;
  data->bin = new_bin;

  g_free(old_bin);
}



/**
 * Checks that the header line contains the appropriate number of
 * columns
 */
static void check_header_line(char *hdr_line, int file_type) {
  char **tok;
  int n;
  
  tok = g_strsplit(hdr_line, "\t", MAX_LINE_TOK);
  n = 0;
  while((tok[n] != NULL) && (n < MAX_LINE_TOK)) { 
    n++; 
  }
 
  if(file_type == BKGD_COUNT_FILE_TYPE_HCM) {
    if(n != 9) {
      g_error("get_file_type: expected HCM header line to "
	      "have 9 tokens");
    }
  }
  else if(file_type == BKGD_COUNT_FILE_TYPE_HCGOM) {
    if(n != 16) {
      g_error("get_file_type: expected HCGOM header line to "
	      "have 16 tokens");
    }
  }
  else {
    g_error("get_file_type: unknown file type");
  }
}



/**
 * Retrieves the filetype from the config, and
 * checks that header line conforms to format
 */
static int get_file_type(Config *config) {
  char *type_str;

  type_str = config_get_str(config, "BKGD_COUNT_FILE_TYPE");

  if(strcmp(type_str, "HCM")==0) {
    return BKGD_COUNT_FILE_TYPE_HCM;
  }
  else if(strcmp(type_str, "HCGOM")==0) {
    return BKGD_COUNT_FILE_TYPE_HCGOM;
  }

  g_error("get_file_type: unknown file type '%s'", type_str);

  return 0;
}




static BkgdEvoMdlData *read_data_file(const char *file, const double cat_scale,
				      const int file_type) {
  int i;
  char line[MAX_LINE];
  char *str_next, *str_cur;
  long n_lines, val;
  gzFile *gzf;
  BkgdEvoMdlData *data;

  fprintf(stderr, "Reading data\n");

  gzf = gzopen(file, "rb");

  /* count number of lines in file */
  n_lines = 0;
  while((gzgets(gzf, line, MAX_LINE)) != NULL) {
    n_lines++;
  }

  /* move back to beginning of file */
  if(gzrewind(gzf) != 0) {
    g_error("bkgd_evo_mdl_read_data: gzrewind of file '%s' failed", file);
  }

  data = g_new(BkgdEvoMdlData,1);

  if(n_lines > 0) {
    data->n_bin = n_lines - 1;
    data->bin = g_new(BkgdBin, data->n_bin);
  } else {
    data->n_bin = 0;
    data->bin = NULL;
    return data;
  }

  /* get header line */
  if(gzgets(gzf, line, MAX_LINE) == NULL) {
    g_error("bkgd_evo_mdl_read_data: expected %ld lines in file '%s', "
	    "but only got 0",  n_lines, file);
  }

  check_header_line(line, file_type);

    
  for(i = 0; i < data->n_bin; i++) {
    /* read line from file */
    if(gzgets(gzf, line, MAX_LINE) == NULL) {
      g_error("bkgd_evo_mdl_read_data: expected %ld lines in file '%s', "
	      "but only got %d", n_lines, file, i+1);
    }
    
    /* read through tokens on line (all tokens are integers) */

    /* B_ex */
    str_cur = line;
    val = strtol(str_cur, &str_next, 10);
    if(val <= 0) {
      data->bin[i].B_ex = 0.5 / (double)BKGD_BIN_SCALE;
    } else {
      data->bin[i].B_ex = (double)val / (double)BKGD_BIN_SCALE;
    }
    data->bin[i].lB_ex = log(data->bin[i].B_ex);
    
    /* B_nex */
    str_cur = str_next;
    val = strtol(str_cur, &str_next, 10);
    if(val <= 0) {
      data->bin[i].B_nex = 0.5 / (double)BKGD_BIN_SCALE;
    } else {
      data->bin[i].B_nex = (double)val / (double)BKGD_BIN_SCALE;
    }

    data->bin[i].lB_nex = log(data->bin[i].B_nex);

    /* category value */
    str_cur = str_next;
    val = strtol(str_cur, &str_next, 10);
    data->bin[i].cat = (double)val / cat_scale;

    /* N (not used) */
    str_cur = str_next;
    val = strtol(str_cur, &str_next, 10);

    /* human-branch diffs */
    str_cur = str_next;
    data->bin[i].h = strtol(str_cur, &str_next, 10);

    /* chimp-branch diffs */
    str_cur = str_next;
    data->bin[i].c = strtol(str_cur, &str_next, 10);

    if(file_type == BKGD_COUNT_FILE_TYPE_HCM) {
      /* macaque-branch diffs */
      str_cur = str_next;
      data->bin[i].m = strtol(str_cur, &str_next, 10);

      /* all diff (not used) */
      str_cur = str_next;
      val = strtol(str_cur, &str_next, 10);

      /* cons */
      str_cur = str_next;
      data->bin[i].cons = strtol(str_cur, &str_next, 10);

      data->bin[i].g = 0;
      data->bin[i].o = 0;
      data->bin[i].hc = 0;
      data->bin[i].hg = 0;
      data->bin[i].cg = 0;
      data->bin[i].ho = 0;
      data->bin[i].co = 0;
      data->bin[i].hcg = 0;
    } 
    else if(file_type == BKGD_COUNT_FILE_TYPE_HCGOM) {
      /* gorilla-branch diffs */
      str_cur = str_next;
      data->bin[i].g = strtol(str_cur, &str_next, 10);

      /* orang-branch diffs */
      str_cur = str_next;
      data->bin[i].o = strtol(str_cur, &str_next, 10);

      /* macaque-branch diffs */
      str_cur = str_next;
      data->bin[i].m = strtol(str_cur, &str_next, 10);

      /* HC diffs */
      str_cur = str_next;
      data->bin[i].hc = strtol(str_cur, &str_next, 10);

      /* HG diffs */
      str_cur = str_next;
      data->bin[i].hg = strtol(str_cur, &str_next, 10);

      /* CG diffs */
      str_cur = str_next;
      data->bin[i].cg = strtol(str_cur, &str_next, 10);
      
      /* HO diffs */
      str_cur = str_next;
      data->bin[i].ho = strtol(str_cur, &str_next, 10);

      /* CO diffs */
      str_cur = str_next;
      data->bin[i].co = strtol(str_cur, &str_next, 10);

      /* HCG diffs */
      str_cur = str_next;
      data->bin[i].hcg = strtol(str_cur, &str_next, 10);

      /* conserved count */
      str_cur = str_next;
      data->bin[i].cons = strtol(str_cur, &str_next, 10);
    }
    else {
      g_error("bkgd_evo_mdl_read_data: unknown file type");
    }
  }

  gzclose(gzf);

  return data;
}

/**
 * Creates background bins, reading from a file that groups by
 * exonic/non-exonic B-values as well as an additional category such
 * as GC content.
 */
BkgdEvoMdlData *bkgd_data_read_data(Config *config) {
  char *file;
  double cat_scale;
  int file_type;

  file = config_get_str(config, "BKGD_COUNT_FILE");

  cat_scale = config_get_double(config, "CAT_SCALE");

  file_type = get_file_type(config);

  fprintf(stderr, "reading data from file '%s'\n", file);

  return read_data_file(file, cat_scale, file_type);
}



/**
 * Reads from separate files containing counts of transitions
 * and transversions and unifies into single data structure
 */
BkgdEvoMdlData *bkgd_data_read_i_v_data(Config *config) {
  char *i_file, *v_file;
  double cat_scale, cur_b_ex, cur_b_nex, cur_cat;
  BkgdEvoMdlData *i_data, *v_data, *combined_data;
  BkgdBin *new_bin;
  long n_new_bin, n_combined_bin, i, j;
  int file_type;
  
  /* get transition file */
  i_file = config_get_str(config, "BKGD_I_COUNT_FILE");

  /* get transversion file */
  v_file = config_get_str(config, "BKGD_V_COUNT_FILE");

  cat_scale = config_get_double(config, "CAT_SCALE");

  file_type = get_file_type(config);
  
  i_data = read_data_file(i_file, cat_scale, file_type);
  v_data = read_data_file(v_file, cat_scale, file_type);
  
  /* place all bins into single array */
  n_new_bin = i_data->n_bin + v_data->n_bin;
  new_bin = g_new(BkgdBin, n_new_bin);

  fprintf(stderr, "combining i and v counts\n");
  
  /* add transitions counts to array */
  j = 0;
  for(i = 0; i < i_data->n_bin; i++) {
    new_bin[j].B_ex = i_data->bin[i].B_ex;
    new_bin[j].lB_ex = i_data->bin[i].lB_ex;
    new_bin[j].B_nex = i_data->bin[i].B_nex;
    new_bin[j].lB_nex = i_data->bin[i].lB_nex;
    new_bin[j].cat = i_data->bin[i].cat;

    new_bin[j].h_i = i_data->bin[i].h;
    new_bin[j].c_i = i_data->bin[i].c;
    new_bin[j].g_i = i_data->bin[i].g;
    new_bin[j].o_i = i_data->bin[i].o;
    new_bin[j].m_i = i_data->bin[i].m;
    new_bin[j].hc_i = i_data->bin[i].hc;
    new_bin[j].hg_i = i_data->bin[i].hg;
    new_bin[j].cg_i = i_data->bin[i].cg;
    new_bin[j].hcg_i = i_data->bin[i].hcg;
    new_bin[j].ho_i = i_data->bin[i].ho;
    new_bin[j].co_i = i_data->bin[i].co;

    new_bin[j].h_v = 0;
    new_bin[j].c_v = 0;
    new_bin[j].g_v = 0;
    new_bin[j].o_v = 0;
    new_bin[j].m_v = 0;
    new_bin[j].hc_v = 0;
    new_bin[j].hg_v = 0;
    new_bin[j].cg_v = 0;
    new_bin[j].hcg_v = 0;
    new_bin[j].ho_v = 0;
    new_bin[j].co_v = 0;

    new_bin[j].h = i_data->bin[i].h;
    new_bin[j].c = i_data->bin[i].c;
    new_bin[j].g = i_data->bin[i].g;
    new_bin[j].o = i_data->bin[i].o;
    new_bin[j].m = i_data->bin[i].m;
    new_bin[j].hc = i_data->bin[i].hc;
    new_bin[j].hg = i_data->bin[i].hg;
    new_bin[j].cg = i_data->bin[i].cg;
    new_bin[j].hcg = i_data->bin[i].hcg;
    new_bin[j].ho = i_data->bin[i].ho;
    new_bin[j].co = i_data->bin[i].co;

    new_bin[j].cons = i_data->bin[i].cons;

    j++;
  }

  g_free(i_data->bin);
  g_free(i_data);


  /* add transversion counts to array */
  for(i = 0; i < v_data->n_bin; i++) {
    new_bin[j].B_ex = v_data->bin[i].B_ex;
    new_bin[j].lB_ex = v_data->bin[i].lB_ex;
    new_bin[j].B_nex = v_data->bin[i].B_nex;
    new_bin[j].lB_nex = v_data->bin[i].lB_nex;
    new_bin[j].cat = v_data->bin[i].cat;

    new_bin[j].h_v = v_data->bin[i].h;
    new_bin[j].c_v = v_data->bin[i].c;
    new_bin[j].g_v = v_data->bin[i].g;
    new_bin[j].o_v = v_data->bin[i].o;
    new_bin[j].m_v = v_data->bin[i].m;
    new_bin[j].hc_v = v_data->bin[i].hc;
    new_bin[j].hg_v = v_data->bin[i].hg;
    new_bin[j].cg_v = v_data->bin[i].cg;
    new_bin[j].hcg_v = v_data->bin[i].hcg;
    new_bin[j].ho_v = v_data->bin[i].ho;
    new_bin[j].co_v = v_data->bin[i].co;

    new_bin[j].h_i = 0;
    new_bin[j].c_i = 0;
    new_bin[j].g_i = 0;
    new_bin[j].o_i = 0;
    new_bin[j].m_i = 0;
    new_bin[j].hc_i = 0;
    new_bin[j].hg_i = 0;
    new_bin[j].cg_i = 0;
    new_bin[j].hcg_i = 0;
    new_bin[j].ho_i = 0;
    new_bin[j].co_i = 0;

    new_bin[j].h += v_data->bin[i].h;
    new_bin[j].c += v_data->bin[i].c;
    new_bin[j].g += v_data->bin[i].g;
    new_bin[j].o += v_data->bin[i].o;
    new_bin[j].m += v_data->bin[i].m;
    new_bin[j].hc += v_data->bin[i].hc;
    new_bin[j].hg += v_data->bin[i].hg;
    new_bin[j].cg += v_data->bin[i].cg;
    new_bin[j].hcg += v_data->bin[i].hcg;
    new_bin[j].ho += v_data->bin[i].ho;
    new_bin[j].co += v_data->bin[i].co;

    new_bin[j].cons = v_data->bin[i].cons;

    j++;
  }

  g_free(v_data->bin);
  g_free(v_data);

  if(j != n_new_bin) {
     g_error("expected number of new bins to be %ld but got %ld", j, n_new_bin);
  }

  fprintf(stderr, "merging i and v bins\n");
  /* sort bins by their B_ex, B_nex and cat values */
  qsort(new_bin, n_new_bin, sizeof(BkgdBin), &cmp_bin);

  n_combined_bin = 0;
  cur_b_ex = cur_b_nex = cur_cat = -1.0;

  /* first pass: count number of new bins that will be required */
  fprintf(stderr, "counting number of bins required\n");
  for(i = 0; i < n_new_bin; i++) {
    if((dbl_cmp(new_bin[i].B_ex, cur_b_ex) != 0) ||
       (dbl_cmp(new_bin[i].B_nex, cur_b_nex) != 0) ||
       (dbl_cmp(new_bin[i].cat, cur_cat) != 0)) {

      /* last bin did not have same <B_ex, B_nex, cat>, will need new bin */
      n_combined_bin++;
      cur_b_ex = new_bin[i].B_ex;
      cur_b_nex = new_bin[i].B_nex;
      cur_cat = new_bin[i].cat;
    }
  }

  /* second pass: combine bins with same B_ex, B_nex and cat values */
  combined_data = g_new(BkgdEvoMdlData, 1);
  combined_data->n_bin = n_combined_bin;
  combined_data->bin = g_new(BkgdBin, n_combined_bin);

  j = -1;
  cur_b_ex = cur_b_nex = cur_cat = -1.0;

  for(i = 0; i < n_new_bin; i++) {
    if((dbl_cmp(new_bin[i].B_ex, cur_b_ex) != 0) ||
       (dbl_cmp(new_bin[i].B_nex, cur_b_nex) != 0) ||
       (dbl_cmp(new_bin[i].cat, cur_cat) != 0)) {

	/* different B_ex, B_nex, or cat than cur bin, start new bin */
	j++;

	if(j >= n_combined_bin) {
	  g_error("bkgd_evo_mdl_read_i_v_data: expected %ld new bins, "
		  "but there are at least %ld",
		  n_combined_bin, j+1);
	}

	bin_copy(&combined_data->bin[j], &new_bin[i]);

	cur_b_ex = combined_data->bin[j].B_ex;
	cur_b_nex = combined_data->bin[j].B_nex;
	cur_cat = combined_data->bin[j].cat;
    } else {

      /* add counts to current bin */
      combined_data->bin[j].h += new_bin[i].h;
      combined_data->bin[j].c += new_bin[i].c;
      combined_data->bin[j].g += new_bin[i].g;
      combined_data->bin[j].o += new_bin[i].o;
      combined_data->bin[j].m += new_bin[i].m;
      combined_data->bin[j].hc += new_bin[i].hc;
      combined_data->bin[j].hg += new_bin[i].hg;
      combined_data->bin[j].cg += new_bin[i].cg;
      combined_data->bin[j].ho += new_bin[i].ho;
      combined_data->bin[j].co += new_bin[i].co;
      combined_data->bin[j].hcg += new_bin[i].hcg;

      combined_data->bin[j].h_i += new_bin[i].h_i;
      combined_data->bin[j].c_i += new_bin[i].c_i;
      combined_data->bin[j].g_i += new_bin[i].g_i;
      combined_data->bin[j].o_i += new_bin[i].o_i;
      combined_data->bin[j].m_i += new_bin[i].m_i;
      combined_data->bin[j].hc_i += new_bin[i].hc_i;
      combined_data->bin[j].hg_i += new_bin[i].hg_i;
      combined_data->bin[j].cg_i += new_bin[i].cg_i;
      combined_data->bin[j].ho_i += new_bin[i].ho_i;
      combined_data->bin[j].co_i += new_bin[i].co_i;
      combined_data->bin[j].hcg_i += new_bin[i].hcg_i;

      combined_data->bin[j].h_v += new_bin[i].h_v;
      combined_data->bin[j].c_v += new_bin[i].c_v;
      combined_data->bin[j].g_v += new_bin[i].g_v;
      combined_data->bin[j].o_v += new_bin[i].o_v;
      combined_data->bin[j].m_v += new_bin[i].m_v;
      combined_data->bin[j].hc_v += new_bin[i].hc_v;
      combined_data->bin[j].hg_v += new_bin[i].hg_v;
      combined_data->bin[j].cg_v += new_bin[i].cg_v;
      combined_data->bin[j].ho_v += new_bin[i].ho_v;
      combined_data->bin[j].co_v += new_bin[i].co_v;
      combined_data->bin[j].hcg_v += new_bin[i].hcg_v;

      if(combined_data->bin[j].cons != new_bin[i].cons) {

	fprintf(stderr, "new: <B_ex=%g, B_nex=%g, cat=%g> "
		"old: <B_ex=%g, B_nex=%g, cat=%g>\n",
		cur_b_ex, cur_b_nex, cur_cat, 
		new_bin[i].B_ex, new_bin[i].B_nex, 
		new_bin[i].cat);

	g_error("bkgd_data_read_i_v_data: expected I and V bins "
		"to have same number of conserved sites");
      }
    }
  }

  /* sanity check */
  if(j != n_combined_bin-1) {
    g_error("bkgd_data_read_i_v_data: expected %ld combined bins, "
	    "but there were only %ld", n_combined_bin, j+1);
  }

  g_free(new_bin);

  return combined_data;
}





/**
 * Writes out column count totals from provided data structure,
 * useful for debuggging
 */
void bkgd_data_write_site_counts(FILE *fh, BkgdEvoMdlData *data) {
  long i, n_h, n_c, n_g, n_o, n_m, n_hc, n_hg, n_cg, n_hcg, n_ho, n_co,
    n_h_i, n_c_i, n_g_i, n_o_i, n_m_i, n_hc_i, n_hg_i, n_cg_i, n_hcg_i, 
    n_ho_i, n_co_i, n_h_v, n_c_v, n_g_v, n_o_v, n_m_v, n_hc_v, n_hg_v, 
    n_cg_v, n_hcg_v, n_ho_v, n_co_v, n_cons, n;

  n_h = n_c = n_g = n_o = n_m = n_hc = n_hg = n_cg = n_hcg = 
    n_ho = n_co = n_cons = n = 0;

  /* transition totals */
  n_h_i = n_c_i = n_g_i = n_o_i = n_m_i = n_hc_i = n_hg_i = n_cg_i = 
    n_hcg_i = n_ho_i = n_co_i = 0;

  /* transversion totals */
  n_h_v = n_c_v = n_g_v = n_o_v = n_m_v = n_hc_v = n_hg_v = n_cg_v = 
    n_hcg_v = n_ho_v = n_co_v = 0;

  for(i = 0; i < data->n_bin; i++) {
    n_h += data->bin[i].h;
    n_c += data->bin[i].c;
    n_g += data->bin[i].g;
    n_o += data->bin[i].o;
    n_m += data->bin[i].m;
    n_hc += data->bin[i].hc;
    n_hg += data->bin[i].hg;
    n_cg += data->bin[i].cg;
    n_hcg += data->bin[i].hcg;
    n_ho += data->bin[i].ho;
    n_co += data->bin[i].co;

    n_h_i += data->bin[i].h_i;
    n_c_i += data->bin[i].c_i;
    n_g_i += data->bin[i].g_i;
    n_o_i += data->bin[i].o_i;
    n_m_i += data->bin[i].m_i;
    n_hc_i += data->bin[i].hc_i;
    n_hg_i += data->bin[i].hg_i;
    n_cg_i += data->bin[i].cg_i;
    n_hcg_i += data->bin[i].hcg_i;
    n_ho_i += data->bin[i].ho_i;
    n_co_i += data->bin[i].co_i;

    n_h_v += data->bin[i].h_v;
    n_c_v += data->bin[i].c_v;
    n_g_v += data->bin[i].g_v;
    n_o_v += data->bin[i].o_v;
    n_m_v += data->bin[i].m_v;
    n_hc_v += data->bin[i].hc_v;
    n_hg_v += data->bin[i].hg_v;
    n_cg_v += data->bin[i].cg_v;
    n_hcg_v += data->bin[i].hcg_v;
    n_ho_v += data->bin[i].ho_v;
    n_co_v += data->bin[i].co_v;

    n_cons += data->bin[i].cons;
  }

  n = n_h + n_c + n_g + n_o + n_m + n_hc + n_hg + n_cg + n_hcg + n_co 
    + n_ho + n_cons;

  fprintf(fh, "site counts:\n");
  fprintf(fh, "   H: %ld [i=%ld; v=%ld] (%% %.4f)\n", n_h, 
	  n_h_i, n_h_v, 100.0*(double)n_h/(double)n);
  fprintf(fh, "   C: %ld [i=%ld; v=%ld] (%% %.4f)\n", n_c, 
	  n_c_i, n_c_v, 100.0*(double)n_c/(double)n);
  fprintf(fh, "   G: %ld [i=%ld; v=%ld] (%% %.4f)\n", n_g, 
	  n_g_i, n_g_v, 100.0*(double)n_g/(double)n);
  fprintf(fh, "   O: %ld [i=%ld; v=%ld] (%% %.4f)\n", n_o, 
	  n_o_i, n_o_v, 100.0*(double)n_o/(double)n);
  fprintf(fh, "   M: %ld [i=%ld; v=%ld] (%% %.4f)\n", n_m, 
	  n_m_i, n_m_v, 100.0*(double)n_m/(double)n);
  fprintf(fh, "  HC: %ld [i=%ld; v=%ld] (%% %.4f)\n", n_hc, 
	  n_hc_i, n_hc_v, 100.0*(double)n_hc/(double)n);
  fprintf(fh, "  HG: %ld [i=%ld; v=%ld] (%% %.4f)\n", n_hg, 
	  n_hg_i, n_hg_v, 100.0*(double)n_hg/(double)n);
  fprintf(fh, "  CG: %ld [i=%ld; v=%ld] (%% %.4f)\n", n_cg, 
	  n_cg_i, n_cg_v, 100.0*(double)n_cg/(double)n);
  fprintf(fh, " HCG: %ld [i=%ld; v=%ld] (%% %.4f)\n", n_hcg, 
	  n_hcg_i, n_hcg_v, 100.0*(double)n_hcg/(double)n);
  fprintf(fh, "  HO: %ld [i=%ld; v=%ld] (%% %.4f)\n", n_ho, 
	  n_ho_i, n_ho_v, 100.0*(double)n_ho/(double)n);
  fprintf(fh, "  CO: %ld [i=%ld; v=%ld] (%% %.4f)\n", n_co, 
	  n_co_i, n_co_v, 100.0*(double)n_co/(double)n);
  fprintf(fh, "CONS: %ld (%% %.4f)\n", n_cons, 
	  100.0*(double)n_cons/(double)n);
  fprintf(fh, "TTL: %ld (%% %.4f)\n", n, 
	  100.0*(double)n/(double)n);
}




/*
 * Performs a jukes-cantor correction to convert divergence values
 * into approximate branch_length * mutation_rate. This function
 * modifies the category values of the data bins, and it only makes
 * sense to call this function if the cat values are interspecies
 * nucleotide divergence.
 */
void bkgd_data_jukes_cantor_correct_cats(BkgdEvoMdlData *data) {
  long i;

  for(i = 0; i < data->n_bin; i++) {
    /* fprintf(stderr, "old=%g ", data->bin[i].cat); */
    data->bin[i].cat = -0.75 * log(1.0 - (4.0/3.0)*data->bin[i].cat);
    /* fprintf(stderr, "new=%g\n", data->bin[i].cat); */
  }
}



/**
 * Sets category values to equal macaque divergence (proportion of
 * M sites)
 */

void bkgd_data_set_mdiv_cats(BkgdEvoMdlData *data) {
  long i, n;

  for(i = 0; i < data->n_bin; i++) {
    n = 0;
    n += data->bin[i].h;
    n += data->bin[i].c;
    n += data->bin[i].g;
    n += data->bin[i].o;
    n += data->bin[i].m;
    n += data->bin[i].hc;
    n += data->bin[i].hg;
    n += data->bin[i].cg;
    n += data->bin[i].hcg;
    n += data->bin[i].ho;
    n += data->bin[i].co;
    n += data->bin[i].cons;

    data->bin[i].cat = (n == 0) ? 0 : (double)data->bin[i].m / (double)n;
//	data->bin[i].cat = (n == 0) ? : (double)data->bin[i].m / (double)n;
  }
}
