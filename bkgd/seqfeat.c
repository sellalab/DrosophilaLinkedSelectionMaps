#include <stdio.h>
#include <stdlib.h>
#include <glib.h>
#include <string.h>

#include "seqcoord.h"
#include "seqfeat.h"
#include "nuc.h"

#include "util.h"



/**
 * Initializes an empty seqfeature 
 */
void seqfeat_init(SeqFeature *feat) {
  feat->id        = 0;
  feat->c.seqname = NULL;
  feat->c.chr     = NULL;
  feat->c.start   = 0;
  feat->c.end     = 0;
  feat->c.strand  = STRAND_NONE;
  feat->name    = NULL;
  feat->score   = 0.0;

  feat->attrib = NULL;

  feat->sub_feats = NULL;
  feat->n_sub_feat = 0;
}


/**
 * Allocates memory for and initializes a new SeqFeature structure.
 */
SeqFeature *seqfeat_new(const char *seqname, const long start, const long end,
			const short strand, const char *name, 
                        const float score) {
  SeqFeature *feat = g_new(SeqFeature, 1);
  feat->c.seqname = g_strdup(seqname);
  feat->c.start   = start;
  feat->c.end     = end;
  feat->c.strand  = strand;
  feat->name    = g_strdup(name);
  feat->score   = score;

  feat->attrib = NULL;

  feat->sub_feats = NULL;
  feat->n_sub_feat = 0;

  return feat;
}


/**
 * Copies the attributes of one seqfeature to another
 */
void seqfeat_copy_attrib(const SeqFeature *src, SeqFeature *dest) {
  SeqFeatAttrib *attr, *cur, *next;

  /* copy attributes */
  attr = src->attrib;
  dest->attrib = cur = NULL;
  while(attr != NULL) {
    next = g_new(SeqFeatAttrib, 1);
    next->next = NULL;
    if(cur == NULL) {
      cur = next;
      dest->attrib = cur;
    } else {
      cur->next = next;
      cur = next;
    }
    cur->name = g_strdup(attr->name);
    cur->value = g_strdup(attr->value);

    attr = attr->next;
  }
}


/**
 * Copies the features of one feature to another.
 */
void seqfeat_copy_subfeat(const SeqFeature *src, SeqFeature *dest) {
  long i;

  if(src->n_sub_feat == 0) {
    dest->n_sub_feat = 0;
    dest->sub_feats = NULL;
    return;
  }

  dest->sub_feats = g_new(SeqFeature, dest->n_sub_feat);
  for(i = 0; i < dest->n_sub_feat; i++) {
    seqfeat_copy(&src->sub_feats[i], &dest->sub_feats[i]);
  }
  dest->n_sub_feat = src->n_sub_feat;
}


/**
 * Copies contents of one seqfeature into another
 */
void seqfeat_copy(const SeqFeature *src, SeqFeature *dest) {
  /* copy name, id, score */
  dest->id        = src->id;
  dest->name = (src->name) ? g_strdup(src->name) : NULL;
  dest->score = src->score;

  /* copy coordinate */
  seq_coord_copy(&src->c, &dest->c);

  /* copy attributes and sub features */
  seqfeat_copy_attrib(src, dest);
  seqfeat_copy_subfeat(src, dest);

  return;
}




/* frees attributes associated with this seq feature */
void seqfeat_attrib_free(SeqFeature *sf) {
  SeqFeatAttrib *attr, *next_attr;

  attr = sf->attrib;

  while(attr != NULL) {
    g_free(attr->name);
    g_free(attr->value);
    next_attr = attr->next;
    g_free(attr);
    attr = next_attr;
  }

  return;
}



/**
 * Frees a singly-linked list of SeqFeature structures
 * and frees the list itself.
 */
void seqfeat_slist_free(GSList *segments) {
  GSList *cur;
  SeqFeature *feat;

  cur = segments;
  while(cur != NULL) {
    feat = (SeqFeature *)cur->data;
    g_free(feat->c.seqname);
    if(feat->name != NULL) {
      g_free(feat->name);
    }
    if(feat->c.seqname != NULL) {
      g_free(feat->c.seqname);
    }

    seqfeat_attrib_free(feat);

    if(feat->n_sub_feat != 0) {
      seqfeat_array_free(feat->sub_feats, feat->n_sub_feat);
    }

    g_free(feat);

    cur = g_slist_next(cur);
  }
 
  g_slist_free(segments);
}



/**
 * Frees an entire array of allocated SeqFeatures.
 */
void seqfeat_array_free(SeqFeature *feats, const long n) {
  long i;

  if(n == 0) {
    return;
  }

  for(i = 0; i < n; i++) {
    if(feats[i].name != NULL) {
      g_free(feats[i].name);
    }
    if(feats[i].c.seqname != NULL) {
      g_free(feats[i].c.seqname);
    }
    if(feats[i].n_sub_feat != 0) {
      seqfeat_array_free(feats[i].sub_feats, feats[i].n_sub_feat);
    }

    seqfeat_attrib_free(&feats[i]);
  }

  g_free(feats);
}



/**
 * Writes a list of SeqFeatures in BED format to a provided filehandle
 */
void seqfeat_slist_write_bed(FILE *fh, GSList *seqfeats) {
  SeqFeature *feat;
  GSList *cur;
  char *seqname;
  
  cur = seqfeats;
  while(cur != NULL) {
    feat = (SeqFeature *)cur->data;

    if(feat->c.chr != NULL && feat->c.chr->name != NULL) {
      seqname = feat->c.chr->name;
    } else if(feat->c.seqname != NULL) {
      seqname = feat->c.seqname;
    } else {
      seqname = "";
    }

    fprintf(fh, "%s\t%lu\t%lu\t%s\t%f\t%c\n", seqname, 
	    feat->c.start - 1, feat->c.end, 
	    (feat->name == NULL) ? "" : feat->name, 
	    feat->score,
	    strand_to_char(feat->c.strand));

    cur = g_slist_next(cur);
  }
}


/**
 * Writes a representation of seqfeatures in a flatfile
 * format that is compatible with the way genes and transcripts
 * are written.
 */
void seqfeat_slist_write_flatfile(FILE *fh, GSList *seqfeats) {
  SeqFeature *feat;
  GSList *cur;
  char *seqname;

  cur = seqfeats;


  while(cur != NULL) {
    feat = (SeqFeature *)cur->data;

    if(feat->c.chr != NULL && feat->c.chr->name != NULL) {
      seqname = feat->c.chr->name;
    } else if(feat->c.seqname != NULL) {
      seqname = feat->c.seqname;
    } else {
      seqname = "";
    }

    fprintf(fh, "%s\t%ld\t%s\t%s\t%lu\t%lu\t%d\t%g\t\n", 
	    SEQ_FEAT_LINE_PREFIX,
	    feat->id, 
	    (feat->name == NULL) ? "" : feat->name, 
	    seqname, 
	    feat->c.start, feat->c.end,
	    feat->c.strand,
	    feat->score);

    cur = g_slist_next(cur);
  }
}




/**
 * Writes a flatfile line representation of the provided feature to a
 * filehandle.
 */
void seqfeat_write_flatfile_line(FILE *fh, const SeqFeature *feat) {
  char *seqname;

  if(feat->c.chr != NULL && feat->c.chr->name != NULL) {
    seqname = feat->c.chr->name;
  } else if(feat->c.seqname != NULL) {
      seqname = feat->c.seqname;
  } else {
      seqname = "";
  }
  
  fprintf(fh, "%s\t%ld\t%s\t%s\t%lu\t%lu\t%d\t%g\n", 
	  SEQ_FEAT_LINE_PREFIX,
	  feat->id, 
	  (feat->name == NULL) ? "" : feat->name, 
	  seqname, 
	  feat->c.start, feat->c.end, feat->c.strand, feat->score);

}


/**
 * Writes a representation of seqfeatures in a flatfile
 * format that is compatible with the way genes and transcripts
 * are written.
 */
void seqfeat_write_flatfile(FILE *fh, SeqFeature *feats, const long n) {
  long i;
  for(i = 0; i < n; i++) {
    seqfeat_write_flatfile_line(fh, &feats[i]);
  }
}




void seqfeat_write(FILE *fh, SeqFeature *feat, 
		    char **attrib_names, const int n_attrib) {
  
  char *seqname, *str;
  int i;
  
  if(feat->c.chr != NULL && feat->c.chr->name != NULL) {
    seqname = feat->c.chr->name;
  } else if(feat->c.seqname != NULL) {
    seqname = feat->c.seqname;
  } else {
    seqname = "";
  }

  fprintf(fh, "%s\t%ld\t%ld\t%d\t%g",
	  seqname, feat->c.start, feat->c.end, feat->c.strand,
	  feat->score);

  if(feat->name) {
    fprintf(fh, "\t%s", feat->name);
  }

  for(i = 0; i < n_attrib; i++) {
    str = seqfeat_get_attrib_str(feat, attrib_names[i]);
    fprintf(fh, "\t%s", str);
  }
  fprintf(fh, "\n");  
}





void seqfeat_write_header(FILE *fh, char **attrib_names, const int n_attrib) {
  int i;
  char *hdr;

  fprintf(fh, "CHR\tSTART\tEND\tSTRAND\tSCORE\tNAME");

  for(i = 0; i < n_attrib; i++) {
    hdr = g_strdup(attrib_names[i]);
    util_str_replace(hdr, '_', '.');
    util_str_uc(hdr);
    fprintf(fh, "\t%s", hdr);
    g_free(hdr);
  }
  
  fprintf(fh, "\n");
}


/**
 *  Reads an array of SeqFeature structures from a file. The value
 *  pointed to the by the n_read argument is set to the number of
 *  SeqFeatures read.
 * 
 * The lines should be formatted as follows: 
 *   <chr>\t<start>\t<end>\t[<strand>\t[<score>\t[<name>\t[<attrib1>\t[<attrib2>\t[...]]]]]]\n
 *
 * For example, the following line could define a repeat feature 
 * with two extra attributes (with space instead of tab delimiters):
 *   chr12 1324151 1324222 0 0.0 AluSx Alu SINE 
 *
 * At a the very minimum, a feature must specify a chromosome, start,
 * and end, but can also specify an arbitrary number of named
 * attributes. If named attributes are to be specified, the strand, score and
 * name must also be specified. The name can be an empty string, in
 * which case the name is set to NULL. The strand can be 0 (indicating
 * no strand).
 */
SeqFeature *seqfeat_read_file(const char *filename, char **attrib_names,
			      const int n_attrib, long *n_read) {
  char *line;
  char **tok;
  int n_tok, min_tok, max_tok, strand, j, more_tok_warning;
  long n_line, i;
  SeqFeature *feats;

  FILE *fh;
  
  fh = fopen(filename, "r");
  if(fh == NULL) {
    g_error("%s:%d: could not read from file '%s'", __FILE__, __LINE__, 
	    filename);
  }
  
  /* count number of lines in file  */
  n_line = util_fcount_lines(fh);
  feats = g_new(SeqFeature, n_line);
  
  /* determine maximum and minimum number of tokens per line */
  if(n_attrib > 0) {
    min_tok = 6 + n_attrib;
    max_tok = 6 + n_attrib;
  } else {
    min_tok = 3;
    max_tok = 6;
  }

  more_tok_warning = FALSE;

  i = 0;
  while(i < n_line) {
    if((line = util_fgets_line(fh)) == NULL) {
      g_error("%s:%d: Expected %lu lines but got only %ld", __FILE__, 
	      __LINE__, n_line, i);
    }

    tok = g_strsplit(line, "\t", max_tok+1);

    n_tok = 0;
    while(tok[n_tok] != NULL) {
      n_tok++;
    }
    
    if(n_tok < min_tok) {
      g_error("%s:%d: Expected between %d and %d tokens per line, got %d",
	      __FILE__, __LINE__, min_tok, max_tok, n_tok);
    }
    if((n_tok > max_tok) && !more_tok_warning) {
      more_tok_warning = TRUE;
      g_warning("%s:%d: Some lines have more tokens than expected. "
		"Some feature attributes may be ignored because they "
		"were not specified",
		__FILE__, __LINE__);
    }

    /* read chromosome, start, end */
    feats[i].c.seqname = g_strdup(tok[0]);
    feats[i].c.chr = NULL;
    feats[i].c.start   = strtol(tok[1], NULL, 10);
    feats[i].c.end     = strtol(tok[2], NULL, 10);

    if(n_tok > 3) {
      /* read strand */
      strand = (tok[3][0] == '\0') ? STRAND_NONE : strtol(tok[3], NULL, 10);
      if((strand != STRAND_FWD) && (strand != STRAND_REV) && 
	 (strand != STRAND_NONE)) {
	g_error("%s:%d: invalid strand (%s)", __FILE__, __LINE__, tok[3]);
      }
      feats[i].c.strand = strand;
    }
    
    if(n_tok > 4) {
      /* read score */
      feats[i].score = (tok[4][0] == '\0') ? 0.0 : strtod(tok[4], NULL);
    } else {
      feats[i].score = 0.0;
    }
    
    if(n_tok > 5) {
      /* read name */
      feats[i].name = (tok[5][0] == '\0') ?  NULL : g_strdup(tok[5]);
    } else {
      feats[i].name = NULL;
    }

    feats[i].n_sub_feat = 0;
    feats[i].sub_feats = NULL;
    feats[i].attrib = NULL;

    /* add attributes */
    for(j = 0; j < n_attrib; j++) {
      seqfeat_add_attrib(&feats[i], attrib_names[j], tok[6 + j]);
    }

    g_strfreev(tok);
    g_free(line);

    i++;
  }

  *n_read = i;

  fclose(fh);

  return feats;
}




/**
 *  Reads an array of SeqFeature structures from a flatfile.
 *  The value pointed to the num_red argument is set to the number
 *  of SeqFeauture's read.
 */
SeqFeature *seqfeat_read_flatfile(const char *filename, long *n_read) {
  char *line;
  char **toks;
  int n_toks, prefix_len;
  long n_lines, i;
  SeqFeature *feats;

  FILE *fh;
  
  fh = fopen(filename, "r");
  if(fh == NULL) {
    g_error("%s:%d: could not read from file '%s'", 
	    __FILE__, __LINE__, filename);
  }

  /* count number of SEQFEAT lines in file  */
  n_lines = util_fcount_lines_match(fh, SEQ_FEAT_LINE_PREFIX);
  feats = g_new(SeqFeature, n_lines);

  prefix_len = strlen(SEQ_FEAT_LINE_PREFIX);

  i = 0;
  while(i < n_lines) {
    if((line = util_fgets_line(fh)) == NULL) {
      g_error("seqfeat_read_flatfile: Expected %lu %s lines, "
	      "got %lu", n_lines, SEQ_FEAT_LINE_PREFIX, i);
    }
    if(strncmp(SEQ_FEAT_LINE_PREFIX,line, prefix_len) != 0) {
      continue;
    }

    toks = g_strsplit(line, "\t", 8);

    n_toks = 0;
    while(toks[n_toks] != NULL) {
      n_toks++;
    }
    
    if(n_toks < 8) {
      g_error("gene_read_flatfile: Expected %d toks per gene line, got %d", 
	      8, n_toks);
    }

    feats[i].id      = strtol(toks[1], NULL, 10);
    feats[i].name    = g_strdup(toks[2]);
    feats[i].c.seqname = g_strdup(toks[3]);

    feats[i].c.start   = strtol(toks[4], NULL, 10);
    feats[i].c.end     = strtol(toks[5], NULL, 10);
    feats[i].c.strand  = strtol(toks[6], NULL, 10);
    feats[i].score = strtod(toks[7], NULL);

    feats[i].n_sub_feat = 0;
    feats[i].sub_feats = NULL;
    feats[i].attrib = NULL;

    g_strfreev(toks);
    g_free(line);

    i++;
  }

  *n_read = i;

  fclose(fh);

  return feats;
}




/**
 * Reads an array of SeqFeature structures from a BED file.
 * The value pointed to the num_red argument is set to the number
 * of SeqFeauture's read.
 */
SeqFeature *seqfeat_read_bed(const char *filename, long *n_read) {
  char *line;
  char **toks;
  int n_toks;
  SeqFeature *feats;
  long n_lines, i;  

  FILE *fh;

  fh = fopen(filename, "r");
  if(fh == NULL) {
    g_error("Could not read from BED file '%s'", filename);
  }

  n_lines = util_fcount_lines(fh);
  *n_read = 0;
  
  /* allocating feature for every line may be excessive because of
   * comments at the beginning of the file. However, it is much faster
   * to allocate all at once -- and not much space will be wasted
   * unless there is an unusually large number of comment lines
   */
  feats = g_new(SeqFeature, n_lines);

  for(i = 0; i < n_lines; i++) {
    if((line = util_fgets_line(fh)) == NULL) {
      /* Less lines than expected. File changed in length? */
      break;
    }

    /* skip comment lines beginning with ';'*/
    if(line[0] == ';') {
      g_free(line);
      continue;
    }

    toks = g_strsplit(line, "\t", 12);

    n_toks = 0;
    while(toks[n_toks] != NULL) {
      n_toks++;
    }
    
    if(n_toks < 3) {
      g_error("Expected at least 3 toks per BED file line, got %d", 
	      n_toks);
    }

    /* first three fields (chr, start, end) are required */
    feats[*n_read].c.seqname = g_strdup(toks[0]);
    feats[*n_read].c.start = strtol(toks[1], NULL, 10) + 1;
    feats[*n_read].c.end   = strtol(toks[2], NULL, 10);

    /* optional fourth field is feature name */
    if(n_toks >= 4) {
      feats[*n_read].name = g_strdup(toks[3]);
    } else {
      feats[*n_read].name = NULL;
    }

    /* optional fifth field is score */
    if(n_toks >= 5) {
      feats[*n_read].score = strtod(toks[4], NULL);
    } else {
      feats[*n_read].score = 0.0;
    }

    /* optional sixth field is strand */
    if(n_toks >= 6) {
      feats[*n_read].c.strand = char_to_strand(toks[5][0]);
    } else {
      feats[*n_read].c.strand = STRAND_NONE;
    }

    feats[*n_read].attrib = NULL;

    *n_read += 1;

    g_strfreev(toks);
    g_free(line);
  }

  fclose(fh);
	  
  return feats;
}


/**
 * Given two sorted feature arrays, looks for features that are
 * overlapping. An array is returned that is the same length as first
 * array. Each element is a singly-linked list containing ptrs to
 * overlapping SeqFeatures from the second array.  If no features
 * overlapped a SeqFeature the corresponding array element is NULL.
 *
 * If the cmp_strand argument is TRUE features that are on opposite
 * strands are not considered as overlapping. If the cmp_strand
 * argument is FALSE, features on opposite strands are considered
 * as overlapping. In this case the SeqFeatures must be pre-sorted
 * using the seqfeat_cmp_nostrand function rather than the seqfeat_cmp
 * function.
 * 
 */
GSList **seqfeat_overlaps(SeqFeature *feats1, const long feats1_sz, 
			  SeqFeature *feats2, const long feats2_sz,
			  const int cmp_strand) {

  long i, j, j_ovlp;
  GSList **overlaps;
  
  overlaps = g_new(GSList *, feats1_sz);
  for(i = 0; i < feats1_sz; i++) {
    overlaps[i] = NULL;
  }
  
  i = j = 0;
  while(i < feats1_sz && j < feats2_sz) {
    /* move along in both feature arrays until reaching the end of one
     * of the arrays or until features potentially overlap
     */

    if(cmp_strand) {
      while(i < feats1_sz && 
	    !seq_coord_ovlp(&feats1[i].c, &feats2[j].c, cmp_strand) && 
	    (seqfeat_cmp(&feats1[i], &feats2[j]) == -1)) {
	i++;
      }
      
      while(i < feats1_sz && j < feats2_sz && 
	    !seq_coord_ovlp(&feats1[i].c, &feats2[j].c, cmp_strand) && 
	    (seqfeat_cmp(&feats1[i], &feats2[j]) == 1)) {
	j++;
      }
    } else {
      while(i < feats1_sz && 
	    !seq_coord_ovlp(&feats1[i].c, &feats2[j].c, cmp_strand) && 
	    (seqfeat_cmp_nostrand(&feats1[i], &feats2[j]) == -1)) {

/* 	fprintf(stderr, "  %ld-%ld does not overlap %ld-%ld\n", */
/* 		feats1[i].c.start, feats1[i].c.end,  */
/* 		feats2[j].c.start, feats2[j].c.end); */

	i++;
      }
      
      while(i < feats1_sz && j < feats2_sz && 
	    !seq_coord_ovlp(&feats1[i].c, &feats2[j].c, cmp_strand) && 
	    (seqfeat_cmp_nostrand(&feats1[i], &feats2[j]) == 1)) {
	
/* 	fprintf(stderr, "  %ld-%ld does not overlap %ld-%ld\n", */
/* 		feats1[i].c.start, feats1[i].c.end,  */
/* 		feats2[j].c.start, feats2[j].c.end); */

	j++;
      }
    }

    /* keep adding features while they overlap with the feature
     * from the first list
     */
    j_ovlp = j;
    while(j_ovlp < feats2_sz && i < feats1_sz &&
	  seq_coord_ovlp(&feats1[i].c, &feats2[j_ovlp].c, cmp_strand)) {

/*       fprintf(stderr, "%ld-%ld overlaps %ld-%ld\n", */
/* 	      feats1[i].c.start, feats1[i].c.end,  */
/* 	      feats2[j_ovlp].c.start, feats2[j_ovlp].c.end); */

      overlaps[i] = g_slist_append(overlaps[i], &feats2[j_ovlp]);
      j_ovlp++;
    }
    /* advance ptr to first list, since we have found all features
     * that overlapped with this feature from second list
     */
    i++; 
  }

  return overlaps;
}



/*
 * Comparison function for SeqFeatures used for sorting.  Uses the
 * seq_coord_cmp function to do a coordinate comparison
 */
int seqfeat_cmp(const void *p1,const void *p2) {
  return seq_coord_cmp(&((SeqFeature *)p1)->c, &((SeqFeature *)p2)->c);
}



/*
 * Comparison function for SeqFeatures used for sorting. Uses the
 * seq_coord_cmp_nostrand function to do a coordinate comparison.
 */
int seqfeat_cmp_nostrand(const void *p1,const void *p2) {
  return seq_coord_cmp_nostrand(&((SeqFeature *)p1)->c,
				&((SeqFeature *)p2)->c);
}




void seqfeat_write_attribs(const SeqFeature *sf, FILE *fh) {
  SeqFeatAttrib *attr;

  attr = sf->attrib;
  
  while(attr != NULL) {
    fprintf(fh, "  %s=%s\n", attr->name, attr->value);
    
    attr = attr->next;
  }
}


/*
 * Adds a new named attribute with the provided name and value. The
 * name and value are duplicated.
 */
void seqfeat_add_attrib(SeqFeature *sf, const char *name, const char *val) {
  SeqFeatAttrib *attr;

  /*** TODO: use GTree instead of list ***/

  /* create new attribute */
  attr = g_new(SeqFeatAttrib, 1);
  attr->name = g_strdup(name);
  attr->value = g_strdup(val);

  /* add to front of attribute list */
  attr->next = sf->attrib;
  sf->attrib = attr;
}


/*
 * Converts the provided double into a string and adds it as an
 * attribute to this feature.
 */
void seqfeat_add_attrib_double(SeqFeature *sf, const char *name,
			       const double val) {
  char buf[100];

// FIXFIX
    sprintf(buf, "%g", val);
//  snprintf(buf, 100, "%g", val);
  seqfeat_add_attrib(sf, name, buf);
}


/*
 * Converts the provided integer into a string and adds it as an
 * attribute to this feature.
 */
void seqfeat_add_attrib_long(SeqFeature *sf, const char *name,
			       const long val) {
  char buf[100];

// FIXFIX
  sprintf(buf, "%ld", val);
//  snprintf(buf, 100, "%ld", val);
  seqfeat_add_attrib(sf, name, buf);
}


/**
 * Retrieve the named attribute from the provided sequence feature
 * or NULL if the attribute does not exist.
 */
SeqFeatAttrib *seqfeat_get_attrib(const SeqFeature *sf, const char *name) {
  SeqFeatAttrib *attr;

  attr = sf->attrib;

  /*** TODO: use GTree instead of list ***/

  while(attr != NULL) {
    if(strcmp(attr->name, name) == 0) {
      /* found the requested attribute */
      return attr;
    }

    attr = attr->next;
  }

  /* did not find attribute */
  return NULL;
}



/**
 * Returns true if the seqfeature has an attribute
 * of the provided name, returns false otherwise.
 */
int seqfeat_has_attrib(const SeqFeature *sf, const char *name) {
  SeqFeatAttrib *attr;
  attr = seqfeat_get_attrib(sf, name);

  if(attr == NULL) {
    return FALSE;
  }
  return TRUE;
}


/**
 * Gets the value of the provided seq features named
 * attribute. Returns NULL if attribute with name is not found.
 */
char *seqfeat_get_attrib_str(const SeqFeature *sf, const char *name) {
  SeqFeatAttrib *attr;

  attr = seqfeat_get_attrib(sf, name);

  if(attr == NULL) {
    fprintf(stderr, "Attributes are:\n");
    seqfeat_write_attribs(sf, stderr);
    g_error("%s:%d: no attribute with name '%s' exists.",
	    __FILE__, __LINE__, name);
  }
  
  return attr->value;
}


/**
 * Returns the value of the named attribute as a long integer
 */
long seqfeat_get_attrib_long(const SeqFeature *sf, const char *name) {
  char *str;

  str = seqfeat_get_attrib_str(sf, name);
  return strtol(str, NULL, 10);
}


/**
 * Returns the value of the named attribute as a double
 */
double seqfeat_get_attrib_double(const SeqFeature *sf, const char *name) {
  char *str;

  str = seqfeat_get_attrib_str(sf, name);
  return strtod(str, NULL);
}



/**
 * Returns the value of the named attribute as nucleotide identified
 */
unsigned char seqfeat_get_attrib_nuc(const SeqFeature *sf, const char *name) {
  char *str;

  str = seqfeat_get_attrib_str(sf, name);

  if(strlen(str) != 1) {
    g_error("%s:%d: attribute '%s' value '%s' does not look like nucleotide",
	    __FILE__, __LINE__, name, str);
  }

  return nuc_char_to_id(str[0]);
}
