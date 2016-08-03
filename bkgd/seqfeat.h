#ifndef __SEQ_FEAT_H__
#define __SEQ_FEAT_H__

#include <math.h>
#include <glib.h>
#include <stdlib.h>
#include <stdio.h>

#include "seqcoord.h"

#define SEQ_FEAT_LINE_PREFIX "SEQFEAT"

/* flag for undefined scores */
#define SEQFEAT_SCORE_NA -1.0


typedef struct SeqFeatAttrib_struct SeqFeatAttrib;

struct SeqFeatAttrib_struct {
  char *name;
  char *value;

  SeqFeatAttrib *next;
};


typedef struct SeqFeature_struct SeqFeature;

struct SeqFeature_struct {
  SeqCoord c;  /* genomic coordinate */
  long id;     /* database identifier */
  char *name;  /* name of feature */
  float score; /* score of feature */

  SeqFeatAttrib *attrib; /* ptr to list of other named attributes */

  SeqFeature *sub_feats; /* array of sub-features */
  long n_sub_feat;
};



void seqfeat_init(SeqFeature *feat);

SeqFeature *seqfeat_new(const char *seqname, long start, long end, 
			short strand, const char *name, float score);


void seqfeat_copy(const SeqFeature *src, SeqFeature *dest);
void seqfeat_copy_attrib(const SeqFeature *src, SeqFeature *dest);
void seqfeat_copy_subfeat(const SeqFeature *src, SeqFeature *dest);


	     
void seqfeat_slist_free(GSList *seqfeats);
void seqfeat_array_free(SeqFeature *seqfeats, long n);
void seqfeat_slist_write_bed(FILE *fh, GSList *seqfeats);

void seqfeat_slist_write_flatfile(FILE *fh, GSList *seqfeats);
void seqfeat_write_flatfile(FILE *fh, SeqFeature *seqfeats, long n);
void seqfeat_write_flatfile_line(FILE *fh, const SeqFeature *feat);

SeqFeature *seqfeat_read_file(const char *filename, char **attrib_names,
			      const int n_attrib, long *n_read);

SeqFeature *seqfeat_read_flatfile(const char *filename, long *n_read);

SeqFeature *seqfeat_read_bed(const char *filename, long *n_read);

GSList **seqfeat_overlaps(SeqFeature *feats1, const long feats1_sz, 
			  SeqFeature *feats2, const long feats2_sz,
			  const int cmp_strand);

int seqfeat_cmp(const void *p1,const void *p2);
int seqfeat_cmp_nostrand(const void *p1,const void *p2);


void seqfeat_add_attrib(SeqFeature *sf, const char *name, const char *value);
void seqfeat_add_attrib_double(SeqFeature *sf, const char *name, 
			       const double val);
void seqfeat_add_attrib_long(SeqFeature *sf, const char *name, const long value);
SeqFeatAttrib *seqfeat_get_attrib(const SeqFeature *sf, const char *name);
int seqfeat_has_attrib(const SeqFeature *sf, const char *name);
char *seqfeat_get_attrib_str(const SeqFeature *sf, const char *name);
long seqfeat_get_attrib_long(const SeqFeature *sf, const char *name);
double seqfeat_get_attrib_double(const SeqFeature *sf, const char *name);
unsigned char seqfeat_get_attrib_nuc(const SeqFeature *sf, const char *name);


#endif
