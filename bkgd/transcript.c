#include <stdio.h>
#include <glib.h>
#include <string.h>

#include "transcript.h"
#include "seqfeat.h"
#include "seqcoord.h"
#include "aa.h"

#include "util.h"

#define TRANSCRIPT_GENE_FIELDS 12


/*
 * internal function to facilitate sorting of a transcript's exons
 */
int ex_cmp(const void *v1, const void *v2) {
  const SeqCoord *ex1, *ex2;
  
  ex1 = v1;
  ex2 = v2;

  /* if exons on reverse strand, then sort coordinates in descending order */
  if(ex1->strand == STRAND_REV) {
    if(ex1->start > ex2->start) {
      return -1;
    }
    if(ex1->start < ex2->start) {
      return 1;
    }
    return 0;
  }

  /* if exons on fwd strand, then sort coordinates in ascending order */
  if(ex1->start > ex2->start) {
    return 1;
  }
  if(ex1->start < ex2->start) {
    return -1;
  }
  return 0;
}


/**
 * Parses a set of comma-delimited exon start and end strings. This
 * function allows for the presence or absence of a trailing comma. If
 * the coordinates are UCSC-style zero-based coordinates the
 * zero_based flag should be set to TRUE. This function sets the exons
 * and num_exons value of the provided transcript.
 *
 * The exons are always stored in 5' to 3' order. This means that for
 * reverse strand transcripts the first exon (index 0) has the highest
 * start coordinate (coordinates are always relative to the forward
 * strand of the chromosome).
 * 
 */
void transcript_parse_exons(Transcript *tr, char *ex_starts, char *ex_ends,
			    int zero_based) {
  int i, n_exons, ordered;
  char *cur_ptr, *end_ptr;
  long prev_start, cur_start;

  i = 0;

  /* count the number of exons */
  if(ex_starts[i] != '\0' && ex_starts[i] != ',') {
    n_exons = 1;
  } else {
    n_exons = 0;
  }
  while(ex_starts[i] != '\0') {
    if(ex_starts[i] == ',') {
      if(ex_starts[i+1] != '\0') { /* can be trailing comma */
	n_exons++;
      }
    }
    i++;
  }

  tr->num_exons = n_exons;
  tr->exons = g_new(SeqCoord, n_exons);


  /* parse exon starts from comma-delimited string */
  prev_start = -1;
  cur_start = -1;
  ordered = TRUE;

  cur_ptr = ex_starts;
  for(i = 0; i < n_exons; i++) {
    cur_start = strtol(cur_ptr, &end_ptr, 10);
    if(zero_based) {
      cur_start += 1;
    }

    /* check if exons are ordered */
    if(prev_start != -1) {
      if(tr->c.strand == STRAND_REV) {
	if(prev_start <= cur_start) {
	  ordered = FALSE;
	}
      }
      else if(tr->c.strand == STRAND_FWD) {
	if(prev_start >= cur_start) {
	  ordered = FALSE;
	}
      }
      else {
	g_error("transcript_parse_exons: transcript strand "
		"must be defined\n");
      }
    }
    if(end_ptr[0] == '\0') {
      if(i != n_exons-1) {
	g_error("transcript_parse_exons: exon_start string terminated "
		"earlier than expected");
      }
    }
    else if(end_ptr[0] != ',') {
      g_error("transcript_parse_exons: expected comma delimiters, got '%c'",
	      end_ptr[0]);
    }

    tr->exons[i].start = cur_start;
    prev_start = cur_start;

    cur_ptr = &end_ptr[1];
  }

  /* parse exon ends */
  cur_ptr = ex_ends;
  for(i = 0; i < n_exons; i++) {
    tr->exons[i].end = strtol(cur_ptr, &end_ptr, 10);
    cur_ptr = &end_ptr[1];
  }

  /* set exon strands chrs, and seqnames to match transcripts */
  for(i = 0; i < n_exons; i++) {
    tr->exons[i].strand = tr->c.strand;

    if(tr->c.seqname == NULL) {
      tr->exons[i].seqname = NULL;
    }
    else {
      tr->exons[i].seqname = g_strdup(tr->c.seqname);
    }
    tr->exons[i].chr = tr->c.chr;
  }


  /* re-order exons if they were not ordered */
  if(!ordered) {
    qsort(tr->exons, n_exons, sizeof(SeqCoord), ex_cmp);
  }

}


/**
 * Reads an array of transcripts from the UCSC knownGene mySQL table
 * file. The value pointed to by the num_read argument is set to
 * number of transcripts read.  Offset can be set to 1 if there is an
 * extra "bin" column as the first column (recently added field to their
 * database), otherwise this should be 0.
 * 
 */
Transcript *transcript_read_gene_file(char *table_file, long *num_read) {
  char *line;
  char **toks;
  int num_toks, n_exons;
  Transcript *trs;
  long num_lines, i;

  FILE *fh;

  fh = fopen(table_file, "r");
  if(fh == NULL) {
    g_error("Could not read from table file '%s'", table_file);
  }

  num_lines = util_fcount_lines(fh);
  trs = g_new(Transcript, num_lines);

  for(i = 0; i < num_lines; i++) {
    if((line = util_fgets_line(fh)) == NULL) {
      /* Less lines than expected */
      g_error("transcript_read_gene_file: expected %ld lines but got only "
	      "%ld\n", num_lines, i);
    }

    toks = g_strsplit(line, "\t", TRANSCRIPT_GENE_FIELDS);

    num_toks = 0;
    while(toks[num_toks] != NULL) {
      num_toks++;
    }
    
    if(num_toks < 10) {
      g_error("transcript_read_gene_file: expected %d toks per "
	      "known gene line, got %d", num_toks, TRANSCRIPT_GENE_FIELDS);
    }

    trs[i].name    = g_strdup(toks[0]);
    trs[i].c.seqname = g_strdup(toks[1]);
    trs[i].c.strand  = char_to_strand(toks[2][0]);
    trs[i].c.start   = strtol(toks[3], NULL, 10) + 1;
    trs[i].c.end     = strtol(toks[4], NULL, 10);
    trs[i].c.chr     = NULL;

    trs[i].cds.start = strtol(toks[5], NULL, 10) + 1;
    trs[i].cds.end   = strtol(toks[6], NULL, 10);
    trs[i].cds.seqname = g_strdup(trs[i].c.seqname);
    trs[i].cds.strand  = trs[i].c.strand;
    trs[i].cds.chr = NULL;

    n_exons = strtol(toks[7], NULL, 10);
    
    trs[i].method = NULL;

    transcript_parse_exons(&trs[i], toks[8], toks[9], TRUE);

    if(n_exons != trs[i].num_exons) {
      g_error("transcript_read_gene_file: expected %d exons, but got "
	      "%d (line %ld)", n_exons, trs[i].num_exons, i+1);
    }
    
    g_strfreev(toks);
    g_free(line);
  }

  *num_read = i;

  fclose(fh);
	  
  return trs;
}


/**
 * Frees the memory allocated for a transcript. Does not free attached
 * Chromosome structures.
 */
void transcript_free(Transcript *tr) {
  int i;

  /* free exons */
  for(i = 0; i < tr->num_exons; i++) {
    if(tr->exons[i].seqname != NULL) {
      g_free(tr->exons[i].seqname);
    }
  }
  g_free(tr->exons);

  /* free CDS */
  if(tr->cds.seqname != NULL) {
    g_free(tr->cds.seqname);
  }

  /* free transcript */
  if(tr->c.seqname != NULL) {
    g_free(tr->c.seqname);
  }

  g_free(tr->name);
  g_free(tr);

}


/**
 * Frees all of the transcripts in an array and the array itself.
 */
void transcript_array_free(Transcript *trs, long num) {
  long i;
  int j;

  if(num == 0) {
    return;
  }

  for(i = 0; i < num; i++) {
    /* free exons */
    for(j = 0; j < trs[i].num_exons; j++) {
      g_free(trs[i].exons[j].seqname);
    }
    g_free(trs[i].exons);

    /* free CDS */
    if(trs[i].cds.seqname != NULL) {
      g_free(trs[i].cds.seqname);
    }

    /* free transcript */
    if(trs[i].c.seqname != NULL) {
      g_free(trs[i].c.seqname);
    }
    g_free(trs[i].name);
  }
  
  g_free(trs);
}


/**
 * Writes a BED line representation of a transcript to a provided
 * filehandle. BED format is defined here: 
 *    http://genome.ucsc.edu/goldenPath/help/customTrack.html#BED
 * 
 * The first three required BED fields are:
 *
 * 1. chrom - The name of the chromosome 
 * 2. chromStart - The starting position of the feature in the chromosome or 
 *         scaffold. The first base in a chromosome is numbered 0.
 * 3. chromEnd - The ending position of the feature in the chromosome or 
 *         scaffold. The chromEnd base is not included in the display 
 *         of the feature. For example, the first 100 bases of a 
 *         chromosome are defined as chromStart=0, chromEnd=100, 
 *         and span the bases numbered 0-99. 
 *
 * The 9 additional optional BED fields are:
 * 
 * 4. name - Defines the name of the BED line. This label is displayed to 
 *         the left of the BED line in the Genome Browser window when the 
 *         track is open to full display mode or directly to the left of the 
 *         item in pack mode.
 * 5. score - A score between 0 and 1000. If the track line useScore attribute
 *         is set to 1 for this annotation data set, the score value will 
 *         determine the level of gray in which this feature is displayed 
 *         (higher numbers = darker gray).
 * 6. strand - Defines the strand - either '+' or '-'.
 * 7. thickStart - The starting position at which the feature is drawn thickly 
 *         (for example, the start codon in gene displays).
 * 8. thickEnd - The ending position at which the feature is drawn thickly 
 *         (for example, the stop codon in gene displays).
 * 9. itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track 
 *         line itemRgb attribute is set to "On", this RBG value will 
 *         determine the display color of the data contained in this BED 
 *         line. NOTE: It is recommended that a simple color scheme (eight 
 *         colors or less) be used with this attribute to avoid overwhelming 
 *         the color resources of the Genome Browser and your Internet 
 *         browser. At the present time, this option is functional only on 
 *         BED tracks with 9 or less fields.
 * 10. blockCount - The number of blotocks (exons) in the BED line.
 * 11. blockSizes - A comma-separated list of the block sizes. The number of 
 *         items in this list should correspond to blockCount.
 * 12. blockStarts - A comma-separated list of block starts. All of the 
 *         blockStart positions should be calculated relative to chromStart. 
 *         The number of items in this list should correspond to blockCount. 
 *
 */
inline void transcript_write_bed(FILE *fh, Transcript *tr) {
  char *rgb = "255,0,0"; /* required track color field */
  int i;
  char *seqname;

  if(tr->c.chr != NULL && tr->c.chr->name != NULL) {
    seqname = tr->c.chr->name;
  } else if(tr->c.seqname != NULL) {
    seqname = tr->c.seqname;
  } else {
    seqname = "";
  }

  fprintf(fh, "%s\t%lu\t%lu\t%s\t%d\t%c\t%lu\t%lu\t%s\t%d\t", seqname, 
	  tr->c.start - 1, tr->c.end, tr->name, 1000,
	    strand_to_char(tr->c.strand),
	  tr->cds.start - 1, tr->cds.end, rgb,
	  tr->num_exons);

  /* write out block sizes for exons */
  for(i = 0; i < tr->num_exons; i++) {
    fprintf(fh, "%lu,", tr->exons[i].end - tr->exons[i].start + 1);
  }

  fprintf(fh,"\t");
  
  /* write out block starts for exons */
  for(i = 0; i < tr->num_exons; i++) {
    fprintf(fh, "%lu,", tr->exons[i].start - tr->c.start);
  }

  fprintf(fh, "\n");
}



/**
 * Outputs a tab-delimited flat file representation of an array of
 * transcripts to the provided file handle.
 */
void transcript_write_flatfile(FILE *fh, Transcript *tr_array, 
			       long num_trs) {
  long i;
  
  /* assign transcripts and exons identifiers and print them out */
  for(i = 0; i < num_trs; i++) {
    transcript_write_flatfile_line(fh, &tr_array[i]);
  }
}



/**
 * Outputs a representation of a single transcript with tab-delimited
 * fields the the provided file handle.
 */
inline void transcript_write_flatfile_line(FILE *fh, Transcript *tr) {
  long i;
  char *seqname;
  int first;

  if(tr->c.chr != NULL && tr->c.chr->name != NULL) {
    seqname = tr->c.chr->name;
  } else if(tr->c.seqname != NULL) {
    seqname = tr->c.seqname;
  } else {
    seqname = "";
  }

  fprintf(fh, "TRANSCRIPT\t%ld\t%s\t%s\t%lu\t%lu\t%d\t%lu\t%lu\t%d\t", 
	  tr->id, tr->name, seqname, tr->c.start, tr->c.end, 
	  tr->c.strand, tr->cds.start, tr->cds.end, tr->num_exons);
  

  /* write out starts for exons */
  first = TRUE;
  for(i = 0; i < tr->num_exons; i++) {
    if(first) {
      first = FALSE;
    } else {
      fprintf(fh,",");
    }
    fprintf(fh, "%lu", tr->exons[i].start);
  }
  
  fprintf(fh,"\t");
      
  /* write out ends for exons */
  first = TRUE;
  for(i = 0; i < tr->num_exons; i++) {
    if(first) {
      first = FALSE;
    } else {
      fprintf(fh,",");
    }
    fprintf(fh, "%lu", tr->exons[i].end);
  }
  
  fprintf(fh, "\n");
}


/**
 * Reads a flatfile formatted file of transcripts. The variable
 * pointed to be the num_read argument is set to the number
 * of transcripts that were read.
 */
Transcript *transcript_read_flatfile(char *filename, long *num_read) {
  char *line;
  char **toks;
  int num_toks, prefix_len;
  Transcript *trs;
  long num_lines;
  long i;

  FILE *fh;

  fh = fopen(filename, "r");
  if(fh == NULL) {
    g_error("transcript_read_flatfile: could not read from file '%s'", 
	    filename);
  }

  /* count number of TRANSCRIPT lines in file */
  num_lines = util_fcount_lines_match(fh, TRANSCRIPT_LINE_PREFIX);
  trs = g_new(Transcript, num_lines);

  prefix_len = strlen(TRANSCRIPT_LINE_PREFIX);

  i = 0;
  while(i < num_lines) {
    if((line = util_fgets_line(fh)) == NULL) {
      g_error("transcript_read_flatfile: Expected %lu %s lines, "
	      "got %lu", num_lines, TRANSCRIPT_LINE_PREFIX, i);
    }
    if(strncmp(TRANSCRIPT_LINE_PREFIX,line, prefix_len) != 0) {
      continue;
    }

    toks = g_strsplit(line, "\t", 12);

    num_toks = 0;
    while(toks[num_toks] != NULL) {
      num_toks++;
    }
    
    if(num_toks < 12) {
      g_error("transcript_read_flatfile: Expected %d toks per transcript line,"
	      " got %d", 12, num_toks);
    }

    trs[i].id      = strtol(toks[1], NULL, 10);
    trs[i].name    = g_strdup(toks[2]);
    trs[i].c.seqname = g_strdup(toks[3]);
    trs[i].c.chr = NULL;

    trs[i].c.start   = strtol(toks[4], NULL, 10);
    trs[i].c.end     = strtol(toks[5], NULL, 10);
    trs[i].c.strand  = strtol(toks[6], NULL, 10);

    trs[i].cds.start = strtol(toks[7], NULL, 10);
    trs[i].cds.end   = strtol(toks[8], NULL, 10);
    trs[i].cds.seqname = g_strdup(trs[i].c.seqname);
    trs[i].cds.strand  = trs[i].c.strand;
    trs[i].cds.chr = NULL;

    trs[i].num_exons = strtol(toks[9], NULL, 10);
    
    trs[i].method = NULL;
    
    transcript_parse_exons(&trs[i], toks[10], toks[11], FALSE);

    g_strfreev(toks);
    g_free(line);

    i++;
  }

  *num_read = i;

  fclose(fh);

  return trs;
}



/*
 * Comparison function of transcripts used for sorting.  Uses the
 * seq_coord_cmp function to do a coordinate comparison.
 */
int transcript_cmp(const void *p1,const void *p2) {
  return seq_coord_cmp(&((Transcript *)p1)->c, &((Transcript *)p2)->c);
}


/*
 * Comparison function of transcripts used for sorting.  Uses the
 * seq_coord_cmp_nostrand function to do a coordinate comparison.
 */
int transcript_cmp_nostrand(const void *p1,const void *p2) {
  return seq_coord_cmp_nostrand(&((Transcript *)p1)->c,
				&((Transcript *)p2)->c);
}


/*
 * Comparison function of transcripts used for sorting. Uses the
 * end of the transcript for comparison, rather than the start.
 */
int transcript_cmp_end(const void *p1, const void *p2) {
  return seq_coord_cmp_end(&((Transcript *)p1)->c,
			   &((Transcript *)p2)->c);
}


/**
 * Returns the cDNA of the provided transcript, which is the sequence
 * of the exons spliced together.
 */
Seq *transcript_cdna_seq(Transcript *tr, Seq *seq) {
  return seq_subseq_coords(seq, tr->exons, tr->num_exons);
}



/**
 * Returns TRUE if transcript has protein coding sequence, 
 * FALSE if this is a non-coding transcript.
 */
int transcript_has_cds(Transcript *tr) {
  return (tr->cds.start < tr->cds.end);
}


/**
 * Returns the length of the coding sequence (CDS) of the provided
 * transcript. Returns 0 if transcript is non-coding.
 */
long transcript_cds_len(Transcript *tr) {
  long i;
  long cds_len, start, end;
  SeqCoord *exon;

  cds_len = 0;

  if(tr->cds.start > tr->cds.end) {
    /* transcript has no CDS */
    return 0;
  }

  /* figure out how long CDS is by summing lengths of exons that fall into
   * CDS
   */  
  for(i = 0; i < tr->num_exons; i++) {
    exon = &tr->exons[i];
    if(tr->cds.start > exon->end) {
      /* exon is before translation start */
      continue;
    }
    if(tr->cds.end < exon->start) {
      /* exon is after translation end */
      continue;
    }

    /* exon may only be partially translated */
    start = (tr->cds.start > exon->start) ? tr->cds.start : exon->start;
    end   = (tr->cds.end   < exon->end)   ? tr->cds.end   : exon->end;
    
    cds_len += end - start + 1;
  }

  return cds_len;
}



/**
 * Returns the coordinates of the CDS of this transcript. The
 * coordinate array should be freed when it is no longer needed using
 * seq_coord_array_free function. The returned coordinates are on the
 * same sequence strand as the transcript, and are ordered so that
 * coordinates from the 5'-most exon are first and 3'-most exon are
 * last.
 *
 * Returns NULL if this transcript is non-coding
 */
SeqCoord *transcript_cds(Transcript *tr, int *n_coords) {
  SeqCoord *cds, *exon;
  int i, j;

  *n_coords = 0;

  if(tr->cds.start > tr->cds.end) {
    /* transcript has no CDS */
    return NULL;
  }


  /* first pass: count exons that contain CDS */
  for(i = 0; i < tr->num_exons; i++) {
    if(tr->cds.start > tr->exons[i].end) {
      /* exon is before translation start */
      continue;
    }
    if(tr->cds.end < tr->exons[i].start) {
      /* exon is after CDS end */
      continue;
    }
    *n_coords += 1;    
  }

  if(*n_coords == 0) {
    g_error("transcript_cds: transcript %s does not have any CDS exons",
	    tr->name);
  }

  cds = g_new(SeqCoord, *n_coords);

  j = 0;

  for(i = 0; i < tr->num_exons; i++) {
    exon = &tr->exons[i];

    if(tr->cds.start > exon->end) {
      /* exon is past end of CDS  */
      continue;
    }
    if(tr->cds.end < exon->start) {
      /* exon is before start of CDS */
      continue;
    }    
    
    /* exon may only be partially translated */
    cds[j].start = (tr->cds.start > exon->start) ? tr->cds.start : exon->start;
    cds[j].end   = (tr->cds.end   < exon->end)   ? tr->cds.end   : exon->end;
    cds[j].strand = tr->c.strand;
    cds[j].seqname = NULL;
    cds[j].chr = tr->c.chr;
    j++;
  }

  return cds;
}



/**
 * Returns an array of long integers the same length as the CDS of
 * this transcript (or NULL if this is a non-coding transcript).
 * The returned array is a mapping of CDS coordinates back to
 * genomic coordinates (indexed by the CDS position). The returned
 * array should be freed when no longer required.
 */
long *transcript_cds2genome_map(Transcript *tr) {
  long *cds2gn, cds_len, cds_pos, gn_pos;
  int strand, cur_coord, n_cds_coord;
  SeqCoord *cds;

  strand = tr->c.strand;
  cds_len = transcript_cds_len(tr);
  cds = transcript_cds(tr, &n_cds_coord);

  if(cds_len < 1) {
    return NULL;
  }

  cds2gn = g_new(long, cds_len);

  cur_coord = 0;
  if(strand == STRAND_FWD) {
    gn_pos = cds[cur_coord].start;
  } else {
    gn_pos = cds[cur_coord].end;
  }

  for(cds_pos = 1; cds_pos <= cds_len; cds_pos++) {    
    cds2gn[cds_pos-1] = gn_pos;

    /* update genomic position */
    if(strand == STRAND_FWD) {
      gn_pos++;
      if(gn_pos > cds[cur_coord].end) {
	/* move ahead to next exon */
	cur_coord++;
	if(cur_coord >= n_cds_coord) {
	  /* terminate, at end of CDS */
	  if(cds_pos != cds_len) {
	    /* sanity check: we should be at last cds pos */
	    g_error("%s:%d: expected to be at end of CDS. "
		    "cds_pos=%ld, cds_len=%ld, cur_coord=%d, n_cds_coord=%d",
		    __FILE__, __LINE__, cds_pos, cds_len, 
		    cur_coord, n_cds_coord);
	  }
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
	  if(cds_pos != cds_len) {
	    g_error("%s:%d: expected to be at end of CDS. "
		    "cds_pos=%ld, cds_len=%ld, cur_coord=%d, n_cds_coord=%d",
		    __FILE__, __LINE__, cds_pos, cds_len, 
		    cur_coord, n_cds_coord);
	  }
	  break;
	}
	gn_pos = cds[cur_coord].end;
      }
    }
  }

  return cds2gn;
}


/**
 * Returns the coding sequence (CDS) of the provided transcript, given
 * a stretch of genomic sequence that spans (at least) the region of
 * the provided transcript. The returned sequence should be freed when
 * it is no longer needed using the seq_free function.
 *
 * Returns NULL if this transcript is non-coding.
 */
Seq *transcript_cds_seq(Transcript *tr, Seq *seq) {
  Seq *cds_seq;
  SeqCoord *cds_coord;
  int n;

  /* get genomic coordinates of CDS */
  cds_coord = transcript_cds(tr, &n);

  if(n == 0) {
    return NULL;
  }

  /* get sequence corresponding to coordinates */
  cds_seq = seq_subseq_coords(seq, cds_coord, n);

  seq_coord_array_free(cds_coord, n);

  return cds_seq;
}




/**
 * Returns the length of the 5'UTR for this transcript.
 */
long transcript_utr5_len(Transcript *tr) {
  SeqCoord *utr5_c;
  long len;
  int num_c;

  utr5_c = transcript_utr5(tr, &num_c);

  if(num_c == 0) {
    len = 0;
  } else {
    len = seq_coord_array_len(utr5_c, num_c);
    seq_coord_array_free(utr5_c, num_c);
  }

  return len;
}






/**
 * Returns the sequence of the 5'UTR for this transcript.
 * The returned Seq should be freed when it is no longer needed.
 * If there is no 5' UTR defined for this transcript NULL is
 * returned instead.
 *
 * The provided sequence should span (at least) the region of the
 * provided transcript.
 */
Seq *transcript_utr5_seq(Transcript *tr, Seq *seq) {
  SeqCoord *utr5_c;
  Seq *utr5_seq;
  int num_c;

  utr5_c = transcript_utr5(tr, &num_c);

  if(num_c == 0) {
    return NULL;
  }
  
  utr5_seq = seq_subseq_coords(seq, utr5_c, num_c);
  seq_coord_array_free(utr5_c, num_c);

  return utr5_seq;
}



/**
 * Returns an array of coordinates that indicate the genomic positions
 * of this transcript's 5' UTR. The value pointed to by the num_coords
 * argument is set to indicate the number of coordinates returned.
 * If this transcript has no 5' UTR NULL is returned.
 */
SeqCoord *transcript_utr5(Transcript *tr, int *num_coords) {
  SeqCoord *utr5_c, *exon;
  long i, j;

  /* first count the number exons that have 5'UTR sequence */
  *num_coords = 0;

  for(i = 0; i < tr->num_exons; i++) {
    exon = &tr->exons[i];
    if(tr->c.strand == STRAND_FWD) {  
      if(exon->start >= tr->cds.start) {
	break;
      }

      /* this exon is all or partially 5' UTR */
      *num_coords += 1;
    } else {
      if(exon->end <= tr->cds.end) {
	continue;
      }

      /* exon is all or partially 5' UTR */
      *num_coords += 1;
    }
  }

  if(*num_coords == 0) {
    return NULL;
  }
  
  utr5_c = g_new(SeqCoord, *num_coords);

  j = 0;
  for(i = 0; i < tr->num_exons; i++) {
    exon = &tr->exons[i];

    if(tr->c.strand == STRAND_FWD) {  
      if(exon->start >= tr->cds.start) {
	break;
      }
      
      /* this exon is all or partially 5' UTR */
      utr5_c[j].start = exon->start;
      utr5_c[j].end = 
	(exon->end >= tr->cds.start) ? tr->cds.start-1 : exon->end;
      utr5_c[j].strand = STRAND_FWD;
      utr5_c[j].seqname = NULL;
      utr5_c[j].chr = tr->c.chr;

      j++;
    } else {
      if(exon->end <= tr->cds.end) {
	break;
      }

      /* exon is all or partially 5' UTR */
      utr5_c[j].start = 
	(exon->start <= tr->cds.end) ? tr->cds.end+1 : exon->start;
      utr5_c[j].end = exon->end;
      utr5_c[j].strand = STRAND_REV;
      utr5_c[j].seqname = NULL;
      utr5_c[j].chr = tr->c.chr;
      
      j++;
    }
  }

  return utr5_c;  
}





/**
 * Returns the length of the 3'UTR for this transcript.
 */
long transcript_utr3_len(Transcript *tr) {
  SeqCoord *utr3_c;
  long len;
  int num_c;

  utr3_c = transcript_utr3(tr, &num_c);

  if(num_c == 0) {
    len = 0;
  } else {
    len = seq_coord_array_len(utr3_c, num_c);
    seq_coord_array_free(utr3_c, num_c);
  }

  return len;
}



/**
 * Returns the sequence of the 3' UTR for this transcript. The 
 * Seq should be freed when it is no longer needed. If there is
 * no 3' UTR defined for this transcript NULL is returned instead.
 *
 * The provided sequence should span (at least) the region of the
 * provided transcript.
 */
Seq *transcript_utr3_seq(Transcript *tr, Seq *seq) {
  SeqCoord *utr3_c;
  Seq *utr3_seq;
  int num_c;

  utr3_c = transcript_utr3(tr, &num_c);

  if(num_c == 0) {
    return NULL;
  }
  
  utr3_seq = seq_subseq_coords(seq, utr3_c, num_c);
  seq_coord_array_free(utr3_c, num_c);

  return utr3_seq;
}



/* Retrieves a genomic sequence coordinate representing a region
 * around the transcription start site (TSS). (num_bp-1)/2 flanking
 * base pairs are retrieved on each side of the transcription start
 * site as well, unless the transcript is less than num_bp basepairs
 * from the edge of the sequence (then as many bases as possible will
 * be included). The returned SeqCoord should be freed when it is no
 * longer needed.
 */
SeqCoord *transcript_tss(Transcript *tr, Seq *seq, 
			 long upstream, long dnstream) {
  SeqCoord *c;  

  c = g_new(SeqCoord,1);
  c->seqname = NULL;
  
  if(tr->c.strand == STRAND_FWD) {
    c->start  = tr->c.start - upstream;
    c->end    = tr->c.start + dnstream;
    c->strand = STRAND_FWD;
  } else {
    c->start  = tr->c.end - dnstream;
    c->end    = tr->c.end + upstream;
    c->strand = STRAND_REV;
  }

  if(c->start < 1) {
    c->start = 1;
  }   
  if(c->end > seq->len) {
    c->end = seq->len;
  }
  
  return c;
}




/* Retrieves a genomic sequence coordinate representing a region
 * upstream from the provided transcript's transcription start site of
 * length num_bp.  If the transcript is less than num_bp basepairs
 * from the edge of the sequence, the returned coordinate will consist
 * of less than num_bp bases. If the transcript is on the very edge of
 * the sequence, NULL is returned instead. The returned SeqCoord
 * should be freed when it is no longer needed.
 */
SeqCoord *transcript_upstream(Transcript *tr, Seq *seq, long num_bp) {
  SeqCoord *c;  

  c = g_new(SeqCoord,1);
  c->seqname = NULL;
  
  if(tr->c.strand == STRAND_FWD) {
    if(tr->c.start == 1) {
      g_free(c);
      return NULL;
    }

    c->start  = tr->c.start - num_bp;
    c->end    = tr->c.start - 1;
    c->strand = STRAND_FWD;

    if(c->start < 1) {
      c->start = 1;
    } 
  } else {
    if(tr->c.end == seq->len) {
      g_free(c);
      return NULL;
    }

    c->start  = tr->c.end + 1;
    c->end    = tr->c.end + num_bp;
    c->strand = STRAND_REV;
    
    if(c->end > seq->len) {
      c->end = seq->len;
    }
  }
  
  return c;
}



/**
 * Retrieves num_bp bases of upstream sequence from a transcript's
 * transcription start site. The sequence that is retrieved is from
 * the same strand that the transcript is on. If the transcript is
 * less than num_bp basepairs from the edge of the sequence, less than
 * num_bp of upstream sequence is returned. If the transcript is on
 * the very edge of the sequence, NULL is returned instead.
 *
 * The returned sequence should be freed when it is no longer needed.
 */
Seq *transcript_upstream_seq(Transcript *tr, Seq *seq, long num_bp) {
  SeqCoord *c;
  Seq *upstrm_seq;

  c = transcript_upstream(tr, seq, num_bp);
  if(c == NULL) {
    return NULL;
  }

  upstrm_seq = seq_subseq(seq, c);

  g_free(c);

  return upstrm_seq;
}


/* Retrieves genomic sequence from a region around the transcription
 * start site (TSS) of the provided transcript. The number of upstream
 * and downstream bases to retrieve are specified by the upstream and
 * dnstream bases arguments. If both of these are zero a single base
 * at the TSS is returned.
 */
Seq *transcript_tss_seq(Transcript *tr, Seq *seq, long upstream,
			long dnstream) {
  SeqCoord *c;
  Seq *tss_seq;

  c = transcript_tss(tr, seq, upstream, dnstream);
  if(c == NULL) {
    return NULL;
  }

  tss_seq = seq_subseq(seq,c);
  g_free(c);
  return tss_seq;
}




/**
 * Returns an array of genomic coordinates that represents the 3'UTR
 * region of this transcript. The number of coordinates returned is
 * assigned to the value pointed to by the num_coords argument.
 * If this transcript does not have a defined 3'UTR, NULL is returned
 * instead of an array of SeqCoords.  
 */
SeqCoord *transcript_utr3(Transcript *tr, int *num_coords) {
  SeqCoord *utr3_c, *exon;
  int i, j;

  /* first count the number exons that have 3'UTR sequence */
  *num_coords = 0;

  for(i = 0; i < tr->num_exons; i++) {
    exon = &tr->exons[i];
    if(tr->c.strand == STRAND_FWD) {  
      if(exon->end <= tr->cds.end) {
	continue;
      }

      /* this exon is all or partially 3' UTR */
      *num_coords += 1;
    } else {
      if(exon->start >= tr->cds.start) {
	continue;
      }

      /* exon is all or partially 3' UTR */
      *num_coords += 1;
    }
  }

  if(*num_coords == 0) {
    return NULL;
  }
  
  utr3_c = g_new(SeqCoord, *num_coords);

  j = 0;
  for(i = 0; i < tr->num_exons; i++) {
    exon = &tr->exons[i];

    if(tr->c.strand == STRAND_FWD) {  
      if(exon->end <= tr->cds.end) {
	continue;
      }
      
      /* this exon is all or partially 3' UTR */
      utr3_c[j].start = 
	(exon->start <= tr->cds.end) ? tr->cds.end+1 : exon->start;
      utr3_c[j].end = exon->end;
      utr3_c[j].strand = STRAND_FWD;
      utr3_c[j].seqname = NULL;
      utr3_c[j].chr = tr->c.chr;

      j++;
    } else {
      if(exon->start >= tr->cds.start) {
	continue;
      }

      /* this exon is all or partially 3' UTR */
      utr3_c[j].start = exon->start;
      utr3_c[j].end = 
	(exon->end >= tr->cds.start) ? tr->cds.start-1 : exon->end;
      utr3_c[j].strand = STRAND_REV;
      utr3_c[j].seqname = NULL;
      utr3_c[j].chr = tr->c.chr;
      
      j++;
    }
  }

  return utr3_c;
}


/**
 * Returns an array of SeqCoords that represent the coordinates of the
 * introns of the provided transcript. If the provided transcript is
 * intronless NULL is returned instead of an array.  The returned
 * array should be freed when it is no longer needed using
 * seq_coord_array_free.  The number of intron coordinates in the
 * array is assigned to the value pointed to by the num_introns
 * argument.
 */
SeqCoord *transcript_introns(Transcript *tr, int *num_introns) {
  int i;
  SeqCoord *introns;

  if(tr->num_exons < 2) {
    *num_introns = 0;
    return NULL;
  } 

  *num_introns = tr->num_exons - 1;
  introns = g_new(SeqCoord, *num_introns);

  for(i = 0; i < *num_introns; i++) {
    if(tr->c.strand == STRAND_REV) {
      introns[i].start = tr->exons[i+1].end + 1;
      introns[i].end   = tr->exons[i].start - 1;
    } else {
      introns[i].start = tr->exons[i].end + 1;
      introns[i].end   = tr->exons[i+1].start - 1;
    }

    /* ugly hack: Some UCSC genes define 0bp introns for unknown
     * reasons.  We eliminate them here because they are an annoyance
     * for some analyses.
     */
    if(introns[i].start > introns[i].end) {
      i--;
      *num_introns -= 1;
      continue;
    }

    introns[i].strand = tr->c.strand;
    introns[i].chr = tr->c.chr;

    if(tr->c.seqname == NULL) {
      introns[i].seqname = NULL;
    } else {
      introns[i].seqname = g_strdup(tr->c.seqname);
    }
  }

  if(*num_introns == 0) {
    g_free(introns);
    return NULL;
  }


  return introns;
}
