#include <glib.h>
#include <string.h>

#include "transcript.h"
#include "gene.h"
#include "aa.h"

#include "util.h"

 
/**
 * Looks at the frequency of synonymous and non-synonymous
 * substitutions in the coding regions of the transcripts of the
 * provided gene and calculates N,S,dN,dS and dN/dS. Makes sure
 * that each coding base is only considered once.
 */
void gene_count_coding_substs(Gene *gene, Seq *q_seq, Seq *t_seq) {
  long i;
  GSList *cur_tr;
  Transcript *tr;
  int codon_pos, gene_pos, q_aa_id, sub_aa_id;
  Seq *q_cds, *t_cds;
  unsigned char *q_codon, *t_codon, sub_codon[3];
  long nonsyn_count, syn_count, ttl_sites, ttl_degeneracy, gene_len;
  double dN, dS, dNdS;
  int *already_examined;

  /* keep track of which bases have already been looked at
   * so we don't count them twice for different transcripts
   */
  gene_len = gene->c.end - gene->c.start + 1;
  already_examined = g_new(int, gene_len);
  for(i = 0; i < gene_len; i++) {
    already_examined[i] = FALSE;
  }

  /* move along coding sequences, count substitutions as synonymous or
   * non-synonymous
   */
  nonsyn_count = syn_count = 0;
  ttl_sites = 0;
  ttl_degeneracy = 0;

  cur_tr = gene->transcripts;
  while(cur_tr != NULL) {
    tr = (Transcript *)cur_tr->data;

    codon_pos = 0;

    q_cds = transcript_cds_seq(tr, q_seq);
    t_cds = transcript_cds_seq(tr, t_seq);

    q_codon = q_cds->sym;
    t_codon = t_cds->sym;

    for(i = 0; i < q_cds->len; i++) {
      if(codon_pos == 0) {
	/* skip over codon if it is not a simple substitution, i.e. if
	 * there is a gap in either sequence.
	 */
	if(!aa_is_codon(q_codon) || !aa_is_codon(t_codon)) {
	  i += 2;
	  continue;
	}
      }

      /* only examine each position once, despite overlapping transcripts */
      if(tr->c.strand == STRAND_REV) {
	gene_pos = tr->cds.end   - i - gene->c.start;
      } else {
	gene_pos = tr->cds.start + i - gene->c.start;
      }
      if(already_examined[gene_pos]) {
	continue;
      } else {
	already_examined[gene_pos] = TRUE;
      }

      /* Count the total number of unambiguously aligned sites and the
       * total degeneracy of sites (total number of possible synonymous
       * substitutions)
       */
      ttl_sites++;
      ttl_degeneracy += aa_codon_degeneracy(q_codon, codon_pos);

      if(q_cds->sym[i] != t_cds->sym[i]) {
	/* there is a substitution in the CDS,
	 * classify as synonymous or non-synonmous
	 */
	sub_codon[0] = q_codon[0];
	sub_codon[1] = q_codon[1];
	sub_codon[2] = q_codon[2];
	sub_codon[codon_pos] = t_codon[codon_pos];
      
	q_aa_id   = aa_codon_to_id(q_codon);
	sub_aa_id = aa_codon_to_id(sub_codon);

	if(q_aa_id != sub_aa_id) {
	  nonsyn_count++;
	} else {
	  syn_count++;
	}
      }

      /* update codon position */
      if(codon_pos == 2) {
	codon_pos = 0;
	q_codon = &q_cds->sym[i];
	t_codon = &t_cds->sym[i];
      } else {
	codon_pos++;
      }
    }

    seq_free(q_cds);
    seq_free(t_cds);

    cur_tr = g_slist_next(cur_tr);
  }

  if(syn_count == 0) {
    dS = 0.0;
    dN = 0.0;
    dNdS = 0.0;
  } else {
    dS = (double)syn_count / (double)ttl_degeneracy;
    dN = (double)nonsyn_count / (double)((3 * ttl_sites) - ttl_degeneracy);
    dNdS = dN/dS;
  }

  g_free(already_examined);

  fprintf(stderr, "%s: N=%lu, S=%lu, dN=%f, dS=%f, dN/dS=%f\n",
	  gene->name, nonsyn_count, syn_count, dN, dS, dNdS);

  return;
}




/*
 * Takes an array of transcripts and returns an array of Gene
 * structures each containing a set of transcripts that overlap. It's
 * currently assumed that all transcripts are on the same sequence.
 * The value pointed to by num_genes is the number of genes in the
 * returned array of genes. The gene array should be freed when it is
 * no longer needed.
 */
Gene *gene_group_transcripts(Transcript *trs, long num_tr, long *num_genes) {
  long gene_num, i;
  Gene *genes, *gene;
  gchar *new_name;

  if(num_tr < 1) {
    return NULL;
  }

  /* sort all transcripts by their strand and start coordinate */
  fprintf(stderr, "sorting transcripts\n");
  qsort(trs, num_tr, sizeof(Transcript), transcript_cmp);
    
  /* Allocate enough mem for one gene for every transcript. Reallocate
   * later to reflect smaller num actually needed.
   */
  genes = g_new(Gene, num_tr);

  gene_num = 0;

  gene = &genes[gene_num];
  gene->c.start   = trs[0].c.start;
  gene->c.end     = trs[0].c.end;
  gene->c.strand  = trs[0].c.strand;
  gene->c.chr = trs[0].c.chr;

  if(trs[0].c.seqname == NULL) {
    gene->c.seqname = NULL;
  } else {
    gene->c.seqname = g_strdup(trs[0].c.seqname);
  }
  
  if(trs[0].name == NULL) {
    gene->name = g_strdup("");
  } else {
    gene->name = g_strdup(trs[0].name);
  }

  gene->transcripts = g_slist_append(NULL, &trs[0]);

  for(i = 1; i < num_tr; i++) {
    /* check for overlap with current gene */
    if(seq_coord_ovlp(&trs[i].c, &gene->c, TRUE)) {
      /* extend coordinates of gene if needed */
      if(trs[i].c.end > gene->c.end) {
	gene->c.end = trs[i].c.end;
      }
      if(trs[i].c.start < gene->c.start) {
	gene->c.start = trs[i].c.start;
      }

      /* sanity check */
      if(trs[i].c.seqname != NULL && gene->c.seqname != NULL) {
	if(strcmp(trs[i].c.seqname, gene->c.seqname) != 0) {
	  
	  g_error("gene_group_transcripts: attempt to combine transcripts on "
		  "sequences '%s' and '%s'", trs[i].c.seqname, 
		  gene->c.seqname);
	}
      }

      /* append name of transcript */
      new_name = g_strconcat(gene->name, ",", trs[i].name, NULL);
      g_free(gene->name);
      gene->name = new_name;

      /* add transcript to list for this gene */
      gene->transcripts = g_slist_append(gene->transcripts, &trs[i]);
    } else {
      /* initialize new gene */
      gene_num++;
      gene = &genes[gene_num];
      gene->c.start   = trs[i].c.start;
      gene->c.end     = trs[i].c.end;
      gene->c.strand  = trs[i].c.strand;
      gene->c.chr     = trs[i].c.chr;

      if(trs[i].c.seqname == NULL) {
	gene->c.seqname = NULL;
      } else {
	gene->c.seqname = g_strdup(trs[i].c.seqname);
      }

      if(trs[i].name == NULL) {
	gene->name = g_strdup("");
      } else {
	gene->name = g_strdup(trs[i].name);
      }

      gene->transcripts = g_slist_append(NULL, &trs[i]);
    }
  }

  *num_genes = gene_num + 1;

  /* shrink array to just accomodate num merged transcripts */
  genes = g_realloc(genes, sizeof(Gene) * *num_genes);  

  return genes;
}


/**
 * Returns the longest transcript associated with a gene.
 */
Transcript *gene_longest_transcript(Gene *gene) {
  Transcript *tr, *longest_tr;
  GSList *cur;

  cur = gene->transcripts;

  longest_tr = NULL;
  while(cur != NULL) {
    tr = (Transcript *)cur->data;
    if(longest_tr == NULL || 
       seq_coord_len(&tr->c) > seq_coord_len(&longest_tr->c)){
      longest_tr = tr;
    }
    cur = g_slist_next(cur);
  }
  
  if(longest_tr == NULL) {
    g_error("Gene has no associated transcripts");
  }

  return longest_tr;
}




/**
 * Frees memory allocated for an array of genes. Frees memory
 * associated with referenced Transcripts only if the free_transcripts
 * flag is TRUE. Does not free attached Chromosome structures.
 */
void gene_array_free(Gene *genes, long num_genes, int free_transcripts) {
  long i;
  GSList *cur;

  if(num_genes == 0) {
    return;
  }

  for(i = 0; i < num_genes; i++) {
    g_free(genes[i].name);
    
    if(free_transcripts) {
      cur = genes[i].transcripts;
      while(cur != NULL) {
	transcript_free(cur->data);
	cur = g_slist_next(cur);
      }
    }
    g_slist_free(genes[i].transcripts);

    if(genes[i].c.seqname != NULL) {
      g_free(genes[i].c.seqname);
    }
  }

  g_free(genes);
}




/**
 * Writes out a flat file representation of an array of genes the the
 * provided file handle.
 */
void gene_write_flatfile(FILE *fh, Gene *gene_array, long num_genes) {
  long i;

  for(i = 0; i < num_genes; i++) {
    gene_write_flatfile_line(fh, &gene_array[i]);
  }
}



/**
 * Writes out a single gene line in a flat file representation to
 * the provided file handle.
 */
inline void gene_write_flatfile_line(FILE *fh, Gene *gene) {  
  GSList *cur;
  Transcript *tr;
  int first;
  char *seqname;

  if(gene->c.chr != NULL && gene->c.chr->name != NULL) {
    seqname = gene->c.chr->name;
  } else if(gene->c.seqname != NULL) {
    seqname = gene->c.seqname;
  } else {
    seqname = "";
  }

  fprintf(fh, "GENE\t%ld\t%s\t%s\t%lu\t%lu\t%d\t", gene->id, gene->name, 
	  seqname, gene->c.start, gene->c.end,  gene->c.strand);


  /* write out names of associated transcripts */
  cur = gene->transcripts;
  first = TRUE;
  while(cur != NULL) {
    tr = (Transcript *)cur->data;
    if(first) {
      first = FALSE;
    } else {
      fprintf(fh, ",");
    }


    fprintf(fh, "%ld", tr->id);
    cur = g_slist_next(cur);
  }
  fprintf(fh,"\n");
}




/**
 * Helper function, places transcripts from an array into a hashtable
 * indexed on their identifiers.
 */
static GHashTable *gene_transcript_array_to_hash(Transcript *trs,
						 long num_trs) {
  long i;
  GHashTable *htab;  

  htab = g_hash_table_new(util_long_hash, util_long_equal);

  for(i = 0; i < num_trs; i++) {
    g_hash_table_insert(htab, &trs[i].id, &trs[i]);
  }

  return htab;
}


/**
 * Helper function, reads a string of transcript identifiers and adds
 * the associated transcripts to the provided gene.
 */
static inline void gene_read_flatfile_transcripts(Gene *gene, 
						  GHashTable *tr_htab,
						  gchar *tr_str) {
  guint i;
  gchar **toks;
  glong tr_id;
  Transcript *tr;

  gene->transcripts = NULL;
  
  /* parse transcript identifiers from comma-delimited string */
  i = 0;
  toks = g_strsplit(tr_str, ",", 0);
  while(toks[i] != NULL) {
    tr_id = strtol(toks[i], NULL, 10);

    /* use identifier to lookup transcript from table */
    tr = g_hash_table_lookup(tr_htab, &tr_id);
    if(tr == NULL) {
      g_error("gene_read_flatfile_transcripts: "
	      "Gene %ld references unknown transcript %ld", gene->id, tr_id);
    }

    gene->transcripts = g_slist_append(gene->transcripts, tr);

    i++;
  }

  g_strfreev(toks);
}


/**
 * Reads a flatfile representation of genes into an array and returns
 * the array. An array of transcripts must be provided to this
 * function so that the transcripts can be associated with the
 * appropriate genes.  The number of genes read from the file is
 * assigned to the variable pointed to be the num_read argument.
 */
Gene *gene_read_flatfile(gchar *filename, Transcript *trs, long num_trs,
			 long *num_read) {
  gchar *line;
  gchar **toks;
  int num_toks, prefix_len;
  Gene *genes;
  GHashTable *tr_htab;
  long num_lines;
  long i;

  FILE *fh;

  fh = fopen(filename, "r");
  if(fh == NULL) {
    g_error("gene_read_flatfile: could not read from file '%s'", filename);
  }

  /* convert array of transcripts to hash table for easy lookup */
  tr_htab = gene_transcript_array_to_hash(trs, num_trs);


  /* count number of GENE lines in file  */
  num_lines = util_fcount_lines_match(fh, GENE_LINE_PREFIX);
  genes = g_new(Gene, num_lines);

  prefix_len = strlen(GENE_LINE_PREFIX);

  i = 0;
  while(i < num_lines) {
    if((line = util_fgets_line(fh)) == NULL) {
      g_error("gene_read_flatfile: Expected %lu %s lines, "
	      "got %lu", num_lines, GENE_LINE_PREFIX, i);
    }
    if(strncmp(GENE_LINE_PREFIX,line, prefix_len) != 0) {
      continue;
    }

    toks = g_strsplit(line, "\t", 8);

    num_toks = 0;
    while(toks[num_toks] != NULL) {
      num_toks++;
    }
    
    if(num_toks < 8) {
      g_error("gene_read_flatfile: Expected %d toks per gene line, got %d", 
	      8, num_toks);
    }


    genes[i].id      = strtol(toks[1], NULL, 10);
    genes[i].name    = g_strdup(toks[2]);
    genes[i].c.seqname = g_strdup(toks[3]);
    genes[i].c.chr = NULL;

    genes[i].c.start   = strtol(toks[4], NULL, 10);
    genes[i].c.end     = strtol(toks[5], NULL, 10);
    genes[i].c.strand  = strtol(toks[6], NULL, 10);

    gene_read_flatfile_transcripts(&genes[i], tr_htab, toks[7]);

    g_strfreev(toks);
    g_free(line);

    i++;
  }

  g_hash_table_destroy(tr_htab);

  *num_read = i;

  fclose(fh);

  return genes;
}


/*
 * Comparison function of genes used for sorting.  Uses the
 * seq_coord_cmp function to do a coordinate comparison.
 */
int gene_cmp(const void *p1,const void *p2) {
  return seq_coord_cmp(&((Gene *)p1)->c, &((Gene *)p2)->c);
}


/*
 * Comparison function of genes used for sorting.  Uses the
 * seq_coord_cmp function to do a coordinate comparison.
 */
int gene_cmp_nostrand(const void *p1, const void *p2) {
  return seq_coord_cmp_nostrand(&((Gene *)p1)->c, &((Gene *)p2)->c);
}


