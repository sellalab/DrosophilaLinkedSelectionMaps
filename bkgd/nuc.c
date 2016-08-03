#include <glib.h>
#include <stdio.h>
#include <string.h>

#include "nuc.h"
#include "seqfeat.h"

#include "util.h"


static const char NUC_SYMBOL[NUM_NUCS] = {'A', 'C', 'G', 'T', '-', 'N'};


/**
 * Converts a nucleotide ID into a character.
 */
char nuc_id_to_char(const unsigned char id) {
  if(id > NUM_NUCS) {
    g_error("Invalid nucleotide identifier %d", id);
  }
  return NUC_SYMBOL[id];
}



/**
 * Converts a nucleotide to a unique nucleotide integer identifier.
 */
unsigned char nuc_char_to_id(const char nuc) {
  switch(nuc) {
  case('A'): case('a'): return NUC_A;
  case('C'): case('c'): return NUC_C;
  case('T'): case('t'): return NUC_T;
  case('G'): case('g'): return NUC_G;
  case('.'): case('-'): case('*'): return NUC_GAP;
  }
  return NUC_N;
}



/**
 * Fills the provided buffer with a null-terminated string representation 
 * of the provided array of nucleotide arrays. The provided buffer
 * must be at least len+1 bytes long, and the provided nucleotide
 * id array must be at least len bytes long.
 *
 * If the provided buffer is NULL, a new buffer of length len+1 is
 * allocated and returned.
 */
char *nuc_ids_to_str(char *buf, const unsigned char *ids, const long len) {
  long i;
  
  if(buf == NULL) {
    buf = g_new(char, len+1);
  }

  for(i = 0; i < len; i++) {
    buf[i] = nuc_id_to_char(ids[i]);
  }
  buf[len] = '\0';
  
  return buf;
}



/**
 * Fills the provided buffer with nucleotide ids corresponding to the
 * string representation of DNA sequence provided. The returned array
 * is NOT null terminated. The len argument should correspond to the
 * length of the provided DNA character string. The provided buf must
 * be at least len bytes long. If the provided buf is NULL a new buffer
 * of length len is allocated and returned.
 *
 * The return value is just the ptr to the provided buf.
 */
unsigned char *nuc_str_to_ids(unsigned char *buf, const char *str, 
			      const long len) {
  long i;
  
  if(buf == NULL) {
    buf = g_new(unsigned char, len);
  }

  for(i = 0; i < len; i++) {
    buf[i] = nuc_char_to_id(str[i]);
  }
  return buf;
}


/**
 * Reads a nucleotide matrix (2D array) from a file. The values of the
 * matrix could represent conditional probabilities, scores, etc. The
 * first row of the matrix must contain tab-delimited
 * nucleotide column labels. Each subsequent row should contain a row
 * label followed by tab-delimited double values.
 * 
 * The returned array should be freed when it is no longer needed.
 */
gdouble **nuc_read_matrix(const char *filename) {
  gdouble **matrix;
  gint col_nuc_ids[NUM_NUCS];
  gint row_nuc_id;
  gint i, j;
  gint num_toks;
  char *line;
  char **tokens;
  FILE *fh;

  /* allocate matrix memory, initialize scores to 0.0 */
  matrix = g_new(gdouble *, NUM_NUCS);
  for(i = 0; i < NUM_NUCS; i++) {
    matrix[i] = g_new(gdouble, NUM_NUCS);
    for(j = 0; j < NUM_NUCS; j++) {
      matrix[i][j] = 0.0;
    }
  }

  fh = fopen(filename, "r");

  if(fh == NULL) {
    g_error("Could not open nucleotide matrix file '%s'", filename);
  }

  /* read the header line containing nucleotide col labels */
  line = util_fgets_line(fh);
  tokens = g_strsplit(line, "\t", -1);
  i = 0;
  num_toks = 0;
  while(tokens[i] != NULL) {
    if(num_toks > NUM_NUCS) {
      g_error("Too many matrix column labels, Max = %d", NUM_NUCS+1);
    }

    /* skip blank tokens from leading/trailing tab */
    if(tokens[i][0] == '\0') {
      i++;
      continue;
    }

    if(strlen(tokens[i]) > 1) {
      g_error("Nucleotide matrix column labels must be one character.");
    }

    col_nuc_ids[num_toks] = nuc_char_to_id(tokens[i][0]);
    num_toks++;
    i++;
  }
  g_strfreev(tokens);
  g_free(line);

  i = 0;
  while((line = util_fgets_line(fh)) != NULL) {
    if(i == NUM_NUCS) {
      g_error("Too many matrix rows, Max = %d", NUM_NUCS);      
    }

    tokens = g_strsplit(line, "\t", -1);

    /* first column should be nucleotide row label */    
    if(tokens[0] == NULL) {
      g_error("Matrix line should not be blank.");
    }
    row_nuc_id = nuc_char_to_id(tokens[0][0]);

    /* subsequent columns have actual values */
    j = 1;
    while(tokens[j] != NULL) {
      if(j > NUM_NUCS) {
	g_error("Too many matrix cols in row %d, Max = %d", i, NUM_NUCS);      
      }
      

      matrix[row_nuc_id][col_nuc_ids[j-1]] = strtod(tokens[j], NULL);
      j++;
    }
    
    g_strfreev(tokens);
    g_free(line);

    i++;
  }

  fclose(fh);
  
  return matrix;
}

/**
 * Creates a new unsigned-long NUM_NUCSxNUM_NUCS nucleotide matrix
 * with all counts initialized to 0. The matrix should be freed when
 * it is no longer needed.
 */
long **nuc_matrix_ul_new() {
  long **matrix;
  gint i;

  matrix = g_new(long *, NUM_NUCS);
  for(i = 0; i < NUM_NUCS; i++) {
    matrix[i] = g_new0(long, NUM_NUCS);
  }

  return matrix;
}


/**
 * Creates a new unsigned-long NUM_NUCSxNUM_NUCS double floating point
 * nucleotide matrix with all values initialized to 0.0. The matrix
 * should be freed when it is no longer needed.
 */
gdouble **nuc_matrix_dbl_new() {
  gdouble **matrix;
  gint i, j;

  matrix = g_new(gdouble *, NUM_NUCS);
  for(i = 0; i < NUM_NUCS; i++) {
    matrix[i] = g_new(gdouble, NUM_NUCS);
    for(j = 0; j < NUM_NUCS; j++) {
      matrix[i][j] = 0.0;
    }
  }

  return matrix;
}




/**
 * Sets the elements of an usigned-long NUM_NUCS x NUM_NUCS  nucleotide
 * count matrix to 0.
 */
void nuc_matrix_ul_set_zero(long **matrix) {
  gint i, j;

  for(i = 0; i < NUM_NUCS; i++) {
    for(j = 0; j < NUM_NUCS; j++) {
      matrix[i][j] = 0;
    }
  }
}



/**
 * Writes a representation of a (double floating point) nucleotide
 * matrix (2D array) to a provided filehandle.
 */
void nuc_write_matrix_dbl(FILE *fh, gdouble **matrix) {
  gshort i,j;

  for(i = 0; i < NUM_NUCS; i++) {
    fprintf(fh, "\t%c", nuc_id_to_char(i)); 
  }
  fprintf(fh, "\n");

  for(i = 0; i < NUM_NUCS; i++) {
    fprintf(fh, "%c", nuc_id_to_char(i));
    for(j = 0; j < NUM_NUCS; j++) {
      fprintf(fh, "\t%f", matrix[i][j]);
    }
    fprintf(fh, "\n");
  }
}


/**
 * Writes a representation of an (unsigned long integer) nucleotide
 * matrix (2D array) to a provided filehandle.
 */
void nuc_write_matrix_ul(FILE *fh, long **matrix) {
  gshort i,j;

  for(i = 0; i < NUM_NUCS; i++) {
    fprintf(fh, "\t%c", nuc_id_to_char(i)); 
  }
  fprintf(fh, "\n");

  for(i = 0; i < NUM_NUCS; i++) {
    fprintf(fh, "%c", nuc_id_to_char(i));
    for(j = 0; j < NUM_NUCS; j++) {
      fprintf(fh, "\t%lu", matrix[i][j]);
    }
    fprintf(fh, "\n");
  }
}
		      

/*
 * Returns an array representing a table of nucleotide
 * composition for the provided coordinates, optionally ignoring any
 * bases that are soft-masked (lower-cased). If the region is on the
 * reverse strand the reverse complement of the sequence is used
 * instead.
 *
 * The nucleotides of the comp_array correspond to nucleotides (as
 * defined by the NUC ids defined in nuc.h). If an array of existing
 * counts is provided, the counts in the existing array are
 * incremented with the counts from the provided region.  If,on the
 * other hand, the comp_array argument is NULL a new array is created
 * and returned with counts for the provided region (and the array
 * should be freed when it is no longer needed).
 *
 * An array of bitmasks can be provided which specify what type of
 * features should be ignored for the composition analysis.
 */
long *nuc_composition(Seq *seq, SeqCoord *coord, SeqMask *ignore_mask, 
		      unsigned char ignore_ids, long *comp_array) {
  long i;

  if(comp_array == NULL) {
    comp_array = g_new0(long, NUM_NUCS);
  }
 
  if(coord->strand == STRAND_REV) {
    for(i = coord->start-1; i < coord->end; i++) {
      /* ignore bases that overlap w/ flagged regions in seqmask */
      if(ignore_mask != NULL && (ignore_mask->mask[i] & ignore_ids)) {
	continue;
      }
      comp_array[nuc_comp(seq->sym[i])] += 1;
    }
  } else {
    for(i = coord->start-1; i < coord->end; i++) {
      /* ignore bases that overlap w/ flagged regions in seqmask */
      if(ignore_mask != NULL && (ignore_mask->mask[i] & ignore_ids)) {
	continue;
      }
      comp_array[seq->sym[i]] += 1;
    }
  }

  return comp_array;
}


/*
 * Counts nucleotides in an array of coords. For details
 * see the nuc_composition function.
 */
long *nuc_composition_coords(Seq *seq, SeqCoord *coords, long num_coords,
			     SeqMask *ignore_mask, unsigned char ignore_ids,
			     long *comp_array) {
  
  long i;

  if(num_coords == 0) {
    return NULL;
  }
  
  for(i = 0; i < num_coords; i++) {
    comp_array = nuc_composition(seq, &coords[i], ignore_mask, 
				 ignore_ids, comp_array);
  }
  
  return comp_array;
}
