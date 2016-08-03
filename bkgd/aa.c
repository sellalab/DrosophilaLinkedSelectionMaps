
#include <glib.h>

#include "aa.h"
#include "seq.h"
#include "nuc.h"


static unsigned char ***__AA_CODON_TABLE = NULL;


static void init_codon_table() {
  if(__AA_CODON_TABLE == NULL) {
    __AA_CODON_TABLE = aa_codon_table();
  }
}


/**
 * Constructs a standard eukaryotic codon table that can be used to
 * lookup the amino acid translation of a codon. The first index
 * corresponds to the first nucleotide of the codon, the second index
 * corresponds to the second nucleotide of the codon and the third
 * index corresponds to the third nucleotide of the codon.  For
 * simplicity, Uracil is represented by the identifier NUC_T.
 * 
 * The table should be freed when it is no longer needed.
 */
unsigned char ***aa_codon_table() {
  unsigned char i, j;
  unsigned char ***codon_table;

  codon_table = g_new(unsigned char **, NUM_NUCS);
  for(i = 0; i < NUM_NUCS; i++) {
    codon_table[i] = g_new(unsigned char *, NUM_NUCS);
    for(j = 0; j < NUM_NUCS; j++) {
      codon_table[i][j] = g_new(unsigned char, NUM_NUCS);
    }
  }

  /* first base A */
  codon_table[NUC_A][NUC_A][NUC_A] = AA_K;
  codon_table[NUC_A][NUC_A][NUC_G] = AA_K;
  codon_table[NUC_A][NUC_A][NUC_C] = AA_N;
  codon_table[NUC_A][NUC_A][NUC_T] = AA_N;
 
  codon_table[NUC_A][NUC_G][NUC_A] = AA_R;
  codon_table[NUC_A][NUC_G][NUC_G] = AA_R;
  codon_table[NUC_A][NUC_G][NUC_C] = AA_S;
  codon_table[NUC_A][NUC_G][NUC_T] = AA_S;
 
  codon_table[NUC_A][NUC_C][NUC_A] = AA_T;
  codon_table[NUC_A][NUC_C][NUC_G] = AA_T;
  codon_table[NUC_A][NUC_C][NUC_C] = AA_T;
  codon_table[NUC_A][NUC_C][NUC_T] = AA_T;

  codon_table[NUC_A][NUC_T][NUC_A] = AA_I;
  codon_table[NUC_A][NUC_T][NUC_G] = AA_M;
  codon_table[NUC_A][NUC_T][NUC_C] = AA_I;
  codon_table[NUC_A][NUC_T][NUC_T] = AA_I;


  /* first base G */
  codon_table[NUC_G][NUC_A][NUC_A] = AA_E;
  codon_table[NUC_G][NUC_A][NUC_G] = AA_E;
  codon_table[NUC_G][NUC_A][NUC_C] = AA_D;
  codon_table[NUC_G][NUC_A][NUC_T] = AA_D;
 
  codon_table[NUC_G][NUC_G][NUC_A] = AA_G;
  codon_table[NUC_G][NUC_G][NUC_G] = AA_G;
  codon_table[NUC_G][NUC_G][NUC_C] = AA_G;
  codon_table[NUC_G][NUC_G][NUC_T] = AA_G;
 
  codon_table[NUC_G][NUC_C][NUC_A] = AA_A;
  codon_table[NUC_G][NUC_C][NUC_G] = AA_A;
  codon_table[NUC_G][NUC_C][NUC_C] = AA_A;
  codon_table[NUC_G][NUC_C][NUC_T] = AA_A;

  codon_table[NUC_G][NUC_T][NUC_A] = AA_V;
  codon_table[NUC_G][NUC_T][NUC_G] = AA_V;
  codon_table[NUC_G][NUC_T][NUC_C] = AA_V;
  codon_table[NUC_G][NUC_T][NUC_T] = AA_V;


  /* first base C */
  codon_table[NUC_C][NUC_A][NUC_A] = AA_Q;
  codon_table[NUC_C][NUC_A][NUC_G] = AA_Q;
  codon_table[NUC_C][NUC_A][NUC_C] = AA_H;
  codon_table[NUC_C][NUC_A][NUC_T] = AA_H;
 
  codon_table[NUC_C][NUC_G][NUC_A] = AA_R;
  codon_table[NUC_C][NUC_G][NUC_G] = AA_R;
  codon_table[NUC_C][NUC_G][NUC_C] = AA_R;
  codon_table[NUC_C][NUC_G][NUC_T] = AA_R;
 
  codon_table[NUC_C][NUC_C][NUC_A] = AA_P;
  codon_table[NUC_C][NUC_C][NUC_G] = AA_P;
  codon_table[NUC_C][NUC_C][NUC_C] = AA_P;
  codon_table[NUC_C][NUC_C][NUC_T] = AA_P;

  codon_table[NUC_C][NUC_T][NUC_A] = AA_L;
  codon_table[NUC_C][NUC_T][NUC_G] = AA_L;
  codon_table[NUC_C][NUC_T][NUC_C] = AA_L;
  codon_table[NUC_C][NUC_T][NUC_T] = AA_L;


  /* first base T */
  codon_table[NUC_T][NUC_A][NUC_A] = AA_STOP; /* ochre */
  codon_table[NUC_T][NUC_A][NUC_G] = AA_STOP; /* amber */
  codon_table[NUC_T][NUC_A][NUC_C] = AA_Y;
  codon_table[NUC_T][NUC_A][NUC_T] = AA_Y;
 
  codon_table[NUC_T][NUC_G][NUC_A] = AA_STOP; /* opal */
  codon_table[NUC_T][NUC_G][NUC_G] = AA_W;
  codon_table[NUC_T][NUC_G][NUC_C] = AA_C;
  codon_table[NUC_T][NUC_G][NUC_T] = AA_C;
 
  codon_table[NUC_T][NUC_C][NUC_A] = AA_S;
  codon_table[NUC_T][NUC_C][NUC_G] = AA_S;
  codon_table[NUC_T][NUC_C][NUC_C] = AA_S;
  codon_table[NUC_T][NUC_C][NUC_T] = AA_S;

  codon_table[NUC_T][NUC_T][NUC_A] = AA_L;
  codon_table[NUC_T][NUC_T][NUC_G] = AA_L;
  codon_table[NUC_T][NUC_T][NUC_C] = AA_F;
  codon_table[NUC_T][NUC_T][NUC_T] = AA_F;


  /* codons with unknown (N) nucleotides */
  codon_table[NUC_A][NUC_A][NUC_N] = AA_X;
  codon_table[NUC_A][NUC_G][NUC_N] = AA_X;
  codon_table[NUC_A][NUC_C][NUC_N] = AA_T;
  codon_table[NUC_A][NUC_T][NUC_N] = AA_X;

  codon_table[NUC_G][NUC_A][NUC_N] = AA_X;
  codon_table[NUC_G][NUC_G][NUC_N] = AA_G;
  codon_table[NUC_G][NUC_C][NUC_N] = AA_N;
  codon_table[NUC_G][NUC_T][NUC_N] = AA_V;

  codon_table[NUC_C][NUC_A][NUC_N] = AA_X;
  codon_table[NUC_C][NUC_G][NUC_N] = AA_R;
  codon_table[NUC_C][NUC_C][NUC_N] = AA_P;
  codon_table[NUC_C][NUC_T][NUC_N] = AA_L;

  codon_table[NUC_T][NUC_A][NUC_N] = AA_X;
  codon_table[NUC_T][NUC_G][NUC_N] = AA_X;
  codon_table[NUC_T][NUC_C][NUC_N] = AA_S;
  codon_table[NUC_T][NUC_T][NUC_N] = AA_X;
  
  for(i = 0; i < NUM_NUCS; i++) {
    for(j = 0; j < NUM_NUCS; j++) {
      codon_table[i][NUC_N][j] = AA_X;
      codon_table[NUC_N][i][j] = AA_X;
    }
  }

  return codon_table;
}


/**
 * Frees the memory allocated for a codon table.
 */
void aa_codon_table_free(unsigned char ***codon_table) {
  int i,j;

  for(i = 0; i < NUM_NUCS; i++) {

    for(j = 0; j < NUM_NUCS; j++) {
      g_free(codon_table[i][j]);
    }
    g_free(codon_table[i]);
  }

  g_free(codon_table);
}



/**
 * Converts an integer amino acid id to a single letter representation
 * of the amino acid.
 */
char aa_id_to_char(const unsigned char aa_id) {
  switch(aa_id) {
  case(AA_A):
    return 'A';
  case(AA_R):
    return 'R';
  case(AA_N):
    return 'N';
  case(AA_D):
    return 'D';
  case(AA_C):
    return 'C';
  case(AA_E):
    return 'E';
  case(AA_Q):
    return 'Q';
  case(AA_G):
    return 'G';
  case(AA_H):
    return 'H';
  case(AA_I):
    return 'I';
  case(AA_L):
    return 'L';
  case(AA_K):
    return 'K';
  case(AA_M):
    return 'M';
  case(AA_F):
    return 'F';
  case(AA_P):
    return 'P';
  case(AA_S):
    return 'S';
  case(AA_T):
    return 'T';
  case(AA_W):
    return 'W';
  case(AA_Y):
    return 'Y';
  case(AA_V):
    return 'V';
  case(AA_STOP):
    return '*';
  case(AA_X):
    return 'X';
  default:
    g_error("Unknown amino acid id '%d'", aa_id);	 
  }

  return '\0';
}


/**
 * Converts an amino acid character to an integer identifier.
 */
unsigned char aa_char_to_id(const char aa) {
  switch(aa) {
  case('A'):
  case('a'):
    return AA_A;
  case('R'):
  case('r'):
    return AA_R;
  case('N'):
  case('n'):
    return AA_N;
  case('D'):
  case('d'):
    return AA_D;
  case('C'):
  case('c'):
    return AA_C;
  case('E'):
  case('e'):
    return AA_E;
  case('Q'):
  case('q'):
    return AA_Q;
  case('G'):
  case('g'):
    return AA_G;
  case('H'):
  case('h'):
    return AA_H;
  case('I'):
  case('i'):
    return AA_I;
  case('L'):
  case('l'):
    return AA_L;
  case('K'):
  case('k'):
    return AA_K;
  case('M'):
  case('m'):
    return AA_M;
  case('F'):
  case('f'):
    return AA_F;
  case('P'):
  case('p'):
    return AA_P;
  case('S'):
  case('s'):
    return AA_S;
  case('T'):
  case('t'):
    return AA_T;
  case('W'):
  case('w'):
    return AA_W;
  case('Y'):
  case('y'):
    return AA_Y;
  case('V'):
  case('v'):
    return AA_V;
  case('*'):
    return AA_STOP;
  default:
    g_error("Unknown amino acid '%c'", aa);	 
  }

  return -1;
}



/**
 * Given a ptr to an array of three NUC ids (as defined in nuc.h)
 * representing the sequence of a 3 base-pair codon, returns the
 * integer id representing the amino acid translation of the codon.
 */
unsigned char aa_codon_to_id(const unsigned char *codon) {
  unsigned char aa;
  init_codon_table();

  aa = __AA_CODON_TABLE[codon[0]][codon[1]][codon[2]];

/*   fprintf(stderr, "looking up aa for codon: %c%c%c: %c\n", */
/* 	  nuc_id_to_char(codon[0]), nuc_id_to_char(codon[1]), */
/* 	  nuc_id_to_char(codon[2]), aa_id_to_char(aa)); */

  return aa;
}


/**
 * Returns single character amino-acid code that is the translation of
 * the codon (provided as array of three NUC ids).
 */
char aa_codon_to_aa(const unsigned char *codon) {
  return aa_id_to_char(aa_codon_to_id(codon));
}


/**
 * Returns the degeneracy of the provided position in the provided
 * codon. The provided codon should be an array of length 3 containing
 * NUC identifiers representing nucleotides (defined in nuc.h) The idx
 * (codon position) argument must be either 0,1, or 2.  The degeneracy
 * is how many of the 3 possible nucleotide changes at the given
 * position would result in a synonymous change (either 0, 1, 2, or
 * 3). 0 is returned if the codon contains an ambiguity nucleotide
 * (e.g. an N).
 */
int aa_codon_degeneracy(const unsigned char *codon, const int idx) {
  int change_count;
  int nuc_id, aa_id, new_aa_id, i;
  unsigned char new_codon[3];

  if((idx < 0) || (idx > 2)) {
    g_error("%s:%d: codon idx must be 0, 1 or 2", __FILE__, __LINE__);
  }

  nuc_id = nuc_char_to_id(codon[idx]);
  aa_id  = aa_codon_to_id(codon);
  
  /* this is an ambiguous amino acid, just return 0 */
  if(aa_id == AA_X) {
    return 0;
  }

  /* Try every other nucleotide at the specified position
   * Count when amino acid is different.
   */
  new_codon[0] = codon[0];
  new_codon[1] = codon[1];
  new_codon[2] = codon[2];
  change_count = 0;
  for(i = 0; i < NUM_REAL_NUCS; i++) {
    if(i == nuc_id) {
      continue;
    }
    new_codon[idx] = i;
    new_aa_id = aa_codon_to_id(new_codon);

    if(new_aa_id == AA_X) {
      /* too complicated to handle codons with ambiguity symbols */
      return 0;
    }

    if(aa_id != new_aa_id) {
      change_count++;
    }
  }

  return 3 - change_count;
}



/**
 * Given a nucleotide sequence, returns a sequence representing the
 * amino acid translation of the sequence.
 */
Seq *aa_translate_seq(const Seq *rna) {
  Seq *peptide;
  long rna_pos, pep_pos;

  if((rna->len % 3) != 0) {
    g_error("Cannot translate sequence that is not divisible by 3.");
  }

  peptide = g_new(Seq, 1);
  peptide->len = rna->len / 3;
  peptide->sym = g_new(unsigned char, peptide->len +  1);
  peptide->c.strand = STRAND_NONE;
  peptide->c.start = 1;
  peptide->c.end = peptide->len;
  peptide->c.seqname = NULL;
  peptide->c.chr = NULL;
  peptide->name = NULL;

  rna_pos = 0;
  for(pep_pos = 0; pep_pos < peptide->len; pep_pos++) {
    peptide->sym[pep_pos] = aa_codon_to_aa(&rna->sym[rna_pos]);
    rna_pos += 3;
  }

  peptide->sym[peptide->len] = '\0';
  
  return peptide;
}



/**
 * Utility method. Returns FALSE if any of the three NUC ids in the
 * provided codon array are NUC_N or NUC_GAP. Returns TRUE otherwise.
 */
int aa_is_codon(const unsigned char *codon) {
  gshort i;

  for(i = 0; i < 3; i++) {
    if(codon[i] == NUC_N || codon[i] == NUC_GAP) {
      return FALSE;
    }
  }

  return TRUE;
}

