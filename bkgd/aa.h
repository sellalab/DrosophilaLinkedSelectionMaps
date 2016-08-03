#ifndef __AA_H__
#define __AA_H__
#include <math.h>
#include <glib.h>
#include "seq.h"

enum amino_acid {AA_A=0,  /* alanine  */
		 AA_R,    /* arginine */
		 AA_N,    /* asparagine */
		 AA_D,    /* aspartic acid */
		 AA_C,    /* cysteine */
		 AA_E,    /* glutamic acid */
		 AA_Q,    /* glutamine */
		 AA_G,    /* glycine */
		 AA_H,    /* histidine */
		 AA_I,    /* isoleucine */
		 AA_L,    /* leucine */
		 AA_K,    /* lysine */
		 AA_M,    /* methionine */
		 AA_F,    /* phenylalanine */
		 AA_P,    /* proline */
		 AA_S,    /* serine */
		 AA_T,    /* threonine */
		 AA_W,    /* tryptophan */
		 AA_Y,    /* tyrosine */
		 AA_V,    /* valine */
		 AA_STOP, /* stop codon */
		 AA_X,    /* unknown */
		 NUM_AA};


unsigned char ***aa_codon_table(void);
void aa_codon_table_free(unsigned char ***codon_table);
char aa_id_to_char(unsigned char aa_id);
unsigned char aa_char_to_id(const char aa);
unsigned char aa_codon_to_id(const unsigned char *codon);
char aa_codon_to_aa(const unsigned char *codon);
Seq  *aa_translate_seq(const Seq *rna);
int aa_codon_degeneracy(const unsigned char *codon, int pos);
int aa_is_codon(const unsigned char *codon);


#endif
