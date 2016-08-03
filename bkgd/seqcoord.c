#include <glib.h>
#include <string.h>

#include "seqcoord.h"




/**
 * Copies the attributes of one sequence coordinate
 * into another
 */
void seq_coord_copy(const SeqCoord *src, SeqCoord *dst) {
  if(!memcpy(dst, src, sizeof(SeqCoord))) {
    g_error("%s:%d: sequence coordinate copy failed", __FILE__, __LINE__);
  }
  if(src->seqname) {
    dst->seqname = g_strdup(src->seqname);
  }

}


static int seq_coord_cmp_helper(SeqCoord *sc1, SeqCoord *sc2,
				 int cmp_strand, int cmp_start) {
  /* should consider comparing assembly too */
  int cmp_val;
  
  /* first need to order by sequence name */
  if(sc1->chr != NULL && sc2->chr != NULL) {
    cmp_val = strcmp(sc1->chr->name, sc2->chr->name);
    if(cmp_val != 0) {
      return cmp_val;
    }
  }

  if(sc1->seqname != NULL && sc2->seqname != NULL) {
    cmp_val = strcmp(sc1->seqname, sc2->seqname);
    if(cmp_val != 0) {
      return cmp_val;
    }
  }

  if(cmp_strand) {
    /* order with rev strand before fwd strands next*/
    if(sc1->strand < sc2->strand) {
      return -1;
    }
    if(sc1->strand > sc2->strand) {
      return 1;
    }
  }

  if(cmp_start) {
    /* sort by start position */
    if(sc1->start == sc2->start) {
      return 0;
    }
    if(sc1->start < sc2->start) {
      return -1;
    }
  } else {
    /* sort by end position */
    if(sc1->end == sc2->end) {
      return 0;
    }
    if(sc1->end < sc2->end) {
      return -1;
    }
  }

  return 1;  
}



/**
 * Comparison function for SeqCoords used for sorting. First
 * chromosome names are compared (if both names are non-NULL) and
 * coordinates are ordered alphanumerically by their seqname. Then
 * strands are compared and fwd strand coordinates are taken to be
 * "higher" than unknown strand coordinates, which are "higher" than
 * negative strand coordinates. If the seqnames and strands of the two
 * coordinates to compare are equal, the START of each of the
 * coordinates is compared.
 *
 * Currently the assembly version of the chromosomes are not taken
 * into account.
 */
int seq_coord_cmp(const void *p1, const void *p2) {
  SeqCoord *sc1, *sc2;

  sc1 = (SeqCoord *)p1;
  sc2 = (SeqCoord *)p2;

  return seq_coord_cmp_helper(sc1,sc2, TRUE, TRUE);
}



/**
 * Comparison function for SeqCoords used for sorting. First
 * chromosome names are compared (if both names are non-NULL) and
 * coordinates are ordered alphanumerically by their seqname. Then
 * strands are compared and fwd strand coordinates are taken to be
 * "higher" than unknown strand coordinates, which are "higher" than
 * negative strand coordinates. If the seqnames and strands of the two
 * coordinates to compare are equal, the END of each of the
 * coordinates is compared.
 */
int seq_coord_cmp_end(const void *p1, const void *p2) {
  SeqCoord *sc1, *sc2;

  sc1 = (SeqCoord *)p1;
  sc2 = (SeqCoord *)p2;

  return seq_coord_cmp_helper(sc1,sc2, TRUE, FALSE);
}



/**
 * Returns the cumulative length of all of the coords in
 * the provided array
 */
long seq_coord_array_len(SeqCoord *c, long num_coords) {
  long len,i;

  len = 0;
  for(i = 0; i < num_coords; i++) {
    len += seq_coord_len(&c[i]);
  }

  return len;
}


/**
 * Comparison function for SeqCoords used for sorting. First chr names
 * are compared (if both names are non-NULL) and coordinates are
 * ordered alphanumerically by their chromosome name. If the
 * chromosome names of the two coordinates are equal the start of each
 * of the coordinates is compared.
 *
 * This is the same as the seq_coord_cmp function but strand
 * information is not used.
 */
int seq_coord_cmp_nostrand(const void *p1, const void *p2) {
  SeqCoord *sc1, *sc2;

  sc1 = (SeqCoord *)p1;
  sc2 = (SeqCoord *)p2;

  return seq_coord_cmp_helper(sc1,sc2, FALSE, TRUE);  
}




void seq_coord_write(FILE *fh, SeqCoord *sc) {
  fprintf(fh, "%s:%ld-%ld(%c)", (sc->chr == NULL) ? "" : sc->chr->name,
	  sc->start, sc->end, strand_to_char(sc->strand));
}



/**
 * Returns a string representation of the sequence coordinate.  The
 * returned string should be freed when it is no longer needed.
 */
char *seq_coord_str(SeqCoord *sc) {
  char *c_str;
  GString *str;

  if(sc->chr != NULL) {
    str = g_string_new(sc->chr->name);
  }
  else if (sc->seqname != NULL) {
    str = g_string_new(sc->seqname);
  }
  else {
    str = g_string_new("");
  }

  g_string_sprintfa(str, ":%ld-%ld(%c)", sc->start, sc->end, 
		    strand_to_char(sc->strand));

  c_str = str->str;

  g_string_free(str, FALSE);

  return c_str;
}


/**
 * Returns true if the provided coordinates overlap, false otherwise.
 * If the cmp_strand argument is true, coordinates are only considered
 * overlapping if their strands match. If both coordinates have
 * non-NULL chromosomes, the chromosome name is used for comparison as
 * well.
 */
int seq_coord_ovlp(SeqCoord *sc1, SeqCoord *sc2, int cmp_strand) {
  if(sc1->chr != NULL && sc2->chr != NULL) {
    if(strcmp(sc1->chr->name, sc2->chr->name) != 0) {
      return FALSE;
    }
  } else if(sc1->seqname != NULL && sc2->seqname != NULL) {
    if(strcmp(sc1->seqname, sc2->seqname) != 0) {
      return FALSE;
    }
  }

  if(cmp_strand) {
    if(sc1->strand != sc2->strand) {
      return FALSE;
    }
  }

  if(sc1->start <= sc2->end && sc1->end >= sc2->start) {
    return TRUE;
  }

  return FALSE;
}



/**
 * Frees an array of SeqCoords. Does not free associated chromosomes.
 */
void seq_coord_array_free(SeqCoord *scs, long num) {
  int i;

  if(num > 0) {
    for(i = 0; i < num; i++) {
      if(scs[i].seqname != NULL) {
	g_free(scs[i].seqname);
      }
    }

    g_free(scs);
  }
}
