#include <glib.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "util.h"

#include "bkgd_reader.h"


#define BKGD_READER_MAX_B 1000
#define BKGD_READER_MAX_LINE 1024


/**
 * Reads a B value line from the specified filehandle. Sets the length
 * in the int pointed to by the len argument, and returns the B value.
 */
static int bkgd_read_line(/*gz*/FILE *gzf, int *len) {
  char line[BKGD_READER_MAX_LINE];
  char *l_ptr;
  int b;
  
  /* skip over blank and comment lines */
  line[0] = '\0';
  while(line[0] == '\0' || line[0] == '#') {
	  if(fgets(line, BKGD_READER_MAX_LINE, gzf) == NULL) {
//	  if(gzgets(gzf, line, BKGD_READER_MAX_LINE) == NULL) {
      g_error("%s:%d: Could not read bkgd line. At EOF?", __FILE__, __LINE__);
    }    
  }

  /* read b-value and len from line */
  b = strtol(line, &l_ptr, 10);
  *len = strtol(l_ptr, NULL, 10);

  if(*len < 1) {
    g_error("%s:%d: len (%d) should not be < 1", __FILE__, __LINE__, *len);
  }

  if(b > BKGD_READER_MAX_B) {
    g_error("%s:%d: b (%d) should not be > %d", __FILE__, __LINE__, b,
	    BKGD_READER_MAX_B);
  }

  /* fprintf(stderr, "%s\n", line);*/

  return b;
}



/**
 * Creates a new bkgd reader that can be used to read through
 * B-values from a file sequentially without using much memory.
 */
BkgdReader *bkgd_reader_new(char *filename) {
  BkgdReader *br;
  size_t len;
  
  br = g_new(BkgdReader, 1);

  len = strlen(filename);
  
  br->gzf = /*gz*/fopen(filename, "rt"/*"rb"*/);
  if(br->gzf == NULL) {
    g_error("bkgd_reader_new: could not open file '%s'", filename);
  }  

  br->b = bkgd_read_line(br->gzf, &br->len);  
  br->pos = 0;

  return br;
}


/**
 * Frees memory allocated for BkgdReader and closes associated
 * filehandle.
 */
void bkgd_reader_free(BkgdReader *br) {
  /*gz*/fclose(br->gzf);
  g_free(br);
}


/**
 * Retrieves a B value at a specified position. Positions must be
 * requested in increasing order.
 */
int bkgd_reader_get_pos(BkgdReader *br, long pos) {
  if(br->pos > pos) {
    g_error("bkgd_reader_get_pos: position (%ld) must be >= "
	    "last requested position (%ld)", pos, br->pos);
  }

  while(br->pos < pos) {
    if(br->len == 0) {
      /* read next line from file */
      br->b = bkgd_read_line(br->gzf, &br->len);
    }

    br->pos += 1;
    br->len -= 1;    
  }

  return br->b;
}


