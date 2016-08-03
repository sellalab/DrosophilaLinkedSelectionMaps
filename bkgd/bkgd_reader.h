#ifndef __BKGD_READER_H__
#define __BKGD_READER_H__

#include <math.h>
#include <stdio.h>
#include <zlib.h>


typedef struct {
	FILE* gzf; /* handle to file containing b values for chr */
//  gzFile *gzf; /* handle to file containing b values for chr */

  long pos; /* current position on chr */
  int b;    /* current b value */
  int len;  /* remaining bases with same b value ahead of this position */

} BkgdReader;



BkgdReader *bkgd_reader_new(char *filename);
void bkgd_reader_free(BkgdReader *br);
int bkgd_reader_get_pos(BkgdReader *br, long pos);


#endif
