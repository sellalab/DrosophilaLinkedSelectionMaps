#ifndef __UTIL_H__
#define __UTIL_H__

#include <math.h>
#include <glib.h>
#include <stdio.h>
#include <zlib.h>

#define UTIL_FGETS_BUF_SZ 16348
#define UTIL_STR_BUF_SZ 20

// PATCH ON PATCH - GUY
/*
// PATCH - EYAL
#define isnan(x) _isnan(x) 
#define isinf(x) (!_finite(x))
*/ 

int util_has_gz_ext(char *filename);

char *util_read_entire_file(char *filename);
long util_fcount_lines(FILE *fh);
long util_gzcount_lines(gzFile *gzf);
long util_fcount_lines_match(FILE *fh, const char *starts_with);
long util_gzcount_lines_match(gzFile *gzf, const char *starts_with);
char *util_fgets_line(FILE *fh);
char *util_gzgets_line(gzFile *fh);

char *util_long_to_comma_str(const long x);

void util_reverse(void *base, const size_t nmemb, const size_t size);
void util_breverse(void *base, const size_t nmemb);

void util_str_replace(char *str, char from, char to);
void util_str_reverse(char *str);
void util_str_uc(char *str);
void util_str_lc(char *str);
void util_str_remove_char(char *str, const char c);
void util_str_remove_whitespace(char *str);
void util_str_lstrip(char *str);
void util_str_rstrip(char *str);
void util_str_strip(char *str);
int util_str_starts_with(const char *str, const char *start);
int util_str_ends_with(const char *str, const char *end);


char **util_hash_table_keys(GHashTable *hash_table, int *n_keys);
void util_hash_table_flush(GHashTable *hash_table);
void util_hash_table_free(GHashTable *hash_table);

unsigned int util_long_hash(gconstpointer v);
int util_long_equal(gconstpointer v1, gconstpointer v2);
int util_dbl_cmp(const void *x, const void *y);

#endif
