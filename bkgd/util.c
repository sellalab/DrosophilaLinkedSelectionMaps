
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <ctype.h>
#include <zlib.h>

#include "util.h"




/**
 * Returns TRUE if the provided filename ends with ".gz", returns
 * FALSE otherwise.
 */
int util_has_gz_ext(char *filename) {
  size_t len;

  len = strlen(filename);

  if((len > 3) && (strncmp(&filename[len-3], ".gz",3)==0)) {
    return TRUE;
  }
  return FALSE;
}




/**
 * Helper function, reads an entire file into memory and returns the
 * result as a null-terminated string. The returned string should be
 * freed when it is no longer needed.
 */
char *util_read_entire_file(char *filename) {
  FILE *fh;
  struct stat fs;
  char *buf;

  fh = fopen(filename, "r");

  if(fh == NULL) {
    g_error("%s:%d: Could not open tree file '%s'", __FILE__, 
	    __LINE__, filename);
  }

  /* determine file size */
  fstat(fileno(fh), &fs);

  /* read entire file at once */
  buf = g_new(char, fs.st_size + 1);
  if(fread(buf, fs.st_size, 1, fh) == 0) {
    g_error("%s:%d: Could not read entire file '%s'", __FILE__,
	    __LINE__, filename);
  }

  buf[fs.st_size] = '\0';

  return buf;
}



/**
 * Counts the number of newline characters in the file pointed to by
 * filehandle. The filehandle is rewound to the beginning of the file
 * before and after the count of newlines.
 */
long util_fcount_lines(FILE *fh) {
  char buf[UTIL_FGETS_BUF_SZ];
  long line_count;
  int len;

  if(fseek(fh, 0L, SEEK_SET) != 0) {
    g_error("%s:%d: could not rewind filehandle", __FILE__, __LINE__);
  }

  line_count = 0;
  while(fgets(buf, UTIL_FGETS_BUF_SZ, fh) != NULL) {
    len = strlen(buf);
    if(buf[len-1] == '\n') {
      line_count++;
    }
  }

  if(fseek(fh, 0L, SEEK_SET) != 0) {
    g_error("%s:%d: could not rewind filehandle", __FILE__, __LINE__);
  }

  return line_count;
}



/**
 * Counts the number of lines in a gzipped fileThe gzFile is rewound
 * to the beginning of the file before and after the count of
 * newlines.
 */
long util_gzcount_lines(gzFile *gzf) {
  char buf[UTIL_FGETS_BUF_SZ];
  long line_count;
  size_t len;

  if(gzseek(gzf, 0L, SEEK_SET) != 0) {
    g_error("%s:%d: could not rewind filehandle", __FILE__, __LINE__);
  }

  line_count = 0;

  while(gzgets(gzf, buf, UTIL_FGETS_BUF_SZ) != NULL) {
    len = strlen(buf);
    if(buf[len-1] == '\n') {
      line_count++;
    }
  }

  if(gzseek(gzf, 0L, SEEK_SET) != 0) {
    g_error("%s:%d: could not rewind filehandle", __FILE__, __LINE__);
  }

  return line_count;
}




/**
 * Counts the number of lines in a file that begin with a provided
 * string. The filehandle is rewound to the beginning of the file
 * before and after the count of newlines.
 */
long util_fcount_lines_match(FILE *fh, const char *starts_with) {
  char buf[UTIL_FGETS_BUF_SZ];
  long line_count;
  gboolean started_with, at_line_start;
  gint match_len, len;

  if(fseek(fh, 0L, SEEK_SET) != 0) {
    g_error("%s:%d: could not rewind filehandle", __FILE__, __LINE__);
  }

  match_len = strlen(starts_with);
  if(match_len == 0) {
    return 0;
  }

  if(match_len > UTIL_FGETS_BUF_SZ) {
    g_error("%s:%d: length of string to match must be"
	    "<= %d bytes", __FILE__, __LINE__, UTIL_FGETS_BUF_SZ);
  }

  line_count = 0;
  started_with  = FALSE;
  at_line_start = TRUE;

  while(fgets(buf, UTIL_FGETS_BUF_SZ, fh) != NULL) {
    if(at_line_start) {
      /* we are at the beginning of a line, does it match the string? */
      if(strncmp(starts_with, buf, match_len) == 0) {
	started_with = TRUE;
      }
    }

    len = strlen(buf);
    if(buf[len-1] == '\n') {
      /* we are at the end of the line, did the beginning match? */
      if(started_with) {
	line_count++;
	started_with = FALSE;
      }
      at_line_start = TRUE;
    } else {
      at_line_start = FALSE;
    }
  }

  if(fseek(fh, 0L, SEEK_SET) != 0) {
    g_error("%s:%d: could not rewind filehandle", __FILE__, __LINE__);
  }

  return line_count;
}



/**
 * Counts the number of lines in a gzipped file that begin with a
 * provided string. The gzFile is rewound to the beginning of the
 * file before and after the count of newlines.
 */
long util_gzcount_lines_match(gzFile *gzf, const char *starts_with) {
  char buf[UTIL_FGETS_BUF_SZ];
  long line_count;
  int started_with, at_line_start;
  int match_len, len;

  if(gzseek(gzf, 0L, SEEK_SET) != 0) {
    g_error("%s:%d: could not rewind filehandle", __FILE__, __LINE__);
  }

  match_len = strlen(starts_with);
  if(match_len == 0) {
    return 0;
  }

  if(match_len > UTIL_FGETS_BUF_SZ) {
    g_error("%s:%d: length of string to match must be"
	    "<= %d bytes", __FILE__, __LINE__, UTIL_FGETS_BUF_SZ);
  }

  line_count = 0;
  started_with  = FALSE;
  at_line_start = TRUE;

  while(gzgets(gzf, buf, UTIL_FGETS_BUF_SZ) != NULL) {
    if(at_line_start) {
      /* we are at the beginning of a line, does it match the string? */
      if(strncmp(starts_with, buf, match_len) == 0) {
	started_with = TRUE;
      }
    }

    len = strlen(buf);
    if(buf[len-1] == '\n') {
      /* we are at the end of the line, did the beginning match? */
      if(started_with) {
	line_count++;
	started_with = FALSE;
      }
      at_line_start = TRUE;
    } else {
      at_line_start = FALSE;
    }
  }

  if(gzseek(gzf, 0L, SEEK_SET) != 0) {
    g_error("%s:%d: could not rewind filehandle", __FILE__, __LINE__);
  }

  return line_count;
}




/**
 * same as util_fgets_line but works with gz-compressed files.
 */
char *util_gzgets_line(gzFile *gzf) {
  char buf[UTIL_FGETS_BUF_SZ];
  char *line, *old_line;
  long len, ttl_len;

  if(gzgets(gzf, buf, UTIL_FGETS_BUF_SZ) == NULL) {
    return NULL;
  }

  len = strlen(buf);
  ttl_len = len;
  line = g_new(char, len+1);
  strcpy(line, buf);

  /* keep extending the line until a newline or EOF is reached */
  while(buf[len-1] != '\n') {
    if(gzgets(gzf, buf, UTIL_FGETS_BUF_SZ) == NULL) {
      return line;
    }

    /* extend the line */
    len = strlen(buf);
    old_line = line;

    line = g_strconcat(old_line, buf, NULL);
    ttl_len += len;

    g_free(old_line);
  }

  /* remove trailing newline char */
  if(line[ttl_len-1] == '\n') {
    line[ttl_len-1] = '\0';
  }

  return line;
}


/**
 * Like fgets, but reads an entire line, without the
 * need to provide a buffer of a fixed size.
 * The returned line should be freed when it is no
 * longer needed.
 * The newline character is replaced by a '\0'
 */
char *util_fgets_line(FILE *fh) {
  char buf[UTIL_FGETS_BUF_SZ];
  char *line, *old_line;
  long len, ttl_len;

  if(fgets(buf, UTIL_FGETS_BUF_SZ, fh) == NULL) {
    return NULL;
  }

  len = strlen(buf);
  ttl_len = len;
  line = g_new(char, len+1);
  strcpy(line, buf);

  /* keep extending the line until a newline or EOF is reached */
  while(buf[len-1] != '\n') {
    if(fgets(buf, UTIL_FGETS_BUF_SZ, fh) == NULL) {
      return line;
    }

    /* extend the line */
    len = strlen(buf);
    old_line = line;

    line = g_strconcat(old_line, buf, NULL);
    ttl_len += len;

    g_free(old_line);
  }

  /* remove trailing newline char */
  if(line[ttl_len-1] == '\n') {
    line[ttl_len-1] = '\0';
  }

  return line;
}


/**
 * Converts a singly-linked list structure into an array and returns
 * the array. Does not free the singly-linked list or the elements in
 * the single linked list (all of which shallow copied).  The array
 * should be freed when no longer needed. The value pointed to by
 * num_elem is set to the number of elements in the array.
 */
void *util_slist_to_array(GSList *slist, size_t sz, long *num_elem) {
  GSList *cur;
  long count;
  size_t mempos;
  char *array;

  /* count the number of elements in the list */
  count = 0;
  cur = slist;
  while(cur != NULL) {
    count++;
    cur = g_slist_next(cur);
  }

  /* allocate memory for the array */
  array = g_malloc(sz * count);

  /* copy the elements from the linked list */
  cur = slist;
  mempos = 0;
  while(cur != NULL) {
    memcpy(&array[mempos], cur->data, sz);
    mempos += sz;
    cur = g_slist_next(cur);
  }

  *num_elem = count;

  return (void *)array;
}


/**
 * Does an in-place uppercase of all characters in a string
 */
void util_str_uc(char *str) {
  size_t i;

  i = 0;
  while(str[i] != '\0') {
    str[i] = toupper(str[i]);
    i++;
  }

  return;
}



/**
 * Does an in-place lowercase of all characters in a string
 */
void util_str_lc(char *str) {
  size_t i;

  i = 0;
  while(str[i] != '\0') {
    str[i] = tolower(str[i]);
    i++;
  }

  return;
}


/*
 * Does an in-place replacement of one character by another
 * in a string.
 */
void util_str_replace(char *str, char from, char to) {
  size_t i;

  i = 0;
  while(str[i] != '\0') {
    if(str[i] == from) {
      str[i] = to;
    }
    i++;
  }
}


/*
 * Does an in-place removal of all instances of
 * a specified character in a string.
 */
void util_str_remove_char(char *str, const char c) {
  size_t i,j;

  i = j = 0;
  while(str[i] != '\0') {
    if(str[i] != c) {
      str[j] = str[i];
      j++;
    }
    i++;
  }
  str[j] = '\0';
}






/**
 * Converts a long integer to a string with commas delimited
 * the hundreds thousands, millions etc. The returned string should be
 * freed when it is no longer needed.
 */
char *util_long_to_comma_str(const long x) {
  char *str, *comma_str;
  gint len,i,j,k;
  gint max_commas;

  str = g_new(char, UTIL_STR_BUF_SZ);
  max_commas = (UTIL_STR_BUF_SZ / 3);

// FIXFIX
  len = sprintf(str, "%lu", x);
//  len = snprintf(str, UTIL_STR_BUF_SZ-max_commas, "%lu", x);

  if(len == 0) {
    g_error("%s:%d: could not write number %lu to string", 
	    __FILE__, __LINE__, x);
  }
  
  max_commas = len / 3;
  comma_str = g_new(char, len + max_commas + 1);

  /* copy the string in reverse, adding commas every three */
  j = 0;
  k = 0;
  for(i = len-1; i >= 0; i--) {
    comma_str[j] = str[i];
    k++;
    
    if(i > 0 && ((k %3) == 0)) {
      /* insert comma */
      j++;
      comma_str[j] = ',';
    }
    j++;
  }

  comma_str[j] = '\0';
  g_free(str);

  util_breverse(comma_str, j);

  return comma_str;
}



/**
 * Reverses the contents of a provided array. Arguments
 * are similar to those for qsort():
 *   base - ptr to the beginning of the array
 *   nmemb - number of elements in array
 *   size - size of each element in bytes
 */
void util_reverse(void *base, const size_t nmemb, const size_t size) {
  void *tmp;
  unsigned char *fwd_ptr, *rev_ptr;

  /* allocate mem for holding swap */
  tmp = malloc(size);
  if(tmp == NULL) {
    g_error("%s:%d: could not allocate mem", __FILE__, __LINE__);
  }
  
  /* initiailize pointers so that they point to beginning and end of array */
  fwd_ptr = (unsigned char *)base;
  rev_ptr = &((unsigned char *)base)[nmemb*size - size];

  /* move ptrs towards each other from either end of array */
  while(fwd_ptr < rev_ptr) {
    /* swap memory at each pointer */
    tmp = memcpy(tmp, fwd_ptr, size);
    memcpy(fwd_ptr, rev_ptr, size);
    memcpy(rev_ptr, tmp, size);

    /* advance ptrs */
    fwd_ptr += size;
    rev_ptr -= size;
  }

  free(tmp);
}


/**
 * Reverses a byte array of a given length
 */
void util_breverse(void *base, const size_t nmemb) {
  util_reverse(base, nmemb, sizeof(unsigned char));
}

/*
 * Does an in-place reversal of the order of the characters within a
 * string. This function calls strlen to determine the length of the
 * string. If this is already known, or this function will be called
 * repeatedly it is more efficient to call util_breverse instead.
 *
 */
void util_str_reverse(char *str) {
  util_breverse(str, strlen(str));
}


/*
 * Performs an in-place removal of whitespace characters from the
 * provided string.
 */
void util_str_remove_whitespace(char *str) {
  long i,j;

  i = j = 0;
  while(str[i] != '\0') {
    if(!isspace((unsigned char)str[i])) {
      str[j] = str[i];
      j++;
    }
    i++;
  }

  str[j] = '\0';
}


/*
 * Performs in-place removal of leading whitespace
 * characters from the provided string.
 */
void util_str_lstrip(char *str) {
  long leading_ws, len;
  
  leading_ws = 0;
  while(str[leading_ws] != '\0' && isspace((unsigned char)str[leading_ws])) {
    leading_ws++;
  }

  if(leading_ws > 0) {
    len = strlen(str);
    memmove(str, &str[leading_ws], len-leading_ws+1);
  }
}


/*
 * Performs in-place removal of trailing whitespace
 * characters from the provided string.
 */
void util_str_rstrip(char *str) {
  long len,i;

  len = strlen(str);
  
  i = len-1;
  while(i >= 0 && isspace((unsigned char)str[i])) {
    str[i] = '\0';
    i--;
  }
}


/**
 * Performs in-place removal of leading and trailing whitespace
 * from provided string.
 */
void util_str_strip(char *str) {
  util_str_lstrip(str);
  util_str_rstrip(str);
}


/**
 * Returns TRUE if the start of str exactly matches the characters in
 * start
 */
int util_str_starts_with(const char *str, const char *start) {
  size_t i;

  i = 0;

  while(start[i] != '\0') {
    if(str[i] != start[i]) {
      return FALSE;
    }
    i++;
  }

  return TRUE;
}


/**
 * Returns FALSE if the end of str exactly matches the characters in
 * end
 */
int util_str_ends_with(const char *str, const char *end) {
  size_t i, j;

  i = strlen(str);
  j = strlen(end);

  if(j > i) {
    /* string to match is longer than provided string */
    return FALSE;
  }

  while(j > 0) {
    if(str[i-1] != end[j-1]) {
      /* character does not match */
      return FALSE;
    }
    i--;
    j--;
  }

  return TRUE;
}



void __util_hash_table_keys(gpointer key, gpointer val, gpointer user_data) {
  char **keys;
  int *idx_ptr;
  void **data;

  data = user_data;

  keys = (char **)data[0];
  idx_ptr = (int *)data[1];

  keys[*idx_ptr] = g_strdup(key);

  *idx_ptr += 1;

  return;
}


/**
 * Returns an array of the keys for the provided hash table. The
 * strings are all copies (i.e. do not point to the actual keys within
 * the hash table data structureq), and should be freed when they are
 * no longer needed.
 * 
 */
char **util_hash_table_keys(GHashTable *hash_tab, int *n_keys) {
  char **keys;
  int i;
  void *data[2];

  *n_keys = g_hash_table_size(hash_tab);
  if(*n_keys == 0) {
    return NULL;
  } 

  keys = g_new(char *, *n_keys);

  i = 0;
  data[0] = keys;
  data[1] = &i;

  g_hash_table_foreach(hash_tab, __util_hash_table_keys, data);

  return keys;
}


void __util_hash_table_free(gpointer key, gpointer val, gpointer user_data) {
  g_free(key);
  g_free(val);
}


int __util_hash_table_flush(gpointer key, gpointer val, gpointer hash_tab) {
  g_free(key);
  g_free(val);
  return TRUE;
}


/**
 * Frees all allocated keys and values in a hash table. Assumes that
 * it is ok to deallocate key and value memory using free().
 */
void util_hash_table_flush(GHashTable *hash_tab) {
  g_hash_table_foreach_remove(hash_tab, __util_hash_table_flush, hash_tab);
}


/**
 * Frees all allocated keys and values in a hash table and then frees
 * the hash table itself. Only safe to use if it is ok to free all of
 * the keys and values.
 */
void util_hash_table_free(GHashTable *hash_tab) {
  g_hash_table_foreach(hash_tab, __util_hash_table_free, NULL);
  g_hash_table_destroy(hash_tab);
}





/**
 * Takes a pointer to a long int and returns the long value cast as
 * an unsigned int so that it can be used as a hash value in a
 * GHashTable.
 */
unsigned int util_long_hash(gconstpointer v) {
  long *l;
  l = (long *)v;
  return (guint)*l;
}

/**
 * Takes pointers to two long ints and returns TRUE if the values they
 * point to are equal.  Used for comparison of long int hash keys in a
 * GHashTable
 */
int util_long_equal(gconstpointer v1, gconstpointer v2) {
  return *(long *)v1 == *(long *)v2;
}


/**
 * Comparison function used to compare two doubles Useful for sorting.
 * Returns -1 if x < y, returns 0 if x == y, returns 1 if x > y
 */
int util_dbl_cmp(const void *x, const void *y) {
  if(*(double *)x < *(double *)y) {
    return -1;
  }
  if(*(double *)x > *(double *)y) {
    return 1;
  }
  return 0;
}
