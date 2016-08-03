#ifndef __METHOD_H__
#define __METHOD_H__

#include <math.h>
#include <glib.h>

typedef struct {
  gint id;
  gchar *name;
  gchar *desc;
} Method;


void method_free(Method *method);

#endif
