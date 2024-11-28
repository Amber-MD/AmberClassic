dnl
dnl
include(cblas.m4)dnl
include(test-common.m4)dnl
dnl
dnl
dnl Usage: PRINT_VECTOR($1)
dnl        PRINT_VECTOR_NAME($1)
dnl        PRINT_VECTOR_HEAD($1)
dnl        PRINT_VECTOR_BODY($1)
dnl
define(`PRINT_VECTOR_NAME', `$1print_vector')dnl
dnl
dnl
define(`PRINT_VECTOR_HEAD', 
  `void PRINT_VECTOR_NAME($1) (const $1_array x, int n, int inc, const char *name)')dnl
dnl
define(`PRINT_VECTOR_BODY', 
  `PTR_CAST(x, $1_type, `const')
   int i, xi;

   INC_ADJUST(inc, $1_type)
   xi = (inc >= 0) ? 0 : (1-n)*inc;

   if (name) { printf("%s = ", name); }
   printf("[\n");
   for (i = 0; i < n; i++, xi += inc) {
     printf("  ");
     PRINT_ARRAY_ELEM(x_i, xi, $1_type)
     printf("\n");
   }
   printf("];\n");
')dnl
dnl
dnl
define(`PRINT_VECTOR', 
  `PRINT_VECTOR_HEAD($1) {
   PRINT_VECTOR_BODY($1)
  }')dnl
dnl
dnl
define(`PROTOTYPES', `
PRINT_VECTOR_HEAD(s);
PRINT_VECTOR_HEAD(d);
PRINT_VECTOR_HEAD(c);
PRINT_VECTOR_HEAD(z);
')dnl
dnl
dnl
define(`SOURCE', `dnl
#include <stdio.h>
#include "blas_extended.h"
#include "blas_extended_test.h"

PRINT_VECTOR(s)
PRINT_VECTOR(d)
PRINT_VECTOR(c)
PRINT_VECTOR(z)
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
