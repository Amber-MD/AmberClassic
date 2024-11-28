dnl
dnl
include(cblas.m4)dnl
include(test-common.m4)dnl
dnl
dnl
dnl Usage: COPY_VECTOR($1)
dnl        COPY_VECTOR_NAME($1)
dnl        COPY_VECTOR_HEAD($1)
dnl        COPY_VECTOR_BODY($1)
dnl
define(`COPY_VECTOR_NAME', `$1copy_vector')dnl
dnl
dnl
define(`COPY_VECTOR_HEAD', 
  `void COPY_VECTOR_NAME($1) (const $1_array x, int n, int incx, dnl
                              $1_array y, int incy)')dnl
dnl
define(`COPY_VECTOR_BODY', 
  `PTR_CAST(x, $1_type, `const')
   PTR_CAST(y, $1_type)
   int i, xi, yi;

   INC_ADJUST(incx, $1_type)
   INC_ADJUST(incy, $1_type)
   xi = (incx >= 0) ? 0 : (1-n)*incx;
   yi = (incy >= 0) ? 0 : (1-n)*incy;

   for (i = 0; i < n; i++, xi += incx, yi += incy) {
     COPY_VECTOR_ELEMENT(y_i, yi, x_i, xi, $1_type)
   }
')dnl
dnl
define(`COPY_VECTOR', 
  `COPY_VECTOR_HEAD($1) {
   COPY_VECTOR_BODY($1)
  }')dnl
dnl
dnl
dnl
dnl
define(`PROTOTYPES', `
COPY_VECTOR_HEAD(s);
COPY_VECTOR_HEAD(d);
COPY_VECTOR_HEAD(c);
COPY_VECTOR_HEAD(z);
')dnl
dnl
dnl
define(`SOURCE', `dnl
#include <stdio.h>
#include "blas_extended.h"
#include "blas_extended_test.h"

COPY_VECTOR(s)
COPY_VECTOR(d)
COPY_VECTOR(c)
COPY_VECTOR(z)
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
