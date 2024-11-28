include(cblas.m4)dnl
dnl
dnl
dnl Usage: GBMV_PREPARE(ap_typeltr)
dnl        Creates the gbmv-prepare routine.
dnl
define(`GBMV_PREPARE_HEAD', 
  `void $1gbmv_prepare(enum blas_order_type order, dnl
       enum blas_trans_type trans, int m, int n, int kl, int ku, dnl
       $1_array AB, int lda, $1_array y, int row, int *nfix2, dnl
       int *nmix, int *ysize)')dnl
dnl
dnl
define(`GBMV_PREPARE',
`GBMV_PREPARE_HEAD($1)
{
  int ra, la, i, iy;
  int n_i;
  int inc = 1;
  PTR_CAST(y,$1_type)
  if ((order == blas_colmajor) && (trans == blas_no_trans)) {
    ra = ku;
    la = kl;
  } else if ((order == blas_colmajor) && (trans == blas_trans)) {
    ra = kl;
    la = ku;
  } else if ((order == blas_rowmajor) && (trans == blas_no_trans)) {
    ra = ku;
    la = kl;
  } else {    /* rowmajor and blas_trans */
    ra = kl;
    la = ku;
  }

  if (trans == blas_no_trans)
    n_i=n;
  else
    n_i=m;

  /* DETERMINE SIZE OF VECTOR TO BE SUBMITTED TO TESTGEN */  

  if (row + ra + 1< n_i) {
    *ysize = ra + row + 1;
  } else {
    *ysize = n_i;
  }

  /* SET NMIX AND NFIX */
  if (row == 0) {
    *nmix = 0;
    *nfix2 = 0;
  } else {
    if (row > la) {
      *nfix2 = row - la;
    } else {
      *nfix2 = 0;
    }
    if (*nfix2 > n_i) *nfix2 = n_i;
    if (row  <= n_i-ra-1) {
      *nmix = (*ysize) - (*nfix2) -1;
    } else {
      *nmix = (*ysize) - (*nfix2);
    }
  }
   
  /* SET ALL VALUES OF Y TO ZERO */
  INC_ADJUST(inc, $1_type)
  for (i = 0, iy = 0; i < *ysize; i++, iy += inc) {
    SET_ZERO_VECTOR_ELEMENT(y_i, iy, $1_type)
  }
}')
dnl
dnl
define(`GBMV_COMMIT_HEAD', 
  `void $1gbmv_commit(enum blas_order_type order, dnl
       enum blas_trans_type trans, int m, int n, int kl, int ku, dnl
       $1_array AB, int lda, $1_array y, int row)')dnl
dnl
dnl
define(`GBMV_COMMIT',
`GBMV_COMMIT_HEAD($1)
{ 
  int i,j;
  PTR_CAST(AB, $1_type)
  PTR_CAST(y, $1_type)
  DECLARE(ytemp, $1_type)
  int inc = 1;

  row++;
  INC_ADJUST(inc, $1_type)

  if (trans == blas_no_trans) {
    if (order == blas_colmajor) {
      for (j = MAX(1, row-kl); j <= MIN(n, row + ku); j++) {
        GET_VECTOR_ELEMENT(ytemp, y_i, (j-1)*inc, $1_type)
        SET_VECTOR_ELEMENT(AB_i, (ku+row-j + lda*(j-1))*inc, ytemp, $1_type)
      }
    } else {
      for (j = MAX(1, row-kl); j <= MIN(n, row + ku); j++) {
        GET_VECTOR_ELEMENT(ytemp, y_i, (j-1)*inc, $1_type)
        SET_VECTOR_ELEMENT(AB_i, ((row-1)*lda + kl+j-row)*inc, ytemp, $1_type)
      }
    }
  } else {
    if (order == blas_colmajor) {
      for (i = MAX(1,row-ku); i <= MIN(m, row + kl); i++){
        GET_VECTOR_ELEMENT(ytemp, y_i, (i-1)*inc, $1_type)
        SET_VECTOR_ELEMENT(AB_i, (ku+i-row + lda*(row-1))*inc, ytemp, $1_type)
      }
    } else {
      for (i = MAX(1,row-ku); i <= MIN(m, row + kl); i++){
        GET_VECTOR_ELEMENT(ytemp, y_i, (i-1)*inc, $1_type)
        SET_VECTOR_ELEMENT(AB_i, ((i-1)*lda + kl - i+row)*inc, ytemp, $1_type)
      }
    }
  }
}')dnl
dnl
dnl
define(`GBMV_COPY_HEAD', 
  `void $1gbmv_copy(enum blas_order_type order, dnl
       enum blas_trans_type trans, int m, int n, int kl, int ku, dnl
       const $1_array AB, int lda, $1_array y, int row)')dnl
dnl
dnl
define(`GBMV_COPY',
`GBMV_COPY_HEAD($1)
{       
  int i ,j, iy;
  int max_mn;
  PTR_CAST(AB, $1_type, `const')
  PTR_CAST(y, $1_type)
  DECLARE(ytemp, $1_type)
  int inc = 1;

  row++;
  INC_ADJUST(inc, $1_type)
  max_mn=MAX(m, n);

  for (i = 0, iy = 0; i < max_mn; i++, iy += inc) {
    SET_ZERO_VECTOR_ELEMENT(y_i, iy, $1_type)
  }

  if (trans == blas_no_trans) {
    if (order == blas_colmajor) {
      for (j = MAX(1, row-kl); j <= MIN(n, row + ku); j++) {
        GET_VECTOR_ELEMENT(ytemp, AB_i, (ku+row-j + lda*(j-1))*inc , $1_type)
        SET_VECTOR_ELEMENT(y_i, (j-1)*inc, ytemp, $1_type)
      }
    } else {
      for (j = MAX(1, row-kl); j <= MIN(n, row + ku); j++) {
        GET_VECTOR_ELEMENT(ytemp, AB_i, ((row-1)*lda + kl+j-row)*inc, $1_type)
        SET_VECTOR_ELEMENT(y_i, (j-1)*inc, ytemp, $1_type)
      }
    }
  } else {
    if (order == blas_colmajor) {
      for (i = MAX(1,row-ku); i <= MIN(m, row + kl); i++) {
        GET_VECTOR_ELEMENT(ytemp, AB_i, ((ku+i-row) + lda*(row-1))*inc , $1_type)
        SET_VECTOR_ELEMENT(y_i, (i-1)*inc, ytemp, $1_type)
      }
    } else {
      for (i = MAX(1,row-ku); i <= MIN(m, row + kl); i++) {
        GET_VECTOR_ELEMENT(ytemp, AB_i, ((i-1)*lda + kl -i+row)*inc, $1_type)
        SET_VECTOR_ELEMENT(y_i, (i-1)*inc, ytemp, $1_type)
      }
    }
  }
}')dnl
dnl
dnl
define(`PROTOTYPES', `dnl
FOREACH(`s, c, d, z', `
GBMV_PREPARE_HEAD(arg);')

FOREACH(`s, c, d, z', `
GBMV_COMMIT_HEAD(arg);')

FOREACH(`s, c, d, z', `
GBMV_COPY_HEAD(arg);')
')dnl
dnl
dnl
define(`SOURCE', `dnl
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"

FOREACH(`s, c, d, z', `
GBMV_PREPARE(arg)
GBMV_COMMIT(arg)
GBMV_COPY(arg)
')
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
