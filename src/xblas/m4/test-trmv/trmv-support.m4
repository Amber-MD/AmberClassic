dnl **********************************************************************
dnl * Generate tr_copy_row, and tr_commit_row                                *
dnl **********************************************************************
dnl
include(cblas.m4)dnl
dnl
dnl Usage: TR_COPY_ROW(T_typeltr)
dnl        copy a row to the vector y
dnl ----------------------------------------------------------------------
define(`TR_COPY_ROW_HEAD', 
  `void $1tr_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, dnl
       enum blas_trans_type trans, int n, const $1_array T, int ldt, dnl
       $1_array y, int row)')dnl
dnl
dnl
define(`TR_COPY_ROW', 
`TR_COPY_ROW_HEAD($1)
/*
 * Purpose
 * =======
 *
 * Copy a row from T to y
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of T; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether T is upper or lower
 *
 * trans        (input) blas_trans_type
 *              Whether T is no trans, trans, or conj trans 
 *
 * n            (input) int
 *              Dimension of AP and the length of vector y
 *
 * T            (input) $1_array
 *              The triangular matrix T
 *
 * ldt          (input) int
 *              Leading dimension of matrix T
 *
 * y            (output) $1_array
 *              The vector y
 *
 * row          (input) int
 *              The row to be copyied to y
 */
{
  int i, inc=1;
  DECLARE(tmp, $1_type)
  PTR_CAST(T, $1_type, `const')
  PTR_CAST(y, $1_type)
  INC_ADJUST(inc, $1_type)

  for (i = 0; i < n*inc; i += inc ) {
    SET_ZERO_VECTOR_ELEMENT(y_i, i, $1_type)
  }

  if ((order == blas_rowmajor && 
       trans == blas_no_trans && uplo == blas_upper) ||
      (order == blas_colmajor && 
       trans != blas_no_trans && uplo == blas_lower)) {
    for(i = row; i < n; i++) {
      GET_VECTOR_ELEMENT(tmp, T_i, (row*ldt+i)*inc, $1_type)
      SET_VECTOR_ELEMENT(y_i, i*inc, tmp, $1_type)
    }
  } else if ((order == blas_rowmajor && 
              trans == blas_no_trans && uplo == blas_lower) ||
             (order == blas_colmajor && 
              trans != blas_no_trans && uplo == blas_upper)) {
    for(i = 0; i <= row; i++) {
      GET_VECTOR_ELEMENT(tmp, T_i, (row*ldt+i)*inc, $1_type)
      SET_VECTOR_ELEMENT(y_i, i*inc, tmp, $1_type)
    }
  } else if ((order == blas_rowmajor && 
              trans != blas_no_trans && uplo == blas_lower) ||
             (order == blas_colmajor && 
              trans == blas_no_trans && uplo == blas_upper)) {
    for(i = row; i < n; i++) {
      GET_VECTOR_ELEMENT(tmp, T_i, (i*ldt+row)*inc, $1_type)
      SET_VECTOR_ELEMENT(y_i, i*inc, tmp, $1_type)
    }
  } else if ((order == blas_rowmajor && 
              trans != blas_no_trans && uplo == blas_upper) ||
             (order == blas_colmajor && 
              trans == blas_no_trans && uplo == blas_lower)) {
    for(i = 0; i <= row; i++) {
      GET_VECTOR_ELEMENT(tmp, T_i, (i*ldt+row)*inc, $1_type)
      SET_VECTOR_ELEMENT(y_i, i*inc, tmp, $1_type)
    }
  } 
}
') dnl
dnl
dnl
define(`PROTOTYPES', `dnl
TR_COPY_ROW_HEAD(s);
TR_COPY_ROW_HEAD(d);
TR_COPY_ROW_HEAD(c);
TR_COPY_ROW_HEAD(z);
')dnl
dnl
dnl
define(`SOURCE', `dnl
#include "blas_extended.h"
#include "blas_extended_test.h"

TR_COPY_ROW(s)
TR_COPY_ROW(d)
TR_COPY_ROW(c)
TR_COPY_ROW(z)
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
