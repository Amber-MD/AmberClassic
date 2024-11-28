dnl
dnl
include(cblas.m4)dnl
include(test-common.m4)dnl
dnl
dnl
dnl  Usage: 
dnl    SKEW_COMMIT_ROW($1)
dnl    SKEW_COMMIT_ROW_HEAD($1)
dnl    SKEW_COMMIT_ROW_NAME($1)
dnl    SKEW_COMMIT_ROW_PARAMS($1)
dnl    SKEW_COMMIT_ROW_COMMENT($1)
dnl    SKEW_COMMIT_ROW_BODY($1)
dnl
dnl    $1 is the type of the matrix.
dnl
define(`SKEW_COMMIT_ROW_NAME', `$1skew_commit_row')dnl
dnl
dnl
define(`SKEW_COMMIT_ROW_PARAMS', 
  `enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_side_type side, int n, $1_array a, int lda, 
   $1_array a_vec, int row')dnl
dnl
dnl
define(`SKEW_COMMIT_ROW_HEAD', 
  `void SKEW_COMMIT_ROW_NAME($1)(SKEW_COMMIT_ROW_PARAMS($1))')dnl
dnl
dnl
define(`SKEW_COMMIT_ROW_COMMENT', `
/*
 *  Copies the given vector into the given row of symmetric matrix a.
 */')dnl
dnl
dnl
define(`SKEW_COMMIT_ROW_BODY', `{
  int conj_flag;
  int i;
  int ai, incai1, incai2;
  int vi, incvi = 1;

  DECLARE(a_elem, $1_type)
  PTR_CAST(a, $1_type)
  PTR_CAST(a_vec, $1_type, `const')

  if ((side == blas_left_side  && uplo == blas_upper) ||
      (side == blas_right_side && uplo == blas_lower))
    conj_flag = 1;
  else
    conj_flag = 0;

  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai1 = 1;
    incai2 = lda;
    ai = lda * row;
  } else {
    incai1 = lda;
    incai2 = 1;
    ai = row;
  }

  INC_ADJUST(incai1, $1_type)
  INC_ADJUST(incai2, $1_type)
  INC_ADJUST(ai, $1_type)
  INC_ADJUST(incvi, $1_type)

  for (i = 0, vi = 0; i < row; i++, ai += incai1, vi += incvi) {
    GET_VECTOR_ELEMENT(a_elem, a_vec_i, vi, $1_type)
    if (conj_flag == 1) {
      NEGATE(a_elem, $1_type)
    }
    SET_VECTOR_ELEMENT(a_i, ai, a_elem, $1_type)
  }

  for (; i < n; i++, ai += incai2, vi += incvi) {
    GET_VECTOR_ELEMENT(a_elem, a_vec_i, vi, $1_type)
    if (conj_flag == 0) {
      NEGATE(a_elem, $1_type)
    }
    SET_VECTOR_ELEMENT(a_i, ai, a_elem, $1_type)
  }
}')dnl
dnl
dnl
define(`SKEW_COMMIT_ROW', 
  `SKEW_COMMIT_ROW_HEAD($1)
   SKEW_COMMIT_ROW_COMMENT($1)
   SKEW_COMMIT_ROW_BODY($1)')dnl
dnl
dnl
dnl
dnl
dnl
dnl  SKEW_COPY_ROW
dnl     |
dnl     |-- SKEW_COPY_ROW_HEAD
dnl     |      |
dnl     |      |-- SKEW_COPY_ROW_NAME
dnl     |      |
dnl     |      |-- SKEW_COPY_ROW_PARAMS
dnl     |
dnl     |-- SKEW_COPY_ROW_BODY
dnl
dnl
dnl  Usage: 
dnl    SKEW_COPY_ROW($1)
dnl    SKEW_COPY_ROW_HEAD($1)
dnl    SKEW_COPY_ROW_NAME($1)
dnl    SKEW_COPY_ROW_PARAMS($1)
dnl    SKEW_COPY_ROW_COMMENT($1)
dnl    SKEW_COPY_ROW_BODY($1)
dnl
dnl    $1 is the type of the matrix.
dnl
dnl
define(`SKEW_COPY_ROW_NAME', `$1skew_copy_row')dnl
dnl
dnl
define(`SKEW_COPY_ROW_PARAMS', 
  `enum blas_order_type order, enum blas_uplo_type uplo, 
   enum blas_side_type side, int n, $1_array a, int lda, 
   $1_array a_vec, int row')dnl
dnl
dnl
define(`SKEW_COPY_ROW_HEAD', 
  `void SKEW_COPY_ROW_NAME($1)(SKEW_COPY_ROW_PARAMS($1))')dnl
dnl
dnl
define(`SKEW_COPY_ROW_COMMENT', `
/*
 *  Copies the given row of skew matrix a into the supplied vector.
 */')dnl
dnl
dnl
define(`SKEW_COPY_ROW_BODY', `{
  int conj_flag;
  int i;
  int ai, incai1, incai2;
  int vi, incvi = 1;

  DECLARE(a_elem, $1_type)
  PTR_CAST(a, $1_type, `const')
  PTR_CAST(a_vec, $1_type)

  if ((side == blas_left_side  && uplo == blas_upper) ||
      (side == blas_right_side && uplo == blas_lower))
    conj_flag = 1;
  else
    conj_flag = 0;

  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai1 = 1;
    incai2 = lda;
    ai = lda * row;
  } else {
    incai1 = lda;
    incai2 = 1;
    ai = row;
  }

  INC_ADJUST(incai1, $1_type)
  INC_ADJUST(incai2, $1_type)
  INC_ADJUST(ai, $1_type)
  INC_ADJUST(incvi, $1_type)

  for (i = 0, vi = 0; i < row; i++, ai += incai1, vi += incvi) {
    GET_VECTOR_ELEMENT(a_elem, a_i, ai, $1_type)
    if (conj_flag == 1) {
      NEGATE(a_elem, $1_type)
    }
    SET_VECTOR_ELEMENT(a_vec_i, vi, a_elem, $1_type)
  }

  for (; i < n; i++, ai += incai2, vi += incvi) {
    GET_VECTOR_ELEMENT(a_elem, a_i, ai, $1_type)
    if (conj_flag == 0) {
      NEGATE(a_elem, $1_type)
    }
    SET_VECTOR_ELEMENT(a_vec_i, vi, a_elem, $1_type)
  }
}')dnl
dnl
dnl
define(`SKEW_COPY_ROW', 
  `SKEW_COPY_ROW_HEAD($1)
   SKEW_COPY_ROW_COMMENT($1)
   SKEW_COPY_ROW_BODY($1)')dnl
dnl
dnl
dnl
dnl  Usage: 
dnl    HE_COPY_ROW($1)
dnl    HE_COPY_ROW_HEAD($1)
dnl    HE_COPY_ROW_NAME($1)
dnl    HE_COPY_ROW_PARAMS($1)
dnl    HE_COPY_ROW_COMMENT($1)
dnl    HE_COPY_ROW_BODY($1)
dnl
dnl    $1 is the type of the matrix.
dnl
dnl
define(`HE_COPY_ROW_NAME', `$1he_copy_row')dnl
dnl
dnl
define(`HE_COPY_ROW_HEAD', 
  `void HE_COPY_ROW_NAME($1)(enum blas_order_type order, dnl
       enum blas_uplo_type uplo, enum blas_side_type side, dnl
       int n, $1_array a, int lda, $1_array a_vec, int row)')dnl
dnl
dnl
define(`HE_COPY_ROW_COMMENT', `
/* Copies the given row of matrix a into the supplied vector. */')dnl
dnl
dnl
define(`HE_COPY_ROW_BODY', `{
  int conj_flag;
  int i;
  int ai, incai1, incai2;
  int vi, incvi = 1;

  DECLARE(a_elem, $1_type)
  PTR_CAST(a, $1_type, `const')
  PTR_CAST(a_vec, $1_type)

  if ((side == blas_left_side  && uplo == blas_upper) ||
      (side == blas_right_side && uplo == blas_lower))
    conj_flag = 1;
  else
    conj_flag = 0;

  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai1 = 1;
    incai2 = lda;
    ai = lda * row;
  } else {
    incai1 = lda;
    incai2 = 1;
    ai = row;
  }

  INC_ADJUST(incai1, $1_type)
  INC_ADJUST(incai2, $1_type)
  INC_ADJUST(ai, $1_type)
  INC_ADJUST(incvi, $1_type)

  for (i = 0, vi = 0; i < row; i++, ai += incai1, vi += incvi) {
    GET_VECTOR_ELEMENT(a_elem, a_i, ai, $1_type)
    if (conj_flag == 1) {
        CONJ_AUX(a_elem, $1_type)
    }
    SET_VECTOR_ELEMENT(a_vec_i, vi, a_elem, $1_type)
  }
  for (; i < n; i++, ai += incai2, vi += incvi) {
    GET_VECTOR_ELEMENT(a_elem, a_i, ai, $1_type)
    if (conj_flag == 0) {
        CONJ_AUX(a_elem, $1_type)
    }
    SET_VECTOR_ELEMENT(a_vec_i, vi, a_elem, $1_type)
  }
}')dnl
dnl
dnl
define(`HE_COPY_ROW', 
  `HE_COPY_ROW_HEAD($1)
   HE_COPY_ROW_COMMENT($1)
   HE_COPY_ROW_BODY($1)')dnl
dnl
dnl
define(`HE_PRINT_MATRIX_NAME', `$1he_print_matrix')dnl
dnl
define(`HE_PRINT_MATRIX_HEAD', 
  `void HE_PRINT_MATRIX_NAME($1) ($1_array a, int n, int lda, 
     enum blas_order_type order, enum blas_uplo_type uplo)')dnl
dnl
define(`HE_PRINT_MATRIX_BODY', 
  `int ai, aij;
   int incai, incaij1, incaij2;
   int i, j;

   DECLARE(a_elem, $1_type)
   PTR_CAST(a, $1_type, `const')
 
   if ((order == blas_colmajor && uplo == blas_upper) ||
       (order == blas_rowmajor && uplo == blas_lower)) {
     incai = lda;
     incaij1 = 1;
     incaij2 = lda;
   } else {
     incai = 1;
     incaij1 = lda;
     incaij2 = 1;
   }

   INC_ADJUST(incai, $1_type)
   INC_ADJUST(incaij1, $1_type)
   INC_ADJUST(incaij2, $1_type)

   for (i = 0, ai = 0; i < n; i++, ai += incai) {

     for (j = 0, aij = ai; j < i; j++, aij += incaij1) {
       GET_VECTOR_ELEMENT(a_elem, a_i, aij, $1_type)
       if (uplo == blas_upper) {
           CONJ_AUX(a_elem, $1_type)
       }
       PRINT_NUMBER(a_elem, $1_type)
     }
     GET_VECTOR_ELEMENT(a_elem[0], a_i, ai, REAL_TYPE($1_type))
     a_elem[1] = 0.0;
     PRINT_NUMBER(a_elem, $1_type)
     j++;
     aij += incaij2;
     for (; j < n; j++, aij += incaij2) {
       GET_VECTOR_ELEMENT(a_elem, a_i, aij, $1_type)
       if (uplo == blas_lower) {
           CONJ_AUX(a_elem, $1_type)
       }
       PRINT_NUMBER(a_elem, $1_type)
     }
     printf("\n");
  }
')dnl
dnl
define(`HE_PRINT_MATRIX', 
  `HE_PRINT_MATRIX_HEAD($1) {
   HE_PRINT_MATRIX_BODY($1)
  }')dnl
dnl
dnl
define(`PROTOTYPES', `dnl
HE_COPY_ROW_HEAD(c);
HE_COPY_ROW_HEAD(z);

HE_PRINT_MATRIX_HEAD(c);
HE_PRINT_MATRIX_HEAD(z);

SKEW_COMMIT_ROW_HEAD(s);
SKEW_COMMIT_ROW_HEAD(d);

SKEW_COPY_ROW_HEAD(s);
SKEW_COPY_ROW_HEAD(d);
')dnl
dnl
dnl
define(`SOURCE', `dnl
#include <stdio.h>
#include "blas_extended.h"

HE_COPY_ROW(c)
HE_COPY_ROW(z)

HE_PRINT_MATRIX(c)
HE_PRINT_MATRIX(z)

SKEW_COMMIT_ROW(s)
SKEW_COMMIT_ROW(d)

SKEW_COPY_ROW(s)
SKEW_COPY_ROW(d)
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
