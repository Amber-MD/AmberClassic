dnl
dnl
include(cblas.m4)dnl
include(test-common.m4)dnl
dnl
dnl
dnl  Usage: 
dnl    SY_COMMIT_ROW($1)
dnl    SY_COMMIT_ROW_HEAD($1)
dnl    SY_COMMIT_ROW_NAME($1)
dnl    SY_COMMIT_ROW_COMMENT($1)
dnl    SY_COMMIT_ROW_BODY($1)
dnl
dnl    $1 is the type of the matrix.
dnl
define(`SY_COMMIT_ROW_NAME', `$1sy_commit_row')dnl
dnl
dnl
define(`SY_COMMIT_ROW_HEAD', 
  `void SY_COMMIT_ROW_NAME($1)(enum blas_order_type order, dnl
       enum blas_uplo_type uplo, int n, $1_array a, int lda, dnl
       $1_array a_vec, int row)')dnl
dnl
dnl
define(`SY_COMMIT_ROW_COMMENT', `
/*
 *  Copies the given vector into the given row of symmetric matrix a.
 */')dnl
dnl
dnl
define(`SY_COMMIT_ROW_BODY', `{
  int ai;
  int i;
  int incai1;
  int incai2;
  int vi;
  int incvi = 1;

  DECLARE(a_elem, $1_type)
  PTR_CAST(a, $1_type)
  PTR_CAST(a_vec, $1_type, `const')

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
    SET_VECTOR_ELEMENT(a_i, ai, a_elem, $1_type)
  }

  for (; i < n; i++, ai += incai2, vi += incvi) {
    GET_VECTOR_ELEMENT(a_elem, a_vec_i, vi, $1_type)
    SET_VECTOR_ELEMENT(a_i, ai, a_elem, $1_type)
  }
}')dnl
dnl
dnl
define(`SY_COMMIT_ROW', 
  `SY_COMMIT_ROW_HEAD($1)
   SY_COMMIT_ROW_COMMENT($1)
   SY_COMMIT_ROW_BODY($1)')dnl
dnl
dnl
dnl  Usage: 
dnl    SY_COPY_ROW($1)
dnl    SY_COPY_ROW_HEAD($1)
dnl    SY_COPY_ROW_NAME($1)
dnl    SY_COPY_ROW_COMMENT($1)
dnl    SY_COPY_ROW_BODY($1)
dnl
dnl    $1 is the type of the matrix.
dnl
define(`SY_COPY_ROW_NAME', `$1sy_copy_row')dnl
dnl
dnl
define(`SY_COPY_ROW_HEAD', 
  `void SY_COPY_ROW_NAME($1)(enum blas_order_type order, dnl
       enum blas_uplo_type uplo, int n, $1_array a, int lda, dnl
       $1_array a_vec, int row)')dnl
dnl
dnl
define(`SY_COPY_ROW_COMMENT', `
/*
 *  Copies the given row of matrix a into the supplied vector.
 */')dnl
dnl
dnl
define(`SY_COPY_ROW_BODY', `{
  int ai;
  int i;
  int incai1;
  int incai2;
  int vi;
  int incvi = 1;

  DECLARE(a_elem, $1_type)
  PTR_CAST(a, $1_type, `const')
  PTR_CAST(a_vec, $1_type)

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
    SET_VECTOR_ELEMENT(a_vec_i, vi, a_elem, $1_type)
  }

  for (; i < n; i++, ai += incai2, vi += incvi) {
    GET_VECTOR_ELEMENT(a_elem, a_i, ai, $1_type)
    SET_VECTOR_ELEMENT(a_vec_i, vi, a_elem, $1_type)
  }
}')dnl
dnl
dnl
define(`SY_COPY_ROW', 
  `SY_COPY_ROW_HEAD($1)
   SY_COPY_ROW_COMMENT($1)
   SY_COPY_ROW_BODY($1)')dnl
dnl
dnl
dnl
dnl
define(`SY_PRINT_MATRIX_NAME', `$1sy_print_matrix')dnl
dnl
dnl
define(`SY_PRINT_MATRIX_HEAD', 
  `void SY_PRINT_MATRIX_NAME($1) ($1_array a, int n, int lda, 
     enum blas_order_type order, enum blas_uplo_type uplo)')dnl
dnl
dnl
define(`SY_PRINT_MATRIX_BODY', 
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
       PRINT_NUMBER(a_elem, $1_type)
     }
     for (; j < n; j++, aij += incaij2) {
       GET_VECTOR_ELEMENT(a_elem, a_i, aij, $1_type)
       PRINT_NUMBER(a_elem, $1_type)
     }
     printf("\n");
  }
')dnl
dnl
dnl
define(`SY_PRINT_MATRIX', 
  `SY_PRINT_MATRIX_HEAD($1) {
   SY_PRINT_MATRIX_BODY($1)
  }')dnl
dnl
dnl
define(`PROTOTYPES', `
SY_COMMIT_ROW_HEAD(s);
SY_COMMIT_ROW_HEAD(d);
SY_COMMIT_ROW_HEAD(c);
SY_COMMIT_ROW_HEAD(z);

SY_COPY_ROW_HEAD(s);
SY_COPY_ROW_HEAD(d);
SY_COPY_ROW_HEAD(c);
SY_COPY_ROW_HEAD(z);

SY_PRINT_MATRIX_HEAD(s);
SY_PRINT_MATRIX_HEAD(d);
SY_PRINT_MATRIX_HEAD(c);
SY_PRINT_MATRIX_HEAD(z);
')dnl
dnl
dnl
define(`SOURCE', `dnl
#include <stdio.h>
#include "blas_extended.h"
SY_COMMIT_ROW(s)
SY_COMMIT_ROW(d)
SY_COMMIT_ROW(c)
SY_COMMIT_ROW(z)

SY_COPY_ROW(s)
SY_COPY_ROW(d)
SY_COPY_ROW(c)
SY_COPY_ROW(z)

SY_PRINT_MATRIX(s)
SY_PRINT_MATRIX(d)
SY_PRINT_MATRIX(c)
SY_PRINT_MATRIX(z)
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
