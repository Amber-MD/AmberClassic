include(cblas.m4)dnl
include(test-common.m4)dnl
dnl
dnl
dnl
dnl Usage:  GE_COMMIT_ROW_NAME($1)
dnl         GE_COMMIT_ROW_BODY($1)
dnl         GE_COMMIT_ROW_HEADER($1)
dnl         GE_COMMIT_ROW($1)
dnl     $1 -- type of vector/array
dnl
define(`GE_COMMIT_ROW_NAME', 
  `$1ge_commit_row')dnl
dnl
dnl
define(`GE_COMMIT_ROW_HEAD',
  `void GE_COMMIT_ROW_NAME($1)(enum blas_order_type order, dnl
       enum blas_trans_type trans, int m, int n, $1_array a, dnl
       int lda, $1_array a_vec, int row)')dnl
dnl
dnl
define(`GE_COMMIT_ROW_BODY', `{
  int ai;
  int i;
  int incai;
  int incvi = 1;
  int vi;

  DECLARE(a_elem, $1_type)
  PTR_CAST(a, $1_type)
  PTR_CAST(a_vec, $1_type, `const')

  if (order == blas_colmajor) {
    if (trans == blas_no_trans) {
      incai = lda;
      ai = row;
    } else {
      incai = 1;
      ai = row * lda;
    }
  } else {
    if (trans == blas_no_trans) {
      incai = 1;
      ai = row * lda;
    } else {
      incai = lda;
      ai = row;
    }
  }

  INC_ADJUST(ai, $1_type)
  INC_ADJUST(incai, $1_type)
  INC_ADJUST(incvi, $1_type)

  for (i = 0, vi = 0; i < n; i++, vi += incvi, ai += incai) {
    GET_VECTOR_ELEMENT(a_elem, a_vec_i, vi, $1_type)
    IF_COMPLEX($1_type, `
      if (trans == blas_conj_trans) {
        CONJ_AUX(a_elem, $1_type)
      }
    ')dnl
    SET_VECTOR_ELEMENT(a_i, ai, a_elem, $1_type)
  }
}')dnl
dnl
dnl
define(`GE_COMMIT_ROW',
  `GE_COMMIT_ROW_HEAD($1) GE_COMMIT_ROW_BODY($1)')dnl
dnl
dnl
dnl
dnl
dnl Usage:  GE_COMMIT_COL_NAME($1)
dnl         GE_COMMIT_COL_BODY($1)
dnl         GE_COMMIT_COL_HEAD($1)
dnl         GE_COMMIT_COL($1)
dnl     $1 -- type of vector/array
dnl
define(`GE_COMMIT_COL_NAME', `$1ge_commit_col')dnl
dnl
dnl
define(`GE_COMMIT_COL_HEAD',
  `void GE_COMMIT_COL_NAME($1)(enum blas_order_type order, dnl
       enum blas_trans_type trans, int m, int n, $1_array a, dnl
       int lda, $1_array a_vec, int col)')dnl
dnl
dnl
define(`GE_COMMIT_COL_BODY', `{
  int ai;
  int i;
  int incai;
  int incvi = 1;
  int vi;

  DECLARE(a_elem, $1_type)
  PTR_CAST(a, $1_type)
  PTR_CAST(a_vec, $1_type, `const')

  if (order == blas_colmajor) {
    if (trans == blas_no_trans) {
      incai = 1;
      ai = col * lda;
    } else {
      incai = lda;
      ai = col;
    }
  } else {
    if (trans == blas_no_trans) {
      incai = lda;
      ai = col;
    } else {
      incai = 1;
      ai = col * lda;
    }
  }

  INC_ADJUST(ai, $1_type)
  INC_ADJUST(incai, $1_type)
  INC_ADJUST(incvi, $1_type)

  for (i = 0, vi = 0; i < m; i++, vi += incvi, ai += incai) {
    GET_VECTOR_ELEMENT(a_elem, a_vec_i, vi, $1_type)
    IF_COMPLEX($1_type, `
      if (trans == blas_conj_trans) {
        CONJ_AUX(a_elem, $1_type)
      }
    ')dnl
    SET_VECTOR_ELEMENT(a_i, ai, a_elem, $1_type)
  }
}')dnl
dnl
dnl
define(`GE_COMMIT_COL',
  `GE_COMMIT_COL_HEAD($1) GE_COMMIT_COL_BODY($1)')dnl
dnl
dnl
dnl
dnl
dnl Usage: GE_COPY_ROW($1)
dnl        GE_COPY_ROW_NAME($1)
dnl        GE_COPY_ROW_HEAD($1)
dnl        GE_COPY_ROW_BODY($1)
dnl
define(`GE_COPY_ROW_NAME', $1ge_copy_row)dnl
dnl
dnl
define(`GE_COPY_ROW_HEAD', 
  `void GE_COPY_ROW_NAME($1)(enum blas_order_type order, dnl
       enum blas_trans_type trans, int m, int n, $1_array a, dnl
       int lda, $1_array a_vec, int row)')dnl
dnl
dnl
define(`GE_COPY_ROW', 
  `GE_COPY_ROW_HEAD($1) { 
    GE_COPY_ROW_BODY($1)
  }')dnl
dnl
dnl
define(`GE_COPY_ROW_BODY', `
  int ai;
  int incai;
  int i;
  int incvi = 1;        
  int vi;

  DECLARE(a_elem, $1_type)
  PTR_CAST(a, $1_type, `const')
  PTR_CAST(a_vec, $1_type)

  if (order == blas_colmajor) {
    if (trans == blas_no_trans) {
      ai = row;
      incai = lda;
    } else {
      ai = row * lda;
      incai = 1;
    }
  } else {
    if (trans == blas_no_trans) {
      ai = row * lda;
      incai = 1;
    } else {
      ai = row;
      incai = lda;
    }

  }

  INC_ADJUST(ai, $1_type)
  INC_ADJUST(incai, $1_type)
  INC_ADJUST(incvi, $1_type)

  for (i = 0, vi = 0; i < n; i++, vi += incvi, ai += incai) {
    GET_VECTOR_ELEMENT(a_elem, a_i, ai, $1_type)
    IF_COMPLEX($1_type, `
      if (trans == blas_conj_trans) {
        CONJ_AUX(a_elem, $1_type)
      }
    ')dnl
    SET_VECTOR_ELEMENT(a_vec_i, vi, a_elem, $1_type)
  }
')dnl
dnl
dnl
dnl Usage: GE_COPY_COL($1)
dnl        GE_COPY_COL_NAME($1)
dnl        GE_COPY_COL_HEAD($1)
dnl        GE_COPY_COL_BODY($1)
dnl
define(`GE_COPY_COL_NAME', `$1ge_copy_col')dnl
dnl
dnl
define(`GE_COPY_COL_HEAD', 
  `void GE_COPY_COL_NAME($1)(enum blas_order_type order, dnl
       enum blas_trans_type trans, int m, int n, $1_array a, dnl
       int lda, $1_array a_vec, int col)')dnl
dnl
define(`GE_COPY_COL', 
  `GE_COPY_COL_HEAD($1) { 
    GE_COPY_COL_BODY($1)
  }')dnl
dnl
dnl
define(`GE_COPY_COL_BODY', `
  int ai;
  int incai;
  int i;
  int incvi = 1;
  int vi;

  DECLARE(a_elem, $1_type)
  PTR_CAST(a, $1_type, `const')
  PTR_CAST(a_vec, $1_type)

  if (order == blas_rowmajor) {
    if (trans == blas_no_trans) {
      ai = col;
      incai = lda;
    } else {
      ai = col * lda;
      incai = 1;
    }
  } else {
    if (trans == blas_no_trans) {
      ai = col * lda;
      incai = 1;
    } else {
      ai = col;
      incai = lda;
    }
  }

  INC_ADJUST(ai, $1_type)
  INC_ADJUST(incai, $1_type)
  INC_ADJUST(incvi, $1_type)

  for (i = 0, vi = 0; i < m; i++, vi += incvi, ai += incai) {
    GET_VECTOR_ELEMENT(a_elem, a_i, ai, $1_type)
    IF_COMPLEX($1_type, `
      if (trans == blas_conj_trans) {
        CONJ_AUX(a_elem, $1_type)
      }
    ')dnl
    SET_VECTOR_ELEMENT(a_vec_i, vi, a_elem, $1_type)
  }
')dnl
dnl
dnl
dnl Performs a <-- b
dnl Usage: GE_COPY_MATRIX($1)
dnl        GE_COPY_MATRIX_NAME($1)
dnl        GE_COPY_MATRIX_HEAD($1)
dnl        GE_COPY_MATRIX_BODY($1)
dnl
define(`GE_COPY_MATRIX_NAME', `$1ge_copy_matrix')
dnl
dnl
define(`GE_COPY_MATRIX_HEAD', 
  `void GE_COPY_MATRIX_NAME($1)(enum blas_order_type order, int m, dnl
       int n, $1_array a, int lda, $1_array b, int ldb)')
dnl
dnl
define(`GE_COPY_MATRIX_BODY', `
  int ai, aij;
  int bi, bij;
  int incai, incaij;
  int incbi, incbij;

  int i, j;

  DECLARE(elem, $1_type)
  PTR_CAST(a, $1_type)
  PTR_CAST(b, $1_type, `const')

  if (order == blas_colmajor) {
    incai = 1;
    incbi = 1;
    incaij = lda;
    incbij = ldb;
  } else {
    incai = lda;
    incbi = ldb;
    incaij = 1;
    incbij = 1;
  }

  INC_ADJUST(incai, $1_type)
  INC_ADJUST(incaij, $1_type)
  INC_ADJUST(incbi, $1_type)
  INC_ADJUST(incbij, $1_type)

  for (i = 0, ai = 0, bi = 0; i < m; i++, ai += incai, bi += incbi) {
    for (j = 0, aij = ai, bij = bi; j < n; j++, aij += incaij, bij += incbij) {
      GET_VECTOR_ELEMENT(elem, b_i, bij, $1_type)
      SET_VECTOR_ELEMENT(a_i, aij, elem, $1_type)
    }
  }
')
dnl
dnl
define(`GE_COPY_MATRIX', 
  `GE_COPY_MATRIX_HEAD($1) {
    GE_COPY_MATRIX_BODY($1)
  }')dnl
dnl
dnl
dnl Usage: GE_PRINT_MATRIX($1)
dnl        GE_PRINT_MATRIX_NAME($1)
dnl        GE_PRINT_MATRIX_HEAD($1)
dnl        GE_PRINT_MATRIX_BODY($1)
dnl
define(`GE_PRINT_MATRIX_NAME', `$1ge_print_matrix')dnl
dnl
define(`GE_PRINT_MATRIX_HEAD', 
  `void GE_PRINT_MATRIX_NAME($1) ($1_array a, int m, int n, dnl
       int lda, enum blas_order_type order, const char *name)')dnl
dnl
define(`GE_PRINT_MATRIX_BODY', 
  `int ai, aij;
   int incai, incaij;
   int i, j;

   PTR_CAST(a, $1_type, `const')
 
   if (order == blas_rowmajor) {
     incai = lda;
     incaij = 1;
   } else {
     incai = 1;
     incaij = lda;
   }

   INC_ADJUST(incai, $1_type)
   INC_ADJUST(incaij, $1_type)

   if (name) { printf("%s = ", name); }
   printf("[\n");
   for (i = 0, ai = 0; i < m; i++, ai += incai) {
     printf("  ");
     for (j = 0, aij = ai; j < n; j++, aij += incaij) {
       PRINT_ARRAY_ELEM(a_i, aij, $1_type)
     }
     printf("\n");
   }
   printf("];\n");
')dnl
dnl
dnl
define(`GE_PRINT_MATRIX', 
  `GE_PRINT_MATRIX_HEAD($1) {
   GE_PRINT_MATRIX_BODY($1)
  }')dnl
dnl
dnl
define(`PROTOTYPES', `
GE_COMMIT_ROW_HEAD(s);
GE_COMMIT_ROW_HEAD(d);
GE_COMMIT_ROW_HEAD(c);
GE_COMMIT_ROW_HEAD(z);

GE_COMMIT_COL_HEAD(s);
GE_COMMIT_COL_HEAD(d);
GE_COMMIT_COL_HEAD(c);
GE_COMMIT_COL_HEAD(z);

GE_COPY_ROW_HEAD(s);
GE_COPY_ROW_HEAD(d);
GE_COPY_ROW_HEAD(c);
GE_COPY_ROW_HEAD(z);

GE_COPY_COL_HEAD(s);
GE_COPY_COL_HEAD(d);
GE_COPY_COL_HEAD(c);
GE_COPY_COL_HEAD(z);

GE_COPY_MATRIX_HEAD(s);
GE_COPY_MATRIX_HEAD(d);
GE_COPY_MATRIX_HEAD(c);
GE_COPY_MATRIX_HEAD(z);

GE_PRINT_MATRIX_HEAD(s);
GE_PRINT_MATRIX_HEAD(d);
GE_PRINT_MATRIX_HEAD(c);
GE_PRINT_MATRIX_HEAD(z);
')dnl
dnl
dnl
define(`SOURCE', `dnl
#include <stdio.h>
#include "blas_extended.h"
#include "blas_extended_test.h"

GE_COMMIT_ROW(s)
GE_COMMIT_ROW(d)
GE_COMMIT_ROW(c)
GE_COMMIT_ROW(z)

GE_COMMIT_COL(s)
GE_COMMIT_COL(d)
GE_COMMIT_COL(c)
GE_COMMIT_COL(z)

GE_COPY_ROW(s)
GE_COPY_ROW(d)
GE_COPY_ROW(c)
GE_COPY_ROW(z)

GE_COPY_COL(s)
GE_COPY_COL(d)
GE_COPY_COL(c)
GE_COPY_COL(z)

GE_COPY_MATRIX(s)
GE_COPY_MATRIX(d)
GE_COPY_MATRIX(c)
GE_COPY_MATRIX(z)

GE_PRINT_MATRIX(s)
GE_PRINT_MATRIX(d)
GE_PRINT_MATRIX(c)
GE_PRINT_MATRIX(z)
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
