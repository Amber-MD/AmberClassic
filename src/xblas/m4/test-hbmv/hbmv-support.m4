dnl
dnl
include(cblas.m4)dnl
include(test-common.m4)dnl
dnl
dnl
dnl  HBMV_COMMIT_ROW
dnl     |
dnl     |-- HBMV_COMMIT_ROW_HEAD
dnl     |      |
dnl     |      |-- HBMV_COMMIT_ROW_NAME
dnl     |      |
dnl     |      |-- HBMV_COMMIT_ROW_PARAMS
dnl     |
dnl     |-- HBMV_COMMIT_ROW_BODY
dnl
dnl
dnl  Usage: 
dnl    HBMV_COMMIT_ROW($1)
dnl    HBMV_COMMIT_ROW_HEAD($1)
dnl    HBMV_COMMIT_ROW_NAME($1)
dnl    HBMV_COMMIT_ROW_PARAMS($1)
dnl    HBMV_COMMIT_ROW_COMMENT($1)
dnl    HBMV_COMMIT_ROW_BODY($1)
dnl
dnl    $1 is the type of the matrix.
dnl
define(`HBMV_COMMIT_ROW_NAME', `$1hbmv_commit_row')dnl
dnl
dnl
define(`HBMV_COMMIT_ROW_PARAMS', 
  `enum blas_order_type order, enum blas_uplo_type uplo, dnl
   int n, $1_array a, int k, int lda, dnl
   $1_array a_vec, int row')dnl
dnl
dnl
define(`HBMV_COMMIT_ROW_HEAD', 
  `void HBMV_COMMIT_ROW_NAME($1)(HBMV_COMMIT_ROW_PARAMS($1))')dnl
dnl
dnl
define(`HBMV_COMMIT_ROW_COMMENT', `
/*
 *  Copies the given vector into the given row of symmetric matrix a.
 */')dnl
dnl
dnl
define(`HBMV_COMMIT_ROW_BODY', `{

  int ai;
  int i;
  int incai1;
  int incai2;
  int vi;
  int incvi = 1;
  int loopmax;

  DECLARE(a_elem, $1_type)
  PTR_CAST(a, $1_type)
  PTR_CAST(a_vec, $1_type, `const')

  if (order == blas_colmajor) {
        if (uplo == blas_upper) {
                incai1 = 1;
                incai2 = lda - 1; 
        } else {
                incai1 = lda - 1;
                incai2 = 1;
        }
  } else {
        if (uplo == blas_upper) {
                incai1 = lda - 1;
                incai2 = 1; 
        } else {
                incai1 = 1;
                incai2 = lda - 1;
        }
  }

  ai = 0;
  if ( (uplo == blas_upper && order == blas_colmajor) ||
         (uplo == blas_lower && order == blas_rowmajor)) {
                /* starting place */
    ai = lda * row + ( (row<k) ? (k - row) : 0); 
  } else {
                /* starting place */
    ai = (row>k) ? ( k + lda*(row-k)) : row; 
  }

  INC_ADJUST(incai1, $1_type)
  INC_ADJUST(incai2, $1_type)
  INC_ADJUST(ai, $1_type)
  INC_ADJUST(incvi, $1_type)

  for (i = 0, vi = 0; i < row-k; i++, vi += incvi) {
        /* this is a wasteful loop but important */
  }

  for (; i < row; i++, ai += incai1, vi += incvi) {
    GET_VECTOR_ELEMENT(a_elem, a_vec_i, vi, $1_type)
    if (uplo == blas_upper) {
        CONJ_AUX(a_elem, $1_type);
    }
    SET_VECTOR_ELEMENT(a_i, ai, a_elem, $1_type)
  }

  loopmax = MIN(row+k+1, n);
  for (; i < loopmax; i++, ai += incai2, vi += incvi) {
    GET_VECTOR_ELEMENT(a_elem, a_vec_i, vi, $1_type)
    if (uplo == blas_lower) {
        CONJ_AUX(a_elem, $1_type);
    }
    SET_VECTOR_ELEMENT(a_i, ai, a_elem, $1_type)
  }
}')dnl
dnl
dnl
define(`HBMV_COMMIT_ROW', 
  `HBMV_COMMIT_ROW_HEAD($1)
   HBMV_COMMIT_ROW_COMMENT($1)
   HBMV_COMMIT_ROW_BODY($1)')dnl
dnl
dnl
dnl
dnl
dnl
dnl  HBMV_COPY_ROW
dnl     |
dnl     |-- HBMV_COPY_ROW_HEAD
dnl     |      |
dnl     |      |-- HBMV_COPY_ROW_NAME
dnl     |      |
dnl     |      |-- HBMV_COPY_ROW_PARAMS
dnl     |
dnl     |-- HBMV_COPY_ROW_BODY
dnl
dnl
dnl  Usage: 
dnl    HBMV_COPY_ROW($1)
dnl    HBMV_COPY_ROW_HEAD($1)
dnl    HBMV_COPY_ROW_NAME($1)
dnl    HBMV_COPY_ROW_PARAMS($1)
dnl    HBMV_COPY_ROW_COMMENT($1)
dnl    HBMV_COPY_ROW_BODY($1)
dnl
dnl    $1 is the type of the matrix.
dnl
dnl
define(`HBMV_COPY_ROW_NAME', `$1hbmv_copy_row')dnl
dnl
dnl
define(`HBMV_COPY_ROW_PARAMS', 
  `enum blas_order_type order, enum blas_uplo_type uplo, dnl
   int n, $1_array a, int k, int lda, dnl
   $1_array a_vec, int row')dnl
dnl
dnl
define(`HBMV_COPY_ROW_HEAD', 
  `void HBMV_COPY_ROW_NAME($1)(HBMV_COPY_ROW_PARAMS($1))')dnl
dnl
dnl
define(`HBMV_COPY_ROW_COMMENT', `
/*
 *  Copies the given row of matrix a into the supplied vector.
 */')dnl
dnl
dnl
define(`HBMV_COPY_ROW_BODY', `{

  int ai;
  int i;
  int incai1;
  int incai2;
  int vi;
  int incvi = 1;
  int loopmax;

  DECLARE(a_elem, $1_type)
  PTR_CAST(a, $1_type, `const')
  PTR_CAST(a_vec, $1_type)


  if (order == blas_colmajor) {
        if (uplo == blas_upper) {
                incai1 = 1;
                incai2 = lda - 1; 
        } else {
                incai1 = lda - 1;
                incai2 = 1;
        }
  } else {
        if (uplo == blas_upper) {
                incai1 = lda - 1;
                incai2 = 1; 
        } else {
                incai1 = 1;
                incai2 = lda - 1;
        }
  }

  ai = 0;
  if ( (uplo == blas_upper && order == blas_colmajor) ||
         (uplo == blas_lower && order == blas_rowmajor)) {
                /* starting place */
    ai = lda * row + ( (row<k) ? (k - row) : 0); 
  } else {
                /* starting place */
    ai = (row>k) ? ( k + lda*(row-k)) : row; 
  }

  INC_ADJUST(incai1, $1_type)
  INC_ADJUST(incai2, $1_type)
  INC_ADJUST(ai, $1_type)
  INC_ADJUST(incvi, $1_type)

  for (i = 0, vi = 0; i < row-k; i++, vi += incvi) {
        SET_ZERO_VECTOR_ELEMENT(a_vec_i, vi, $1_type)
  }

  for (; i < row; i++, ai += incai1, vi += incvi) {
    GET_VECTOR_ELEMENT(a_elem, a_i, ai, $1_type)
    if (uplo == blas_upper) {
        CONJ_AUX(a_elem, $1_type);
    }
    SET_VECTOR_ELEMENT(a_vec_i, vi, a_elem, $1_type)
  }

  loopmax = MIN(row+k+1, n);
  for (; i < loopmax; i++, ai += incai2, vi += incvi) {
    GET_VECTOR_ELEMENT(a_elem, a_i, ai, $1_type)
    if (i == row) {
      ZERO_IMAG_PART(a_elem, $1_type)
    }
    if (uplo == blas_lower) {
        CONJ_AUX(a_elem, $1_type);
    }
    SET_VECTOR_ELEMENT(a_vec_i, vi, a_elem, $1_type)
  }
  for (; i < n; i++, vi += incvi) {
        SET_ZERO_VECTOR_ELEMENT(a_vec_i, vi, $1_type)
  }
}')dnl
dnl
dnl
define(`HBMV_COPY_ROW', 
  `HBMV_COPY_ROW_HEAD($1)
   HBMV_COPY_ROW_COMMENT($1)
   HBMV_COPY_ROW_BODY($1)')dnl
dnl
dnl
define(`PRINT_HBMV_MATRIX_NAME', `$1print_hbmv_matrix')dnl
dnl
define(`PRINT_HBMV_MATRIX_HEAD', 
  `void PRINT_HBMV_MATRIX_NAME($1)($1_array a, int n, int k, dnl
     int lda, enum blas_order_type order, enum blas_uplo_type uplo)')
dnl
define(`PRINT_HBMV_MATRIX_BODY', 
  `
  int row;
  DECLARE_VECTOR(x, $1_type)
  MALLOC_VECTOR(x, $1_type, n)
  
  for (row = 0; row < n; row++ ) {
        $1hbmv_copy_row(order, uplo, n, a, k, lda, x, row);
        $1print_vector(x, n, 1, NULL);
  }
  printf("\n");
  FREE_VECTOR(x, $1_type)
')dnl
dnl
define(`PRINT_HBMV_MATRIX', 
  `PRINT_HBMV_MATRIX_HEAD($1) {
   PRINT_HBMV_MATRIX_BODY($1)
  }')dnl
dnl
dnl
dnl
dnl  SKEW_COMMIT_ROW_HBMV
dnl     |
dnl     |-- SKEW_COMMIT_ROW_HBMV_HEAD
dnl     |      |
dnl     |      |-- SKEW_COMMIT_ROW_HBMV_NAME
dnl     |      |
dnl     |      |-- SKEW_COMMIT_ROW_HBMV_PARAMS
dnl     |
dnl     |-- SKEW_COMMIT_ROW_HBMV_BODY
dnl
dnl
dnl  Usage: 
dnl    SKEW_COMMIT_ROW_HBMV($1)
dnl    SKEW_COMMIT_ROW_HBMV_HEAD($1)
dnl    SKEW_COMMIT_ROW_HBMV_NAME($1)
dnl    SKEW_COMMIT_ROW_HBMV_PARAMS($1)
dnl    SKEW_COMMIT_ROW_HBMV_COMMENT($1)
dnl    SKEW_COMMIT_ROW_HBMV_BODY($1)
dnl
dnl    $1 is the type of the matrix.
dnl
define(`SKEW_COMMIT_ROW_HBMV_NAME', `$1skew_commit_row_hbmv')dnl
dnl
dnl
define(`SKEW_COMMIT_ROW_HBMV_PARAMS', 
  `enum blas_order_type order, enum blas_uplo_type uplo, dnl
   int n, $1_array a, int k, int lda, dnl
   $1_array a_vec, int row')dnl
dnl
dnl
define(`SKEW_COMMIT_ROW_HBMV_HEAD', 
  `void SKEW_COMMIT_ROW_HBMV_NAME($1)(SKEW_COMMIT_ROW_HBMV_PARAMS($1))')dnl
dnl
dnl
define(`SKEW_COMMIT_ROW_HBMV_COMMENT', `
/*
 *  Copies the given vector into the given row of symmetric matrix a.
 */')dnl
dnl
dnl
define(`SKEW_COMMIT_ROW_HBMV_BODY', `{

  int ai;
  int i;
  int incai1;
  int incai2;
  int vi;
  int incvi = 1;
  int loopmax;

  DECLARE(a_elem, $1_type)
  PTR_CAST(a, $1_type)
  PTR_CAST(a_vec, $1_type, `const')

  if (order == blas_colmajor) {
        if (uplo == blas_upper) {
                incai1 = 1;
                incai2 = lda - 1; 
        } else {
                incai1 = lda - 1;
                incai2 = 1;
        }
  } else {
        if (uplo == blas_upper) {
                incai1 = lda - 1;
                incai2 = 1; 
        } else {
                incai1 = 1;
                incai2 = lda - 1;
        }
  }

  ai = 0;
  if ( (uplo == blas_upper && order == blas_colmajor) ||
         (uplo == blas_lower && order == blas_rowmajor)) {
                /* starting place */
    ai = lda * row + ( (row<k) ? (k - row) : 0); 
  } else {
                /* starting place */
    ai = (row>k) ? ( k + lda*(row-k)) : row; 
  }

  INC_ADJUST(incai1, $1_type)
  INC_ADJUST(incai2, $1_type)
  INC_ADJUST(ai, $1_type)
  INC_ADJUST(incvi, $1_type)

  for (i = 0, vi = 0; i < row-k; i++, vi += incvi) {
        /* this is a wasteful loop but important */
  }

  for (; i < row; i++, ai += incai1, vi += incvi) {
    GET_VECTOR_ELEMENT(a_elem, a_vec_i, vi, $1_type)
    if (uplo == blas_upper) {
        NEGATE(a_elem, $1_type);
    }
    SET_VECTOR_ELEMENT(a_i, ai, a_elem, $1_type)
  }

  loopmax = MIN(row+k+1, n);
  for (; i < loopmax; i++, ai += incai2, vi += incvi) {
    GET_VECTOR_ELEMENT(a_elem, a_vec_i, vi, $1_type)
    if (uplo == blas_lower) {
        NEGATE(a_elem, $1_type);
    }
    SET_VECTOR_ELEMENT(a_i, ai, a_elem, $1_type)
  }
}')dnl
dnl
dnl
define(`SKEW_COMMIT_ROW_HBMV', 
  `SKEW_COMMIT_ROW_HBMV_HEAD($1)
   SKEW_COMMIT_ROW_HBMV_COMMENT($1)
   SKEW_COMMIT_ROW_HBMV_BODY($1)')dnl
dnl
dnl
dnl
dnl
dnl
dnl  SKEW_COPY_ROW_HBMV
dnl     |
dnl     |-- SKEW_COPY_ROW_HBMV_HEAD
dnl     |      |
dnl     |      |-- SKEW_COPY_ROW_HBMV_NAME
dnl     |      |
dnl     |      |-- SKEW_COPY_ROW_HBMV_PARAMS
dnl     |
dnl     |-- SKEW_COPY_ROW_HBMV_BODY
dnl
dnl
dnl  Usage: 
dnl    SKEW_COPY_ROW_HBMV($1)
dnl    SKEW_COPY_ROW_HBMV_HEAD($1)
dnl    SKEW_COPY_ROW_HBMV_NAME($1)
dnl    SKEW_COPY_ROW_HBMV_PARAMS($1)
dnl    SKEW_COPY_ROW_HBMV_COMMENT($1)
dnl    SKEW_COPY_ROW_HBMV_BODY($1)
dnl
dnl    $1 is the type of the matrix.
dnl
dnl
define(`SKEW_COPY_ROW_HBMV_NAME', `$1skew_copy_row_hbmv')dnl
dnl
dnl
define(`SKEW_COPY_ROW_HBMV_PARAMS', 
  `enum blas_order_type order, enum blas_uplo_type uplo, 
   int n, $1_array a, int k, int lda, 
   $1_array a_vec, int row')dnl
dnl
dnl
define(`SKEW_COPY_ROW_HBMV_HEAD', 
  `void SKEW_COPY_ROW_HBMV_NAME($1)(SKEW_COPY_ROW_HBMV_PARAMS($1))')dnl
dnl
dnl
define(`SKEW_COPY_ROW_HBMV_COMMENT', `
/*
 *  Copies the given row of matrix a into the supplied vector.
 */')dnl
dnl
dnl
define(`SKEW_COPY_ROW_HBMV_BODY', `{

  int ai;
  int i;
  int incai1;
  int incai2;
  int vi;
  int incvi = 1;
  int loopmax;

  DECLARE(a_elem, $1_type)
  PTR_CAST(a, $1_type, `const')
  PTR_CAST(a_vec, $1_type)


  if (order == blas_colmajor) {
        if (uplo == blas_upper) {
                incai1 = 1;
                incai2 = lda - 1; 
        } else {
                incai1 = lda - 1;
                incai2 = 1;
        }
  } else {
        if (uplo == blas_upper) {
                incai1 = lda - 1;
                incai2 = 1; 
        } else {
                incai1 = 1;
                incai2 = lda - 1;
        }
  }

  ai = 0;
  if ( (uplo == blas_upper && order == blas_colmajor) ||
         (uplo == blas_lower && order == blas_rowmajor)) {
                /* starting place */
    ai = lda * row + ( (row<k) ? (k - row) : 0); 
  } else {
                /* starting place */
    ai = (row>k) ? ( k + lda*(row-k)) : row; 
  }

  INC_ADJUST(incai1, $1_type)
  INC_ADJUST(incai2, $1_type)
  INC_ADJUST(ai, $1_type)
  INC_ADJUST(incvi, $1_type)

  for (i = 0, vi = 0; i < row-k; i++, vi += incvi) {
        SET_ZERO_VECTOR_ELEMENT(a_vec_i, vi, $1_type)
  }

  for (; i < row; i++, ai += incai1, vi += incvi) {
    GET_VECTOR_ELEMENT(a_elem, a_i, ai, $1_type)
    if (uplo == blas_upper) {
        NEGATE(a_elem, $1_type);
    }
    SET_VECTOR_ELEMENT(a_vec_i, vi, a_elem, $1_type)
  }

  loopmax = MIN(row+k+1, n);
  for (; i < loopmax; i++, ai += incai2, vi += incvi) {
    GET_VECTOR_ELEMENT(a_elem, a_i, ai, $1_type)
    if (i == row) {
      ZERO_IMAG_PART(a_elem, $1_type)
    }
    if (uplo == blas_lower) {
        NEGATE(a_elem, $1_type);
    }
    SET_VECTOR_ELEMENT(a_vec_i, vi, a_elem, $1_type)
  }
  for (; i < n; i++, vi += incvi) {
        SET_ZERO_VECTOR_ELEMENT(a_vec_i, vi, $1_type)
  }


}')dnl
dnl
dnl
define(`SKEW_COPY_ROW_HBMV', 
  `SKEW_COPY_ROW_HBMV_HEAD($1)
   SKEW_COPY_ROW_HBMV_COMMENT($1)
    SKEW_COPY_ROW_HBMV_BODY($1)')dnl
dnl
dnl
define(`PROTOTYPES', `dnl
SKEW_COMMIT_ROW_HBMV_HEAD(s);
SKEW_COMMIT_ROW_HBMV_HEAD(d);

HBMV_COMMIT_ROW_HEAD(c);
HBMV_COMMIT_ROW_HEAD(z);

SKEW_COPY_ROW_HBMV_HEAD(s);
SKEW_COPY_ROW_HBMV_HEAD(d);

HBMV_COPY_ROW_HEAD(c);
HBMV_COPY_ROW_HEAD(z);

PRINT_HBMV_MATRIX_HEAD(c);
PRINT_HBMV_MATRIX_HEAD(z);
')dnl
dnl
dnl
define(`SOURCE', `dnl
#include <stdio.h>
#include <stdlib.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"

SKEW_COMMIT_ROW_HBMV(s)
SKEW_COMMIT_ROW_HBMV(d)

HBMV_COMMIT_ROW(c)
HBMV_COMMIT_ROW(z)

SKEW_COPY_ROW_HBMV(s)
SKEW_COPY_ROW_HBMV(d)

HBMV_COPY_ROW(c)
HBMV_COPY_ROW(z)

PRINT_HBMV_MATRIX(c)
PRINT_HBMV_MATRIX(z)
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
