include(cblas.m4)dnl
include(packed.m4)dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: SPMV_COPY_ROW(ap_typeltr)
dnl        copy a row to the vector a_vec
dnl ----------------------------------------------------------------------
define(`SPMV_COPY_ROW_HEAD', `
void $1spmv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo, dnl
    int n, const $1_array a, $1_array a_vec, int row)')dnl
dnl
dnl
define(`SPMV_COPY_ROW',
`SPMV_COPY_ROW_HEAD($1)
/*
 * Purpose
 * =======
 *
 * Copy a row from a to a_vec
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of a; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether a is upper or lower
 *
 * n            (input) int
 *              Dimension of  and the length of vector x
 *
 * a           (input) $1_array
 *
 * a_vec            (input) $1_array
 *
 * row          (input) int
 *
 */
{
        int i, ind, inc = 1;
        PTR_CAST(a, $1_type, `const')
        PTR_CAST(a_vec, $1_type)
        DECLARE(tmp, $1_type)

        INC_ADJUST(inc, $1_type)

        if(((order == blas_rowmajor)&&(uplo == blas_upper))||
           ((order == blas_colmajor)&&(uplo == blas_lower)))
        {
                /* Pretend it is colmajor/lower.  We can do this in the
                   symmetric case. */
                ind = row * inc;
                for(i = 0; i < row; i++)
                {
                        GET_VECTOR_ELEMENT(tmp, a_i, ind, $1_type)
                        SET_VECTOR_ELEMENT(a_vec_i, i * inc, tmp, $1_type)
                        COLUMN_STRIDE(ind, n, i, blas_colmajor, blas_lower, inc)
                }

                ind = (row + ((2 * n - row - 1) * row) / 2) * inc;
                for(i = row; i < n; i++)
                {
                        GET_VECTOR_ELEMENT(tmp, a_i, ind, $1_type)
                        SET_VECTOR_ELEMENT(a_vec_i, i * inc, tmp, $1_type)
                        COLUMN_STRIDE(ind, n, i, blas_rowmajor, blas_upper, inc)
                }
        }
        else
        {
                /* Pretend it is rowmajor/lower.  We can do this in the
                   symmetric case. */
                ind = row * (row + 1) * inc / 2;
                for(i = 0; i < row; i++)
                {
                        GET_VECTOR_ELEMENT(tmp, a_i, ind, $1_type)
                        SET_VECTOR_ELEMENT(a_vec_i, i * inc, tmp, $1_type)
                        COLUMN_STRIDE(ind, n, i, blas_rowmajor, blas_lower, inc)
                }

                ind = (row + (row * (row + 1)) / 2) * inc;
                for(i = row; i < n; i++)
                {
                        GET_VECTOR_ELEMENT(tmp, a_i, ind, $1_type)
                        SET_VECTOR_ELEMENT(a_vec_i, i * inc, tmp, $1_type)
                        COLUMN_STRIDE(ind, n, i, blas_colmajor, blas_upper, inc)
                }
        }
}
')dnl
dnl
dnl
dnl
dnl -------------------------------------------------------------------------
dnl Usage: SPMV_COMMIT_ROW(ap_typeltr)
dnl        copy the vector a_vec to a
dnl -------------------------------------------------------------------------
define(`SPMV_COMMIT_ROW_HEAD', `
void $1spmv_commit_row(enum blas_order_type order, dnl
    enum blas_uplo_type uplo, int n, $1_array a, dnl
    const $1_array a_vec, int row)')dnl
dnl
dnl
define(`SPMV_COMMIT_ROW',
`SPMV_COMMIT_ROW_HEAD($1)
/*
 * Purpose
 * =======
 *
 * Copy a_vec to a
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of a; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether a is upper or lower
 *
 * n            (input) int
 *              Dimension of a and the length of vector a_vec
 *
 * a            (output) $1_array
 *
 * a_vec            (input) $1_array
 *
 * row          (input) int
 *
 */
{
        int i, ind, inc = 1;
        PTR_CAST(a, $1_type)
        PTR_CAST(a_vec, $1_type, `const')
        DECLARE(tmp, $1_type)

        INC_ADJUST(inc, $1_type)

        if(((order == blas_rowmajor)&&(uplo == blas_upper))||
           ((order == blas_colmajor)&&(uplo == blas_lower)))
        {
                /* Pretend it is rowmajor/upper.  We can do this in the
                   symmetric case. */
                ind = (row + ((2 * n - row - 1) * row) / 2) * inc;
                for(i = row; i < n; i++)
                {
                        GET_VECTOR_ELEMENT(tmp, a_vec_i, i * inc, $1_type)
                        SET_VECTOR_ELEMENT(a_i, ind, tmp, $1_type)
                        COLUMN_STRIDE(ind, n, i, blas_rowmajor, blas_upper, inc)
                }
        }
        else
        {
                /* Pretend it is colmajor/upper.  We can do this in the
                   symmetric case. */
                ind = (row + (row * (row + 1)) / 2) * inc;
                for(i = row; i < n; i++)
                {
                        GET_VECTOR_ELEMENT(tmp, a_vec_i, i * inc, $1_type)
                        SET_VECTOR_ELEMENT(a_i, ind, tmp, $1_type)
                        COLUMN_STRIDE(ind, n, i, blas_colmajor, blas_upper, inc)
                }
        }
}
')dnl
dnl  SPMV_PACK_MATRIX
dnl     |
dnl     |-- SPMV_PACK_MATRIX_HEAD
dnl     |      |
dnl     |      |-- SPMV_PACK_MATRIX_NAME
dnl     |      |
dnl     |      |-- SPMV_PACK_MATRIX_PARAMS
dnl     |
dnl     |-- SPMV_PACK_MATRIX_BODY
dnl
dnl
dnl  Usage: 
dnl    SPMV_PACK_MATRIX($1)
dnl    SPMV_PACK_MATRIX_HEAD($1)
dnl    SPMV_PACK_MATRIX_NAME($1)
dnl    SPMV_PACK_MATRIX_PARAMS($1)
dnl    SPMV_PACK_MATRIX_COMMENT($1)
dnl    SPMV_PACK_MATRIX_BODY($1)
dnl
dnl    $1 is the type of the matrix.
dnl
dnl
define(`SPMV_PACK_MATRIX_NAME', `$1spmv_pack_matrix')dnl
dnl
dnl
define(`SPMV_PACK_MATRIX_PARAMS', 
  `enum blas_order_type order, enum blas_uplo_type uplo, dnl 
   int n, $1_array a_packed, $1_array a_full, int lda')dnl
dnl
dnl
define(`SPMV_PACK_MATRIX_HEAD', 
  `void SPMV_PACK_MATRIX_NAME($1)(SPMV_PACK_MATRIX_PARAMS($1))')dnl
dnl
dnl
define(`SPMV_PACK_MATRIX_COMMENT', `
/*
 *  Packs the he matrix a_full into packed form a.
 */')dnl
dnl
dnl
define(`SPMV_PACK_MATRIX_BODY',
`{
        int row;
        DECLARE_VECTOR(a_row, $1_type);
        
        MALLOC_VECTOR(a_row, $1_type, n);       
        for (row = 0; row < n; row ++) {
                $1sy_copy_row(order, uplo, n, a_full, lda,
                        a_row, row);
                $1spmv_commit_row(order, uplo, n, a_packed,
                        a_row, row); 
        }

        FREE_VECTOR(a_row, $1_type)
}')dnl
dnl
dnl
define(`SPMV_PACK_MATRIX', 
  `SPMV_PACK_MATRIX_HEAD($1)
   SPMV_PACK_MATRIX_COMMENT($1)
   SPMV_PACK_MATRIX_BODY($1)')dnl
dnl
dnl
dnl
define(`PRINT_SPMV_MATRIX_NAME', `$1print_spmv_matrix')dnl
dnl
define(`PRINT_SPMV_MATRIX_HEAD', 
  `void PRINT_SPMV_MATRIX_NAME($1) ($1_array a, int n, 
     enum blas_order_type order, enum blas_uplo_type uplo)')dnl
dnl
define(`PRINT_SPMV_MATRIX_BODY', 
  `
{
  int row;
  DECLARE_VECTOR(x, $1_type)
  MALLOC_VECTOR(x, $1_type, n)
  
  for (row = 0; row < n; row++ ) {
        $1spmv_copy_row(order, uplo, n, a, x, row);
        $1print_vector(x, n, 1, NULL);
  }
  printf("\n");
  FREE_VECTOR(x, $1_type)
}
')dnl
dnl
define(`PRINT_SPMV_MATRIX', 
  `PRINT_SPMV_MATRIX_HEAD($1) {
   PRINT_SPMV_MATRIX_BODY($1)
  }')dnl
dnl
dnl
define(`PROTOTYPES', `dnl
SPMV_COPY_ROW_HEAD(s);
SPMV_COPY_ROW_HEAD(d);
SPMV_COPY_ROW_HEAD(c);
SPMV_COPY_ROW_HEAD(z);

SPMV_COMMIT_ROW_HEAD(s);
SPMV_COMMIT_ROW_HEAD(d);
SPMV_COMMIT_ROW_HEAD(c);
SPMV_COMMIT_ROW_HEAD(z);

SPMV_PACK_MATRIX_HEAD(s);
SPMV_PACK_MATRIX_HEAD(d);
SPMV_PACK_MATRIX_HEAD(c);
SPMV_PACK_MATRIX_HEAD(z);

PRINT_SPMV_MATRIX_HEAD(s);
PRINT_SPMV_MATRIX_HEAD(d);
PRINT_SPMV_MATRIX_HEAD(c);
PRINT_SPMV_MATRIX_HEAD(z);
')dnl
dnl
dnl
define(`SOURCE', `dnl
#include <stdio.h>
#include "blas_extended.h"
#include "blas_extended_test.h"

SPMV_COPY_ROW(s)
SPMV_COPY_ROW(d)
SPMV_COPY_ROW(c)
SPMV_COPY_ROW(z)

SPMV_COMMIT_ROW(s)
SPMV_COMMIT_ROW(d)
SPMV_COMMIT_ROW(c)
SPMV_COMMIT_ROW(z)

SPMV_PACK_MATRIX(s)
SPMV_PACK_MATRIX(d)
SPMV_PACK_MATRIX(c)
SPMV_PACK_MATRIX(z)

PRINT_SPMV_MATRIX(s)
PRINT_SPMV_MATRIX(d)
PRINT_SPMV_MATRIX(c)
PRINT_SPMV_MATRIX(z)
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
