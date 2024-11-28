dnl
dnl
include(cblas.m4)dnl
include(packed.m4)dnl
include(test-common.m4)dnl
dnl
dnl
dnl  HPMV_COPY_ROW
dnl     |
dnl     |-- HPMV_COPY_ROW_HEAD
dnl     |      |
dnl     |      |-- HPMV_COPY_ROW_NAME
dnl     |      |
dnl     |      |-- HPMV_COPY_ROW_PARAMS
dnl     |
dnl     |-- HPMV_COPY_ROW_BODY
dnl
dnl
dnl  Usage: 
dnl    HPMV_COPY_ROW($1)
dnl    HPMV_COPY_ROW_HEAD($1)
dnl    HPMV_COPY_ROW_NAME($1)
dnl    HPMV_COPY_ROW_PARAMS($1)
dnl    HPMV_COPY_ROW_COMMENT($1)
dnl    HPMV_COPY_ROW_BODY($1)
dnl
dnl    $1 is the type of the matrix.
dnl
dnl
define(`HPMV_COPY_ROW_NAME', `$1hpmv_copy_row')dnl
dnl
dnl
define(`HPMV_COPY_ROW_PARAMS', 
  `enum blas_order_type order, enum blas_uplo_type uplo, dnl
   int n, $1_array a, $1_array a_vec, int row')dnl
dnl
dnl
define(`HPMV_COPY_ROW_HEAD', 
  `void HPMV_COPY_ROW_NAME($1)(HPMV_COPY_ROW_PARAMS($1))')dnl
dnl
dnl
define(`HPMV_COPY_ROW_COMMENT', `
/*
 *  Copies the given row of packed hermitianmatrix a into the supplied vector.
 */')dnl
dnl
dnl
define(`HPMV_COPY_ROW_BODY',
`
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
                        if ( uplo == blas_upper ) {
                                CONJ_AUX(tmp, $1_type)
                        }
                        SET_VECTOR_ELEMENT(a_vec_i, i * inc, tmp, $1_type)
                        COLUMN_STRIDE(ind, n, i, 
                                blas_colmajor, blas_lower, inc)
                }

                ind = (row + ((2 * n - row - 1) * row) / 2) * inc;
                for(i = row; i < n; i++)
                {
                        GET_VECTOR_ELEMENT(tmp, a_i, ind, $1_type)
                        if ( uplo == blas_lower ) {
                                CONJ_AUX(tmp, $1_type)
                        }
                        if ( i == row ) {
                          ZERO_IMAG_PART(tmp, $1_type)
                        }
                        SET_VECTOR_ELEMENT(a_vec_i, i * inc, tmp, $1_type)
                        COLUMN_STRIDE(ind, n, i, 
                                blas_rowmajor, blas_upper, inc)
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
                        if ( uplo == blas_upper ) {
                                CONJ_AUX(tmp, $1_type)
                        }
                        SET_VECTOR_ELEMENT(a_vec_i, i * inc, tmp, $1_type)
                        COLUMN_STRIDE(ind, n, i, blas_rowmajor, blas_lower, inc)
                }

                ind = (row + (row * (row + 1)) / 2) * inc;
                for(i = row; i < n; i++) {
                        GET_VECTOR_ELEMENT(tmp, a_i, ind, $1_type)
                        if ( uplo == blas_lower ) {
                                CONJ_AUX(tmp, $1_type)
                        }
                        if ( i == row ) {
                          ZERO_IMAG_PART(tmp, $1_type);
                        }
                        SET_VECTOR_ELEMENT(a_vec_i, i * inc, tmp, $1_type)
                        COLUMN_STRIDE(ind, n, i, 
                                blas_colmajor, blas_upper, inc)
                }
        }
} 
')dnl
dnl
dnl
define(`HPMV_COPY_ROW', 
  `HPMV_COPY_ROW_HEAD($1)
   HPMV_COPY_ROW_COMMENT($1)
   HPMV_COPY_ROW_BODY($1)')dnl
dnl
dnl
dnl
dnl
dnl
dnl  HPMV_COMMIT_ROW
dnl     |
dnl     |-- HPMV_COMMIT_ROW_HEAD
dnl     |      |
dnl     |      |-- HPMV_COMMIT_ROW_NAME
dnl     |      |
dnl     |      |-- HPMV_COMMIT_ROW_PARAMS
dnl     |
dnl     |-- HPMV_COMMIT_ROW_BODY
dnl
dnl
dnl  Usage: 
dnl    HPMV_COMMIT_ROW($1)
dnl    HPMV_COMMIT_ROW_HEAD($1)
dnl    HPMV_COMMIT_ROW_NAME($1)
dnl    HPMV_COMMIT_ROW_PARAMS($1)
dnl    HPMV_COMMIT_ROW_COMMENT($1)
dnl    HPMV_COMMIT_ROW_BODY($1)
dnl
dnl    $1 is the type of the matrix.
dnl
dnl
define(`HPMV_COMMIT_ROW_NAME', `$1hpmv_commit_row')dnl
dnl
dnl
define(`HPMV_COMMIT_ROW_PARAMS', 
  `enum blas_order_type order, enum blas_uplo_type uplo, dnl
   int n, $1_array a, $1_array a_vec, int row')dnl
dnl
dnl
define(`HPMV_COMMIT_ROW_HEAD', 
  `void HPMV_COMMIT_ROW_NAME($1)(HPMV_COMMIT_ROW_PARAMS($1))')dnl
dnl
dnl
define(`HPMV_COMMIT_ROW_COMMENT', `
/*
 *  Commits the supplied vector a_vec to the given row of 
 *      packed hermitian matrix a.
 */')dnl
dnl
dnl
define(`HPMV_COMMIT_ROW_BODY',
`
{
        int i, ind, inc = 1;
        PTR_CAST(a, $1_type)
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
                        GET_VECTOR_ELEMENT(tmp, a_vec_i, i * inc, $1_type)
                        if ( uplo == blas_upper ) {
                                CONJ_AUX(tmp, $1_type)
                        }
                        SET_VECTOR_ELEMENT(a_i, ind, tmp, $1_type)
                        COLUMN_STRIDE(ind, n, i, 
                                blas_colmajor, blas_lower, inc)
                }

                ind = (row + ((2 * n - row - 1) * row) / 2) * inc;
                for(i = row; i < n; i++)
                {
                        GET_VECTOR_ELEMENT(tmp, a_vec_i, i * inc, $1_type)
                        if ( uplo == blas_lower ) {
                                CONJ_AUX(tmp, $1_type)
                        }
                        SET_VECTOR_ELEMENT(a_i, ind, tmp, $1_type)
                        COLUMN_STRIDE(ind, n, i, 
                                blas_rowmajor, blas_upper, inc)
                }
        }
        else
        {
                /* Pretend it is rowmajor/lower.  We can do this in the
                   symmetric case. */
                ind = row * (row + 1) * inc / 2;
                for(i = 0; i < row; i++)
                {
                        GET_VECTOR_ELEMENT(tmp, a_vec_i, i * inc, $1_type)
                        if ( uplo == blas_upper ) {
                                CONJ_AUX(tmp, $1_type)
                        }
                        SET_VECTOR_ELEMENT(a_i, ind, tmp, $1_type)
                        COLUMN_STRIDE(ind, n, i, 
                                blas_rowmajor, blas_lower, inc)
                }

                ind = (row + (row * (row + 1)) / 2) * inc;
                for(i = row; i < n; i++)
                {
                        GET_VECTOR_ELEMENT(tmp, a_vec_i, i * inc, $1_type)
                        if ( uplo == blas_lower ) {
                                CONJ_AUX(tmp, $1_type)
                        }
                        SET_VECTOR_ELEMENT(a_i, ind, tmp, $1_type)
                        COLUMN_STRIDE(ind, n, i, 
                                blas_colmajor, blas_upper, inc)
                }
        }
} 
')dnl
dnl
dnl
define(`HPMV_COMMIT_ROW', 
  `HPMV_COMMIT_ROW_HEAD($1)
   HPMV_COMMIT_ROW_COMMENT($1)
   HPMV_COMMIT_ROW_BODY($1)')dnl
dnl
dnl
define(`PRINT_HPMV_MATRIX_NAME', `$1print_hpmv_matrix')dnl
dnl
define(`PRINT_HPMV_MATRIX_HEAD', 
  `void PRINT_HPMV_MATRIX_NAME($1) ($1_array a, int n, 
     enum blas_order_type order, enum blas_uplo_type uplo)')dnl
dnl
define(`PRINT_HPMV_MATRIX_BODY', 
  `
{
  int row;
  DECLARE_VECTOR(x, $1_type)
  MALLOC_VECTOR(x, $1_type, n)
  
  for (row = 0; row < n; row++ ) {
        $1hpmv_copy_row(order, uplo, n, a, x, row);
        $1print_vector(x, n, 1, NULL);
  }
  printf("\n");
  FREE_VECTOR(x, $1_type)
}
')dnl
dnl
define(`PRINT_HPMV_MATRIX', 
  `PRINT_HPMV_MATRIX_HEAD($1) {
   PRINT_HPMV_MATRIX_BODY($1)
  }')dnl
dnl
dnl
dnl
dnl  HPMV_PACK_MATRIX
dnl     |
dnl     |-- HPMV_PACK_MATRIX_HEAD
dnl     |      |
dnl     |      |-- HPMV_PACK_MATRIX_NAME
dnl     |      |
dnl     |      |-- HPMV_PACK_MATRIX_PARAMS
dnl     |
dnl     |-- HPMV_PACK_MATRIX_BODY
dnl
dnl
dnl  Usage: 
dnl    HPMV_PACK_MATRIX($1)
dnl    HPMV_PACK_MATRIX_HEAD($1)
dnl    HPMV_PACK_MATRIX_NAME($1)
dnl    HPMV_PACK_MATRIX_PARAMS($1)
dnl    HPMV_PACK_MATRIX_COMMENT($1)
dnl    HPMV_PACK_MATRIX_BODY($1)
dnl
dnl    $1 is the type of the matrix.
dnl
dnl
define(`HPMV_PACK_MATRIX_NAME', `$1hpmv_pack_matrix')dnl
dnl
dnl
define(`HPMV_PACK_MATRIX_PARAMS', 
  `enum blas_order_type order, enum blas_uplo_type uplo, dnl 
   int n, $1_array a_packed, $1_array a_full, int lda')dnl
dnl
dnl
define(`HPMV_PACK_MATRIX_HEAD', 
  `void HPMV_PACK_MATRIX_NAME($1)(HPMV_PACK_MATRIX_PARAMS($1))')dnl
dnl
dnl
define(`HPMV_PACK_MATRIX_COMMENT', `
/*
 *  Packs the he matrix a_full into packed form a.
 */')dnl
dnl
dnl
define(`HPMV_PACK_MATRIX_BODY',
`{
        int row;
        DECLARE_VECTOR(a_row, $1_type);
        
        MALLOC_VECTOR(a_row, $1_type, n);       
        for (row = 0; row < n; row ++) {
                $1he_copy_row(order, uplo, blas_left_side, n, a_full, lda, a_row, row);
                $1hpmv_commit_row(order, uplo, n, a_packed,
                        a_row, row); 
        }

        FREE_VECTOR(a_row, $1_type)
}')dnl
dnl
dnl
define(`HPMV_PACK_MATRIX', 
  `HPMV_PACK_MATRIX_HEAD($1)
   HPMV_PACK_MATRIX_COMMENT($1)
   HPMV_PACK_MATRIX_BODY($1)')dnl
dnl
dnl
define(`PROTOTYPES', `dnl
HPMV_COPY_ROW_HEAD(c);
HPMV_COPY_ROW_HEAD(z);

HPMV_COMMIT_ROW_HEAD(c);
HPMV_COMMIT_ROW_HEAD(z);

HPMV_PACK_MATRIX_HEAD(c);
HPMV_PACK_MATRIX_HEAD(z);

PRINT_HPMV_MATRIX_HEAD(c);
PRINT_HPMV_MATRIX_HEAD(z);
')dnl
dnl
dnl
define(`SOURCE', `dnl
#include <stdio.h>
#include "blas_extended.h"
#include "blas_extended_test.h"

HPMV_COPY_ROW(c)
HPMV_COPY_ROW(z)

HPMV_COMMIT_ROW(c)
HPMV_COMMIT_ROW(z)

HPMV_PACK_MATRIX(c)
HPMV_PACK_MATRIX(z)

PRINT_HPMV_MATRIX(c)
PRINT_HPMV_MATRIX(z)
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
