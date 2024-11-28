dnl **********************************************************************
dnl * Generate trsv_prepare, trsv_copy, and trsv_commit                  *
dnl **********************************************************************
dnl
dnl
include(cblas.m4)dnl
dnl
dnl
dnl -------------------------------------------------------------------------
dnl Usage: TRSV_COMMIT(T_typeltr)
dnl        copy the vector y to T
dnl -------------------------------------------------------------------------
define(`TRSV_COMMIT_HEAD', 
  `void $1trsv_commit(enum blas_order_type order, enum blas_uplo_type uplo, dnl
       enum blas_trans_type trans, int length, $1_array T, int lda, dnl
       const $1_array y, int row)')dnl
dnl
dnl
define(`TRSV_COMMIT',
`TRSV_COMMIT_HEAD($1)
/*
 * Purpose
 * =======
 *
 * Copy y to T
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
 *              Dimension of AP and the length of vector x
 *
 * T            (input) $1_array
 *              The triangular matrix T
 *
 * lda          (input) int
 *              Leading dimension
 *
 * y            (output) $1_array
 *              The vector y
 *
 * row          (input) int
 *              The row to be copyied to y
 */
{
        int i; 
        DECLARE(tmp, $1_type)

        if ((order == blas_rowmajor && uplo == blas_upper && trans == blas_no_trans) ||
            (order == blas_colmajor && uplo == blas_lower && trans != blas_no_trans)) { 
                for(i=0; i<length; i++){
                        GET_VECTOR_ELEMENT(tmp, y, i, $1_type)
                        SET_VECTOR_ELEMENT(T, (row*lda) + (row+i+1), tmp, $1_type)
                }
        }
        else
        if ((order == blas_rowmajor && uplo == blas_lower && trans == blas_no_trans) ||
            (order == blas_colmajor && uplo == blas_upper && trans != blas_no_trans)) { 
                for(i=0; i<length; i++){
                        GET_VECTOR_ELEMENT(tmp, y, i, $1_type)
                        SET_VECTOR_ELEMENT(T, row*lda + i, tmp, $1_type)
                }
        }
        else
        if ((order == blas_rowmajor && uplo == blas_lower && trans != blas_no_trans) ||
            (order == blas_colmajor && uplo == blas_upper && trans == blas_no_trans)) {
                for(i=0; i<length; i++){
                        GET_VECTOR_ELEMENT(tmp, y, i, $1_type)
                        SET_VECTOR_ELEMENT(T, (row+i+1)*lda + row, tmp, $1_type)
                }
        }
        else
        if ((order == blas_rowmajor && uplo == blas_upper && trans != blas_no_trans) ||
            (order == blas_colmajor && uplo == blas_lower && trans == blas_no_trans)) {
                for(i=0; i<length; i++){
                        GET_VECTOR_ELEMENT(tmp, y, i, $1_type)
                        SET_VECTOR_ELEMENT(T, i*lda + row, tmp, $1_type)
                }
        }
}
')dnl
dnl
dnl
define(`PROTOTYPES', `dnl
TRSV_COMMIT_HEAD(s);
TRSV_COMMIT_HEAD(d);
')dnl
dnl
dnl
define(`SOURCE', `dnl
#include "blas_extended.h"
#include "blas_extended_test.h"

TRSV_COMMIT(s)
TRSV_COMMIT(d)
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
