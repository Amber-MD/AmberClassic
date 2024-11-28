dnl **********************************************************************
dnl * Generate tbsv_copy, and tbsv_commit                  *
dnl **********************************************************************
dnl
dnl
include(cblas.m4)dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: TBSV_COPY(T_typeltr)
dnl        copy a row to the vector y
dnl ----------------------------------------------------------------------
define(`TBSV_COPY_HEAD', 
  `void $1tbsv_copy(enum blas_order_type order, enum blas_uplo_type uplo, dnl
       enum blas_trans_type trans, int n, int k, const $1_array T, dnl
       int ldt, $1_array y, int row)')dnl
dnl
dnl
define(`TBSV_COPY',
`TBSV_COPY_HEAD($1)
/*
 * Purpose
 * =======
 *
 * Copy a row from T to y
 *
 *      NOTE: this function just encapsulates $1gbmv_copy.
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
 * k            (input) int
 *              Number of super(sub)diagonals of T. 
 *
 * T            (input) $1_array
 *              The triangular matrix T
 *
 * ldt          (input) int
 *              Leading dimension
 *
 * y            (output) $1_array
 *              The vector y
 *
 * row          (input) int
 *              The row to be copyied to y
 */
{
        enum blas_trans_type new_trans;
        int conj=0;
        int kl, ku;
        IF_COMPLEX($1_type, `PTR_CAST(y, $1_type)')

        if(uplo == blas_upper) {
                ku = k;
                kl = 0; 
        } else {
                ku = 0;
                kl = k;
        }
        

        if(trans == blas_no_trans) {
                new_trans = blas_no_trans;
        } else if(trans == blas_conj_trans) {
                new_trans = blas_trans;
                conj = 1;
        } else if (trans == blas_trans) {
                new_trans = blas_trans;
        } else {
                /* conj, no trans */
                new_trans = blas_no_trans;
                conj = 1;
        }
        $1gbmv_copy(order, new_trans, n, n, kl, ku,
                T, ldt, y, row);

        IF_COMPLEX($1_type, 
        `if(conj) {
                DECLARE(y_elem, $1_type)
                int i, incyi, ni;

                incyi = 1;
                ni = n;
                INC_ADJUST(ni, $1_type)
                INC_ADJUST(incyi, $1_type)
                for(i=0; i < ni; i += incyi) {
                        GET_VECTOR_ELEMENT(y_elem, y_i, `i', $1_type)
                        CONJ_AUX(y_elem, $1_type)
                        SET_VECTOR_ELEMENT(y_i, `i', y_elem, $1_type)
                } 
        }')             
        
}
')dnl
dnl
dnl
dnl -------------------------------------------------------------------------
dnl Usage: TBSV_COMMIT(T_typeltr)
dnl        copy the vector y to T
dnl -------------------------------------------------------------------------
define(`TBSV_COMMIT_HEAD', 
  `void $1tbsv_commit(enum blas_order_type order, dnl
       enum blas_uplo_type uplo, enum blas_trans_type trans, dnl
       int n, int k, $1_array T, int ldt, $1_array y, int row)')dnl
dnl
dnl
define(`TBSV_COMMIT',
`TBSV_COMMIT_HEAD($1)
/*
 * Purpose
 * =======
 *
 * Copy y to T
 *
 *      NOTE: this function just encapsulates $1gbmv_commit.
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
 * k            (input) int
 *              Number of super(sub)diagonals of T
 *
 * T            (input) $1_array
 *              The triangular matrix T
 *
 * ldt          (input) int
 *              Leading dimension
 *
 * y            (output) $1_array
 *              The vector y
 *
 * row          (input) int
 *              The row to be copyied to y
 */
{
        enum blas_trans_type new_trans;
        int conj=0;
        int kl, ku;
        IF_COMPLEX($1_type, `PTR_CAST(y, $1_type)')

        if(uplo == blas_upper) {
                ku = k;
                kl = 0; 
        } else {
                ku = 0;
                kl = k;
        }
        

        if(trans == blas_no_trans) {
                new_trans = blas_no_trans;
        } else if(trans == blas_conj_trans) {
                new_trans = blas_trans;
                conj = 1;
        } else if (trans == blas_trans) {
                new_trans = blas_trans;
        } else {
                /* conj, no trans */
                new_trans = blas_no_trans;
                conj = 1;
        }

        IF_COMPLEX($1_type, 
        `if(conj) {
                DECLARE(y_elem, $1_type)
                int i, incyi, ni;

                incyi = 1;
                ni = n;
                INC_ADJUST(ni, $1_type)
                INC_ADJUST(incyi, $1_type)
                for(i=0; i < ni; i += incyi) {
                        GET_VECTOR_ELEMENT(y_elem, y_i, `i', $1_type)
                        CONJ_AUX(y_elem, $1_type)
                        SET_VECTOR_ELEMENT(y_i, `i', y_elem, $1_type)
                } 
        }')             

        $1gbmv_commit(order, new_trans, n, n, kl, ku,
                T, ldt, y, row);


        IF_COMPLEX($1_type,
        `if(conj) {
                /* now conjugate back - leave original */
                DECLARE(y_elem, $1_type)
                int i, incyi, ni;

                incyi = 1;
                ni = n;
                INC_ADJUST(ni, $1_type)
                INC_ADJUST(incyi, $1_type)
                for(i=0; i < ni; i += incyi) {
                        GET_VECTOR_ELEMENT(y_elem, y_i, `i', $1_type)
                        CONJ_AUX(y_elem, $1_type)
                        SET_VECTOR_ELEMENT(y_i, `i', y_elem, $1_type)
                } 
        }')             

}
')dnl
dnl 
dnl
dnl -------------------------------------------------------------------------
dnl Usage: TBSV_PRINT_MATRIX(T_typeltr)
dnl        print T
dnl -------------------------------------------------------------------------
define(`TBSV_PRINT_MATRIX_HEAD', 
  `void $1print_tbsv_matrix($1_array T, int n, int k, int ldt, dnl
       enum blas_order_type order, enum blas_uplo_type uplo, dnl
       enum blas_trans_type trans)')dnl
dnl
dnl
define(`TBSV_PRINT_MATRIX',
`TBSV_PRINT_MATRIX_HEAD($1)
/*
 * Purpose
 * =======
 *
 * Print T
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
 * k            (input) int
 *              Number of super(sub)diagonals of T
 *
 * T            (input) $1_array
 *              The triangular matrix T
 *
 * ldt          (input) int
 *              Leading dimension
 *
 */
{
        int i;
        DECLARE_VECTOR(T_row, $1_type)
        MALLOC_VECTOR(T_row, $1_type, n)

        for(i=0; i< n; i++) {
                $1tbsv_copy(order, uplo, trans, n, k, T, ldt, T_row, i);
                $1print_vector(T_row, n, 1, NULL);
        }
        printf("\n");
        FREE_VECTOR(T_row, $1_type)
}
')dnl
dnl
dnl
define(`PROTOTYPES', `dnl
TBSV_COPY_HEAD(s);
TBSV_COPY_HEAD(d);
TBSV_COPY_HEAD(c);
TBSV_COPY_HEAD(z);

TBSV_COMMIT_HEAD(s);
TBSV_COMMIT_HEAD(d);
TBSV_COMMIT_HEAD(c);
TBSV_COMMIT_HEAD(z);

TBSV_PRINT_MATRIX_HEAD(s);
TBSV_PRINT_MATRIX_HEAD(d);
TBSV_PRINT_MATRIX_HEAD(c);
TBSV_PRINT_MATRIX_HEAD(z);
')dnl
dnl
dnl
define(`SOURCE', `dnl
#include <stdio.h>
#include "blas_extended.h"
#include "blas_extended_test.h"

TBSV_COPY(s)
TBSV_COMMIT(s)
TBSV_COPY(d)
TBSV_COMMIT(d)
TBSV_COPY(c)
TBSV_COMMIT(c)
TBSV_COPY(z)
TBSV_COMMIT(z)

TBSV_PRINT_MATRIX(s)
TBSV_PRINT_MATRIX(d)
TBSV_PRINT_MATRIX(c)
TBSV_PRINT_MATRIX(z)
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
