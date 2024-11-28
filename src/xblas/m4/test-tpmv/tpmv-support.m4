include(cblas.m4)dnl
include(packed.m4)dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: TPMV_COPY_ROW(tp_typeltr)
dnl        copy a row to the vector a_vec
dnl ----------------------------------------------------------------------
define(`TPMV_COPY_ROW_HEAD', 
  `void $1tpmv_copy_row(enum blas_order_type order, dnl
       enum blas_uplo_type uplo, enum blas_trans_type trans, dnl
       int n, const $1_array a, $1_array a_vec, int row)')dnl
dnl
dnl
define(`TPMV_COPY_ROW',
`TPMV_COPY_ROW_HEAD($1)
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
 *              Dimension of tp and the length of vector x
 *
 * a           (input) $1_array
 *
 * a_vec            (input) $1_array
 *
 * row          (input) int
 *
 */
{
        int i,j, ind, stride, inc = 1;
        PTR_CAST(a, $1_type, `const')
        PTR_CAST(a_vec, $1_type)
        DECLARE(tmp, $1_type)
        DECLARE(tmp2, $1_type)
        INC_ADJUST(inc, $1_type)
        ZERO(tmp2, $1_type)
        for (j = 0; j < n; j++){
                SET_VECTOR_ELEMENT(a_vec_i, j*inc, tmp2, $1_type)
        }               

        if(((order == blas_rowmajor)&&(uplo == blas_upper)&&(trans != blas_no_trans))||
           ((order == blas_colmajor)&&(uplo == blas_lower)&&(trans == blas_no_trans)))
        {
                /* colmajor/lower.*/
                ind = row * inc;
                stride = (n-1)*inc;
                for(i = 0; i <= row; i++)
                {
                        GET_VECTOR_ELEMENT(tmp, a_i, ind, $1_type)
                        SET_VECTOR_ELEMENT(a_vec_i, i * inc, tmp, $1_type)
                        ind += stride;
                        stride-=inc;                    

                }

        }
        else if(((order == blas_rowmajor)&&(uplo == blas_lower)&&(trans == blas_no_trans))||
                ((order == blas_colmajor)&&(uplo == blas_upper)&&(trans != blas_no_trans)))
        {
                /* Pretend it is rowmajor/lower.*/

                ind = row * (row + 1) * inc / 2;
                for(i = 0; i <= row; i++)
                {
                        GET_VECTOR_ELEMENT(tmp, a_i, ind, $1_type)
                        SET_VECTOR_ELEMENT(a_vec_i, i * inc, tmp, $1_type)
                        ind+=inc;

                }


        }
        else if(((order == blas_rowmajor)&&(uplo == blas_upper)&&(trans == blas_no_trans))||
                ((order == blas_colmajor)&&(uplo == blas_lower)&&(trans != blas_no_trans)))
        {
                /* Pretend it is rowmajor/upper.*/


                ind = (row + ((2 * n - row - 1) * row) / 2) * inc;
                for(i = row; i < n; i++)
                {
                        GET_VECTOR_ELEMENT(tmp, a_i, ind, $1_type)
                        SET_VECTOR_ELEMENT(a_vec_i, i * inc, tmp, $1_type)
                        ind+=inc;

                }
        }
        else
        {
                /* Pretend it is colmajor/upper.*/

                ind = (row + (row * (row + 1)) / 2) * inc;
                stride = (row + 1)*inc;
                for(i = row; i < n; i++)
                {
                        GET_VECTOR_ELEMENT(tmp, a_i, ind, $1_type)
                        SET_VECTOR_ELEMENT(a_vec_i, i * inc, tmp, $1_type)
                        ind +=stride;
                        stride+=inc;

                }
        }


}
')dnl
dnl
dnl
dnl -------------------------------------------------------------------------
dnl Usage: TPMV_COMMIT_ROW(tp_typeltr)
dnl        copy the vector a_vec to a
dnl -------------------------------------------------------------------------
define(`TPMV_COMMIT_ROW_HEAD', 
  `void $1tpmv_commit_row(enum blas_order_type order, dnl
       enum blas_uplo_type uplo, enum blas_trans_type trans, dnl
       int n, $1_array a, const $1_array a_vec, int row)')dnl
dnl
dnl
define(`TPMV_COMMIT_ROW',
`TPMV_COMMIT_ROW_HEAD($1)
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
        int i, stride, ind, inc = 1;
        PTR_CAST(a, $1_type)
        PTR_CAST(a_vec, $1_type, `const')
        DECLARE(tmp, $1_type)

        INC_ADJUST(inc, $1_type)



        if(((order == blas_rowmajor)&&(uplo == blas_upper)&&(trans != blas_no_trans))||
           ((order == blas_colmajor)&&(uplo == blas_lower)&&(trans == blas_no_trans)))
        {
                /* colmajor/lower.*/
                stride = (n-1)*inc;
                ind = row * inc;
                for(i = 0; i <= row; i++)
                {
                        GET_VECTOR_ELEMENT(tmp, a_vec_i, i * inc, $1_type)
                        SET_VECTOR_ELEMENT(a_i, ind, tmp, $1_type)
                        ind += stride;
                        stride-=inc;


                }

        }
        else if(((order == blas_rowmajor)&&(uplo == blas_lower)&&(trans == blas_no_trans))||
                ((order == blas_colmajor)&&(uplo == blas_upper)&&(trans != blas_no_trans)))
        {
                /* Pretend it is rowmajor/lower.*/

                ind = row * (row + 1) * inc / 2;
                for(i = 0; i <= row; i++)
                {
                        GET_VECTOR_ELEMENT(tmp, a_vec_i, i * inc, $1_type)
                        SET_VECTOR_ELEMENT(a_i, ind, tmp, $1_type)
                        ind+=inc;

                }


        }
        else if(((order == blas_rowmajor)&&(uplo == blas_upper)&&(trans == blas_no_trans))||
                ((order == blas_colmajor)&&(uplo == blas_lower)&&(trans != blas_no_trans)))
        {
                /* Pretend it is rowmajor/upper.*/


                ind = (row + ((2 * n - row - 1) * row) / 2) * inc;
                for(i = row; i < n; i++)
                {
                        GET_VECTOR_ELEMENT(tmp, a_vec_i, i * inc, $1_type)
                        SET_VECTOR_ELEMENT(a_i, ind, tmp, $1_type)
                        ind+=inc;

                }
        }
        else
        {
                /* Pretend it is colmajor/upper.*/

                ind = (row + (row * (row + 1)) / 2) * inc;
                stride = (row+1)*inc;
                for(i = row; i < n; i++)
                {
                        GET_VECTOR_ELEMENT(tmp, a_vec_i, i * inc, $1_type)
                        SET_VECTOR_ELEMENT(a_i, ind, tmp, $1_type)
                        ind += stride;
                        stride+=inc;

                }
        }
}
')dnl
dnl
dnl
define(`PROTOTYPES', `dnl
TPMV_COPY_ROW_HEAD(s);
TPMV_COPY_ROW_HEAD(d);
TPMV_COPY_ROW_HEAD(c);
TPMV_COPY_ROW_HEAD(z);

TPMV_COMMIT_ROW_HEAD(s);
TPMV_COMMIT_ROW_HEAD(d);
TPMV_COMMIT_ROW_HEAD(c);
TPMV_COMMIT_ROW_HEAD(z);
')dnl
dnl
dnl
define(`SOURCE', `dnl
#include "blas_extended.h"
#include "blas_extended_test.h"

TPMV_COPY_ROW(s)
TPMV_COPY_ROW(d)
TPMV_COPY_ROW(c)
TPMV_COPY_ROW(z)

TPMV_COMMIT_ROW(s)
TPMV_COMMIT_ROW(d)
TPMV_COMMIT_ROW(c)
TPMV_COMMIT_ROW(z)
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
