dnl **********************************************************************
dnl * Generates alpha, tp, and x, where tp is a triangular packed matrix;* 
dnl * and computes r_true.                                               *
dnl **********************************************************************
dnl
dnl
include(cblas.m4)dnl
include(test-common.m4)dnl
dnl
dnl
define(`TPMV_TESTGEN_COMMENT', `
/*
 * Purpose
 * =======
 *
 * Generates alpha, tp and x, where tp is a triangular matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of tp; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether tp is upper or lower
 *
 * trans        (input) blas_trans_type
 *              No trans, trans, conj trans
 *
 * diag         (input) blas_diag_type
 *              non unit, unit
 *
 * n            (input) int
 *              Dimension of tp and the length of vector x
 *
 * alpha        (input/output) $1_array
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * tp            (output) $2_array
 *
 * x            (input/output) $1_array
 *
 * seed         (input/output) int
 *
 * HEAD(r_true)     (output) double*
 *              The leading part of the truth in double-double.
 *
 * TAIL(r_true)     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */')
dnl
dnl
dnl
dnl ---------------------------------------------------------------------
dnl Usage: TPMV_TESTGEN(ax_typeltr, tp_typeltr)
dnl        produce tpmv_prepare signature
dnl ---------------------------------------------------------------------
define(`TPMV_TESTGEN_NAME', `BLAS_$1tpmv`'ifelse(`$1', `$2', `', `_$2')_testgen')dnl
dnl
dnl
define(`TPMV_TESTGEN_HEAD', 
  `void TPMV_TESTGEN_NAME($1, $2)(int norm, enum blas_order_type order, dnl
       enum blas_uplo_type uplo, enum blas_trans_type trans, dnl
       enum blas_diag_type diag, int n, $1_array alpha, int alpha_flag, dnl
       $2_array tp, $1_array x,  int *seed, dnl
       double *HEAD(r_true), double *TAIL(r_true))')dnl
dnl
dnl
define(`TPMV_TESTGEN', `
TPMV_TESTGEN_HEAD($1, $2)
TPMV_TESTGEN_COMMENT($1, $2)
TPMV_TESTGEN_BODY($1, $2)')dnl
dnl
dnl
dnl ---------------------------------------------------------------------
dnl Usage: TPMV_TESTGEN_BODY(ax_typeltr, tp_typeltr)
dnl        
dnl ---------------------------------------------------------------------
define(`TPMV_TESTGEN_BODY',
`{
  int n_fix2 = 0;
  int n_mix = 0;
  int i;
  int j;
  int length;   
  DECLARE_VECTOR(temp, $2_type)
  DECLARE(beta_zero_fake, $1_type)
  int beta_flag = 1;
  DECLARE(r_zero_fake, $1_type)
  DECLARE(one_this, $2_type)
  DECLARE(tp_swapTemp1, $2_type)
  DECLARE(tp_swapTemp2, $2_type)
  DECLARE(x_swapTemp1, $1_type)
  DECLARE(x_swapTemp2, $1_type)
  DECLARE(r_true_swapTemp1, R_TRUE_TYPE($1))    
  DECLARE(r_true_swapTemp2, R_TRUE_TYPE($1))    

  PTR_CAST(x, $1_type)
  PTR_CAST(alpha, $1_type)
  int inctp = 1;
  int inc_index = 1;

  ZERO(beta_zero_fake, $1_type)
  ZERO(r_zero_fake, $1_type)
  ONE(one_this, $2_type)

  INC_ADJUST(inctp, $2_type)
  INC_ADJUST(inc_index, $1_type)

  MALLOC_VECTOR(temp, $2_type, n*2*inc_index)
  for (i=0; i<n*2*inc_index; i+=inc_index ) {
        SET_ZERO_VECTOR_ELEMENT(temp,i,$2_type)
  }
  for(i=0; i<(n-1+n-1+1)*n*2*inctp; i++){
    ifelse(`$2', `c', `((float*)tp)[i]=0.0;',
           `$2', `z', `((double*)tp)[i]=0.0;',
                      `tp[i]=0.0;')
  }

  /* calling dot_testgen n time. in each iteration, one row of tp, and one 
     element of x are produced. */

  if (alpha_flag == 0 && diag == blas_unit_diag) {
    RAND_PTR(alpha_i, $1_type)
    alpha_flag = 1;
  }
  
  for(i=0; i<n; i++){

    if (i == 0){
      n_fix2 = 0;
      n_mix  = 0;
    }
    else {
       n_mix  = i;

      /* from now on, fix alpha */
      alpha_flag = 1;
    }

    for (j=0; j<n*inc_index; j+=inc_index ) {
        SET_ZERO_VECTOR_ELEMENT(temp,j,$2_type)
    }
    length = i + 1;
        
    if (diag == blas_unit_diag){
        length = i;

                        
    DOT_TESTGEN_NAME($1, $1, $2)(length, n_fix2, n_mix, norm, blas_no_conj,
                        alpha, alpha_flag,
                        alpha, beta_flag, x, temp, seed, PASS_BY_REF(r_zero_fake, $1_type),
                        &((double*)HEAD(r_true))[i*inc_index], &((double*)TAIL(r_true))[i*inc_index]);
   }else{       
    DOT_TESTGEN_NAME($1, $1, $2)(length, n_fix2, n_mix, norm, blas_no_conj,
                        alpha, alpha_flag,
                        PASS_BY_REF(beta_zero_fake, $1_type), beta_flag, x, temp, seed, PASS_BY_REF(r_zero_fake, $1_type),
                        &((double*)HEAD(r_true))[i*inc_index], &((double*)TAIL(r_true))[i*inc_index]);
  }

/* one type is $1 */
  if (diag == blas_unit_diag){
        ASSIGN_THIS(temp, one_this, length*inctp, $2)
        ASSIGN_THIS(x_i, r_zero_fake, length*inc_index, $1)
  }
  ifelse(`$2_type', `c', 
        `if (trans == blas_conj_trans){
           for(j=0; j<n_mix; j++){      
                temp[j*inctp+1] = -temp[j*inctp+1];
           }
         }',
         `$2_type', `z',
        `if (trans == blas_conj_trans){
           for(j=0; j<n_mix; j++){      
                temp[j*inctp+1] = -temp[j*inctp+1];
           }
         }', `')

    /* copy temp to tp */
    if (((order == blas_rowmajor)&&(trans == blas_no_trans)&&(uplo == blas_upper))
        ||
        ((order == blas_colmajor)&&(trans != blas_no_trans)&&(uplo == blas_lower))
        ||
        ((order == blas_rowmajor)&&(trans != blas_no_trans)&&(uplo == blas_lower))
        ||
        ((order == blas_colmajor)&&(trans == blas_no_trans)&&(uplo == blas_upper)))
                {
          for(j = 0; j < n/2; j++) {
                GET_VECTOR_ELEMENT(tp_swapTemp1, temp, (n-1-j)*inctp, $2_type)
                GET_VECTOR_ELEMENT(tp_swapTemp2, temp, j*inctp, $2_type)
                SET_VECTOR_ELEMENT(temp, (n-1-j)*inctp, tp_swapTemp2, $2_type)
                SET_VECTOR_ELEMENT(temp, j*inctp, tp_swapTemp1, $2_type)
          }

          $2tpmv_commit_row(order, uplo, trans, n, tp, temp, (n-1-i));
        }
        else if(((order == blas_rowmajor)&&(trans == blas_no_trans)&&(uplo == blas_lower))
                ||
                ((order == blas_colmajor)&&(trans != blas_no_trans)&&(uplo == blas_upper))
                ||
                ((order == blas_rowmajor)&&(trans != blas_no_trans)&&(uplo == blas_upper))
                ||
                ((order == blas_colmajor)&&(trans == blas_no_trans)&&(uplo == blas_lower)))
        {
          $2tpmv_commit_row(order, uplo, trans, n, tp, temp, i);
        } 

  }  

         /* The upper triangular matrices generate x backwards, so need to rearrange the elements */
         if (((order == blas_rowmajor)&&(trans == blas_no_trans)&&(uplo == blas_upper))
                ||
                ((order == blas_colmajor)&&(trans != blas_no_trans)&&(uplo == blas_lower))
                ||
                ((order == blas_rowmajor)&&(trans != blas_no_trans)&&(uplo == blas_lower))
                ||
                ((order == blas_colmajor)&&(trans == blas_no_trans)&&(uplo == blas_upper)))
            {
              for(j = 0; j < n/2; j++) {
                GET_VECTOR_ELEMENT(x_swapTemp1, x_i, (n-1-j)*inc_index, $1_type)
                GET_VECTOR_ELEMENT(x_swapTemp2, x_i, j*inc_index, $1_type)
                SET_VECTOR_ELEMENT(x_i, (n-1-j)*inc_index, x_swapTemp2, $1_type)
                SET_VECTOR_ELEMENT(x_i, j*inc_index, x_swapTemp1, $1_type)

                GET_ARRAY_ELEMENT(r_true_swapTemp1, R_TRUE_TYPE($1), r_true, R_TRUE_TYPE($1), (n-1-j)*inc_index)
                GET_ARRAY_ELEMENT(r_true_swapTemp2, R_TRUE_TYPE($1), r_true, R_TRUE_TYPE($1), j*inc_index)
                SET_ARRAY_ELEMENT(r_true_swapTemp2,R_TRUE_TYPE($1),r_true, R_TRUE_TYPE($1),(n-1-j)*inc_index) 
                SET_ARRAY_ELEMENT(r_true_swapTemp1,R_TRUE_TYPE($1),r_true, R_TRUE_TYPE($1), j*inc_index)
              }
            }

  FREE_VECTOR(temp, $2_type)
}')dnl
dnl
dnl
define(`R_TRUE_TYPE', `ifelse(`$1', `s', real_E,
        `$1', `d', real_E,
        complex_E)')dnl
dnl
dnl
define(`ASSIGN_THIS',
        `ifelse(`$4', `s', `$1[$3] = $2;',
                `$4', `d', `$1[$3] = $2;',
                `$1[$3] = $2[0]; $1[$3+1] = $2[1];')')dnl
dnl
dnl
define(`PROTOTYPES', `dnl
TPMV_TESTGEN_HEAD(s, s);
TPMV_TESTGEN_HEAD(d, d);
TPMV_TESTGEN_HEAD(c, c);
TPMV_TESTGEN_HEAD(z, z);
TPMV_TESTGEN_HEAD(d, s);
TPMV_TESTGEN_HEAD(z, c);
TPMV_TESTGEN_HEAD(c, s);
TPMV_TESTGEN_HEAD(z, d);
')dnl
dnl
dnl
define(`SOURCE', `dnl
#include "blas_extended.h"
#include "blas_extended_test.h"

TPMV_TESTGEN(s, s)
TPMV_TESTGEN(d, d)
TPMV_TESTGEN(c, c)
TPMV_TESTGEN(z, z)
TPMV_TESTGEN(d, s)
TPMV_TESTGEN(z, c)
TPMV_TESTGEN(c, s)
TPMV_TESTGEN(z, d)
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
