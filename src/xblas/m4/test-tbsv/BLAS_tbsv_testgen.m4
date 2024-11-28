dnl *************************************************************************
dnl * Generates alpha, T, x, and y, where T is a triangular Banded matrix;  *
dnl * and computes r_true.                                                  *
dnl *************************************************************************
dnl
dnl
include(cblas.m4)dnl
include(test-common.m4)dnl
dnl
dnl
define(`TBSV_TESTGEN_COMMENT', `
/*
 * Purpose
 * =======
 *
 * Generates alpha, x and T, where T is a triangular Banded matrix; 
 *      and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of T; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether T is upper or lower
 *
 * trans        (input) blas_trans_type
 *              No trans, trans, conj trans
 *
 * diag         (input) blas_diag_type
 *              non unit, unit
 *
 * n            (input) int
 *              Dimension of T and the length of vector x
 *
 * k            (input) int
 *              Number of super(sub)diagonals of T
 *
 * alpha        (input/output) $1_array
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * T            (output) $2_array
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
 * row          (input) int
 *              The true row being generated
 *
 * prec         (input) blas_prec_type
 *              single, double, or extra precision   
 *
 */')dnl
dnl
dnl
dnl
define(`TBSV_TESTGEN_NAME', `ifelse(`$1&&$2', `$1&&$1', 
  `BLAS_$1tbsv_testgen', `BLAS_$1tbsv_$2_testgen')')dnl
dnl
dnl
define(`TBSV_TESTGEN_HEAD', 
  `void TBSV_TESTGEN_NAME($1, $2)(int norm, enum blas_order_type order, dnl
       enum blas_uplo_type uplo, enum blas_trans_type trans, dnl
       enum blas_diag_type diag, int n, int k, int randomize, dnl
       $1_array alpha, int alpha_flag, $2_array T, int ldt, dnl
       $1_array x,int *seed, double *HEAD(r_true), double *TAIL(r_true), dnl
       int row, enum blas_prec_type prec)')dnl
dnl
dnl
dnl ---------------------------------------------------------------------
dnl Usage: TBSV_TESTGEN(ax_typeltr, T_typeltr)
dnl        produce tbsv_prepare signature
dnl ---------------------------------------------------------------------
define(`TBSV_TESTGEN', `dnl
TBSV_TESTGEN_HEAD($1, $2)
TBSV_TESTGEN_COMMENT($1, $2)
IF_REAL($1_type, `TBSV_TESTGEN_BODY($1, $2)', 
  `IF_REAL($2_type, 
    `TBSV_TESTGEN_MIX_COMPLEX_BODY($1, $2)', 
    `TBSV_TESTGEN_PURE_COMPLEX_BODY($1, $2)')')
')dnl
dnl
dnl
dnl ---------------------------------------------------------------------
dnl Usage: TBSV_TESTGEN_BODY(ax_typeltr, T_typeltr)
dnl        
dnl ---------------------------------------------------------------------
define(`TBSV_TESTGEN_BODY',
`{
  int start;
  int length;
  int i, j;
  int tempi, inc_tempi;
  DECLARE(alpha_i, real_S)
  DECLARE(minus_one, $1_type)
  DECLARE(Tii, $2_type)
  DECLARE_VECTOR(temp, $2_type)
  DECLARE_VECTOR(xtemp2, $1_type)         

  MALLOC_VECTOR(temp, $2_type, n*2)
  
  /* always allocate, not always needed */
  MALLOC_VECTOR(xtemp2, $1_type, n*2)

  minus_one = -1.0;  

  /* if alpha_flag=0, gives a random value to alpha */
  if (alpha_flag == 0){
    alpha_i = xrand(seed);
    *alpha = alpha_i; 
    alpha_flag = 1;
  }
  
  for(i=0; i<4*n*ldt; i++){
    ZERO(T[i], $1_type)
  }

  inc_tempi = 1;
  INC_ADJUST(inc_tempi, $2_type)

  for(i=0; i<n; i++){
  
    if (i!=row){
      if (diag == blas_non_unit_diag){
        Tii   = xrand(seed);
       }
      else {
        ONE(Tii, $2_type)
      }
      for(j=0, tempi=0; j<n; j++, tempi += inc_tempi) {
        if(j != i) {
          SET_ZERO_VECTOR_ELEMENT(temp, tempi, $2_type)
        } else {
          SET_VECTOR_ELEMENT(temp, tempi, Tii, $2_type) 
        }
      }
      $2tbsv_commit(order, uplo, trans, n, k, T, ldt, temp, i); 


      x[i] = xrand(seed);

      switch(prec){
      case blas_prec_single:     
      {
          DECLARE(multemp, $1_type)
          DECLARE(divtemp, $1_type)

          MUL(multemp, $1_type, x[i], $1_type, *alpha, $1_type)
          DIV(divtemp, $1_type, multemp, $1_type, Tii, $2_type)
          ASSIGN(HEAD(r_true)[i], real_D, divtemp, $1_type)
          TAIL(r_true)[i]=0.0;
          ASSIGN(xtemp2[i], $1_type, divtemp, $1_type)
          break;
      }
      case blas_prec_indigenous:
      case blas_prec_double:
      {
  ifelse(`$1', `s',
          `DECLARE(multemp, real_D)
          DECLARE(divtemp, real_D)

          MUL(multemp, real_D, x[i], $1_type, *alpha, $1_type)
          DIV(divtemp, real_D, multemp, $1_type, Tii, $2_type)
          ASSIGN(HEAD(r_true)[i], real_D, divtemp, real_D)
          TAIL(r_true)[i]=0.0;',
         `$1', `d',
         `DECLARE(multemp, $1_type)
          DECLARE(divtemp, $1_type)

          MUL(multemp, $1_type, x[i], $1_type, *alpha, $1_type)
          DIV(divtemp, $1_type, multemp, $1_type, Tii, $2_type)
          ASSIGN(HEAD(r_true)[i], real_D, divtemp, $1_type)
          TAIL(r_true)[i]=0.0;
          ASSIGN(xtemp2[i], $1_type, divtemp, $1_type)')
          break;
      }
      case blas_prec_extra:
      {
          DECLARE(multemp, real_E)
          DECLARE(divtemp, real_E)

          MUL(multemp, real_E, x[i], $1_type, *alpha, $1_type)
          DIV(divtemp, real_E, multemp, real_E, Tii, $2_type)
          ASSIGN(HEAD(r_true)[i], real_D, HEAD(divtemp), real_D)
          ASSIGN(TAIL(r_true)[i], real_D, TAIL(divtemp), real_D)
          break;
      }
      } /* case */
    } /* if */
  } /* for */

  for(j=0; j<n; j++){
    SET_ZERO_VECTOR_ELEMENT(temp, j, $2_type)
  }

  /* now set the important row - */
  ONE(Tii, $2_type)
  for(j=0, tempi=0; j<n; j++, tempi += inc_tempi) {
        if(j != row) {
                SET_ZERO_VECTOR_ELEMENT(temp, tempi, $2_type)
        } else {
                SET_VECTOR_ELEMENT(temp, tempi, Tii, $2_type)   
        }
  }
        /* this is extra, it will really be committed later */
  $2tbsv_commit(order, uplo, trans, n, k, T, ldt, temp, row);

 
  if ((uplo==blas_lower && (trans==blas_no_trans || trans == blas_conj)) ||
      (uplo==blas_upper && (trans==blas_trans || trans == blas_conj_trans))){
     length = MIN(row, k);
     start  = MAX(0, row - k);
  }
  else{
     length = MIN(n-row-1, k);
     start  = row+1;
  }

  if (length != 0){

    ifelse(`$1', `$2', 
`    switch (prec){
     case blas_prec_single: BLAS_$1dot_testgen(length, 0, length, norm, 
                                blas_no_conj, &minus_one, 1, alpha, 1, 
                                &xtemp2[start], &temp[start], seed, &x[row], 
                                &HEAD(r_true)[row], &TAIL(r_true)[row]);
                            break;
     case blas_prec_indigenous:
     case blas_prec_double:
     case blas_prec_extra:
                        BLAS_$1dot_x_testgen(length,0, length, norm, 
                                blas_no_conj, &minus_one, 1, alpha, 1, 
                                &HEAD(r_true)[start], &TAIL(r_true)[start], 
                                &temp[start], 
                                seed, &x[row], &HEAD(r_true)[row], &TAIL(r_true)[row]);
                            break;
     }',
       `$1&&$2', `d&&s',        
`    switch (prec){
     case blas_prec_single:
     case blas_prec_indigenous:
     case blas_prec_double: /*BLAS_ddot_s_x_testgen(length, 0, length, norm, 
                                blas_no_conj, &minus_one, 1, alpha, 1, 
                                &HEAD(r_true)[start], &TAIL(r_true)[start], 
                                &temp[start], 
                                seed, &x[row], &HEAD(r_true)[row], &TAIL(r_true)[row]);
                            break;*/
     case blas_prec_extra:
                        BLAS_ddot_s_x_testgen(length,0, length, norm, 
                                blas_no_conj, &minus_one, 1, alpha, 1, 
                                &HEAD(r_true)[start], &TAIL(r_true)[start], 
                                &temp[start], 
                                seed, &x[row], &HEAD(r_true)[row], &TAIL(r_true)[row]);
                            break; 
     }') 
    $2tbsv_commit(order, uplo, trans, n, k, T, ldt, temp, row);
  }
  else{
        /* probably low k- must use tricks to test */
        /* tricks are not yet implemented */
    x[row] = xrand(seed);
    
    switch(prec){
    case blas_prec_single:
    { 
      DECLARE(multemp, $1_type)
      
      MUL(multemp, $1_type, x[row], $1_type, *alpha, $1_type)
      HEAD(r_true)[row]=multemp;
      TAIL(r_true)[row]=0.0;
      break;
    }
    case blas_prec_indigenous:
    case blas_prec_double:
    {
ifelse(`$1', `s', 
     `DECLARE(multemp, real_D)
      
      MUL(multemp, real_D, x[row], $1_type, *alpha, $1_type)
      HEAD(r_true)[row]=multemp;
      TAIL(r_true)[row]=0.0;',
      `$1', `d',
      `DECLARE(multemp, $1_type)
      
      MUL(multemp, $1_type, x[row], $1_type, *alpha, $1_type)
      HEAD(r_true)[row]=multemp;
      TAIL(r_true)[row]=0.0;')
      break;
    } 
    case blas_prec_extra:
    {
      DECLARE(multemp, real_E)
      
      MUL(multemp, real_E, x[row], $1_type, *alpha, $1_type)    
      ASSIGN(HEAD(r_true)[row], real_D, HEAD(multemp), real_D)
      ASSIGN(TAIL(r_true)[row], real_D, TAIL(multemp), real_D) 
      break;
    }
    }   
  }      

  FREE_VECTOR(temp, $2_type)
  FREE_VECTOR(xtemp2, $1_type)
}')dnl
dnl
dnl
dnl ---------------------------------------------------------------------
dnl Usage: DOT_TBSV_NAME(abr_typeltr, x_typeltr, y_typeltr)
dnl        produce dot_testgen name
dnl ---------------------------------------------------------------------
define(`DOT_TESTGEN_X_NAME', `ifelse(
        `$1&&$2', `$1&&$1', `BLAS_$1dot_x_testgen',
        `BLAS_$1dot_$2_x_testgen')')dnl
dnl
dnl
define(`TBSV_TESTGEN_PURE_COMPLEX_BODY',
`{  
  PTR_CAST(x, $1_type)
  PTR_CAST(alpha, $1_type)
  PTR_CAST(T, $2_type)
  DECLARE(alpha_r, REAL_TYPE($1_type))
  DECLARE_VECTOR(T_r, REAL_TYPE($2_type))
  DECLARE_VECTOR(x_r, REAL_TYPE($1_type))
  DECLARE_VECTOR(T_temp_r, REAL_TYPE($2_type))
  DECLARE_VECTOR(T_temp_c, $2_type)
  DECLARE_VECTOR(r_true_r, real_E)
  int i, inc=2, length, j;

  MALLOC_VECTOR(T_r, REAL_TYPE($2_type), 8*n*ldt)
  MALLOC_VECTOR(x_r, REAL_TYPE($1_type), n)
  MALLOC_VECTOR(r_true_r, real_E, n)
  MALLOC_VECTOR(T_temp_c, $2_type, n)
  MALLOC_VECTOR(T_temp_r, REAL_TYPE($2_type), n)

  if (alpha_flag == 1){
    alpha_r = alpha_i[0];
  }

  if ((uplo==blas_lower && trans==blas_no_trans) ||
      (uplo==blas_upper && trans!=blas_no_trans)){
     length = row;
  }
  else{
     length = n-row-1;
  }
        
  REAL_TBSV_NAME($1, $2)(norm, order, uplo, trans, diag, n, k, randomize, 
         &alpha_r, 
        alpha_flag, T_r, ldt, x_r, seed, HEAD(r_true_r), TAIL(r_true_r), 
        row, prec);

  /* multiply alpha_r by 1+i */
  alpha_i[0] = alpha_r;
  alpha_i[1] = alpha_r;

  if (diag == blas_non_unit_diag){
    for(i=0; i<n; i++){
      /* multiply x_r by i */
      x_i[i*inc] = 0.0;
      x_i[i*inc+1] = x_r[i];
 
         /* multiply non major rows truth by i*/
      if (i != row){
        HEAD(r_true)[i*inc] = 0.0;
        HEAD(r_true)[i*inc+1] = HEAD(r_true_r)[i];
        TAIL(r_true)[i*inc] = 0.0;
        TAIL(r_true)[i*inc+1] = TAIL(r_true_r)[i];
      }
      else{
        /* multiply major rows truth by (-1+i) */
        HEAD(r_true)[i*inc] = -HEAD(r_true_r)[i];
        HEAD(r_true)[i*inc+1] = HEAD(r_true_r)[i];
        TAIL(r_true)[i*inc] = -TAIL(r_true_r)[i];
        TAIL(r_true)[i*inc+1] = TAIL(r_true_r)[i];
      }
    }

        /* copy T real to complex (and multiply by (1+i)) */
    for(i=0; i < n; i++) {
        REAL_ABBREV($2)tbsv_copy(order, uplo, trans, n, k, T_r, ldt, T_temp_r, i);
        for(j=0; j<n; j++) {
                T_temp_c[j*inc] = T_temp_r[j];
                /* conjugation is handled by the commit later - */
                T_temp_c[j*inc+1] = T_temp_r[j];
        }
        if (i == row) {
                /* zero out imaginary part of diagonal on important row. */
                T_temp_c[row*inc + 1] = 0.0;
        }
        $2tbsv_commit(order, uplo, trans, n, k, T_i, ldt, T_temp_c, i);   
    } 
  }
  else{
    for(i=0; i<n; i++){
        /* multiply x by i */ 
      x_i[i*inc] = 0.0;
      x_i[i*inc+1] = x_r[i];

        /* multiply non-major rows by -1+i */
      if (i != row || length == 0){
        HEAD(r_true)[i*inc] = -HEAD(r_true_r)[i];
        HEAD(r_true)[i*inc+1] = HEAD(r_true_r)[i];
        TAIL(r_true)[i*inc] = -TAIL(r_true_r)[i];
        TAIL(r_true)[i*inc+1] = TAIL(r_true_r)[i];
      }
      else{
        /* multiply x by 1+i (overwrite multiply by i earlier) */
        x_i[i*inc] = x_r[i];    
        x_i[i*inc+1] = x_r[i];  

        HEAD(r_true)[i*inc] = 0.0;
        HEAD(r_true)[i*inc+1] = 2*HEAD(r_true_r)[i];
        TAIL(r_true)[i*inc] = 0.0;
        TAIL(r_true)[i*inc+1] = 2*TAIL(r_true_r)[i];
      }
    }

        /* copy T real to complex (and multiply by (1-i)) */
    for(i=0; i < n; i++) {
        REAL_ABBREV($2)tbsv_copy(order, uplo, trans, n, k, T_r, ldt, T_temp_r, i);
        for(j=0; j<n; j++) {
                T_temp_c[j*inc] = T_temp_r[j];
                /* conjugation is handled by the commit later - */
                T_temp_c[j*inc+1] = -T_temp_r[j];
        }
        /* zero out imaginary part of diagonal on every row. */
        T_temp_c[i*inc + 1] = 0.0;
        $2tbsv_commit(order, uplo, trans, n, k, T_i, ldt, T_temp_c, i);   
    } 
  }

  FREE_VECTOR(T_temp_c, $2_type)
  FREE_VECTOR(T_temp_r, REAL_TYPE($2_type))
  FREE_VECTOR(T_r, REAL_TYPE($2_type))
  FREE_VECTOR(x_r, REAL_TYPE($1_type))
  FREE_VECTOR(r_true_r, real_E)
}')dnl
dnl
dnl
dnl
define(`TBSV_TESTGEN_MIX_COMPLEX_BODY',
`{  
  PTR_CAST(x, $1_type)
  PTR_CAST(alpha, $1_type)
  PTR_CAST(T, $2_type)
  DECLARE(alpha_r, REAL_TYPE($1_type))
  DECLARE_VECTOR(x_r, REAL_TYPE($1_type))
  double *HEAD(r_true_r), *TAIL(r_true_r);
  int i, inc=2;

  MALLOC_VECTOR(x_r, REAL_TYPE($1_type), n)
  MALLOC_VECTOR(r_true_r, real_E, n)

  if (alpha_flag == 1){
    alpha_r = alpha_i[0];
  }
        
  REAL_TBSV_NAME($1, $2)(norm, order, uplo, trans, diag, n, k, 
        randomize, &alpha_r, 
        alpha_flag, T_i, ldt, x_r, seed, HEAD(r_true_r), TAIL(r_true_r), 
        row, prec);

  alpha_i[0] = alpha_r;
  alpha_i[1] = alpha_r;
  
  for(i=0; i<n; i++){
    x_i[i*inc] = 0.0;
    x_i[i*inc+1] = x_r[i];

    HEAD(r_true)[i*inc] = -HEAD(r_true_r)[i];
    HEAD(r_true)[i*inc+1] = HEAD(r_true_r)[i];
    TAIL(r_true)[i*inc] = -TAIL(r_true_r)[i];
    TAIL(r_true)[i*inc+1] = TAIL(r_true_r)[i];
  }

  FREE_VECTOR(x_r, REAL_TYPE($1_type))
  FREE_VECTOR(r_true_r, real_E)
}')dnl
dnl
dnl
define(`REAL_TBSV_NAME', `ifelse(
        `$1&&$2', `c&&c', `BLAS_stbsv_testgen',
        `$1&&$2', `z&&z', `BLAS_dtbsv_testgen',
        `$1&&$2', `z&&c', `BLAS_dtbsv_s_testgen',
        `$1&&$2', `c&&s', `BLAS_stbsv_testgen',
        `$1&&$2', `z&&d', `BLAS_dtbsv_testgen')')dnl
dnl
dnl
define(`PROTOTYPES', `dnl
TBSV_TESTGEN_HEAD(s, s);
TBSV_TESTGEN_HEAD(d, d);
TBSV_TESTGEN_HEAD(d, s);
TBSV_TESTGEN_HEAD(c, c);
TBSV_TESTGEN_HEAD(z, c);
TBSV_TESTGEN_HEAD(z, z);
TBSV_TESTGEN_HEAD(c, s);
TBSV_TESTGEN_HEAD(z, d);
')dnl
dnl
dnl
define(`SOURCE', `dnl
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"

TBSV_TESTGEN(s, s)
TBSV_TESTGEN(d, d)
TBSV_TESTGEN(d, s)
TBSV_TESTGEN(c, c)
TBSV_TESTGEN(z, c)
TBSV_TESTGEN(z, z)
TBSV_TESTGEN(c, s)
TBSV_TESTGEN(z, d)
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
