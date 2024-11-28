dnl ---------------------------------------------------------
dnl  HEMM ---- Hermitian Matrix Multiply
dnl
dnl    C  <---   alpha * A * B + beta * C
dnl    C  <---   alpha * B * A + beta * C
dnl
dnl    where matrix A is Hermitian.
dnl ---------------------------------------------------------
dnl
#include "blas_extended.h"
#include "blas_extended_private.h"
include(cblas.m4)dnl
include(hemm-common.m4)dnl
dnl
dnl
dnl  Usage:
dnl    HEMM        ($1, $2, $3, $4)
dnl    HEMM_HEAD   ($1, $2, $3, $4)
dnl    HEMM_NAME   ($1, $2, $3, $4)
dnl    HEMM_PARAMS ($1, $2, $3, $4)
dnl    HEMM_COMMENT($1, $2, $3, $4)
dnl
dnl    $1 -- type of alpha, beta, c.
dnl    $2 -- type of b
dnl    $3 -- type of c
dnl    $4 -- set to `_x' for _x routines.
dnl          Otherwise, `'.
dnl
dnl
define(`HEMM_COMMENT', `
/* 
 * Purpose
 * =======
 *
 * This routines computes one of the matrix product:
 *
 *     C  <-  alpha * A * B  +  beta * C
 *     C  <-  alpha * B * A  +  beta * C
 * 
 * where A is a hermitian matrix.
 *
 * Arguments
 * =========
 *
 * order   (input) enum blas_order_type
 *         Storage format of input matrices A, B, and C.
 * 
 * side    (input) enum blas_side_type
 *         Determines which side of matrix B is matrix A is multiplied.
 *
 * uplo    (input) enum blas_uplo_type
 *         Determines which half of matrix A (upper or lower triangle)
 *         is accessed.
 *
 * m n     (input) int
 *         Size of matrices A, B, and C.
 *         Matrix A is m-by-m if it is multiplied on the left, 
 *                     n-by-n otherwise.
 *         Matrices B and C are m-by-n.
 *
 * alpha   (input) $1_scalar
 * 
 * a       (input) const $2_array
 *         Matrix A.
 *
 * lda     (input) int
 *         Leading dimension of matrix A.
 *
 * b       (input) const $3_array
 *         Matrix B.
 *   
 * ldb     (input) int
 *         Leading dimension of matrix B.
 *
 * beta    (input) $1_scalar
 * 
 * c       (input/output) $1_array
 *         Matrix C.
 * 
 * ldc     (input) int
 *         Leading dimension of matrix C.
 *
PREC_COMMENT($4)dnl
 */')dnl
dnl
dnl
dnl Usage:  HEMM_BODY($1, $2, $3, $4, $5, [$6])
dnl     $1 - type of matrix a
dnl     $2 - type of matrix b
dnl     $3 - type of alpha, beta, and matrix c
dnl     $4 - type of sum/prod
dnl     $5 - type of temp
dnl     $6 - [optional] String `FPU' is passed if FPU fix
dnl          is needed.  The empty string is passed otherwise.
dnl
define(`HEMM_BODY', `
  /* Integer Index Variables */
  int i, j, k;

  int ai, bj, ci;
  int aik, bkj, cij;

  int incai, incbj, incci;
  int incaik1, incaik2, incbkj, inccij;

  int m_i, n_i;

  int conj_flag;

  /* Input Matrices */
  PTR_CAST(a, $2, `const')
  PTR_CAST(b, $3, `const')

  /* Output Matrix */
  PTR_CAST(c, $1)

  /* Input Scalars */
  SCALAR_CAST(alpha, $1)
  SCALAR_CAST(beta, $1)

  /* Temporary Floating-Point Variables */
  DECLARE(a_elem, $2)
  DECLARE(b_elem, $3)
  DECLARE(c_elem, $1)
  DECLARE(prod, $4)
  DECLARE(sum, $4)
  DECLARE(tmp1, $5)
  DECLARE(tmp2, $5)

  ifelse(`$6', `FPU', `FPU_FIX_DECL;')

  /* Check for error conditions. */
  if (m <= 0 || n <= 0) {
    return;
  }

  if (order == blas_colmajor && (ldb < m || ldc < m)) {
    return;
  }
  if (order == blas_rowmajor && (ldb < n || ldc < n)) {
    return;
  }

  if (side == blas_left_side && lda < m) {
    return;
  }

  if (side == blas_right_side && lda < n) {
    return;
  }

  /* Test for no-op */
  if (TEST_0(alpha_i, $1) && TEST_1(beta_i, $1)) {
    return;
  }

  /* Set Index Parameters */
  if (side == blas_left_side) {
    m_i = m; n_i = n;
  } else {
    m_i = n; n_i = m;
  }

  if ((order == blas_colmajor && side == blas_left_side) ||
      (order == blas_rowmajor && side == blas_right_side)) {
    incbj = ldb;
    incbkj = 1;
    incci = 1;
    inccij = ldc;
  } else {
    incbj = 1;
    incbkj = ldb;
    incci = ldc;
    inccij = 1;
  }

  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai = lda;
    incaik1 = 1;
    incaik2 = lda;
  } else {
    incai = 1;
    incaik1 = lda;
    incaik2 = 1;
  }

  if ((side == blas_left_side  && uplo == blas_upper) ||
      (side == blas_right_side && uplo == blas_lower))
    conj_flag = 1;
  else
    conj_flag = 0;

  ifelse(`$6', `FPU', `FPU_FIX_START;')

  /* Adjustment to increments (if any) */
  INC_ADJUST(incci, $1)
  INC_ADJUST(inccij, $1)
  INC_ADJUST(incai, $2)
  INC_ADJUST(incaik1, $2)
  INC_ADJUST(incaik2, $2)
  INC_ADJUST(incbj, $3)
  INC_ADJUST(incbkj, $3)
 
  /* alpha = 0.  In this case, just return beta * C */
  if (TEST_0(alpha_i, $1)) {
    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {
        GET_VECTOR_ELEMENT(c_elem, c_i, cij, $1)
        MUL(tmp1, $5, c_elem, $1, beta_i, $1)
        SET_ROUND_VECTOR_ELEMENT(c_i, cij, tmp1, $5)
      }
    }
  } else if (TEST_1(alpha_i, $1)) {
    
    /* Case alpha == 1. */
    
    if (TEST_0(beta_i, $1)) {
      /* Case alpha = 1, beta = 0.  We compute  C <--- A * B   or  B * A */
      for (i = 0, ci = 0, ai = 0; i < m_i; i++, ci += incci, ai += incai) {
        for (j = 0, cij = ci, bj = 0; j < n_i; 
             j++, cij += inccij, bj += incbj) {
          
          ZERO(sum, $4)
          
          for (k = 0, aik = ai, bkj = bj; k < i; 
               k++, aik += incaik1, bkj += incbkj) {
            GET_VECTOR_ELEMENT(a_elem, a_i, aik, $2)
            GET_VECTOR_ELEMENT(b_elem, b_i, bkj, $3)
            if (conj_flag == 1) {
              CONJ_AUX(a_elem, $2)
            }
            MUL(prod, $4, a_elem, $2, b_elem, $3)
            ADD(sum, $4, sum, $4, prod, $4)
          }
          for (; k < m_i; k++, aik += incaik2, bkj += incbkj) {
            GET_VECTOR_ELEMENT(a_elem, a_i, aik, $2)
            GET_VECTOR_ELEMENT(b_elem, b_i, bkj, $3)
            if (conj_flag == 0) {
              CONJ_AUX(a_elem, $2)
            }
            MUL(prod, $4, a_elem, $2, b_elem, $3)
            ADD(sum, $4, sum, $4, prod, $4)
          }
          SET_ROUND_VECTOR_ELEMENT(c_i, cij, sum, $4)
        }
      }
    } else {
      /* Case alpha = 1, but beta != 0. 
         We compute  C  <--- A * B + beta * C 
                 or  C  <--- B * A + beta * C  */
      
      for (i = 0, ci = 0, ai = 0; i < m_i; i++, ci += incci, ai += incai) {
        for (j = 0, cij = ci, bj = 0; j < n_i; 
             j++, cij += inccij, bj += incbj) {
          
          ZERO(sum, $4)
          
          for (k = 0, aik = ai, bkj = bj; k < i; 
               k++, aik += incaik1, bkj += incbkj) {
            GET_VECTOR_ELEMENT(a_elem, a_i, aik, $2)
            GET_VECTOR_ELEMENT(b_elem, b_i, bkj, $3)
            if (conj_flag == 1) {
              CONJ_AUX(a_elem, $2)
            }
            MUL(prod, $4, a_elem, $2, b_elem, $3)
            ADD(sum, $4, sum, $4, prod, $4)
          }
          for (; k < m_i; k++, aik += incaik2, bkj += incbkj) {
            GET_VECTOR_ELEMENT(a_elem, a_i, aik, $2)
            GET_VECTOR_ELEMENT(b_elem, b_i, bkj, $3)
            if (conj_flag == 0) {
              CONJ_AUX(a_elem, $2)
            }
            MUL(prod, $4, a_elem, $2, b_elem, $3)
            ADD(sum, $4, sum, $4, prod, $4)
          }
          GET_VECTOR_ELEMENT(c_elem, c_i, cij, $1)
          MUL(tmp2, $5, c_elem, $1, beta_i, $1)
          ASSIGN(tmp1, $5, sum, $4)
          ADD(tmp1, $5, tmp2, $5, tmp1, $5)
          SET_ROUND_VECTOR_ELEMENT(c_i, cij, tmp1, $5)
        }
      }
    }
    
  } else {
    /* The most general form,   C <--- alpha * A * B + beta * C  
                           or   C <--- alpha * B * A + beta * C  */
    
    for (i = 0, ci = 0, ai = 0; i < m_i; i++, ci += incci, ai += incai) {
      for (j = 0, cij = ci, bj = 0; j < n_i; 
           j++, cij += inccij, bj += incbj) {
        
        ZERO(sum, $4)
          
        for (k = 0, aik = ai, bkj = bj; k < i; 
             k++, aik += incaik1, bkj += incbkj) {
          GET_VECTOR_ELEMENT(a_elem, a_i, aik, $2)
          GET_VECTOR_ELEMENT(b_elem, b_i, bkj, $3)
          if (conj_flag == 1) {
            CONJ_AUX(a_elem, $2)
          }
          MUL(prod, $4, a_elem, $2, b_elem, $3)
          ADD(sum, $4, sum, $4, prod, $4)
        }
        for (; k < m_i; k++, aik += incaik2, bkj += incbkj) {
          GET_VECTOR_ELEMENT(a_elem, a_i, aik, $2)
          GET_VECTOR_ELEMENT(b_elem, b_i, bkj, $3)
          if (conj_flag == 0) {
            CONJ_AUX(a_elem, $2)
          }
          MUL(prod, $4, a_elem, $2, b_elem, $3)
          ADD(sum, $4, sum, $4, prod, $4)
        }
        MUL(tmp1, $5, sum, $4, alpha_i, $1)
        GET_VECTOR_ELEMENT(c_elem, c_i, cij, $1)
        MUL(tmp2, $5, c_elem, $1, beta_i, $1)
        ADD(tmp1, $5, tmp1, $5, tmp2, $5)
        SET_ROUND_VECTOR_ELEMENT(c_i, cij, tmp1, $5)
      }
    }
  }

  ifelse(`$6', `FPU', `FPU_FIX_STOP;')
')dnl
dnl
dnl
dnl Usage: SWITCH_prec($1, $2, $3, $4, $5, $6, $7, $8, $9)
dnl Generates a 3-way switch statement based on prec.
dnl Arguments
dnl    $1      --  type of matrix C
dnl    $2      --  type of matrix A
dnl    $3      --  type of matrix B
dnl    $4, $5  --  type of `sum' and `tmp' in single case
dnl    $6, $7  --  type of `sum' and `tmp' in double/indigenous case
dnl    $8, $9  --  type of `sum' and `tmp' in extra case
dnl
define(`SWITCH_prec', `switch (prec) {

    case blas_prec_single: ifelse(`$4&&$5', `$6&&$7', `', `{
      HEMM_BODY($1, $2, $3, $4, $5)
      break;
    }
    ')dnl
    case blas_prec_double:
    case blas_prec_indigenous: {
      HEMM_BODY($1, $2, $3, $6, $7)
      break;
    }

    case blas_prec_extra: { 
      HEMM_BODY($1, $2, $3, $8, $9, FPU)
      break;
    }
  }')dnl
dnl
dnl
dnl
dnl
dnl
define(`HEMM_X_BODY', 
  `SWITCH_prec($1, $2, $3, 
    SUM_TYPE_X($1, $2, $3, S), TMP_TYPE_X($1, S),
    SUM_TYPE_X($1, $2, $3, D), TMP_TYPE_X($1, D),
    SUM_TYPE_X($1, $2, $3, E), TMP_TYPE_X($1, E))')dnl
dnl
dnl
define(`HEMM', 
  `HEMM_HEAD($1, $2, $3, $4) 
   HEMM_COMMENT($1, $2, $3, $4)
  {
    ifelse($4, `_x', `HEMM_X_BODY($1_type, $2_type, $3_type)', 
      `HEMM_BODY($1_type, $2_type, $3_type, 
        SUM_TYPE($1_type, $2_type, $3_type), $1_type)')
  }')dnl
dnl
dnl
