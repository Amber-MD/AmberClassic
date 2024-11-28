dnl ------------------------------------------------------
dnl  GEMM --- General Matrix Multiply
dnl     C  <---  alpha * op(A) * op(B) + beta * C
dnl
dnl     where op can be no-op, transpose, conj. transpose.
dnl ------------------------------------------------------
dnl
dnl
#include "blas_extended.h"
#include "blas_extended_private.h"
include(cblas.m4)dnl
include(gemm-common.m4)dnl
dnl
dnl
dnl Usage:  GEMM(c_type, a_type, b_type, sum_type, tmp_type, FPU)
dnl
dnl  $1  c_type   -- type of alpha, beta, and matrix C
dnl  $2  a_type   -- type of matrix A
dnl  $3  b_type   -- type of matrix B
dnl  $4  sum_type -- type of sum/product above
dnl  $5  tmp_type -- type of temporary variable
dnl  $6  FPU      -- [optional] String `FPU' is passed if FPU fix
dnl                  is needed.  The empty string is passed otherwise.
dnl
dnl Each type can be one of
dnl
dnl    real_S    -- real and single
dnl    real_D    -- real and double
dnl    real_I    -- real and indigenous
dnl    real_E    -- real and extra
dnl    complex_S -- complex and single
dnl    complex_D -- complex and double
dnl    complex_I -- complex and indigenous
dnl    complex_E -- complex and extra
dnl
dnl
define(`GEMM_COMMENT', `
/* 
 * Purpose
 * =======
 *
 * This routine computes the matrix product:
 *
 *      C   <-  alpha * op(A) * op(B)  +  beta * C .
 * 
 * where op(M) represents either M, M transpose, 
 * or M conjugate transpose.
 *
 * Arguments
 * =========
 *
 * order   (input) enum blas_order_type
 *         Storage format of input matrices A, B, and C.
 *
 * transa  (input) enum blas_trans_type
 *         Operation to be done on matrix A before multiplication.
 *         Can be no operation, transposition, or conjugate transposition.
 *
 * transb  (input) enum blas_trans_type
 *         Operation to be done on matrix B before multiplication.
 *         Can be no operation, transposition, or conjugate transposition.
 * 
 * m n k   (input) int
 *         The dimensions of matrices A, B, and C.
 *         Matrix C is m-by-n matrix.
 *         Matrix A is m-by-k if A is not transposed, 
 *                     k-by-m otherwise.
 *         Matrix B is k-by-n if B is not transposed, 
 *                     n-by-k otherwise.
 *      
 * alpha   (input) $1_scalar
 *
 * a       (input) const $2_array
 *         matrix A.
 * 
 * lda     (input) int
 *         leading dimension of A.
 * 
 * b       (input) const $3_array
 *         matrix B
 *
 * ldb     (input) int
 *         leading dimension of B.
 *
 * beta    (input) $1_scalar
 *
 * c       (input/output) $1_array
 *         matrix C
 *
 * ldc     (input) int
 *         leading dimension of C.
 *
PREC_COMMENT($4)dnl
 */')dnl
dnl
dnl
define(`GEMM_BODY', `

  /* Integer Index Variables */
  int i, j, h;

  int ai, bj, ci;
  int aih, bhj, cij;  /* Index into matrices a, b, c during multiply */

  int incai, incaih;  /* Index increments for matrix a */
  int incbj, incbhj;  /* Index increments for matrix b */
  int incci, inccij;  /* Index increments for matrix c */

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

  /* Test for error conditions */
  if (m < 0)
    BLAS_error(routine_name, -4, m, NULL);
  if (n < 0)
    BLAS_error(routine_name, -5, n, NULL);
  if (k < 0)
    BLAS_error(routine_name, -6, k, NULL);

  if (order == blas_colmajor) {

    if (ldc < m)
      BLAS_error(routine_name, -14, ldc, NULL);

    if (transa == blas_no_trans) {
      if (lda < m)
        BLAS_error(routine_name, -9, lda, NULL);
    } else {
      if (lda < k)
        BLAS_error(routine_name, -9, lda, NULL);
    }
 
    if (transb == blas_no_trans) {
      if (ldb < k)
        BLAS_error(routine_name, -11, ldb, NULL);
    } else {
      if (ldb < n)
        BLAS_error(routine_name, -11, ldb, NULL);
    }

  } else {
    /* row major */
    if (ldc < n)
      BLAS_error(routine_name, -14, ldc, NULL);

    if (transa == blas_no_trans) {
      if (lda < k)
        BLAS_error(routine_name, -9, lda, NULL);
    } else {
      if (lda < m)
        BLAS_error(routine_name, -9, lda, NULL);
    }
 
    if (transb == blas_no_trans) {
      if (ldb < n)
        BLAS_error(routine_name, -11, ldb, NULL);
    } else {
      if (ldb < k)
        BLAS_error(routine_name, -11, ldb, NULL);
    }
  }
  
  /* Test for no-op */
  if (n == 0 || m == 0 || k == 0)
    return;
  if (TEST_0(alpha_i, $1) && TEST_1(beta_i, $1)) {
    return;
  }

  /* Set Index Parameters */
  if (order == blas_colmajor) {
    incci = 1;
    inccij = ldc;

    if (transa == blas_no_trans) {
      incai = 1;
      incaih = lda; 
    } else {
      incai = lda;
      incaih = 1;
    }

    if (transb == blas_no_trans) {
      incbj = ldb;
      incbhj = 1;
    } else {
      incbj = 1;
      incbhj = ldb;
    }

  } else {
    /* row major */
    incci = ldc;
    inccij = 1;

    if (transa == blas_no_trans) {
      incai = lda;
      incaih = 1; 
    } else {
      incai = 1;
      incaih = lda;
    }

    if (transb == blas_no_trans) {
      incbj = 1;
      incbhj = ldb;
    } else {
      incbj = ldb;
      incbhj = 1;
    }

  }

  ifelse(`$6', `FPU', `FPU_FIX_START;')

  /* Ajustment to increments */
  INC_ADJUST(incci, $1)
  INC_ADJUST(inccij, $1)
  INC_ADJUST(incai, $2)
  INC_ADJUST(incaih, $2)
  INC_ADJUST(incbj, $3)
  INC_ADJUST(incbhj, $3) 

  /* alpha = 0.  In this case, just return beta * C */
  if (TEST_0(alpha_i, $1)) {

    ci = 0;
    for (i = 0; i < m; i++, ci += incci) {
      cij = ci;
      for (j = 0; j < n; j++, cij += inccij) {
        GET_VECTOR_ELEMENT(c_elem, c_i, cij, $1)
        MUL(tmp1, $5, c_elem, $1, beta_i, $1)
        SET_ROUND_VECTOR_ELEMENT(c_i, cij, tmp1, $5)
      }
    }

  } else if (TEST_1(alpha_i, $1)) {

    /* Case alpha == 1. */

    if (TEST_0(beta_i, $1)) {
      /* Case alpha == 1, beta == 0.   We compute  C <--- A * B */

      ci = 0;
      ai = 0;
      for (i = 0; i < m; i++, ci += incci, ai += incai) {

        cij = ci;
        bj = 0;

        for (j = 0; j < n; j++, cij += inccij, bj += incbj) {

          aih = ai;
          bhj = bj;

          ZERO(sum, $4)

          for (h = 0; h < k; h++, aih += incaih, bhj += incbhj) {
            GET_VECTOR_ELEMENT(a_elem, a_i, aih, $2)
            GET_VECTOR_ELEMENT(b_elem, b_i, bhj, $3)
            if (transa == blas_conj_trans) {
              CONJ_AUX(a_elem, $2)
            }
            if (transb == blas_conj_trans) {
              CONJ_AUX(b_elem, $3)
            }
            MUL(prod, $4, a_elem, $2, b_elem, $3)
            ADD(sum, $4, sum, $4, prod, $4)
          }
          ifelse(`$5', `$4', 
                  `SET_ROUND_VECTOR_ELEMENT(c_i, cij, sum, $4)',
                `ASSIGN(tmp1, $5, sum, $4)
                SET_ROUND_VECTOR_ELEMENT(c_i, cij, tmp1, $5)')
        } 
      }

    } else {
      /* Case alpha == 1, but beta != 0.
         We compute   C <--- A * B + beta * C   */

      ci = 0;
      ai = 0;
      for (i = 0; i < m; i++, ci += incci, ai += incai) {

        cij = ci;
        bj = 0;
  
        for (j = 0; j < n; j++, cij += inccij, bj += incbj) {
  
          aih = ai;
          bhj = bj;
  
          ZERO(sum, $4)
  
          for (h = 0; h < k; h++, aih += incaih, bhj += incbhj) {
            GET_VECTOR_ELEMENT(a_elem, a_i, aih, $2)
            GET_VECTOR_ELEMENT(b_elem, b_i, bhj, $3)
            if (transa == blas_conj_trans) {
              CONJ_AUX(a_elem, $2)
            }
            if (transb == blas_conj_trans) {
              CONJ_AUX(b_elem, $3)
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

    /* The most general form,   C <-- alpha * A * B + beta * C  */
    ci = 0;
    ai = 0;
    for (i = 0; i < m; i++, ci += incci, ai += incai) {

      cij = ci;
      bj = 0;

      for (j = 0; j < n; j++, cij += inccij, bj += incbj) {

        aih = ai;
        bhj = bj;

        ZERO(sum, $4)

        for (h = 0; h < k; h++, aih += incaih, bhj += incbhj) {
          GET_VECTOR_ELEMENT(a_elem, a_i, aih, $2)
          GET_VECTOR_ELEMENT(b_elem, b_i, bhj, $3)
          if (transa == blas_conj_trans) {
            CONJ_AUX(a_elem, $2)
          }
          if (transb == blas_conj_trans) {
            CONJ_AUX(b_elem, $3)
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
      GEMM_BODY($1, $2, $3, $4, $5)
      break;
    }
    ')dnl
    case blas_prec_double:
    case blas_prec_indigenous: {
      GEMM_BODY($1, $2, $3, $6, $7)
      break;
    }

    case blas_prec_extra: { 
      GEMM_BODY($1, $2, $3, $8, $9, FPU) 
      break;
    }
  }')dnl
dnl
dnl
define(`GEMM_X_BODY', 
  `SWITCH_prec($1, $2, $3, 
    SUM_TYPE_X($1, $2, $3, S), TMP_TYPE_X($1, S), 
    SUM_TYPE_X($1, $2, $3, D), TMP_TYPE_X($1, D), 
    SUM_TYPE_X($1, $2, $3, E), TMP_TYPE_X($1, E))')dnl
dnl
dnl
define(`GEMM', 
 `GEMM_HEAD($1, $2, $3, $4)
  GEMM_COMMENT($1, $2, $3, $4)
{
  static const char routine_name[] = "GEMM_NAME($1, $2, $3, $4)";
  ifelse($4, `_x', `GEMM_X_BODY($1_type, $2_type, $3_type)', 
    `GEMM_BODY($1_type, $2_type, $3_type, 
      SUM_TYPE($1_type, $2_type, $3_type), $1_type)')
}')dnl
dnl
dnl
