dnl ***************************************************************************
dnl *                Generate trmv functions                                  *
dnl *                x <-- alpha * T * x                                      *
dnl ***************************************************************************
dnl
#include "blas_extended.h"
#include "blas_extended_private.h"
include(cblas.m4)dnl
include(trmv-common.m4)dnl
dnl -----------------------------------------------------------------------
dnl Usage: TRMV_COMMENT(ax_type, T_type, x_type, extended)
dnl
dnl    ax_type   : type and precision of alpha and x.
dnl    T_type    : type and precision of T
dnl    extended  : `_x' if extended precision required.
dnl 
dnl -----------------------------------------------------------------------
define(`TRMV_COMMENT', `
/*
 * Purpose
 * =======
 *
 * Computes x <-- alpha * T * x, where T is a triangular matrix.
 *
 * Arguments
 * =========
 * 
 * order  (input) enum blas_order_type
 *        column major, row major
 *
 * uplo   (input) enum blas_uplo_type
 *        upper, lower
 *
 * trans  (input) enum blas_trans_type
 *        no trans, trans, conj trans
 * 
 * diag   (input) enum blas_diag_type
 *        unit, non unit
 *
 * n      (input) int
 *        the dimension of T
 * 
 * alpha  (input) $1_scalar
 * 
 * T      (input) $2_array
 *        Triangular matrix
 *
 * ldt    (input) int 
 *        Leading dimension of T
 *
 * x      (input) const $1_array
 *    Array of length n.
 * 
 * incx   (input) int
 *     The stride used to access components x[i].
 *
PREC_COMMENT($3)dnl
 */')dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: TRMV_LOOP(ax_type, T_type, sum_type, tmp_type, FPU_flag)
dnl
dnl                             n-1
dnl            x[i] <-- alpha * SUM  (T[i][j] * x[j]) 
dnl                             j=0
dnl
dnl   $1   ax_type  : type and precision of alpha and x
dnl   $2   T_type   : type and precision of T
dnl   $3   sum_type : the type and precision of auxiliary variables prod/sum
dnl   $4   tmp_type : the type and precision of auxiliary variable tmp
dnl   $5   FPU      : [optional] The string `FPU' is passed if FPU
dnl                   fix is needed.  Empty string is passed otherwise.
dnl
dnl Each type and precision specifier can be one of
dnl     real_S   ... real and single
dnl     real_D   ... real and double
dnl     real_I   ... real and indigenous
dnl     real_E   ... real and extra
dnl     complex_S  ... complex and single
dnl     complex_D  ... complex and double
dnl     complex_I  ... complex and indigeneous
dnl     complex_E  ... complex and extra
dnl ----------------------------------------------------------------------
define(`TRMV_BODY', `
  int i, j;    /* used to idx matrix */
  int xj, xj0;
  int ti, tij, tij0;

  int inc_ti, inc_tij;
  int inc_x;

  PTR_CAST(T, $2, `const') /* internal matrix T */
  PTR_CAST(x, $1, `')      /* internal x */
  SCALAR_CAST(alpha, $1)   /* internal alpha */

  DECLARE(t_elem, $2)
  DECLARE(x_elem, $1)
  DECLARE(prod, $3)
  DECLARE(sum, $3)
  DECLARE(tmp, $4)

  ifelse(`$5', `FPU', `FPU_FIX_DECL; FPU_FIX_START;')
   
  /* all error calls */
  if ((order != blas_rowmajor && order != blas_colmajor) ||
      (uplo != blas_upper && uplo != blas_lower) ||
      (trans != blas_trans && 
       trans != blas_no_trans && 
       trans != blas_conj_trans) ||
      (diag != blas_non_unit_diag && diag != blas_unit_diag) ||
      (ldt < n) ||
      (incx == 0)) {
    BLAS_error(routine_name, 0, 0, NULL);
  } else if (n <= 0) {
    BLAS_error(routine_name, -4, n, NULL);
  } else if (incx == 0) {
    BLAS_error(routine_name, -9, incx, NULL);
  }

  if (trans == blas_no_trans) {
    if (uplo == blas_upper) {
      inc_x = -incx;
      if (order == blas_rowmajor) {
        inc_ti = ldt;
        inc_tij = -1;
      } else {
        inc_ti = 1;
        inc_tij = -ldt;
      }
    } else {
      inc_x = incx;
      if (order == blas_rowmajor) {
        inc_ti = -ldt;
        inc_tij = 1;
      } else {
        inc_ti = -1;
        inc_tij = ldt;
      }
    }
  } else {
    if (uplo == blas_upper) {
      inc_x = incx;
      if (order == blas_rowmajor) {
        inc_ti = -1;
        inc_tij = ldt;
      } else {
        inc_ti = -ldt;
        inc_tij = 1;
      }
    } else {
      inc_x = -incx;
      if (order == blas_rowmajor) {
        inc_ti = 1;
        inc_tij = -ldt;
      } else {
        inc_ti = ldt;
        inc_tij = -1;
      }
    }
  }

  INC_ADJUST(inc_ti, $2)
  INC_ADJUST(inc_tij, $2)
  INC_ADJUST(inc_x, $1)

  xj0 = (inc_x > 0 ? 0 : -(n-1) * inc_x);
  if (TEST_0(alpha_i, $1)) {
    xj = xj0;
    for (j = 0; j < n; j++) {
      SET_ZERO_VECTOR_ELEMENT(x_i, xj, $1)
      xj += inc_x;
    }
  } else {

    if (diag == blas_unit_diag) {
      TRMV_CONJ($1, $2, $3, $4, blas_unit_diag)
    } else {
      TRMV_CONJ($1, $2, $3, $4, blas_non_unit_diag)
    }

  }

  ifelse(`$5', `FPU', `FPU_FIX_STOP;')
')dnl
dnl
dnl
define(`TRMV_CONJ', 
 `IF_COMPLEX($2, 
   `if (trans == blas_conj_trans) {
      TRMV_LOOP($1, $2, $3, $4, blas_conj, $5)
    } else {
      TRMV_LOOP($1, $2, $3, $4, blas_no_conj, $5)
    }', 
    `TRMV_LOOP($1, $2, $3, $4, blas_no_conj, $5)')')dnl
dnl
dnl
dnl ---------------------------------------------------------------------
dnl Usage: TRMV_LOOP(ax_type, T_type, sum_type, tmp_type, conj, diag)
dnl                  x = alpha * T * x
dnl ---------------------------------------------------------------------
define(`TRMV_LOOP',`
  
  ti = (inc_ti > 0 ? 0 : -(n-1) * inc_ti);
  tij0 = (inc_tij > 0 ? 0 : -(n-1) * inc_tij);
  for (i = 0; i < n; i++) {

    ZERO(sum, $3)

    xj = xj0;
    tij = ti + tij0;
    ifelse(`$6', `blas_unit_diag', 
      `for (j = i; j < (n-1); j++) {',
      `for (j = i; j < n; j++) {')

      GET_VECTOR_ELEMENT(t_elem, T_i, tij, $2)
      CONJ(t_elem, $2, $5)
      GET_VECTOR_ELEMENT(x_elem, x_i, xj, $1)
      MUL(prod, $3, x_elem, $1, t_elem, $2)
      ADD(sum, $3, sum, $3, prod, $3)

      xj += inc_x;
      tij += inc_tij;
    }

    ifelse(`$6', `blas_unit_diag', 
     `GET_VECTOR_ELEMENT(x_elem, x_i, xj, $1)
      ADD(sum, $3, sum, $3, x_elem, $1)

      if (TEST_1(alpha_i, $1)) {
        SET_ROUND_VECTOR_ELEMENT(x_i, xj, sum, $3)
      } else {
        MUL(tmp, $4, sum, $3, alpha_i, $1)
        SET_ROUND_VECTOR_ELEMENT(x_i, xj, tmp, $3)
      }', 
      `if (TEST_1(alpha_i, $1)) {
        SET_ROUND_VECTOR_ELEMENT(x_i, xj - inc_x, sum, $3)
      } else {
        MUL(tmp, $4, sum, $3, alpha_i, $1)
        SET_ROUND_VECTOR_ELEMENT(x_i, xj - inc_x, tmp, $3)
      }')

    ti += inc_ti;
  }
')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: SWITCH_prec($1, $2, $3, $4, $5, $6, $7, $8)
dnl        
dnl        Generates a 3-way switch statement based on prec.
dnl        $1 is the type of alpha and x 
dnl        $2 is the type of T 
dnl        $3, $4 are types of sum and tmp in single case
dnl        $5, $6 are types of sum and tmp in double / indigenous case
dnl        $7, $8 are types of sum and tmp in extra case
dnl ----------------------------------------------------------------------
define(`SWITCH_prec', `
  switch (prec) {
    case blas_prec_single: ifelse(`$3&&$4', `$5&&$6', `', `{
      TRMV_BODY($1, $2, $3, $4)
      break;
    }
    ')dnl
    case blas_prec_double:
    case blas_prec_indigenous: {
      TRMV_BODY($1, $2, $5, $6)
      break;
    }

    case blas_prec_extra: { 
      TRMV_BODY($1, $2, $7, $8, FPU)
      break;
    }
  }')dnl
dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: TRMV_X_BODY(ax_type, T_type) ... dispatches
dnl        TRMV with appropriate type and precision info of
dnl        the specified internal prec.       
dnl Each type specifier can be one of
dnl     s   ... real and single
dnl     d   ... real and double
dnl     c   ... complex and single
dnl     z   ... complex and double
dnl ----------------------------------------------------------------------
define(`TRMV_X_BODY', 
 `SWITCH_prec($1, $2, 
   TMP_TYPE_X($1, S), TMP_TYPE_X($1, S), 
   TMP_TYPE_X($1, D), TMP_TYPE_X($1, D), 
   TMP_TYPE_X($1, E), TMP_TYPE_X($1, E))')dnl
dnl
dnl
dnl --------------------------------------------------------------------
dnl Usage: TRMV(ax_type, T_type, extended) 
dnl        ... generates the calling sequence for the standard case.
dnl --------------------------------------------------------------------
define(`TRMV', `
TRMV_HEAD($1, $2, $3)
TRMV_COMMENT($1, $2, $3)
{
  static const char routine_name[] = "TRMV_NAME($1, $2, $3)";
  ifelse(`$3', `_x', `TRMV_X_BODY($1_type, $2_type)', 
    `TRMV_BODY($1_type, $2_type, HIGHER_TYPE($1_type, $2_type), $1_type)')
}')dnl
dnl
dnl
