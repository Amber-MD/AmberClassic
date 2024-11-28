dnl ----------------------------------------------------------------------
dnl TPMV -- Triangular packed matrix vector multiply
dnl   x <--- alpha * op(T) * x
dnl
dnl   where op can be no-op, tranpose, or conjugate transpose
dnl ----------------------------------------------------------------------
dnl
dnl tpmv.m4 contains the macros for producing the TPMV (triangular-
dnl packed matrix-vector multiply) routines.  See the documentation
dnl for packed.m4 for data on packed storage.
dnl
dnl
#include "blas_extended.h"
#include "blas_extended_private.h"
include(cblas.m4)dnl
include(packed.m4)dnl
include(tpmv-common.m4)dnl
dnl
dnl
dnl ---------------------------------------------------------------------
dnl Usage: TPMV_COMMENT(ax_typeltr, tp_typeltr, extended) ...
dnl        generates the leading comments for TPMV.
dnl ---------------------------------------------------------------------
define(`TPMV_COMMENT', `
/*
 * Purpose
 * =======
 *
 * Computes x = alpha * tp * x, x = alpha * tp_transpose * x,
 * or x = alpha * tp_conjugate_transpose where tp is a triangular
 * packed matrix.
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of tp; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether tp is upper or lower
 *
 * trans        (input) blas_trans_type
 *
 * diag         (input) blas_diag_type
 *              Whether the diagonal entries of tp are 1
 *
 * n            (input) int
 *              Dimension of tp and the length of vector x
 *
 * alpha        (input) $1_scalar
 *
 * tp           (input) $2_array
 *
 * x            (input) $1_array
 *
 * incx         (input) int
 *              The stride for vector x.
 *
PREC_COMMENT($3)dnl
 */')dnl
dnl
dnl ---------------------------------------------------------------------
dnl Usage: TPMV_LOOP(ax_type, tp_type, sum_type, tmp_type)
dnl ---------------------------------------------------------------------
define(`TPMV_LOOP', `
{
  if ((uplo == blas_upper &&
       trans == blas_no_trans && order == blas_rowmajor) ||
      (uplo == blas_lower &&
       trans != blas_no_trans && order == blas_colmajor)) {
dnl ---------------------------------------------------------------------
dnl 1) upper triangular and row major with no trans
dnl 2) lower triangular and column major with trans
dnl ---------------------------------------------------------------------
    tp_start = 0;
    tp_index = tp_start;
    for (matrix_row = 0; matrix_row < n; matrix_row++) {
      x_index = x_start + incx*matrix_row;
      x_index2 = x_index;
      col_index = matrix_row;
      ZERO(rowsum, $3)
      ZERO(rowtmp, $4)
      ZERO(result, $3)
      while (col_index < n) {
        GET_VECTOR_ELEMENT(vecval, x_i, x_index, $1)
        if ((diag == blas_unit_diag) && (col_index == matrix_row)) {
          MUL(rowtmp, $4, vecval, $1, one, $2)
        } else {
          GET_VECTOR_ELEMENT(matval, tp_i, tp_index, $2)
          MUL(rowtmp, $4, matval, $2, vecval, $1)
        }
        ADD(rowsum, $3, rowsum, $3, rowtmp, $4)
        x_index += incx;
        tp_index+= inctp;
        col_index++;
      }
      MUL(result, $3, rowsum, $3, alpha_i, $1)
      SET_ROUND_VECTOR_ELEMENT(x_i, x_index2, result, $3)
    }
  } else if ((uplo == blas_upper &&
              trans == blas_no_trans && order == blas_colmajor) ||
             (uplo == blas_lower &&
              trans != blas_no_trans && order == blas_rowmajor)) {
dnl ---------------------------------------------------------------------
dnl 1) upper triangular and column major with no trans
dnl 2) lower triangular and row major with trans
dnl ---------------------------------------------------------------------
    tp_start = ((n-1)*n)/ 2;
    inctp2 = n-1;
    x_index2 = x_start;
    for (matrix_row = 0; matrix_row < n; matrix_row++, inctp2 = n-1) {
      x_index = x_start + incx*(n-1);
      tp_index = (tp_start + matrix_row)*inctp;
      col_index = (n-1) - matrix_row;
      ZERO(rowsum, $3)
      ZERO(rowtmp, $4)
      ZERO(result, $3)
      while (col_index >= 0) {
        GET_VECTOR_ELEMENT(vecval, x_i, x_index, $1)
        if ((diag == blas_unit_diag) && (col_index == 0)) {
          MUL(rowtmp, $4, vecval, $1, one, $2)
        } else {
          GET_VECTOR_ELEMENT(matval, tp_i, tp_index, $2)
          MUL(rowtmp, $4, matval, $2, vecval, $1)
        }
        ADD(rowsum, $3, rowsum, $3, rowtmp, $4)
        x_index -= incx;
        tp_index -= inctp2*inctp;
        inctp2--;
        col_index--;
      }
      MUL(result, $3, rowsum, $3, alpha_i, $1)
      SET_ROUND_VECTOR_ELEMENT(x_i, x_index2, result, $3)
      x_index2 += incx;
    }
  } else if ((uplo == blas_lower &&
              trans == blas_no_trans && order == blas_rowmajor) ||
             (uplo == blas_upper &&
              trans != blas_no_trans && order == blas_colmajor)) {
dnl ---------------------------------------------------------------------
dnl 1) lower triangular and row major with no trans
dnl 2) upper triangular and column major with trans
dnl ---------------------------------------------------------------------
    tp_start = (n-1) + ((n-1) * n)/ 2;
    tp_index = tp_start * inctp;
    x_index = x_start + (n-1) * incx;

    for(matrix_row = n-1; matrix_row >=0; matrix_row--) {
      x_index2 = x_index;
      ZERO(rowsum, $3)
      ZERO(rowtmp, $4)
      ZERO(result, $3)
      for(step = 0; step <= matrix_row; step++) {
        GET_VECTOR_ELEMENT(vecval, x_i, x_index2, $1)
        if ((diag == blas_unit_diag) && (step == 0)) {
          MUL(rowtmp, $4, vecval, $1, one, $2)
        } else {
          GET_VECTOR_ELEMENT(matval, tp_i, tp_index, $2)
          MUL(rowtmp, $4, matval, $2, vecval, $1)
        }
        ADD(rowsum, $3, rowsum, $3, rowtmp, $4)
        x_index2 -= incx;
        tp_index-=inctp;
      }
      MUL(result, $3, rowsum, $3, alpha_i, $1)
      SET_ROUND_VECTOR_ELEMENT(x_i, x_index, result, $3)
      x_index -= incx;
    }
  } else {
dnl ---------------------------------------------------------------------
dnl 1) lower triangular and column major with no trans
dnl 2) upper triangular and row major with trans
dnl ---------------------------------------------------------------------
    tp_start = 0;
    x_index = x_start + (n-1) * incx;
    for (matrix_row = n-1; matrix_row >= 0; matrix_row--) {
      tp_index = matrix_row*inctp;
      x_index2 = x_start;
      ZERO(rowsum, $3)
      ZERO(rowtmp, $4)
      ZERO(result, $3)
      stride = n;
      for (step = 0; step <= matrix_row; step++) {
        GET_VECTOR_ELEMENT(vecval, x_i, x_index2, $1)
        if ((diag == blas_unit_diag) && (step == matrix_row)) {
          MUL(rowtmp, $4, vecval, $1, one, $2)
        } else {
          GET_VECTOR_ELEMENT(matval, tp_i, tp_index, $2)
          MUL(rowtmp, $4, matval, $2, vecval, $1)
        }
        ADD(rowsum, $3, rowsum, $3, rowtmp, $4)
        stride--;
        tp_index += stride*inctp;
        x_index2 += incx;
      }
      MUL(result, $3, rowsum, $3, alpha_i, $1)
      SET_ROUND_VECTOR_ELEMENT(x_i, x_index, result, $3)
      x_index -= incx;
    }
  }
}')dnl
dnl
dnl
dnl ---------------------------------------------------------------------
dnl Usage: TPMV_BODY(ax_type, tp_type, sum_type, tmp_type, [FPU])
dnl        perform the top level dispatch for TPMV.  Branches on order
dnl        and on whether alpha is equal to zero.
dnl The last parameter ($5 FPU) is optional.  The string `FPU' is passed
dnl if FPU fix is needed.  Otherwise, an empty string is passed.
dnl ---------------------------------------------------------------------
define(`TPMV_BODY', `
{
  int matrix_row, step, tp_index, tp_start, x_index, x_start;
  int inctp, x_index2, stride, col_index, inctp2;

  SCALAR_CAST(alpha, $1)

  PTR_CAST(tp, $2, `const')
  PTR_CAST(x, $1)
  DECLARE(rowsum, $3)
  DECLARE(rowtmp, $4)
  DECLARE(result, $3)
  DECLARE(matval, $2)
  DECLARE(vecval, $1)
  DECLARE(one, $2)

  ifelse(`$5', `FPU', `FPU_FIX_DECL;')
  ONE(one, $2)

  inctp = 1;
  INC_ADJUST(inctp, $2)
  INC_ADJUST(incx, $1)

  if(incx<0) x_start = (-n+1)*incx; else x_start = 0;

  if (n < 1) {
    return;
  }

  /* Check for error conditions. */
  if (order != blas_colmajor && order != blas_rowmajor) {
    BLAS_error(routine_name, -1, order, NULL);
  }
  if (uplo != blas_upper && uplo != blas_lower) {
    BLAS_error(routine_name, -2, uplo, NULL);
  }
  if (incx == 0) {
    BLAS_error(routine_name, -9, incx, NULL);
  }
  ifelse(`$5', `FPU', `FPU_FIX_START;')

  TPMV_LOOP($1, $2, $3, $4)

  ifelse(`$6', `FPU', `FPU_FIX_STOP;')
}')dnl
dnl
dnl
define(`SWITCH_prec', `
  switch(prec) {
    case blas_prec_single: ifelse(`$3&&$4', `$5&&$6', `', `{
      TPMV_BODY($1, $2, $3, $4)
      break;
    }
    ')dnl
    case blas_prec_double:
    case blas_prec_indigenous: {
      TPMV_BODY($1, $2, $5, $6)
      break;
    }

    case blas_prec_extra: {
      TPMV_BODY($1, $2, $7, $8, FPU)
      break;
    }
  }
')dnl
dnl
dnl
define(`TPMV_X_BODY', `
  SWITCH_prec($1, $2, 
    TMP_TYPE_X($1, S), TMP_TYPE_X($1, S), 
    TMP_TYPE_X($1, D), TMP_TYPE_X($1, D), 
    TMP_TYPE_X($1, E), TMP_TYPE_X($1, E))')dnl
dnl
dnl
define(`TPMV', `
  TPMV_HEAD($1, $2, $3)
  TPMV_COMMENT($1, $2, $3)
  {
    static const char routine_name[] = "TPMV_NAME($1, $2, $3)";
    ifelse($3, _x, `TPMV_X_BODY($1_type, $2_type)', 
      `TPMV_BODY($1_type, $2_type, HIGHER_TYPE($1_type, $2_type), $1_type)')
  }')dnl
dnl
dnl
