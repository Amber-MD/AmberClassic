dnl **********************************************************************
dnl * Generates alpha, T, and x, where T is a triangular  matrix;        *
dnl * and computes r_true.                                               *
dnl **********************************************************************
dnl
dnl
include(cblas.m4)dnl
include(test-common.m4)dnl
dnl
dnl
define(`TRMV_TESTGEN_COMMENT', `
/*
 * Purpose
 * =======
 *
 * Generates alpha, T and x, where T is a triangular matrix; and 
 * computes r_true.
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
 *              Dimension of AP and the length of vector x
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
 */')dnl
dnl
dnl
dnl
define(`TRMV_TESTGEN_NAME', `ifelse(`$1', `$2', 
  `BLAS_$1trmv_testgen', `BLAS_$1trmv_$2_testgen')')dnl
dnl
dnl
define(`TRMV_TESTGEN_PARAMS', 
  `int norm, enum blas_order_type order, dnl
   enum blas_uplo_type uplo, enum blas_trans_type trans, dnl
   enum blas_diag_type diag, int n, $1_array alpha, dnl
   int alpha_flag, $2_array T, int ldt, $1_array x, dnl
   int *seed, double *HEAD(r_true), double *TAIL(r_true)')dnl
dnl
dnl
define(`TRMV_TESTGEN_HEAD', 
`void TRMV_TESTGEN_NAME($1, $2)(TRMV_TESTGEN_PARAMS($1, $2))')dnl
dnl
dnl
define(`TRMV_TESTGEN', `
TRMV_TESTGEN_HEAD($1, $2) 
TRMV_TESTGEN_COMMENT($1, $2)
TRMV_TESTGEN_BODY($1, $2)  /* end of TRMV_TESTGEN_NAME($1, $2) */
')dnl
dnl
dnl
define(`TRMV_TESTGEN_BODY', `
{
  PTR_CAST(x, $1_type)
  PTR_CAST(T, $2_type)
  PTR_CAST(alpha, $1_type)
  DECLARE_VECTOR(x_vec, $1_type)
  DECLARE_VECTOR(t_vec, $2_type)
  DECLARE(beta, $1_type)
  DECLARE(r, $1_type)
  DECLARE(r_true_elem, EXTRA_TYPE($1_type))
  DECLARE(x_elem, $1_type)
  DECLARE(t_elem, $2_type)

  int inc_tvec = 1, inc_xvec = 1;
  int xvec_i, tvec_j;
  int xi;
  int ti, tij;
  int inc_ti, inc_tij;
  int inc_xi;
  int i, j;

  ZERO(r, $1_type)
  ZERO(beta, $1_type)

  INC_ADJUST(inc_tvec, $2_type)
  INC_ADJUST(inc_xvec, $1_type)

  MALLOC_VECTOR(t_vec, $2_type, n);
  MALLOC_VECTOR(x_vec, $1_type, n);

  if (trans == blas_no_trans) {
    if (uplo == blas_upper) {
      inc_xi = -1;
      if (order == blas_rowmajor) {
        inc_ti = -ldt;
        inc_tij = -1;
      } else {
        inc_ti = -1;
        inc_tij = -ldt;
      }
    } else {
      inc_xi = 1;
      if (order == blas_rowmajor) {
        inc_ti = ldt;
        inc_tij = 1;
      } else {
        inc_ti = 1;
        inc_tij = ldt;
      }
    }
  } else {
    if (uplo == blas_upper) {
      inc_xi = 1;
      if (order == blas_rowmajor) {
        inc_ti = 1;
        inc_tij = ldt;
      } else {
        inc_ti = ldt;
        inc_tij = 1;
      }
    } else {
      inc_xi = -1;
      if (order == blas_rowmajor) {
        inc_ti = -1;
        inc_tij = -ldt;
      } else {
        inc_ti = -ldt;
        inc_tij = -1;
      }
    }
  }

  INC_ADJUST(inc_xi, $1_type)

  INC_ADJUST(inc_ti, $2_type)
  INC_ADJUST(inc_tij, $2_type)

  /* Call dot_testgen n times.  Each call will generate
   * one row of T and one element of x.                   */
  ti = (inc_ti > 0 ? 0 : -(n-1) * inc_ti);
  xi = (inc_xi > 0 ? 0 : -(n-1) * inc_xi);
  xvec_i = 0;
  for (i = 0; i < n; i++) {
    
    /* Generate the i-th element of x_vec and all of t_vec. */
    if (diag == blas_unit_diag) {
      /* Since we need alpha = beta, we fix alpha if alpha_flag = 0. */
      if (i == 0 && alpha_flag == 0) {
        RAND_PTR(alpha_i, $1_type)
      }
      DOT_TESTGEN_NAME($1, $1, $2)(i, 0, i, norm, blas_no_conj, alpha_i, 
          1, alpha_i, 1, x_vec, t_vec, 
          seed, PASS_BY_REF(r, $1_type), 
          PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));

      /* Copy generated t_vec to T. */
      tvec_j = 0;
      tij = (inc_tij > 0 ? ti : ti - (n-1) * inc_tij);
      for (j = 0; j < i; j++) {
        GET_VECTOR_ELEMENT(t_elem, t_vec, tvec_j, $2_type)
        IF_COMPLEX($2_type, `
          if (trans == blas_conj_trans) {
            CONJ_AUX(t_elem, $2_type)
          }', `')
        SET_VECTOR_ELEMENT(T_i, tij, t_elem, $2_type)
        tvec_j += inc_tvec;
        tij += inc_tij;
      }

      /* Set the diagonal element to 1. */
      ONE(t_elem, $2_type)
      SET_VECTOR_ELEMENT(T_i, tij, t_elem, $2_type)

      /* Set x[i] to be r. */
      SET_VECTOR_ELEMENT(x_i, xi, r, $1_type)
      SET_VECTOR_ELEMENT(x_vec, xvec_i, r, $1_type)

    } else {
      DOT_TESTGEN_NAME($1, $1, $2)(i+1, 0, i, norm, blas_no_conj, alpha, 
          (i == 0 ? alpha_flag : 1), PASS_BY_REF(beta, $1_type), 1, x_vec, t_vec, 
          seed, PASS_BY_REF(r, $1_type), 
          PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));

      /* Copy generated t_vec to T. */
      tvec_j = 0;
      tij = (inc_tij > 0 ? ti : ti - (n-1) * inc_tij);
      for (j = 0; j <= i; j++) {
        GET_VECTOR_ELEMENT(t_elem, t_vec, tvec_j, $2_type)
        IF_COMPLEX($2_type, `
          if (trans == blas_conj_trans) {
            CONJ_AUX(t_elem, $2_type)
          }', `')
        SET_VECTOR_ELEMENT(T_i, tij, t_elem, $2_type)
        tvec_j += inc_tvec;
        tij += inc_tij;
      }

      /* Copy generated x_vec[i] to appropriate position in x. */
      GET_VECTOR_ELEMENT(x_elem, x_vec, xvec_i, $1_type)
      SET_VECTOR_ELEMENT(x_i, xi, x_elem, $1_type)
    }

    /* Copy r_true */
    SET_VECTOR_ELEMENT(r_true, xi, r_true_elem, EXTRA_TYPE($1_type))

    xvec_i += inc_xvec;
    xi += inc_xi;

    ti += inc_ti;
  }

  FREE_VECTOR(x_vec, $1_type)
  FREE_VECTOR(t_vec, $2_type)
}
')dnl
dnl
dnl
define(`PROTOTYPES', `dnl
TRMV_TESTGEN_HEAD(s, s);
TRMV_TESTGEN_HEAD(d, d);
TRMV_TESTGEN_HEAD(d, s);
TRMV_TESTGEN_HEAD(c, c);
TRMV_TESTGEN_HEAD(z, c);
TRMV_TESTGEN_HEAD(z, z);
TRMV_TESTGEN_HEAD(c, s);
TRMV_TESTGEN_HEAD(z, d);
')dnl
dnl
dnl
define(`SOURCE', `dnl
#include "blas_extended.h"
#include "blas_extended_test.h"

TRMV_TESTGEN(s, s)
TRMV_TESTGEN(d, d)
TRMV_TESTGEN(d, s)
TRMV_TESTGEN(c, c)
TRMV_TESTGEN(z, c)
TRMV_TESTGEN(z, z)
TRMV_TESTGEN(c, s)
TRMV_TESTGEN(z, d)
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
