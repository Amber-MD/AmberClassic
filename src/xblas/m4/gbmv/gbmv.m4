dnl ----------------------------------------------------------------------
dnl GBMV --- General Banded Matrix Vector Product
dnl   y <--- alpha * op(A) * x + beta * y
dnl   
dnl   where op can be no-op, transpose, or conjugate transpose.
dnl ----------------------------------------------------------------------
dnl
dnl
#include "blas_extended.h" 
#include "blas_extended_private.h"
include(cblas.m4)dnl
include(gbmv-common.m4)dnl
dnl
dnl
dnl  Usage:
dnl    GBMV        ($1, $2, $3, $4)
dnl    GBMV_HEAD   ($1, $2, $3, $4)
dnl    GBMV_NAME   ($1, $2, $3, $4)
dnl    GBMV_PARAMS ($1, $2, $3, $4)
dnl    GBMV_COMMENT($1, $2, $3, $4)
dnl
dnl      $1 -- type of alpha, beta, y.
dnl      $2 -- type of A
dnl      $3 -- type of x
dnl      $4 -- Set to `_x' for _x routines.  Otherwise set to `'.
dnl
dnl
define(`GBMV_COMMENT',`
/*           
 * Purpose
 * =======
 *
 *  gbmv computes y = alpha * A * x + beta * y, where 
 *
 *  A is a m x n banded matrix
 *  x is a n x 1 vector
 *  y is a m x 1 vector
 *  alpha and beta are scalars 
 *
 *   
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of AP; row or column major
 *
 * trans        (input) blas_trans_type
 *              Transpose of AB; no trans, 
 *              trans, or conjugate trans
 *
 * m            (input) int
 *              Dimension of AB
 *
 * n            (input) int
 *              Dimension of AB and the length of vector x
 *
 * kl           (input) int 
 *              Number of lower diagnols of AB
 *
 * ku           (input) int
 *              Number of upper diagnols of AB
 *
 * alpha        (input) $1_scalar
 *              
 * AB           (input) $2_array
 *
 * lda          (input) int 
 *              Leading dimension of AB
 *              lda >= ku + kl + 1
 *
 * x            (input) $3_array
 * 
 * incx         (input) int
 *              The stride for vector x.
 *
 * beta         (input) $1_scalar
 *
 * y            (input/output) $1_array
 *
 * incy         (input) int
 *              The stride for vector y.
 * 
PREC_COMMENT($4)dnl
 *
 * LOCAL VARIABLES 
 * ===============
 * 
 *  As an example, these variables are described on the mxn, column 
 *  major, banded matrix described in section 2.2.3 of the specification  
 *
 *  astart      indexes first element in A where computation begins
 *
 *  incai1      indexes first element in row where row is less than lbound
 * 
 *  incai2      indexes first element in row where row exceeds lbound
 *   
 *  lbound      denotes the number of rows before  first element shifts 
 *
 *  rbound      denotes the columns where there is blank space
 *   
 *  ra          index of the rightmost element for a given row
 *  
 *  la          index of leftmost  elements for a given row
 *
 *  ra - la     width of a row
 *
 *                        rbound 
 *            la   ra    ____|_____ 
 *             |    |   |          |
 *         |  a00  a01   *    *   *
 * lbound -|  a10  a11  a12   *   *
 *         |  a20  a21  a22  a23  *
 *             *   a31  a32  a33 a34
 *             *    *   a42  a43 a44
 *
 *  Varations on order and transpose have been implemented by modifying these
 *  local variables. 
 *
 */')dnl
dnl
dnl
dnl  Usage: GBMV_BODY($1, $2, $3, $4, $5, $6)
dnl  Generates the main body of the product code.
dnl    $1 - type of alpha, beta, r
dnl    $2 - type of a
dnl    $3 - type of b
dnl    $4 - type of sum/prod
dnl    $5 - type of temp
dnl    $6 - [optional] String `FPU' is passed if FPU fix 
dnl         is needed.  Empty string is passed otherwise.
dnl
define(`GBMV_BODY', `
  int ky, iy, kx, jx, j, i, rbound, lbound, ra, la, lenx, leny;
  int incaij, aij, incai1,incai2,astart,ai;
  PTR_CAST(y, $1)
  PTR_CAST(a, $2,`const')
  PTR_CAST(x, $3,`const')
  SCALAR_CAST(alpha, $1)
  SCALAR_CAST(beta, $1)
  DECLARE(tmp1, $5)
  DECLARE(tmp2, $5)
  DECLARE(result, $1)
  DECLARE(sum, $4)
  DECLARE(prod, $4)
  DECLARE(a_elem, $2)
  DECLARE(x_elem, $3)
  DECLARE(y_elem, $1)
  ifelse(`$6', `FPU', `FPU_FIX_DECL;')

  if (order != blas_colmajor && order != blas_rowmajor)
    BLAS_error(routine_name, -1, order, NULL);
  if (trans != blas_no_trans &&
      trans != blas_trans &&
      trans != blas_conj_trans) {
    BLAS_error(routine_name, -2, trans, NULL);
  }
  if (m < 0)
    BLAS_error(routine_name, -3, m, NULL);
  if (n < 0)
    BLAS_error(routine_name, -4, n, NULL);
  if (kl < 0 || kl >= m)
    BLAS_error(routine_name, -5, kl, NULL);
  if (ku < 0 || ku >= n)
    BLAS_error(routine_name, -6, ku, NULL);
  if (lda < kl + ku + 1)
    BLAS_error(routine_name, -9, lda, NULL);
  if (incx == 0)
    BLAS_error(routine_name, -11, incx, NULL);
  if (incy == 0)
    BLAS_error(routine_name, -14, incy, NULL);

  if ((m == 0) || (n == 0) || 
      (((TEST_0(alpha_i, $1)) && (TEST_1(beta_i, $1)))))
    return;

  if (trans == blas_no_trans) {
    lenx = n;
    leny = m;
  } else {      /* change back */
    lenx = m;
    leny = n;
  }

  if ( incx < 0 ) {
    kx = - ( lenx - 1 ) * incx;
  } else {
    kx = 0;
  }
  
  if ( incy < 0 ) {
    ky = - ( leny - 1 ) * incy;
  } else {
    ky = 0;
  }

  ifelse(`$6', `FPU', `FPU_FIX_START;') 

  /* if alpha = 0, return y = y*beta */
  if ((order == blas_colmajor) && (trans == blas_no_trans)) {
    astart = ku;
    incai1 = 1;
    incai2 = lda;
    incaij = lda-1;
    lbound = kl;
    rbound = n-ku-1;
    ra = ku;
  } else if ((order == blas_colmajor) && (trans != blas_no_trans)) {
    astart = ku;
    incai1 = lda-1;
    incai2 = lda;
    incaij = 1;
    lbound = ku;
    rbound = m-kl-1;
    ra = kl;
  } else if ((order == blas_rowmajor) && (trans == blas_no_trans)) {
    astart = kl; 
    incai1 = lda-1; 
    incai2 = lda; 
    incaij = 1; 
    lbound = kl; 
    rbound = n-ku-1; 
    ra = ku;
  } else {    /* rowmajor and blas_trans */
    astart = kl; 
    incai1 = 1; 
    incai2 = lda; 
    incaij = lda-1; 
    lbound = ku; 
    rbound = m-kl-1; 
    ra = kl;
  }
  INC_ADJUST(incx, $3)
  INC_ADJUST(incy, $1)
  INC_ADJUST(incaij, $2)
  INC_ADJUST(incai1, $2)
  INC_ADJUST(incai2, $2)
  INC_ADJUST(astart,$2)
  INC_ADJUST(ky,$1)
  INC_ADJUST(kx, $3)

  la = 0;
  ai = astart;
  iy = ky;
  for (i = 0; i < leny; i++) {
    ZERO(sum, $4)
    aij = ai;
    jx = kx;
    ifelse(
      `$2', `complex_S', `if(trans != blas_conj_trans) {',
      `$2', `complex_D', `if(trans != blas_conj_trans) {')
    for (j = ra - la; j >= 0; j--) {
      GET_VECTOR_ELEMENT(x_elem, x_i, jx, $3)
      GET_VECTOR_ELEMENT(a_elem, a_i, aij, $2)
      MUL(prod, $4, x_elem, $3, a_elem, $2)
      ADD(sum, $4, sum, $4, prod, $4)
      aij += incaij;
      jx += incx;
    }       
    ifelse(
      `$2', `complex_S',  `
      } else {
        for (j = ra - la; j >= 0; j--) {
          GET_VECTOR_ELEMENT(x_elem, x_i, jx, $3)
          GET_VECTOR_ELEMENT(a_elem, a_i, aij, $2)
          CONJ_AUX(a_elem, $2)
          MUL(prod, $4, x_elem, $3, a_elem, $2)
          ADD(sum, $4, sum, $4, prod, $4)
          aij += incaij;
          jx += incx;
        }
      }',
      `$2', `complex_D', `
      } else {
        for (j = ra - la; j >= 0; j--) {
          GET_VECTOR_ELEMENT(x_elem, x_i, jx, $3)
          GET_VECTOR_ELEMENT(a_elem, a_i, aij, $2)
          CONJ_AUX(a_elem, $2)
          MUL(prod, $4, x_elem, $3, a_elem, $2)
          ADD(sum, $4, sum, $4, prod, $4)
          aij += incaij;
          jx += incx;
      }
    }')
    
    MUL(tmp1, $5, sum, $4, alpha_i, $1)
    GET_VECTOR_ELEMENT(y_elem, y_i, iy, $1)
    MUL(tmp2, $5, beta_i, $1, y_elem, $1)
    ADD(result, `$1', tmp1, `$5', tmp2, `$5')
    SET_ROUND_VECTOR_ELEMENT(y_i, iy, result, $1)
    iy += incy;
    if (i >=lbound) {
      kx +=incx;
      ai += incai2;
      la++;
    } else {
      ai += incai1;
    }
    if (i < rbound) {
      ra++;
    }
  }

  ifelse(`$6', `FPU', `FPU_FIX_STOP;')
')dnl
dnl
dnl
dnl  Usage: GBMV_X_BODY($1, $2, $3)
dnl  Generates the main body of the extended version of gbmv code.
dnl    $1 -- type of alpha, beta, r
dnl    $2 -- type of a
dnl    $3 -- type of b
dnl
define(`GBMV_X_BODY', 
  `SWITCH_prec($1, $2, $3, 
    SUM_TYPE_X($1, $2, $3, S), TMP_TYPE_X($1, S), 
    SUM_TYPE_X($1, $2, $3, D), TMP_TYPE_X($1, D), 
    SUM_TYPE_X($1, $2, $3, E), TMP_TYPE_X($1, E))')dnl
dnl
dnl
dnl  Usage: SWITCH_prec($1, $2, $3, $4, $5, $6, $7, $8, $9)
dnl  Generates a 3-way switch statement based on prec.
dnl    $1      --  type of alpha, beta, r
dnl    $2      --  type of a
dnl    $3      --  type of b
dnl    $4, $5  --  type of `sum' and `tmp' in single case
dnl    $6, $7  --  type of `sum' and `tmp' in double/indigenous case
dnl    $8, $9  --  type of `sum' and `tmp' in extra case
dnl
define(`SWITCH_prec', `
  switch ( prec ) {
  case blas_prec_single: ifelse(`$4&&$5', `$6&&$7', `', `{
    GBMV_BODY($1, $2, $3, $4, $5)
    break;
  }
  ')dnl
  case blas_prec_double:
  case blas_prec_indigenous:
    { GBMV_BODY($1, $2, $3, $6, $7) }
    break;
  case blas_prec_extra:
    { GBMV_BODY($1, $2, $3, $8, $9, FPU) }
    break;
 }')dnl
dnl
dnl
define(`GBMV', 
  `GBMV_HEAD($1, $2, $3, $4)
   GBMV_COMMENT($1, $2, $3, $4)
   {
     static const char routine_name[] = "GBMV_NAME($1, $2, $3, $4)";
     ifelse($4, _x, `GBMV_X_BODY($1_type, $2_type, $3_type)', 
       `GBMV_BODY($1_type, $2_type, $3_type, 
         SUM_TYPE($1_type, $2_type, $3_type), $1_type)')
   } /* end GEMV_NAME($1, $2, $3, $4) */
')dnl
