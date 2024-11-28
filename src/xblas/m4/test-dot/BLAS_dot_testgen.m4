dnl **********************************************************************
dnl Generate routines to test DOT.
dnl **********************************************************************
dnl
include(cblas.m4)dnl
include(test-common.m4)dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: DOT_TESTGEN_BODY(abr_type, x_type, y_type, int_typeltr)
dnl        abr_type    : the type and precision of alpha, beta and r
dnl      x_type      : the type and precision of x
dnl      y_type      : the type and precision of y
dnl        int_typeltr : the internal precision used to generate test vectors;
dnl                      it can be one of s, d, c, z.
dnl ----------------------------------------------------------------------
define(`DOT_TESTGEN_BODY',
`{  
  int i, inc = 1;
  PTR_CAST(alpha, $1)
  PTR_CAST(beta, $1)
  PTR_CAST(r, $1)
  PTR_CAST(x, $2)
  PTR_CAST(y, $3)
  DECLARE(alpha_tmp, $4_type)
  DECLARE(beta_tmp, $4_type)
  DECLARE(r_tmp, $4_type)
  DECLARE_VECTOR(x_vec, $4_type)
  DECLARE_VECTOR(y_vec, $4_type)

  ASSIGN_PTR_TO_SCALAR(alpha_tmp, $4_type, alpha_i, $1)
  ASSIGN_PTR_TO_SCALAR(beta_tmp, $4_type, beta_i, $1)
  INC_ADJUST(inc, $4_type)
  MALLOC_VECTOR(x_vec, $4_type, 2*n)
  y_vec = x_vec + inc*n;
  for (i = 0; i < inc*n_fix2; i += inc) {
    COPY_VECTOR_ELEMENT2(x_vec, $4_type, x_i, $2, i)
    COPY_VECTOR_ELEMENT2(y_vec, $4_type, y_i, $3, i)
  }
  for (; i < inc*(n_fix2 + n_mix); i += inc) {
    COPY_VECTOR_ELEMENT2(x_vec, $4_type, x_i, $2, i)
  }

  /* Call generator now. */
  testgen_BLAS_$4dot(n, n_fix2, n_mix, norm,  conj,  
      PASS_BY_REF(alpha_tmp, $4_type), alpha_flag,
      PASS_BY_REF(beta_tmp, $4_type), beta_flag,
      x_vec, y_vec, seed, PASS_BY_REF(r_tmp, $4_type), r_true_l, r_true_t);

  ASSIGN_SCALAR_TO_PTR(alpha_i, $1, alpha_tmp, $4_type)
  ASSIGN_SCALAR_TO_PTR(beta_i, $1, beta_tmp, $4_type)
  ASSIGN_SCALAR_TO_PTR(r_i, $1, r_tmp, $4_type)
  for (i = 0; i < inc*n; i += inc) {
    COPY_VECTOR_ELEMENT2(x_i, $2, x_vec, $4_type, i)
    COPY_VECTOR_ELEMENT2(y_i, $3, y_vec, $4_type, i)
  }

  FREE_VECTOR(x_vec, $4_type) /* also y_i */
}')dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: MIXED_DOT_TESTGEN_BODY(abr_type, x_type, y_type, int_typeltr)
dnl        ... mixed real and complex types
dnl        abr_type    : the type and precision of alpha, beta and r
dnl        x_type      : the type and precision of x
dnl      y_type      : the type and precision of y
dnl        int_typeltr : the internal precision used to generate test vectors;
dnl                      it can be one of s and d.
dnl
dnl Algorithm description:
dnl     We always call the real test generator to get real test inputs first.
dnl     Then we multiply them by appropriate complex numbers of the
dnl     form (2^j) + (2^k)i in a way that the final result is the product of
dnl     this complex number and the real result.
dnl     We should have small magnitude in both real and imaginary parts.
dnl  
dnl Example 1
dnl     (abr_type, x_type, y_type) = (c, s, s)
dnl     After getting the real alpha, beta and r, we can apply one of the
dnl     following complex scalings:
dnl             alpha   beta    r
dnl             -----   ----    -
dnl             1+i     1-i     i
dnl            -1+i     1+i     i
dnl             2i      1+i     1+i
dnl             2       1+i     1-i
dnl     The first (or second) one is the best.  
dnl
dnl Example 2
dnl     (abr_type, x_type, y_type) = (c, c, s)
dnl     After getting the real alpha, beta, r, x and y, we can apply the
dnl     following complex scalings:
dnl             alpha   x       y       beta    r
dnl             -----   ----    ----    ----    ----
dnl             1+i     1+i     1       1+i     1+i  ==> product is 2i
dnl ----------------------------------------------------------------------
define(`MIXED_DOT_TESTGEN_BODY',
`{  
  int i;
  PTR_CAST(alpha, $1)
  PTR_CAST(beta, $1)
  PTR_CAST(r, $1)
  PTR_CAST(x, $2)
  PTR_CAST(y, $3)
  DECLARE(alpha_i_r, $4_type)
  DECLARE(alpha_i_i, $4_type)
  DECLARE(beta_i_r, $4_type)
  DECLARE(beta_i_i, $4_type)
  DECLARE(r_tmp, $4_type)
  DECLARE_VECTOR(x_vec, $4_type)
  DECLARE_VECTOR(y_vec, $4_type)

  GET_REAL_PART(alpha_i_r, $4_type, alpha_i, $1)
  GET_IMAG_PART(alpha_i_i, $4_type, alpha_i, $1)
  GET_REAL_PART(beta_i_r, $4_type, beta_i, $1)
  GET_IMAG_PART(beta_i_i, $4_type, beta_i, $1)
  MALLOC_VECTOR(x_vec, $4_type, 2*n)
  y_vec = x_vec + n;
  for (i = 0; i < n_fix2; i++) {
    COPY_VECTOR_ELEMENT2(x_vec, $4_type, x_i, $2, i)
    COPY_VECTOR_ELEMENT2(y_vec, $4_type, y_i, $3, i)
  }
  for (; i < n_fix2 + n_mix; i++) {
    COPY_VECTOR_ELEMENT2(x_vec, $4_type, x_i, $2, i)
  }

  /* Call generator now. */
  testgen_BLAS_$4dot(n, n_fix2, n_mix, norm,  conj,  
      PASS_BY_REF(alpha_i_r, $4_type), alpha_flag,
      PASS_BY_REF(beta_i_r, $4_type), beta_flag,
      x_vec, y_vec, seed, PASS_BY_REF(r_tmp, $4_type),
      &r_true_l[0], &r_true_t[0]);

  if ( alpha_flag == 1 ) { /* alpha_i is fixed */
    if ( alpha_i_r == 1.0 && alpha_i_i == 0. ) { /* alpha_i == 1.0 */
      ifelse(`$2', `$3',  dnl ... both x_i and y_i real
        `if (beta_flag == 1 && 
             ((beta_i_r == 0. && beta_i_i == 0.) ||
              (beta_i_r == 1. && beta_i_i == 0.))) { /* beta_i == 0 or 1 */
           PUT_REAL_PART(r_i, $1, r_tmp, $4_type)
           PUT_IMAG_PART(r_i, $1, 0.0, $4_type)
         } else { /* beta_i *= (1-i), r_i *= (1+i)/2 --> prod = 1 */
           PUT_REAL_PART(beta_i, $1, beta_i_r, $4_type)
           PUT_IMAG_PART(beta_i, $1, -beta_i_r, $4_type)
           PUT_REAL_PART(r_i, $1, r_tmp/2., $4_type)
           PUT_IMAG_PART(r_i, $1, r_tmp/2., $4_type)
         }
         r_true_l[1] = r_true_t[1] = 0.0;',
        `                 dnl ... one of x_i and y_i is complex
         if (beta_flag == 1 && 
            ((beta_i_r == 0. && beta_i_i == 0.) || 
             (beta_i_r == 1. && beta_i_i == 0.))) {
           /* beta_i == 0 or 1 --> r_i *= (1+i) */
           PUT_REAL_PART(r_i, $1, r_tmp, $4_type)
           PUT_IMAG_PART(r_i, $1, r_tmp, $4_type)
         } else { /* beta_i *= (1-i), r_i *= i --> prod = 1+i */
           PUT_REAL_PART(beta_i, $1, beta_i_r, $4_type)
           PUT_IMAG_PART(beta_i, $1, -beta_i_r, $4_type)
           PUT_REAL_PART(r_i, $1, 0.0, $4_type)
           PUT_IMAG_PART(r_i, $1, r_tmp, $4_type)
         }
         r_true_l[1] = r_true_l[0];
         r_true_t[1] = r_true_t[0];')dnl
    } else if (alpha_i_r == 0. && alpha_i_i == 0.) { /* alpha_i == 0.0 */
      if (beta_flag == 1 && 
          ((beta_i_r == 0. && beta_i_i == 0.) ||
           (beta_i_r == 1. && beta_i_i == 0.))) {
        /* beta_i == 0 or 1 --> r_i *= (1+i) */
        PUT_REAL_PART(r_i, $1, r_tmp, $4_type)
        PUT_IMAG_PART(r_i, $1, r_tmp, $4_type)
      } else { /* beta_i *= (1-i), r_i *= i --> prod = 1+i */
        PUT_REAL_PART(beta_i, $1, beta_i_r, $4_type)
        PUT_IMAG_PART(beta_i, $1, -beta_i_r, $4_type)
        PUT_REAL_PART(r_i, $1, 0.0, $4_type)
        PUT_IMAG_PART(r_i, $1, r_tmp, $4_type)
      }
      r_true_l[1] = r_true_l[0];
      r_true_t[1] = r_true_t[0];
    } else { /* alpha_i is a fixed multiple of (1+i) */
      PUT_REAL_PART(alpha_i, $1, alpha_i_r, $4_type)
      PUT_IMAG_PART(alpha_i, $1, alpha_i_r, $4_type)
      ifelse(`$2', `$3',  dnl ... both x_i and y_i real
        `if (beta_flag == 1 &&   
             ((beta_i_r == 0. && beta_i_i == 0.) ||
              (beta_i_r == 1. && beta_i_i == 0.))) {
           /* beta_i == 0 or 1 --> r_i *= (1+i) */
           PUT_REAL_PART(r_i, $1, r_tmp, $4_type)
           PUT_IMAG_PART(r_i, $1, r_tmp, $4_type)
         } else { /* beta_i *= (1-i), r_i *= i --> prod = 1+i */
           PUT_REAL_PART(beta_i, $1, beta_i_r, $4_type)
           PUT_IMAG_PART(beta_i, $1, -beta_i_r, $4_type)
           PUT_REAL_PART(r_i, $1, 0.0, $4_type)
           PUT_IMAG_PART(r_i, $1, r_tmp, $4_type)
         }
         r_true_l[1] = r_true_l[0];
         r_true_t[1] = r_true_t[0];',
        `                dnl ... one of x_i and y_i is complex
         if (beta_flag == 1 && 
             ((beta_i_r == 0. && beta_i_i == 0.) ||
              (beta_i_r == 1. && beta_i_i == 0.))) {
           /* beta_i is 0 or 1 --> r_i *= 2i --> prod = 2i */
           PUT_REAL_PART(r_i, $1, 0.0, $4_type)
           PUT_IMAG_PART(r_i, $1, 2.0*r_tmp, $4_type)
         } else { /* beta_i *= (1+i), r_i *= (1+i) --> prod = 2i */
           PUT_REAL_PART(beta_i, $1, beta_i_r, $4_type)
           PUT_IMAG_PART(beta_i, $1, beta_i_r, $4_type)
           PUT_REAL_PART(r_i, $1, r_tmp, $4_type)
           PUT_IMAG_PART(r_i, $1, r_tmp, $4_type)
         }
         r_true_l[1] = 2.0 * r_true_l[0];
         r_true_t[1] = 2.0 * r_true_t[0];
         r_true_l[0] = r_true_t[0] = 0.0;')dnl
    }
  } else if ( beta_flag == 1 ) { /* alpha_i is free, beta_i is fixed */
    /* alpha_i *= (1+i) */
    PUT_REAL_PART(alpha_i, $1, alpha_i_r, $4_type)
    PUT_IMAG_PART(alpha_i, $1, alpha_i_r, $4_type)
    ifelse(`$2', `$3',  dnl ... both x_i and y_i real
      `if ((beta_i_r == 0. && beta_i_i == 0.) ||
           (beta_i_r == 1. && beta_i_i == 0.)) { /* r_i *= (1+i) */
         PUT_REAL_PART(r_i, $1, r_tmp, $4_type)
         PUT_IMAG_PART(r_i, $1, r_tmp, $4_type)
       } else { /* beta_i *= (1-i), r_i *= i */
         PUT_REAL_PART(beta_i, $1, beta_i_r, $4_type)
         PUT_IMAG_PART(beta_i, $1, -beta_i_r, $4_type)
         PUT_REAL_PART(r_i, $1, 0., $4_type)
         PUT_IMAG_PART(r_i, $1, r_tmp, $4_type)
       }
       r_true_l[1] = r_true_l[0];
       r_true_t[1] = r_true_t[0];',
      `                 dnl ... one of x_i and y_i is complex
       if ((beta_i_r == 0. && beta_i_i == 0.) ||
           (beta_i_r == 1. && beta_i_i == 0.)) { /* r_i*=2i --> prod = 2i */
         PUT_REAL_PART(r_i, $1, 0.0, $4_type)
         PUT_IMAG_PART(r_i, $1, 2.0*r_tmp, $4_type)
       } else { /* beta_i *= (1+i), r_i *= (1+i) */
         PUT_REAL_PART(beta_i, $1, beta_i_r, $4_type)
         PUT_IMAG_PART(beta_i, $1, beta_i_r, $4_type)
         PUT_REAL_PART(r_i, $1, r_tmp, $4_type)
         PUT_IMAG_PART(r_i, $1, r_tmp, $4_type)
       }
       r_true_l[1] = 2.0 * r_true_l[0];
       r_true_t[1] = 2.0 * r_true_t[0];
       r_true_l[0] = r_true_t[0] = 0.0;')dnl
  } else { /* both alpha_i and beta_i are free */
    assert(alpha_flag==0 && beta_flag==0);
    PUT_REAL_PART(alpha_i, $1, alpha_i_r, $4_type)
    PUT_IMAG_PART(alpha_i, $1, alpha_i_r, $4_type)
    ifelse(`$2', `$3',           dnl ... both x_i and y_i real
      `PUT_REAL_PART(beta_i, $1, beta_i_r, $4_type)
       PUT_IMAG_PART(beta_i, $1, -beta_i_r, $4_type)
       PUT_REAL_PART(r_i, $1, 0, $4_type)
       PUT_IMAG_PART(r_i, $1, r_tmp, $4_type)
       /* imaginary part of r_true */
       r_true_l[1] = r_true_l[0];
       r_true_t[1] = r_true_t[0];',
      `                              dnl ... one of x_i and y_i is complex
       PUT_REAL_PART(beta_i, $1, beta_i_r, $4_type)
       PUT_IMAG_PART(beta_i, $1, beta_i_r, $4_type)
       PUT_REAL_PART(r_i, $1, r_tmp, $4_type)
       PUT_IMAG_PART(r_i, $1, r_tmp, $4_type)
       /* imaginary part of r_true */
       ddmuld(r_true_l[0], r_true_t[0], 2.0, &r_true_l[1], &r_true_t[1]);
       /* real part of r_true */
       r_true_l[0] = 0.;
       r_true_t[0] = 0.;')dnl ...
  }
  for (i = 0; i < n; ++i) {
    COPY_VECTOR_ELEMENT2(x_i, $2, x_vec, $4_type, i)
    COPY_VECTOR_ELEMENT2(y_i, $3, y_vec, $4_type, i)
  }
  IF_COMPLEX($2, 
    `if (conj == blas_conj) {
        for (i = 0; i < n; ++i) x_i[2*i+1] = -x_i[2*i+1];
    }')

  FREE_VECTOR(x_vec, $4_type) /* also y_vec */
}')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: COPY_VECTOR_ELEMENT2(x, x_type, y, y_type, i)  ... x[i] = y[i]
dnl Scales by (1 + i) when copying real into complex, copies only the real
dnl part when copying complex into real.
dnl ----------------------------------------------------------------------
define(`COPY_VECTOR_ELEMENT2', `ifelse( 
  REAL_COMPLEX($2)&&REAL_COMPLEX($4), `complex&&complex', 
    `$1[$5] = $3[$5]; $1[$5+1] = $3[$5+1];',
  REAL_COMPLEX($2)&&REAL_COMPLEX($4), `complex&&real', 
    `$1[2*$5] = $3[$5]; $1[2*$5+1] = $3[$5];', 
  REAL_COMPLEX($2)&&REAL_COMPLEX($4), `real&&complex', 
    `$1[$5] = $3[2*$5];', 
  REAL_COMPLEX($2)&&REAL_COMPLEX($4), `real&&real', 
    `$1[$5] = $3[$5];', `
#error Unimplemented case in COPY_VECTOR_ELEMENT2')')dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: DOT_TESTGEN_COMMENT(abr_typeltr, x_typeltr, y_typeltr)
dnl        ... generate the leading comment for the dot product
dnl            testing routines.
dnl Each _typeltr specifier can be one of
dnl        s   ... real and single
dnl     d   ... real and double
dnl     c  ... complex and single
dnl     z  ... complex and double
dnl ----------------------------------------------------------------------
dnl
define(`DOT_TESTGEN_COMMENT',`
/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to DOT_NAME_PREFIX($1,$2,$3){_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) $1_array
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) $1_array
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) $2_array
 *
 * y       (input/output) $3_array
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) $1_array
 *         The generated scalar r that will be used as an input to DOT.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: DOT_NAME_PREFIX(abr_typeltr, x_typeltr, y_typeltr)
dnl        ... generate the prefix of the dot product routines.
dnl Each type letter can be one of
dnl        s   ... real and single
dnl     d   ... real and double
dnl     c  ... complex and single
dnl     z  ... complex and double
dnl ----------------------------------------------------------------------
dnl
define(`DOT_NAME_PREFIX', `ifelse(
     `$2&&$3', `$1&&$1',`BLAS_$1dot',
  `BLAS_$1dot_$2_$3')')dnl
dnl
dnl 
define(`DOT_TESTGEN_HEAD', 
  `void DOT_NAME_PREFIX($1, $2, $3)_testgen(int n, int n_fix2, int n_mix, dnl
      int norm, enum blas_conj_type conj, $1_array alpha, int alpha_flag, dnl
      $1_array beta,  int beta_flag, $2_array x, $3_array y, int *seed, dnl
      $1_array r, double *r_true_l, double *r_true_t)')dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: DOT_TESTGEN(abr_typeltr, x_typeltr, y_typeltr)
dnl        ... generate the top level testing routines for dot product.
dnl Each type letter can be one of
dnl        s   ... real and single
dnl     d   ... real and double
dnl     c  ... complex and single
dnl     z  ... complex and double
dnl ----------------------------------------------------------------------
dnl
define(`DOT_TESTGEN',
 `DOT_TESTGEN_HEAD($1, $2, $3)
  DOT_TESTGEN_COMMENT($1, $2, $3)
  ifelse(IS_MIXED($1, $2, $3), `t', 
    `MIXED_DOT_TESTGEN_BODY($1_type, $2_type, $3_type, REAL_ABBREV($1))', 
    `DOT_TESTGEN_BODY($1_type, $2_type, $3_type, 
      ABBREV(REAL_COMPLEX($1_type)_`'LOWER_PREC(PREC($1_type), PREC($2_type), PREC($3_type))))')dnl
  /* end DOT_NAME_PREFIX($1,$2,$3)_testgen */
')dnl
dnl
dnl
define(`PROTOTYPES', `FOREACH(`PREC_ARGS', `
DOT_TESTGEN_HEAD(arg);')')dnl
dnl
dnl
define(`SOURCE', `dnl
#include <stdio.h>
#include <assert.h>
#include "blas_extended.h"
#include "blas_extended_test.h"

FOREACH(`PREC_ARGS', `
DOT_TESTGEN(arg)')')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
