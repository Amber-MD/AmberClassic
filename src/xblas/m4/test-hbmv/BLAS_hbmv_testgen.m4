dnl
dnl hbmv_testgen.m4
dnl
dnl Test case generator for hbmv routines.
dnl Generates test cases for alpha, A, beta, x and y and
dnl computes r_true in double-double precision.
dnl
dnl
dnl C Interface has the form 
dnl
dnl   BLAS_xhbmv{_a_b}{_x}_testgen(int norm, 
dnl                             enum blas_order_type order, 
dnl                             enum blas_uplo_type uplo, 
dnl                             int n, int randomize, 
dnl                             SCALAR *alpha, int alpha_flag, 
dnl                             SCALAR *beta, int beta_flag, 
dnl                             ARRAY a, int k, int lda 
dnl                             ARRAY x, int incx, 
dnl                             ARRAY y, int incy, 
dnl                             int *seed, 
dnl                             double *HEAD(r_true), double *TAIL(r_true))
dnl
dnl Arguments
dnl   norm      (in) int
dnl
dnl   order     (in) blas_order_type
dnl                    determines the storage format for the A matrix
dnl   
dnl   uplo      (in) blas_uplo_type
dnl                    determines whether the upper triangular portion
dnl                    or the lower triangular portion of the hermitian
dnl                    matrix A is used.
dnl
dnl   n         (in) int
dnl                    the size of the
dnl                    vectors x and y is n
dnl                    matrix A is n-by-n.
dnl   
dnl   randomize (in) int
dnl                  if 0, entries in matrices A, x will be chosen for
dnl                        maximum cancellation, but with less randomness.
dnl                  if 1, every entry in the matrix A, x will be 
dnl                        random.
dnl
dnl   alpha     (in/out) SCALAR
dnl                    if alpha_flag = 1, alpha is input
dnl                    if alpha_flag = 0, alpha is output
dnl
dnl   alpha_flag (in) int
dnl                    see above
dnl
dnl   beta      (in/out) SCALAR
dnl                    if beta_flag = 1, beta is input
dnl                    if beta_flag = 0, beta is output
dnl
dnl   beta_flag (in) int
dnl                    see above
dnl
dnl   a         (out) matrix a.
dnl
dnl   k         (in) number of super/sub diagonals of matrix a.
dnl
dnl   lda       (in) leading dimensions of matrix a.
dnl
dnl   x, y      (out) vectors, x, y.
dnl   incx, incy (in) strides for vectors x, y.
dnl
dnl   seed      (in/out) int
dnl
dnl   HEAD(r_true)  (out) double *  (these are vectors of size n)
dnl   TAIL(r_true)  (out) double *
dnl               the leading/trailing part of the true in double-double
dnl
include(cblas.m4)dnl
include(test-common.m4)dnl
dnl
dnl
define(`HBMV_NAME', `ifelse(
  `$2&&$3', `$1&&$1', `BLAS_$1hbmv$4', 
  `BLAS_$1hbmv_$2_$3$4')')dnl
dnl
dnl
define(`SBMV_TESTGEN_NAME', `ifelse(`$2&&$3', `$1&&$1', 
  `BLAS_$1sbmv_testgen', `BLAS_$1sbmv_$2_$3_testgen')')dnl
dnl
dnl
dnl  HBMV_TESTGEN
dnl     |
dnl     |-- HBMV_TESTGEN_HEAD
dnl     |      |
dnl     |      |-- HBMV_TESTGEN_NAME
dnl     |      |
dnl     |      |-- HBMV_TESTGEN_PARAMS
dnl     |      
dnl     |-- HBMV_TESTGEN_COMMENT
dnl     |
dnl     |-- HBMV_TESTGEN_BODY
dnl
dnl  Usage:
dnl    HBMV_TESTGEN        ($1, $2, $3)
dnl    HBMV_TESTGEN_HEAD   ($1, $2, $3)
dnl    HBMV_TESTGEN_NAME   ($1, $2, $3)
dnl    HBMV_TESTGEN_PARAMS ($1, $2, $3)
dnl    HBMV_TESTGEN_COMMENT($1, $2, $3)
dnl    HBMV_TESTGEN_BODY   ($1, $2, $3)
dnl
dnl    $1 -- type of alpha, beta, y
dnl    $2 -- type of a
dnl    $3 -- type of x
dnl
define(`HBMV_TESTGEN_NAME', `ifelse(`$2&&$3', `$1&&$1', 
  `BLAS_$1hbmv_testgen', `BLAS_$1hbmv_$2_$3_testgen')')dnl
dnl
dnl
define(`HBMV_TESTGEN_PARAMS',
  `int norm, enum blas_order_type order, enum blas_uplo_type uplo, dnl
   int n,  int randomize, $1_array alpha, int alpha_flag, $1_array beta, dnl
   int beta_flag, $2_array a, int k, int lda, $3_array x, int incx, dnl
   $1_array y, int incy, int *seed, dnl
   double *HEAD(r_true), double *TAIL(r_true)')dnl
dnl
dnl
define(`HBMV_TESTGEN_HEAD', 
  `void HBMV_TESTGEN_NAME($1, $2, $3)(HBMV_TESTGEN_PARAMS($1, $2, $3))')dnl
dnl
dnl
define(`HBMV_TESTGEN_COMMENT', `
/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to HBMV_NAME($1, $2, $3, `'){_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the hermitian matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
 *              random.
 *
 * alpha   (input/output) $1_array
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) $1_array
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) $2_array
 *
 * k       (input) k
 *         number of sub/super diagonals of hermitian matrix A
 * 
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * x       (input/output) $3_array
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) $1_array
 *         generated vector y that will be used as an input to HBMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * HEAD(r_true)  (output) double *
 *         the leading part of the truth in double-double.
 *
 * TAIL(r_true)  (output) double *
 *         the trailing part of the truth in double-double
 *
 */')dnl
dnl
dnl
dnl
dnl
define(`HBMV_TESTGEN_BODY', `{

  /* Strategy:  (applies to the randomize == 0 case) :
         R1 = alpha * A1 * x + beta * y1
         R2 = alpha * A2 * x + beta * y2
     where all the matrices and vectors are real.  Then let R = R1 + i R2, 
     A = A1 + i A2, y = y1 + i y2.  To make A hermitian, A1 is 
     symmetric, and A2 is a skew matrix (trans(A2) = -A2).
  */

  int i, j;
  int yi;
  int aij, ai;
  int a1ij, a1i;
  int xi;
  int incyi, x_starti, y_starti;
  int incaij, incai;
  int inca_real_ij, inca_real_i;
  int incxi;
  int inca_vec, incx_vec, a_veci;
  int n_i;
  int ab;
  int ri, incri;

  DECLARE_VECTOR(a1, REAL_TYPE($2_type))
  DECLARE_VECTOR(a2, REAL_TYPE($2_type))
  DECLARE_VECTOR(y1, REAL_TYPE($1_type))
  DECLARE_VECTOR(y2, REAL_TYPE($1_type))
  DECLARE_VECTOR(x0, REAL_TYPE($3_type))

  DECLARE_VECTOR(r1_true, EXTRA_TYPE(REAL_TYPE($1_type)))
  DECLARE_VECTOR(r2_true, EXTRA_TYPE(REAL_TYPE($1_type)))

  DECLARE(r_elem1, real_E)
  DECLARE(r_elem2, real_E)
  DECLARE(r_elem, real_E)

  DECLARE_VECTOR(a_vec, $2_type)
  DECLARE_VECTOR(x_vec, $3_type)

  PTR_CAST(y, $1_type)
  PTR_CAST(alpha, $1_type)
  PTR_CAST(beta, $1_type)

  PTR_CAST(a, $2_type)
  PTR_CAST(x, $3_type)

  n_i = n;

  inca_real_i = lda; /* distance between rows (columns for colmajor)
                        in real matricies a1, a2*/
  inca_real_ij = 1;  /* distance between colums (rows for colmajor)
                        in real matricies a1, a2 */
  incai = lda;
  incaij = 1;
  INC_ADJUST(incai, $2_type)
  INC_ADJUST(incaij, $2_type)

  incyi = incy;
  incxi = incx;

  if ((0 == incx) || (0 == incy)) {
        BLAS_error(routine_name, 0, 0, NULL);
  }
  if (incx > 0) {
        x_starti = 0;
  } else {
        x_starti = (-n_i + 1) * incx;
  }
  if (incy > 0) {
        y_starti = 0;
  } else {
        y_starti = (-n_i + 1) * incy;
  }
  incri = 1;  
  INC_ADJUST(incri, EXTRA_TYPE($1_type))
  INC_ADJUST(x_starti, $3_type)
  INC_ADJUST(y_starti, $1_type)
  INC_ADJUST(incyi, $1_type)
  INC_ADJUST(incxi, $3_type)

  inca_vec = incx_vec = 1;
  INC_ADJUST(inca_vec, $2_type)
  INC_ADJUST(incx_vec, $3_type)
  MALLOC_VECTOR(a_vec, $2_type, n_i)
  for (i = 0; i < n_i*inca_vec; i += inca_vec) {
        SET_ZERO_VECTOR_ELEMENT(a_vec, i, $2_type)
  }
  MALLOC_VECTOR(x_vec, $3_type, n_i)
  for (i = 0; i < n_i*incx_vec; i += incx_vec) {
        SET_ZERO_VECTOR_ELEMENT(x_vec, i, $3_type)
  }

  if (randomize == 0) {
      int mi, incmi = 1;
      int incx0, incy1, incy2;
      MALLOC_VECTOR(a1, REAL_TYPE($2_type), n_i * lda)
      MALLOC_VECTOR(a2, REAL_TYPE($2_type), n_i * lda)
      for (i = 0; i < n_i * (k+1); ++i) {
        a1[i] = a2[i] = 0.0;
      }
      MALLOC_VECTOR(y1, REAL_TYPE($1_type), n_i)
      incy1 = 1;
      INC_ADJUST(incy1, REAL_TYPE($1_type))
      MALLOC_VECTOR(y2, REAL_TYPE($1_type), n_i)
      incy2 = 1;
      INC_ADJUST(incy2, REAL_TYPE($1_type))
      MALLOC_VECTOR(x0, REAL_TYPE($3_type), n_i)
      incx0 = 1;
      INC_ADJUST(incx0, REAL_TYPE($3_type))
      for (i = 0; i < n_i; ++i) {
        y1[i] = y2[i] = 0.0;
        x0[i] = 0.0;
      }
      MALLOC_VECTOR(r1_true, EXTRA_TYPE(REAL_TYPE($1_type)), n_i);
      MALLOC_VECTOR(r2_true, EXTRA_TYPE(REAL_TYPE($1_type)), n_i);

      /* First generate the real portion of matrix A, and matrix B.
         Note that Re(A) is a symmetric banded matrix. */
        /* a1, x0 is output from this call */
      SBMV_TESTGEN_NAME(REAL_ABBREV($1), REAL_ABBREV($2), REAL_ABBREV($3))
        (norm, order, uplo, n_i, 0/*randomize*/, 
          alpha_i, alpha_flag, 
          beta_i, beta_flag, 
          a1, k, lda, 
          x0, incx0, y1, incy1, 
          seed, HEAD(r1_true), TAIL(r1_true));

        /* x0 is now fixed and is input to this call */
      SKEW_TESTGEN_HBMV_NAME(REAL_ABBREV($1), REAL_ABBREV($2), REAL_ABBREV($3))
        (norm, order, uplo, n_i, 0, alpha_i, 
                beta_i, a2, k, lda, x0,
                incx0, y2, incy2, seed, HEAD(r2_true), TAIL(r2_true));


      ifelse($3_type, `real_S', `HBMV_TESTGEN_REAL_ADJ($1, $2, $3)', 
             $3_type, `real_D', `HBMV_TESTGEN_REAL_ADJ($1, $2, $3)', 
             $3_type, `complex_S', `HBMV_TESTGEN_COMPLEX_ADJ($1, $2, $3)', 
             $3_type, `complex_D', `HBMV_TESTGEN_COMPLEX_ADJ($1, $2, $3)')

      FREE_VECTOR(a1, REAL_TYPE($2_type))
      FREE_VECTOR(a2, REAL_TYPE($2_type))
      FREE_VECTOR(y1, REAL_TYPE($1_type))
      FREE_VECTOR(y2, REAL_TYPE($1_type))
      FREE_VECTOR(x0, REAL_TYPE($3_type))
      FREE_VECTOR(r1_true, EXTRA_TYPE(REAL_TYPE($1_type)))
      FREE_VECTOR(r2_true, EXTRA_TYPE(REAL_TYPE($1_type)))
  } else {
        /* get random A, x, then compute y for some cancellation. */    
    DECLARE(a_elem, $2_type)
    DECLARE(x_elem, $3_type)
    DECLARE(y_elem, $1_type)
    DECLARE(r_true_elem, EXTRA_TYPE($1_type))

    /* Since mixed real/complex test generator for dot
       scales the vectors, we need to used the non-mixed
       version if x is real (since A is always complex). */
 
    ifelse(`$3', `s', `DECLARE_VECTOR(xx_vec, complex_S)
                       MALLOC_VECTOR(xx_vec, complex_S, n_i)', 
           `$3', `d', `DECLARE_VECTOR(xx_vec, complex_D)
                       MALLOC_VECTOR(xx_vec, complex_D, n_i)')


    if (alpha_flag == 0) {
      RANDOM(alpha_i, $1_type, IS_MIXED_PREC($1_type, $2_type, $3_type))
    }
    if (beta_flag == 0) {
      RANDOM(beta_i, $1_type, IS_MIXED_PREC($1_type, $2_type, $3_type))
    }

    /* Fill in matrix A -- Hermitian. */
    inca_vec = 1;
    INC_ADJUST(inca_vec, $2_type)
    for (i = 0; i < n_i; i++) {
      j = a_veci = MAX(0, i - k);
      INC_ADJUST(a_veci, $2_type)
      for (; j < MIN(n_i, i + k + 1); j++, a_veci += inca_vec ) {
        RANDOM(a_elem, $2_type, IS_MIXED_PREC($1_type, $2_type, $3_type))         
        if (j==i) {
          ZERO_IMAG_PART(a_elem, $2_type)
        }
        SET_VECTOR_ELEMENT(a_vec, a_veci, a_elem, $2_type)
      }
      $2hbmv_commit_row(order, uplo, n_i, a_i, k, lda, a_vec, i);
        
    }

    /* Fill in vector x */
    for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
        RANDOM(x_elem, $3_type, IS_MIXED_PREC($1_type, $2_type, $3_type))
        SET_VECTOR_ELEMENT(x_i, xi, x_elem, $3_type)
    }

    $3copy_vector(x_i, n_i, incx, x_vec, 1);
    IF_REAL($3_type, `
        /* copy the real x_vec into complex xx_vec, so that 
           pure complex test case generator can be called. */
        { int k;
          for (k = 0; k < n_i; k++) {
            xx_vec[2*k] = x_vec[k];
            xx_vec[2*k+1] = 0.0;
          }
        }')

    for (i = 0, yi = y_starti, ri = 0; i < n_i; i++, yi += incyi, 
                        ri += incri) {
      $2hbmv_copy_row(order, uplo,  n_i, a, k, lda, a_vec, i);
        
        DOT_TESTGEN_NAME($1, $2, COMPLEX_ABBREV($3))(n_i, n_i, 0, 
                        norm, blas_no_conj, 
                     alpha_i, 1, beta_i, 1, a_vec, 
                     IF_REAL($3_type, `xx_vec', `x_vec'), seed, 
                     PASS_BY_REF(y_elem, $1_type), 
                     PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));

        SET_VECTOR_ELEMENT(y_i, yi, y_elem, $1_type)
        SET_VECTOR_ELEMENT(r_true, ri, r_true_elem, EXTRA_TYPE($1_type))
    }    

    IF_REAL($3_type, `FREE_VECTOR(xx_vec, COMPLEX_TYPE($3_type));')
  }

  FREE_VECTOR(a_vec, $2_type);
  FREE_VECTOR(x_vec, $3_type);

}')dnl
dnl
dnl
define(`HBMV_TESTGEN_REAL_ADJ', 
`  
  /* The case where x is a real vector. 

     There are four cases to consider, depending on the 
     values of alpha and beta.

         values                             scaling
      alpha  beta         alpha    A    x    beta    y     R (truth)
   0    1      1            
   1    1      ?                              -i     i     
   2    ?      1            i                        i     i
   3    ?      ?           1+i                1+i         1+i

   Note that we can afford to scale truth by (1+i) since they
   are computed in double-double.
  */

  if (alpha_i[0] == 1.0 && beta_i[0] == 1.0) {
    ab = 0;
    alpha_i[1] = beta_i[1] = 0.0;  /* set alpha, beta to be 1. */
  } else if (alpha_i[0] == 1.0) {
    ab = 1;
    /* set alpha to 1, multiply beta by -i. */
    alpha_i[1] = 0.0;
    beta_i[1] = - beta_i[0];
    beta_i[0] = 0.0;        
  } else if (beta_i[0] == 1.0) {
    ab = 2;
    /* set beta to 1, multiply alpha by i. */
    beta_i[1] = 0.0;
    alpha_i[1] = alpha_i[0];
    alpha_i[0] = 0.0;    
  } else {
    ab = 3;
    /* multiply alpha, beta by (1 + i). */
    alpha_i[1] = alpha_i[0];
    beta_i[1] = beta_i[0];
  }

  /* Now fill in a */
  for (i = 0, ai = 0, a1i = 0; i < n_i; i++, ai += incai,
                a1i += inca_real_i) {
    for (j = 0, aij = ai, a1ij = a1i; j < lda; 
          j++, aij += incaij, a1ij += inca_real_ij) {
      a_i[aij] = a1[a1ij];
      a_i[aij+1] = a2[a1ij];
    }
  }     

  /* Fill in x */
  for (i = 0, xi = x_starti, mi = 0; i < n_i; i++, xi += incxi, 
                mi += incx0) {
      x_i[xi] = x0[mi];    
  }

  /* Fill in y */
  for (i = 0, yi = y_starti, mi = 0; i < n_i; i++, yi += incyi, 
        mi += incy1) {
      if (ab == 1 || ab == 2) {
        y_i[yi] = -y2[mi];
        y_i[yi+1] = y1[mi];
      } else {
        y_i[yi] = y1[mi];
        y_i[yi+1] = y2[mi];
      }
  }

  /* Fill in truth */
  for (i = 0, ri = 0, mi = 0; i < n_i; i++, ri += incri, 
                mi += incmi) {
      if (ab == 0 || ab == 1) {
        HEAD(r_true)[ri] = HEAD(r1_true)[mi];
        TAIL(r_true)[ri] = TAIL(r1_true)[mi];
        HEAD(r_true)[ri+1] = HEAD(r2_true)[mi];
        TAIL(r_true)[ri+1] = TAIL(r2_true)[mi];
      } else if (ab == 2) {
        HEAD(r_true)[ri] = -HEAD(r2_true)[mi];
        TAIL(r_true)[ri] = -TAIL(r2_true)[mi];
        HEAD(r_true)[ri+1] = HEAD(r1_true)[mi];
        TAIL(r_true)[ri+1] = TAIL(r1_true)[mi];
      } else {
        HEAD(r_elem1) = HEAD(r1_true)[mi];
        TAIL(r_elem1) = TAIL(r1_true)[mi];

        HEAD(r_elem2) = HEAD(r2_true)[mi];
        TAIL(r_elem2) = TAIL(r2_true)[mi];

        ADD(r_elem, real_E, r_elem1, real_E, r_elem2, real_E)

        /* Set the imaginary part to  R1 + R2 */
        TAIL(r_true)[ri+1] = TAIL(r_elem);
        HEAD(r_true)[ri+1] = HEAD(r_elem);

        /* Set the real part to R1 - R2. */
        SUB(r_elem, real_E, r_elem1, real_E, r_elem2, real_E)
        TAIL(r_true)[ri] = TAIL(r_elem);
        HEAD(r_true)[ri] = HEAD(r_elem);
      }
  }
')dnl
dnl
dnl
define(`HBMV_TESTGEN_COMPLEX_ADJ', 
` 
  /* The case where x is a complex vector.  Since x is generated
     as a real vector, we need to perform some scaling.

     There are four cases to consider, depending on the values
     of alpha and beta.

            values                         scaling
         alpha   beta      alpha  A    x       beta    y    R (truth)
      0    1      1                    i               i    i
      1    1      ?                   1+i      1+i         1+i
      2    ?      1         1+i       1+i             2i    2i
      3    ?      ?         1+i       1+i      2i           2i

     Note that we can afford to scale R by 1+i, since they are
     computed in double-double precision.
   */

  if (alpha_i[0] == 1.0 && beta_i[0] == 1.0) {
    ab = 0;
    alpha_i[1] = beta_i[1] = 0.0;  /* set alpha, beta to be 1. */
  } else if (alpha_i[0] == 1.0) {
    ab = 1;
    /* set alpha to 1, multiply beta by 1+i. */
    alpha_i[1] = 0.0;
    beta_i[1] = beta_i[0];
  } else if (beta_i[0] == 1.0) {
    ab = 2;
    /* set beta to 1, multiply alpha by 1+i. */
    beta_i[1] = 0.0;
    alpha_i[1] = alpha_i[0];
  } else {
    ab = 3;
    /* multiply alpha by 1+i, beta by 2i. */
    alpha_i[1] = alpha_i[0];
    beta_i[1] = 2.0 * beta_i[0];
    beta_i[0] = 0.0;
  }

  /* Now fill in a */
  for (i = 0, ai = 0, a1i = 0; i < n_i; i++, ai += incai, 
                a1i += inca_real_i) {
    for (j = 0, aij = ai, a1ij = a1i; j < lda; 
          j++, aij += incaij, a1ij += inca_real_ij) {
      a_i[aij] = a1[a1ij];
      a_i[aij+1] = a2[a1ij];
    }
  }

  /* Fill in x */
  for (i = 0, xi = x_starti, mi = 0; i < n_i; i++, xi += incxi, 
                mi += incx0) {
      if (ab == 0) {
        x_i[xi] = 0.0;
        x_i[xi+1] = x0[mi];
      } else {
        x_i[xi] = x0[mi];
        x_i[xi+1] = x0[mi];
      }
  }

  /* Fill in y */
  for (i = 0, yi = y_starti, mi = 0; i < n_i; 
        i++, yi += incyi, mi += incy1) {
      if (ab == 0) {
        y_i[yi] = -y2[mi];
        y_i[yi+1] = y1[mi];
      } else if (ab == 2) {
        y_i[yi] = -2.0 * y2[mi];
        y_i[yi+1] = 2.0 * y1[mi];
      } else {
        y_i[yi] = y1[mi];
        y_i[yi+1] = y2[mi];
      }
  }

  /* Fill in the truth */
  for (i = 0, ri = 0, mi = 0; i < n_i; i++, ri += incri, 
        mi += incmi) {

        HEAD(r_elem1) = HEAD(r1_true)[mi];
        TAIL(r_elem1) = TAIL(r1_true)[mi];

        HEAD(r_elem2) = HEAD(r2_true)[mi];
        TAIL(r_elem2) = TAIL(r2_true)[mi];

      if (ab == 0) {
        HEAD(r_true)[ri] = -HEAD(r_elem2);
        TAIL(r_true)[ri] = -TAIL(r_elem2);
        HEAD(r_true)[ri+1] = HEAD(r_elem1);
        TAIL(r_true)[ri+1] = TAIL(r_elem1);
      } else if (ab == 1) {
        ADD(r_elem, real_E, r_elem1, real_E, r_elem2, real_E)

        /* Set the imaginary part to  R1 + R2 */
        TAIL(r_true)[ri+1] = TAIL(r_elem);
        HEAD(r_true)[ri+1] = HEAD(r_elem);

        /* Set the real part to R1 - R2. */
        SUB(r_elem, real_E, r_elem1, real_E, r_elem2, real_E)
        TAIL(r_true)[ri] = TAIL(r_elem);
        HEAD(r_true)[ri] = HEAD(r_elem);
      } else {

        /* Real part */
        MUL(r_elem, real_E, r_elem2, real_E, -2.0, real_D)
        HEAD(r_true)[ri] = HEAD(r_elem);
        TAIL(r_true)[ri] = TAIL(r_elem);

        /* Imaginary Part */
        MUL(r_elem, real_E, r_elem1, real_E, 2.0, real_D)
        HEAD(r_true)[ri+1] = HEAD(r_elem);
        TAIL(r_true)[ri+1] = TAIL(r_elem);
      }
  }
')dnl
dnl
dnl
dnl
dnl
define(`HBMV_TESTGEN', 
  `HBMV_TESTGEN_HEAD($1, $2, $3)
   HBMV_TESTGEN_COMMENT($1, $2, $3)
{
   char *routine_name = "HBMV_TESTGEN_NAME($1, $2, $3)";
   HBMV_TESTGEN_BODY($1, $2, $3)
} /* end HBMV_TESTGEN_NAME($1, $2, $3) */')dnl
dnl
dnl
dnl
dnl ------------------------------ Skew Matrix Generator:
dnl
dnl
dnl  SKEW_HBMV_TESTGEN
dnl     |
dnl     |-- SKEW_HBMV_TESTGEN_HEAD
dnl     |      |
dnl     |      |-- SKEW_HBMV_TESTGEN_NAME
dnl     |      |
dnl     |      |-- SKEW_HBMV_TESTGEN_PARAMS
dnl     |      
dnl     |-- SKEW_HBMV_TESTGEN_COMMENT
dnl     |
dnl     |-- SKEW_HBMV_TESTGEN_BODY
dnl
dnl  Usage:
dnl    SKEW_HBMV_TESTGEN        ($1, $2, $3)
dnl    SKEW_HBMV_TESTGEN_HEAD   ($1, $2, $3)
dnl    SKEW_HBMV_TESTGEN_NAME   ($1, $2, $3)
dnl    SKEW_HBMV_TESTGEN_PARAMS ($1, $2, $3)
dnl    SKEW_HBMV_TESTGEN_COMMENT($1, $2, $3)
dnl    SKEW_HBMV_TESTGEN_BODY   ($1, $2, $3)
dnl
dnl    $1 -- type of alpha, beta, y
dnl    $2 -- type of a
dnl    $3 -- type of x
dnl
define(`SKEW_TESTGEN_HBMV_NAME', `ifelse(`$2&&$3', `$1&&$1', 
  `BLAS_$1skew_testgen_hbmv', `BLAS_$1skew_testgen_hbmv_$2_$3')')dnl
dnl
dnl
define(`SKEW_TESTGEN_HBMV_PARAMS',
  `int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n,  int randomize, 
   $1_array alpha, $1_array beta, 
   $2_array a, int k, int lda, $3_array x, int incx, $1_array y, int incy, 
   int *seed, double *HEAD(r_true), double *TAIL(r_true)')dnl
dnl
dnl
define(`SKEW_TESTGEN_HBMV_HEAD', 
  `void SKEW_TESTGEN_HBMV_NAME($1, $2, $3)(SKEW_TESTGEN_HBMV_PARAMS($1, $2, $3))')dnl
dnl
dnl
define(`SKEW_TESTGEN_HBMV_COMMENT', `
/*
 * Purpose
 * =======
 *
 *   Generates a skew Matrix for use with hbmv testing
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the skew symmetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
 *              random.
 *
 * alpha   (input/output) $1_array
 *
 * beta    (input) $1_array
 *
 * a       (input/output) $2_array
 *
 * k       (input) k
 *         Number of sub/super diagonals in A
 * 
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * x       (input) $3_array
 *         note : x is input only. x should be determined before calling
 *         this function.
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) $1_array
 *         generated vector y.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * HEAD(r_true)  (output) double *
 *         the leading part of the truth in double-double.
 *
 * TAIL(r_true)  (output) double *
 *         the trailing part of the truth in double-double
 *
 */')dnl
dnl
dnl
define(`SKEW_TESTGEN_HBMV_BODY', `{

  int i;
  int yi;
  int ri;
  int incyi, incri;
  int incx_veci, y_starti;
  int inca_vec;
  int n_i;

  DECLARE(y_elem, $1_type)
  DECLARE(r_true_elem, EXTRA_TYPE($1_type))

  DECLARE_VECTOR(a_vec, $2_type)
  DECLARE_VECTOR(x_vec, $3_type)

  PTR_CAST(y, $1_type)
  PTR_CAST(alpha, $1_type)
  PTR_CAST(beta, $1_type)
  PTR_CAST(a, $2_type)
  PTR_CAST(x, $3_type)

  n_i = n;

  /*a_vec must have stride of 1*/
  inca_vec = 1;
  INC_ADJUST(inca_vec, $2_type)

  MALLOC_VECTOR(a_vec, $2_type, n_i)
  for (i = 0; i < n_i; i += inca_vec) {
        SET_ZERO_VECTOR_ELEMENT(a_vec, i, $2_type)
  }
  MALLOC_VECTOR(x_vec, $3_type, 2*n_i)
  $3copy_vector(x_i, n_i, incx, x_vec, 1);

  incyi = incy;
  INC_ADJUST(incyi, $1_type)
  if (incyi < 0) {
        y_starti = (-n+1) * incyi;
  }
  else {
        y_starti = 0;
  }

  incri = 1;
  INC_ADJUST(incri, EXTRA_TYPE($1_type))

  incx_veci = 1;
  INC_ADJUST(incx_veci, $3_type)

  if (randomize == 0) {
        int n_on_row_to_do, both_fixed, mixed_fixed;

    /* Fill in skew matrix A */
    for(i = 0, yi = y_starti, ri = 0; i < n_i; i++, ri += incri, 
                yi += incyi) {
        /* x_i has already been copied to x_vec */
      $2skew_copy_row_hbmv(order, uplo, n_i, a, k, lda, a_vec, i);

        /* skew matricies have zeroed diagonals */
      SET_ZERO_VECTOR_ELEMENT(a_vec, i, $2_type)
        
        n_on_row_to_do = MIN(n_i, i+k+1);
        both_fixed = MIN(n_i, i+1);
        mixed_fixed = n_on_row_to_do - both_fixed;

      DOT_TESTGEN_NAME($1, $3, $2)(n_on_row_to_do, both_fixed, mixed_fixed, 
                                   norm, 
                                   blas_no_conj, alpha_i, 1, 
                                   beta_i, 1, x_vec, a_vec, seed, 
                                   PASS_BY_REF(y_elem, $1_type), 
                                   PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));

      $2skew_commit_row_hbmv(order, uplo, n_i, a_i, k, lda, a_vec, i);
        

        /*commits an element to the generated y */
      SET_VECTOR_ELEMENT(y_i, yi, y_elem, $1_type)
      SET_VECTOR_ELEMENT(r_true, ri, r_true_elem, EXTRA_TYPE($1_type))
    }

  } else {
        printf("\nError- should not use random case for skew generator\n");
        fflush(stdout);
  }

  FREE_VECTOR(a_vec, $2_type)
  FREE_VECTOR(x_vec, $3_type)
}')dnl
dnl
dnl
define(`SKEW_TESTGEN_HBMV', 
  `SKEW_TESTGEN_HBMV_HEAD($1, $2, $3)
   SKEW_TESTGEN_HBMV_COMMENT($1, $2, $3)
   SKEW_TESTGEN_HBMV_BODY($1, $2, $3)')dnl
dnl
dnl
define(`PROTOTYPES', `dnl
SKEW_TESTGEN_HBMV_HEAD(s, s, s);
SKEW_TESTGEN_HBMV_HEAD(d, d, d);
SKEW_TESTGEN_HBMV_HEAD(d, d, s);
SKEW_TESTGEN_HBMV_HEAD(d, s, d);
SKEW_TESTGEN_HBMV_HEAD(d, s, s);

HBMV_TESTGEN_HEAD(c, c, c);
HBMV_TESTGEN_HEAD(z, z, z);
HBMV_TESTGEN_HEAD(z, c, z);
HBMV_TESTGEN_HEAD(z, z, c);
HBMV_TESTGEN_HEAD(z, c, c);
HBMV_TESTGEN_HEAD(z, z, d);
HBMV_TESTGEN_HEAD(c, c, s);
')dnl
dnl
dnl
define(`SOURCE', `
#include <stdio.h>
#include <stdlib.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"

SKEW_TESTGEN_HBMV(s, s, s)
SKEW_TESTGEN_HBMV(d, d, d)
SKEW_TESTGEN_HBMV(d, d, s)
SKEW_TESTGEN_HBMV(d, s, d)
SKEW_TESTGEN_HBMV(d, s, s)

HBMV_TESTGEN(c, c, c)
HBMV_TESTGEN(z, z, z)
HBMV_TESTGEN(z, c, z)
HBMV_TESTGEN(z, z, c)
HBMV_TESTGEN(z, c, c)
HBMV_TESTGEN(z, z, d)
HBMV_TESTGEN(c, c, s)
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
