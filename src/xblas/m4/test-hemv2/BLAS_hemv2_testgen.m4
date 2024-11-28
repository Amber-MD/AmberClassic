dnl
dnl
include(cblas.m4)dnl
include(test-common.m4)dnl
dnl
dnl
define(`HEMV2_NAME', `ifelse(`$2&&$3', `$1&&$1',
  `BLAS_$1hemv2$4', `BLAS_$1hemv_$2_$3$4')')dnl
dnl
dnl
define(`SYMV2_TESTGEN_NAME', `ifelse(`$2&&$3', `$1&&$1', 
  `BLAS_$1symv2_testgen', `BLAS_$1symv2_$2_$3_testgen')')dnl
dnl
dnl
dnl  HEMV2_TESTGEN
dnl     |
dnl     +-- HEMV2_TESTGEN_HEAD
dnl     |      |
dnl     |      +-- HEMV2_TESTGEN_NAME
dnl     |      |
dnl     |      +-- HEMV2_TESTGEN_PARAMS
dnl     |      
dnl     +-- HEMV2_TESTGEN_COMMENT
dnl     |
dnl     +-- HEMV2_TESTGEN_BODY
dnl
dnl  Usage:
dnl    HEMV2_TESTGEN        ($1, $2, $3)
dnl    HEMV2_TESTGEN_HEAD   ($1, $2, $3)
dnl    HEMV2_TESTGEN_NAME   ($1, $2, $3)
dnl    HEMV2_TESTGEN_PARAMS ($1, $2, $3)
dnl    HEMV2_TESTGEN_COMMENT($1, $2, $3)
dnl    HEMV2_TESTGEN_BODY   ($1, $2, $3)
dnl
dnl    $1 -- type of alpha, beta, y
dnl    $2 -- type of a
dnl    $3 -- type of x
dnl
define(`HEMV2_TESTGEN_NAME', `ifelse(`$2&&$3', `$1&&$1', 
  `BLAS_$1hemv2_testgen', `BLAS_$1hemv2_$2_$3_testgen')')dnl
dnl
dnl
define(`HEMV2_TESTGEN_PARAMS', `int norm, enum blas_order_type order, dnl
   enum blas_uplo_type uplo, int n, $1_array alpha, dnl
   int alpha_flag, $1_array beta, int beta_flag, $2_array a, int lda, dnl
   $3_array x_head, $3_array x_tail, $1_array y, dnl
   int *seed, double *HEAD(r_true), double *TAIL(r_true)')dnl
dnl
dnl
define(`HEMV2_TESTGEN_HEAD', 
  `void HEMV2_TESTGEN_NAME($1, $2, $3)(HEMV2_TESTGEN_PARAMS($1, $2, $3))')dnl
dnl
dnl
define(`HEMV2_TESTGEN_COMMENT', `
/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to HEMV2_NAME($1, $2, $3){_x}
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
 *           which half of the hemv2etric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
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
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * x_head  (input/output) $3_array
 * x_tail  (input/output) $3_array
 *         vector x = (x_head + x_tail)
 *
 * y       (input/output) $1_array
 *         generated vector y that will be used as an input to HEMV2.
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
 */')
dnl
dnl
define(`HEMV2_TESTGEN_BODY', `{

  /* Strategy:  
         r1 = alpha * A1 * x + beta * y1
         r2 = alpha * A2 * x + beta * y2
     where all the matrices and vectors are real.  Then let R = R1 + i R2, 
     A = A1 + i A2, y = y1 + i y2.  To make A Hermitian, A1 is symmetric,
     and A2 is skew-symmetric.
  */

  int i, j;
  int yi;
  int aij, ai;
  int a1ij, a1i;
  int xi;
  int mi;
  int incx = 1, incy = 1;
  int incyi, xi0, yi0;
  int incaij, incai;
  int inca1ij, inca1i;
  int incxi;
  int inca_vec, incx_vec;
  int ld;
  int ab;
  int ri, incri;

  DECLARE_VECTOR(a1, REAL_TYPE($2_type))
  DECLARE_VECTOR(a2, REAL_TYPE($2_type))
  DECLARE_VECTOR(y1, REAL_TYPE($1_type))
  DECLARE_VECTOR(y2, REAL_TYPE($1_type))
  DECLARE_VECTOR(x_head_0, REAL_TYPE($3_type))
  DECLARE_VECTOR(x_tail_0, REAL_TYPE($3_type))

  DECLARE_VECTOR(r1_true, EXTRA_TYPE(REAL_TYPE($1_type)))
  DECLARE_VECTOR(r2_true, EXTRA_TYPE(REAL_TYPE($1_type)))

  DECLARE(r_elem1, real_E)
  DECLARE(r_elem2, real_E)
  DECLARE(r_elem, real_E)

  DECLARE_VECTOR(a_vec, $2_type)
  DECLARE_VECTOR(x_head_vec, $3_type)
  DECLARE_VECTOR(x_tail_vec, $3_type)

  PTR_CAST(y, $1_type)
  PTR_CAST(alpha, $1_type)
  PTR_CAST(beta, $1_type)

  PTR_CAST(a, $2_type)
  PTR_CAST(x_head, $3_type)
  PTR_CAST(x_tail, $3_type)

  ld = n;
  if (order == blas_colmajor) {
    inca1i = incai = 1;
    incyi = incy;
    incaij = lda;
    incxi = incx;
    inca1ij = n;
  } else {
    incyi = incy;
    incai = lda;
    incxi = incx;
    inca1i = n;
    inca1ij = incaij = 1;
  }
  xi0 = (incx > 0) ? 0 : -(n-1)*incx;
  yi0 = (incy > 0) ? 0 : -(n-1)*incy;
  incri = 1;  
  INC_ADJUST(incri, EXTRA_TYPE($1_type))
  INC_ADJUST(xi0, $3_type)
  INC_ADJUST(yi0, $1_type)
  INC_ADJUST(incyi, $1_type)
  INC_ADJUST(incai, $2_type)
  INC_ADJUST(incaij, $2_type)
  INC_ADJUST(incxi, $3_type)

  inca_vec = incx_vec = 1;
  INC_ADJUST(inca_vec, $2_type)
  INC_ADJUST(incx_vec, $3_type)
  MALLOC_VECTOR(a_vec, $2_type, n)
  MALLOC_VECTOR(x_head_vec, $3_type, n)
  MALLOC_VECTOR(x_tail_vec, $3_type, n)

  int incx0, incy1, incy2, incmi = 1;
  MALLOC_VECTOR(a1, REAL_TYPE($2_type), n * n)
  MALLOC_VECTOR(a2, REAL_TYPE($2_type), n * n)
  MALLOC_VECTOR(y1, REAL_TYPE($1_type), n)
  incy1 = 1;
  INC_ADJUST(incy1, REAL_TYPE($1_type))
  MALLOC_VECTOR(y2, REAL_TYPE($1_type), n)
  incy2 = 1;
  INC_ADJUST(incy2, REAL_TYPE($1_type))
  MALLOC_VECTOR(x_head_0, REAL_TYPE($3_type), n)
  MALLOC_VECTOR(x_tail_0, REAL_TYPE($3_type), n)
  incx0 = 1;
  INC_ADJUST(incx0, REAL_TYPE($3_type))
  MALLOC_VECTOR(r1_true, EXTRA_TYPE(REAL_TYPE($1_type)), n);
  MALLOC_VECTOR(r2_true, EXTRA_TYPE(REAL_TYPE($1_type)), n);

  /* First generate the real portion of A and x. */
  SYMV2_TESTGEN_NAME(REAL_ABBREV($1), REAL_ABBREV($2), REAL_ABBREV($3))( dnl
      norm, order, uplo, n, alpha_i, alpha_flag, a1, n, x_head_0, x_tail_0, dnl
      beta_i, beta_flag, y1, seed, HEAD(r1_true), TAIL(r1_true));

  /* x0 is now fixed and is input to this call */
  SKMV2_TESTGEN_NAME(REAL_ABBREV($1), REAL_ABBREV($2), REAL_ABBREV($3))( dnl
      norm, order, uplo, n, alpha_i, beta_i, a2, n, x_head_0, x_tail_0, dnl
      y2, seed, HEAD(r2_true), TAIL(r2_true));

  IF_REAL($3_type, `HEMV2_TESTGEN_REAL_ADJ($1, $2, $3)',
                  `HEMV2_TESTGEN_COMPLEX_ADJ($1, $2, $3)')

  FREE_VECTOR(a1, REAL_TYPE($2_type))
  FREE_VECTOR(a2, REAL_TYPE($2_type))
  FREE_VECTOR(y1, REAL_TYPE($1_type))
  FREE_VECTOR(y2, REAL_TYPE($1_type))
  FREE_VECTOR(x_head_0, REAL_TYPE($3_type))
  FREE_VECTOR(x_tail_0, REAL_TYPE($3_type))
  FREE_VECTOR(r1_true, EXTRA_TYPE(REAL_TYPE($1_type)))
  FREE_VECTOR(r2_true, EXTRA_TYPE(REAL_TYPE($1_type)))
  FREE_VECTOR(a_vec, $2_type)
  FREE_VECTOR(x_head_vec, $3_type)
  FREE_VECTOR(x_tail_vec, $3_type)
}')dnl
dnl
dnl
define(`HEMV2_TESTGEN_REAL_ADJ', 
`  
  /* The case where x is a real vector. 

     There are four cases to consider, depending on the 
     values of alpha and beta.

         values                             scaling
      alpha  beta         alpha    A    x    beta    y     r (truth)
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
  for (i = 0, ai = 0, a1i = 0; i < n; i++, ai += incai, a1i += inca1i) {
    for (j = 0, aij = ai, a1ij = a1i; j < n; 
          j++, aij += incaij, a1ij += inca1ij) {
      a_i[aij] = a1[a1ij];
      a_i[aij+1] = a2[a1ij];
    }
  }

  /* Fill in x */
  for (i = 0, xi = xi0, mi = 0; i < n; i++, xi += incxi, mi += incx0) {
    x_head_i[xi] = x_head_0[mi];    
    x_tail_i[xi] = x_tail_0[mi];    
  }

  /* Fill in y */
  for (i = 0, yi = yi0, mi = 0; i < n; i++, yi += incyi, mi += incy1) {
    if (ab == 1 || ab == 2) {
      y_i[yi] = -y2[mi];
      y_i[yi+1] = y1[mi];
    } else {
      y_i[yi] = y1[mi];
      y_i[yi+1] = y2[mi];
    }
  }

  /* Fill in truth */
  for (i = 0, ri = 0, mi = 0; i < n; i++, ri += incri, mi += incmi) {
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

      /* Set the imaginary part to  r1 + r2 */
      TAIL(r_true)[ri+1] = TAIL(r_elem);
      HEAD(r_true)[ri+1] = HEAD(r_elem);

      /* Set the real part to r1 - r2. */
      SUB(r_elem, real_E, r_elem1, real_E, r_elem2, real_E)
      TAIL(r_true)[ri] = TAIL(r_elem);
      HEAD(r_true)[ri] = HEAD(r_elem);
    }
  }
')dnl
dnl
dnl
define(`HEMV2_TESTGEN_COMPLEX_ADJ', 
` 
  /* The case where x is a complex vector.  Since x is generated
     as a real vector, we need to perform some scaling.

     There are four cases to consider, depending on the values
     of alpha and beta.

            values                         scaling
         alpha   beta      alpha  A    x       beta    y    r (truth)
      0    1      1                    i               i    i
      1    1      ?                   1+i      1+i         1+i
      2    ?      1         1+i       1+i             2i    2i
      3    ?      ?         1+i       1+i      2i           2i

     Note that we can afford to scale r by 1+i, since they are
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
  for (i = 0, ai = 0, a1i = 0; i < n; i++, ai += incai, a1i += inca1i) {
    for (j = 0, aij = ai, a1ij = a1i; j < n; 
         j++, aij += incaij, a1ij += inca1ij) {
      a_i[aij] = a1[a1ij];
      a_i[aij+1] = a2[a1ij];
    }
  }

  /* Fill in x */
  for (i = 0, xi = xi0, mi = 0; i < n; i++, xi += incxi, mi += incx0) {
    if (ab == 0) {
      x_head_i[xi] = 0.0;
      x_head_i[xi+1] = x_head_0[mi];
      x_tail_i[xi] = 0.0;
      x_tail_i[xi+1] = x_tail_0[mi];
    } else {
      x_head_i[xi] = x_head_0[mi];
      x_head_i[xi+1] = x_head_0[mi];
      x_tail_i[xi] = x_tail_0[mi];
      x_tail_i[xi+1] = x_tail_0[mi];
    }
  }

  /* Fill in y */
  for (i = 0, yi = yi0, mi = 0; i < n; 
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
  for (i = 0, ri = 0, mi = 0; i < n; i++, ri += incri, mi += incmi) {

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
define(`HEMV2_TESTGEN', `
HEMV2_TESTGEN_HEAD($1, $2, $3)
HEMV2_TESTGEN_COMMENT($1, $2, $3)
{
  HEMV2_TESTGEN_BODY($1, $2, $3)
}')dnl
dnl
dnl
dnl
dnl ------------------------------ Skew Matrix Generator:
dnl
dnl
dnl  SKMV2_TESTGEN
dnl     |
dnl     |-- SKMV2_TESTGEN_HEAD
dnl     |      |
dnl     |      |-- SKMV2_TESTGEN_NAME
dnl     |      |
dnl     |      |-- SKMV2_TESTGEN_PARAMS
dnl     |      
dnl     |-- SKMV2_TESTGEN_COMMENT
dnl     |
dnl     |-- SKMV2_TESTGEN_BODY
dnl
dnl  Usage:
dnl    SKMV2_TESTGEN        ($1, $2, $3)
dnl    SKMV2_TESTGEN_HEAD   ($1, $2, $3)
dnl    SKMV2_TESTGEN_NAME   ($1, $2, $3)
dnl    SKMV2_TESTGEN_PARAMS ($1, $2, $3)
dnl    SKMV2_TESTGEN_COMMENT($1, $2, $3)
dnl    SKMV2_TESTGEN_BODY   ($1, $2, $3)
dnl
dnl    $1 -- type of alpha, beta, y
dnl    $2 -- type of a
dnl    $3 -- type of x
dnl
define(`SKMV2_TESTGEN_NAME', `ifelse(`$2&&$3', `$1&&$1', 
  `BLAS_$1skmv2_testgen', `BLAS_$1skmv2_testgen_$2_$3')')dnl
dnl
dnl
define(`SKMV2_TESTGEN_PARAMS', `int norm, enum blas_order_type order, dnl
   enum blas_uplo_type uplo, int n,  $1_array alpha, $1_array beta, dnl
   $2_array a, int lda, $3_array x_head, $3_array x_tail, $1_array y, dnl
   int *seed, double *HEAD(r_true), double *TAIL(r_true)')dnl
dnl
dnl
define(`SKMV2_TESTGEN_HEAD', 
  `void SKMV2_TESTGEN_NAME($1, $2, $3)(SKMV2_TESTGEN_PARAMS($1, $2, $3))')dnl
dnl
dnl
define(`SKMV2_TESTGEN_COMMENT', `
/*
 * Purpose
 * =======
 *
 *   Generates a skew Matrix for use with hemv2 testing
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
 * alpha   (input/output) $1_array
 *
 * beta    (input) $1_array
 *
 * a       (input/output) $2_array
 * 
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * x_head  (input) $3_array
 * x_tail  (input) $3_array
 *         note : x is input only. x should be determined before calling
 *         this function.
 *
 * y       (input/output) $1_array
 *         generated vector y.
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
define(`SKMV2_TESTGEN_BODY', `{

  int i, yi, ri;
  int incx = 1, incy = 1;
  int incyi, incri;
  int incx_veci, yi0;
  int inca_vec;

  DECLARE(y_elem, $1_type)
  DECLARE(r_true_elem, EXTRA_TYPE($1_type))

  DECLARE_VECTOR(a_vec, $2_type)
  DECLARE_VECTOR(x_head_vec, $3_type)
  DECLARE_VECTOR(x_tail_vec, $3_type)

  PTR_CAST(y, $1_type)
  PTR_CAST(alpha, $1_type)
  PTR_CAST(beta, $1_type)
  PTR_CAST(x_head, $3_type)
  PTR_CAST(x_tail, $3_type)

  /*a_vec must have stride of 1*/
  inca_vec = 1;
  INC_ADJUST(inca_vec, $2_type)

  MALLOC_VECTOR(a_vec, $2_type, n)
  for (i = 0; i < n; i += inca_vec) {
        SET_ZERO_VECTOR_ELEMENT(a_vec, i, $2_type)
  }
  MALLOC_VECTOR(x_head_vec, $3_type, n)
  MALLOC_VECTOR(x_tail_vec, $3_type, n)
  $3copy_vector(x_head_i, n, incx, x_head_vec, 1);
  $3copy_vector(x_tail_i, n, incx, x_tail_vec, 1);

  incyi = incy;
  INC_ADJUST(incyi, $1_type)
  yi0 = (incy > 0) ? 0 : -(n-1)*incyi;

  incri = 1;
  INC_ADJUST(incri, EXTRA_TYPE($1_type))

  incx_veci = 1;
  INC_ADJUST(incx_veci, $3_type)

  /* Fill in skew matrix A */
  for(i = 0, yi = yi0, ri = 0; i < n; i++, ri += incri, yi += incyi) {
    /* x_i has already been copied to x_vec */
    $2skew_copy_row(order, uplo, blas_left_side, n, a, lda, a_vec, i);

    /* skew matricies have zeroed diagonals */
    SET_ZERO_VECTOR_ELEMENT(a_vec, i, $2_type)
    DOT2_TESTGEN_NAME($1, $3, $2)(n, i+1, n-i-1, norm, dnl
                                  blas_no_conj, alpha_i, 1, dnl
                                  beta_i, 1, x_head_vec, x_tail_vec, a_vec, dnl
                                  seed, PASS_BY_REF(y_elem, $1_type), dnl
                                  PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));

    $2skew_commit_row(order, uplo, blas_left_side, n, a, lda, a_vec, i);

    /*commits an element to the generated y */
    SET_VECTOR_ELEMENT(y_i, yi, y_elem, $1_type)
    SET_VECTOR_ELEMENT(r_true, ri, r_true_elem, EXTRA_TYPE($1_type))
  }

  FREE_VECTOR(a_vec, $2_type)
  FREE_VECTOR(x_head_vec, $3_type)
  FREE_VECTOR(x_tail_vec, $3_type)
}')dnl
dnl
dnl
define(`SKMV2_TESTGEN', 
  `SKMV2_TESTGEN_HEAD($1, $2, $3)
   SKMV2_TESTGEN_COMMENT($1, $2, $3)
   SKMV2_TESTGEN_BODY($1, $2, $3)')dnl
dnl
dnl
define(`PROTOTYPES', `dnl
SKMV2_TESTGEN_HEAD(s, s, s);
SKMV2_TESTGEN_HEAD(d, d, d);
SKMV2_TESTGEN_HEAD(d, d, s);
SKMV2_TESTGEN_HEAD(d, s, d);
SKMV2_TESTGEN_HEAD(d, s, s);

HEMV2_TESTGEN_HEAD(c, c, c);
HEMV2_TESTGEN_HEAD(z, z, z);
HEMV2_TESTGEN_HEAD(z, c, z);
HEMV2_TESTGEN_HEAD(z, z, c);
HEMV2_TESTGEN_HEAD(z, c, c);
HEMV2_TESTGEN_HEAD(z, z, d);
HEMV2_TESTGEN_HEAD(c, c, s);
')dnl
dnl
dnl
define(`SOURCE', `dnl
#include <stdio.h>
#include <stdlib.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"

SKMV2_TESTGEN(s, s, s)
SKMV2_TESTGEN(d, d, d)
SKMV2_TESTGEN(d, d, s)
SKMV2_TESTGEN(d, s, d)
SKMV2_TESTGEN(d, s, s)

HEMV2_TESTGEN(c, c, c)
HEMV2_TESTGEN(z, z, z)
HEMV2_TESTGEN(z, c, z)
HEMV2_TESTGEN(z, z, c)
HEMV2_TESTGEN(z, c, c)
HEMV2_TESTGEN(z, z, d)
HEMV2_TESTGEN(c, c, s)
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
