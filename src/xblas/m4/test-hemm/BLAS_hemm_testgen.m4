dnl
dnl hemm_testgen.m4
dnl
dnl Test case generator for hemm routines.
dnl Generates test cases for alpha, A, beta, B and C and
dnl computes r_true in double-double precision.
dnl
dnl
dnl C Interface has the form 
dnl
dnl   BLAS_xhemm{_a_b}{_x}_testgen(int norm, 
dnl                             enum blas_order_type order, 
dnl                             enum blas_uplo_type uplo, 
dnl                             enum blas_side_type side, 
dnl                             int m, int n, int randomize, 
dnl                             SCALAR *alpha, int alpha_flag, 
dnl                             SCALAR *beta, int beta_flag, 
dnl                             ARRAY a, int lda, 
dnl                             ARRAY b, int ldb, 
dnl                             ARRAY c, int ldc, 
dnl                             int *seed, 
dnl                             double *HEAD(r_true), double *TAIL(r_true))
dnl
dnl Arguments
dnl   norm      (in) int
dnl
dnl   order     (in) blas_order_type
dnl                    determines the storage format for matrices
dnl   
dnl   uplo      (in) blas_uplo_type
dnl                    determines whether the upper triangular portion
dnl                    or the lower triangular portion of the hermitian
dnl                    matrix A is used.
dnl 
dnl   side      (in) blas_side_type
dnl                    determines on which side the hermitian matrix A
dnl                    is multiplied.
dnl
dnl   m, n      (in) int
dnl                    the dimensions of matrices
dnl                    matrix B and C are m-by-n
dnl                    matrix A is m-by-m if side = left, 
dnl                                n-by-n if side = right.
dnl   randomize (in) int
dnl                    if 1, the test case is highly randomized, 
dnl                    if 0, the test case is made to cancel the most.
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
dnl   a, b, c (out) matrices a, b, c.
dnl   lda, ldb, ldc (in) leading dimensions of matrices above
dnl
dnl   seed      (in/out) int
dnl
dnl   HEAD(r_true)  (out) double *
dnl   TAIL(r_true)  (out) double *
dnl               the leading/trailing part of the true in double-double
dnl
include(cblas.m4)dnl
include(test-common.m4)dnl
dnl
dnl
define(`HEMM_NAME', `ifelse(
  `$2&&$3', `$1&&$1', `BLAS_$1hemm$4', 
  `BLAS_$1hemm_$2_$3$4')')dnl
dnl
dnl
define(`SYMM_TESTGEN_NAME', `ifelse(`$2&&$3', `$1&&$1', 
  `BLAS_$1symm_testgen', `BLAS_$1symm_$2_$3_testgen')')dnl
dnl
dnl
dnl  HEMM_TESTGEN
dnl     |
dnl     |-- HEMM_TESTGEN_HEAD
dnl     |      |
dnl     |      |-- HEMM_TESTGEN_NAME
dnl     |      |
dnl     |      |-- HEMM_TESTGEN_PARAMS
dnl     |      
dnl     |-- HEMM_TESTGEN_COMMENT
dnl     |
dnl     |-- HEMM_TESTGEN_BODY
dnl
dnl  Usage:
dnl    HEMM_TESTGEN        ($1, $2, $3)
dnl    HEMM_TESTGEN_HEAD   ($1, $2, $3)
dnl    HEMM_TESTGEN_NAME   ($1, $2, $3)
dnl    HEMM_TESTGEN_PARAMS ($1, $2, $3)
dnl    HEMM_TESTGEN_COMMENT($1, $2, $3)
dnl    HEMM_TESTGEN_BODY   ($1, $2, $3)
dnl
dnl    $1 -- type of alpha, beta, c
dnl    $2 -- type of a
dnl    $3 -- type of b
dnl
define(`HEMM_TESTGEN_NAME', `ifelse(`$2&&$3', `$1&&$1', 
  `BLAS_$1hemm_testgen', `BLAS_$1hemm_$2_$3_testgen')')dnl
dnl
dnl
define(`HEMM_TESTGEN_PARAMS',
  `int norm, enum blas_order_type order, dnl
   enum blas_uplo_type uplo, enum blas_side_type side, dnl
   int m, int n, int randomize, dnl
   $1_array alpha, int alpha_flag, $1_array beta, int beta_flag, dnl
   $2_array a, int lda, $3_array b, int ldb, $1_array c, int ldc, dnl
   int *seed, double *HEAD(r_true), double *TAIL(r_true)')dnl
dnl
dnl
define(`HEMM_TESTGEN_HEAD', 
  `void HEMM_TESTGEN_NAME($1, $2, $3)(HEMM_TESTGEN_PARAMS($1, $2, $3))')dnl
dnl
dnl
define(`HEMM_TESTGEN_COMMENT', `
/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to HEMM_NAME($1, $2, $3, `'){_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_side_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the hermitian matrix a is to be stored.
 *
 * side    (input) enum blas_side_type
 *           which side of matrix b matrix a is to be multiplied.
 *
 * m n     (input) int
 *           sizes of matrices a, b, c:
 *              matrix a is m-by-m for left multiplication
 *                          n-by-n otherwise, 
 *              matrices b, c are m-by-n.
 * 
 * randomize (input) int
 *           = 0: test case made for maximum cancellation.
 *           = 1: test case made for maximum radomization.
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
 * b       (input/output) $3_array
 *
 * ldb     (input) int
 *         leading dimension of matrix B.
 * 
 * c       (input/output) $1_array
 *         generated matrix C that will be used as an input to HEMM.
 * 
 * ldc     (input) int
 *         leading dimension of matrix C.
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
define(`HEMM_TESTGEN_BODY', `{

  /* Strategy:  
         R1 = alpha * A1 * B + beta * C1
         R2 = alpha * A2 * B + beta * C2
     where all the matrices are real.  Then let R = R1 + i R2, 
     A = A1 + i A2, C = C1 + i C2.  To make A hermitian, A1 is 
     symmetric, and A2 is a skew matrix (trans(A2) = -A2).
  */

  int i, j;
  int cij, ci;
  int aij, ai;
  int a1ij, a1i;
  int bij, bi;
  int mij, mi;
  int inccij, incci;
  int incaij, incai;
  int inca1ij, inca1i;
  int incbij, incbi;
  int inci, incij;
  int inca, incb;
  int m_i, n_i;
  int ld;
  int ab;

  DECLARE_VECTOR(a1, REAL_TYPE($2_type))
  DECLARE_VECTOR(a2, REAL_TYPE($2_type))
  DECLARE_VECTOR(c1, REAL_TYPE($1_type))
  DECLARE_VECTOR(c2, REAL_TYPE($1_type))
  DECLARE_VECTOR(b0, REAL_TYPE($3_type))

  DECLARE_VECTOR(r1_true, EXTRA_TYPE(REAL_TYPE($1_type)))
  DECLARE_VECTOR(r2_true, EXTRA_TYPE(REAL_TYPE($1_type)))

  DECLARE(r_elem1, real_E)
  DECLARE(r_elem2, real_E)
  DECLARE(r_elem, real_E)

  DECLARE_VECTOR(a_vec, REAL_TYPE($2_type))
  DECLARE_VECTOR(b_vec, REAL_TYPE($3_type))

  PTR_CAST(c, $1_type)
  PTR_CAST(alpha, $1_type)
  PTR_CAST(beta, $1_type)
  PTR_CAST(a, REAL_TYPE($2_type))
  PTR_CAST(b, REAL_TYPE($3_type))

  if (side == blas_left_side) {
    m_i = m; n_i = n;
  } else {
    m_i = n; n_i = m;
  }

  if (order == blas_colmajor)
    ld = m;
  else
    ld = n;


  if ((order == blas_colmajor && side == blas_left_side) ||
      (order == blas_rowmajor && side == blas_right_side)) {
    inci = inca1i = incai = incbi = incci = 1;
    inccij = ldc;
    incaij = lda;
    incbij = ldb;
    inca1ij = m_i;
    incij = ld;
  } else {
    incci = ldc;
    incai = lda;
    incbi = ldb;
    inca1i = m_i;
    inci = ld;
    incij = inca1ij = incaij = incbij = inccij = 1;
  }

  INC_ADJUST(incci, $1_type)
  INC_ADJUST(inccij, $1_type)
  INC_ADJUST(incai, $2_type)
  INC_ADJUST(incaij, $2_type)
  INC_ADJUST(incbi, $3_type)
  INC_ADJUST(incbij, $3_type)

  inca = incb = 1;
  INC_ADJUST(inca, $2_type)
  INC_ADJUST(incb, $3_type)
  MALLOC_VECTOR(a_vec, $2_type, m_i)
  for (i = 0; i < m_i*inca; i += inca) {
        SET_ZERO_VECTOR_ELEMENT(a_vec, i, $2_type)
  }
  MALLOC_VECTOR(b_vec, $3_type, m_i)
  for (i = 0; i < m_i*incb; i += incb) {
        SET_ZERO_VECTOR_ELEMENT(b_vec, i, $3_type)
  }

  if (randomize == 0) {
      MALLOC_VECTOR(a1, REAL_TYPE($2_type), m_i * m_i)
      MALLOC_VECTOR(a2, REAL_TYPE($2_type), m_i * m_i)
      for (i = 0; i < m_i * m_i; ++i) {
        a1[i] = a2[i] = 0.0;
      }
      MALLOC_VECTOR(c1, REAL_TYPE($1_type), m_i * n_i)
      MALLOC_VECTOR(c2, REAL_TYPE($1_type), m_i * n_i)
      MALLOC_VECTOR(b0, REAL_TYPE($3_type), m_i * n_i)
      for (i = 0; i < m_i * n_i; ++i) {
        c1[i] = c2[i] = b0[i] = 0.0;
      }
      MALLOC_VECTOR(r1_true, EXTRA_TYPE(REAL_TYPE($1_type)), m_i * n_i);
      MALLOC_VECTOR(r2_true, EXTRA_TYPE(REAL_TYPE($1_type)), m_i * n_i);

      /* First generate the real portion of matrix A, and matrix B.
         Note that Re(A) is a symmetric matrix. */
      SYMM_TESTGEN_NAME(REAL_ABBREV($1), REAL_ABBREV($2), REAL_ABBREV($3))
        (norm, order, uplo, side, m, n, 0, alpha_i, alpha_flag, beta_i, 
         beta_flag, a1, m_i, b0, ld, c1, ld, seed, HEAD(r1_true), TAIL(r1_true));

      SKEW_TESTGEN_NAME(REAL_ABBREV($1), REAL_ABBREV($2), REAL_ABBREV($3))
        (norm, order, uplo, side, m, n, alpha_i, 1, beta_i, 1, a2, m_i, b0, ld, 
         c2, ld, seed, HEAD(r2_true), TAIL(r2_true));

      ifelse($3_type, `real_S', `HEMM_TESTGEN_REAL_ADJ($1, $2, $3)', 
             $3_type, `real_D', `HEMM_TESTGEN_REAL_ADJ($1, $2, $3)', 
             $3_type, `complex_S', `HEMM_TESTGEN_COMPLEX_ADJ($1, $2, $3)', 
             $3_type, `complex_D', `HEMM_TESTGEN_COMPLEX_ADJ($1, $2, $3)')
      FREE_VECTOR(a1, REAL_TYPE($2_type))
      FREE_VECTOR(a2, REAL_TYPE($2_type))
      FREE_VECTOR(c1, REAL_TYPE($1_type))
      FREE_VECTOR(c2, REAL_TYPE($1_type))
      FREE_VECTOR(b0, REAL_TYPE($3_type))
      FREE_VECTOR(r1_true, EXTRA_TYPE(REAL_TYPE($1_type)))
      FREE_VECTOR(r2_true, EXTRA_TYPE(REAL_TYPE($1_type)))

  } else {
    /* random matrix */
    DECLARE(a_elem, $2_type)
    DECLARE(b_elem, $3_type)
    DECLARE(c_elem, $1_type)
    DECLARE(r_true_elem, EXTRA_TYPE($1_type))

    /* Since mixed real/complex test generator for dot
       scales the vectors, we need to used the non-mixed
       version if B is real (since A is always complex). */
    IF_REAL($3_type, 
      `DECLARE_VECTOR(bb_vec, COMPLEX_TYPE($3_type))
       MALLOC_VECTOR(bb_vec, COMPLEX_TYPE($3_type), m_i)')

    if (alpha_flag == 0) {
      RANDOM(alpha_i, $1_type, IS_MIXED_PREC($1_type, $2_type, $3_type))
    }
    if (beta_flag == 0) {
      RANDOM(beta_i, $1_type, IS_MIXED_PREC($1_type, $2_type, $3_type))
    }

    /* Fill in matrix A -- Hermitian. */
    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < m_i; j++, aij += incaij) {
        RANDOM(a_elem, $2_type, IS_MIXED_PREC($1_type, $2_type, $3_type))
        SET_VECTOR_ELEMENT(a_i, aij, a_elem, $2_type)
        if (i == j)
          a_i[aij+1] = 0.0;
      }
    }

    /* Fill in matrix B */
    for (i = 0, bi = 0; i < m_i; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n_i; j++, bij += incbij) {
        RANDOM(b_elem, $3_type, IS_MIXED_PREC($1_type, $2_type, $3_type))
        SET_VECTOR_ELEMENT(b_i, bij, b_elem, $3_type)
      }
    }

    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      $2he_copy_row(order, uplo, side, m_i, a, lda, a_vec, i);
      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {

        if (side == blas_left_side)
          $3ge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
        else
          $3ge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);
       
        /* copy the real b_vec into complex bb_vec, so that 
           pure complex test case generator can be called. */
        ifelse(`$3', `s', `HEMM_RANDOM_MIXED_TRUTH($1, $2, c)', 
               `$3', `d', `HEMM_RANDOM_MIXED_TRUTH($1, $2, z)', 
               `
        DOT_TESTGEN_NAME($1, $2, $3)(m_i, m_i, 0, norm, blas_no_conj, 
                                     alpha, 1, beta, 1, a_vec, b_vec, seed, 
                                     PASS_BY_REF(c_elem, $1_type), 
                                     PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));')

        SET_VECTOR_ELEMENT(c_i, cij, c_elem, $1_type)
        SET_VECTOR_ELEMENT(r_true, cij, r_true_elem, EXTRA_TYPE($1_type))
      }
    }    

    IF_REAL($3_type, `FREE_VECTOR(bb_vec, COMPLEX_TYPE($3_type))')

  }

  FREE_VECTOR(a_vec, REAL_TYPE($2_type))
  FREE_VECTOR(b_vec, REAL_TYPE($3_type))
}')dnl
dnl
dnl
define(`HEMM_RANDOM_MIXED_TRUTH', `
        { int k;
          for (k = 0; k < m_i; k++) {
            bb_vec[2*k] = b_vec[k];
            bb_vec[2*k+1] = 0.0;
          }
        }
        DOT_TESTGEN_NAME($1, $2, $3)(m_i, m_i, 0, norm, blas_no_conj, 
                                     alpha, 1, beta, 1, a_vec, bb_vec, seed, 
                                     PASS_BY_REF(c_elem, $1_type), 
                                     PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));')dnl
dnl
dnl
define(`HEMM_TESTGEN_REAL_ADJ', 
`  
  /* The case where B is a real matrix. 

     There are four cases to consider, depending on the 
     values of alpha and beta.

         values                             scaling
      alpha  beta         alpha    A    B    beta    C     R (truth)
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
  for (i = 0, ai = 0, a1i = 0; i < m_i; i++, ai += incai, a1i += inca1i) {
    for (j = 0, aij = ai, a1ij = a1i; j < m_i; 
          j++, aij += incaij, a1ij += inca1ij) {
      a_i[aij] = a1[a1ij];
      a_i[aij+1] = a2[a1ij];
    }
  }

  /* Fill in b */
  for (i = 0, bi = 0, mi = 0; i < m_i; i++, bi += incbi, mi += inci) {
    for (j = 0, bij = bi, mij = mi; j < n_i; 
         j++, bij += incbij, mij += incij) {
      b_i[bij] = b0[mij];    
    }
  }

  /* Fill in c */
  for (i = 0, ci = 0, mi = 0; i < m_i; i++, ci += incci, mi += inci) {
    for (j = 0, cij = ci, mij = mi; j < n_i; 
         j++, cij += inccij, mij += incij) {
      if (ab == 1 || ab == 2) {
        c_i[cij] = -c2[mij];
        c_i[cij+1] = c1[mij];
      } else {
        c_i[cij] = c1[mij];
        c_i[cij+1] = c2[mij];
      }
    }
  }

  /* Fill in truth */
  for (i = 0, ci = 0, mi = 0; i < m_i; i++, ci += incci, mi += inci) {
    for (j = 0, cij = ci, mij = mi; j < n_i; 
         j++, cij += inccij, mij += incij) {

      if (ab == 0 || ab == 1) {
        HEAD(r_true)[cij] = HEAD(r1_true)[mij];
        TAIL(r_true)[cij] = TAIL(r1_true)[mij];
        HEAD(r_true)[cij+1] = HEAD(r2_true)[mij];
        TAIL(r_true)[cij+1] = TAIL(r2_true)[mij];
      } else if (ab == 2) {
        HEAD(r_true)[cij] = -HEAD(r2_true)[mij];
        TAIL(r_true)[cij] = -TAIL(r2_true)[mij];
        HEAD(r_true)[cij+1] = HEAD(r1_true)[mij];
        TAIL(r_true)[cij+1] = TAIL(r1_true)[mij];
      } else {
        HEAD(r_elem1) = HEAD(r1_true)[mij];
        TAIL(r_elem1) = TAIL(r1_true)[mij];

        HEAD(r_elem2) = HEAD(r2_true)[mij];
        TAIL(r_elem2) = TAIL(r2_true)[mij];

        ADD(r_elem, real_E, r_elem1, real_E, r_elem2, real_E)

        /* Set the imaginary part to  R1 + R2 */
        TAIL(r_true)[cij+1] = TAIL(r_elem);
        HEAD(r_true)[cij+1] = HEAD(r_elem);

        /* Set the real part to R1 - R2. */
        SUB(r_elem, real_E, r_elem1, real_E, r_elem2, real_E)
        TAIL(r_true)[cij] = TAIL(r_elem);
        HEAD(r_true)[cij] = HEAD(r_elem);
      }

    }
  }
')dnl
dnl
dnl
define(`HEMM_TESTGEN_COMPLEX_ADJ', 
` 
  /* The case where B is a complex matrix.  Since B is generated
     as a real matrix, we need to perform some scaling.

     There are four cases to consider, depending on the values
     of alpha and beta.

            values                         scaling
         alpha   beta      alpha  A    B       beta    C    R (truth)
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
  for (i = 0, ai = 0, a1i = 0; i < m_i; i++, ai += incai, a1i += inca1i) {
    for (j = 0, aij = ai, a1ij = a1i; j < m_i; 
          j++, aij += incaij, a1ij += inca1ij) {
      a_i[aij] = a1[a1ij];
      a_i[aij+1] = a2[a1ij];
    }
  }

  /* Fill in b */
  for (i = 0, bi = 0, mi = 0; i < m_i; i++, bi += incbi, mi += inci) {
    for (j = 0, bij = bi, mij = mi; j < n_i; 
         j++, bij += incbij, mij += incij) {
      if (ab == 0) {
        b_i[bij] = 0.0;
        b_i[bij+1] = b0[mij];
      } else {
        b_i[bij] = b0[mij];
        b_i[bij+1] = b0[mij];
      }
    }
  }

  /* Fill in c */
  for (i = 0, ci = 0, mi = 0; i < m_i; i++, ci += incci, mi += inci) {
    for (j = 0, cij = ci, mij = mi; j < n_i; 
         j++, cij += inccij, mij += incij) {
      if (ab == 0) {
        c_i[cij] = -c2[mij];
        c_i[cij+1] = c1[mij];
      } else if (ab == 2) {
        c_i[cij] = -2.0 * c2[mij];
        c_i[cij+1] = 2.0 * c1[mij];
      } else {
        c_i[cij] = c1[mij];
        c_i[cij+1] = c2[mij];
      }
    }
  }

  /* Fill in truth */
  for (i = 0, ci = 0, mi = 0; i < m_i; i++, ci += incci, mi += inci) {
    for (j = 0, cij = ci, mij = mi; j < n_i; 
         j++, cij += inccij, mij += incij) {

        HEAD(r_elem1) = HEAD(r1_true)[mij];
        TAIL(r_elem1) = TAIL(r1_true)[mij];

        HEAD(r_elem2) = HEAD(r2_true)[mij];
        TAIL(r_elem2) = TAIL(r2_true)[mij];

      if (ab == 0) {
        HEAD(r_true)[cij] = -HEAD(r_elem2);
        TAIL(r_true)[cij] = -TAIL(r_elem2);
        HEAD(r_true)[cij+1] = HEAD(r_elem1);
        TAIL(r_true)[cij+1] = TAIL(r_elem1);
      } else if (ab == 1) {
        ADD(r_elem, real_E, r_elem1, real_E, r_elem2, real_E)

        /* Set the imaginary part to  R1 + R2 */
        TAIL(r_true)[cij+1] = TAIL(r_elem);
        HEAD(r_true)[cij+1] = HEAD(r_elem);

        /* Set the real part to R1 - R2. */
        SUB(r_elem, real_E, r_elem1, real_E, r_elem2, real_E)
        TAIL(r_true)[cij] = TAIL(r_elem);
        HEAD(r_true)[cij] = HEAD(r_elem);
      } else {

        /* Real part */
        MUL(r_elem, real_E, r_elem2, real_E, -2.0, real_D)
        HEAD(r_true)[cij] = HEAD(r_elem);
        TAIL(r_true)[cij] = TAIL(r_elem);

        /* Imaginary Part */
        MUL(r_elem, real_E, r_elem1, real_E, 2.0, real_D)
        HEAD(r_true)[cij+1] = HEAD(r_elem);
        TAIL(r_true)[cij+1] = TAIL(r_elem);
      }
    }
  }
')dnl
dnl
dnl
dnl
dnl
define(`HEMM_TESTGEN', 
  `HEMM_TESTGEN_HEAD($1, $2, $3)
   HEMM_TESTGEN_COMMENT($1, $2, $3)
   HEMM_TESTGEN_BODY($1, $2, $3)')dnl
dnl
dnl
dnl
dnl
dnl
define(`SKEW_TESTGEN_NAME', `ifelse(`$2&&$3', `$1&&$1', 
  `BLAS_$1skew_testgen', `BLAS_$1skew_$2_$3_testgen')')dnl
dnl
dnl
define(`SKEW_TESTGEN_PARAMS',
  `int norm, enum blas_order_type order, dnl
   enum blas_uplo_type uplo, enum blas_side_type side, dnl
   int m, int n, dnl
   $1_array alpha, int alpha_flag, $1_array beta, int beta_flag, dnl
   $2_array a, int lda, $3_array b, int ldb, $1_array c, int ldc, dnl
   int *seed, double *HEAD(r_true), double *TAIL(r_true)')dnl
dnl
dnl
define(`SKEW_TESTGEN_HEAD', 
  `void SKEW_TESTGEN_NAME($1, $2, $3)(SKEW_TESTGEN_PARAMS($1, $2, $3))')dnl
dnl
dnl
dnl
dnl
define(`SKEW_TESTGEN_BODY', `{

  int i, j;
  int cij, ci;
  int inccij, incci;
  int inca, incb;
  int m_i, n_i;

  DECLARE(c_elem, $1_type)
  DECLARE(r_true_elem, EXTRA_TYPE($1_type))

  DECLARE_VECTOR(a_vec, $2_type)
  DECLARE_VECTOR(b_vec, $3_type)

  PTR_CAST(c, $1_type)

  if (side == blas_left_side) {
    m_i = m; n_i = n;
  } else {
    m_i = n; n_i = m;
  }

  inca = incb = 1;
  INC_ADJUST(inca, $2_type)
  INC_ADJUST(incb, $3_type)
  MALLOC_VECTOR(a_vec, $2_type, m_i)
  for (i = 0; i < m_i*inca; i += inca) {
        SET_ZERO_VECTOR_ELEMENT(a_vec, i, $2_type)
  }
  MALLOC_VECTOR(b_vec, $3_type, m_i)
  for (i = 0; i < m_i*incb; i += incb) {
        SET_ZERO_VECTOR_ELEMENT(b_vec, i, $3_type)
  }

  if ((order == blas_colmajor && side == blas_left_side) ||
      (order == blas_rowmajor && side == blas_right_side)) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }

  INC_ADJUST(incci, $1_type)
  INC_ADJUST(inccij, $1_type)

  if (side == blas_left_side)
    $3ge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, 0);
  else
    $3ge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, 0); 

  /* Fill in matrix A */
  cij = 0;
  for (i = 0; i < m_i; i++, cij += incci) {
    $2skew_copy_row(order, uplo, side, m_i, a, lda, a_vec, i);
    DOT_TESTGEN_NAME($1, $3, $2)(m_i, i, m_i-i, norm, 
                                 blas_no_conj, alpha, 1, 
                                 beta, 1, b_vec, a_vec, seed, 
                                 PASS_BY_REF(c_elem, $1_type), 
                                 PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));

    $2skew_commit_row(order, uplo, side, m_i, a, lda, a_vec, i);

    SET_VECTOR_ELEMENT(c_i, cij, c_elem, $1_type)
    SET_VECTOR_ELEMENT(r_true, cij, r_true_elem, EXTRA_TYPE($1_type))
  }

  /* Now fill in c and r_true */
  
  for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
    for (j = 1, cij = ci + inccij; j < n_i; j++, cij += inccij) {
      COPY_VECTOR_ELEMENT(c_i, cij, c_i, ci, $1_type)
      COPY_VECTOR_ELEMENT(r_true, cij, r_true, ci, EXTRA_TYPE($1_type))
    } 
  }

  FREE_VECTOR(a_vec, $2_type)
  FREE_VECTOR(b_vec, $3_type)
}')dnl
dnl
dnl
dnl
define(`SKEW_TESTGEN', 
  `SKEW_TESTGEN_HEAD($1, $2, $3)
   SKEW_TESTGEN_BODY($1, $2, $3)')dnl
dnl
dnl
define(`PROTOTYPES', `dnl
HEMM_TESTGEN_HEAD(c, c, c);
HEMM_TESTGEN_HEAD(z, z, z);
HEMM_TESTGEN_HEAD(z, c, z);
HEMM_TESTGEN_HEAD(z, z, c);
HEMM_TESTGEN_HEAD(z, c, c);

HEMM_TESTGEN_HEAD(z, z, d);
HEMM_TESTGEN_HEAD(c, c, s);

SKEW_TESTGEN_HEAD(s, s, s);
SKEW_TESTGEN_HEAD(d, d, d);
SKEW_TESTGEN_HEAD(d, d, s);
SKEW_TESTGEN_HEAD(d, s, d);
SKEW_TESTGEN_HEAD(d, s, s);
')dnl
dnl
dnl
define(`SOURCE', `dnl
#include <stdio.h>
#include <stdlib.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"

HEMM_TESTGEN(c, c, c)
HEMM_TESTGEN(z, z, z)
HEMM_TESTGEN(z, c, z)
HEMM_TESTGEN(z, z, c)
HEMM_TESTGEN(z, c, c)

HEMM_TESTGEN(z, z, d)
HEMM_TESTGEN(c, c, s)

SKEW_TESTGEN(s, s, s)
SKEW_TESTGEN(d, d, d)
SKEW_TESTGEN(d, d, s)
SKEW_TESTGEN(d, s, d)
SKEW_TESTGEN(d, s, s)
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
