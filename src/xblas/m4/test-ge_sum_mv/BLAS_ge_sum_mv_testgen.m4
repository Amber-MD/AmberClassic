dnl
dnl ge_sum_mv_testgen.m4
dnl
dnl Test case generator for ge_sum_mv routines.
dnl Generates test cases for alpha, A, beta, B, and x and
dnl computes r_true in double-double precision.
dnl
dnl
dnl C Interface has the form 
dnl
dnl   BLAS_xge_sum_mv{_a_b}{_x}_testgen(int norm, 
dnl                             enum blas_order_type order, 
dnl                             int m, int n, int randomize, 
dnl                             SCALAR *alpha, int alpha_flag, 
dnl                             SCALAR *beta, int beta_flag, 
dnl                             ARRAY a, int lda, 
dnl                             ARRAY b, int ldb,
dnl                             ARRAY x, int incx, 
dnl                             SCALAR *alpha_use_ptr, 
dnl                             ARRAY a_use, ARRAY b_use
dnl                             int *seed, 
dnl                             double *HEAD(r_true), double *TAIL(r_true))
dnl
dnl Arguments
dnl   norm      (in) int
dnl
dnl   order     (in) blas_order_type
dnl                    determines the storage format for matrices
dnl   
dnl   m, n      (in) int, int
dnl                    the size of the
dnl                    vector x is n
dnl                    matrix A is m-by-n.
dnl   
dnl   randomize (in) int
dnl                  if 0, entries in matrices A, B will be chosen for
dnl                        maximum cancellation, but with less randomness.
dnl                  if 1, every entry in the matrix A, B will be 
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
dnl   lda       (in) leading dimension of matrix a.
dnl
dnl   b         (out) matrix b.
dnl   ldb       (in) leading dimension of matrix b.
dnl
dnl   x         (out) vector x.
dnl   incx      (in) stride for vector x.
dnl
dnl   alpha_use_ptr (out) pointer to the scalar
dnl             this routine puts in the alpha or beta used
dnl             before scaling (see strategy below)
dnl
dnl   a_use     (out) these arrays 
dnl   b_use             get the matricies a, b before any scaling
dnl                     they have the same dimensions and lda/ldb
dnl                     as a, b, respectively
dnl
dnl   seed      (in/out) int
dnl
dnl   HEAD(r_true)  (out) double *  (these are vectors of size n)
dnl   TAIL(r_true)  (out) double *
dnl               the leading/trailing part of the true in double-double
dnl
dnl
dnl
include(cblas.m4)dnl
include(test-common.m4)dnl
dnl
dnl
dnl   DEREF - 
dnl  This macro allows dereferencing of number pointers,
dnl    if the number is real, the name is returned dereferenced
dnl    if the number is complex, the name is returned still referenced,
dnl     so that it can be array subscripted.
define(`DEREF', `IF_COMPLEX($2_type, `$1', `(*$1)')')dnl
dnl
dnl
define(`GE_SUM_MV_NAME', `ifelse(
  `$2&&$3', `$1&&$1', `BLAS_$1ge_sum_mv$4', 
  `BLAS_$1ge_sum_mv_$2_$3$4')')dnl
dnl
define(`GEMV_TESTGEN_NAME', `ifelse(
  `$2&&$3', `$1&&$1', `BLAS_$1gemv_testgen', 
  `BLAS_$1gemv_$2_$3_testgen')')dnl
dnl
dnl
dnl  GE_SUM_MV_TESTGEN
dnl     |
dnl     |-- GE_SUM_MV_TESTGEN_HEAD
dnl     |      |
dnl     |      |-- GE_SUM_MV_TESTGEN_NAME
dnl     |      |
dnl     |      |-- GE_SUM_MV_TESTGEN_PARAMS
dnl     |      
dnl     |-- GE_SUM_MV_TESTGEN_COMMENT
dnl     |
dnl     |-- GE_SUM_MV_TESTGEN_BODY
dnl
dnl  Usage:
dnl    GE_SUM_MV_TESTGEN        ($1, $2, $3)
dnl    GE_SUM_MV_TESTGEN_HEAD   ($1, $2, $3)
dnl    GE_SUM_MV_TESTGEN_NAME   ($1, $2, $3)
dnl    GE_SUM_MV_TESTGEN_PARAMS ($1, $2, $3)
dnl    GE_SUM_MV_TESTGEN_COMMENT($1, $2, $3)
dnl    GE_SUM_MV_TESTGEN_BODY   ($1, $2, $3)
dnl
dnl    $1 -- type of alpha, beta, y (output of ge_sum_mv routine)
dnl    $2 -- type of a, b 
dnl    $3 -- type of x
dnl
define(`GE_SUM_MV_TESTGEN_NAME', `ifelse(`$2&&$3', `$1&&$1', 
  `BLAS_$1ge_sum_mv_testgen', `BLAS_$1ge_sum_mv_$2_$3_testgen')')dnl
dnl
dnl
define(`GE_SUM_MV_TESTGEN_PARAMS',
  `int norm, enum blas_order_type order, 
   int m, int n,  int randomize,
   $1_array alpha, int alpha_flag, $1_array beta, int beta_flag, 
   $2_array a, int lda, $2_array b, int ldb, $3_array x, int incx, 
   $1_array alpha_use_ptr, $2_array a_use, $2_array b_use, 
   int *seed, double *HEAD(r_true), double *TAIL(r_true)')dnl
dnl
dnl
define(`GE_SUM_MV_TESTGEN_HEAD', 
  `void GE_SUM_MV_TESTGEN_NAME($1, $2, $3)(GE_SUM_MV_TESTGEN_PARAMS($1, $2, $3))')dnl
dnl
dnl
define(`GE_SUM_MV_TESTGEN_COMMENT', `
/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to GE_SUM_MV_NAME($1, $2, $3, `'){_x}
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
 * m, n    (input) int
 *              vector x is length n.
 *              Matricies A, B are size m-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, B will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, B will be 
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
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * b       (input/output) $2_array
 * 
 * ldb     (input) ldb
 *         leading dimension of matrix B.
 *
 * x       (input/output) $3_array
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * alpha_use_ptr (output) $1_array
 *              must contain a valid pointer. 
 *              used to return the value of alpha, beta before scaling
 *              (see strategy below)
 *
 * a_use   (output) $2_array
 *              matrix of dimension m by n, leading dimension lda.
 *              a_use will get the a matrix before any scaling.
 *
 * b_use   (output) $2_array
 *              matrix of dimension m by n, leading dimension ldb.
 *              b_use will get the b matrix before any scaling.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * double  (output) *HEAD(r_true)
 *         the leading part of the truth in double-double.
 *
 * double  (output) *TAIL(r_true)
 *         the trailing part of the truth in double-double
 *
 *
 * strategy :
 * the test generation for ge_sum_mv is broken up into cases.
 * first off, x is selected randomly, and put twice into
 * a vector of length 2*n, x_vec.  x_vec will be used in most 
 * cases in the call to the dot test generator.
 *
 * Then, the case is determined, and the type is stored in
 *      case_type.
 *
 * Note that ge_sum_mv is symmetric with respect to matricies a, b.
 *   
 *
 * IF_REAL($1_type, `
 *cases:  alpha, beta are real:
 * case 1: alpha, beta are free:
 *    In this case, we select alpha randomly, and make 
 *      beta = (2^k) * alpha, where k is an 
 *      integer between +- 4.  
 *      The generator is run as if alpha == beta, 
 *      with dot products with length 2*n,
 *      and then afterwards each element in B is scaled
 *      by (2^(-k)).
 * case 2: alpha = 0, beta not 0 (alpha not zero, beta = 0):
 *    This case degrades into the GEMV case, with beta=0.0.
 *    the matrix a_use (b_use) is set to zero, and
 *    a (b) is filled with random numbers. 
 * case 3: alpha = 1, beta free (or alpha free, beta = 1):
 *    This case is treated similar to case 1. alpha (beta) is
 *    held fixed, and beta (alpha) becomes (2^k)*alpha ((2^k)*beta).
 * case 4: alpha = 1, beta = 1
 *    This case is treated as in case 1, with k = 0. no scaling
 *    is done.
 *', `
 *cases:  alpha, beta are complex:
 * case 1: alpha, beta are free:
 *    In this case, we select alpha randomly, and make 
 *      beta = (2^k) * alpha, where k is an 
 *      integer between +- 4.  
 *      The generator is run as if alpha == beta, 
 *      with dot products with length 2*n,
 *      and then afterwards each element in B is scaled
 *      by (2^(-k)).
 * case 2: alpha = 0, beta not 0 (alpha not zero, beta = 0):
 *    This case degrades into the GEMV case, with beta=0.0.
 *    the matrix a_use (b_use) is set to zero, and
 *    a (b) is filled with random numbers. 
 * case 3: alpha = 1, beta free (or alpha free, beta = 1):
 *ifelse(REAL_COMPLEX($1_type)&&REAL_COMPLEX($2_type), `complex&&real', `
 *    THIS CASE IS NOT PROPERLY TESTED.
 *      because of the difficulty in testing this case,
 *      a call with this case and randomize = 0 is
 *      converted into a call with randomize = 1.
 *      THERE IS INSUFFICIENT TESTING OF CANCELLATION IN THIS CASE.
 *      It is suggested that implementors be aware of this
 *      and take caution when working on ge_sum_mv.
',`
 *    This becomes tricky; In this case,
 *      When randomize == 1, treat similar to case 1.
 *      When randomize == 0,
 *        k is determined as usual. 
 *        x_vec is selected real randomly,
 *        then a, b, are generated real for cancellation,
 *          and the truth is obtained (at this point, it is real)
 *        x_vec is scaled by 1+i.
 *        the truth is scaled by 1+i.
 *        b (a) is scaled by (2^-(k+1))*(1+i)
 *        beta (alpha) is scaled by (2^k)*(1-i)
 *        because (1+i)*(1-i) == 2+0i.
')dnl
 * case 4: alpha = 1, beta = 1
 *    This case is treated as in case 1, with k = 0. no scaling
 *    is done.')
 */
')dnl
dnl
dnl
define(`GE_SUM_MV_TESTGEN_BODY', `{

  int i, j, k;
  int xi;
  int aij, ai, ri;
  int incri;
  int incxi, incx_veci, x_starti;
  int incaij, incai;
  int inca_veci;
  int n_i, m_i;
  int case_type;
  int which_free;

  DECLARE(y_elem, $1_type)
  DECLARE(beta_zero_fake, $1_type)
  DECLARE(a_elem, $2_type)
  DECLARE(x_elem, $3_type)
  DECLARE(r_true_elem, EXTRA_TYPE($1_type))
  DECLARE(multiplier, REAL_TYPE($1_type))
  DECLARE(divider, REAL_TYPE($1_type))
  DECLARE(alpha_use, $1_type)

  DECLARE_VECTOR(a_vec, $2_type)
  DECLARE_VECTOR(x_vec, $3_type)

  PTR_CAST(alpha_use_ptr, $1_type)
  PTR_CAST(alpha, $1_type)
  PTR_CAST(beta, $1_type)
  PTR_CAST(a, $2_type)
  PTR_CAST(b, $2_type)
  PTR_CAST(a_use, $2_type)
  PTR_CAST(b_use, $2_type)
  PTR_CAST(x, $3_type)

  n_i = n;
  m_i = m;

  ZERO(beta_zero_fake, $1_type)

  /*x_vec, a_vec must have stride of 1*/
  inca_veci = 1;
  INC_ADJUST(inca_veci, $2_type)

  MALLOC_VECTOR(a_vec, $2_type, 2*n_i)

  incri = 1;
  INC_ADJUST(incri, EXTRA_TYPE($1_type))

  incxi = incx;
  incx_veci = 1;
  INC_ADJUST(incx_veci, $3_type)
  INC_ADJUST(incxi, $3_type)

  if (incxi < 0) {
    x_starti = (-n+1) * incxi;
  } else {
    x_starti = 0;
  }

  MALLOC_VECTOR(x_vec, $3_type, 2*n_i);

        /* choose k */ 
    k = 0;
    while(!k) { k = xrand(seed) * 7 - 3; }

    ONE(multiplier, REAL_TYPE($1_type))
    ONE(divider, REAL_TYPE($1_type))
    for(i =0; i < k; i++ ) {
        MUL(multiplier, REAL_TYPE($1_type), multiplier, 
                REAL_TYPE($1_type), 2.0, REAL_TYPE($1_type))
        MUL(divider, REAL_TYPE($1_type), divider, REAL_TYPE($1_type), 
                0.5, REAL_TYPE($1_type))
    }
    for(i=0; i > k; i--) {
        MUL(multiplier, REAL_TYPE($1_type), 
                multiplier, REAL_TYPE($1_type), 0.5, REAL_TYPE($1_type))
        MUL(divider, REAL_TYPE($1_type), divider, REAL_TYPE($1_type), 
                2.0, REAL_TYPE($1_type))
    }
        /* decide which case */
    if(alpha_flag) {
        if(TEST_0(DEREF(alpha_i, $1), $1_type)) {
                /* case 2*/
                case_type = 2;
                which_free = ALPHA_USE_IS_BETA; /* for use beta */
        } else {
                if(beta_flag) {
                        if(TEST_0(DEREF(beta_i, $1), $1_type)) {
                                /* case 2 */
                                case_type = 2; 
                                which_free = ALPHA_USE_IS_ALPHA; 
                                        /*for use alpha*/
                        } else {
                                /* case 4 */
                                case_type = 4;
                                k = 0;
                                which_free = ALPHA_USE_IS_EITHER;
                        }
                } else {
                        /* case 3 */
                        case_type = 3;
                        which_free = ALPHA_USE_IS_ALPHA; 
                                /* for beta free, use alpha*/
                }
        }
    } else {
        if(beta_flag) {
                if(TEST_0(DEREF(beta_i, $1), $1_type)) {
                        /* case 2 */
                        case_type = 2;
                        which_free = ALPHA_USE_IS_ALPHA; 
                        /*alpha is nonzero */
                } else {
                        /* case 3 */
                        case_type = 3;
                        which_free = ALPHA_USE_IS_BETA; 
                                /* for alpha free, use beta*/
                }
        } else {
                /* case 1 */
                case_type = 1;
                which_free = ALPHA_USE_IS_ALPHA;
        }
    }

    if(which_free == ALPHA_USE_IS_BETA) {
        if(!beta_flag) {
              RANDOM(y_elem, $1_type, IS_MIXED_PREC($1_type, $2_type, $3_type))
              SET_VECTOR_ELEMENT(beta_i, 0, y_elem, $1_type)
        }
        ASSIGN(alpha_use, $1_type, DEREF(beta_i, $1), $1_type)
    } else {
        if(!alpha_flag) {
              RANDOM(y_elem, $1_type, IS_MIXED_PREC($1_type, $2_type, $3_type))
              SET_VECTOR_ELEMENT(alpha_i, 0, y_elem, $1_type)
        }
        ASSIGN(alpha_use, $1_type, DEREF(alpha_i, $1), $1_type)
    }
        /* put in return value */
    ASSIGN(DEREF(alpha_use_ptr_i, $1), $1_type, alpha_use, $1_type)

  if (randomize == 0) {

        /*first pick x randomly*/
    for (i = 0, xi = x_starti; i < n_i; i++, xi ++) {
      RANDOM(x_elem, $3_type, IS_MIXED_PREC($1_type, $2_type, $3_type))
      SET_VECTOR_ELEMENT(x_i, xi*incxi, x_elem, $3_type)
    }
        /*copy new x into x_vec (twice) */
    $3copy_vector(x, n_i, incx, x_vec, 1);
    $3copy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);
        
    if(case_type == 2) {
        /* degenerate case - similar to gemv */
        if(which_free == ALPHA_USE_IS_ALPHA) {
                /* alpha == alpha_use */
                DEGENERATE_TESTGEN($1, $2, $3, alpha, a, lda) 

                /*now fill a, x, and return */
                FILL_MATRIX_RANDOM(b, $2, ldb)
                FILL_WITH_ZERO_MATRIX(b_use, ldb, $2)
                COPY_MATRIX(a_use, a, lda, $2)
                ASSIGN(DEREF(alpha_use_ptr_i, $1), $1_type, 
                        alpha_use, $1_type)
                ASSIGN(DEREF(alpha_i, $1), $1_type, alpha_use, $1_type)
                $3copy_vector(x_vec, n_i, 1, x, incx);
                FREE_VECTOR(a_vec, $2_type)
                FREE_VECTOR(x_vec, $3_type)
                return; 
        } else {
                DEGENERATE_TESTGEN($1, $2, $3, beta, b, ldb) 

                /*now fill b, x, and return */
                FILL_MATRIX_RANDOM(a, $2, lda)
                FILL_WITH_ZERO_MATRIX(a_use, lda, $2)
                COPY_MATRIX(b_use, b, ldb, $2)
                ASSIGN(DEREF(alpha_use_ptr_i, $1), $1_type, 
                        alpha_use, $1_type)
                ASSIGN(DEREF(beta_i, $1), $1_type, alpha_use, $1_type)
                $3copy_vector(x_vec, n_i, 1, x, incx);
                FREE_VECTOR(a_vec, $2_type)
                FREE_VECTOR(x_vec, $3_type)
                return; 
        }
    }

    ifelse(REAL_COMPLEX($1_type)&&REAL_COMPLEX($2_type),
        `complex&&real',
        `       /* case 3, Not fully tested */
        if(case_type == 3) {
                GE_SUM_MV_TESTGEN_NAME($1, $2, $3, $4)(
                        norm, order, m, n, 1/*randomize*/,
                        alpha_i, alpha_flag, beta_i, beta_flag,
                        a, lda, b, ldb, x, incx, 
                        alpha_use_ptr_i, a_use, b_use, 
                        seed, HEAD(r_true), TAIL(r_true));
                FREE_VECTOR(a_vec, $2_type)
                FREE_VECTOR(x_vec, $3_type)
                return;
        }')

    ifelse(REAL_COMPLEX($1_type)&&REAL_COMPLEX($2_type),
        `complex&&complex',
        `       /* case 3, start with real matricies, x */
        if(case_type == 3) {
                DECLARE_VECTOR(a_vec_2, $2_type)
                DECLARE_VECTOR(x_vec_2, $3_type)
                MALLOC_VECTOR(a_vec_2, $2_type, 4*n_i)
                MALLOC_VECTOR(x_vec_2, $3_type, 4*n_i)
                for (i = 0; i < 2*n_i*inca_veci; i+= inca_veci) {
                 SET_ZERO_VECTOR_ELEMENT(a_vec, i, $2_type)
                }
        
                /*first pick x randomly, but real*/
            for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
              RANDOM(x_elem, $3_type, IS_MIXED_PREC($1_type, $2_type, $3_type))
              IF_COMPLEX($3_type, `x_elem[1] = 0.0;')
              SET_VECTOR_ELEMENT(x_i, xi, x_elem, $3_type)
            }
                /*copy new x into x_vec_2 (twice) */
            REAL_ABBREV($3)copy_vector(x, n_i, 2*incx, x_vec_2, 1);
            REAL_ABBREV($3)copy_vector(x_vec_2, n_i, 1, (x_vec_2 + n_i), 1);

                /* Now Fill in matrix A, B real*/
                /*since we have case 3, we know alpha_use == 1.0+0i,
                        so we will force it to be real */
                for(i = 0, ri = 0; i < m_i; i++, ri += incri) {
                ZERO(y_elem, $1_type)
                DOT_TESTGEN_NAME(REAL_ABBREV($1), REAL_ABBREV($3), 
                        REAL_ABBREV($2))(2*n_i, 0, 2*n_i, norm, 
                                   blas_no_conj, 
                                   PASS_BY_REF(alpha_use, $1_type), 1, 
                                   PASS_BY_REF(beta_zero_fake, $1_type), 1, 
                                   x_vec_2, 
                                   a_vec_2, seed, 
                                   PASS_BY_REF(y_elem, $1_type), 
                                   PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));
        
                IF_COMPLEX($3_type, `
                  /*multiply truth by 1+i (we will multiply 1+i to x later)*/
                  HEAD(r_true_elem)[1] = HEAD(r_true_elem)[0];
                  TAIL(r_true_elem)[1] = TAIL(r_true_elem)[0];', 
                 `HEAD(r_true_elem)[1] = TAIL(r_true_elem)[1] = 0.0;')
                for(j=0; j<2*n_i; j++) {
                        a_vec[2*j] = a_vec_2[j];
                        a_vec[2*j+1] = 0.0;
                }
                $2ge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
                $2ge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb, a_vec + inca_veci* n_i, i);

                /*commits an element to the truth */
                SET_VECTOR_ELEMENT(r_true, ri, r_true_elem, EXTRA_TYPE($1_type))
                }
                /* copy to x_vec - will be copied to x_i later */
                IF_COMPLEX($3_type, `
                /* also multiply x by 1+i, to compensate for change in
                        truth above */
                for(j=0; j<n_i; j++) {
                        x_vec[2*j] = x_vec_2[j];
                        x_vec[2*j+1] = x_vec_2[j];
                }',
                `for(j=0; j<n_i; j++) {
                        x_vec[j] = x_vec_2[j];
                }')
                FREE_VECTOR(x_vec_2, $2_type)
                FREE_VECTOR(a_vec_2, $3_type)
        } else {
                /*not case 3 */
                REGULAR_GENERATION($1, $2, $3, $4, $5)
        }', 
        `REGULAR_GENERATION($1, $2, $3, $4, $5)')

  } else {
        /* randomize == 1 */
    ifelse(IS_MIXED($1, $2), `t', 
      `DECLARE_VECTOR(aa_vec, COMPLEX_TYPE($2_type))')
    ifelse(IS_MIXED($1, $3), `t', 
      `DECLARE_VECTOR(xx_vec, COMPLEX_TYPE($3_type))')

    ifelse(IS_MIXED($1, $2), `t', 
      `MALLOC_VECTOR(aa_vec, COMPLEX_TYPE($2_type), 2 * n_i)')
    ifelse(IS_MIXED($1, $3), `t', 
      `MALLOC_VECTOR(xx_vec, COMPLEX_TYPE($3_type), 2 * n_i)')


        /*first pick x randomly*/
    for (i = 0, xi = x_starti; i < n_i; i++, xi++) {
      RANDOM(x_elem, $3_type, IS_MIXED_PREC($1_type, $2_type, $3_type))
      SET_VECTOR_ELEMENT(x_i, xi*incxi, x_elem, $3_type)
    }
    ifelse(REAL_COMPLEX($1_type)&&REAL_COMPLEX($2_type)&&$1,
        `complex&&complex&&z',
        `if(case_type == 3) {
            FILL_MATRIX_RANDOM(a, $2, lda, real)
            FILL_MATRIX_RANDOM(b, $2, ldb, real)
        } else {
            FILL_MATRIX_RANDOM(a, $2, lda)
            FILL_MATRIX_RANDOM(b, $2, ldb)
        }',
        `
            FILL_MATRIX_RANDOM(a, $2, lda)
            FILL_MATRIX_RANDOM(b, $2, ldb)')

    /* now compute appropriate truth */
      
      /* get x */
        /*copy new x into x_vec (twice) */
      $3copy_vector(x, n_i, incx, x_vec, 1);
      $3copy_vector(x, n_i, incx, (x_vec + incx_veci * n_i), 1);
      ifelse(IS_MIXED($1, $3), `t', 
        `{
          /* promote to complex */
         int r;
         for (r = 0; r < 2*n_i; r++) {
           xx_vec[2*r] = x_vec[r];
           xx_vec[2*r+1] = 0.0;
         }
        }')
        
    if(case_type == 2) {
        if(which_free == ALPHA_USE_IS_BETA) {
                DEGENERATE_RANDOM_TESTGEN($1, $2, $3, beta, b, ldb)
                FILL_WITH_ZERO_MATRIX(a_use, lda, $2)
                COPY_MATRIX(b_use, b, ldb, $2)
                ASSIGN(DEREF(alpha_use_ptr_i, $1), $1_type, 
                        alpha_use, $1_type)
                $3copy_vector(x_vec, n_i, 1, x, incx);
                FREE_VECTOR(a_vec, $2_type)
                FREE_VECTOR(x_vec, $3_type)
                ifelse(IS_MIXED($1, $2), `t', `FREE_VECTOR(aa_vec, $2_type)')
                ifelse(IS_MIXED($1, $3), `t', `FREE_VECTOR(xx_vec, $3_type)')
                return; 
                
        } else {
                DEGENERATE_RANDOM_TESTGEN($1, $2, $3, alpha, a, lda)
                FILL_WITH_ZERO_MATRIX(b_use, ldb, $2)
                COPY_MATRIX(a_use, a, lda, $2)
                ASSIGN(DEREF(alpha_use_ptr_i, $1), $1_type, 
                        alpha_use, $1_type)
                $3copy_vector(x_vec, n_i, 1, x, incx);
                FREE_VECTOR(a_vec, $2_type)
                FREE_VECTOR(x_vec, $3_type)
                ifelse(IS_MIXED($1, $2), `t', `FREE_VECTOR(aa_vec, $2_type)')
                ifelse(IS_MIXED($1, $3), `t', `FREE_VECTOR(xx_vec, $3_type)')
                return; 
        }
    } else {    
      for (i = 0, ri = 0; 
          i < m_i; i++, ri += incri) {
        $2ge_copy_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
        $2ge_copy_row(order, blas_no_trans, m_i, n_i, b, ldb, (a_vec + inca_veci*n_i), i);

        ifelse(IS_MIXED($1, $2), `t', 
          `{
                  /* promote to complex */
           int r;
           for (r = 0; r < 2 * n_i; r++) {
           aa_vec[2*r] = a_vec[r];
             aa_vec[2*r+1] = 0.0;
           }
          }')

        ZERO(y_elem, $1_type)
        ifelse(IS_MIXED($1, $2, $3), `t', 
        `DOT_TESTGEN_NAME($1, $1, $1)(2*n_i, 2*n_i, 0, norm, blas_no_conj, 
                       &alpha_use, 1, 
                       &beta_zero_fake, 1, 
                       ifelse(IS_MIXED($1, $3), `t', `xx_vec', `x_vec'), 
                       ifelse(IS_MIXED($1, $2), `t', `aa_vec', `a_vec'), 
                       seed, 
                       PASS_BY_REF(y_elem, $1_type), 
                       PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));',
        `DOT_TESTGEN_NAME($1, $3, $2)(2*n_i, 2*n_i, 0, norm, blas_no_conj, 
                       &alpha_use, 1, 
                       &beta_zero_fake, 1, 
                       x_vec, a_vec, seed, 
                       PASS_BY_REF(y_elem, $1_type), 
                       PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));')

        SET_VECTOR_ELEMENT(r_true, ri, r_true_elem, EXTRA_TYPE($1_type))
      }
    }
    ifelse(IS_MIXED($1, $2), `t', `FREE_VECTOR(aa_vec, $2_type)')
    ifelse(IS_MIXED($1, $3), `t', `FREE_VECTOR(xx_vec, $3_type)')
  }


        COPY_MATRIX(a_use, a, lda, $2)
        COPY_MATRIX(b_use, b, ldb, $2)
        ASSIGN(DEREF(alpha_use_ptr_i, $1), $1_type, alpha_use, $1_type)


        /* now we scale */
    if(which_free == ALPHA_USE_IS_BETA) {
        IF_COMPLEX($1_type, `
          SCALE_MATRIX_COMPLEX(a, lda, $2, $1, 
           `ifelse(REAL_COMPLEX($1_type)&&REAL_COMPLEX($2_type), `complex&&complex', `complex')')
          SCALE_MULTIPLIER_COMPLEX(alpha_i, beta_i, $1,
           `ifelse(REAL_COMPLEX($1_type)&&REAL_COMPLEX($2_type), `complex&&complex', `complex')')', 
         `SCALE_MATRIX_REAL(a, lda, $2, $1)
          SCALE_MULTIPLIER_REAL(DEREF(alpha_i, $1), DEREF(beta_i, $1), $1)')
    } else {
        if(which_free == ALPHA_USE_IS_ALPHA) {
        IF_COMPLEX($1_type, `
          SCALE_MATRIX_COMPLEX(b, ldb, $2, $1,
           `ifelse(REAL_COMPLEX($1_type)&&REAL_COMPLEX($2_type), `complex&&complex', `complex')')
          SCALE_MULTIPLIER_COMPLEX(beta_i, alpha_i, $1,
           `ifelse(REAL_COMPLEX($1_type)&&REAL_COMPLEX($2_type), `complex&&complex', `complex')')', 
         `SCALE_MATRIX_REAL(b, ldb, $2, $1)
          SCALE_MULTIPLIER_REAL(DEREF(beta_i, $1), DEREF(alpha_i, $1), $1)')
        } else {
        /*which_free = ALPHA_USE_IS_EITHER , case 4 */
        }
    } /* which_free if */

        /*copy x_vec into x : it is possible that the generator
                changed x_vec, even though none were free*/
    $3copy_vector(x_vec, n_i, 1, x, incx);
  FREE_VECTOR(a_vec, $2_type)
  FREE_VECTOR(x_vec, $3_type)
}')dnl
dnl
dnl
dnl
define(`FILL_MATRIX_RANDOM', `
    /*set $1 randomly */
    if(order == blas_colmajor)
    {
        incai = 1;
        incaij = $3; dnl lda/ldb
    } else {
        incai = $3;
        incaij = 1;
    }   
    INC_ADJUST(incai, $2_type)
    INC_ADJUST(incaij, $2_type)

    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
        ifelse($4, `real', 
                `RANDOM(a_elem[0], REAL_TYPE($2_type), IS_MIXED_PREC($1_type, $2_type, $3_type))
                a_elem[1] = 0.0;',
                `RANDOM(a_elem, $2_type, IS_MIXED_PREC($1_type, $2_type, $3_type)) ')
        SET_VECTOR_ELEMENT($1_i, aij, a_elem, $2_type) 
      }
    }')dnl
dnl
dnl
define(`COPY_MATRIX', `
    /*set $1 = $2 */
    if(order == blas_colmajor)
    {
        incai = 1;
        incaij = $3; dnl lda/ldb
    } else {
        incai = $3;
        incaij = 1;
    }   
    INC_ADJUST(incai, $4_type)
    INC_ADJUST(incaij, $4_type)

    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
        GET_VECTOR_ELEMENT(a_elem, $2_i, aij, $4_type)
        SET_VECTOR_ELEMENT($1_i, aij, a_elem, $4_type) 
      }
    }')dnl
dnl
define(`FILL_WITH_ZERO_MATRIX', `
    /*set $1 = 0 */
    if(order == blas_colmajor)
    {
        incai = 1;
        incaij = $2; dnl lda/ldb
    } else {
        incai = $2;
        incaij = 1;
    }   
    INC_ADJUST(incai, $3_type)
    INC_ADJUST(incaij, $3_type)

    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
        SET_ZERO_VECTOR_ELEMENT($1_i, aij, $3_type) 
      }
    }')dnl
dnl
define(`GE_SUM_MV_TESTGEN', 
  `GE_SUM_MV_TESTGEN_HEAD($1, $2, $3)
   GE_SUM_MV_TESTGEN_COMMENT($1, $2, $3)
   GE_SUM_MV_TESTGEN_BODY($1, $2, $3)')dnl
dnl
dnl ************************************************************
dnl     SCALING M4 PROCEDURES.
dnl
dnl ************************************************************
dnl  SCALE_MULTIPLIER_REAL - scales beta when beta is real
dnl
dnl $1 is the variable (beta or alpha) to be assigned from $2
define(`SCALE_MULTIPLIER_REAL',
        `
        switch(case_type) {
        case 1:
        case 3:
                ASSIGN($2, $3_type, alpha_use, $3_type)
                MUL($1, $3_type, $2, $3_type, multiplier, REAL_TYPE($3_type))
                break;
        case 2: /*should not happen */
        case 4:
                break;
        }')dnl
dnl
dnl
define(`SCALE_MULTIPLIER_COMPLEX',
        `
        switch(case_type) {
        case 1:
                ASSIGN($2, $3_type, alpha_use, $3_type)
                MUL($1, $3_type, $2, $3_type, multiplier, REAL_TYPE($3_type))
                break;
        case 2: /*should not happen */
                break;
        case 3:
                ASSIGN($2, $3_type, alpha_use, $3_type)
                ifelse(`complex', $4,
                `MUL($1, $3_type, $2, $3_type, multiplier, 
                        REAL_TYPE($3_type))
                $1[1] = $1[0];
                break;',
                `MUL($1, $3_type, $2, $3_type, 
                        multiplier, REAL_TYPE($3_type))
                break;')
        case 4:
                break;
        }')dnl
dnl
dnl
define(`SCALE_MATRIX_REAL',
`{
    if(order == blas_colmajor)
    {
        incai = 1;
        incaij = $2; dnl lda/ldb
    } else {
        incai = $2;
        incaij = 1;
    }   
    INC_ADJUST(incai, $3_type)
    INC_ADJUST(incaij, $3_type)

    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
        GET_VECTOR_ELEMENT(a_elem, $1_i, aij, $3_type)
        switch(case_type) {
        case 1:
        case 3:
                MUL(a_elem, $3_type, a_elem, $3_type, divider, 
                        REAL_TYPE($3_type))
                break;
        case 2: /*should not happen */
        case 4: /*k ==0 */
    break;
        }
        SET_VECTOR_ELEMENT($1_i, aij, a_elem, $3_type) 
      }
    }    
}')dnl
dnl
dnl
define(`SCALE_MATRIX_COMPLEX',
`{
    ifelse($5, `complex', `
                DECLARE(one_minus_i, $3_type)
                DECLARE(a_elem_2, EXTRA_TYPE($3_type))
                DECLARE(a_elem_3, EXTRA_TYPE($3_type))
                one_minus_i[0] = 0.5;
                one_minus_i[1] = -0.5;') dnl

    if(order == blas_colmajor)
    {
        incai = 1;
        incaij = $2; dnl lda/ldb
    } else {
        incai = $2;
        incaij = 1;
    }   
    INC_ADJUST(incai, $3_type)
    INC_ADJUST(incaij, $3_type)

    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
        GET_VECTOR_ELEMENT(a_elem, $1_i, aij, $3_type)
        switch(case_type) {
        case 1:
                MUL(a_elem, $3_type, a_elem, $3_type, divider, REAL_TYPE($3_type))
                break;
        case 2: /*should not happen */
        case 3:
                ifelse($5, `complex', dnl must use a_elem_2, otherwise
                                        dnl multiplication does not work
                                        dnl properly
                `MUL(a_elem_2, EXTRA_TYPE($3_type), a_elem, $3_type, divider, 
                        REAL_TYPE($3_type))
                MUL(a_elem_3, EXTRA_TYPE($3_type), 
                        a_elem_2, EXTRA_TYPE($3_type), one_minus_i,
                        $3_type)
                ROUND(a_elem, $3_type, a_elem_3, EXTRA_TYPE($3_type)) 
                break;',
                `
                MUL(a_elem, $3_type, a_elem, $3_type, divider, REAL_TYPE($3_type))
                break;')
        case 4: /*k ==0 */
      break;
        }
        SET_VECTOR_ELEMENT($1_i, aij, a_elem, $3_type) 
      }
    }    
}')dnl
dnl
dnl
dnl ******************************
dnl     DEGENERATE_TESTGEN
dnl
dnl  Purpose:  Handles the #2 case, where one of alpha, beta is zero
dnl
dnl
dnl  arguments :        $1 = type of y, alpha, beta, 
dnl                     $2 = type of Matrix, 
dnl                     $3 = type of x, 
dnl                     $4 = name of nonzero-scalar (alpha or beta), 
dnl                     $5 = name of important matrix (a or b), 
dnl                     $6 = name of leading dimension 
dnl                             of important matrix (lda/ldb) 
define(`DEGENERATE_TESTGEN',
`
        /* now Fill in matrix $4 only */
    for(i = 0, ri = 0; i < m_i; i++, ri += incri) {

        ZERO(y_elem, $1_type)
      DOT_TESTGEN_NAME($1, $3, $2)(n_i, 0, n_i, norm, 
                                   blas_no_conj, &alpha_use, 1, 
                                   &beta_zero_fake, 1, x_vec, a_vec, seed, 
                                   PASS_BY_REF(y_elem, $1_type), 
                                   PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));
                        
      $2ge_commit_row(order, blas_no_trans, m_i, n_i, $5, $6, a_vec, i);

        /*commits an element to the truth */
      SET_VECTOR_ELEMENT(r_true, ri, r_true_elem, EXTRA_TYPE($1_type))
    }')dnl
dnl
dnl ******************************
dnl     DEGENERATE_RANDOM_TESTGEN
dnl
dnl  Purpose:  Handles the #2 case, randomize == 1, where one of alpha, beta is zero
dnl
dnl
dnl  arguments :        $1 = type of y, alpha, beta, 
dnl                     $2 = type of Matrix, 
dnl                     $3 = type of x, 
dnl                     $4 = name of nonzero-scalar (alpha or beta), 
dnl                     $5 = name of important matrix (a or b), 
dnl                     $6 = name of leading dimension of important matrix (lda/ldb) 
define(`DEGENERATE_RANDOM_TESTGEN',
`
        /* Fill in truth from $5, $4_i only */
      for (i = 0, ri = 0; 
          i < m_i; i++, ri += incri) {
        $2ge_copy_row(order, blas_no_trans, m_i, n_i, $5, $6, a_vec, i);

        ifelse(IS_MIXED($1, $2), `t', 
          `{
                  /* promote to complex */
           int r;
           for (r = 0; r < n_i; r++) {
           aa_vec[2*r] = a_vec[r];
             aa_vec[2*r+1] = 0.0;
           }
          }')
        ZERO(y_elem, $1_type)

        ifelse(IS_MIXED($1, $2, $3), `t', 
        `DOT_TESTGEN_NAME($1, $1, $1)(n_i, n_i, 0, norm, blas_no_conj, 
                       &alpha_use, 1, 
                       &beta_zero_fake, 1, 
                       ifelse(IS_MIXED($1, $3), `t', `xx_vec', `x_vec'), 
                       ifelse(IS_MIXED($1, $2), `t', `aa_vec', `a_vec'), 
                       seed, 
                       PASS_BY_REF(y_elem, $1_type), 
                       PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));',
        `DOT_TESTGEN_NAME($1, $3, $2)(n_i, n_i, 0, norm, blas_no_conj, 
                       &alpha_use, 1, 
                       &beta_zero_fake, 1, 
                       x_vec, a_vec, seed, 
                       PASS_BY_REF(y_elem, $1_type), 
                       PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));')

        SET_VECTOR_ELEMENT(r_true, ri, r_true_elem, EXTRA_TYPE($1_type))
      }')dnl
dnl
dnl
define(`REGULAR_GENERATION',
`
        /* Fill in matrix A, B */
    for(i = 0, ri = 0; i < m_i; i++, ri += incri) {
        ZERO(y_elem, $1_type)
      DOT_TESTGEN_NAME($1, $3, $2)(2*n_i, 0, 2*n_i, norm, 
                                   blas_no_conj, &alpha_use, 1, 
                                   &beta_zero_fake, 1, x_vec, a_vec, seed, 
                                   PASS_BY_REF(y_elem, $1_type), 
                                   PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));
                                
      $2ge_commit_row(order, blas_no_trans, m_i, n_i, a, lda, a_vec, i);
      $2ge_commit_row(order, blas_no_trans, m_i, n_i, b, ldb, (a_vec + inca_veci*n_i), i);

        /*commits an element to the truth */
      SET_VECTOR_ELEMENT(r_true, ri, r_true_elem, EXTRA_TYPE($1_type))
    }
')dnl
dnl
dnl
define(`PROTOTYPES', `FOREACH(`PREC_ARGS', `GE_SUM_MV_TESTGEN_HEAD(arg);
')')dnl
dnl
dnl
define(`SOURCE', `dnl
#include <stdio.h>
#include <stdlib.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"

#ifndef ALPHA_USE_IS_ALPHA
#define ALPHA_USE_IS_ALPHA 1
#define ALPHA_USE_IS_BETA 0
#define ALPHA_USE_IS_EITHER -1
#endif

FOREACH(`PREC_ARGS', `GE_SUM_MV_TESTGEN(arg)
')
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
