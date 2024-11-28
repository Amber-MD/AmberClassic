dnl **********************************************************************
dnl Perform x = alpha * inverse(T) * x  
dnl **********************************************************************
dnl
#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
include(cblas.m4)dnl
include(trsv-common.m4)dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: TRSV_COMMENT(ax_typeltr, T_typeltr, _x)
dnl        ... generate the leading comment for the TRSV routine
dnl ----------------------------------------------------------------------
dnl
define(`TRSV_COMMENT',`
/*
 * Purpose
 * =======
 * 
 * This routine solve :
 * 
 *     x <- alpha * inverse(T) * x
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
 * x      (input) const $1_array
 *           Array of length n.
 * 
 * incx   (input) int
 *           The stride used to access components x[i].
 *
PREC_COMMENT($3)dnl
 */')dnl
dnl
dnl
dnl
dnl --------------------------------------------------------------------
dnl Usage: TRSV_INIT(ax_typeltr, T_typeltr, _x)
dnl        declare trsv local variables, and check inputs
dnl --------------------------------------------------------------------
dnl
define(`TRSV_INIT',
` int i, j; /* used to idx matrix */
   int ix, jx; /* used to idx vector x */
   int start_x; /* used as the starting idx to vector x */
   PTR_CAST(T, $2_type, `const') /* internal matrix T */
   PTR_CAST(x, $1_type, `')      /* internal x */
   SCALAR_CAST(alpha, $1_type) /* internal alpha */
   DECLARE(T_element, $2_type) /* temporary variable for an element of matrix A */
   int incT=1; /* internal ldt */

   if ((order != blas_rowmajor && order != blas_colmajor) ||
       (uplo != blas_upper && uplo != blas_lower) ||
       (trans != blas_trans && trans != 
        blas_no_trans && trans != blas_conj_trans) ||
       (diag != blas_non_unit_diag && diag != blas_unit_diag) ||
       (ldt < n) ||
       (incx == 0)){
     BLAS_error(routine_name, 0, 0, NULL);
   }
   
   if (n <= 0) return;

   INC_ADJUST(incT, $2_type)
   INC_ADJUST(incx, $1_type)
   /* configuring the vector starting idx */ 
   if (incx <= 0) { start_x = -(n-1)*incx; }
   else { start_x = 0; }

   /* if alpha is zero, then return x as a zero vector */
   if (TEST_0(alpha_i, $1_type)) {
     ix = start_x;
     for(i=0; i<n; i++) {
        SET_ZERO_VECTOR_ELEMENT(x_i, ix, $1_type)
        ix += incx;
     }
     return;
   }')dnl
dnl     
dnl
dnl --------------------------------------------------------------------
dnl Usage: TRSV1(ax_typeltr, T_typeltr, internal_prec )
dnl        ... x = alpha * transpose(T) * x
dnl        ... x = alpha * transpose(inverse(T)) * x
dnl
dnl        Each type and precision can be one of
dnl                 real_S ... real, single precision
dnl                 real_D ... real, double precision
dnl                 real_I ... real, indigenous precision
dnl                 real_E ... real, extra precision
dnl              complex_S ... complex, single precision
dnl              complex_D ... complex, double precision
dnl              complex_I ... complex, indigenous precision
dnl              complex_E ... complex, extra precision
dnl 
dnl       temp1, temp2:    temp scalar that holds an entry in array x.
dnl       adjust_row: adjusts the row index into matrix a
dnl       adjust_col: adjusts the coloumn index into matrix a
dnl --------------------------------------------------------------------
dnl
define(`TRSV1', `
  {
    DECLARE(temp1, $3) /* temporary variable for calculations */ 
    DECLARE(temp2, $3) /* temporary variable for calculations */
    DECLARE(temp3, $3) /* temporary variable for calculations */

    if ((order == blas_rowmajor && 
         trans == blas_no_trans && uplo == blas_upper) ||
        (order == blas_colmajor && 
         trans != blas_no_trans && uplo == blas_lower)) {
      TRSV1_CONJ_TRANS_BACKWARD($2, $3, $1, i, j)       
    } else if ((order == blas_rowmajor && 
                trans == blas_no_trans && uplo == blas_lower) ||
               (order == blas_colmajor && 
                trans != blas_no_trans && uplo == blas_upper)) {
      TRSV1_CONJ_TRANS_FORWARD($2, $3, $1, i, j)
    } else if ((order == blas_rowmajor && 
                trans != blas_no_trans && uplo == blas_lower) ||
               (order == blas_colmajor && 
                trans == blas_no_trans && uplo == blas_upper)) {
      TRSV1_CONJ_TRANS_BACKWARD($2, $3, $1, j, i)
    } else if ((order == blas_rowmajor && 
                trans != blas_no_trans && uplo == blas_upper) ||
               (order == blas_colmajor && 
                trans == blas_no_trans && uplo == blas_lower)) {
      TRSV1_CONJ_TRANS_FORWARD($2, $3, $1, j, i)
    }
  }')dnl
dnl
dnl
dnl
dnl --------------------------------------------------------------------
dnl Usage: TRSV2(ax_typeltr, T_typeltr, internal_prec )
dnl      
dnl        *************************************************************
dnl        NOTE: The only difference between TRSV1 and TRSV2 is that
dnl              TRSV2 dynamically allocates a vector of type 
dnl              internal_prec used as an intermediary for computing x,
dnl              so as to keep any intermediate computations in 
dnl              internal_prec accuracy. 
dnl        *************************************************************
dnl
dnl        ... x = alpha * transpose(T) * x
dnl        ... x = alpha * transpose(inverse(T)) * x
dnl
dnl        Each type and precision can be one of
dnl                 real_S ... real, single precision
dnl                 real_D ... real, double precision
dnl                 real_I ... real, indigenous precision
dnl                 real_E ... real, extra precision
dnl              complex_S ... complex, single precision
dnl              complex_D ... complex, double precision
dnl              complex_I ... complex, indigenous precision
dnl              complex_E ... complex, extra precision
dnl 
dnl       temp1, temp2:    temp scalar that holds an entry in array x.
dnl       adjust_row: adjusts the row index into matrix a
dnl       adjust_col: adjusts the coloumn index into matrix a
dnl --------------------------------------------------------------------
dnl
define(`TRSV2', `
  {
    int inc_intx; /* inc for intx */
    DECLARE(temp1, $3) /* temporary variable for calculations */ 
    DECLARE(temp2, $3) /* temporary variable for calculations */
    DECLARE(temp3, $3) /* temporary variable for calculations */
    DECLARE_VECTOR(intx, $3) /* copy of x used for calculations */

    /* allocate space for intx */
    MALLOC_VECTOR(intx, $3, n)  

    /* since intx is for internal usage, set it to 1 and then adjust
       it if necessary */
    inc_intx = 1;
    INC_ADJUST(inc_intx, $1)

    /* copy x to intx */
    ix=start_x; 
    jx=0;
    for (i=0; i<n; i++){
      GET_ARRAY_ELEMENT(temp1, $3, x_i, $1, ix)
      SET_ARRAY_ELEMENT(temp1, $3, intx, $3, jx)
      ix+=incx; 
      jx+=inc_intx;
    }

    if ((order == blas_rowmajor && 
         trans == blas_no_trans && uplo == blas_upper) ||
        (order == blas_colmajor && 
         trans != blas_no_trans && uplo == blas_lower)) {
      TRSV2_CONJ_TRANS_BACKWARD($2, $3, $1, i, j)       
    } else if ((order == blas_rowmajor && 
                trans == blas_no_trans && uplo == blas_lower) ||
               (order == blas_colmajor && 
                trans != blas_no_trans && uplo == blas_upper)) {
      TRSV2_CONJ_TRANS_FORWARD($2, $3, $1, i, j)
    } else if ((order == blas_rowmajor && 
                trans != blas_no_trans && uplo == blas_lower) ||
               (order == blas_colmajor && 
                trans == blas_no_trans && uplo == blas_upper)) {
      TRSV2_CONJ_TRANS_BACKWARD($2, $3, $1, j, i)
    } else if ((order == blas_rowmajor && 
                trans != blas_no_trans && uplo == blas_upper) ||
               (order == blas_colmajor && 
                trans == blas_no_trans && uplo == blas_lower)) {
      TRSV2_CONJ_TRANS_FORWARD($2, $3, $1, j, i)
    }
        
    /* copy the final results from intx to x */ 
    ix=start_x; 
    jx=0;
    for (i=0; i<n; i++){
      GET_ARRAY_ELEMENT(temp1, $3, intx, $3, jx)
      SET_ROUND_VECTOR_ELEMENT(x_i, ix, temp1, $3)
      ix+=incx; 
      jx+=inc_intx;
    }

    FREE_VECTOR(intx, $3)
  }')dnl
dnl
dnl
dnl -----------------------------------------------------------------------
dnl Usage: TRSV1_CONJ_TRANS_BACKWARD(T_typeltr, internal_prec, x_typeltr, col, row)
dnl        If the matrix type is complex, then allow blas_conj_trans option
dnl -----------------------------------------------------------------------
dnl
define(`TRSV1_CONJ_TRANS_BACKWARD',`ifelse(
  `$1', `complex_S', 
  `if (trans == blas_conj_trans) {
     TRSV1_BACKWARD_SUBST($1, $2, $3, $4, $5, blas_conj)
   } else {
     TRSV1_BACKWARD_SUBST($1, $2, $3, $4, $5, blas_no_conj)
   }', 
  `$1', `complex_D',
  `if (trans == blas_conj_trans) {
     TRSV1_BACKWARD_SUBST($1, $2, $3, $4, $5, blas_conj)
   } else {
     TRSV1_BACKWARD_SUBST($1, $2, $3, $4, $5, blas_no_conj)
   }', 
  `TRSV1_BACKWARD_SUBST($1, $2, $3, $4, $5, blas_no_conj)')')dnl
dnl
dnl
dnl -----------------------------------------------------------------------
dnl Usage: TRSV1_CONJ_TRANS_FORWARD(T_typeltr, internal_prec, x_typeltr, col, row)
dnl        If the matrix type is complex, then allow blas_conj_trans option
dnl -----------------------------------------------------------------------
dnl
define(`TRSV1_CONJ_TRANS_FORWARD',`ifelse(
  `$1', `complex_S', 
  `if (trans == blas_conj_trans) {
     TRSV1_FORWARD_SUBST($1, $2, $3, $4, $5, blas_conj)
   } else {
     TRSV1_FORWARD_SUBST($1, $2, $3, $4, $5, blas_no_conj)
   }', 
  `$1', `complex_D',
  `if (trans == blas_conj_trans) {
     TRSV1_FORWARD_SUBST($1, $2, $3, $4, $5, blas_conj)
   } else {
     TRSV1_FORWARD_SUBST($1, $2, $3, $4, $5, blas_no_conj)
   }', 
  `TRSV1_FORWARD_SUBST($1, $2, $3, $4, $5, blas_no_conj)')')dnl
dnl
dnl
dnl -----------------------------------------------------------------------
dnl Usage: TRSV2_CONJ_TRANS_BACKWARD(a_type, internal_prec, x_type, col, row)
dnl        If the matrix type is complex, then allow blas_conj_trans option
dnl -----------------------------------------------------------------------
dnl
define(`TRSV2_CONJ_TRANS_BACKWARD',`ifelse(
  `$1', `complex_S', 
  `if (trans == blas_conj_trans) {
     TRSV2_BACKWARD_SUBST($1, $2, $3, $4, $5, blas_conj)
   } else {
     TRSV2_BACKWARD_SUBST($1, $2, $3, $4, $5, blas_no_conj)
   }', 
  `$1', `complex_D',
  `if (trans == blas_conj_trans) {
     TRSV2_BACKWARD_SUBST($1, $2, $3, $4, $5, blas_conj)
   } else {
     TRSV2_BACKWARD_SUBST($1, $2, $3, $4, $5, blas_no_conj)
   }', 
  `TRSV2_BACKWARD_SUBST($1, $2, $3, $4, $5, blas_no_conj)')')dnl
dnl
dnl
dnl -----------------------------------------------------------------------
dnl Usage: TRSV2_CONJ_TRANS_FORWARD(a_type, internal_prec, x_type, col, row)
dnl        If the matrix type is complex, then allow blas_conj_trans option
dnl -----------------------------------------------------------------------
dnl
define(`TRSV2_CONJ_TRANS_FORWARD',`ifelse(
  `$1', `complex_S', 
  `if (trans == blas_conj_trans) {
     TRSV2_FORWARD_SUBST($1, $2, $3, $4, $5, blas_conj)
   } else {
     TRSV2_FORWARD_SUBST($1, $2, $3, $4, $5, blas_no_conj)
   }', 
  `$1', `complex_D',
  `if (trans == blas_conj_trans) {
     TRSV2_FORWARD_SUBST($1, $2, $3, $4, $5, blas_conj)
   } else {
     TRSV2_FORWARD_SUBST($1, $2, $3, $4, $5, blas_no_conj)
   }', 
  `TRSV2_FORWARD_SUBST($1, $2, $3, $4, $5, blas_no_conj)')')dnl
dnl
dnl
dnl
dnl ---------------------------------------------------------------------
dnl Usage: TRSV1_BACKWARD_SUBST(T_typeltr, internal_prec, x_typeltr, col, row, conj)
dnl        perform backward substitution to compute
dnl                  x = alpha * inverse(T) * x
dnl
dnl        col and row are used to index a matrix element.
dnl        The possible pairs of (col, row) values are (i, j) or (j, i)
dnl
dnl        TRSV1_BACKWARD_SUBST is used by a matrix of type:
dnl        1. blas_rowmajor, blas_no_trans, blas_upper             (i, j) 
dnl        2. blas_colmajor, blas_trans/blas_conj_trans, blas_lower(i, j)
dnl        3. blas_rowmajor, blas_trans/blas_conj_trans, blas_lower(j, i)
dnl        4. blas_colmajor, blas_no_trans, blas_upper             (j, i)
dnl ---------------------------------------------------------------------
define(`TRSV1_BACKWARD_SUBST', `
  jx = start_x+(n-1)*incx;
  for(j=n-1; j>=0; j--) {
   
  /* compute Xj = alpha*Xj - SUM Tij(or Tji) * Xi
                        i=j+1 to n-1           */  
  GET_ARRAY_ELEMENT(temp3, $2, x_i, $3, jx)
  MUL(temp1, $2, temp3, $2, alpha_i, $3)
  
  ix=start_x+(n-1)*incx;
  for (i=n-1; i>=j+1; i--) {
    GET_MATRIX_ELEMENT(T_element, T_i, $4*incT, $5, ldt*incT, $1)
    CONJ(T_element, $1, $6)
    GET_ARRAY_ELEMENT(temp3, $2, x_i, $3, ix)
    MUL(temp2, $2, temp3, $2, T_element, $1)
    SUB(temp1, $2, temp1, $2, temp2, $2)
    ix -= incx;
  } /* for j<n */

  /* if the diagonal entry is not equal to one, then divide Xj by 
     the entry */
  if (diag == blas_non_unit_diag) {
    GET_MATRIX_ELEMENT(T_element, T_i, j*incT, j, ldt*incT, $1)
    CONJ(T_element, $1, $6)

    DIV(temp1, $2, temp1, $2, T_element, $1)

  } /* if (diag == blas_non_unit_diag) */

  SET_ROUND_VECTOR_ELEMENT(x_i, jx, temp1, $2)

  jx -= incx;
} /* for j>=0 */')dnl
dnl
dnl
dnl
dnl ---------------------------------------------------------------------
dnl Usage: TRSV2_BACKWARD_SUBST(T_typeltr, internal_prec, x_typeltr, col, row, conj)
dnl        perform backward substitution to compute
dnl                  x = alpha * inverse(T) * x
dnl
dnl        col and row are used to index a matrix element.
dnl        The possible pairs of (col, row) values are (i, j) or (j, i)
dnl
dnl        TRSV2_BACKWARD_SUBST is used by a matrix of type:
dnl        1. blas_rowmajor, blas_no_trans, blas_upper             (i, j) 
dnl        2. blas_colmajor, blas_trans/blas_conj_trans, blas_lower(i, j)
dnl        3. blas_rowmajor, blas_trans/blas_conj_trans, blas_lower(j, i)
dnl        4. blas_colmajor, blas_no_trans, blas_upper             (j, i)
dnl ---------------------------------------------------------------------
define(`TRSV2_BACKWARD_SUBST', `
  jx = (n-1)*inc_intx;
  for(j=n-1; j>=0; j--) {
   
  /* compute Xj = alpha*Xj - SUM Aij(or Aji) * Xi
                       i=j+1 to n-1           */  
  GET_ARRAY_ELEMENT(temp3, $2, intx, $2, jx)
  /* multiply by alpha */
  MUL(temp1, $2, temp3, $2, alpha_i, $3)

  ix=(n-1)*inc_intx;
  for (i=n-1; i>=j+1; i--) {
    GET_MATRIX_ELEMENT(T_element, T_i, $4*incT, $5, ldt*incT, $1)
    CONJ(T_element, $1, $6)
    GET_ARRAY_ELEMENT(temp3, $2, intx, $2, ix)
    MUL(temp2, $2, temp3, $2, T_element, $1)
    SUB(temp1, $2, temp1, $2, temp2, $2)
    ix -= inc_intx;
  } /* for j<n */

  /* if the diagonal entry is not equal to one, then divide Xj by 
     the entry */
  if (diag == blas_non_unit_diag) {
    GET_MATRIX_ELEMENT(T_element, T_i, j*incT, j, ldt*incT, $1)
    CONJ(T_element, $1, $6)

    DIV(temp1, $2, temp1, $2, T_element, $1)

  } /* if (diag == blas_non_unit_diag) */

  SET_ARRAY_ELEMENT(temp1, $2, intx, $2, jx)

  jx -= inc_intx;
} /* for j>=0 */')dnl
dnl
dnl
dnl ---------------------------------------------------------------------
dnl Usage: TRSV1_FORWARD_SUBST(T_typeltr, internal_prec, x_typeltr, col, row)
dnl        perform forward substitution to compute
dnl                  x = alpha * inverse(T) * x
dnl
dnl        col and row are used to index a matrix element.
dnl        The possible pairs of (col, row) values are (i, j) or (j, i)
dnl
dnl        TRSV1_FORWARD_SUBST is used by a matrix of type:
dnl        1. blas_rowmajor, blas_no_trans, blas_lower             (i, j) 
dnl        2. blas_colmajor, blas_trans/blas_conj_trans, blas_upper(i, j)
dnl        3. blas_rowmajor, blas_trans/blas_conj_trans, blas_upper(j, i)
dnl        4. blas_colmajor, blas_no_trans, blas_lower             (j, i)
dnl ---------------------------------------------------------------------
define(`TRSV1_FORWARD_SUBST', `
  jx = start_x;
  for(j=0; j<n; j++){
   
  /* compute Xj = alpha*Xj - SUM Aij(or Aji) * Xi
                       i=j+1 to n-1           */  
  GET_ARRAY_ELEMENT(temp3, $2, x_i, $3, jx)
  /* multiply by alpha */
  MUL(temp1, $2, temp3, $2, alpha_i, $3)

  ix = start_x; 
  for (i=0; i<j; i++) {
    GET_MATRIX_ELEMENT(T_element, T_i, $4*incT, $5, ldt*incT, $1)
    CONJ(T_element, $1, $6)
    GET_ARRAY_ELEMENT(temp3, $2, x_i, $3, ix)
    MUL(temp2, $2, temp3, $2, T_element, $1)
    SUB(temp1, $2, temp1, $2, temp2, $2)
    ix += incx;
  } /* for i<j */
 
  /* if the diagonal entry is not equal to one, then divide Xj by 
     the entry */
  if (diag == blas_non_unit_diag) {
    GET_MATRIX_ELEMENT(T_element, T_i, j*incT, j, ldt*incT, $1)
    CONJ(T_element, $1, $6)

    DIV(temp1, $2, temp1, $2, T_element, $1)

  } /* if (diag == blas_non_unit_diag) */

  SET_ROUND_VECTOR_ELEMENT(x_i, jx, temp1, $2)
  jx += incx;
} /* for j<n */')dnl
dnl
dnl
dnl
dnl ---------------------------------------------------------------------
dnl Usage: TRSV2_FORWARD_SUBST(T_typeltr, internal_prec, ax_typeltr, col, row)
dnl        perform forward substitution to compute
dnl                  x = alpha * inverse(T) * x
dnl
dnl        col and row are used to index a matrix element.
dnl        The possible pairs of (col, row) values are (i, j) or (j, i)
dnl
dnl        TRSV2_FORWARD_SUBST is used by a matrix of type:
dnl        1. blas_rowmajor, blas_no_trans, blas_lower             (i, j) 
dnl        2. blas_colmajor, blas_trans/blas_conj_trans, blas_upper(i, j)
dnl        3. blas_rowmajor, blas_trans/blas_conj_trans, blas_upper(j, i)
dnl        4. blas_colmajor, blas_no_trans, blas_lower             (j, i)
dnl ---------------------------------------------------------------------
define(`TRSV2_FORWARD_SUBST', `
  jx = 0;
  for(j=0; j<n; j++){
   
    /* compute Xj = Xj - SUM Aij(or Aji) * Xi
                         i=j+1 to n-1           */  
    GET_ARRAY_ELEMENT(temp3, $2, intx, $2, jx)
    /* multiply by alpha */
    MUL(temp1, $2, temp3, $2, alpha_i, $3)

    ix = 0;
    for (i=0; i<j; i++) {
      GET_MATRIX_ELEMENT(T_element, T_i, $4*incT, $5, ldt*incT, $1)
      CONJ(T_element, $1, $6)
      GET_ARRAY_ELEMENT(temp3, $2, intx, $2, ix)
      MUL(temp2, $2, temp3, $2, T_element, $1)
      SUB(temp1, $2, temp1, $2, temp2, $2)
      ix += inc_intx;
    } /* for i<j */

    /* if the diagonal entry is not equal to one, then divide Xj by 
       the entry */
    if (diag == blas_non_unit_diag) {
      GET_MATRIX_ELEMENT(T_element, T_i, j*incT, j, ldt*incT, $1)
      CONJ(T_element, $1, $6)

      DIV(temp1, $2, temp1, $2, T_element, $1)
   
    } /* if (diag == blas_non_unit_diag) */

    SET_ARRAY_ELEMENT(temp1, $2, intx, $2, jx)
    jx += inc_intx;
  } /* for j<n */')dnl
dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: SWITCH_prec($1, $2, $3, $4, $5) ... generate
dnl        $3 is the type of 'temp1' and 'temp2' in single case.
dnl        $4 is the type of 'temp1' and 'temp2' in double/indigenous case.
dnl        $5 is the type of 'temp1' and 'temp2' in extra case.
dnl ----------------------------------------------------------------------
dnl 
define(`SWITCH_prec',
       `switch ( prec ) {
dnl
dnl $3 >= $1_type. Therefore,
dnl if ($1_type == $3) then do not need to allocate memory,
dnl so use TRSV1 instead of TRSV2
dnl
  case blas_prec_single: ifelse(`$3', `$4', `', `
    ifelse(`$1', `$3', `TRSV1($1, $2, $3)', `TRSV2($1, $2, $3)')
    break;
  ')dnl
dnl
dnl likewise if ($1 == $4)
dnl
  case blas_prec_double:
  case blas_prec_indigenous:
    ifelse(`$1', `$4', `TRSV1($1, $2, $4)', `TRSV2($1, $2, $4)')
    break;
  case blas_prec_extra:
dnl
dnl since $5 > $1, use TRSV2
dnl
  { FPU_FIX_DECL; FPU_FIX_START;
    { TRSV2($1, $2, $5) }
    FPU_FIX_STOP; }
  break;
  }')dnl
dnl
dnl
dnl
dnl --------------------------------------------------------------------
dnl Usage: TRSV_X_BODY(ax_prec, T_prec) ... dispatches
dnl        TRSV with appropriate type and precision info of
dnl        the specified internal prec.
dnl --------------------------------------------------------------------
dnl
define(`TRSV_X_BODY', 
  `SWITCH_prec($1, $2, 
    TMP_TYPE_X($1, S), TMP_TYPE_X($1, D), TMP_TYPE_X($1, E))')dnl
dnl
dnl
define(`TRSV_BODY', 
  `TRSV1($1, $2, $3)')dnl
dnl
dnl
define(`TRSV', `
  TRSV_HEAD($1, $2, $3)
  TRSV_COMMENT($1, $2, $3)
  {
    char *routine_name = "TRSV_NAME($1, $2)";

    TRSV_INIT($1, $2, $3)
    ifelse($3, _x, `TRSV_X_BODY($1_type, $2_type)', 
      `TRSV_BODY($1_type, $2_type, $1_type)')
  }
')dnl
dnl
dnl
