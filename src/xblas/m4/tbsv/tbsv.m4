dnl ---------------------------------------------------------
dnl  TBSV ---- Triangular Banded Solve
dnl
dnl    x  <---   alpha * T^-1 * x
dnl
dnl    where matrix T is a triangular band matrix.
dnl ---------------------------------------------------------
dnl
#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
include(cblas.m4)dnl
include(tbsv-common.m4)dnl
dnl
dnl
dnl  Usage:
dnl    TBSV        ($1, $2, $3)
dnl    TBSV_HEAD   ($1, $2, $3)
dnl    TBSV_NAME   ($1, $2, $3)
dnl    TBSV_PARAMS ($1, $2, $3)
dnl    TBSV_COMMENT($1, $2, $3)
dnl
dnl    $1 -- type of alpha, x.
dnl    $2 -- type of T matrix
dnl    $3 -- set to `_x' for _x routines.
dnl          Otherwise, `'.
dnl
define(`TBSV_COMMENT',`
/*
 * Purpose
 * =======
 * 
 * This routine solves :
 * 
 *     x <- alpha * inverse(t) * x
 * 
 * Arguments
 * =========
 * 
 * order  (input) enum blas_order_type
 *        column major, row major (blas_rowmajor, blas_colmajor)
 *
 * uplo   (input) enum blas_uplo_type
 *        upper, lower (blas_upper, blas_lower)
 *
 * trans  (input) enum blas_trans_type
 *        no trans, trans, conj trans
 * 
 * diag   (input) enum blas_diag_type
 *        unit, non unit (blas_unit_diag, blas_non_unit_diag)
 *
 * n      (input) int
 *        the dimension of t
 * 
 * k      (input) int
 *        the number of subdiagonals/superdiagonals of t
 *
 * alpha  (input) $1_scalar
 * 
 * t      (input) $2_array
 *        Triangular Banded matrix
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
dnl --------------------------------------------------------------------
dnl Usage: TBSV_INIT(ax_typeltr, T_typeltr, _x)
dnl        declare TBSV local variables, and check inputs
dnl --------------------------------------------------------------------
dnl
define(`TBSV_INIT',
` 
   int i, j; /* used to keep track of loop counts */
   int xi; /* used to `index' vector x */
   int start_xi; /* used as the starting idx to vector x */
   int incxi;
   int Tij; /* `index' inside of Banded structure */
   int dot_start, dot_start_inc1, dot_start_inc2,  dot_inc;

   PTR_CAST(t, $2_type, `const') /* internal matrix t */
   PTR_CAST(x, $1_type, `')      /* internal x */
   SCALAR_CAST(alpha, $1_type) /* internal alpha */

   if (order != blas_rowmajor && order != blas_colmajor) {
     BLAS_error(routine_name, -1, order, 0);
   }
   if (uplo != blas_upper && uplo != blas_lower) {
     BLAS_error(routine_name, -2, uplo, 0);
   }
   if ((trans != blas_trans) && (trans != blas_no_trans) &&
        (trans != blas_conj) && (trans != blas_conj_trans)) {
     BLAS_error(routine_name, -2, uplo, 0);
   }
   if (diag != blas_non_unit_diag && diag != blas_unit_diag) {
     BLAS_error(routine_name, -4, diag, 0);
   }
   if (n < 0) {
     BLAS_error(routine_name, -5, n, 0);
   }
   if (k >= n) {
     BLAS_error(routine_name, -6, k, 0);
   }
   if ((ldt < 1) || (ldt <= k)) {
     BLAS_error(routine_name, -9, ldt, 0);
   }
   if (incx == 0) {
     BLAS_error(routine_name, -11, incx, 0);
   }

   if (n <= 0) return;
        
   incxi = incx;
   INC_ADJUST(incxi, $1_type)

   /* configuring the vector starting idx */ 
   if (incxi < 0) { 
     start_xi = (1-n)*incxi; 
   } else {
     start_xi = 0;
   }

   /* if alpha is zero, then return x as a zero vector */
   if (TEST_0(alpha_i, $1_type)){
     xi = start_xi;
     for(i=0; i<n; i++){
       SET_ZERO_VECTOR_ELEMENT(x_i, xi, $1_type)
       xi += incxi;
     }
     return;
   }
   /* check to see if k=0.  if so, we can optimize somewhat */
   if (k == 0) {
     if((TEST_1(alpha_i, $1_type)) && (diag == blas_unit_diag)) {
       /* nothing to do */
       return;
     } else {
       /* just run the loops as is. */
       ifelse(`_x', $3,
         `/* must set prec to output. Ignore user input of prec */
          prec = ifelse(
            $1, `s', blas_prec_single,
            $1, `d', blas_prec_double,
            $1, `c', blas_prec_single,
            $1, `z', blas_prec_double);')
     }
   }
')dnl
dnl     
dnl
dnl ------------------------------------------------------------
dnl Usage : TBSV1_MAIN_BODY(x type, T type, internal type)
dnl
dnl   Makes decisions and call TBSV1_SUBST, 
dnl     with what conjugation.
dnl ------------------------------------------------------------
define(`TBSV1_MAIN_BODY', 
`
  {
    DECLARE(temp1, $3) /* temporary variable for calculations */ 
    DECLARE(temp2, $3) /* temporary variable for calculations */
    DECLARE(x_elem, $1) 
    DECLARE(T_element, $2) 
    ifelse($4, `FPU', `FPU_FIX_DECL;')


    ifelse($4, `FPU', `FPU_FIX_START;')
    IF_COMPLEX($2,
      `if ((trans == blas_conj) || (trans == blas_conj_trans)) {
         /* conjugated */
         TBSV1_SUBS($1, $2, $3, blas_conj)
       }  else {
        /* not conjugated */
        TBSV1_SUBS($1, $2, $3, no_conj)
       }',
      ` TBSV1_SUBS($1, $2, $3, all_real)')
    ifelse($4, `FPU', `FPU_FIX_STOP;')
  }')dnl
dnl
dnl
dnl ---------------------------------------------------------------------
dnl Usage: TBSV1_SUBS(x type, T type, internal type, 
dnl     $4 is blas_conj or no_conj)
dnl
dnl        perform substitution to compute
dnl                  x = alpha * inverse(T) * x
dnl
dnl
dnl ---------------------------------------------------------------------
define(`TBSV1_SUBS', `

        /*loop 1*/
        xi = start_xi;
       for(j=0; j<k; j++){

         /* each time through loop, xi lands on next x to compute.*/
         GET_VECTOR_ELEMENT(x_elem, x_i, xi, $1)
         /* preform the multiplication -
                in this implementation we do not seperate the alpha = 1 case */
         MUL(temp1, $3, x_elem, $1, alpha_i, $1)

         xi = start_xi; 

         Tij = dot_start;
         dot_start += dot_start_inc1;

         for (i=j; i>0; i--){
           GET_VECTOR_ELEMENT(T_element, t_i, Tij, $2)
           CONJ(T_element, $2, $4)
           GET_VECTOR_ELEMENT(x_elem, x_i, xi, $1)
           MUL(temp2, $3, x_elem, $1, T_element, $2)
           SUB(temp1, $3, temp1, $3, temp2, $3)
           xi += incxi;
           Tij += dot_inc;
         } /* for across row */
  

         /* if the diagonal entry is not equal to one, then divide Xj by 
            the entry */
         if (diag == blas_non_unit_diag){
           GET_VECTOR_ELEMENT(T_element, t_i, Tij, $2)
           CONJ(T_element, $2, $4)

           DIV(temp1, $3, temp1, $3, T_element, $2)

         } /* if (diag == blas_non_unit_diag) */

         SET_ROUND_VECTOR_ELEMENT(x_i, xi, temp1, $3)
         xi += incxi;
       } /* for j<k */
        /*end loop 1*/

        /*loop 2 continue without changing j to start*/
       for(; j<n; j++){

         /* each time through loop, xi lands on next x to compute.*/
         GET_VECTOR_ELEMENT(x_elem, x_i, xi, $1)
         MUL(temp1, $3, x_elem, $1, alpha_i, $1)

         xi = start_xi; 
         start_xi += incxi;

         Tij = dot_start;
         dot_start += dot_start_inc2;

         for (i = k; i>0; i--){
           GET_VECTOR_ELEMENT(T_element, t_i, Tij, $2)
           CONJ(T_element, $2, $4)
           GET_VECTOR_ELEMENT(x_elem, x_i, xi, $1)
           MUL(temp2, $3, x_elem, $1, T_element, $2)
           SUB(temp1, $3, temp1, $3, temp2, $3)
           xi += incxi;
           Tij += dot_inc;
         } /* for across row */
  

         /* if the diagonal entry is not equal to one, then divide by 
            the entry */
         if (diag == blas_non_unit_diag){
           GET_VECTOR_ELEMENT(T_element, t_i, Tij, $2)
           CONJ(T_element, $2, $4)

           DIV(temp1, $3, temp1, $3, T_element, $2)

         } /* if (diag == blas_non_unit_diag) */

         SET_ROUND_VECTOR_ELEMENT(x_i, xi, temp1, $3)
         xi += incxi;
       } /* for j<n */
')dnl
dnl
dnl
dnl
dnl
dnl ------------------------------------------------------------
dnl  Usage: TBSV(x type, t type, extra)
dnl   where
dnl    x,T types are one of { s,d,c,z }
dnl    extra can be either `_x' or ommited.
dnl
dnl    The types represent single, double, complex, or double complex
dnl     precision and type for the respective arguments.
dnl     Note that y type also specifies the type of alpha, beta
dnl    if extra is set to `_x', then the precision parameter, prec,
dnl     is added to the parameter list of the function. 
dnl  ------------------------------------------------------------
dnl 
define(`TBSV', 
  `TBSV_HEAD($1, $2, $3) 
   TBSV_COMMENT($1, $2, $3)
  {
  /* Routine name */
  static const char routine_name[] = "TBSV_NAME($1, $2, $3)";
  TBSV_INIT($1, $2, $3)
  /* get `index' variables prepared */
  INDEX_LOGIC($1, $2, $3)

    ifelse($3, `_x', `TBSV_X_BODY($1_type, $2_type, $3_type)', 
      `TBSV_SWITCH_BODY($1_type, $2_type, 
      HIGHER_TYPE($1_type, $2_type))')
  } /* end TBSV_NAME($1, $2, $3) */')dnl
dnl
dnl
dnl ------------------------------------------------------------
dnl Usage : TBSV_SWITCH_BODY(x type, T type, internal_precision, FPU)
dnl     Determines whether TBSV1 or TBSV2 is needed (TBSV2 mallocs space
dnl     for high precision internal variables)
dnl ------------------------------------------------------------
define(`TBSV_SWITCH_BODY',
`{
ifelse(`$1&&$1', `$1&&$3',
        `TBSV1_MAIN_BODY($1, $2, $3, $4)',
        `TBSV2_MAIN_BODY($1, $2, $3, $4)')
}')dnl
dnl
dnl ************************************************************
dnl  TBSV2 Stuff --
dnl     For cases fwhen we need to allocate space.
dnl
dnl ------------------------------------------------------------
dnl Usage : TBSV2_MAIN_BODY(x type, T type, internal type)
dnl
dnl   Makes decisions about which TBSV2_***_SUBS to call, 
dnl     and with what conjugation.
dnl ------------------------------------------------------------
define(`TBSV2_MAIN_BODY', 
`
  {
    DECLARE(temp1, $3) /* temporary variable for calculations */ 
    DECLARE(temp2, $3) /* temporary variable for calculations */
    DECLARE(temp3, $3) /* temporary variable for calculations */
    DECLARE(x_elem, $1) 
    DECLARE(T_element, $2) /* temporary variable for an element of matrix T */

    int x_inti=0, inc_x_inti=1;
    int k_compare=k; /*used for comparisons with x_inti*/
    DECLARE_VECTOR(x_internal, $3)
    ifelse($4, `FPU', `FPU_FIX_DECL;')

    INC_ADJUST(k_compare, $3)
    INC_ADJUST(inc_x_inti, $3) 
    MALLOC_VECTOR(x_internal, $3, `k')


    ifelse($4, `FPU', `FPU_FIX_START;')

    IF_COMPLEX($2, `
    if ((trans == blas_conj) || (trans == blas_conj_trans)) {
        /* conjugated */
        TBSV2_SUBS($1, $2, $3, blas_conj)
    }  else {
        /* not conjugated */
        TBSV2_SUBS($1, $2, $3, no_conj)
    }', 
        ` TBSV2_SUBS($1, $2, $3, all_real)')
    ifelse($4, `FPU', `FPU_FIX_STOP;')

    FREE_VECTOR(x_internal, $3)
  }')dnl
dnl
dnl
dnl ---------------------------------------------------------------------
dnl Usage: TBSV2_SUBS(x type, T type, internal type, 
dnl     $4 is blas_conj or no_conj)
dnl
dnl        perform substitution to compute
dnl                  x = alpha * inverse(T) * x
dnl
dnl     TBSV2 preforms internal variable allocation.  TBSV1 does not.
dnl
dnl ---------------------------------------------------------------------
define(`TBSV2_SUBS', 
  `/*loop 1*/
   xi = start_xi;
   /* x_inti already initialized to 0 */ 
   for(j=0; j<k; j++){

     /* each time through loop, xi lands on next x to compute.*/
     GET_VECTOR_ELEMENT(x_elem, x_i, xi, $1)
     /* preform the multiplication -
            in this implementation we do not seperate the alpha = 1 case */
     MUL(temp1, $3, x_elem, $1, alpha_i, $1)

     Tij = dot_start;
     dot_start += dot_start_inc1;

     /*start loop buffer over in loop 1*/
     x_inti = 0;
     for (i=j; i>0; i--){
       GET_VECTOR_ELEMENT(T_element, t_i, Tij, $2)
       CONJ(T_element, $2, $4)
       GET_ARRAY_ELEMENT(temp3, $3, x_internal, $3, x_inti)
       MUL(temp2, $3, temp3, $3, T_element, $2)
       SUB(temp1, $3, temp1, $3, temp2, $3)
       x_inti += inc_x_inti;
       Tij += dot_inc;
     } /* for across row */


     /* if the diagonal entry is not equal to one, then divide Xj by 
        the entry */
     if (diag == blas_non_unit_diag){
       GET_VECTOR_ELEMENT(T_element, t_i, Tij, $2)
       CONJ(T_element, $2, $4)

       DIV(temp1, $3, temp1, $3, T_element, $2)

     } /* if (diag == blas_non_unit_diag) */

     /* place internal precision result in internal buffer */
     SET_ARRAY_ELEMENT(temp1, $3, x_internal, $3,  x_inti)

     /* place result x in same place as got x this loop */
     SET_ROUND_VECTOR_ELEMENT(x_i, xi, temp1, $3)
     xi += incxi;
   } /* for j<k */
    /*end loop 1*/
    

    /* loop2 ******************************/
    x_inti = 0;
    /*loop 2 continue without changing j to start*/
   for(; j<n; j++){

     /* each time through loop, xi lands on next x to compute.*/
     GET_VECTOR_ELEMENT(x_elem, x_i, xi, $1)
     MUL(temp1, $3, x_elem, $1, alpha_i, $1)


     Tij = dot_start;
     dot_start += dot_start_inc2;

     for (i = k; i>0 && (x_inti < k_compare); i--){
       GET_VECTOR_ELEMENT(T_element, t_i, Tij, $2)
       CONJ(T_element, $2, $4)
       GET_ARRAY_ELEMENT(temp3, $3, x_internal, $3, x_inti)
       MUL(temp2, $3, temp3, $3, T_element, $2)
       SUB(temp1, $3, temp1, $3, temp2, $3)
       x_inti += inc_x_inti;
       Tij += dot_inc;
     } /* for across row */
     /*reset `index' to internal storage loop buffer.*/
     x_inti = 0; 
     for (; i>0; i--){
       GET_VECTOR_ELEMENT(T_element, t_i, Tij, $2)
       CONJ(T_element, $2, $4)
       GET_ARRAY_ELEMENT(temp3, $3, x_internal, $3, x_inti)
       MUL(temp2, $3, temp3, $3, T_element, $2)
       SUB(temp1, $3, temp1, $3, temp2, $3)
       x_inti += inc_x_inti;
       Tij += dot_inc;
     } /* for across row */


     /* if the diagonal entry is not equal to one, then divide by 
        the entry */
     if (diag == blas_non_unit_diag){
       GET_VECTOR_ELEMENT(T_element, t_i, Tij, $2)
       CONJ(T_element, $2, $4)

       DIV(temp1, $3, temp1, $3, T_element, $2)

     } /* if (diag == blas_non_unit_diag) */

     /* place internal precision result in internal buffer */
     SET_ARRAY_ELEMENT(temp1, $3, x_internal, $3, x_inti)
     x_inti += inc_x_inti; 
     if(x_inti >= k_compare)
            x_inti = 0;

     /* place result x in same place as got x this loop */
     SET_ROUND_VECTOR_ELEMENT(x_i, xi, temp1, $3)
     xi += incxi;
   } /* for j<n */
')dnl
dnl
dnl
dnl Usage: SWITCH_prec($1, $2, $3, $4, $5)
dnl Generates a 4-way switch statement based on prec.
dnl Arguments
dnl    $1      --  type of vector x
dnl    $2      --  type of matrix T
dnl    $3      --  type of internal variables in single case
dnl    $4      --  type of internal variables in double/indigenous case
dnl    $5      --  type of internal variables in extra case
dnl
define(`SWITCH_prec', 
  `switch (prec) {

    case blas_prec_single: ifelse(`$3', `$4', `', `{
      TBSV_SWITCH_BODY($1, $2, $3)
      break;
    }
    ')dnl
    case blas_prec_indigenous: 
    case blas_prec_double: {
      TBSV_SWITCH_BODY($1, $2, $4)
      break;
    }

    case blas_prec_extra: { 
      TBSV_SWITCH_BODY($1, $2, $5, FPU) 
      break;
    }

    default:
      BLAS_error(routine_name, -13, prec, 0);
      break;
  } /* end prec switch */')dnl
dnl
dnl 
dnl ------------------------------------------------------------ 
dnl  Usage : TBSV_X_BODY(x type, T type)
dnl ------------------------------------------------------------
dnl 
define(`TBSV_X_BODY', 
  `SWITCH_prec($1, $2, 
    TMP_TYPE_X($1, S), TMP_TYPE_X($1, D), TMP_TYPE_X($1, E))')dnl
dnl
dnl
dnl ------------------------------------------------------------
dnl  INDEX_LOGIC:
dnl     prepares indexing variables dot_start, dot_start_inc1,
dnl     dot_start_inc2, dot_inc, start_xi, incxi.
dnl ------------------------------------------------------------
define(`INDEX_LOGIC',
  `if( ((trans == blas_trans) || (trans == blas_conj_trans)) ^ 
       (order == blas_rowmajor)) {
     dot_start = k;
   } else {
     dot_start = 0;
   }

   if( ((trans == blas_trans) || (trans == blas_conj_trans)) ^
       (order == blas_rowmajor)) {
     dot_inc = 1;
     dot_start_inc1 = ldt -1;
     dot_start_inc2 = ldt;
   } else {
     dot_inc = ldt -1;
     dot_start_inc1 = 1;
     dot_start_inc2 = ldt;
   }
        
   if( ((trans == blas_trans) || (trans == blas_conj_trans)) ^
       (uplo == blas_lower)) {
     /*start at the first element of x */
     /* substitution will proceed forwards (forwardsubstitution)*/ 
   } else {
     /*start at the last element of x */
     /* substitution will proceed backwards (backsubstitution)*/ 
     dot_inc = -dot_inc;
     dot_start_inc1 = -dot_start_inc1;
     dot_start_inc2 = -dot_start_inc2;
     dot_start = ldt*(n-1) + k - dot_start;
     /*order of the following 2 statements matters! */
     start_xi = start_xi + (n-1)*incxi;
     incxi = -incxi;
   }

   INC_ADJUST(dot_inc, $2_type)
   INC_ADJUST(dot_start, $2_type)
   INC_ADJUST(dot_start_inc1, $2_type)
   INC_ADJUST(dot_start_inc2, $2_type)
')dnl
dnl
dnl
