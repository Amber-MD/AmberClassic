dnl cblas.m4.h 
dnl
dnl
dnl -----------------------------------------------------------------------
dnl LIBNAMES
dnl List of all the routines currently implemented.
dnl -----------------------------------------------------------------------
define(`LIBNAMES', `dot, sum, axpby, waxpby, gemv, ge_sum_mv, gbmv, symv, 
  spmv, sbmv, hemv, hpmv, hbmv, trmv, tpmv, trsv, tbsv, gemm, symm, hemm,
  gemv2, symv2, hemv2, gbmv2')dnl
dnl -----------------------------------------------------------------------
dnl LIBNAMES_AMB
dnl List of all the routines needed by Amber.
dnl -----------------------------------------------------------------------
define(`LIBNAMES_AMB', `dot, axpby, waxpby, gemv, gemm, gemv2')dnl
dnl
dnl
dnl -----------------------------------------------------------------------
dnl FOREACH($1, $2) 
dnl Iterates through $1 and expands $2 for each dnl item in the list $1.
dnl Within $2, `arg' will be the unquoted argument from $1, while
dnl `argq' will be the quoted argument.
dnl -----------------------------------------------------------------------
dnl
define(`FOREACH', 
  `pushdef(`func', `$2')__FOREACH($1)popdef(`func')')dnl
dnl
dnl
define(`__FOREACH', `ifelse(`$1', `', `', 
  `pushdef(`argq', ``$1'')pushdef(`arg', `$1')func`'__FOREACH(shift($@))popdef(`arg')popdef(`argq')')')dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl if_blas(xxx) and if_noblas(xxx) expands to xxx depending whether
dnl no_plain_blas is defined.  Used to exclude plain BLAS from building.
dnl ----------------------------------------------------------------------
define(`if_blas', `ifdef(`no_plain_blas', `$2', `$1')')dnl
define(`if_noblas', `ifdef(`no_plain_blas', `$1', `$2')')dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl HEAD(x) and TAIL(x) gives access to leading and trailing 
dnl part of double-double type.
dnl ----------------------------------------------------------------------
define(`HEAD', `head_$1')dnl
define(`TAIL', `tail_$1')dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl REAL_COMPLEX(x) returns either `real' or `complex' depending
dnl on the type x.
dnl ----------------------------------------------------------------------
define(`REAL_COMPLEX', `ifelse(
  `$1', `real_S', `real',
  `$1', `real_D', `real',
  `$1', `real_E', `real',
  `complex')')dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl PREC(x) Returns one of S, D, or E depending on the precision 
dnl of the type x.
dnl ----------------------------------------------------------------------
define(`PREC', `ifelse(
  `$1', `real_S', `S',
  `$1', `complex_S', `S',
  `$1', `real_D', `D',
  `$1', `complex_D', `D',
  `E')')dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl HIGHER_PREC(x, y) Returns the higher of the two precisions 
dnl x and y.  The input x and y must be one of S, D, or E.
dnl ----------------------------------------------------------------------
define(`HIGHER_PREC', `ifelse(
  `$1', `S', `$2', 
  `$2', `S', `$1', 
  `$1', `D', `$2', 
  `$2', `D', `$1', 
  `E')')dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl HIGHER_REAL_COMPLEX(x, y) Returns complex if at least one of the two
dnl types is complex, otherwise returning real.  The input x and y must
dnl be one of `real' or `complex'.
dnl ----------------------------------------------------------------------
define(`HIGHER_REAL_COMPLEX', `ifelse(
  `$1&&$2', `real&&real', `real', 
  `complex')')dnl
dnl
dnl
define(`HIGHER_TYPE2', `HIGHER_REAL_COMPLEX(REAL_COMPLEX($1), 
  REAL_COMPLEX($2))`_'HIGHER_PREC(PREC($1), PREC($2))')dnl
dnl
dnl
dnl HIGHER_TYPE(x, y[, z]) returns the higher type of x, y, and z 
dnl (if z exists).
define(`HIGHER_TYPE', `ifelse(`$3', `', 
  `HIGHER_TYPE2($1, $2)', `HIGHER_TYPE2(HIGHER_TYPE2($1, $2), $3)')')dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl SUM_TYPE(x, y, z) For dot product like routines (x = y dot z), this 
dnl gives the type to be used for accumulating the sum.
dnl ----------------------------------------------------------------------
define(`SUM_TYPE', `HIGHER_REAL_COMPLEX(REAL_COMPLEX($2), REAL_COMPLEX($3))`_'PREC($1)')dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl SUM_TYPE_X(x, y, z, w) For dot product like routines (x = y dot z), 
dnl this gives the type to be used for accumulating the sum when
dnl using extra precision specified by w.  Input w must be one of S, 
dnl D, or E.
dnl ----------------------------------------------------------------------
define(`SUM_TYPE_X', `HIGHER_REAL_COMPLEX(REAL_COMPLEX($2), REAL_COMPLEX($3))`_'HIGHER_PREC(PREC($1), $4)')dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl TMP_TYPE_X(x, y) This gives the type of the temporary variable
dnl when using extra precision.  Input x is the output type, while 
dnl input y is the specified extra precision.
dnl ----------------------------------------------------------------------
define(`TMP_TYPE_X', `REAL_COMPLEX($1)`_'HIGHER_PREC(PREC($1), $2)')dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl PREC_COMMENT(a) If a is set to `_x', this produces the comment
dnl for prec parameter.  Otherwise it produces nothing.
dnl ----------------------------------------------------------------------
define(`PREC_COMMENT', `ifelse(
  `$1', `_x',
` * prec   (input) enum blas_prec_type
 *        Specifies the internal precision to be used.
 *        = blas_prec_single: single precision.
 *        = blas_prec_double: double precision.
 *        = blas_prec_extra : anything at least 1.5 times as accurate
 *                            than double, and wider than 80-bits.
 *                            We use double-double in our implementation.
 *
')')dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl IS_EXTRA returns `t' is the given type is a extra type, `f' 
dnl otherwise.  Use to detect whether to turn on IEEE double precision
dnl rounding for intel x86 machines.
dnl ----------------------------------------------------------------------
define(`IS_EXTRA', `ifelse(PREC($1), `E', `t', `f')')dnl
dnl
dnl
dnl IF_COMPLEX(type, expr_1, expr_2) expands to expr_1 if type is a complex 
dnl type, otherwise it expands to expr_2.
define(`IF_COMPLEX', `ifelse(REAL_COMPLEX($1), `complex', `$2', `$3')')dnl
dnl
dnl
dnl IF_COMPLEX(type, expr_1, expr_2) expands to expr_1 if type is a real 
dnl type, otherwise it expands to expr_2.
define(`IF_REAL', `ifelse(REAL_COMPLEX($1), `real', `$2', `$3')')dnl
dnl
dnl
dnl DOUBLE_TYPE(type) returns appropriate double precision type (real or complex).
define(`DOUBLE_TYPE', `REAL_COMPLEX($1)_D')dnl
dnl
dnl
dnl REAL_TYPE(type) returns real version of the given type.
define(`REAL_TYPE', `real_`'PREC($1)')dnl
dnl
dnl
dnl COMPLEX_TYPE(type) returns complex version of the given type.
define(`COMPLEX_TYPE', `complex_`'PREC($1)')dnl
dnl
dnl
dnl EXTRA_TYPE(type) returns appropriate double-double type (real or complex).
define(`EXTRA_TYPE', `REAL_COMPLEX($1)_E')dnl
dnl
dnl
dnl DOUBLE_ABBREV(abbrev) returns appropriate abbreviated type form for double 
dnl precision version of the given type.
define(`DOUBLE_ABBREV', 
  `ifelse(`$1', `s', `d', 
          `$1', `d', `d', 
          `$1', `c', `z', 
          `$1', `z', `z')')dnl
dnl
dnl
dnl REAL_ABBREV(abbrev) returns appropriate abbreviated type form for real
dnl version of the given type.
define(`REAL_ABBREV', 
  `ifelse(`$1', `c', `s', 
          `$1', `z', `d', 
          `$1')')dnl
dnl
dnl
dnl COMPLEX_ABBREV(abbrev) returns appropriate abbreviated type form for complex
dnl version of the given type.
define(`COMPLEX_ABBREV', 
  `ifelse(`$1', `s', `c', 
          `$1', `d', `z', 
          `$1')')dnl
dnl
dnl
dnl ABBREV(type) returns the appropriate abbreviated form for the given type.
define(`ABBREV', `ifelse(
  `$1', `real_S', `s', 
  `$1', `real_D', `d', 
  `$1', `complex_S', `c', 
  `$1', `complex_D', `z', `
#error Unsupported type for ABB`'REV.')')dnl
dnl
dnl
dnl BLAS_PREC(type) returns appropriate blas_prec_type enum
define(`BLAS_PREC', `ifelse(
  PREC($1), `S', `blas_prec_single', 
  PREC($1), `D', `blas_prec_double', 
  PREC($1), `E', `blas_prec_extra', 
  `blas_prec_unknown')')dnl
dnl
dnl
dnl NEGATE(var, type) negates the variable var
define(`NEGATE', `ifelse(
  IS_EXTRA($2), `t', `
#error negation of double-double unsupported.', 
  REAL_COMPLEX($2), `real', `$1 = -$1;', 
  `$1[0] = -$1[0]; $1[1] = -$1[1];')')dnl
dnl
dnl
dnl Usage: COPY_VECTOR_ELEMENT(x, i, y, j, type)  ... x[i] = y[j]
define(`COPY_VECTOR_ELEMENT', `ifelse(
  IS_EXTRA($5), `t', 
    `ifelse(REAL_COMPLEX($5), `real', 
      `HEAD($1)[$2] = HEAD($3)[$4]; TAIL($1)[$2] = TAIL($3)[$4];', 
      `HEAD($1)[$2] = HEAD($3)[$4]; TAIL($1)[$2] = TAIL($3)[$4];
       HEAD($1)[$2+1] = HEAD($3)[$4+1]; TAIL($1)[$2+1] = TAIL($3)[$4+1];')', 
  REAL_COMPLEX($5), `real', `$1[$2] = $3[$4];', 
  `$1[$2] = $3[$4]; $1[$2+1] = $3[$4+1];')')dnl
dnl
dnl
dnl VECTOR_ELEMENT(array, index, type) takes array[index] and massages it 
dnl appropriate for passing to another function.
define(`VECTOR_ELEMENT', `ifelse(
  IS_EXTRA($3), `t', 
    `ifelse(REAL_COMPLEX($3), `real', `HEAD($1)[$2], TAIL($1)[$2]', 
      `&HEAD($1)[$2], &TAIL($1)[$2]')', 
  REAL_COMPLEX($3), `real', `$1[$2]', `&$1[$2]')')dnl
dnl
dnl
dnl PASS_BY_REF(var, type) takes the variable and massages it so that it
dnl is appropriate for passing them by reference.
define(`PASS_BY_REF', `ifelse(
  IS_EXTRA($2), `t', 
    `ifelse(REAL_COMPLEX($2), `real', `&HEAD($1), &TAIL($1)', `HEAD($1), TAIL($1)')', 
  REAL_COMPLEX($2), `real', `&$1', `$1')')dnl
dnl
dnl
dnl GET_REAL_PART(a, a_type, b, b_type) ... set a = real(b)
define(`GET_REAL_PART', `ifelse(
  `$2&&$4', `real_S&&complex_S', `$1 = $3[0];', 
  `$2&&$4', `real_D&&complex_D', `$1 = $3[0];', `
#error GET_REAL_PART: invalid argument.')')dnl
dnl
dnl
dnl GET_IMAG_PART(a, a_type, b, b_type) ... set a = imag(b)
define(`GET_IMAG_PART', `ifelse(
  `$2&&$4',`real_S&&complex_S', `$1 = $3[1];',
  `$2&&$4',`real_D&&complex_D', `$1 = $3[1];', `
#error GET_IMAG_PART: invalid argument.')')dnl
dnl
dnl
dnl PUT_REAL_PART(a, a_type, b, b_type) ... set real(a) = b
define(`PUT_REAL_PART', `ifelse(
  `$2&&$4',`complex_S&&real_S', `$1[0] = $3;',
  `$2&&$4',`complex_D&&real_D', `$1[0] = $3;', `
#error PUT_REAL_PART: invalid argument.')')dnl
dnl
dnl
dnl PUT_IMAG_PART(a, a_type, b, b_type) ... set imag(a) = b
define(`PUT_IMAG_PART', `ifelse(
  `$2&&$4',`complex_S&&real_S', `$1[1] = $3;',
  `$2&&$4',`complex_D&&real_D', `$1[1] = $3;', `
#error PUT_IMAG_PART: invalid argument.')')dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl [typeltr]_scalar returns the C type that corresponds to scalars of
dnl type typeltr in function invocations (so s_scalar becomes float,
dnl z_scalar becomes void*, etc.).
dnl ----------------------------------------------------------------------
define(`s_scalar', `float')dnl
define(`d_scalar', `double')dnl
define(`c_scalar', `const void*')dnl
define(`z_scalar', `const void*')dnl
dnl
dnl ----------------------------------------------------------------------
dnl [typeltr]_array returns the C type that corresponds to the arrays of
dnl type typeltr in function invocations (so s_array becomes float*,
dnl z_array becomes void*, etc.).
dnl ----------------------------------------------------------------------
define(`s_array',  `float*')dnl
define(`d_array',  `double*')dnl
define(`c_array',  `void*')dnl
define(`z_array',  `void*')dnl
dnl
dnl ----------------------------------------------------------------------
dnl [typeltr]_type returns the type and precision that correspond to the
dnl typeltr (so s_type becomes real_S, z_type becomes complex_D, etc.).
dnl ----------------------------------------------------------------------
define(`s_type', `real_S')dnl
define(`d_type', `real_D')dnl
define(`c_type', `complex_S')dnl
define(`z_type', `complex_D')dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: PTR_CAST(var, var_type, `const') ... cast the input pointer
dnl        argument into an internal pointer of the proper type.
dnl        The third argument is a string and takes either the values
dnl        `const' or `'.
dnl ----------------------------------------------------------------------
define(`PTR_CAST', `ifelse(
	`$2',`real_S', `$3 float *$1_i = $1;',
	`$2',`real_D', `$3 double *$1_i = $1;',
	`$2',`complex_S', `$3 float *$1_i = (float*) $1;',
	`$2',`complex_D', `$3 double *$1_i = (double*) $1;',
      )')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: SCALAR_CAST(var, var_type) ... cast the input scalar
dnl        argument into an internal representation.
dnl        The internal variable has the name var_i.
dnl ----------------------------------------------------------------------
define(`SCALAR_CAST', `ifelse(
	`$2',`real_S', `float $1_i = $1;',
	`$2',`real_D', `double $1_i = $1;',
	`$2',`complex_S', `float *$1_i = (float*) $1;',
	`$2',`complex_D', `double *$1_i = (double*) $1;',
      )')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: R_TRUE_TYPE(abr_typeltr)
dnl        ... generate proper type specifier for TAIL(r)rue.
dnl ----------------------------------------------------------------------
define(`R_TRUE_TYPE', `ifelse(
	`$1', `s', `double', `$1', `d', `double',
	`double*')')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: DECLARE(var, var_type) ... declare a variable
dnl ----------------------------------------------------------------------
define(`DECLARE', `ifelse(
	`$2',`real_S', `float $1;',
	`$2',`real_D', `double $1;',
	`$2',`real_E', `double HEAD($1), TAIL($1);',
	`$2',`complex_S', `float $1[2];',
	`$2',`complex_D', `double $1[2];',
	`$2',`complex_E', `double HEAD($1)[2], TAIL($1)[2];',
      )')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: ONE(a, a_type) ... a = 1
dnl ----------------------------------------------------------------------
define(`ONE', `ifelse(
	`$2',`real_E', `HEAD($1) = 1.0; TAIL($1) = 0.0;',
	`$2',`complex_S', `$1[0] = 1.0; $1[1] = 0.0;',
	`$2',`complex_D', `$1[0] = 1.0; $1[1] = 0.0;',
	`$2',`complex_E', `HEAD($1)[0] = 1.0; HEAD($1)[1] = TAIL($1)[0] = TAIL($1)[1] = 0.0;',
        `$1 = 1.0;')')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: ZERO(a, a_type) ... a = 0
dnl ----------------------------------------------------------------------
define(`ZERO', `ifelse(
	`$2',`real_E', `HEAD($1) = TAIL($1) = 0.0;',
	`$2',`complex_S', `$1[0] = $1[1] = 0.0;',
	`$2',`complex_D', `$1[0] = $1[1] = 0.0;',
	`$2',`complex_E', `HEAD($1)[0] = HEAD($1)[1] = TAIL($1)[0] = TAIL($1)[1] = 0.0;',
        `$1 = 0.0;')')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: ZERO_OUT(a, a_type) ... *a = 0
dnl ----------------------------------------------------------------------
define(`ZERO_OUT', `ifelse(
	`$2',`complex_S', `$1[0] = $1[1] = 0.0;',
	`$2',`complex_D', `$1[0] = $1[1] = 0.0;',
	`*$1 = 0.0;')')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: ASSIGN(a, a_type, b, b_type) ... set a = b
dnl        a must have at least as high a type and precision as b.
dnl ----------------------------------------------------------------------
define(`ASSIGN', `ifelse(
        `$2&&$4',`real_E&&real_E', `HEAD($1) = HEAD($3); TAIL($1) = TAIL($3);',
	`$2',`real_E', `HEAD($1) = $3; TAIL($1) = 0.0;',
	`$2&&$4',`complex_S&&complex_S', `$1[0] = $3[0]; $1[1] = $3[1];',
	`$2&&$4',`complex_D&&complex_D', `$1[0] = $3[0]; $1[1] = $3[1];',
	`$2&&$4',`complex_D&&complex_S', `$1[0] = $3[0]; $1[1] = $3[1];',
	`$2&&$4',`complex_S&&real_S', `$1[0] = $3; $1[1] = 0.0;',
	`$2&&$4',`complex_D&&real_D', `$1[0] = $3; $1[1] = 0.0;',
	`$2&&$4',`complex_D&&real_S', `$1[0] = $3; $1[1] = 0.0;',
	`$2&&$4',`complex_E&&complex_E',
	    `HEAD($1)[0] = HEAD($3)[0]; TAIL($1)[0] = TAIL($3)[0];
	    HEAD($1)[1] = HEAD($3)[1]; TAIL($1)[1] = TAIL($3)[1];',
	`$2&&$4',`complex_E&&complex_D',
            `HEAD($1)[0] = $3[0]; TAIL($1)[0] = 0.0;
            HEAD($1)[1] = $3[1]; TAIL($1)[1] = 0.0;',
	`$2&&$4',`complex_E&&complex_S',
            `HEAD($1)[0] = $3[0]; TAIL($1)[0] = 0.0; HEAD($1)[1] = $3[1]; TAIL($1)[1] = 0.0;',
	`$2&&$4',`complex_E&&real_E',
            `HEAD($1)[0] = HEAD($3); TAIL($1)[0] = TAIL($3); HEAD($1)[1] = TAIL($1)[1] = 0.0;',
	`$2&&$4',`complex_E&&real_D',
	    `HEAD($1)[0] = $3; TAIL($1)[0] = HEAD($1)[1] = TAIL($1)[1] = 0.0;',
	`$2&&$4',`complex_E&&real_S',
            `HEAD($1)[0] = $3; TAIL($1)[0] = HEAD($1)[1] = TAIL($1)[1] = 0.0;',
	`$1 = $3;')')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: ASSIGN_PTR_TO_SCALAR(a, a_type, b, b_type) ... set a = *b
dnl ----------------------------------------------------------------------
define(`ASSIGN_PTR_TO_SCALAR', `ifelse(
  `$4',`complex_S',`$1[0] = ((float*)$3)[0]; $1[1] = ((float*)$3)[1];',
  `$4',`complex_D',`$1[0] = ((double*)$3)[0]; $1[1] = ((double*)$3)[1];',
  `$1 = *$3;')')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: ASSIGN_SCALAR_TO_PTR(a, a_type, b, b_type) ... set *a = b
dnl ----------------------------------------------------------------------
define(`ASSIGN_SCALAR_TO_PTR', `ifelse(
  `$2',`complex_S',`((float*)$1)[0] = $3[0]; ((float*)$1)[1] = $3[1];',
  `$2',`complex_D',`((double*)$1)[0] = $3[0]; ((double*)$1)[1] = $3[1];',
  `$2&&$4',`real_E&&real_E',`*HEAD($1) = HEAD($3); *TAIL($1) = TAIL($3);',
  `$2&&$4',`complex_E&&complex_E',
        `HEAD($1)[0] = HEAD($3)[0]; TAIL($1)[0] = TAIL($3)[0];
         HEAD($1)[1] = HEAD($3)[1]; TAIL($1)[1] = TAIL($3)[1];',
  `*$1 = $3;')')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: ROUND(a, a_type, b, b_type) ... *a = (round)b
dnl        *a must have the same type as b, but may have any precision
dnl        but extra.
dnl ----------------------------------------------------------------------
define(`ROUND', `ifelse( 
	`$2&&$4',`complex_D&&complex_E',
            `((double*)$1)[0] = HEAD($3)[0]; ((double*)$1)[1] = HEAD($3)[1];',
	`$2&&$4',`complex_S&&complex_E',
            `((float*)$1)[0] = HEAD($3)[0]; ((float*)$1)[1] = HEAD($3)[1];',
	`$2',`complex_S',`((float*)$1)[0] = $3[0]; ((float*)$1)[1] = $3[1];',
	`$2',`complex_D',`((double*)$1)[0] = $3[0]; ((double*)$1)[1] = $3[1];',
	`$4',`real_E', `*$1 = HEAD($3);',
	`*$1 = $3;')')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: INC_ADJUST(inc, x_type) ... modify incx to correcly deal with
dnl        cases where x is complex.
dnl ----------------------------------------------------------------------
define(`INC_ADJUST', `ifelse(
        `$2', `complex_S', `$1 *= 2;',
        `$2', `complex_D', `$1 *= 2;',
	`$2', `complex_E', `$1 *= 2;',)')dnl
dnl
dnl --------------------------------------------------------------------
dnl Usage: CONJ(a, a_type, conj) ... conjugates a if
dnl        conj == blas_conj; otherwise does nothing.
dnl --------------------------------------------------------------------
define(`CONJ', `ifelse(
	`$3', `blas_conj', `CONJ_AUX($1, $2)', `')')dnl
dnl
dnl --------------------------------------------------------------------
dnl Usage: CONJ_AUX(a, a_type) ... a = conjugate of a
dnl --------------------------------------------------------------------
define(`CONJ_AUX', `ifelse(
	`$2', `complex_S', `$1[1] = - $1[1];',
	`$2', `complex_D', `$1[1] = - $1[1];',
	`$2', `complex_I', `$1[1] = - $1[1];',
	`$2', `complex_E', `HEAD($1)[1] = - HEAD($1)[1]; TAIL($1)[1] = -TAIL($1)[1];'
	`')')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: GET_VECTOR_ELEMENT(xi, x, i, type)   ... xi = x[i]
dnl ----------------------------------------------------------------------
define(`GET_VECTOR_ELEMENT', `ifelse( 
        `$4', `complex_S', `$1[0] = $2[$3]; $1[1] = $2[$3+1];',
        `$4', `complex_D', `$1[0] = $2[$3]; $1[1] = $2[$3+1];',
        `$4', `real_E', `HEAD($1) = HEAD($2)[$3]; TAIL($1) = TAIL($2)[$3];',
        `$4', `complex_E', 
          `HEAD($1)[0] = HEAD($2)[$3]; HEAD($1)[1] = HEAD($2)[$3+1];
           TAIL($1)[0] = TAIL($2)[$3]; TAIL($1)[1] = TAIL($2)[$3+1];',
        `$1 = $2[$3];')')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: GET_REAL_VECTOR_ELEMENT(xi, x, i, type)   ... xi = real(x[i])
dnl ----------------------------------------------------------------------
define(`GET_REAL_VECTOR_ELEMENT', `ifelse( 
        `$4', `complex_S', `$1 = $2[$3];',
        `$4', `complex_D', `$1 = $2[$3];',
        `$4', `real_E', `HEAD($1) = HEAD($2)[$3]; TAIL($1) = TAIL($2)[$3];',
        `$4', `complex_E', `HEAD($1) = HEAD($2)[$3]; TAIL($1) = TAIL($2)[$3];',
        `$1 = $2[$3];')')dnl
dnl
dnl ---------------------------------------------------------------------
dnl Usage: GET_ARRAY_ELEMENT(xi, xi_type, x, x_type, i) ... xi = x[i]
dnl        xi_type must be higher than x_type.
dnl----------------------------------------------------------------------
define(`GET_ARRAY_ELEMENT', `ifelse(
	`$2&&$4', `real_E&&real_E', `HEAD($1) = HEAD($3)[$5]; TAIL($1) = TAIL($3)[$5];',
	`$2',`real_E', `HEAD($1) = $3[$5]; TAIL($1) = 0.0;',
	`$2&&$4',`complex_S&&complex_S', 
	    `$1[0] = $3[$5];  $1[1] = $3[1+$5];',
	`$2&&$4',`complex_D&&complex_D', 
 	    `$1[0] = $3[$5];  $1[1] = $3[1+$5];',
	`$2&&$4',`complex_D&&complex_S', 
	    `$1[0] = $3[$5];  $1[1] = $3[1+$5];',
	`$2&&$4',`complex_S&&real_S', 
	    `$1[0] = $3[$5];  $1[1] = 0.0;',
	`$2&&$4',`complex_D&&real_D', 
  	    `$1[0] = $3[$5];  $1[1] = 0.0;',
	`$2&&$4',`complex_D&&real_S', 
	    `$1[0] = $3[$5];  $1[1] = 0.0;',
	`$2&&$4',`complex_E&&complex_E',
	     `HEAD($1)[0] = HEAD($3)[$5]; HEAD($1)[1] = HEAD($3)[1+$5];
	     TAIL($1)[0] = TAIL($3)[$5]; TAIL($1)[1] = TAIL($3)[1+$5];',
	`$2&&$4',`complex_E&&complex_D', 
	     `HEAD($1)[0] = $3[$5]; TAIL($1)[0] = 0.0;
	     HEAD($1)[1] = $3[1+$5]; TAIL($1)[1] = 0.0;',
	`$2&&$4',`complex_E&&complex_S',
             `HEAD($1)[0] = $3[$5]; TAIL($1)[0] = 0.0;
	     HEAD($1)[1] = $3[1+$5]; TAIL($1)[1] = 0.0;',
	`$2&&$4',`complex_E&&real_E',
            `HEAD($1)[0] = HEAD($3)[$5]; TAIL($1)[0] = TAIL($3)[$5]; HEAD($1)[1] = TAIL($1)[1] = 0.0;',
	`$2&&$4',`complex_E&&real_D',
	    `HEAD($1)[0] = $3[$5]; TAIL($1)[0] = HEAD($1)[1] = TAIL($1)[1] = 0.0;',
	`$2&&$4',`complex_E&&real_S',
            `HEAD($1)[0] = $3[$5]; TAIL($1)[0] = HEAD($1)[1] = TAIL($1)[1] = 0.0;',
	`$1 = $3[$5];')')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: SET_VECTOR_ELEMENT(x, i, val, type)  ... x[i] = val
dnl ----------------------------------------------------------------------
define(`SET_VECTOR_ELEMENT', `ifelse(
  IS_EXTRA($4), `t', 
    `ifelse(REAL_COMPLEX($4), `real', 
      `HEAD($1)[$2] = HEAD($3); TAIL($1)[$2] = TAIL($3);', 
      `HEAD($1)[$2] = HEAD($3)[0]; HEAD($1)[$2+1] = HEAD($3)[1];
       TAIL($1)[$2] = TAIL($3)[0]; TAIL($1)[$2+1] = TAIL($3)[1];')', 
  REAL_COMPLEX($4), `real', `$1[$2] = $3;', 
  `$1[$2] = $3[0]; $1[$2+1] = $3[1];')')dnl
dnl
dnl
dnl ---------------------------------------------------------------------
dnl Usage: SET_ARRAY_ELEMENT(xi, xi_type, x, x_type, i) ... x[i] = xi
dnl        xi_type must be higher than x_type.
dnl----------------------------------------------------------------------
define(`SET_ARRAY_ELEMENT', `ifelse(
	`$2&&$4', `real_E&&real_E', `HEAD($3)[$5] = HEAD($1); TAIL($3)[$5] = TAIL($1);',
	`$2',`real_E', `$3[$5] = HEAD($1);',
	`$2&&$4',`complex_S&&complex_S', 
	    `$3[$5] = $1[0]; $3[1+$5] = $1[1];',
	`$2&&$4',`complex_D&&complex_D', 
	    `$3[$5] = $1[0]; $3[1+$5] = $1[1];',
	`$2&&$4',`complex_D&&complex_S', 
	    `$3[$5] = $1[0]; $3[1+$5] = $1[1];',
	`$2&&$4',`complex_E&&complex_E',
	    `HEAD($3)[$5] = HEAD($1)[0]; TAIL($3)[$5] = TAIL($1)[0];
	    HEAD($3)[1+$5] = HEAD($1)[1]; TAIL($3)[1+$5] = TAIL($1)[1];',
	`$2&&$4',`complex_E&&complex_D',
            `$3[$5] = HEAD($1)[0]; $3[1+$5] = HEAD($1)[1];',
	`$2&&$4',`complex_E&&complex_S',
            `$3[$5] = HEAD($1)[0]; $3[1+$5] = HEAD($1)[1];',
	`$3[$5] = $1;')')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: SET_ZERO_VECTOR_ELEMENT(x, i, x_type)  ... x[i] = 0
dnl ----------------------------------------------------------------------
define(`SET_ZERO_VECTOR_ELEMENT', `ifelse( 
        `$3', `complex_S', `$1[$2] = 0.0; $1[$2+1] = 0.0;',
        `$3', `complex_D', `$1[$2] = 0.0; $1[$2+1] = 0.0;',
        `$1[$2] = 0.0;')')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: SET_ROUND_VECTOR_ELEMENT(x, i, val, val_type) ... x[i] = val
dnl        x[i] may be lower precision than val. x may not be extra
dnl        precision.
dnl ----------------------------------------------------------------------
define(`SET_ROUND_VECTOR_ELEMENT', `ifelse(
	`$4', `real_E', `$1[$2] = HEAD($3);',
	`$4', `complex_S', `$1[$2] = $3[0]; $1[$2+1] = $3[1];',
	`$4', `complex_D', `$1[$2] = $3[0]; $1[$2+1] = $3[1];',
	`$4', `complex_E', `$1[$2] = HEAD($3)[0]; $1[$2+1] = HEAD($3)[1];',
	`$1[$2] = $3;')')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: GET_MATRIX_ELEMENT(aij, A, i, j, lda, type) ... aij = A[i,j]
dnl ----------------------------------------------------------------------
define(`GET_MATRIX_ELEMENT', `ifelse(
        `$6', `complex_S', `$1[0] = $2[$3+$4*$5]; $1[1] = $2[$3+$4*$5+1];',
        `$6', `complex_D', `$1[0] = $2[$3+$4*$5]; $1[1] = $2[$3+$4*$5+1];',
        `$1 = $2[$3+$4*$5];')')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: SET_MATRIX_ELEMENT(A, i, j, lda, val, type) ... A[i,j] = val
dnl ----------------------------------------------------------------------
define(`SET_MATRIX_ELEMENT', `ifelse(
        `$6', `complex_S', `$1[$2+$3*$4] = $5[0]; $1[$2+$3*$4+1] = $5[1];',
        `$6', `complex_D', `$1[$2+$3*$4] = $5[0]; $1[$2+$3*$4+1] = $5[1];',
        `$1[$2] = $5;')')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: SET_ZERO_MATRIX_ELEMENT(A, i, j, lda, type) ... A[i,j] = 0
dnl ----------------------------------------------------------------------
define(`SET_ZERO_MATRIX_ELEMENT', `ifelse(
        `$5', `complex_S', `$1[$2+$3*$4] = 0.0; $1[$2+$3*$4+1] = 0.0;',
        `$5', `complex_D', `$1[$2+$3*$4] = 0.0; $1[$2+$3*$4+1] = 0.0;',
        `$1[$2] = 0.0;')')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: SET_ROUND_MATRIX_ELEMENT(A, i, j, lda, val, val_type)
dnl        ... A[i,j] = (round)val
dnl        A may be lower precision than val.
dnl ----------------------------------------------------------------------
define(`SET_ROUND_MATRIX_ELEMENT', `ifelse(
	`$6', `real_E', `$1[$2] = HEAD($3);',
	`$6', `complex_S', `$1[$2] = $3[0]; $1[$2+1] = $3[1];',
	`$6', `complex_D', `$1[$2] = $3[0]; $1[$2+1] = $3[1];',
	`$6', `complex_E', `$1[$2] = HEAD($3)[0]; $1[$2+1] = HEAD($3)[1];',
	`$1[$2] = $3;')')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: TEST_0(val, val_type)  ... returns a boolean expression that
dnl        evaluates true if val == 0.
dnl ----------------------------------------------------------------------
define(`TEST_0', `ifelse(
        `$2', `complex_S', `$1[0] == 0.0 && $1[1] == 0.0',
        `$2', `complex_D', `$1[0] == 0.0 && $1[1] == 0.0',
        `$1 == 0.0')')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: TEST_1(val, val_type)  ... returns a boolean expression that
dnl        evaluates true if val == 1.
dnl ----------------------------------------------------------------------
define(`TEST_1', `ifelse(
        `$2', `complex_S', `($1[0] == 1.0 && $1[1] == 0.0)',
        `$2', `complex_D', `($1[0] == 1.0 && $1[1] == 0.0)',
        `$1 == 1.0')')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: MUL(c, c_type, a, a_type, b, b_type) ... c = a * b
dnl ----------------------------------------------------------------------
define(`MUL', `ifelse(
	`$2&&$4&&$6', `real_E&&real_D&&real_D',
        `{
  	    /* Compute double_double = double * double. */
            double a1, a2, b1, b2, con;
            
            con = $3 * split;
            a1 = con - $3;
            a1 = con - a1;
            a2 = $3 - a1;
            con = $5 * split;
            b1 = con - $5;
            b1 = con - b1;
            b2 = $5 - b1;

            HEAD($1) = $3 * $5;
            TAIL($1) = (((a1 * b1 - HEAD($1)) + a1 * b2) + a2 * b1) + a2 * b2;
            }',
	`$2&&$4&&$6', `real_E&&real_E&&real_D',
	`{
  	    /* Compute double-double = double-double * double. */
            double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	    con = HEAD($3) * split;
            a11 = con - HEAD($3);
	    a11 = con - a11;
	    a21 = HEAD($3) - a11;
	    con = $5 * split;
	    b1 = con - $5;
	    b1 = con - b1;
	    b2 = $5 - b1;

	    c11 = HEAD($3) * $5;
	    c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

 	    c2 = TAIL($3) * $5;
	    t1 = c11 + c2;
	    t2 = (c2 - (t1 - c11)) + c21;
dnl	    e = t1 - c11;
dnl	    t2 = ((c2 - e) + (c11 - (t1 - e))) + c21;
	   
	    HEAD($1) = t1 + t2;
	    TAIL($1) = t2 - (HEAD($1) - t1);
            }',
	`$2&&$4&&$6', `real_E&&real_S&&real_S',
	    `HEAD($1) = (double) $3 * $5; TAIL($1) = 0.0;',
	`$2&&$4&&$6', `real_E&&real_D&&real_S',
        `{
            double dt = (double) $5;
            MUL($1, $2, $3, $4, dt, real_D)
         }',
	`$2&&$4&&$6', `real_E&&real_S&&real_D',
        `{
            double dt = (double) $3;
            MUL($1, $2, dt, real_D, $5, $6)
         }',
	`$2&&$4&&$6', `real_E&&real_E&&real_S',
        `{
            double dt = (double) $5;
            MUL($1, $2, $3, $4, dt, real_D)
         }',
	`$2&&$4&&$6', `real_D&&real_S&&real_S', `$1 = (double) $3 * $5;',
	`$2&&$4&&$6', `complex_S&&real_S&&real_S', 
		`$1[0] = $3 * $5; $1[1] = 0.0;',
	`$2&&$4&&$6', `complex_D&&real_D&&real_D',
		`$1[0] = $3 * $5; $1[1] = 0.0;',
        `$2&&$4&&$6', `complex_D&&real_S&&real_S',
		`$1[0] = (double) $3 * $5; $1[1] = 0.0;',
        `$2&&$4&&$6', `complex_E&&real_S&&real_S',
		 `HEAD($1)[0] = (double) $3 * $5; HEAD($1)[1] = 0.0;
		  TAIL($1)[0] = TAIL($1)[1] = 0.0;',
        `$2&&$4&&$6', `complex_E&&real_D&&real_D',
        `{
  	    /* Compute complex double_double = double * double. */
            double a1, a2, b1, b2, con;
            
            con = $3 * split;
            a1 = con - $3;
            a1 = con - a1;
            a2 = $3 - a1;
            con = $5 * split;
            b1 = con - $5;
            b1 = con - b1;
            b2 = $5 - b1;

            HEAD($1)[0] = $3 * $5;
            TAIL($1)[0] = (((a1 * b1 - HEAD($1)[0]) + a1 * b2) + a2 * b1) + a2 * b2;
	    HEAD($1)[1] = TAIL($1)[1] = 0.0;
            }',
        `$2',`complex_S', `S_CMPLX_MUL($1, $2, $3, $4, $5, $6)',
        `$2',`complex_D', `D_CMPLX_MUL($1, $2, $3, $4, $5, $6)',
        `$2',`complex_E', `E_CMPLX_MUL($1, $2, $3, $4, $5, $6)',
	`$1 = $3 * $5;')')dnl
dnl 
dnl S_CMPLX_MUL(c, c_type, a, a_type, b, b_type) ... c_type is complex_S
dnl
define(`S_CMPLX_MUL', `ifelse(
        `$4&&$6', `complex_S&&real_S', `C_TIMES_R($1, $3, $5)',
        `$4&&$6', `real_S&&complex_S', `C_TIMES_R($1, $5, $3)',
        `$4&&$6', `complex_S&&complex_S', `C_TIMES_C($1, $2, $3, $5)'
       )')dnl
dnl
dnl D_CMPLX_MUL(c, c_type, a, a_type, b, b_type) ... c_type is complex_D
dnl
define(`D_CMPLX_MUL', `ifelse(
        `$4&&$6', `complex_S&&real_S', `D_C_TIMES_R($1, $3, $5)',
        `$4&&$6', `real_S&&complex_S', `D_C_TIMES_R($1, $5, $3)',
        `$4&&$6', `complex_D&&real_S', `C_TIMES_R($1, $3, $5)',
        `$4&&$6', `real_S&&complex_D', `C_TIMES_R($1, $5, $3)',
        `$4&&$6', `complex_D&&real_D', `C_TIMES_R($1, $3, $5)',
        `$4&&$6', `real_D&&complex_D', `C_TIMES_R($1, $5, $3)',
        `$4&&$6', `real_D&&complex_S', `C_TIMES_R($1, $5, $3)',
        `$4&&$6', `complex_S&&real_D', `C_TIMES_R($1, $3, $5)',
        `C_TIMES_C($1, $2, $3, $5)' )')dnl
dnl
dnl E_CMPLX_MUL(c, c_type, a, a_type, b, b_type) ... c_type is complex_E
dnl
define(`E_CMPLX_MUL', `ifelse(
        `$4&&$6',`complex_D&&real_D', `E_C_TIMES_R($1, $2, $3, $4, $5, $6)',
        `$4&&$6',`real_D&&complex_D', `E_C_TIMES_R($1, $2, $5, $6, $3, $4)',
	`$4&&$6',`real_E&&complex_D', `E_C_TIMES_R($1, $2, $5, $6, $3, $4)',
        `$4&&$6',`complex_E&&real_D', `E_C_TIMES_R($1, $2, $3, $4, $5, $6)',
        `$4&&$6',`complex_D&&complex_D', `E_C_TIMES_C($1, $2, $3, $4, $5, $6)',
        `$4&&$6',`complex_E&&complex_D', `E_C_TIMES_C($1, $2, $3, $4, $5, $6)',
        `$4&&$6',`complex_S&&complex_S', `E_C_TIMES_C($1, $2, $3, $4, $5, $6)',
        `$4&&$6',`complex_D&&complex_S', `E_C_TIMES_C($1, $2, $3, $4, $5, $6)',
        `$4&&$6',`complex_S&&complex_D', `E_C_TIMES_C($1, $2, $5, $6, $3, $4)',
        `$4&&$6',`complex_E&&complex_S', `E_C_TIMES_C($1, $2, $3, $4, $5, $6)',
        `$4&&$6',`complex_S&&real_S', `E_C_TIMES_R($1, $2, $3, $4, $5, $6)',
        `$4&&$6',`real_S&&complex_S', `E_C_TIMES_R($1, $2, $5, $6, $3, $4)',
        `$4&&$6',`complex_S&&real_D', `E_C_TIMES_R($1, $2, $3, $4, $5, $6)',
        `$4&&$6',`real_D&&complex_S', `E_C_TIMES_R($1, $2, $5, $6, $3, $4)',
        `$4&&$6',`complex_D&&real_S', `E_C_TIMES_R($1, $2, $3, $4, $5, $6)',
        `$4&&$6',`complex_E&&real_S', `E_C_TIMES_R($1, $2, $3, $4, $5, $6)',
        `$4&&$6',`real_E&&complex_S', `E_C_TIMES_R($1, $2, $5, $6, $3, $4)',
        `miss($1, $2, $3, $4, $5, $6)' )')dnl
define(`D_C_TIMES_R',
        `{
           $1[0] = (double) $2[0] * $3; $1[1] = (double) $2[1] * $3;
         }')dnl
define(`C_TIMES_R',
        `{
           $1[0] = $2[0] * $3; $1[1] = $2[1] * $3;
         }')dnl
define(`C_TIMES_C', `ifelse(
        `$2', `complex_D', 
        `{
           $1[0] = (double) $3[0] * $4[0] - (double) $3[1] * $4[1];
           $1[1] = (double) $3[0] * $4[1] + (double) $3[1] * $4[0];
         }',
        `{
           $1[0] = $3[0] * $4[0] - $3[1] * $4[1];
           $1[1] = $3[0] * $4[1] + $3[1] * $4[0];
         }')')dnl
define(`E_C_TIMES_R', `ifelse(
        `$4&&$6', `complex_S&&real_S',
        `{
           HEAD($1)[0] = (double) $3[0] * $5; TAIL($1)[0] = 0.0;
           HEAD($1)[1] = (double) $3[1] * $5; TAIL($1)[1] = 0.0;
         }',
        `$4&&$6', `complex_S&&real_D',
        `{
           HEAD($1)[0] = (double) $3[0] * $5; TAIL($1)[0] = 0.0;
           HEAD($1)[1] = (double) $3[1] * $5; TAIL($1)[1] = 0.0;
         }',
        `$4&&$6', `complex_S&&real_E',
        `{
           DECLARE(e1, real_E)
	   DECLARE(dt, real_D)
	   dt = (double) $3[0];
           MUL(e1, real_E, $5, real_E, dt, real_D)
	   HEAD($1)[0] = HEAD(e1); TAIL($1)[0] = TAIL(e1);
	   dt = (double) $3[1];
           MUL(e1, real_E, $5, real_E, dt, real_D)
	   HEAD($1)[1] = HEAD(e1); TAIL($1)[1] = TAIL(e1);
         }',
        `$4&&$6', `complex_D&&real_S',
        `{
	   double dt = (double) $5;
	   MUL($1, $2, $3, $4, dt, real_D)
         }',
        `$4&&$6', `complex_E&&real_S',
        `{
           double dt = (double) $5;
	   MUL($1, $2, $3, $4, dt, real_D)
         }',
	`$4', `complex_D',
        `{
	   /* Compute complex-extra = complex-double * real. */
           DECLARE(t, real_E)
	   MUL(t, real_E, $5, $6, $3[0], real_D)
           HEAD($1)[0] = HEAD(t); TAIL($1)[0] = TAIL(t);
           MUL(t, real_E, $5, $6, $3[1], real_D)
           HEAD($1)[1] = HEAD(t); TAIL($1)[1] = TAIL(t);
	}',
        `$4', `complex_E',
        `{
           /* Compute complex-extra = complex-extra * real. */
           DECLARE(a0, real_E)
           DECLARE(a1, real_E)
           DECLARE(t, real_E)
	   HEAD(a0) = HEAD($3)[0]; TAIL(a0) = TAIL($3)[0];
	   HEAD(a1) = HEAD($3)[1]; TAIL(a1) = TAIL($3)[1];
           MUL(t, real_E, a0, real_E, $5, $6)
           HEAD($1)[0] = HEAD(t); TAIL($1)[0] = TAIL(t);
           MUL(t, real_E, a1, real_E, $5, $6)
           HEAD($1)[1] = HEAD(t); TAIL($1)[1] = TAIL(t);
         }'
       )')dnl
define(`E_C_TIMES_C', `ifelse(
        `$4&&$6', `complex_S&&complex_S', 
	`{
	    DECLARE(e1, real_E)
	    DECLARE(d1, real_D)
	    DECLARE(d2, real_D)
	    /* Real part */
	    d1 = (double) $3[0] * $5[0];
	    d2 = (double) -$3[1] * $5[1];
	    ADD(e1, real_E, d1, real_D, d2, real_D)
	    HEAD($1)[0] = HEAD(e1); TAIL($1)[0] = TAIL(e1);
	    /* imaginary part */
	    d1 = (double) $3[0] * $5[1];
	    d2 = (double) $3[1] * $5[0];
	    ADD(e1, real_E, d1, real_D, d2, real_D)
	    HEAD($1)[1] = HEAD(e1); TAIL($1)[1] = TAIL(e1);
	 }',
        `$4&&$6', `complex_D&&complex_D', 
        `{
           /* Compute complex-extra = complex-double * complex-double. */
           DECLARE(t1, real_E)
           DECLARE(t2, real_E)
           /* Real part */
           MUL(t1, real_E, $3[0], real_D, $5[0], real_D)
           MUL(t2, real_E, $3[1], real_D, $5[1], real_D)
           HEAD(t2) = -HEAD(t2); TAIL(t2) = -TAIL(t2);
           ADD(t1, real_E, t1, real_E, t2, real_E)
           HEAD($1)[0] = HEAD(t1); TAIL($1)[0] = TAIL(t1);
           /* Imaginary part */
           MUL(t1, real_E, $3[1], real_D, $5[0], real_D)
           MUL(t2, real_E, $3[0], real_D, $5[1], real_D)
           ADD(t1, real_E, t1, real_E, t2, real_E)
           HEAD($1)[1] = HEAD(t1); TAIL($1)[1] = TAIL(t1);
         }',
        `$4&&$6', `complex_D&&complex_S',
	`{
	   DECLARE(cd, complex_D)
	   cd[0] = (double) $5[0]; cd[1] = (double) $5[1];
	   MUL($1, $2, $3, $4, cd, complex_D)
	 }',
        `$4&&$6', `complex_E&&complex_S',
	`{
	   DECLARE(cd, complex_D)
	   cd[0] = (double) $5[0]; cd[1] = (double) $5[1];
	   MUL($1, $2, $3, $4, cd, complex_D)
	 }',
        `$4&&$6', `complex_E&&complex_D',
	`{
           /* Compute complex-extra = complex-extra * complex-double. */
           DECLARE(a0, real_E)
           DECLARE(a1, real_E)
           DECLARE(t1, real_E)
           DECLARE(t2, real_E)
	   HEAD(a0) = HEAD($3)[0]; TAIL(a0) = TAIL($3)[0];
	   HEAD(a1) = HEAD($3)[1]; TAIL(a1) = TAIL($3)[1];
	   /* real part */
           MUL(t1, real_E, a0, real_E, $5[0], real_D)
           MUL(t2, real_E, a1, real_E, $5[1], real_D)
           HEAD(t2) = -HEAD(t2); TAIL(t2) = -TAIL(t2);
           ADD(t1, real_E, t1, real_E, t2, real_E)
           HEAD($1)[0] = HEAD(t1); TAIL($1)[0] = TAIL(t1);
           /* imaginary part */
           MUL(t1, real_E, a1, real_E, $5[0], real_D)
           MUL(t2, real_E, a0, real_E, $5[1], real_D)
           ADD(t1, real_E, t1, real_E, t2, real_E)
           HEAD($1)[1] = HEAD(t1); TAIL($1)[1] = TAIL(t1);
         }'
       )')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: DIV(c, c_type, a, a_type, b, b_type) ... c = a / b
dnl ----------------------------------------------------------------------
define(`DIV', `ifelse(
        `$2&&$4&&$6', `real_E&&real_E&&real_E',
        `{
	    double q1, q2, q3;
	    double a1, a2, b1, b2;
	    double p1, p2, c;
	    double s1, s2, v;
	    double t1, t2;
	    double r1, r2;
	    double cona, conb;

            q1 = HEAD($3) / HEAD($5);   /*  approximate quotient */

            /*  Compute  q1 * b  */
	    cona = q1 * split;
	    conb = HEAD($5) * split;
	    a1 = cona - (cona - q1);
	    b1 = conb - (conb - HEAD($5));
	    a2 = q1 - a1;
	    b2 = HEAD($5) - b1;

	    /*  (p1, p2) is the product of high order terms. */
	    p1 = q1 * HEAD($5);
	    p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;
	
	    /*  Compute the low-order term */
            c = q1 * TAIL($5);


	    /*  Compute  (s1, s2) = (p1, p2) + c */
	    s1 = p1 + c;
	    v = s1 - p1;
            s2 = ((c - v) + (p1 - (s1 - v))) + p2;

            /*  Renormalize. */
            p1 = s1 + s2;
	    p2 = s2 - (p1 - s1);

	    
	    /*  Compute  a - (p1, p2)    */
            s1 = HEAD($3) - p1;
	    v = s1 - HEAD($3);
            s2 = (HEAD($3) - (s1 - v)) - (p1 + v);

	    t1 = TAIL($3) - p2;
	    v = t1 - TAIL($3);
	    t2 = (TAIL($3) - (t1 - v)) - (p2 + v);

	    s2 += t1;
            t1 = s1 + s2;
            s2 = s2 - (t1 - s1);

	    t2 += s2;
            r1 = t1 + t2;
            r2 = t2 - (r1 - t1);


	    /*  Compute the next quotient. */
            q2 = r1 / HEAD($5);


	    /*  Compute residual   r1 - q2 * b          */
	    cona = q2 * split;
	    a1 = cona - (cona - q2);
	    a2 = q2 - a1;

	    /*  (p1, p2) is the product of high order terms. */
	    p1 = q2 * HEAD($5);
	    p2 = (((a1 * b1 - p1) + a1 * b2) + a2 * b1) + a2 * b2;
	
	    /*  Compute the low-order term */
            c = q2 * TAIL($5);

	    /*  Compute  (s1, s2) = (p1, p2) + c */
	    s1 = p1 + c;
	    v = s1 - p1;
            s2 = ((c - v) + (p1 - (s1 - v))) + p2;

            /*  Renormalize. */
            p1 = s1 + s2;
	    p2 = s2 - (p1 - s1);


	    /*  Compute  (r1, r2) - (p1, p2)    */
            s1 = r1 - p1;
	    v = s1 - r1;
            s2 = (r1 - (s1 - v)) - (p1 + v);

	    t1 = r2 - p2;
	    v = t1 - r2;
	    t2 = (r2 - (t1 - v)) - (p2 + v);

	    s2 += t1;
            t1 = s1 + s2;
            s2 = s2 - (t1 - s1);

	    t2 += s2;
            s1 = t1 + t2;

	    /*  Compute the last correction. */
	    q3 = s1 / HEAD($5);

	    /* Renormalize q1, q2, q3. */
            s1 = q2 + q3;
            s2 = q3 - (s1 - q2);

            HEAD($1) = q1 + s1;
            t1 = s1 - (HEAD($1) - q1);

	    TAIL($1) = s2 + t1;

         }',
	`$2&&$4&&$6', `real_E&&real_E&&real_D',
        `{
  	    /* Compute double-double = double-double / double,
               using a Newton iteration scheme. */
            double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

            /* Compute a DP approximation to the quotient. */
            t1 = HEAD($3) / $5;

            /* Split t1 and b into two parts with at most 26 bits each,
               using the Dekker-Veltkamp method. */
            con = t1 * split;
            t11 = con - (con - t1);
            t21 = t1 - t11;
            con = $5 * split;
            b1 = con - (con - $5);
            b2 = $5 - b1;

            /* Compute t1 * b using Dekker method. */
            t12 = t1 * $5;
            t22 = (((t11 * b1 - t12) + t11 * b2) + t21 * b1) + t21 * b2;

	    /* Compute dda - (t12, t22) using Knuth trick. */
	    t11 = HEAD($3) - t12;
	    e = t11 - HEAD($3);
	    t21 = ((-t12 - e) + (HEAD($3) - (t11 - e))) + TAIL($3) - t22;

	    /* Compute high-order word of (t11, t21) and divide by b. */
	    t2 = (t11 + t21) / $5;
	    
	    /* The result is t1 + t2, after normalization. */
	    HEAD($1) = t1 + t2;
	    TAIL($1) = t2 - (HEAD($1) - t1);
         }',
	`$2&&$4&&$6', `real_E&&real_E&&real_S',
        `{
            double dt = (double) $5;
            DIV($1, $2, $3, $4, dt, real_D)
         }',
	`$2&&$4&&$6', `real_E&&real_D&&real_D',
        `{
            DECLARE(t, real_E)
            HEAD(t) = $3; TAIL(t) = 0.0;
            DIV($1, $2, t, real_E, $5, $6)
         }',
	`$2&&$4&&$6', `real_D&&real_S&&real_S', `$1 = (double) $3 / $5;',
        `$2',`complex_S', `S_CMPLX_DIV($1, $2, $3, $4, $5, $6)',
        `$2',`complex_D', `D_CMPLX_DIV($1, $2, $3, $4, $5, $6)',
        `$2',`complex_E', `E_CMPLX_DIV($1, $2, $3, $4, $5, $6)',
	`$1 = $3 / $5;')')dnl
dnl
dnl S_CMPLX_DIV(c, c_type, a, a_type, b, b_type) ... c_type is complex_S
dnl
define(`S_CMPLX_DIV', `ifelse(
        `$4&&$6', `complex_S&&real_S',
	`  DIV($1[0], real_S, $3[0], real_S, $5, $6)
	   DIV($1[1], real_S, $3[1], real_S, $5, $6)',
        `$4&&$6', `complex_D&&real_S', 
	`  DIV($1[0], real_S, $3[0], real_D, $5, $6)
	   DIV($1[1], real_S, $3[1], real_D, $5, $6)',
        `$4&&$6', `complex_E&&real_S', 
	`{
            DECLARE(a, real_E)
            DECLARE(b, real_E)
            HEAD(a) = HEAD($3)[0]; TAIL(a) = TAIL($3)[0];
            DIV(b, real_E, a, real_E, $5, $6)
            $1[0] = TAIL(b);
	    HEAD(a) = HEAD($3)[1]; TAIL(a) = TAIL($3)[1];
            DIV(b, real_E, a, real_E, $5, $6)
            $1[1] = TAIL(b);
         }',
        `$4&&$6', `complex_S&&complex_S',
	`{
            double S = 1.0, eps, ov, un;
            double abs_a, abs_b, abs_c, abs_d, ab, cd;
	    DECLARE(r, real_D)
	    DECLARE(t, real_D)
            DECLARE(q, complex_D)

            eps = pow(2.0, -24.0);
            un = pow(2.0, -126.0);
            ov = pow(2.0, 128.0) * (1 - eps);
            abs_a = fabs((double) $3[0]);
            abs_b = fabs((double) $3[1]);
            abs_c = fabs((double) $5[0]);
            abs_d = fabs((double) $5[1]);
            ab = MAX(abs_a, abs_b);
            cd = MAX(abs_c, abs_d);

            /* Scaling */
            if ( ab > ov/16 ) { /* scale down a, b */
                $3[0] /= 16; $3[1] /= 16;  S = S * 16;
            }
            if ( cd > ov/16 ) { /* scale down c, d */
                $5[0] /= 16; $5[1] /= 16;  S = S / 16;
            }
            if ( ab < un/eps*2 ) { /* scale up a, b */
                t = 2.0 / (eps*eps);
                $3[0] *= t; $3[1] *= t;  S = S / t;
            }
            if ( cd < un/eps*2) { /* scale up c, d */
                t = 2.0 / (eps*eps);
                $5[0] *= t; $5[1] *= t;  S = S * t;
            }
     
            /* Now un/eps*2 <= (a, b, c, d) >= ov/16 */
            if ( abs_c > abs_d ) {
               r = $5[1] / $5[0];
               t = 1 / ($5[0] + $5[1] * r);
               q[0] = ($3[0] + $3[1] * r) * t;
               q[1] = ($3[1] - $3[0] * r) * t;
            } else {
               r = $5[0] / $5[1];
               t = 1 / ($5[1] + $5[0] * r);
               q[0] = ($3[1] + $3[0] * r) * t;
               q[1] = (-$3[0] + $3[1] * r) * t;
            }
            /* Scale back */
            $1[0] = q[0] * S;  $1[1] = q[1] * S;
         }',
        `$4&&$6', `complex_D&&complex_S',
	`{
            double S = 1.0, eps, ov, un, eps1, ov1, un1;
            double abs_a, abs_b, abs_c, abs_d, ab, cd;
	    DECLARE(r, real_D)
	    DECLARE(t, real_D)
            DECLARE(q, complex_D)

            eps = pow(2.0, -24.0);
            un = pow(2.0, -126.0);
            ov = pow(2.0, 128.0) * (1 - eps);
            eps1 = pow(2.0, -53.0);
            un1 = pow(2.0, -1022.0);
            ov1 = 1.79769313486231571e+308;
              /* = (pow(2.0, 1023.0) * (1 - eps1)) * 2.0; */
            abs_a = fabs((double) $3[0]);
            abs_b = fabs((double) $3[1]);
            abs_c = fabs((double) $5[0]);
            abs_d = fabs((double) $5[1]);
            ab = MAX(abs_a, abs_b);
            cd = MAX(abs_c, abs_d);

            /* Scaling */
            if ( ab > ov1/16 ) { /* scale down a, b */
                $3[0] /= 16; $3[1] /= 16;  S = S * 16;
            }
            if ( cd > ov/16 ) { /* scale down c, d */
                $5[0] /= 16; $5[1] /= 16;  S = S / 16;
            }
            if ( ab < un1/eps1*2 ) { /* scale up a, b */
                t = 2.0 / (eps1*eps1);
                $3[0] *= t; $3[1] *= t;  S = S / t;
            }
            if ( cd < un/eps*2) { /* scale up c, d */
                t = 2.0 / (eps*eps);
                $5[0] *= t; $5[1] *= t;  S = S * t;
            }
     
            /* Now un/eps*2 <= (a, b, c, d) >= ov/16 */
            if ( abs_c > abs_d ) {
               r = $5[1] / $5[0];
               t = 1 / ($5[0] + $5[1] * r);
               q[0] = ($3[0] + $3[1] * r) * t;
               q[1] = ($3[1] - $3[0] * r) * t;
            } else {
               r = $5[0] / $5[1];
               t = 1 / ($5[1] + $5[0] * r);
               q[0] = ($3[1] + $3[0] * r) * t;
               q[1] = (-$3[0] + $3[1] * r) * t;
            }
            /* Scale back */
            $1[0] = q[0] * S;  $1[1] = q[1] * S;
         }',
        `$4&&$6', `complex_E&&complex_S', 
	`{
            wrong type $2&&$4&&$6
	 }')')dnl ... end S_CMPLX_DIV
dnl
dnl D_CMPLX_DIV(c, c_type, a, a_type, b, b_type) ... c_type is complex_D
dnl
define(`D_CMPLX_DIV', `ifelse(
        `$4&&$6', `complex_D&&real_S', 
	`{   double dt = $5;
            DIV($1[0], real_D, $3[0], real_D, dt, real_D)
	    DIV($1[1], real_D, $3[1], real_D, dt, real_D)}',			      
	`$4&&$6', `complex_D&&real_D', 
	`   DIV($1[0], real_D, $3[0], real_D, $5, $6)
	    DIV($1[1], real_D, $3[1], real_D, $5, $6)',
	`$4&&$6', `complex_E&&real_D', 
	`{
            DECLARE(a, real_E)
            DECLARE(b, real_E)
            HEAD(a) = HEAD($3)[0]; TAIL(a) = TAIL($3)[0];
            DIV(b, real_E, a, real_E, $5, $6)
            $1[0] = HEAD(b);
            HEAD(a) = HEAD($3)[1]; TAIL(a) = TAIL($3)[1];
            DIV(b, real_E, a, real_E, $5, $6)
            $1[1] = HEAD(b);
         }',
        `$4&&$6', `complex_D&&complex_S',
	`{
            double S = 1.0, eps, ov, un, eps1, ov1, un1;
            double abs_a, abs_b, abs_c, abs_d, ab, cd;
	    DECLARE(r, real_D)
	    DECLARE(t, real_D)
            DECLARE(q, complex_D)

            eps = pow(2.0, -24.0);
            un = pow(2.0, -126.0);
            ov = pow(2.0, 128.0) * (1 - eps);
            eps1 = pow(2.0, -53.0);
            un1 = pow(2.0, -1022.0);
            ov1 = 1.79769313486231571e+308; 
              /* = (pow(2.0, 1023.0) * (1 - eps1)) * 2.0; */
            abs_a = fabs($3[0]);
            abs_b = fabs($3[1]);
            abs_c = fabs((double) $5[0]);
            abs_d = fabs((double) $5[1]);
            ab = MAX(abs_a, abs_b);
            cd = MAX(abs_c, abs_d);

            /* Scaling */
            if ( ab > ov1/16 ) { /* scale down a, b */
                $3[0] /= 16; $3[1] /= 16;  S = S * 16;
            }
            if ( cd > ov/16 ) { /* scale down c, d */
                $5[0] /= 16; $5[1] /= 16;  S = S / 16;
            }
            if ( ab < un1/eps1*2 ) { /* scale up a, b */
                t = 2.0 / (eps1*eps1);
                $3[0] *= t; $3[1] *= t;  S = S / t;
            }
            if ( cd < un/eps*2) { /* scale up c, d */
                t = 2.0 / (eps*eps);
                $5[0] *= t; $5[1] *= t;  S = S * t;
            }
     
            /* Now un/eps*2 <= (a, b, c, d) >= ov/16 */
            if ( abs_c > abs_d ) {
               r = $5[1] / $5[0];
               t = 1 / ($5[0] + $5[1] * r);
               q[0] = ($3[0] + $3[1] * r) * t;
               q[1] = ($3[1] - $3[0] * r) * t;
            } else {
               r = $5[0] / $5[1];
               t = 1 / ($5[1] + $5[0] * r);
               q[0] = ($3[1] + $3[0] * r) * t;
               q[1] = (-$3[0] + $3[1] * r) * t;
            }
            /* Scale back */
            $1[0] = q[0] * S;  $1[1] = q[1] * S;
	 }',
        `$4&&$6', `complex_D&&complex_D', 
	`{
            double S = 1.0, eps, ov, un;
            double abs_a, abs_b, abs_c, abs_d, ab, cd;
	    DECLARE(r, real_D)
	    DECLARE(t, real_D)
            DECLARE(q, complex_D)

            eps = pow(2.0, -24.0);
            un = pow(2.0, -126.0);
            ov = pow(2.0, 128.0) * (1 - eps);
            abs_a = fabs($3[0]);
            abs_b = fabs($3[1]);
            abs_c = fabs((double) $5[0]);
            abs_d = fabs((double) $5[1]);
            ab = MAX(abs_a, abs_b);
            cd = MAX(abs_c, abs_d);

            /* Scaling */
            if ( ab > ov/16 ) { /* scale down a, b */
                $3[0] /= 16; $3[1] /= 16;  S = S * 16;
            }
            if ( cd > ov/16 ) { /* scale down c, d */
                $5[0] /= 16; $5[1] /= 16;  S = S / 16;
            }
            if ( ab < un/eps*2 ) { /* scale up a, b */
                t = 2.0 / (eps*eps);
                $3[0] *= t; $3[1] *= t;  S = S / t;
            }
            if ( cd < un/eps*2) { /* scale up c, d */
                t = 2.0 / (eps*eps);
                $5[0] *= t; $5[1] *= t;  S = S * t;
            }
     
            /* Now un/eps*2 <= (a, b, c, d) >= ov/16 */
            if ( abs_c > abs_d ) {
               r = $5[1] / $5[0];
               t = 1 / ($5[0] + $5[1] * r);
               q[0] = ($3[0] + $3[1] * r) * t;
               q[1] = ($3[1] - $3[0] * r) * t;
            } else {
               r = $5[0] / $5[1];
               t = 1 / ($5[1] + $5[0] * r);
               q[0] = ($3[1] + $3[0] * r) * t;
               q[1] = (-$3[0] + $3[1] * r) * t;
            }
            /* Scale back */
            $1[0] = q[0] * S;  $1[1] = q[1] * S;
	 }',
        `$4&&$6', `complex_E&&complex_S', 
	`{
            wrong type $2&&$4&&$6
	 }',
        `$4&&$6', `complex_E&&complex_D', 
	`{
            wrong type $2&&$4&&$6
	 }')')dnl ... end D_CMPLX_DIV
dnl
dnl E_CMPLX_DIV(c, c_type, a, a_type, b, b_type) ... c_type is complex_E
dnl
define(`E_CMPLX_DIV', `ifelse(
	`$4&&$6', `complex_D&&real_D', 
        `{
            DECLARE(t, real_E)
            DIV(t, real_E, $3[0], real_D, $5, $6)
            HEAD($1)[0] = HEAD(t); TAIL($1)[0] = TAIL(t);
	    DIV(t, real_E, $3[1], real_D, $5, $6)
            HEAD($1)[1] = HEAD(t); TAIL($1)[1] = TAIL(t);
         }',
	`$4&&$6', `complex_E&&real_S', 
        `{  DECLARE(a, real_E)
            DECLARE(b, real_E)
            HEAD(a) = HEAD($3)[0]; TAIL(a) = TAIL($3)[0];
            DIV(b, real_E, a, real_E, $5, real_D)
            HEAD($1)[0] = HEAD(b); TAIL($1)[0] = TAIL(b);
            HEAD(a) = HEAD($3)[1]; TAIL(a) = TAIL($3)[1];
            DIV(b, real_E, a, real_E, $5, real_D)
            HEAD($1)[1] = HEAD(b); TAIL($1)[1] = TAIL(b);
         }',
	`$4&&$6', `complex_E&&real_D', 
	`{
            DECLARE(a, real_E)
            DECLARE(b, real_E)
            HEAD(a) = HEAD($3)[0]; TAIL(a) = TAIL($3)[0];
            DIV(b, real_E, a, real_E, $5, $6)
            HEAD($1)[0] = HEAD(b); TAIL($1)[0] = TAIL(b);
            HEAD(a) = HEAD($3)[1]; TAIL(a) = TAIL($3)[1];
            DIV(b, real_E, a, real_E, $5, $6)
            HEAD($1)[1] = HEAD(b); TAIL($1)[1] = TAIL(b);
         }',
        `$4&&$6', `complex_D&&complex_S', 
	`{
	    DECLARE(t, complex_D)
	    t[0] = (double) $5[0]; t[1] = (double) $5[1];
	    DIV($1, $2, $3, $4, t, complex_D)
	 }',
        `$4&&$6', `complex_D&&complex_D', 
	`{
            double S = 1.0, eps, ov, un;
            double abs_a, abs_b, abs_c, abs_d, ab, cd;
            DECLARE(s, real_D)
	    DECLARE(r, real_D)
	    DECLARE(t, real_E)
	    DECLARE(t1, real_E)
	    DECLARE(t2, real_E)
            DECLARE(q, complex_E)

            eps = pow(2.0, -53.0);  /* double precision */
            un = pow(2.0, -1022.0);
            ov = 1.79769313486231571e+308;
              /* = (pow(2.0, 1023.0) * (1 - eps)) * 2.0; */
            abs_a = fabs($3[0]);
            abs_b = fabs($3[1]);
            abs_c = fabs((double)$5[0]);
            abs_d = fabs((double)$5[1]);
            ab = MAX(abs_a, abs_b);
            cd = MAX(abs_c, abs_d);

            /* Scaling */
            if ( ab > ov/16 ) { /* scale down a, b */
                $3[0] /= 16; $3[1] /= 16;  S = S * 16;
            }
            if ( cd > ov/16 ) { /* scale down c, d */
                $5[0] /= 16; $5[1] /= 16;  S = S / 16;
            }
            if ( ab < un/eps*2 ) { /* scale up a, b */
                s = 2.0 / (eps*eps);
                $3[0] *= s; $3[1] *= s;  S = S / s;
            }
            if ( cd < un/eps*2) { /* scale up c, d */
                s = 2.0 / (eps*eps);
                $5[0] *= s; $5[1] *= s;  S = S * s;
            }
     
            /* Now un/eps*2 <= (a,b,c,d) >= ov/16 */
            if ( abs_c > abs_d ) {
               r = $5[1] / $5[0];
               MUL(t, real_E, r, real_D, $5[1], real_D)
               ADD(t, real_E, t, real_E, $5[0], real_D)
               MUL(t1, real_E, $3[1], real_D, r, real_D)
               ADD(t1, real_E, t1, real_E, $3[0], real_D) /* a+b*r */
               DIV(t2, real_E, t1, real_E, t, real_E)
               HEAD(q)[0] = HEAD(t2); TAIL(q)[0] = TAIL(t2);
               MUL(t1, real_E, $3[0], real_D, r, real_D)
               HEAD(t1) = -HEAD(t1); TAIL(t1) = -TAIL(t1);
               ADD(t1, real_E, t1, real_E, $3[1], real_D) /* b-a*r */
               DIV(t2, real_E, t1, real_E, t, real_E)
               HEAD(q)[1] = HEAD(t2); TAIL(q)[1] = TAIL(t2);
            } else {
               r = $5[0] / $5[1];
               MUL(t, real_E, r, real_D, $5[0], real_D)
               ADD(t, real_E, t, real_E, $5[1], real_D)
               MUL(t1, real_E, $3[0], real_D, r, real_D)
               ADD(t1, real_E, t1, real_E, $3[1], real_D) /* b+a*r */
               DIV(t2, real_E, t1, real_E, t, real_E)
               HEAD(q)[0] = HEAD(t2); TAIL(q)[0] = TAIL(t2);
               MUL(t1, real_E, $3[1], real_D, r, real_D)
               SUB(t1, real_E, t1, real_E, $3[0], real_D) /* -a+b*r */
               DIV(t2, real_E, t1, real_E, t, real_E)
               HEAD(q)[1] = HEAD(t2); TAIL(q)[1] = TAIL(t2);
            }
            /* Scale back */
            if ( S == 1.0 ) {
                HEAD($1)[0] = HEAD(q)[0]; TAIL($1)[0] = TAIL(q)[0];
                HEAD($1)[1] = HEAD(q)[1]; TAIL($1)[1] = TAIL(q)[1];
            } else
                MUL($1, $2, q, complex_E, S, real_D)
	 }',
        `$4&&$6', `complex_E&&complex_S', 
	`{
            double S = 1.0, eps, ov, un, eps1, ov1, un1;
            double abs_a, abs_b, abs_c, abs_d, ab, cd;
            DECLARE(s, real_D)
	    DECLARE(r, real_D)
	    DECLARE(t, real_E)
	    DECLARE(t1, real_E)
	    DECLARE(t2, real_E)
            DECLARE(q, complex_E)

            eps = pow(2.0, -24.0);  /* single precision */
            un = pow(2.0, -126.0);
            ov = pow(2.0, 128.0) * (1 - eps);
            eps1 = pow(2.0, -104.0); /* extra precision */
            un1 = pow(2.0, -1022.0);
            ov1 = 1.79769313486231571e+308;
              /* = (pow(2.0, 1023.0) * (1 - eps1)) * 2.0 */
            abs_a = fabs(HEAD($3)[0]);
            abs_b = fabs(HEAD($3)[1]);
            abs_c = fabs((double)$5[0]);
            abs_d = fabs((double)$5[1]);
            ab = MAX(abs_a, abs_b);
            cd = MAX(abs_c, abs_d);

            /* Scaling */
            if ( ab > ov1/16 ) { /* scale down a, b */
                DIV($3, $4, $3, $4, 16.0, real_D)
                S = S * 16;
            }
            if ( cd > ov/16 ) { /* scale down c, d */
                $5[0] /= 16; $5[1] /= 16;  S = S / 16;
            }
            if ( ab < un1/eps1*2 ) { /* scale up a, b */
                s = 2.0 / (eps1*eps1);
                MUL($3, $4, $3, $4, s, real_D)
                S = S / s;
            }
            if ( cd < un/eps*2) { /* scale up c, d */
                s = 2.0 / (eps*eps);
                $5[0] *= s; $5[1] *= s;  S = S * s;
            }
     
            /* Now un1/eps1*2 <= (a,b) >= ov1/16, un/eps*2 <= (c,d) >= ov/16 */
            if ( abs_c > abs_d ) {
               r = $5[1] / $5[0];
               MUL(t, real_E, r, real_D, $5[1], real_S)
               ADD(t, real_E, t, real_E, $5[0], real_S)
               HEAD(t1) = HEAD($3)[1]; TAIL(t1) = TAIL($3)[1]; /* b */
               MUL(t2, real_E, t1, real_E, r, real_D)
               HEAD(t1) = HEAD($3)[0]; TAIL(t1) = TAIL($3)[0]; /* a */
               ADD(t2, real_E, t2, real_E, t1, real_E) /* a + b*r */
               DIV(t2, real_E, t2, real_E, t, real_E)
               HEAD(q)[0] = HEAD(t2); TAIL(q)[0] = TAIL(t2);
               HEAD(t1) = HEAD($3)[0]; TAIL(t1) = TAIL($3)[0]; /* a */
               MUL(t2, real_E, t1, real_E, r, real_D)
               HEAD(t1) = HEAD($3)[1]; TAIL(t1) = TAIL($3)[1]; /* b */
               SUB(t2, real_E, t1, real_E, t2, real_E) /* b - a*r */
               DIV(t2, real_E, t2, real_E, t, real_E)
               HEAD(q)[1] = HEAD(t2); TAIL(q)[1] = TAIL(t2);
            } else {
               r = $5[0] / $5[1];
               MUL(t, real_E, r, real_D, $5[0], real_S)
               ADD(t, real_E, t, real_E, $5[1], real_S)
               HEAD(t1) = HEAD($3)[0]; TAIL(t1) = TAIL($3)[0]; /* a */
               MUL(t2, real_E, t1, real_E, r, real_D)
               HEAD(t1) = HEAD($3)[1]; TAIL(t1) = TAIL($3)[1]; /* b */
               ADD(t2, real_E, t2, real_E, t1, real_E) /* b + a*r */
               DIV(t2, real_E, t2, real_E, t, real_E)
               HEAD(q)[0] = HEAD(t2); TAIL(q)[0] = TAIL(t2);
               HEAD(t1) = HEAD($3)[1]; TAIL(t1) = TAIL($3)[1]; /* b */
               MUL(t2, real_E, t1, real_E, r, real_D)
               HEAD(t1) = HEAD($3)[0]; TAIL(t1) = TAIL($3)[0]; /* a */
               SUB(t2, real_E, t2, real_E, t1, real_E) /* -a + b*r */
               DIV(t2, real_E, t2, real_E, t, real_E)
               HEAD(q)[1] = HEAD(t2); TAIL(q)[1] = TAIL(t2);
            }
            /* Scale back */
            if ( S == 1.0 ) {
                HEAD($1)[0] = HEAD(q)[0]; TAIL($1)[0] = TAIL(q)[0];
                HEAD($1)[1] = HEAD(q)[1]; TAIL($1)[1] = TAIL(q)[1];
            } else
                MUL($1, $2, q, complex_E, S, real_D)
	 }',
        `$4&&$6', `complex_E&&complex_D', 
	`{
            double S = 1.0, eps, ov, un, eps1, ov1, un1;
            double abs_a, abs_b, abs_c, abs_d, ab, cd;
            DECLARE(s, real_D)
	    DECLARE(r, real_D)
	    DECLARE(t, real_E)
	    DECLARE(t1, real_E)
	    DECLARE(t2, real_E)
            DECLARE(q, complex_E)

            eps = pow(2.0, -53.0);  /* double precision */
            un = pow(2.0, -1022.0);
            ov = 1.79769313486231571e+308;
              /* = (pow(2.0, 1023.0) * (1 - eps)) * 2.0 */
            eps1 = pow(2.0, -104.0); /* extra precision */
            un1 = pow(2.0, -1022.0);
            ov1 = 1.79769313486231571e+308;
              /* = (pow(2.0, 1023.0) * (1 - eps1)) * 2.0 */
            abs_a = fabs(HEAD($3)[0]);
            abs_b = fabs(HEAD($3)[1]);
            abs_c = fabs((double)$5[0]);
            abs_d = fabs((double)$5[1]);
            ab = MAX(abs_a, abs_b);
            cd = MAX(abs_c, abs_d);

            /* Scaling */
            if ( ab > ov1/16 ) { /* scale down a, b */
                DIV($3, $4, $3, $4, 16.0, real_D)
                S = S * 16;
            }
            if ( cd > ov/16 ) { /* scale down c, d */
                $5[0] /= 16; $5[1] /= 16;  S = S / 16;
            }
            if ( ab < un1/eps1*2 ) { /* scale up a, b */
                s = 2.0 / (eps1*eps1);
                MUL($3, $4, $3, $4, s, real_D)
                S = S / s;
            }
            if ( cd < un/eps*2) { /* scale up c, d */
                s = 2.0 / (eps*eps);
                $5[0] *= s; $5[1] *= s;  S = S * s;
            }
     
            /* Now un1/eps1*2 <= (a,b) >= ov1/16, un/eps*2 <= (c,d) >= ov/16 */
            if ( abs_c > abs_d ) {
               r = $5[1] / $5[0];
               MUL(t, real_E, r, real_D, $5[1], real_D)
               ADD(t, real_E, t, real_E, $5[0], real_D)
               HEAD(t1) = HEAD($3)[1]; TAIL(t1) = TAIL($3)[1]; /* b */
               MUL(t2, real_E, t1, real_E, r, real_D)
               HEAD(t1) = HEAD($3)[0]; TAIL(t1) = TAIL($3)[0]; /* a */
               ADD(t2, real_E, t2, real_E, t1, real_E) /* a + b*r */
               DIV(t2, real_E, t2, real_E, t, real_E)
               HEAD(q)[0] = HEAD(t2); TAIL(q)[0] = TAIL(t2);
               HEAD(t1) = HEAD($3)[0]; TAIL(t1) = TAIL($3)[0]; /* a */
               MUL(t2, real_E, t1, real_E, r, real_D)
               HEAD(t1) = HEAD($3)[1]; TAIL(t1) = TAIL($3)[1]; /* b */
               SUB(t2, real_E, t1, real_E, t2, real_E) /* b - a*r */
               DIV(t2, real_E, t2, real_E, t, real_E)
               HEAD(q)[1] = HEAD(t2); TAIL(q)[1] = TAIL(t2);
            } else {
               r = $5[0] / $5[1];
               MUL(t, real_E, r, real_D, $5[0], real_D)
               ADD(t, real_E, t, real_E, $5[1], real_D)
               HEAD(t1) = HEAD($3)[0]; TAIL(t1) = TAIL($3)[0]; /* a */
               MUL(t2, real_E, t1, real_E, r, real_D)
               HEAD(t1) = HEAD($3)[1]; TAIL(t1) = TAIL($3)[1]; /* b */
               ADD(t2, real_E, t2, real_E, t1, real_E) /* b + a*r */
               DIV(t2, real_E, t2, real_E, t, real_E)
               HEAD(q)[0] = HEAD(t2); TAIL(q)[0] = TAIL(t2);
               HEAD(t1) = HEAD($3)[1]; TAIL(t1) = TAIL($3)[1]; /* b */
               MUL(t2, real_E, t1, real_E, r, real_D)
               HEAD(t1) = HEAD($3)[0]; TAIL(t1) = TAIL($3)[0]; /* a */
               SUB(t2, real_E, t2, real_E, t1, real_E) /* -a + b*r */
               DIV(t2, real_E, t2, real_E, t, real_E)
               HEAD(q)[1] = HEAD(t2); TAIL(q)[1] = TAIL(t2);
            }
            /* Scale back */
            if ( S == 1.0 ) {
                HEAD($1)[0] = HEAD(q)[0]; TAIL($1)[0] = TAIL(q)[0];
                HEAD($1)[1] = HEAD(q)[1]; TAIL($1)[1] = TAIL(q)[1];
            } else
                MUL($1, $2, q, complex_E, S, real_D)
         }')')dnl ... end E_CMPLX_DIV
dnl
dnl ----------------------------------------------------------------------
dnl Usage: ADD(c, c_type, a, a_type, b, b_type) ... c = a + b
dnl ----------------------------------------------------------------------
define(`ADD', `ifelse(
	`$2&&$4&&$6',`real_E&&real_D&&real_D',
        `{
  	    /* Compute double-double = double + double. */
  	    double e, t1, t2;

	    /* Knuth trick. */
	    t1 = $3 + $5;
      	    e = t1 - $3;
	    t2 = (($5 - e) + ($3 - (t1 - e)));

	    /* The result is t1 + t2, after normalization. */
	    HEAD($1) = t1 + t2;
	    TAIL($1) = t2 - (HEAD($1) - t1);
            }',
	`$2&&$4&&$6',`real_E&&real_E&&real_D',
        `{
  	    /* Compute double-double = double-double + double. */
  	    double e, t1, t2;

	    /* Knuth trick. */
	    t1 = HEAD($3) + $5;
      	    e = t1 - HEAD($3);
	    t2 = (($5 - e) + (HEAD($3) - (t1 - e))) + TAIL($3);

	    /* The result is t1 + t2, after normalization. */
	    HEAD($1) = t1 + t2;
	    TAIL($1) = t2 - (HEAD($1) - t1);
            }',
	`$2&&$4&&$6',`real_E&&real_E&&real_E',
        `{ 
	    /* Compute double-double = double-double + double-double. */
	    double bv;
	    double s1, s2, t1, t2;

	    /* Add two hi words. */
            s1 = HEAD($3) + HEAD($5);
	    bv = s1 - HEAD($3);
            s2 = ((HEAD($5) - bv) + (HEAD($3) - (s1 - bv)));

	    /* Add two lo words. */
            t1 = TAIL($3) + TAIL($5);
            bv = t1 - TAIL($3);
            t2 = ((TAIL($5) - bv) + (TAIL($3) - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    HEAD($1) = t1 + t2;
            TAIL($1) = t2 - (HEAD($1) - t1);
         }',
	`$2&&$4&&$6', `real_D&&real_S&&real_S', `$1 = (double) $3 + $5;',
	`$2&&$4&&$6', `real_E&&real_E&&real_S', 
        `{
            double dt = (double) $5;
            ADD($1, $2, $3, $4, dt, real_D);
         }',
  `$2&&$4&&$6', `complex_D&&complex_E&&complex_E',
	`{
	    DECLARE(a, real_E)
	    DECLARE(b, real_E)
	    /* Real part */
	    HEAD(a) = HEAD($3)[0]; TAIL(a) = TAIL($3)[0];
	    HEAD(b) = HEAD($5)[0]; TAIL(b) = TAIL($5)[0];
	    ADD($1[0], real_D, a, real_E, b, real_E)
	    /* Imaginary part */
	    HEAD(a) = HEAD($3)[1]; TAIL(a) = TAIL($3)[1];
	    HEAD(b) = HEAD($5)[1]; TAIL(b) = TAIL($5)[1];
	    ADD($1[1], real_D, a, real_E, b, real_E)
	}',
  `$2&&$4&&$6', `complex_S&&complex_E&&complex_E',
	`{
	    DECLARE(a, real_E)
	    DECLARE(b, real_E)
	    /* Real part */
	    HEAD(a) = HEAD($3)[0]; TAIL(a) = TAIL($3)[0];
	    HEAD(b) = HEAD($5)[0]; TAIL(b) = TAIL($5)[0];
	    ADD($1[0], real_S, a, real_E, b, real_E)
	    /* Imaginary part */
	    HEAD(a) = HEAD($3)[1]; TAIL(a) = TAIL($3)[1];
	    HEAD(b) = HEAD($5)[1]; TAIL(b) = TAIL($5)[1];
	    ADD($1[1], real_S, a, real_E, b, real_E)
	}',
        `$2',`complex_S', `$1[0] = $3[0] + $5[0]; $1[1] = $3[1] + $5[1];',
        `$2',`complex_D', `$1[0] = $3[0] + $5[0]; $1[1] = $3[1] + $5[1];',
        `$2&&$4&&$6',`complex_E&&complex_E&&complex_S',
	`{
	    DECLARE(cd, complex_D)
	    cd[0] = (double) $5[0]; cd[1] = (double) $5[1];
	    ADD($1, $2, $3, $4, cd, complex_D)
	}',
        `$2&&$4&&$6',`complex_E&&complex_E&&complex_D',
        `{
            DECLARE(t, real_E)
            DECLARE(a, real_E)
            HEAD(a) = HEAD($3)[0]; TAIL(a) = TAIL($3)[0];
            ADD(t, real_E, a, real_E, $5[0], real_D)
            HEAD($1)[0] = HEAD(t); TAIL($1)[0] = TAIL(t);
            HEAD(a) = HEAD($3)[1]; TAIL(a) = TAIL($3)[1];
            ADD(t, real_E, a, real_E, $5[1], real_D)
            HEAD($1)[1] = HEAD(t); TAIL($1)[1] = TAIL(t);
         }',
        `$2&&$4&&$6',`complex_E&&complex_E&&complex_E',
	`{
	    DECLARE(t, real_E)
	    DECLARE(a, real_E)
	    DECLARE(b, real_E)
	    /* Real part */
	    HEAD(a) = HEAD($3)[0]; TAIL(a) = TAIL($3)[0];
	    HEAD(b) = HEAD($5)[0]; TAIL(b) = TAIL($5)[0];
	    ADD(t, real_E, a, real_E, b, real_E)
            HEAD($1)[0] = HEAD(t); TAIL($1)[0] = TAIL(t);
	    /* Imaginary part */
	    HEAD(a) = HEAD($3)[1]; TAIL(a) = TAIL($3)[1];
	    HEAD(b) = HEAD($5)[1]; TAIL(b) = TAIL($5)[1];
	    ADD(t, real_E, a, real_E, b, real_E)
            HEAD($1)[1] = HEAD(t); TAIL($1)[1] = TAIL(t);
	}',
  `$2&&$4&&$6', `real_D&&real_E&&real_E', 
  ` {
    /* Compute double-double = double-double + double-double. */
    double bv;
    double s1, s2, t1, t2;

    /* Add two hi words. */
    s1 = HEAD($3) + HEAD($5);
    bv = s1 - HEAD($3);
    s2 = ((HEAD($5) - bv) + (HEAD($3) - (s1 - bv)));

    /* Add two lo words. */
    t1 = TAIL($3) + TAIL($5);
    bv = t1 - TAIL($3);
    t2 = ((TAIL($5) - bv) + (TAIL($3) - (t1 - bv)));

    s2 += t1;

    /* Renormalize (s1, s2)  to  (t1, s2) */
    t1 = s1 + s2;
    s2 = s2 - (t1 - s1);

    t2 += s2;

    /* Renormalize (t1, t2)  */
    $1 = t1 + t2;
    }',
  `$2&&$4&&$6', `real_S&&real_E&&real_E', 
  ` {
    /* Compute double-double = double-double + double-double. */
    double bv;
    double s1, s2, t1, t2;

    /* Add two hi words. */
    s1 = HEAD($3) + HEAD($5);
    bv = s1 - HEAD($3);
    s2 = ((HEAD($5) - bv) + (HEAD($3) - (s1 - bv)));

    /* Add two lo words. */
    t1 = TAIL($3) + TAIL($5);
    bv = t1 - TAIL($3);
    t2 = ((TAIL($5) - bv) + (TAIL($3) - (t1 - bv)));

    s2 += t1;

    /* Renormalize (s1, s2)  to  (t1, s2) */
    t1 = s1 + s2;
    s2 = s2 - (t1 - s1);

    t2 += s2;

    /* Renormalize (t1, t2)  */
    $1 = t1 + t2;
    }',
	`$1 = $3 + $5;')')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: SUB(c, c_type, a, a_type, b, b_type) ... c = a - b
dnl ----------------------------------------------------------------------
define(`SUB', `ifelse(
	`$2&&$4&&$6',`real_E&&real_D&&real_D',
        `ADD($1, $2, $3, $4, (-$5), $6)',
	`$2&&$4&&$6',`real_E&&real_E&&real_D',
        `ADD($1, $2, $3, $4, (-$5), $6)',
	`$2&&$4&&$6',`real_E&&real_E&&real_E',
        `{
            DECLARE(bt, $6)
            HEAD(bt) = -HEAD($5); TAIL(bt) = -TAIL($5);
            ADD($1, $2, $3, $4, bt, $6)
         }',
	`$2&&$4&&$6', `real_D&&real_S&&real_S', `$1 = (double) $3 + (-$5);',
	`$2&&$4&&$6', `real_E&&real_E&&real_S', 
        `{
            double bt = -(double) $5;
            ADD($1, $2, $3, $4, bt, real_D);
         }',
        `$2',`complex_S', `$1[0] = $3[0] - $5[0]; $1[1] = $3[1] - $5[1];',
        `$2',`complex_D', `$1[0] = $3[0] - $5[0]; $1[1] = $3[1] - $5[1];',
        `$2&&$4&&$6',`complex_E&&complex_E&&complex_S',
	`{
	    DECLARE(bt, complex_D)
	    bt[0] = -(double) $5[0]; bt[1] = -(double) $5[1];
	    ADD($1, $2, $3, $4, bt, complex_D)
	}',
        `$2&&$4&&$6',`complex_E&&complex_E&&complex_D',
        `{
	    DECLARE(bt, complex_D)
	    bt[0] = -$5[0]; bt[1] = -$5[1];
	    ADD($1, $2, $3, $4, bt, $6)
         }',
        `$2&&$4&&$6',`complex_E&&complex_E&&complex_E',
	`{
	    DECLARE(at, real_E)
	    DECLARE(bt, real_E)
	    DECLARE(ct, real_E)
	    
	    /* Real part */
	    HEAD(at) = HEAD($3)[0]; TAIL(at) = TAIL($3)[0];
	    HEAD(bt) = -HEAD($5)[0]; TAIL(bt) = -TAIL($5)[0];
	    ADD(ct, real_E, at, real_E, bt, real_E)
            HEAD($1)[0] = HEAD(ct); TAIL($1)[0] = TAIL(ct);
	    /* Imaginary part */
	    HEAD(at) = HEAD($3)[1]; TAIL(at) = TAIL($3)[1];
	    HEAD(bt) = -HEAD($5)[1]; TAIL(bt) = -TAIL($5)[1];
	    ADD(ct, real_E, at, real_E, bt, real_E)
            HEAD($1)[1] = HEAD(ct); TAIL($1)[1] = TAIL(ct);
	}',
	`$1 = $3 + (-$5);')')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: DECLARE_VECTOR(var, var_type) ... declare a vector
dnl ----------------------------------------------------------------------
define(`DECLARE_VECTOR', `ifelse(
	`$2',`real_S', `float *$1;',
	`$2',`real_D', `double *$1;',
	`$2',`real_E', `double *HEAD($1), *TAIL($1);',
	`$2',`complex_S', `float *$1;',
	`$2',`complex_D', `double *$1;',
	`$2',`complex_E', `double *HEAD($1), *TAIL($1);'
      )')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: MALLOC_VECTOR(var, var_type, length) ... malloc a vector
dnl ----------------------------------------------------------------------
define(`MALLOC_VECTOR', `ifelse(
	`$2',`real_S', 
	     `$1 = (float*)blas_malloc($3 * sizeof(float));
	      if($3 > 0 && $1 == NULL){
                BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
              }',
	`$2',`real_D', 
             `$1 = (double*)blas_malloc($3 * sizeof(double));
              if($3 > 0 && $1 == NULL){
                BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
              }',
	`$2',`real_E', 
             `HEAD($1) = (double*)blas_malloc($3 * sizeof(double));
              TAIL($1) = (double*)blas_malloc($3 * sizeof(double));
              if($3 > 0 && (HEAD($1) == NULL || TAIL($1) == NULL)){
                BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
              }',
	`$2',`complex_S', 
             `$1 = (float*)blas_malloc($3 * sizeof(float)*2);
	      if($3 > 0 && $1 == NULL){
                BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
              }',
	`$2',`complex_D', 
             `$1 = (double*)blas_malloc($3 * sizeof(double)*2);
              if($3 > 0 && $1 == NULL){
                BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
              }',
	`$2',`complex_E', 
             `HEAD($1) = (double*)blas_malloc($3 * sizeof(double)*2); 
              TAIL($1) = (double*)blas_malloc($3 * sizeof(double)*2);
              if($3 > 0 && (HEAD($1) == NULL || TAIL($1) == NULL)){
                BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
              }',
      )')dnl
dnl
dnl
dnl -----------------------------------------------------------------------
dnl Usage: FREE_VECTOR(var, type)
dnl        free memory used by var
dnl -----------------------------------------------------------------------
dnl
define(`FREE_VECTOR', `ifelse(
  `$2', `real_E', `blas_free(HEAD($1)); blas_free(TAIL($1));',
  `$2', `complex_E', `blas_free(HEAD($1)); blas_free(TAIL($1));',
  `blas_free($1);')')dnl
dnl
dnl
