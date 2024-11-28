dnl
dnl Contains common m4 macros used in testing routines.
dnl
dnl
dnl standard precision combination to generate various testing procedures
define(`PREC_ARGS', 
  ``s, s, s', `d, d, d', `c, c, c', `z, z, z', 
   `c, s, s', `c, s, c', `c, c, s',
   `z, d, d', `z, d, z', `z, z, d', 
   `d, s, s', `d, s, d', `d, d, s', 
   `z, c, c', `z, c, z', `z, z, c'')dnl
dnl
dnl
dnl CABS(x) ... compute the absolute value of complex number x.
define(`CABS',`sqrt($1[0]*$1[0] + $1[1]*$1[1])')dnl
dnl
dnl
dnl ASSIGN_PTR_TO_SCALAR(a, a_type, b, b_type) ... set a = *b
define(`ASSIGN_PTR_TO_SCALAR', 
  `IF_REAL($2, `$1 = *$3;', `$1[0] = $3[0]; $1[1] = $3[1];')')dnl
dnl
dnl ASSIGN_SCALAR_TO_PTR(a, a_type, b, b_type) ... set *a = b
define(`ASSIGN_SCALAR_TO_PTR', 
  `IF_REAL($4, `*$1 = $3;', `$1[0] = $3[0]; $1[1] = $3[1];')')dnl
dnl
dnl
dnl ZERO_IMAG_PART(x, type) ... imag(x) = 0
define(`ZERO_IMAG_PART', `ifelse(
  REAL_COMPLEX($2), `complex', `$1[1] = 0.0;')')dnl
dnl
dnl
define(`IS_MIXED', `ifelse(
  __IS_MIXED_ABBREV($1, $2), `t', `t', 
  `$3', `', `f', 
  __IS_MIXED_ABBREV($2, $3))')dnl
dnl
dnl
define(`__IS_MIXED_ABBREV', 
  `ifelse(REAL_COMPLEX($1_type), REAL_COMPLEX($2_type), `f', `t')')dnl
dnl
dnl
define(`IS_MIXED_PREC', `ifelse(
  `$#', `2', `ifelse(PREC($1)PREC($2), PREC($1)PREC($1), `f', `t')', 
  `ifelse(PREC($1)PREC($2)PREC($3), PREC($1)PREC($1)PREC($1), `f', `t')')')dnl
dnl
dnl
dnl RANDOM(x, x_type, round_to_single)
dnl Set x to a random number.  If round_to_single is t, round the number to single.
define(`RANDOM', `ifelse(
  REAL_COMPLEX($2), `real', `$1 = ifelse(`$3', `t', `(float)') xrand(seed);', 
  `$1[0] = ifelse(`$3', `t', `(float)') xrand(seed);
   $1[1] = ifelse(`$3', `t', `(float)') xrand(seed);')')dnl
dnl
dnl
define(`LOWER_PREC', `ifelse(
  `$#', `1', `$1', 
  `$#', `2', `ifelse(
    `$1', `D', `$2',
    `$2', `D', `$1', `S')', 
  `LOWER_PREC($1, LOWER_PREC(shift($@)))')')dnl
dnl
dnl
define(`SET_EPS', 
  `$2 = power(2, -BITS_`'PREC($1));')dnl
dnl
dnl
define(`SET_UN', 
  `$2 = pow((double) BLAS_fpinfo_x(blas_base, BLAS_PREC($1)), 
            (double) BLAS_fpinfo_x(blas_emin, BLAS_PREC($1)));')dnl
dnl
dnl
dnl SET_INTERNAL_PARAMS(type, _x)
dnl If called with one argument, sets underflow and precision of given 
dnl type.  If _x is present, then sets the appropriate underflow and 
dnl precision values based on prec_val variable.
define(`SET_INTERNAL_PARAMS', 
  `ifelse(`$2', `_x', 
    `switch(prec_val){
       case 0:
         SET_EPS($1, eps_int)
         SET_UN($1, un_int)
         prec = BLAS_PREC($1);
         break;
       case 1: 
         SET_EPS(real_D, eps_int)
         SET_UN(real_D, un_int)
         prec = BLAS_PREC(real_D);
         break;
       case 2: 
       default:
         SET_EPS(real_E, eps_int)
         SET_UN(real_E, un_int)
         prec = BLAS_PREC(real_E);
         break;
     }',
    `SET_EPS($1, eps_int)
     SET_UN($1, un_int)
     prec = BLAS_PREC($1);')')dnl
dnl
dnl
dnl PRINT_NUMBER(var, type) ... print value of var
define(`PRINT_NUMBER', `ifelse(
  `$2', `real_S', `printf("%16.8e", $1);', 
  `$2', `real_D', `printf("%24.16e", $1);', 
  `$2', `complex_S', `printf("(%16.8e, %16.8e)", $1[0], $1[1]);',
  `$2', `complex_D', `printf("(%24.16e, %24.16e)", $1[0], $1[1]);',
  `$2', `real_E', `printf("[%24.16e %24.16e]", HEAD($1), TAIL($1));',
  `$2', `complex_E', 
    `printf("([%24.16e  %24.16e], [%24.16e %24.16e])", 
      HEAD($1)[0], TAIL($1)[0], HEAD($1)[1], TAIL($1)[1]);', `
#error Unknown type for PRINT_NUMBER')')dnl
dnl
dnl
dnl ---------------------------------------------------------------------
dnl Usage: PRINT_VAR(var, type, name)
dnl        print var
dnl
dnl        type : the type and precision var
dnl The type and precision specifier can be one of
dnl        real_S       ... real and single
dnl        real_D       ... real and double
dnl        real_I       ... real and indigenous
dnl        real_E       ... real and extra
dnl        complex_S    ... complex and single
dnl        complex_D    ... complex and double
dnl        complex_I    ... complex and indigeneous
dnl        complex_E    ... complex and extra
dnl ----------------------------------------------------------------------
define(`PRINT_VAR', 
  `ifelse(`$3', `', `printf("$1 = "); ', `printf("$3 = "); ')PRINT_NUMBER($1, $2)')dnl
dnl
dnl
dnl ---------------------------------------------------------------------
dnl Usage: PRINT_ARRAY_ELEM(var, index, type, name)
dnl        print var[index]
dnl
dnl        type : the type and precision var
dnl The type and precision specifier can be one of
dnl        real_S       ... real and single
dnl        real_D       ... real and double
dnl        real_I       ... real and indigenous
dnl        real_E       ... real and extra
dnl        complex_S    ... complex and single
dnl        complex_D    ... complex and double
dnl        complex_I    ... complex and indigeneous
dnl        complex_E    ... complex and extra
dnl ----------------------------------------------------------------------
define(`PRINT_ARRAY_ELEM',`ifelse(`$4', `', `', `printf("$4[%d] = ", $2); ')ifelse(
  `$3', `real_S', `printf("%16.8e", $1[$2]);', 
  `$3', `real_D', `printf("%24.16e", $1[$2]);', 
  `$3', `complex_S', `printf("(%16.8e, %16.8e)", $1[$2], $1[$2+1]);',
  `$3', `complex_D', `printf("(%24.16e, %24.16e)", $1[$2], $1[$2+1]);',
  `$3', `real_E', `printf("[%24.16e, %24.16e]", HEAD($1)[$2], TAIL($1)[$2]);', 
  `$3', `complex_E', 
    `printf("([%24.16e  %24.16e], [%24.16e %24.16e])", 
      HEAD($1)[$2], TAIL($1)[$2], HEAD($1)[$2+1], TAIL($1)[$2+1]);', `
#error Unknown type for PRINT_ARRAY_ELEM')')dnl
dnl
dnl
dnl PRINT_PREC(prec)
define(`PRINT_PREC', 
  `switch($1) {
   case blas_prec_single:     printf("single "); break;
   case blas_prec_double:     printf("double "); break;
   case blas_prec_indigenous: printf("indigenous "); break;
   case blas_prec_extra:      printf("extra "); break;
}')dnl
dnl
dnl
dnl PRINT_TRANS(name, var)
define(`PRINT_TRANS', 
  `ifelse(`$2', `', `', `printf("$2:");
   ')dnl
   switch($1) {
   case blas_no_trans:   printf("no_trans "); break;
   case blas_trans:      printf("trans "); break;
   case blas_conj_trans: printf("conj_trans "); break;
}')dnl
dnl
dnl
dnl PRINT_NORM(var)
define(`PRINT_NORM', 
 `switch ($1) {
  case -1: printf("near_underflow "); break;
  case  0: printf("near_one "); break;
  case  1: printf("near_overflow "); break;
}')dnl
dnl
dnl
dnl PRINT_ORDER(var)
define(`PRINT_ORDER', 
 `switch ($1) {
  case blas_rowmajor:
       printf("row_major "); break;
  case blas_colmajor:
       printf("col_major "); break;
}')dnl
dnl
dnl
dnl PRINT_UPLO(var)
define(`PRINT_UPLO', 
 `switch($1) {
    case blas_upper: printf("upper "); break;
    case blas_lower: printf("lower "); break;
}')dnl
dnl
dnl
dnl PRINT_DIAG(var)
define(`PRINT_DIAG', 
 `switch($1) {
    case blas_non_unit_diag: printf("non_unit_diag "); break;
    case blas_unit_diag:     printf("unit_diag "); break;
}')dnl
dnl
dnl
dnl PRINT_CONJ(var)
define(`PRINT_CONJ', 
 `switch($1){
  case blas_no_conj: printf("no_conj "); break;
  case blas_conj: printf("conj "); break;
}')dnl
dnl
dnl
dnl SET_ALPHA(type)
dnl Sets alpha to zero, one, or leave it alone depending on alpha_val
define(`SET_ALPHA', 
 `alpha_flag = 0;
  switch (alpha_val) {
  case 0: 
    ZERO(alpha, $1)
    alpha_flag = 1;
    break;
  case 1: 
    ONE(alpha, $1)
    alpha_flag = 1;
    break;
  }')dnl
dnl
dnl
dnl SET_BETA(type)
dnl Sets beta to zero, one, or leave it alone depending on beta_val
define(`SET_BETA', 
 `beta_flag = 0;
  switch (beta_val) {
  case 0: 
    ZERO(beta, $1)
    beta_flag = 1;
    break;
  case 1: 
    ONE(beta, $1)
    beta_flag = 1;
    break;
  }')dnl
dnl
dnl
dnl DOT_TESTGEN_NAME(abr_typeltr, x_typeltr, y_typeltr)
define(`DOT_TESTGEN_NAME', `ifelse(
  `$2&&$3', `$1&&$1', `BLAS_$1dot_testgen', `BLAS_$1dot_$2_$3_testgen')')dnl
dnl
dnl
dnl DOT2_TESTGEN_NAME(abr_typeltr, x_typeltr, y_typeltr)
define(`DOT2_TESTGEN_NAME', `ifelse(
  `$2&&$3', `$1&&$1', `BLAS_$1dot2_testgen', `BLAS_$1dot2_$2_$3_testgen')')dnl
dnl
dnl
dnl --------------------------------------------------------------------
dnl Usage: TEST_DOT_NAME(abr_typeltr, x_typeltr, y_typeltr, extended) 
dnl        create a test_dot name 
dnl --------------------------------------------------------------------
define(`TEST_DOT_NAME', `ifelse(
        `$2&&$3', `$1&&$1',
        `test_BLAS_$1dot',
        `test_BLAS_$1dot_$2_$3')')dnl
dnl
dnl
dnl --------------------------------------------------------------------
dnl Usage: TEST_DOT2_NAME(aby_typeltr, A_typeltr, x_typeltr, extended) 
dnl        create a test_dot2 name 
dnl --------------------------------------------------------------------
define(`TEST_DOT2_NAME', `ifelse(
        `$2&&$3', `$1&&$1',
        `test_BLAS_$1dot2',
        `test_BLAS_$1dot2_$2_$3')')dnl
dnl
dnl
dnl RAND_PTR(a, a_type)  ... *a = random
define(`RAND_PTR', `ifelse(
  `$2', `complex_S', `$1[0] = xrand(seed);  $1[1] = xrand(seed);', 
  `$2', `complex_D', `$1[0] = xrand(seed);  $1[1] = xrand(seed);', 
  `*$1 = xrand(seed);')')dnl
dnl
dnl
define(`PRINT_TEST_RESULT', `
  printf("\n");
  if (nr_failed_routines)
    printf("FAILED ");
  else
    printf("PASSED ");
  printf("%-10s: FAIL/TOTAL = %d/%d\n", 
         base_routine, nr_failed_routines, nr_routines);
')dnl
dnl
dnl
define(`MAIN_TOP', `
int main(int argc, char **argv) {
  int nsizes, ntests, debug;
  double thresh, test_prob;
  double total_min_ratio, total_max_ratio;
  int total_bad_ratios;
  int seed, num_bad_ratio, num_tests;
  int total_tests, nr_failed_routines = 0, nr_routines = 0;
  double min_ratio, max_ratio;
  const char *base_routine = "BASE_ROUTINE()";
  char *fname;
  int n;
  $1

  if (argc != 6) {
    printf("Usage:\n");
    printf("do_test_`'BASE_ROUTINE() <nsizes> <ntests> <thresh> <debug> <test_prob>\n");
    printf("   <nsizes>: number of sizes to be run.\n");
    printf("   <ntests>: the number of tests performed for each set of attributes\n");
    printf("   <thresh>: to catch bad ratios if it is greater than <thresh>\n");
    printf("    <debug>: 0, 1, 2, or 3; \n");
    printf("        if 0, no printing \n");
    printf("        if 1, print error summary only if tests fail\n");
    printf("        if 2, print error summary for each n\n");
    printf("        if 3, print complete info each test fails \n");
    printf("<test_prob>: probability of preforming a given \n");
    printf("           test case: 0.0 does no tests, 1.0 does all tests\n");
    return -1;
  } else {
    nsizes= atoi(argv[1]);
    ntests = atoi(argv[2]);
    thresh = atof(argv[3]);
    debug = atoi(argv[4]);
    test_prob = atof(argv[5]);
  }

  seed = 1999;

  if (nsizes<0 || ntests<0 || debug<0 || debug>3)
    BLAS_error("Testing BASE_ROUTINE()", 0, 0, NULL);

  printf("Testing %s...\n", base_routine);
  printf("INPUT: nsizes = %d, ntests = %d, thresh = %4.2f, debug = %d\n\n",
        nsizes, ntests, thresh, debug);  

  ')dnl
dnl
dnl
define(`MAIN_BOTTOM', `
  PRINT_TEST_RESULT
  return 0;
}
')dnl
dnl
dnl
dnl MAIN(preamble, body)
define(`MAIN', 
`MAIN_TOP(`$1')
  $2
MAIN_BOTTOM()')dnl
