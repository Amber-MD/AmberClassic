dnl
dnl
define(`AXPBY_NAME', 
  `ifelse( `$1', `$2', `BLAS_$1axpby$3', `BLAS_$1axpby_$2$3')')dnl
dnl
dnl
define(`AXPBY_PARAMS',
  `int n, $1_scalar alpha, const $2_array x, int incx, 
   $1_scalar beta, $1_array y, 
   int incy`'ifelse(`$3', `_x', `, enum blas_prec_type prec')')dnl
dnl
dnl
define(`AXPBY_HEAD', 
  `void AXPBY_NAME($1, $2, $3)(AXPBY_PARAMS($1, $2, $3))')dnl
dnl
dnl
define(`AXPBY_ARGS', 
 `if_blas(``s, s', `d, d', `c, c', `z, z', ')`d, s', `c, s', `z, c', `z, d',
  `s, s, _x', `d, d, _x', `c, c, _x', `z, z, _x',
  `d, s, _x', `z, c, _x', `c, s, _x', `z, d, _x'')dnl
define(`AXPBY_AMB_ARGS', 
 `if_blas(``d, d', ')`d, d, _x'')dnl
dnl
dnl
