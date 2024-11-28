dnl
dnl
define(`SUM_NAME', `BLAS_$1sum$2')dnl
dnl
dnl
define(`SUM_PARAMS', 
  `int n, const $1_array x, int incx, 
   $1_array sum`'ifelse(`$2', `_x', `, enum blas_prec_type prec')')dnl
dnl
dnl
define(`SUM_HEAD', 
  `void SUM_NAME($1, $2)(SUM_PARAMS($1, $2))')dnl
dnl
dnl
define(`SUM_ARGS', 
`if_blas(``s', `d', `c', `z',')`s, _x', `d, _x', `c, _x', `z, _x'')dnl
dnl
dnl
