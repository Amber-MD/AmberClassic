dnl
dnl
define(`TRSV_NAME', 
  `ifelse(`$2&&$2', `$1&&$1', `BLAS_$1trsv$3', `BLAS_$1trsv_$2$3')')dnl
dnl
dnl
define(`TRSV_PARAMS', 
  `enum blas_order_type order, enum blas_uplo_type uplo,
   enum blas_trans_type trans, enum blas_diag_type diag, 
   int n, $1_scalar alpha, const $2_array T, int ldt, 
   $1_array x, int incx`'ifelse($3, _x, `, enum blas_prec_type prec')')dnl
dnl
dnl
define(`TRSV_HEAD', 
  `void TRSV_NAME($1, $2, $3)(TRSV_PARAMS($1, $2, $3))')dnl
dnl
dnl
define(`TRSV_ARGS', 
`if_blas(``s, s', `d, d', `c, c', `z, z', ')`d, s', `z, c', `c, s', `z, d', 
 `s, s, _x', `d, d, _x', `d, s, _x', `c, c, _x', 
 `z, z, _x', `z, c, _x', `c, s, _x', `z, d, _x'')dnl
dnl
dnl
