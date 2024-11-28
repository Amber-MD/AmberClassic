dnl
dnl
define(`HEMV2_NAME', `ifelse(
  `$2&&$3', `$1&&$1', `BLAS_$1hemv2$4', `BLAS_$1hemv2_$2_$3$4')')dnl
dnl
dnl
define(`HEMV2_PARAMS', 
  `enum blas_order_type order, enum blas_uplo_type uplo, int n,
   $1_scalar alpha, const $2_array a, int lda, const $3_array x_head,
   const $3_array x_tail, int incx, $1_scalar beta, const $3_array y, 
   int incy`'ifelse(`$4', `_x', `, enum blas_prec_type prec')')dnl
dnl
dnl
define(`HEMV2_HEAD', 
  `void HEMV2_NAME($1, $2, $3, $4)(HEMV2_PARAMS($1, $2, $3, $4))')dnl
dnl
dnl
define(`HEMV2_ARGS', 
  `if_blas(``c, c, c', `z, z, z', ')`z, z, c', `z, c, z',
  `z, c, c', `c, c, s', `z, z, d',
  `c, c, c, _x', `z, z, z, _x',
  `z, z, c, _x', `z, c, z, _x',
  `z, c, c, _x', `c, c, s, _x', `z, z, d, _x'')dnl
dnl
dnl
