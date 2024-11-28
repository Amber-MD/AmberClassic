dnl
dnl
define(`GEMM_NAME', 
  `ifelse(`$2&&$3', `$1&&$1', `BLAS_$1gemm$4', `BLAS_$1gemm_$2_$3$4')')dnl
dnl
dnl
define(`GEMM_PARAMS', 
  `enum blas_order_type order, enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   $1_scalar alpha, const $2_array a, int lda, const $3_array b, int ldb,
   $1_scalar beta, $1_array c, int ldc`'ifelse(`$4', `_x', `, enum blas_prec_type prec')')dnl
dnl
dnl
define(`GEMM_HEAD', 
  `void GEMM_NAME($1, $2, $3, $4)(GEMM_PARAMS($1, $2, $3, $4))')dnl
dnl
dnl
define(`GEMM_ARGS', 
`if_blas(``s, s, s', `d, d, d', 
 `c, c, c', `z, z, z',')dnl
`d, d, s', `d, s, d', 
 `d, s, s', `z, z, c', 
 `z, c, z', `z, c, c', 
 `c, c, s', `c, s, c', 
 `c, s, s', `z, z, d', 
 `z, d, z', `z, d, d', 
 `s, s, s, _x', `d, d, d, _x', 
 `c, c, c, _x', `z, z, z, _x', 
 `d, d, s, _x', `d, s, d, _x', 
 `d, s, s, _x', `z, z, c, _x', 
 `z, c, z, _x', `z, c, c, _x', 
 `c, c, s, _x', `c, s, c, _x', 
 `c, s, s, _x', `z, z, d, _x', 
 `z, d, z, _x', `z, d, d, _x'')dnl
define(`GEMM_AMB_ARGS', 
`if_blas(``d, d, d',')dnl
`d, d, d, _x'')dnl
dnl
dnl
