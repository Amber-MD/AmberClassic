define(`ROUTINE', `translit(routine(), `abcdefghijklmnopqrstuvwxyz', `ABCDEFGHIJKLMNOPQRSTUVWXYZ')')dnl
dnl
dnl
include(cblas.m4)dnl
include(routine()-common.m4)
#include "f2c-bridge.h"
#include "blas_enum.h"
ROUTINE()_HEAD(RARG);

define(`enum',`int')dnl
define(`s_scalar', `float*')dnl
define(`d_scalar', `double*')dnl
define(`blas_order_type',`')dnl
define(`blas_trans_type',`')dnl
define(`blas_uplo_type',`')dnl
define(`blas_diag_type',`')dnl
define(`blas_side_type',`')dnl
define(`blas_cmach_type',`')dnl
define(`blas_norm_type',`')dnl
define(`blas_sort_type',`')dnl
define(`blas_conj_type',`')dnl
define(`blas_jrot_type',`')dnl
define(`blas_prec_type',`')dnl
define(`blas_base_type',`')dnl
define(`blas_symmetry_type',`')dnl
define(`blas_field_type',`')dnl
define(`blas_size_type',`')dnl
define(`blas_handle_type',`')dnl
define(`blas_sparsity_optimization_type',`')dnl
define(`int',`in`'t*')dnl
define(`plonkorder', `ifelse($#, 0, , `$1', `int*  order', `plonkorder(shift($@))', $#, 1, ``$1'', ``$1', plonkorder(shift($@))')')
extern void FC_FUNC_(translit(ROUTINE()_NAME(RARG),'A-Z','a-z'),translit(ROUTINE()_NAME(RARG),'a-z','A-Z'))
  ( plonkorder(ROUTINE()_PARAMS(RARG)) )
{
define(`order',`blas_colmajor')dnl
define(`blas_order_type',`')dnl
define(`enum',`')dnl
define(`void',`')dnl
define(`const',`')dnl
define(`blas_trans_type',`(en`'um blas_trans`'_type)*')dnl
define(`blas_uplo_type',`(en`'um blas_uplo`'_type)*')dnl
define(`blas_diag_type',`(en`'um blas_diag`'_type)*')dnl
define(`blas_side_type',`(en`'um blas_side`'_type)*')dnl
define(`blas_cmach_type',`(en`'um blas_cmach`'_type)*')dnl
define(`blas_norm_type',`(en`'um blas_norm`'_type)*')dnl
define(`blas_sort_type',`(en`'um blas_sort`'_type)*')dnl
define(`blas_conj_type',`(en`'um blas_conj`'_type)*')dnl
define(`blas_jrot_type',`(en`'um blas_jrot`'_type)*')dnl
define(`blas_prec_type',`(en`'um blas_prec`'_type)*')dnl
define(`blas_base_type',`(en`'um blas_base`'_type)*')dnl
define(`blas_symmetry_type',`(en`'um blas_symmetry`'_type)*')dnl
define(`blas_field_type',`(en`'um blas_field`'_type)*')dnl
define(`blas_size_type',`(en`'um blas_size`'_type)*')dnl
define(`blas_handle_type',`(en`'um blas_handle`'_type)*')dnl
define(`blas_sparsity_optimization_type',`(en`'um blas_sparsity_optimization`'_type)*')dnl
define(`int',`*')dnl
define(`s_scalar',`*')dnl
define(`d_scalar',`*')dnl
define(`c_scalar',`')dnl
define(`z_scalar',`')dnl
define(`s_array',`')dnl
define(`d_array',`')dnl
define(`c_array',`')dnl
define(`z_array',`')dnl
	ROUTINE()_HEAD(RARG);
}
