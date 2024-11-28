dnl
dnl --------------------------------------------------------------------
dnl Usage: PACKED_ADDRESS(order, uplo, AP, row, col, size) ... compute 
dnl        the address of AP[row, column].  Does not check that you are
dnl        in the right part the matrix.  This routine is very expensive 
dnl        and dangerous and should not be used.  It is a function, 
dnl        returning a mathematical expression.
dnl --------------------------------------------------------------------
define(`PACKED_ADDRESS', `ifelse(
	`$1&&$2', `blas_rowmajor&&blas_upper', 
	`($3 + ($5 + $4 * ($4 - 1) / 2))',
	`$1&&$2', `blas_colmajor&&blas_upper', 
	`($3 + ($4 + $5 * ($5 - 1) / 2))',
	`$1&&$2', `blas_rowmajor&&blas_lower', 
	`($3 + ($5 + (2 * $6 - $4) * ($4 - 1) / 2))',
	`$1&&$2', `blas_colmajor&&blas_lower', 
	`($3 + ($4 + (2 * $6 - $5)*($5 - 1) / 2))')')dnl
dnl
dnl --------------------------------------------------------------------
dnl Usage: COLUMN_STRIDE(index, size, column, order, uplo, inc) ... 
dnl        Advances the index one column to the right in a packed array.
dnl        In the column major case, this is a function of the present
dnl        column.
dnl --------------------------------------------------------------------
dnl
define(`COLUMN_STRIDE', `ifelse(
	`$4&&$5', `blas_colmajor&&blas_lower', `$1 += ($2 - $3 - 1) * $6;',
	`$4&&$5', `blas_colmajor&&blas_upper', `$1 += ($3 + 1) * $6;',
	`$1 += $6;')')dnl
dnl
dnl --------------------------------------------------------------------
dnl Usage: ROW_STRIDE(index, size, row, order, uplo, inc) ... Advances
dnl        the index one row down in a packed array.  In the
dnl        row major case, this is a function of the present row.
dnl --------------------------------------------------------------------
dnl
define(`ROW_STRIDE', `ifelse(
	`$4&&$5', `blas_rowmajor&&blas_upper', `$1 += ($2 - $3 - 1) * $6;',
	`$4&&$5', `blas_rowmajor&&blas_lower', `$1 += ($3 + 1) * $6;',
	`$1 += $6;')')dnl
dnl
dnl The macros below are general usage and might be useful in cblas.m4.h
dnl
dnl
dnl --------------------------------------------------------------------
dnl Usage: FLIP(flip_enum, variable) ... computes the reversed value of 
dnl        any two-value operator argument.  flip_enum is one of the
dnl        enumerated types in cblas.h.  
dnl --------------------------------------------------------------------
define(`FLIP', `ifelse(
	`$1', `blas_sort_type', 
    `($2 == blas_increasing_order) ? blas_decreasing_order : blas_increasing_order',
	`$1', `blas_side_type', 
    `($2 == blas_left_side) ? blas_right_side : blas_left_side',
	`$1', `blas_uplo_type',
    `($2 == blas_upper) ? blas_lower : blas_upper',
	`$1', `blas_conj_type', 
    `($2 == blas_conj) ? blas_no_conj : blas_conj',
	`$1', `blas_diag_type', 
    `($2 == blas_unit_diag) ? blas_non_unit_diag : blas_unit_diag', 
	`$1', `blas_direct_type', 
    `($2 == blas_forward_seq) ? blas_backward_seq : blas_forward_seq',
	`$1', `blas_order_type', 
    `($2 == blas_rowmajor) ? blas_colmajor : blas_rowmajor')')dnl
dnl
dnl
