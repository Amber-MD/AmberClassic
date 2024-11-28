dnl
dnl
include(cblas.m4)dnl
dnl
dnl
#ifndef BLAS_EXTENDED_PROTO_H
#define BLAS_EXTENDED_PROTO_H

FOREACH(`LIBNAMES', `
include(arg/arg-common.m4)dnl
define(`ARG', translit(arg, `abcdefghijklmnopqrstuvwxyz', 
`ABCDEFGHIJKLMNOPQRSTUVWXYZ'))dnl
define(`ROUTINE', ARG)dnl
define(`routine', arg)dnl
FOREACH(`ROUTINE()_ARGS', `ROUTINE()_HEAD(arg);
')
')

int BLAS_fpinfo_x(enum blas_cmach_type cmach, enum blas_prec_type prec);
void BLAS_error(const char *rname, int iflag, int ival, char *form, ...);

#endif /* BLAS_EXTENDED_PROTO_H */

