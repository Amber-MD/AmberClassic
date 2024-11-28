dnl M4 file for generating m4/Makefile.
dnl
include(cblas.m4)dnl
dnl
dnl
`include' ../make.conf
`include' ../$(MAKEINC)
dnl
dnl
dnl Can't use 'MAKEFILES' since it is used by GNU make.
MKFILES = Makefile \
  FOREACH(`LIBNAMES', ` \
  arg/Makefile test-arg/Makefile') \
  test-dot2/Makefile

HEADERS = ../src/blas_extended_proto.h ../testing/blas_extended_test.h

makefiles: $(MKFILES)

header: $(HEADERS)

sources: FOREACH(`LIBNAMES', `
	cd arg && $(MAKE) arg-source')

test-sources: FOREACH(`LIBNAMES', `
	cd test-arg && $(MAKE) test-source')
	cd test-dot2 && $(MAKE) test-source

Makefile: Makefile.m4 cblas.m4
	$(M4) $(M4_OPTS) Makefile.m4 >Makefile

../src/blas_extended_proto.h: Makefile.m4 proto.m4`'FOREACH(`LIBNAMES', `\
	arg/arg-common.m4')
	$(M4) $(M4_OPTS) proto.m4 > blas_extended_proto.tmp.h && \
        $(INDENT) $(INDENT_OPTS) blas_extended_proto.tmp.h && \
        mv blas_extended_proto.tmp.h $`'@ && rm -f blas_extended_proto.tmp.h*

../testing/blas_extended_test.h: Makefile.m4 test-proto.m4
	$(M4) $(M4_OPTS) test-proto.m4 > blas_extended_test.tmp.h && \
        $(INDENT) $(INDENT_OPTS) blas_extended_test.tmp.h && \
        mv blas_extended_test.tmp.h $`'@ && rm -f blas_extended_test.tmp.h*

dnl
dnl LIB_MKFILE_RULE(routine) generates rule for creating 
dnl routine/Makefile from lib-makefile.template.
define(`LIB_MKFILE_RULE', `
$1/Makefile: lib-makefile.template cblas.m4 arg/$1-common.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=$1 lib-makefile.template >$`'@
')dnl
dnl
dnl TEST_MKFILE_RULE(routine, extra_srcs) generates rule for creating
dnl test-routine/Makefile from test-makefile.template.  extra_srcs is 
dnl a comma separated list of any extra C source files generated from 
dnl M4 files with suffixes stripped (e.g., BLAS_gemv_testgen,gemv-support).
define(`TEST_MKFILE_RULE', `
test-$1/Makefile: test-makefile.template cblas.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS)ifelse(`$2', `', `', ` -D EXTRA_SRCS=$2') -Droutine=$1 test-makefile.template >$`'@
')dnl
dnl
FOREACH(`LIBNAMES', `LIB_MKFILE_RULE(arg)')dnl

FOREACH(`sum, gbmv, gemv, hbmv, hemv, 
         hpmv, sbmv, tbsv, spmv, symv, tpmv, trmv, trsv', 
        `TEST_MKFILE_RULE(arg, `BLAS_`'arg`'_testgen,arg-support')')dnl

FOREACH(`axpby, waxpby', `TEST_MKFILE_RULE(arg)')dnl

FOREACH(`ge_sum_mv, gemv2, symv2, hemv2, gbmv2, gemm, symm, hemm', `TEST_MKFILE_RULE(arg, `BLAS_`'arg`'_testgen')')dnl

TEST_MKFILE_RULE(dot, `BLAS_dot_testgen,test_dot,print_vector,copy_vector')dnl

TEST_MKFILE_RULE(dot2, `BLAS_dot2_testgen,r_truth2,dot2,test_dot2')dnl

clean:
	rm -f temp *~

source-clean: clean FOREACH(`LIBNAMES', `
	cd arg && $(MAKE) source-clean
	cd test-arg && $(MAKE) source-clean')
	cd test-dot2 && $(MAKE) source-clean
	rm -f $(HEADERS)

maintainer-clean: source-clean
	rm -f $(MKFILES)

.PHONY: clean source-clean maintainer-clean makefiles header

