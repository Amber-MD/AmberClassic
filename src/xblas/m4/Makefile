include ../make.conf
include ../$(MAKEINC)
MKFILES = Makefile \
   \
  dot/Makefile test-dot/Makefile \
  sum/Makefile test-sum/Makefile \
  axpby/Makefile test-axpby/Makefile \
  waxpby/Makefile test-waxpby/Makefile \
  gemv/Makefile test-gemv/Makefile \
  ge_sum_mv/Makefile test-ge_sum_mv/Makefile \
  gbmv/Makefile test-gbmv/Makefile \
  symv/Makefile test-symv/Makefile \
  spmv/Makefile test-spmv/Makefile \
  sbmv/Makefile test-sbmv/Makefile \
  hemv/Makefile test-hemv/Makefile \
  hpmv/Makefile test-hpmv/Makefile \
  hbmv/Makefile test-hbmv/Makefile \
  trmv/Makefile test-trmv/Makefile \
  tpmv/Makefile test-tpmv/Makefile \
  trsv/Makefile test-trsv/Makefile \
  tbsv/Makefile test-tbsv/Makefile \
  gemm/Makefile test-gemm/Makefile \
  symm/Makefile test-symm/Makefile \
  hemm/Makefile test-hemm/Makefile \
  gemv2/Makefile test-gemv2/Makefile \
  symv2/Makefile test-symv2/Makefile \
  hemv2/Makefile test-hemv2/Makefile \
  gbmv2/Makefile test-gbmv2/Makefile \
  test-dot2/Makefile

HEADERS = ../src/blas_extended_proto.h ../testing/blas_extended_test.h

makefiles: $(MKFILES)

header: $(HEADERS)

sources: 
	cd dot && $(MAKE) dot-source
	cd sum && $(MAKE) sum-source
	cd axpby && $(MAKE) axpby-source
	cd waxpby && $(MAKE) waxpby-source
	cd gemv && $(MAKE) gemv-source
	cd ge_sum_mv && $(MAKE) ge_sum_mv-source
	cd gbmv && $(MAKE) gbmv-source
	cd symv && $(MAKE) symv-source
	cd spmv && $(MAKE) spmv-source
	cd sbmv && $(MAKE) sbmv-source
	cd hemv && $(MAKE) hemv-source
	cd hpmv && $(MAKE) hpmv-source
	cd hbmv && $(MAKE) hbmv-source
	cd trmv && $(MAKE) trmv-source
	cd tpmv && $(MAKE) tpmv-source
	cd trsv && $(MAKE) trsv-source
	cd tbsv && $(MAKE) tbsv-source
	cd gemm && $(MAKE) gemm-source
	cd symm && $(MAKE) symm-source
	cd hemm && $(MAKE) hemm-source
	cd gemv2 && $(MAKE) gemv2-source
	cd symv2 && $(MAKE) symv2-source
	cd hemv2 && $(MAKE) hemv2-source
	cd gbmv2 && $(MAKE) gbmv2-source

test-sources: 
	cd test-dot && $(MAKE) test-source
	cd test-sum && $(MAKE) test-source
	cd test-axpby && $(MAKE) test-source
	cd test-waxpby && $(MAKE) test-source
	cd test-gemv && $(MAKE) test-source
	cd test-ge_sum_mv && $(MAKE) test-source
	cd test-gbmv && $(MAKE) test-source
	cd test-symv && $(MAKE) test-source
	cd test-spmv && $(MAKE) test-source
	cd test-sbmv && $(MAKE) test-source
	cd test-hemv && $(MAKE) test-source
	cd test-hpmv && $(MAKE) test-source
	cd test-hbmv && $(MAKE) test-source
	cd test-trmv && $(MAKE) test-source
	cd test-tpmv && $(MAKE) test-source
	cd test-trsv && $(MAKE) test-source
	cd test-tbsv && $(MAKE) test-source
	cd test-gemm && $(MAKE) test-source
	cd test-symm && $(MAKE) test-source
	cd test-hemm && $(MAKE) test-source
	cd test-gemv2 && $(MAKE) test-source
	cd test-symv2 && $(MAKE) test-source
	cd test-hemv2 && $(MAKE) test-source
	cd test-gbmv2 && $(MAKE) test-source
	cd test-dot2 && $(MAKE) test-source

Makefile: Makefile.m4 cblas.m4
	$(M4) $(M4_OPTS) Makefile.m4 >Makefile

../src/blas_extended_proto.h: Makefile.m4 proto.m4\
	dot/dot-common.m4\
	sum/sum-common.m4\
	axpby/axpby-common.m4\
	waxpby/waxpby-common.m4\
	gemv/gemv-common.m4\
	ge_sum_mv/ge_sum_mv-common.m4\
	gbmv/gbmv-common.m4\
	symv/symv-common.m4\
	spmv/spmv-common.m4\
	sbmv/sbmv-common.m4\
	hemv/hemv-common.m4\
	hpmv/hpmv-common.m4\
	hbmv/hbmv-common.m4\
	trmv/trmv-common.m4\
	tpmv/tpmv-common.m4\
	trsv/trsv-common.m4\
	tbsv/tbsv-common.m4\
	gemm/gemm-common.m4\
	symm/symm-common.m4\
	hemm/hemm-common.m4\
	gemv2/gemv2-common.m4\
	symv2/symv2-common.m4\
	hemv2/hemv2-common.m4\
	gbmv2/gbmv2-common.m4
	$(M4) $(M4_OPTS) proto.m4 > blas_extended_proto.tmp.h && \
        $(INDENT) $(INDENT_OPTS) blas_extended_proto.tmp.h && \
        mv blas_extended_proto.tmp.h $@ && rm -f blas_extended_proto.tmp.h*

../testing/blas_extended_test.h: Makefile.m4 test-proto.m4
	$(M4) $(M4_OPTS) test-proto.m4 > blas_extended_test.tmp.h && \
        $(INDENT) $(INDENT_OPTS) blas_extended_test.tmp.h && \
        mv blas_extended_test.tmp.h $@ && rm -f blas_extended_test.tmp.h*


dot/Makefile: lib-makefile.template cblas.m4 dot/dot-common.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=dot lib-makefile.template >$@

sum/Makefile: lib-makefile.template cblas.m4 sum/sum-common.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=sum lib-makefile.template >$@

axpby/Makefile: lib-makefile.template cblas.m4 axpby/axpby-common.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=axpby lib-makefile.template >$@

waxpby/Makefile: lib-makefile.template cblas.m4 waxpby/waxpby-common.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=waxpby lib-makefile.template >$@

gemv/Makefile: lib-makefile.template cblas.m4 gemv/gemv-common.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=gemv lib-makefile.template >$@

ge_sum_mv/Makefile: lib-makefile.template cblas.m4 ge_sum_mv/ge_sum_mv-common.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=ge_sum_mv lib-makefile.template >$@

gbmv/Makefile: lib-makefile.template cblas.m4 gbmv/gbmv-common.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=gbmv lib-makefile.template >$@

symv/Makefile: lib-makefile.template cblas.m4 symv/symv-common.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=symv lib-makefile.template >$@

spmv/Makefile: lib-makefile.template cblas.m4 spmv/spmv-common.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=spmv lib-makefile.template >$@

sbmv/Makefile: lib-makefile.template cblas.m4 sbmv/sbmv-common.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=sbmv lib-makefile.template >$@

hemv/Makefile: lib-makefile.template cblas.m4 hemv/hemv-common.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=hemv lib-makefile.template >$@

hpmv/Makefile: lib-makefile.template cblas.m4 hpmv/hpmv-common.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=hpmv lib-makefile.template >$@

hbmv/Makefile: lib-makefile.template cblas.m4 hbmv/hbmv-common.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=hbmv lib-makefile.template >$@

trmv/Makefile: lib-makefile.template cblas.m4 trmv/trmv-common.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=trmv lib-makefile.template >$@

tpmv/Makefile: lib-makefile.template cblas.m4 tpmv/tpmv-common.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=tpmv lib-makefile.template >$@

trsv/Makefile: lib-makefile.template cblas.m4 trsv/trsv-common.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=trsv lib-makefile.template >$@

tbsv/Makefile: lib-makefile.template cblas.m4 tbsv/tbsv-common.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=tbsv lib-makefile.template >$@

gemm/Makefile: lib-makefile.template cblas.m4 gemm/gemm-common.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=gemm lib-makefile.template >$@

symm/Makefile: lib-makefile.template cblas.m4 symm/symm-common.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=symm lib-makefile.template >$@

hemm/Makefile: lib-makefile.template cblas.m4 hemm/hemm-common.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=hemm lib-makefile.template >$@

gemv2/Makefile: lib-makefile.template cblas.m4 gemv2/gemv2-common.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=gemv2 lib-makefile.template >$@

symv2/Makefile: lib-makefile.template cblas.m4 symv2/symv2-common.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=symv2 lib-makefile.template >$@

hemv2/Makefile: lib-makefile.template cblas.m4 hemv2/hemv2-common.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=hemv2 lib-makefile.template >$@

gbmv2/Makefile: lib-makefile.template cblas.m4 gbmv2/gbmv2-common.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=gbmv2 lib-makefile.template >$@


test-sum/Makefile: test-makefile.template cblas.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -D EXTRA_SRCS=BLAS_sum_testgen,sum-support -Droutine=sum test-makefile.template >$@

test-gbmv/Makefile: test-makefile.template cblas.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -D EXTRA_SRCS=BLAS_gbmv_testgen,gbmv-support -Droutine=gbmv test-makefile.template >$@

test-gemv/Makefile: test-makefile.template cblas.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -D EXTRA_SRCS=BLAS_gemv_testgen,gemv-support -Droutine=gemv test-makefile.template >$@

test-hbmv/Makefile: test-makefile.template cblas.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -D EXTRA_SRCS=BLAS_hbmv_testgen,hbmv-support -Droutine=hbmv test-makefile.template >$@

test-hemv/Makefile: test-makefile.template cblas.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -D EXTRA_SRCS=BLAS_hemv_testgen,hemv-support -Droutine=hemv test-makefile.template >$@

test-hpmv/Makefile: test-makefile.template cblas.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -D EXTRA_SRCS=BLAS_hpmv_testgen,hpmv-support -Droutine=hpmv test-makefile.template >$@

test-sbmv/Makefile: test-makefile.template cblas.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -D EXTRA_SRCS=BLAS_sbmv_testgen,sbmv-support -Droutine=sbmv test-makefile.template >$@

test-tbsv/Makefile: test-makefile.template cblas.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -D EXTRA_SRCS=BLAS_tbsv_testgen,tbsv-support -Droutine=tbsv test-makefile.template >$@

test-spmv/Makefile: test-makefile.template cblas.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -D EXTRA_SRCS=BLAS_spmv_testgen,spmv-support -Droutine=spmv test-makefile.template >$@

test-symv/Makefile: test-makefile.template cblas.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -D EXTRA_SRCS=BLAS_symv_testgen,symv-support -Droutine=symv test-makefile.template >$@

test-tpmv/Makefile: test-makefile.template cblas.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -D EXTRA_SRCS=BLAS_tpmv_testgen,tpmv-support -Droutine=tpmv test-makefile.template >$@

test-trmv/Makefile: test-makefile.template cblas.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -D EXTRA_SRCS=BLAS_trmv_testgen,trmv-support -Droutine=trmv test-makefile.template >$@

test-trsv/Makefile: test-makefile.template cblas.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -D EXTRA_SRCS=BLAS_trsv_testgen,trsv-support -Droutine=trsv test-makefile.template >$@


test-axpby/Makefile: test-makefile.template cblas.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=axpby test-makefile.template >$@

test-waxpby/Makefile: test-makefile.template cblas.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=waxpby test-makefile.template >$@


test-ge_sum_mv/Makefile: test-makefile.template cblas.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -D EXTRA_SRCS=BLAS_ge_sum_mv_testgen -Droutine=ge_sum_mv test-makefile.template >$@

test-gemv2/Makefile: test-makefile.template cblas.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -D EXTRA_SRCS=BLAS_gemv2_testgen -Droutine=gemv2 test-makefile.template >$@

test-symv2/Makefile: test-makefile.template cblas.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -D EXTRA_SRCS=BLAS_symv2_testgen -Droutine=symv2 test-makefile.template >$@

test-hemv2/Makefile: test-makefile.template cblas.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -D EXTRA_SRCS=BLAS_hemv2_testgen -Droutine=hemv2 test-makefile.template >$@

test-gbmv2/Makefile: test-makefile.template cblas.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -D EXTRA_SRCS=BLAS_gbmv2_testgen -Droutine=gbmv2 test-makefile.template >$@

test-gemm/Makefile: test-makefile.template cblas.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -D EXTRA_SRCS=BLAS_gemm_testgen -Droutine=gemm test-makefile.template >$@

test-symm/Makefile: test-makefile.template cblas.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -D EXTRA_SRCS=BLAS_symm_testgen -Droutine=symm test-makefile.template >$@

test-hemm/Makefile: test-makefile.template cblas.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -D EXTRA_SRCS=BLAS_hemm_testgen -Droutine=hemm test-makefile.template >$@


test-dot/Makefile: test-makefile.template cblas.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -D EXTRA_SRCS=BLAS_dot_testgen,test_dot,print_vector,copy_vector -Droutine=dot test-makefile.template >$@


test-dot2/Makefile: test-makefile.template cblas.m4 Makefile ../$(MAKEINC)
	$(M4) $(M4_OPTS) -D EXTRA_SRCS=BLAS_dot2_testgen,r_truth2,dot2,test_dot2 -Droutine=dot2 test-makefile.template >$@

clean:
	rm -f temp *~

source-clean: clean 
	cd dot && $(MAKE) source-clean
	cd test-dot && $(MAKE) source-clean
	cd sum && $(MAKE) source-clean
	cd test-sum && $(MAKE) source-clean
	cd axpby && $(MAKE) source-clean
	cd test-axpby && $(MAKE) source-clean
	cd waxpby && $(MAKE) source-clean
	cd test-waxpby && $(MAKE) source-clean
	cd gemv && $(MAKE) source-clean
	cd test-gemv && $(MAKE) source-clean
	cd ge_sum_mv && $(MAKE) source-clean
	cd test-ge_sum_mv && $(MAKE) source-clean
	cd gbmv && $(MAKE) source-clean
	cd test-gbmv && $(MAKE) source-clean
	cd symv && $(MAKE) source-clean
	cd test-symv && $(MAKE) source-clean
	cd spmv && $(MAKE) source-clean
	cd test-spmv && $(MAKE) source-clean
	cd sbmv && $(MAKE) source-clean
	cd test-sbmv && $(MAKE) source-clean
	cd hemv && $(MAKE) source-clean
	cd test-hemv && $(MAKE) source-clean
	cd hpmv && $(MAKE) source-clean
	cd test-hpmv && $(MAKE) source-clean
	cd hbmv && $(MAKE) source-clean
	cd test-hbmv && $(MAKE) source-clean
	cd trmv && $(MAKE) source-clean
	cd test-trmv && $(MAKE) source-clean
	cd tpmv && $(MAKE) source-clean
	cd test-tpmv && $(MAKE) source-clean
	cd trsv && $(MAKE) source-clean
	cd test-trsv && $(MAKE) source-clean
	cd tbsv && $(MAKE) source-clean
	cd test-tbsv && $(MAKE) source-clean
	cd gemm && $(MAKE) source-clean
	cd test-gemm && $(MAKE) source-clean
	cd symm && $(MAKE) source-clean
	cd test-symm && $(MAKE) source-clean
	cd hemm && $(MAKE) source-clean
	cd test-hemm && $(MAKE) source-clean
	cd gemv2 && $(MAKE) source-clean
	cd test-gemv2 && $(MAKE) source-clean
	cd symv2 && $(MAKE) source-clean
	cd test-symv2 && $(MAKE) source-clean
	cd hemv2 && $(MAKE) source-clean
	cd test-hemv2 && $(MAKE) source-clean
	cd gbmv2 && $(MAKE) source-clean
	cd test-gbmv2 && $(MAKE) source-clean
	cd test-dot2 && $(MAKE) source-clean
	rm -f $(HEADERS)

maintainer-clean: source-clean
	rm -f $(MKFILES)

.PHONY: clean source-clean maintainer-clean makefiles header

