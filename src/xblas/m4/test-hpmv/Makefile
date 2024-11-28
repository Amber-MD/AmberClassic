include ../../make.conf
include ../../$(MAKEINC)

SRC_PATH = ../../testing/test-hpmv
M4_PATH = -I. -I..

TEST_SRCS = \
  $(SRC_PATH)/do_test_hpmv.c\
  $(SRC_PATH)/BLAS_hpmv_testgen.c\
  $(SRC_PATH)/hpmv-support.c

test-source: $(TEST_SRCS)


$(SRC_PATH)/do_test_hpmv.c: do_test_hpmv.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=hpmv do_test_hpmv.m4 > do_test_hpmv.tmp.c && \
	$(INDENT) $(INDENT_OPTS) do_test_hpmv.tmp.c && \
	mv do_test_hpmv.tmp.c $@ && rm -f do_test_hpmv.tmp.c*

$(SRC_PATH)/BLAS_hpmv_testgen.c: BLAS_hpmv_testgen.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=hpmv BLAS_hpmv_testgen.m4 > BLAS_hpmv_testgen.tmp.c && \
	$(INDENT) $(INDENT_OPTS) BLAS_hpmv_testgen.tmp.c && \
	mv BLAS_hpmv_testgen.tmp.c $@ && rm -f BLAS_hpmv_testgen.tmp.c*

$(SRC_PATH)/hpmv-support.c: hpmv-support.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=hpmv hpmv-support.m4 > hpmv-support.tmp.c && \
	$(INDENT) $(INDENT_OPTS) hpmv-support.tmp.c && \
	mv hpmv-support.tmp.c $@ && rm -f hpmv-support.tmp.c*


source-clean:
	rm -f $(TEST_SRCS)

