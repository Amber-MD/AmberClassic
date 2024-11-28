include ../../make.conf
include ../../$(MAKEINC)

SRC_PATH = ../../testing/test-spmv
M4_PATH = -I. -I..

TEST_SRCS = \
  $(SRC_PATH)/do_test_spmv.c\
  $(SRC_PATH)/BLAS_spmv_testgen.c\
  $(SRC_PATH)/spmv-support.c

test-source: $(TEST_SRCS)


$(SRC_PATH)/do_test_spmv.c: do_test_spmv.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=spmv do_test_spmv.m4 > do_test_spmv.tmp.c && \
	$(INDENT) $(INDENT_OPTS) do_test_spmv.tmp.c && \
	mv do_test_spmv.tmp.c $@ && rm -f do_test_spmv.tmp.c*

$(SRC_PATH)/BLAS_spmv_testgen.c: BLAS_spmv_testgen.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=spmv BLAS_spmv_testgen.m4 > BLAS_spmv_testgen.tmp.c && \
	$(INDENT) $(INDENT_OPTS) BLAS_spmv_testgen.tmp.c && \
	mv BLAS_spmv_testgen.tmp.c $@ && rm -f BLAS_spmv_testgen.tmp.c*

$(SRC_PATH)/spmv-support.c: spmv-support.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=spmv spmv-support.m4 > spmv-support.tmp.c && \
	$(INDENT) $(INDENT_OPTS) spmv-support.tmp.c && \
	mv spmv-support.tmp.c $@ && rm -f spmv-support.tmp.c*


source-clean:
	rm -f $(TEST_SRCS)

