include ../../make.conf
include ../../$(MAKEINC)

SRC_PATH = ../../testing/test-gemm
M4_PATH = -I. -I..

TEST_SRCS = \
  $(SRC_PATH)/do_test_gemm.c\
  $(SRC_PATH)/BLAS_gemm_testgen.c

test-source: $(TEST_SRCS)


$(SRC_PATH)/do_test_gemm.c: do_test_gemm.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=gemm do_test_gemm.m4 > do_test_gemm.tmp.c && \
	$(INDENT) $(INDENT_OPTS) do_test_gemm.tmp.c && \
	mv do_test_gemm.tmp.c $@ && rm -f do_test_gemm.tmp.c*

$(SRC_PATH)/BLAS_gemm_testgen.c: BLAS_gemm_testgen.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=gemm BLAS_gemm_testgen.m4 > BLAS_gemm_testgen.tmp.c && \
	$(INDENT) $(INDENT_OPTS) BLAS_gemm_testgen.tmp.c && \
	mv BLAS_gemm_testgen.tmp.c $@ && rm -f BLAS_gemm_testgen.tmp.c*


source-clean:
	rm -f $(TEST_SRCS)

