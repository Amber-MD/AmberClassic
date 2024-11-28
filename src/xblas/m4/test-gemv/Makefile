include ../../make.conf
include ../../$(MAKEINC)

SRC_PATH = ../../testing/test-gemv
M4_PATH = -I. -I..

TEST_SRCS = \
  $(SRC_PATH)/do_test_gemv.c\
  $(SRC_PATH)/BLAS_gemv_testgen.c\
  $(SRC_PATH)/gemv-support.c

test-source: $(TEST_SRCS)


$(SRC_PATH)/do_test_gemv.c: do_test_gemv.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=gemv do_test_gemv.m4 > do_test_gemv.tmp.c && \
	$(INDENT) $(INDENT_OPTS) do_test_gemv.tmp.c && \
	mv do_test_gemv.tmp.c $@ && rm -f do_test_gemv.tmp.c*

$(SRC_PATH)/BLAS_gemv_testgen.c: BLAS_gemv_testgen.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=gemv BLAS_gemv_testgen.m4 > BLAS_gemv_testgen.tmp.c && \
	$(INDENT) $(INDENT_OPTS) BLAS_gemv_testgen.tmp.c && \
	mv BLAS_gemv_testgen.tmp.c $@ && rm -f BLAS_gemv_testgen.tmp.c*

$(SRC_PATH)/gemv-support.c: gemv-support.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=gemv gemv-support.m4 > gemv-support.tmp.c && \
	$(INDENT) $(INDENT_OPTS) gemv-support.tmp.c && \
	mv gemv-support.tmp.c $@ && rm -f gemv-support.tmp.c*


source-clean:
	rm -f $(TEST_SRCS)

