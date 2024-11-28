include ../../make.conf
include ../../$(MAKEINC)

SRC_PATH = ../../testing/test-ge_sum_mv
M4_PATH = -I. -I..

TEST_SRCS = \
  $(SRC_PATH)/do_test_ge_sum_mv.c\
  $(SRC_PATH)/BLAS_ge_sum_mv_testgen.c

test-source: $(TEST_SRCS)


$(SRC_PATH)/do_test_ge_sum_mv.c: do_test_ge_sum_mv.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=ge_sum_mv do_test_ge_sum_mv.m4 > do_test_ge_sum_mv.tmp.c && \
	$(INDENT) $(INDENT_OPTS) do_test_ge_sum_mv.tmp.c && \
	mv do_test_ge_sum_mv.tmp.c $@ && rm -f do_test_ge_sum_mv.tmp.c*

$(SRC_PATH)/BLAS_ge_sum_mv_testgen.c: BLAS_ge_sum_mv_testgen.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=ge_sum_mv BLAS_ge_sum_mv_testgen.m4 > BLAS_ge_sum_mv_testgen.tmp.c && \
	$(INDENT) $(INDENT_OPTS) BLAS_ge_sum_mv_testgen.tmp.c && \
	mv BLAS_ge_sum_mv_testgen.tmp.c $@ && rm -f BLAS_ge_sum_mv_testgen.tmp.c*


source-clean:
	rm -f $(TEST_SRCS)

