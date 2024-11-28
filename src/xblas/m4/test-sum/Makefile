include ../../make.conf
include ../../$(MAKEINC)

SRC_PATH = ../../testing/test-sum
M4_PATH = -I. -I..

TEST_SRCS = \
  $(SRC_PATH)/do_test_sum.c\
  $(SRC_PATH)/BLAS_sum_testgen.c\
  $(SRC_PATH)/sum-support.c

test-source: $(TEST_SRCS)


$(SRC_PATH)/do_test_sum.c: do_test_sum.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=sum do_test_sum.m4 > do_test_sum.tmp.c && \
	$(INDENT) $(INDENT_OPTS) do_test_sum.tmp.c && \
	mv do_test_sum.tmp.c $@ && rm -f do_test_sum.tmp.c*

$(SRC_PATH)/BLAS_sum_testgen.c: BLAS_sum_testgen.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=sum BLAS_sum_testgen.m4 > BLAS_sum_testgen.tmp.c && \
	$(INDENT) $(INDENT_OPTS) BLAS_sum_testgen.tmp.c && \
	mv BLAS_sum_testgen.tmp.c $@ && rm -f BLAS_sum_testgen.tmp.c*

$(SRC_PATH)/sum-support.c: sum-support.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=sum sum-support.m4 > sum-support.tmp.c && \
	$(INDENT) $(INDENT_OPTS) sum-support.tmp.c && \
	mv sum-support.tmp.c $@ && rm -f sum-support.tmp.c*


source-clean:
	rm -f $(TEST_SRCS)

