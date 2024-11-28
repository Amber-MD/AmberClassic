include ../../make.conf
include ../../$(MAKEINC)

SRC_PATH = ../../testing/test-dot2
M4_PATH = -I. -I..

TEST_SRCS = \
  $(SRC_PATH)/do_test_dot2.c\
  $(SRC_PATH)/BLAS_dot2_testgen.c\
  $(SRC_PATH)/r_truth2.c\
  $(SRC_PATH)/dot2.c\
  $(SRC_PATH)/test_dot2.c

test-source: $(TEST_SRCS)


$(SRC_PATH)/do_test_dot2.c: do_test_dot2.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=dot2 do_test_dot2.m4 > do_test_dot2.tmp.c && \
	$(INDENT) $(INDENT_OPTS) do_test_dot2.tmp.c && \
	mv do_test_dot2.tmp.c $@ && rm -f do_test_dot2.tmp.c*

$(SRC_PATH)/BLAS_dot2_testgen.c: BLAS_dot2_testgen.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=dot2 BLAS_dot2_testgen.m4 > BLAS_dot2_testgen.tmp.c && \
	$(INDENT) $(INDENT_OPTS) BLAS_dot2_testgen.tmp.c && \
	mv BLAS_dot2_testgen.tmp.c $@ && rm -f BLAS_dot2_testgen.tmp.c*

$(SRC_PATH)/r_truth2.c: r_truth2.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=dot2 r_truth2.m4 > r_truth2.tmp.c && \
	$(INDENT) $(INDENT_OPTS) r_truth2.tmp.c && \
	mv r_truth2.tmp.c $@ && rm -f r_truth2.tmp.c*

$(SRC_PATH)/dot2.c: dot2.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=dot2 dot2.m4 > dot2.tmp.c && \
	$(INDENT) $(INDENT_OPTS) dot2.tmp.c && \
	mv dot2.tmp.c $@ && rm -f dot2.tmp.c*

$(SRC_PATH)/test_dot2.c: test_dot2.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=dot2 test_dot2.m4 > test_dot2.tmp.c && \
	$(INDENT) $(INDENT_OPTS) test_dot2.tmp.c && \
	mv test_dot2.tmp.c $@ && rm -f test_dot2.tmp.c*


source-clean:
	rm -f $(TEST_SRCS)

