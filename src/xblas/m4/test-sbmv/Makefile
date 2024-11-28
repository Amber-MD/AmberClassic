include ../../make.conf
include ../../$(MAKEINC)

SRC_PATH = ../../testing/test-sbmv
M4_PATH = -I. -I..

TEST_SRCS = \
  $(SRC_PATH)/do_test_sbmv.c\
  $(SRC_PATH)/BLAS_sbmv_testgen.c\
  $(SRC_PATH)/sbmv-support.c

test-source: $(TEST_SRCS)


$(SRC_PATH)/do_test_sbmv.c: do_test_sbmv.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=sbmv do_test_sbmv.m4 > do_test_sbmv.tmp.c && \
	$(INDENT) $(INDENT_OPTS) do_test_sbmv.tmp.c && \
	mv do_test_sbmv.tmp.c $@ && rm -f do_test_sbmv.tmp.c*

$(SRC_PATH)/BLAS_sbmv_testgen.c: BLAS_sbmv_testgen.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=sbmv BLAS_sbmv_testgen.m4 > BLAS_sbmv_testgen.tmp.c && \
	$(INDENT) $(INDENT_OPTS) BLAS_sbmv_testgen.tmp.c && \
	mv BLAS_sbmv_testgen.tmp.c $@ && rm -f BLAS_sbmv_testgen.tmp.c*

$(SRC_PATH)/sbmv-support.c: sbmv-support.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=sbmv sbmv-support.m4 > sbmv-support.tmp.c && \
	$(INDENT) $(INDENT_OPTS) sbmv-support.tmp.c && \
	mv sbmv-support.tmp.c $@ && rm -f sbmv-support.tmp.c*


source-clean:
	rm -f $(TEST_SRCS)

