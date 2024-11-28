include ../../make.conf
include ../../$(MAKEINC)

SRC_PATH = ../../testing/test-hemv
M4_PATH = -I. -I..

TEST_SRCS = \
  $(SRC_PATH)/do_test_hemv.c\
  $(SRC_PATH)/BLAS_hemv_testgen.c\
  $(SRC_PATH)/hemv-support.c

test-source: $(TEST_SRCS)


$(SRC_PATH)/do_test_hemv.c: do_test_hemv.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=hemv do_test_hemv.m4 > do_test_hemv.tmp.c && \
	$(INDENT) $(INDENT_OPTS) do_test_hemv.tmp.c && \
	mv do_test_hemv.tmp.c $@ && rm -f do_test_hemv.tmp.c*

$(SRC_PATH)/BLAS_hemv_testgen.c: BLAS_hemv_testgen.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=hemv BLAS_hemv_testgen.m4 > BLAS_hemv_testgen.tmp.c && \
	$(INDENT) $(INDENT_OPTS) BLAS_hemv_testgen.tmp.c && \
	mv BLAS_hemv_testgen.tmp.c $@ && rm -f BLAS_hemv_testgen.tmp.c*

$(SRC_PATH)/hemv-support.c: hemv-support.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=hemv hemv-support.m4 > hemv-support.tmp.c && \
	$(INDENT) $(INDENT_OPTS) hemv-support.tmp.c && \
	mv hemv-support.tmp.c $@ && rm -f hemv-support.tmp.c*


source-clean:
	rm -f $(TEST_SRCS)

