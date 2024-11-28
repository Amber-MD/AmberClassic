include ../../make.conf
include ../../$(MAKEINC)

SRC_PATH = ../../testing/test-gbmv
M4_PATH = -I. -I..

TEST_SRCS = \
  $(SRC_PATH)/do_test_gbmv.c\
  $(SRC_PATH)/BLAS_gbmv_testgen.c\
  $(SRC_PATH)/gbmv-support.c

test-source: $(TEST_SRCS)


$(SRC_PATH)/do_test_gbmv.c: do_test_gbmv.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=gbmv do_test_gbmv.m4 > do_test_gbmv.tmp.c && \
	$(INDENT) $(INDENT_OPTS) do_test_gbmv.tmp.c && \
	mv do_test_gbmv.tmp.c $@ && rm -f do_test_gbmv.tmp.c*

$(SRC_PATH)/BLAS_gbmv_testgen.c: BLAS_gbmv_testgen.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=gbmv BLAS_gbmv_testgen.m4 > BLAS_gbmv_testgen.tmp.c && \
	$(INDENT) $(INDENT_OPTS) BLAS_gbmv_testgen.tmp.c && \
	mv BLAS_gbmv_testgen.tmp.c $@ && rm -f BLAS_gbmv_testgen.tmp.c*

$(SRC_PATH)/gbmv-support.c: gbmv-support.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=gbmv gbmv-support.m4 > gbmv-support.tmp.c && \
	$(INDENT) $(INDENT_OPTS) gbmv-support.tmp.c && \
	mv gbmv-support.tmp.c $@ && rm -f gbmv-support.tmp.c*


source-clean:
	rm -f $(TEST_SRCS)

