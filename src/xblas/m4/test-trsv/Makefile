include ../../make.conf
include ../../$(MAKEINC)

SRC_PATH = ../../testing/test-trsv
M4_PATH = -I. -I..

TEST_SRCS = \
  $(SRC_PATH)/do_test_trsv.c\
  $(SRC_PATH)/BLAS_trsv_testgen.c\
  $(SRC_PATH)/trsv-support.c

test-source: $(TEST_SRCS)


$(SRC_PATH)/do_test_trsv.c: do_test_trsv.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=trsv do_test_trsv.m4 > do_test_trsv.tmp.c && \
	$(INDENT) $(INDENT_OPTS) do_test_trsv.tmp.c && \
	mv do_test_trsv.tmp.c $@ && rm -f do_test_trsv.tmp.c*

$(SRC_PATH)/BLAS_trsv_testgen.c: BLAS_trsv_testgen.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=trsv BLAS_trsv_testgen.m4 > BLAS_trsv_testgen.tmp.c && \
	$(INDENT) $(INDENT_OPTS) BLAS_trsv_testgen.tmp.c && \
	mv BLAS_trsv_testgen.tmp.c $@ && rm -f BLAS_trsv_testgen.tmp.c*

$(SRC_PATH)/trsv-support.c: trsv-support.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=trsv trsv-support.m4 > trsv-support.tmp.c && \
	$(INDENT) $(INDENT_OPTS) trsv-support.tmp.c && \
	mv trsv-support.tmp.c $@ && rm -f trsv-support.tmp.c*


source-clean:
	rm -f $(TEST_SRCS)

