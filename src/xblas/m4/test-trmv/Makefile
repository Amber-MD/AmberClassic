include ../../make.conf
include ../../$(MAKEINC)

SRC_PATH = ../../testing/test-trmv
M4_PATH = -I. -I..

TEST_SRCS = \
  $(SRC_PATH)/do_test_trmv.c\
  $(SRC_PATH)/BLAS_trmv_testgen.c\
  $(SRC_PATH)/trmv-support.c

test-source: $(TEST_SRCS)


$(SRC_PATH)/do_test_trmv.c: do_test_trmv.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=trmv do_test_trmv.m4 > do_test_trmv.tmp.c && \
	$(INDENT) $(INDENT_OPTS) do_test_trmv.tmp.c && \
	mv do_test_trmv.tmp.c $@ && rm -f do_test_trmv.tmp.c*

$(SRC_PATH)/BLAS_trmv_testgen.c: BLAS_trmv_testgen.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=trmv BLAS_trmv_testgen.m4 > BLAS_trmv_testgen.tmp.c && \
	$(INDENT) $(INDENT_OPTS) BLAS_trmv_testgen.tmp.c && \
	mv BLAS_trmv_testgen.tmp.c $@ && rm -f BLAS_trmv_testgen.tmp.c*

$(SRC_PATH)/trmv-support.c: trmv-support.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=trmv trmv-support.m4 > trmv-support.tmp.c && \
	$(INDENT) $(INDENT_OPTS) trmv-support.tmp.c && \
	mv trmv-support.tmp.c $@ && rm -f trmv-support.tmp.c*


source-clean:
	rm -f $(TEST_SRCS)

