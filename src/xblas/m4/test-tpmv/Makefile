include ../../make.conf
include ../../$(MAKEINC)

SRC_PATH = ../../testing/test-tpmv
M4_PATH = -I. -I..

TEST_SRCS = \
  $(SRC_PATH)/do_test_tpmv.c\
  $(SRC_PATH)/BLAS_tpmv_testgen.c\
  $(SRC_PATH)/tpmv-support.c

test-source: $(TEST_SRCS)


$(SRC_PATH)/do_test_tpmv.c: do_test_tpmv.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=tpmv do_test_tpmv.m4 > do_test_tpmv.tmp.c && \
	$(INDENT) $(INDENT_OPTS) do_test_tpmv.tmp.c && \
	mv do_test_tpmv.tmp.c $@ && rm -f do_test_tpmv.tmp.c*

$(SRC_PATH)/BLAS_tpmv_testgen.c: BLAS_tpmv_testgen.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=tpmv BLAS_tpmv_testgen.m4 > BLAS_tpmv_testgen.tmp.c && \
	$(INDENT) $(INDENT_OPTS) BLAS_tpmv_testgen.tmp.c && \
	mv BLAS_tpmv_testgen.tmp.c $@ && rm -f BLAS_tpmv_testgen.tmp.c*

$(SRC_PATH)/tpmv-support.c: tpmv-support.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=tpmv tpmv-support.m4 > tpmv-support.tmp.c && \
	$(INDENT) $(INDENT_OPTS) tpmv-support.tmp.c && \
	mv tpmv-support.tmp.c $@ && rm -f tpmv-support.tmp.c*


source-clean:
	rm -f $(TEST_SRCS)

