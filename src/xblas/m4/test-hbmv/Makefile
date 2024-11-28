include ../../make.conf
include ../../$(MAKEINC)

SRC_PATH = ../../testing/test-hbmv
M4_PATH = -I. -I..

TEST_SRCS = \
  $(SRC_PATH)/do_test_hbmv.c\
  $(SRC_PATH)/BLAS_hbmv_testgen.c\
  $(SRC_PATH)/hbmv-support.c

test-source: $(TEST_SRCS)


$(SRC_PATH)/do_test_hbmv.c: do_test_hbmv.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=hbmv do_test_hbmv.m4 > do_test_hbmv.tmp.c && \
	$(INDENT) $(INDENT_OPTS) do_test_hbmv.tmp.c && \
	mv do_test_hbmv.tmp.c $@ && rm -f do_test_hbmv.tmp.c*

$(SRC_PATH)/BLAS_hbmv_testgen.c: BLAS_hbmv_testgen.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=hbmv BLAS_hbmv_testgen.m4 > BLAS_hbmv_testgen.tmp.c && \
	$(INDENT) $(INDENT_OPTS) BLAS_hbmv_testgen.tmp.c && \
	mv BLAS_hbmv_testgen.tmp.c $@ && rm -f BLAS_hbmv_testgen.tmp.c*

$(SRC_PATH)/hbmv-support.c: hbmv-support.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=hbmv hbmv-support.m4 > hbmv-support.tmp.c && \
	$(INDENT) $(INDENT_OPTS) hbmv-support.tmp.c && \
	mv hbmv-support.tmp.c $@ && rm -f hbmv-support.tmp.c*


source-clean:
	rm -f $(TEST_SRCS)

