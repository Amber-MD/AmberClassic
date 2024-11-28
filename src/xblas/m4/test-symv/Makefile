include ../../make.conf
include ../../$(MAKEINC)

SRC_PATH = ../../testing/test-symv
M4_PATH = -I. -I..

TEST_SRCS = \
  $(SRC_PATH)/do_test_symv.c\
  $(SRC_PATH)/BLAS_symv_testgen.c\
  $(SRC_PATH)/symv-support.c

test-source: $(TEST_SRCS)


$(SRC_PATH)/do_test_symv.c: do_test_symv.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=symv do_test_symv.m4 > do_test_symv.tmp.c && \
	$(INDENT) $(INDENT_OPTS) do_test_symv.tmp.c && \
	mv do_test_symv.tmp.c $@ && rm -f do_test_symv.tmp.c*

$(SRC_PATH)/BLAS_symv_testgen.c: BLAS_symv_testgen.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=symv BLAS_symv_testgen.m4 > BLAS_symv_testgen.tmp.c && \
	$(INDENT) $(INDENT_OPTS) BLAS_symv_testgen.tmp.c && \
	mv BLAS_symv_testgen.tmp.c $@ && rm -f BLAS_symv_testgen.tmp.c*

$(SRC_PATH)/symv-support.c: symv-support.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=symv symv-support.m4 > symv-support.tmp.c && \
	$(INDENT) $(INDENT_OPTS) symv-support.tmp.c && \
	mv symv-support.tmp.c $@ && rm -f symv-support.tmp.c*


source-clean:
	rm -f $(TEST_SRCS)

