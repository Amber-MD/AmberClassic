include ../../make.conf
include ../../$(MAKEINC)

SRC_PATH = ../../testing/test-symm
M4_PATH = -I. -I..

TEST_SRCS = \
  $(SRC_PATH)/do_test_symm.c\
  $(SRC_PATH)/BLAS_symm_testgen.c

test-source: $(TEST_SRCS)


$(SRC_PATH)/do_test_symm.c: do_test_symm.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=symm do_test_symm.m4 > do_test_symm.tmp.c && \
	$(INDENT) $(INDENT_OPTS) do_test_symm.tmp.c && \
	mv do_test_symm.tmp.c $@ && rm -f do_test_symm.tmp.c*

$(SRC_PATH)/BLAS_symm_testgen.c: BLAS_symm_testgen.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=symm BLAS_symm_testgen.m4 > BLAS_symm_testgen.tmp.c && \
	$(INDENT) $(INDENT_OPTS) BLAS_symm_testgen.tmp.c && \
	mv BLAS_symm_testgen.tmp.c $@ && rm -f BLAS_symm_testgen.tmp.c*


source-clean:
	rm -f $(TEST_SRCS)

