include ../../make.conf
include ../../$(MAKEINC)

SRC_PATH = ../../testing/test-symv2
M4_PATH = -I. -I..

TEST_SRCS = \
  $(SRC_PATH)/do_test_symv2.c\
  $(SRC_PATH)/BLAS_symv2_testgen.c

test-source: $(TEST_SRCS)


$(SRC_PATH)/do_test_symv2.c: do_test_symv2.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=symv2 do_test_symv2.m4 > do_test_symv2.tmp.c && \
	$(INDENT) $(INDENT_OPTS) do_test_symv2.tmp.c && \
	mv do_test_symv2.tmp.c $@ && rm -f do_test_symv2.tmp.c*

$(SRC_PATH)/BLAS_symv2_testgen.c: BLAS_symv2_testgen.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=symv2 BLAS_symv2_testgen.m4 > BLAS_symv2_testgen.tmp.c && \
	$(INDENT) $(INDENT_OPTS) BLAS_symv2_testgen.tmp.c && \
	mv BLAS_symv2_testgen.tmp.c $@ && rm -f BLAS_symv2_testgen.tmp.c*


source-clean:
	rm -f $(TEST_SRCS)

