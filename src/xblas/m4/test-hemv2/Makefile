include ../../make.conf
include ../../$(MAKEINC)

SRC_PATH = ../../testing/test-hemv2
M4_PATH = -I. -I..

TEST_SRCS = \
  $(SRC_PATH)/do_test_hemv2.c\
  $(SRC_PATH)/BLAS_hemv2_testgen.c

test-source: $(TEST_SRCS)


$(SRC_PATH)/do_test_hemv2.c: do_test_hemv2.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=hemv2 do_test_hemv2.m4 > do_test_hemv2.tmp.c && \
	$(INDENT) $(INDENT_OPTS) do_test_hemv2.tmp.c && \
	mv do_test_hemv2.tmp.c $@ && rm -f do_test_hemv2.tmp.c*

$(SRC_PATH)/BLAS_hemv2_testgen.c: BLAS_hemv2_testgen.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=hemv2 BLAS_hemv2_testgen.m4 > BLAS_hemv2_testgen.tmp.c && \
	$(INDENT) $(INDENT_OPTS) BLAS_hemv2_testgen.tmp.c && \
	mv BLAS_hemv2_testgen.tmp.c $@ && rm -f BLAS_hemv2_testgen.tmp.c*


source-clean:
	rm -f $(TEST_SRCS)

