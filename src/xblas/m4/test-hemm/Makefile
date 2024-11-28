include ../../make.conf
include ../../$(MAKEINC)

SRC_PATH = ../../testing/test-hemm
M4_PATH = -I. -I..

TEST_SRCS = \
  $(SRC_PATH)/do_test_hemm.c\
  $(SRC_PATH)/BLAS_hemm_testgen.c

test-source: $(TEST_SRCS)


$(SRC_PATH)/do_test_hemm.c: do_test_hemm.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=hemm do_test_hemm.m4 > do_test_hemm.tmp.c && \
	$(INDENT) $(INDENT_OPTS) do_test_hemm.tmp.c && \
	mv do_test_hemm.tmp.c $@ && rm -f do_test_hemm.tmp.c*

$(SRC_PATH)/BLAS_hemm_testgen.c: BLAS_hemm_testgen.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=hemm BLAS_hemm_testgen.m4 > BLAS_hemm_testgen.tmp.c && \
	$(INDENT) $(INDENT_OPTS) BLAS_hemm_testgen.tmp.c && \
	mv BLAS_hemm_testgen.tmp.c $@ && rm -f BLAS_hemm_testgen.tmp.c*


source-clean:
	rm -f $(TEST_SRCS)

