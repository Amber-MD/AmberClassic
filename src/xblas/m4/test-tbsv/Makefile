include ../../make.conf
include ../../$(MAKEINC)

SRC_PATH = ../../testing/test-tbsv
M4_PATH = -I. -I..

TEST_SRCS = \
  $(SRC_PATH)/do_test_tbsv.c\
  $(SRC_PATH)/BLAS_tbsv_testgen.c\
  $(SRC_PATH)/tbsv-support.c

test-source: $(TEST_SRCS)


$(SRC_PATH)/do_test_tbsv.c: do_test_tbsv.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=tbsv do_test_tbsv.m4 > do_test_tbsv.tmp.c && \
	$(INDENT) $(INDENT_OPTS) do_test_tbsv.tmp.c && \
	mv do_test_tbsv.tmp.c $@ && rm -f do_test_tbsv.tmp.c*

$(SRC_PATH)/BLAS_tbsv_testgen.c: BLAS_tbsv_testgen.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=tbsv BLAS_tbsv_testgen.m4 > BLAS_tbsv_testgen.tmp.c && \
	$(INDENT) $(INDENT_OPTS) BLAS_tbsv_testgen.tmp.c && \
	mv BLAS_tbsv_testgen.tmp.c $@ && rm -f BLAS_tbsv_testgen.tmp.c*

$(SRC_PATH)/tbsv-support.c: tbsv-support.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=tbsv tbsv-support.m4 > tbsv-support.tmp.c && \
	$(INDENT) $(INDENT_OPTS) tbsv-support.tmp.c && \
	mv tbsv-support.tmp.c $@ && rm -f tbsv-support.tmp.c*


source-clean:
	rm -f $(TEST_SRCS)

