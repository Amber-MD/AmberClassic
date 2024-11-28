include ../../make.conf
include ../../$(MAKEINC)

SRC_PATH = ../../testing/test-dot
M4_PATH = -I. -I..

TEST_SRCS = \
  $(SRC_PATH)/do_test_dot.c\
  $(SRC_PATH)/BLAS_dot_testgen.c\
  $(SRC_PATH)/test_dot.c\
  $(SRC_PATH)/print_vector.c\
  $(SRC_PATH)/copy_vector.c

test-source: $(TEST_SRCS)


$(SRC_PATH)/do_test_dot.c: do_test_dot.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=dot do_test_dot.m4 > do_test_dot.tmp.c && \
	$(INDENT) $(INDENT_OPTS) do_test_dot.tmp.c && \
	mv do_test_dot.tmp.c $@ && rm -f do_test_dot.tmp.c*

$(SRC_PATH)/BLAS_dot_testgen.c: BLAS_dot_testgen.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=dot BLAS_dot_testgen.m4 > BLAS_dot_testgen.tmp.c && \
	$(INDENT) $(INDENT_OPTS) BLAS_dot_testgen.tmp.c && \
	mv BLAS_dot_testgen.tmp.c $@ && rm -f BLAS_dot_testgen.tmp.c*

$(SRC_PATH)/test_dot.c: test_dot.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=dot test_dot.m4 > test_dot.tmp.c && \
	$(INDENT) $(INDENT_OPTS) test_dot.tmp.c && \
	mv test_dot.tmp.c $@ && rm -f test_dot.tmp.c*

$(SRC_PATH)/print_vector.c: print_vector.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=dot print_vector.m4 > print_vector.tmp.c && \
	$(INDENT) $(INDENT_OPTS) print_vector.tmp.c && \
	mv print_vector.tmp.c $@ && rm -f print_vector.tmp.c*

$(SRC_PATH)/copy_vector.c: copy_vector.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=dot copy_vector.m4 > copy_vector.tmp.c && \
	$(INDENT) $(INDENT_OPTS) copy_vector.tmp.c && \
	mv copy_vector.tmp.c $@ && rm -f copy_vector.tmp.c*


source-clean:
	rm -f $(TEST_SRCS)

