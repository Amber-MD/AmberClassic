include ../../make.conf
include ../../$(MAKEINC)

SRC_PATH = ../../testing/test-axpby
M4_PATH = -I. -I..

TEST_SRCS = \
  $(SRC_PATH)/do_test_axpby.c

test-source: $(TEST_SRCS)


$(SRC_PATH)/do_test_axpby.c: do_test_axpby.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=axpby do_test_axpby.m4 > do_test_axpby.tmp.c && \
	$(INDENT) $(INDENT_OPTS) do_test_axpby.tmp.c && \
	mv do_test_axpby.tmp.c $@ && rm -f do_test_axpby.tmp.c*


source-clean:
	rm -f $(TEST_SRCS)

