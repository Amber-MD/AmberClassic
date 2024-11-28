include ../../make.conf
include ../../$(MAKEINC)

SRC_PATH = ../../testing/test-waxpby
M4_PATH = -I. -I..

TEST_SRCS = \
  $(SRC_PATH)/do_test_waxpby.c

test-source: $(TEST_SRCS)


$(SRC_PATH)/do_test_waxpby.c: do_test_waxpby.m4 ../cblas.m4 ../test-common.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -D BASE_ROUTINE=waxpby do_test_waxpby.m4 > do_test_waxpby.tmp.c && \
	$(INDENT) $(INDENT_OPTS) do_test_waxpby.tmp.c && \
	mv do_test_waxpby.tmp.c $@ && rm -f do_test_waxpby.tmp.c*


source-clean:
	rm -f $(TEST_SRCS)

