include ../../make.conf
include ../../$(MAKEINC)


SRC_PATH = ../../src/axpby
M4_PATH = -I. -I..

AXPBY_SRCS = \
	$(SRC_PATH)/BLAS_saxpby.c\
	$(SRC_PATH)/BLAS_daxpby.c\
	$(SRC_PATH)/BLAS_caxpby.c\
	$(SRC_PATH)/BLAS_zaxpby.c\
	$(SRC_PATH)/BLAS_daxpby_s.c\
	$(SRC_PATH)/BLAS_caxpby_s.c\
	$(SRC_PATH)/BLAS_zaxpby_c.c\
	$(SRC_PATH)/BLAS_zaxpby_d.c\
	$(SRC_PATH)/BLAS_saxpby_x.c\
	$(SRC_PATH)/BLAS_daxpby_x.c\
	$(SRC_PATH)/BLAS_caxpby_x.c\
	$(SRC_PATH)/BLAS_zaxpby_x.c\
	$(SRC_PATH)/BLAS_daxpby_s_x.c\
	$(SRC_PATH)/BLAS_zaxpby_c_x.c\
	$(SRC_PATH)/BLAS_caxpby_s_x.c\
	$(SRC_PATH)/BLAS_zaxpby_d_x.c \
	\
	$(SRC_PATH)/BLAS_saxpby-f2c.c\
	$(SRC_PATH)/BLAS_daxpby-f2c.c\
	$(SRC_PATH)/BLAS_caxpby-f2c.c\
	$(SRC_PATH)/BLAS_zaxpby-f2c.c\
	$(SRC_PATH)/BLAS_daxpby_s-f2c.c\
	$(SRC_PATH)/BLAS_caxpby_s-f2c.c\
	$(SRC_PATH)/BLAS_zaxpby_c-f2c.c\
	$(SRC_PATH)/BLAS_zaxpby_d-f2c.c\
	$(SRC_PATH)/BLAS_saxpby_x-f2c.c\
	$(SRC_PATH)/BLAS_daxpby_x-f2c.c\
	$(SRC_PATH)/BLAS_caxpby_x-f2c.c\
	$(SRC_PATH)/BLAS_zaxpby_x-f2c.c\
	$(SRC_PATH)/BLAS_daxpby_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zaxpby_c_x-f2c.c\
	$(SRC_PATH)/BLAS_caxpby_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zaxpby_d_x-f2c.c

axpby-source: $(AXPBY_SRCS)

axpby-header: axpby.h

$(SRC_PATH)/BLAS_saxpby.c: axpby.m4 axpby-common.m4 ../cblas.m4
	echo "AXPBY(s, s)" | cat axpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_saxpby.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_saxpby.c && \
    mv tmp.BLAS_saxpby.c $@ && rm -f tmp.BLAS_saxpby.c*

$(SRC_PATH)/BLAS_daxpby.c: axpby.m4 axpby-common.m4 ../cblas.m4
	echo "AXPBY(d, d)" | cat axpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_daxpby.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_daxpby.c && \
    mv tmp.BLAS_daxpby.c $@ && rm -f tmp.BLAS_daxpby.c*

$(SRC_PATH)/BLAS_caxpby.c: axpby.m4 axpby-common.m4 ../cblas.m4
	echo "AXPBY(c, c)" | cat axpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_caxpby.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_caxpby.c && \
    mv tmp.BLAS_caxpby.c $@ && rm -f tmp.BLAS_caxpby.c*

$(SRC_PATH)/BLAS_zaxpby.c: axpby.m4 axpby-common.m4 ../cblas.m4
	echo "AXPBY(z, z)" | cat axpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zaxpby.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zaxpby.c && \
    mv tmp.BLAS_zaxpby.c $@ && rm -f tmp.BLAS_zaxpby.c*

$(SRC_PATH)/BLAS_daxpby_s.c: axpby.m4 axpby-common.m4 ../cblas.m4
	echo "AXPBY(d, s)" | cat axpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_daxpby_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_daxpby_s.c && \
    mv tmp.BLAS_daxpby_s.c $@ && rm -f tmp.BLAS_daxpby_s.c*

$(SRC_PATH)/BLAS_caxpby_s.c: axpby.m4 axpby-common.m4 ../cblas.m4
	echo "AXPBY(c, s)" | cat axpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_caxpby_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_caxpby_s.c && \
    mv tmp.BLAS_caxpby_s.c $@ && rm -f tmp.BLAS_caxpby_s.c*

$(SRC_PATH)/BLAS_zaxpby_c.c: axpby.m4 axpby-common.m4 ../cblas.m4
	echo "AXPBY(z, c)" | cat axpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zaxpby_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zaxpby_c.c && \
    mv tmp.BLAS_zaxpby_c.c $@ && rm -f tmp.BLAS_zaxpby_c.c*

$(SRC_PATH)/BLAS_zaxpby_d.c: axpby.m4 axpby-common.m4 ../cblas.m4
	echo "AXPBY(z, d)" | cat axpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zaxpby_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zaxpby_d.c && \
    mv tmp.BLAS_zaxpby_d.c $@ && rm -f tmp.BLAS_zaxpby_d.c*

$(SRC_PATH)/BLAS_saxpby_x.c: axpby.m4 axpby-common.m4 ../cblas.m4
	echo "AXPBY(s, s, _x)" | cat axpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_saxpby_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_saxpby_x.c && \
    mv tmp.BLAS_saxpby_x.c $@ && rm -f tmp.BLAS_saxpby_x.c*

$(SRC_PATH)/BLAS_daxpby_x.c: axpby.m4 axpby-common.m4 ../cblas.m4
	echo "AXPBY(d, d, _x)" | cat axpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_daxpby_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_daxpby_x.c && \
    mv tmp.BLAS_daxpby_x.c $@ && rm -f tmp.BLAS_daxpby_x.c*

$(SRC_PATH)/BLAS_caxpby_x.c: axpby.m4 axpby-common.m4 ../cblas.m4
	echo "AXPBY(c, c, _x)" | cat axpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_caxpby_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_caxpby_x.c && \
    mv tmp.BLAS_caxpby_x.c $@ && rm -f tmp.BLAS_caxpby_x.c*

$(SRC_PATH)/BLAS_zaxpby_x.c: axpby.m4 axpby-common.m4 ../cblas.m4
	echo "AXPBY(z, z, _x)" | cat axpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zaxpby_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zaxpby_x.c && \
    mv tmp.BLAS_zaxpby_x.c $@ && rm -f tmp.BLAS_zaxpby_x.c*

$(SRC_PATH)/BLAS_daxpby_s_x.c: axpby.m4 axpby-common.m4 ../cblas.m4
	echo "AXPBY(d, s, _x)" | cat axpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_daxpby_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_daxpby_s_x.c && \
    mv tmp.BLAS_daxpby_s_x.c $@ && rm -f tmp.BLAS_daxpby_s_x.c*

$(SRC_PATH)/BLAS_zaxpby_c_x.c: axpby.m4 axpby-common.m4 ../cblas.m4
	echo "AXPBY(z, c, _x)" | cat axpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zaxpby_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zaxpby_c_x.c && \
    mv tmp.BLAS_zaxpby_c_x.c $@ && rm -f tmp.BLAS_zaxpby_c_x.c*

$(SRC_PATH)/BLAS_caxpby_s_x.c: axpby.m4 axpby-common.m4 ../cblas.m4
	echo "AXPBY(c, s, _x)" | cat axpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_caxpby_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_caxpby_s_x.c && \
    mv tmp.BLAS_caxpby_s_x.c $@ && rm -f tmp.BLAS_caxpby_s_x.c*

$(SRC_PATH)/BLAS_zaxpby_d_x.c: axpby.m4 axpby-common.m4 ../cblas.m4
	echo "AXPBY(z, d, _x)" | cat axpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zaxpby_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zaxpby_d_x.c && \
    mv tmp.BLAS_zaxpby_d_x.c $@ && rm -f tmp.BLAS_zaxpby_d_x.c*



$(SRC_PATH)/BLAS_saxpby-f2c.c: ../f2c-bridge.m4 axpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=axpby -DRARG="s, s" ../f2c-bridge.m4 > tmp.BLAS_saxpby-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_saxpby-f2c.c && \
    mv tmp.BLAS_saxpby-f2c.c $@ && rm -f tmp.BLAS_saxpby-f2c.c*

$(SRC_PATH)/BLAS_daxpby-f2c.c: ../f2c-bridge.m4 axpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=axpby -DRARG="d, d" ../f2c-bridge.m4 > tmp.BLAS_daxpby-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_daxpby-f2c.c && \
    mv tmp.BLAS_daxpby-f2c.c $@ && rm -f tmp.BLAS_daxpby-f2c.c*

$(SRC_PATH)/BLAS_caxpby-f2c.c: ../f2c-bridge.m4 axpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=axpby -DRARG="c, c" ../f2c-bridge.m4 > tmp.BLAS_caxpby-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_caxpby-f2c.c && \
    mv tmp.BLAS_caxpby-f2c.c $@ && rm -f tmp.BLAS_caxpby-f2c.c*

$(SRC_PATH)/BLAS_zaxpby-f2c.c: ../f2c-bridge.m4 axpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=axpby -DRARG="z, z" ../f2c-bridge.m4 > tmp.BLAS_zaxpby-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zaxpby-f2c.c && \
    mv tmp.BLAS_zaxpby-f2c.c $@ && rm -f tmp.BLAS_zaxpby-f2c.c*

$(SRC_PATH)/BLAS_daxpby_s-f2c.c: ../f2c-bridge.m4 axpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=axpby -DRARG="d, s" ../f2c-bridge.m4 > tmp.BLAS_daxpby_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_daxpby_s-f2c.c && \
    mv tmp.BLAS_daxpby_s-f2c.c $@ && rm -f tmp.BLAS_daxpby_s-f2c.c*

$(SRC_PATH)/BLAS_caxpby_s-f2c.c: ../f2c-bridge.m4 axpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=axpby -DRARG="c, s" ../f2c-bridge.m4 > tmp.BLAS_caxpby_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_caxpby_s-f2c.c && \
    mv tmp.BLAS_caxpby_s-f2c.c $@ && rm -f tmp.BLAS_caxpby_s-f2c.c*

$(SRC_PATH)/BLAS_zaxpby_c-f2c.c: ../f2c-bridge.m4 axpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=axpby -DRARG="z, c" ../f2c-bridge.m4 > tmp.BLAS_zaxpby_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zaxpby_c-f2c.c && \
    mv tmp.BLAS_zaxpby_c-f2c.c $@ && rm -f tmp.BLAS_zaxpby_c-f2c.c*

$(SRC_PATH)/BLAS_zaxpby_d-f2c.c: ../f2c-bridge.m4 axpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=axpby -DRARG="z, d" ../f2c-bridge.m4 > tmp.BLAS_zaxpby_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zaxpby_d-f2c.c && \
    mv tmp.BLAS_zaxpby_d-f2c.c $@ && rm -f tmp.BLAS_zaxpby_d-f2c.c*

$(SRC_PATH)/BLAS_saxpby_x-f2c.c: ../f2c-bridge.m4 axpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=axpby -DRARG="s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_saxpby_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_saxpby_x-f2c.c && \
    mv tmp.BLAS_saxpby_x-f2c.c $@ && rm -f tmp.BLAS_saxpby_x-f2c.c*

$(SRC_PATH)/BLAS_daxpby_x-f2c.c: ../f2c-bridge.m4 axpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=axpby -DRARG="d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_daxpby_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_daxpby_x-f2c.c && \
    mv tmp.BLAS_daxpby_x-f2c.c $@ && rm -f tmp.BLAS_daxpby_x-f2c.c*

$(SRC_PATH)/BLAS_caxpby_x-f2c.c: ../f2c-bridge.m4 axpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=axpby -DRARG="c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_caxpby_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_caxpby_x-f2c.c && \
    mv tmp.BLAS_caxpby_x-f2c.c $@ && rm -f tmp.BLAS_caxpby_x-f2c.c*

$(SRC_PATH)/BLAS_zaxpby_x-f2c.c: ../f2c-bridge.m4 axpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=axpby -DRARG="z, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zaxpby_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zaxpby_x-f2c.c && \
    mv tmp.BLAS_zaxpby_x-f2c.c $@ && rm -f tmp.BLAS_zaxpby_x-f2c.c*

$(SRC_PATH)/BLAS_daxpby_s_x-f2c.c: ../f2c-bridge.m4 axpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=axpby -DRARG="d, s, _x" ../f2c-bridge.m4 > tmp.BLAS_daxpby_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_daxpby_s_x-f2c.c && \
    mv tmp.BLAS_daxpby_s_x-f2c.c $@ && rm -f tmp.BLAS_daxpby_s_x-f2c.c*

$(SRC_PATH)/BLAS_zaxpby_c_x-f2c.c: ../f2c-bridge.m4 axpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=axpby -DRARG="z, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zaxpby_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zaxpby_c_x-f2c.c && \
    mv tmp.BLAS_zaxpby_c_x-f2c.c $@ && rm -f tmp.BLAS_zaxpby_c_x-f2c.c*

$(SRC_PATH)/BLAS_caxpby_s_x-f2c.c: ../f2c-bridge.m4 axpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=axpby -DRARG="c, s, _x" ../f2c-bridge.m4 > tmp.BLAS_caxpby_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_caxpby_s_x-f2c.c && \
    mv tmp.BLAS_caxpby_s_x-f2c.c $@ && rm -f tmp.BLAS_caxpby_s_x-f2c.c*

$(SRC_PATH)/BLAS_zaxpby_d_x-f2c.c: ../f2c-bridge.m4 axpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=axpby -DRARG="z, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zaxpby_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zaxpby_d_x-f2c.c && \
    mv tmp.BLAS_zaxpby_d_x-f2c.c $@ && rm -f tmp.BLAS_zaxpby_d_x-f2c.c*


axpby.h: ../cblas.m4 ../header.template axpby-common.m4
	sed -e "s/routine/axpby/g" -e "s/ROUTINE/AXPBY/g" \
  ../header.template | $(M4) $(M4_OPTS) $(M4_PATH) - >$@

source-clean:
	rm -f $(AXPBY_SRCS) axpby.h

.PHONY: axpby-source axpby-header source-clean

