include ../../make.conf
include ../../$(MAKEINC)


SRC_PATH = ../../src/trmv
M4_PATH = -I. -I..

TRMV_SRCS = \
	$(SRC_PATH)/BLAS_strmv.c\
	$(SRC_PATH)/BLAS_dtrmv.c\
	$(SRC_PATH)/BLAS_ctrmv.c\
	$(SRC_PATH)/BLAS_ztrmv.c\
	$(SRC_PATH)/BLAS_dtrmv_s.c\
	$(SRC_PATH)/BLAS_ztrmv_c.c\
	$(SRC_PATH)/BLAS_ctrmv_s.c\
	$(SRC_PATH)/BLAS_ztrmv_d.c\
	$(SRC_PATH)/BLAS_strmv_x.c\
	$(SRC_PATH)/BLAS_dtrmv_x.c\
	$(SRC_PATH)/BLAS_ctrmv_x.c\
	$(SRC_PATH)/BLAS_ztrmv_x.c\
	$(SRC_PATH)/BLAS_dtrmv_s_x.c\
	$(SRC_PATH)/BLAS_ztrmv_c_x.c\
	$(SRC_PATH)/BLAS_ctrmv_s_x.c\
	$(SRC_PATH)/BLAS_ztrmv_d_x.c \
	\
	$(SRC_PATH)/BLAS_strmv-f2c.c\
	$(SRC_PATH)/BLAS_dtrmv-f2c.c\
	$(SRC_PATH)/BLAS_ctrmv-f2c.c\
	$(SRC_PATH)/BLAS_ztrmv-f2c.c\
	$(SRC_PATH)/BLAS_dtrmv_s-f2c.c\
	$(SRC_PATH)/BLAS_ztrmv_c-f2c.c\
	$(SRC_PATH)/BLAS_ctrmv_s-f2c.c\
	$(SRC_PATH)/BLAS_ztrmv_d-f2c.c\
	$(SRC_PATH)/BLAS_strmv_x-f2c.c\
	$(SRC_PATH)/BLAS_dtrmv_x-f2c.c\
	$(SRC_PATH)/BLAS_ctrmv_x-f2c.c\
	$(SRC_PATH)/BLAS_ztrmv_x-f2c.c\
	$(SRC_PATH)/BLAS_dtrmv_s_x-f2c.c\
	$(SRC_PATH)/BLAS_ztrmv_c_x-f2c.c\
	$(SRC_PATH)/BLAS_ctrmv_s_x-f2c.c\
	$(SRC_PATH)/BLAS_ztrmv_d_x-f2c.c

trmv-source: $(TRMV_SRCS)

trmv-header: trmv.h

$(SRC_PATH)/BLAS_strmv.c: trmv.m4 trmv-common.m4 ../cblas.m4
	echo "TRMV(s, s)" | cat trmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_strmv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_strmv.c && \
    mv tmp.BLAS_strmv.c $@ && rm -f tmp.BLAS_strmv.c*

$(SRC_PATH)/BLAS_dtrmv.c: trmv.m4 trmv-common.m4 ../cblas.m4
	echo "TRMV(d, d)" | cat trmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dtrmv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtrmv.c && \
    mv tmp.BLAS_dtrmv.c $@ && rm -f tmp.BLAS_dtrmv.c*

$(SRC_PATH)/BLAS_ctrmv.c: trmv.m4 trmv-common.m4 ../cblas.m4
	echo "TRMV(c, c)" | cat trmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ctrmv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctrmv.c && \
    mv tmp.BLAS_ctrmv.c $@ && rm -f tmp.BLAS_ctrmv.c*

$(SRC_PATH)/BLAS_ztrmv.c: trmv.m4 trmv-common.m4 ../cblas.m4
	echo "TRMV(z, z)" | cat trmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ztrmv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztrmv.c && \
    mv tmp.BLAS_ztrmv.c $@ && rm -f tmp.BLAS_ztrmv.c*

$(SRC_PATH)/BLAS_dtrmv_s.c: trmv.m4 trmv-common.m4 ../cblas.m4
	echo "TRMV(d, s)" | cat trmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dtrmv_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtrmv_s.c && \
    mv tmp.BLAS_dtrmv_s.c $@ && rm -f tmp.BLAS_dtrmv_s.c*

$(SRC_PATH)/BLAS_ztrmv_c.c: trmv.m4 trmv-common.m4 ../cblas.m4
	echo "TRMV(z, c)" | cat trmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ztrmv_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztrmv_c.c && \
    mv tmp.BLAS_ztrmv_c.c $@ && rm -f tmp.BLAS_ztrmv_c.c*

$(SRC_PATH)/BLAS_ctrmv_s.c: trmv.m4 trmv-common.m4 ../cblas.m4
	echo "TRMV(c, s)" | cat trmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ctrmv_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctrmv_s.c && \
    mv tmp.BLAS_ctrmv_s.c $@ && rm -f tmp.BLAS_ctrmv_s.c*

$(SRC_PATH)/BLAS_ztrmv_d.c: trmv.m4 trmv-common.m4 ../cblas.m4
	echo "TRMV(z, d)" | cat trmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ztrmv_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztrmv_d.c && \
    mv tmp.BLAS_ztrmv_d.c $@ && rm -f tmp.BLAS_ztrmv_d.c*

$(SRC_PATH)/BLAS_strmv_x.c: trmv.m4 trmv-common.m4 ../cblas.m4
	echo "TRMV(s, s, _x)" | cat trmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_strmv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_strmv_x.c && \
    mv tmp.BLAS_strmv_x.c $@ && rm -f tmp.BLAS_strmv_x.c*

$(SRC_PATH)/BLAS_dtrmv_x.c: trmv.m4 trmv-common.m4 ../cblas.m4
	echo "TRMV(d, d, _x)" | cat trmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dtrmv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtrmv_x.c && \
    mv tmp.BLAS_dtrmv_x.c $@ && rm -f tmp.BLAS_dtrmv_x.c*

$(SRC_PATH)/BLAS_ctrmv_x.c: trmv.m4 trmv-common.m4 ../cblas.m4
	echo "TRMV(c, c, _x)" | cat trmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ctrmv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctrmv_x.c && \
    mv tmp.BLAS_ctrmv_x.c $@ && rm -f tmp.BLAS_ctrmv_x.c*

$(SRC_PATH)/BLAS_ztrmv_x.c: trmv.m4 trmv-common.m4 ../cblas.m4
	echo "TRMV(z, z, _x)" | cat trmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ztrmv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztrmv_x.c && \
    mv tmp.BLAS_ztrmv_x.c $@ && rm -f tmp.BLAS_ztrmv_x.c*

$(SRC_PATH)/BLAS_dtrmv_s_x.c: trmv.m4 trmv-common.m4 ../cblas.m4
	echo "TRMV(d, s, _x)" | cat trmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dtrmv_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtrmv_s_x.c && \
    mv tmp.BLAS_dtrmv_s_x.c $@ && rm -f tmp.BLAS_dtrmv_s_x.c*

$(SRC_PATH)/BLAS_ztrmv_c_x.c: trmv.m4 trmv-common.m4 ../cblas.m4
	echo "TRMV(z, c, _x)" | cat trmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ztrmv_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztrmv_c_x.c && \
    mv tmp.BLAS_ztrmv_c_x.c $@ && rm -f tmp.BLAS_ztrmv_c_x.c*

$(SRC_PATH)/BLAS_ctrmv_s_x.c: trmv.m4 trmv-common.m4 ../cblas.m4
	echo "TRMV(c, s, _x)" | cat trmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ctrmv_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctrmv_s_x.c && \
    mv tmp.BLAS_ctrmv_s_x.c $@ && rm -f tmp.BLAS_ctrmv_s_x.c*

$(SRC_PATH)/BLAS_ztrmv_d_x.c: trmv.m4 trmv-common.m4 ../cblas.m4
	echo "TRMV(z, d, _x)" | cat trmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ztrmv_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztrmv_d_x.c && \
    mv tmp.BLAS_ztrmv_d_x.c $@ && rm -f tmp.BLAS_ztrmv_d_x.c*



$(SRC_PATH)/BLAS_strmv-f2c.c: ../f2c-bridge.m4 trmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trmv -DRARG="s, s" ../f2c-bridge.m4 > tmp.BLAS_strmv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_strmv-f2c.c && \
    mv tmp.BLAS_strmv-f2c.c $@ && rm -f tmp.BLAS_strmv-f2c.c*

$(SRC_PATH)/BLAS_dtrmv-f2c.c: ../f2c-bridge.m4 trmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trmv -DRARG="d, d" ../f2c-bridge.m4 > tmp.BLAS_dtrmv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtrmv-f2c.c && \
    mv tmp.BLAS_dtrmv-f2c.c $@ && rm -f tmp.BLAS_dtrmv-f2c.c*

$(SRC_PATH)/BLAS_ctrmv-f2c.c: ../f2c-bridge.m4 trmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trmv -DRARG="c, c" ../f2c-bridge.m4 > tmp.BLAS_ctrmv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctrmv-f2c.c && \
    mv tmp.BLAS_ctrmv-f2c.c $@ && rm -f tmp.BLAS_ctrmv-f2c.c*

$(SRC_PATH)/BLAS_ztrmv-f2c.c: ../f2c-bridge.m4 trmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trmv -DRARG="z, z" ../f2c-bridge.m4 > tmp.BLAS_ztrmv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztrmv-f2c.c && \
    mv tmp.BLAS_ztrmv-f2c.c $@ && rm -f tmp.BLAS_ztrmv-f2c.c*

$(SRC_PATH)/BLAS_dtrmv_s-f2c.c: ../f2c-bridge.m4 trmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trmv -DRARG="d, s" ../f2c-bridge.m4 > tmp.BLAS_dtrmv_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtrmv_s-f2c.c && \
    mv tmp.BLAS_dtrmv_s-f2c.c $@ && rm -f tmp.BLAS_dtrmv_s-f2c.c*

$(SRC_PATH)/BLAS_ztrmv_c-f2c.c: ../f2c-bridge.m4 trmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trmv -DRARG="z, c" ../f2c-bridge.m4 > tmp.BLAS_ztrmv_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztrmv_c-f2c.c && \
    mv tmp.BLAS_ztrmv_c-f2c.c $@ && rm -f tmp.BLAS_ztrmv_c-f2c.c*

$(SRC_PATH)/BLAS_ctrmv_s-f2c.c: ../f2c-bridge.m4 trmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trmv -DRARG="c, s" ../f2c-bridge.m4 > tmp.BLAS_ctrmv_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctrmv_s-f2c.c && \
    mv tmp.BLAS_ctrmv_s-f2c.c $@ && rm -f tmp.BLAS_ctrmv_s-f2c.c*

$(SRC_PATH)/BLAS_ztrmv_d-f2c.c: ../f2c-bridge.m4 trmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trmv -DRARG="z, d" ../f2c-bridge.m4 > tmp.BLAS_ztrmv_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztrmv_d-f2c.c && \
    mv tmp.BLAS_ztrmv_d-f2c.c $@ && rm -f tmp.BLAS_ztrmv_d-f2c.c*

$(SRC_PATH)/BLAS_strmv_x-f2c.c: ../f2c-bridge.m4 trmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trmv -DRARG="s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_strmv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_strmv_x-f2c.c && \
    mv tmp.BLAS_strmv_x-f2c.c $@ && rm -f tmp.BLAS_strmv_x-f2c.c*

$(SRC_PATH)/BLAS_dtrmv_x-f2c.c: ../f2c-bridge.m4 trmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trmv -DRARG="d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dtrmv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtrmv_x-f2c.c && \
    mv tmp.BLAS_dtrmv_x-f2c.c $@ && rm -f tmp.BLAS_dtrmv_x-f2c.c*

$(SRC_PATH)/BLAS_ctrmv_x-f2c.c: ../f2c-bridge.m4 trmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trmv -DRARG="c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_ctrmv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctrmv_x-f2c.c && \
    mv tmp.BLAS_ctrmv_x-f2c.c $@ && rm -f tmp.BLAS_ctrmv_x-f2c.c*

$(SRC_PATH)/BLAS_ztrmv_x-f2c.c: ../f2c-bridge.m4 trmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trmv -DRARG="z, z, _x" ../f2c-bridge.m4 > tmp.BLAS_ztrmv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztrmv_x-f2c.c && \
    mv tmp.BLAS_ztrmv_x-f2c.c $@ && rm -f tmp.BLAS_ztrmv_x-f2c.c*

$(SRC_PATH)/BLAS_dtrmv_s_x-f2c.c: ../f2c-bridge.m4 trmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trmv -DRARG="d, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dtrmv_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtrmv_s_x-f2c.c && \
    mv tmp.BLAS_dtrmv_s_x-f2c.c $@ && rm -f tmp.BLAS_dtrmv_s_x-f2c.c*

$(SRC_PATH)/BLAS_ztrmv_c_x-f2c.c: ../f2c-bridge.m4 trmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trmv -DRARG="z, c, _x" ../f2c-bridge.m4 > tmp.BLAS_ztrmv_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztrmv_c_x-f2c.c && \
    mv tmp.BLAS_ztrmv_c_x-f2c.c $@ && rm -f tmp.BLAS_ztrmv_c_x-f2c.c*

$(SRC_PATH)/BLAS_ctrmv_s_x-f2c.c: ../f2c-bridge.m4 trmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trmv -DRARG="c, s, _x" ../f2c-bridge.m4 > tmp.BLAS_ctrmv_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctrmv_s_x-f2c.c && \
    mv tmp.BLAS_ctrmv_s_x-f2c.c $@ && rm -f tmp.BLAS_ctrmv_s_x-f2c.c*

$(SRC_PATH)/BLAS_ztrmv_d_x-f2c.c: ../f2c-bridge.m4 trmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trmv -DRARG="z, d, _x" ../f2c-bridge.m4 > tmp.BLAS_ztrmv_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztrmv_d_x-f2c.c && \
    mv tmp.BLAS_ztrmv_d_x-f2c.c $@ && rm -f tmp.BLAS_ztrmv_d_x-f2c.c*


trmv.h: ../cblas.m4 ../header.template trmv-common.m4
	sed -e "s/routine/trmv/g" -e "s/ROUTINE/TRMV/g" \
  ../header.template | $(M4) $(M4_OPTS) $(M4_PATH) - >$@

source-clean:
	rm -f $(TRMV_SRCS) trmv.h

.PHONY: trmv-source trmv-header source-clean

