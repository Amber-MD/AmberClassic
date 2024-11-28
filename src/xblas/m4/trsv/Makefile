include ../../make.conf
include ../../$(MAKEINC)


SRC_PATH = ../../src/trsv
M4_PATH = -I. -I..

TRSV_SRCS = \
	$(SRC_PATH)/BLAS_strsv.c\
	$(SRC_PATH)/BLAS_dtrsv.c\
	$(SRC_PATH)/BLAS_ctrsv.c\
	$(SRC_PATH)/BLAS_ztrsv.c\
	$(SRC_PATH)/BLAS_dtrsv_s.c\
	$(SRC_PATH)/BLAS_ztrsv_c.c\
	$(SRC_PATH)/BLAS_ctrsv_s.c\
	$(SRC_PATH)/BLAS_ztrsv_d.c\
	$(SRC_PATH)/BLAS_strsv_x.c\
	$(SRC_PATH)/BLAS_dtrsv_x.c\
	$(SRC_PATH)/BLAS_dtrsv_s_x.c\
	$(SRC_PATH)/BLAS_ctrsv_x.c\
	$(SRC_PATH)/BLAS_ztrsv_x.c\
	$(SRC_PATH)/BLAS_ztrsv_c_x.c\
	$(SRC_PATH)/BLAS_ctrsv_s_x.c\
	$(SRC_PATH)/BLAS_ztrsv_d_x.c \
	\
	$(SRC_PATH)/BLAS_strsv-f2c.c\
	$(SRC_PATH)/BLAS_dtrsv-f2c.c\
	$(SRC_PATH)/BLAS_ctrsv-f2c.c\
	$(SRC_PATH)/BLAS_ztrsv-f2c.c\
	$(SRC_PATH)/BLAS_dtrsv_s-f2c.c\
	$(SRC_PATH)/BLAS_ztrsv_c-f2c.c\
	$(SRC_PATH)/BLAS_ctrsv_s-f2c.c\
	$(SRC_PATH)/BLAS_ztrsv_d-f2c.c\
	$(SRC_PATH)/BLAS_strsv_x-f2c.c\
	$(SRC_PATH)/BLAS_dtrsv_x-f2c.c\
	$(SRC_PATH)/BLAS_dtrsv_s_x-f2c.c\
	$(SRC_PATH)/BLAS_ctrsv_x-f2c.c\
	$(SRC_PATH)/BLAS_ztrsv_x-f2c.c\
	$(SRC_PATH)/BLAS_ztrsv_c_x-f2c.c\
	$(SRC_PATH)/BLAS_ctrsv_s_x-f2c.c\
	$(SRC_PATH)/BLAS_ztrsv_d_x-f2c.c

trsv-source: $(TRSV_SRCS)

trsv-header: trsv.h

$(SRC_PATH)/BLAS_strsv.c: trsv.m4 trsv-common.m4 ../cblas.m4
	echo "TRSV(s, s)" | cat trsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_strsv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_strsv.c && \
    mv tmp.BLAS_strsv.c $@ && rm -f tmp.BLAS_strsv.c*

$(SRC_PATH)/BLAS_dtrsv.c: trsv.m4 trsv-common.m4 ../cblas.m4
	echo "TRSV(d, d)" | cat trsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dtrsv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtrsv.c && \
    mv tmp.BLAS_dtrsv.c $@ && rm -f tmp.BLAS_dtrsv.c*

$(SRC_PATH)/BLAS_ctrsv.c: trsv.m4 trsv-common.m4 ../cblas.m4
	echo "TRSV(c, c)" | cat trsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ctrsv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctrsv.c && \
    mv tmp.BLAS_ctrsv.c $@ && rm -f tmp.BLAS_ctrsv.c*

$(SRC_PATH)/BLAS_ztrsv.c: trsv.m4 trsv-common.m4 ../cblas.m4
	echo "TRSV(z, z)" | cat trsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ztrsv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztrsv.c && \
    mv tmp.BLAS_ztrsv.c $@ && rm -f tmp.BLAS_ztrsv.c*

$(SRC_PATH)/BLAS_dtrsv_s.c: trsv.m4 trsv-common.m4 ../cblas.m4
	echo "TRSV(d, s)" | cat trsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dtrsv_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtrsv_s.c && \
    mv tmp.BLAS_dtrsv_s.c $@ && rm -f tmp.BLAS_dtrsv_s.c*

$(SRC_PATH)/BLAS_ztrsv_c.c: trsv.m4 trsv-common.m4 ../cblas.m4
	echo "TRSV(z, c)" | cat trsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ztrsv_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztrsv_c.c && \
    mv tmp.BLAS_ztrsv_c.c $@ && rm -f tmp.BLAS_ztrsv_c.c*

$(SRC_PATH)/BLAS_ctrsv_s.c: trsv.m4 trsv-common.m4 ../cblas.m4
	echo "TRSV(c, s)" | cat trsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ctrsv_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctrsv_s.c && \
    mv tmp.BLAS_ctrsv_s.c $@ && rm -f tmp.BLAS_ctrsv_s.c*

$(SRC_PATH)/BLAS_ztrsv_d.c: trsv.m4 trsv-common.m4 ../cblas.m4
	echo "TRSV(z, d)" | cat trsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ztrsv_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztrsv_d.c && \
    mv tmp.BLAS_ztrsv_d.c $@ && rm -f tmp.BLAS_ztrsv_d.c*

$(SRC_PATH)/BLAS_strsv_x.c: trsv.m4 trsv-common.m4 ../cblas.m4
	echo "TRSV(s, s, _x)" | cat trsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_strsv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_strsv_x.c && \
    mv tmp.BLAS_strsv_x.c $@ && rm -f tmp.BLAS_strsv_x.c*

$(SRC_PATH)/BLAS_dtrsv_x.c: trsv.m4 trsv-common.m4 ../cblas.m4
	echo "TRSV(d, d, _x)" | cat trsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dtrsv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtrsv_x.c && \
    mv tmp.BLAS_dtrsv_x.c $@ && rm -f tmp.BLAS_dtrsv_x.c*

$(SRC_PATH)/BLAS_dtrsv_s_x.c: trsv.m4 trsv-common.m4 ../cblas.m4
	echo "TRSV(d, s, _x)" | cat trsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dtrsv_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtrsv_s_x.c && \
    mv tmp.BLAS_dtrsv_s_x.c $@ && rm -f tmp.BLAS_dtrsv_s_x.c*

$(SRC_PATH)/BLAS_ctrsv_x.c: trsv.m4 trsv-common.m4 ../cblas.m4
	echo "TRSV(c, c, _x)" | cat trsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ctrsv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctrsv_x.c && \
    mv tmp.BLAS_ctrsv_x.c $@ && rm -f tmp.BLAS_ctrsv_x.c*

$(SRC_PATH)/BLAS_ztrsv_x.c: trsv.m4 trsv-common.m4 ../cblas.m4
	echo "TRSV(z, z, _x)" | cat trsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ztrsv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztrsv_x.c && \
    mv tmp.BLAS_ztrsv_x.c $@ && rm -f tmp.BLAS_ztrsv_x.c*

$(SRC_PATH)/BLAS_ztrsv_c_x.c: trsv.m4 trsv-common.m4 ../cblas.m4
	echo "TRSV(z, c, _x)" | cat trsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ztrsv_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztrsv_c_x.c && \
    mv tmp.BLAS_ztrsv_c_x.c $@ && rm -f tmp.BLAS_ztrsv_c_x.c*

$(SRC_PATH)/BLAS_ctrsv_s_x.c: trsv.m4 trsv-common.m4 ../cblas.m4
	echo "TRSV(c, s, _x)" | cat trsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ctrsv_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctrsv_s_x.c && \
    mv tmp.BLAS_ctrsv_s_x.c $@ && rm -f tmp.BLAS_ctrsv_s_x.c*

$(SRC_PATH)/BLAS_ztrsv_d_x.c: trsv.m4 trsv-common.m4 ../cblas.m4
	echo "TRSV(z, d, _x)" | cat trsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ztrsv_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztrsv_d_x.c && \
    mv tmp.BLAS_ztrsv_d_x.c $@ && rm -f tmp.BLAS_ztrsv_d_x.c*



$(SRC_PATH)/BLAS_strsv-f2c.c: ../f2c-bridge.m4 trsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trsv -DRARG="s, s" ../f2c-bridge.m4 > tmp.BLAS_strsv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_strsv-f2c.c && \
    mv tmp.BLAS_strsv-f2c.c $@ && rm -f tmp.BLAS_strsv-f2c.c*

$(SRC_PATH)/BLAS_dtrsv-f2c.c: ../f2c-bridge.m4 trsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trsv -DRARG="d, d" ../f2c-bridge.m4 > tmp.BLAS_dtrsv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtrsv-f2c.c && \
    mv tmp.BLAS_dtrsv-f2c.c $@ && rm -f tmp.BLAS_dtrsv-f2c.c*

$(SRC_PATH)/BLAS_ctrsv-f2c.c: ../f2c-bridge.m4 trsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trsv -DRARG="c, c" ../f2c-bridge.m4 > tmp.BLAS_ctrsv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctrsv-f2c.c && \
    mv tmp.BLAS_ctrsv-f2c.c $@ && rm -f tmp.BLAS_ctrsv-f2c.c*

$(SRC_PATH)/BLAS_ztrsv-f2c.c: ../f2c-bridge.m4 trsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trsv -DRARG="z, z" ../f2c-bridge.m4 > tmp.BLAS_ztrsv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztrsv-f2c.c && \
    mv tmp.BLAS_ztrsv-f2c.c $@ && rm -f tmp.BLAS_ztrsv-f2c.c*

$(SRC_PATH)/BLAS_dtrsv_s-f2c.c: ../f2c-bridge.m4 trsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trsv -DRARG="d, s" ../f2c-bridge.m4 > tmp.BLAS_dtrsv_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtrsv_s-f2c.c && \
    mv tmp.BLAS_dtrsv_s-f2c.c $@ && rm -f tmp.BLAS_dtrsv_s-f2c.c*

$(SRC_PATH)/BLAS_ztrsv_c-f2c.c: ../f2c-bridge.m4 trsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trsv -DRARG="z, c" ../f2c-bridge.m4 > tmp.BLAS_ztrsv_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztrsv_c-f2c.c && \
    mv tmp.BLAS_ztrsv_c-f2c.c $@ && rm -f tmp.BLAS_ztrsv_c-f2c.c*

$(SRC_PATH)/BLAS_ctrsv_s-f2c.c: ../f2c-bridge.m4 trsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trsv -DRARG="c, s" ../f2c-bridge.m4 > tmp.BLAS_ctrsv_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctrsv_s-f2c.c && \
    mv tmp.BLAS_ctrsv_s-f2c.c $@ && rm -f tmp.BLAS_ctrsv_s-f2c.c*

$(SRC_PATH)/BLAS_ztrsv_d-f2c.c: ../f2c-bridge.m4 trsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trsv -DRARG="z, d" ../f2c-bridge.m4 > tmp.BLAS_ztrsv_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztrsv_d-f2c.c && \
    mv tmp.BLAS_ztrsv_d-f2c.c $@ && rm -f tmp.BLAS_ztrsv_d-f2c.c*

$(SRC_PATH)/BLAS_strsv_x-f2c.c: ../f2c-bridge.m4 trsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trsv -DRARG="s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_strsv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_strsv_x-f2c.c && \
    mv tmp.BLAS_strsv_x-f2c.c $@ && rm -f tmp.BLAS_strsv_x-f2c.c*

$(SRC_PATH)/BLAS_dtrsv_x-f2c.c: ../f2c-bridge.m4 trsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trsv -DRARG="d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dtrsv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtrsv_x-f2c.c && \
    mv tmp.BLAS_dtrsv_x-f2c.c $@ && rm -f tmp.BLAS_dtrsv_x-f2c.c*

$(SRC_PATH)/BLAS_dtrsv_s_x-f2c.c: ../f2c-bridge.m4 trsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trsv -DRARG="d, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dtrsv_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtrsv_s_x-f2c.c && \
    mv tmp.BLAS_dtrsv_s_x-f2c.c $@ && rm -f tmp.BLAS_dtrsv_s_x-f2c.c*

$(SRC_PATH)/BLAS_ctrsv_x-f2c.c: ../f2c-bridge.m4 trsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trsv -DRARG="c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_ctrsv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctrsv_x-f2c.c && \
    mv tmp.BLAS_ctrsv_x-f2c.c $@ && rm -f tmp.BLAS_ctrsv_x-f2c.c*

$(SRC_PATH)/BLAS_ztrsv_x-f2c.c: ../f2c-bridge.m4 trsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trsv -DRARG="z, z, _x" ../f2c-bridge.m4 > tmp.BLAS_ztrsv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztrsv_x-f2c.c && \
    mv tmp.BLAS_ztrsv_x-f2c.c $@ && rm -f tmp.BLAS_ztrsv_x-f2c.c*

$(SRC_PATH)/BLAS_ztrsv_c_x-f2c.c: ../f2c-bridge.m4 trsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trsv -DRARG="z, c, _x" ../f2c-bridge.m4 > tmp.BLAS_ztrsv_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztrsv_c_x-f2c.c && \
    mv tmp.BLAS_ztrsv_c_x-f2c.c $@ && rm -f tmp.BLAS_ztrsv_c_x-f2c.c*

$(SRC_PATH)/BLAS_ctrsv_s_x-f2c.c: ../f2c-bridge.m4 trsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trsv -DRARG="c, s, _x" ../f2c-bridge.m4 > tmp.BLAS_ctrsv_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctrsv_s_x-f2c.c && \
    mv tmp.BLAS_ctrsv_s_x-f2c.c $@ && rm -f tmp.BLAS_ctrsv_s_x-f2c.c*

$(SRC_PATH)/BLAS_ztrsv_d_x-f2c.c: ../f2c-bridge.m4 trsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=trsv -DRARG="z, d, _x" ../f2c-bridge.m4 > tmp.BLAS_ztrsv_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztrsv_d_x-f2c.c && \
    mv tmp.BLAS_ztrsv_d_x-f2c.c $@ && rm -f tmp.BLAS_ztrsv_d_x-f2c.c*


trsv.h: ../cblas.m4 ../header.template trsv-common.m4
	sed -e "s/routine/trsv/g" -e "s/ROUTINE/TRSV/g" \
  ../header.template | $(M4) $(M4_OPTS) $(M4_PATH) - >$@

source-clean:
	rm -f $(TRSV_SRCS) trsv.h

.PHONY: trsv-source trsv-header source-clean

