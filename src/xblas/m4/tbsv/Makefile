include ../../make.conf
include ../../$(MAKEINC)


SRC_PATH = ../../src/tbsv
M4_PATH = -I. -I..

TBSV_SRCS = \
	$(SRC_PATH)/BLAS_stbsv.c\
	$(SRC_PATH)/BLAS_dtbsv.c\
	$(SRC_PATH)/BLAS_ctbsv.c\
	$(SRC_PATH)/BLAS_ztbsv.c\
	$(SRC_PATH)/BLAS_dtbsv_s.c\
	$(SRC_PATH)/BLAS_ztbsv_c.c\
	$(SRC_PATH)/BLAS_ctbsv_s.c\
	$(SRC_PATH)/BLAS_ztbsv_d.c\
	$(SRC_PATH)/BLAS_stbsv_x.c\
	$(SRC_PATH)/BLAS_dtbsv_x.c\
	$(SRC_PATH)/BLAS_dtbsv_s_x.c\
	$(SRC_PATH)/BLAS_ctbsv_x.c\
	$(SRC_PATH)/BLAS_ztbsv_x.c\
	$(SRC_PATH)/BLAS_ztbsv_c_x.c\
	$(SRC_PATH)/BLAS_ctbsv_s_x.c\
	$(SRC_PATH)/BLAS_ztbsv_d_x.c \
	\
	$(SRC_PATH)/BLAS_stbsv-f2c.c\
	$(SRC_PATH)/BLAS_dtbsv-f2c.c\
	$(SRC_PATH)/BLAS_ctbsv-f2c.c\
	$(SRC_PATH)/BLAS_ztbsv-f2c.c\
	$(SRC_PATH)/BLAS_dtbsv_s-f2c.c\
	$(SRC_PATH)/BLAS_ztbsv_c-f2c.c\
	$(SRC_PATH)/BLAS_ctbsv_s-f2c.c\
	$(SRC_PATH)/BLAS_ztbsv_d-f2c.c\
	$(SRC_PATH)/BLAS_stbsv_x-f2c.c\
	$(SRC_PATH)/BLAS_dtbsv_x-f2c.c\
	$(SRC_PATH)/BLAS_dtbsv_s_x-f2c.c\
	$(SRC_PATH)/BLAS_ctbsv_x-f2c.c\
	$(SRC_PATH)/BLAS_ztbsv_x-f2c.c\
	$(SRC_PATH)/BLAS_ztbsv_c_x-f2c.c\
	$(SRC_PATH)/BLAS_ctbsv_s_x-f2c.c\
	$(SRC_PATH)/BLAS_ztbsv_d_x-f2c.c

tbsv-source: $(TBSV_SRCS)

tbsv-header: tbsv.h

$(SRC_PATH)/BLAS_stbsv.c: tbsv.m4 tbsv-common.m4 ../cblas.m4
	echo "TBSV(s, s)" | cat tbsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_stbsv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_stbsv.c && \
    mv tmp.BLAS_stbsv.c $@ && rm -f tmp.BLAS_stbsv.c*

$(SRC_PATH)/BLAS_dtbsv.c: tbsv.m4 tbsv-common.m4 ../cblas.m4
	echo "TBSV(d, d)" | cat tbsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dtbsv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtbsv.c && \
    mv tmp.BLAS_dtbsv.c $@ && rm -f tmp.BLAS_dtbsv.c*

$(SRC_PATH)/BLAS_ctbsv.c: tbsv.m4 tbsv-common.m4 ../cblas.m4
	echo "TBSV(c, c)" | cat tbsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ctbsv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctbsv.c && \
    mv tmp.BLAS_ctbsv.c $@ && rm -f tmp.BLAS_ctbsv.c*

$(SRC_PATH)/BLAS_ztbsv.c: tbsv.m4 tbsv-common.m4 ../cblas.m4
	echo "TBSV(z, z)" | cat tbsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ztbsv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztbsv.c && \
    mv tmp.BLAS_ztbsv.c $@ && rm -f tmp.BLAS_ztbsv.c*

$(SRC_PATH)/BLAS_dtbsv_s.c: tbsv.m4 tbsv-common.m4 ../cblas.m4
	echo "TBSV(d, s)" | cat tbsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dtbsv_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtbsv_s.c && \
    mv tmp.BLAS_dtbsv_s.c $@ && rm -f tmp.BLAS_dtbsv_s.c*

$(SRC_PATH)/BLAS_ztbsv_c.c: tbsv.m4 tbsv-common.m4 ../cblas.m4
	echo "TBSV(z, c)" | cat tbsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ztbsv_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztbsv_c.c && \
    mv tmp.BLAS_ztbsv_c.c $@ && rm -f tmp.BLAS_ztbsv_c.c*

$(SRC_PATH)/BLAS_ctbsv_s.c: tbsv.m4 tbsv-common.m4 ../cblas.m4
	echo "TBSV(c, s)" | cat tbsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ctbsv_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctbsv_s.c && \
    mv tmp.BLAS_ctbsv_s.c $@ && rm -f tmp.BLAS_ctbsv_s.c*

$(SRC_PATH)/BLAS_ztbsv_d.c: tbsv.m4 tbsv-common.m4 ../cblas.m4
	echo "TBSV(z, d)" | cat tbsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ztbsv_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztbsv_d.c && \
    mv tmp.BLAS_ztbsv_d.c $@ && rm -f tmp.BLAS_ztbsv_d.c*

$(SRC_PATH)/BLAS_stbsv_x.c: tbsv.m4 tbsv-common.m4 ../cblas.m4
	echo "TBSV(s, s, _x)" | cat tbsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_stbsv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_stbsv_x.c && \
    mv tmp.BLAS_stbsv_x.c $@ && rm -f tmp.BLAS_stbsv_x.c*

$(SRC_PATH)/BLAS_dtbsv_x.c: tbsv.m4 tbsv-common.m4 ../cblas.m4
	echo "TBSV(d, d, _x)" | cat tbsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dtbsv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtbsv_x.c && \
    mv tmp.BLAS_dtbsv_x.c $@ && rm -f tmp.BLAS_dtbsv_x.c*

$(SRC_PATH)/BLAS_dtbsv_s_x.c: tbsv.m4 tbsv-common.m4 ../cblas.m4
	echo "TBSV(d, s, _x)" | cat tbsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dtbsv_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtbsv_s_x.c && \
    mv tmp.BLAS_dtbsv_s_x.c $@ && rm -f tmp.BLAS_dtbsv_s_x.c*

$(SRC_PATH)/BLAS_ctbsv_x.c: tbsv.m4 tbsv-common.m4 ../cblas.m4
	echo "TBSV(c, c, _x)" | cat tbsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ctbsv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctbsv_x.c && \
    mv tmp.BLAS_ctbsv_x.c $@ && rm -f tmp.BLAS_ctbsv_x.c*

$(SRC_PATH)/BLAS_ztbsv_x.c: tbsv.m4 tbsv-common.m4 ../cblas.m4
	echo "TBSV(z, z, _x)" | cat tbsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ztbsv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztbsv_x.c && \
    mv tmp.BLAS_ztbsv_x.c $@ && rm -f tmp.BLAS_ztbsv_x.c*

$(SRC_PATH)/BLAS_ztbsv_c_x.c: tbsv.m4 tbsv-common.m4 ../cblas.m4
	echo "TBSV(z, c, _x)" | cat tbsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ztbsv_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztbsv_c_x.c && \
    mv tmp.BLAS_ztbsv_c_x.c $@ && rm -f tmp.BLAS_ztbsv_c_x.c*

$(SRC_PATH)/BLAS_ctbsv_s_x.c: tbsv.m4 tbsv-common.m4 ../cblas.m4
	echo "TBSV(c, s, _x)" | cat tbsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ctbsv_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctbsv_s_x.c && \
    mv tmp.BLAS_ctbsv_s_x.c $@ && rm -f tmp.BLAS_ctbsv_s_x.c*

$(SRC_PATH)/BLAS_ztbsv_d_x.c: tbsv.m4 tbsv-common.m4 ../cblas.m4
	echo "TBSV(z, d, _x)" | cat tbsv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ztbsv_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztbsv_d_x.c && \
    mv tmp.BLAS_ztbsv_d_x.c $@ && rm -f tmp.BLAS_ztbsv_d_x.c*



$(SRC_PATH)/BLAS_stbsv-f2c.c: ../f2c-bridge.m4 tbsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tbsv -DRARG="s, s" ../f2c-bridge.m4 > tmp.BLAS_stbsv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_stbsv-f2c.c && \
    mv tmp.BLAS_stbsv-f2c.c $@ && rm -f tmp.BLAS_stbsv-f2c.c*

$(SRC_PATH)/BLAS_dtbsv-f2c.c: ../f2c-bridge.m4 tbsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tbsv -DRARG="d, d" ../f2c-bridge.m4 > tmp.BLAS_dtbsv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtbsv-f2c.c && \
    mv tmp.BLAS_dtbsv-f2c.c $@ && rm -f tmp.BLAS_dtbsv-f2c.c*

$(SRC_PATH)/BLAS_ctbsv-f2c.c: ../f2c-bridge.m4 tbsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tbsv -DRARG="c, c" ../f2c-bridge.m4 > tmp.BLAS_ctbsv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctbsv-f2c.c && \
    mv tmp.BLAS_ctbsv-f2c.c $@ && rm -f tmp.BLAS_ctbsv-f2c.c*

$(SRC_PATH)/BLAS_ztbsv-f2c.c: ../f2c-bridge.m4 tbsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tbsv -DRARG="z, z" ../f2c-bridge.m4 > tmp.BLAS_ztbsv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztbsv-f2c.c && \
    mv tmp.BLAS_ztbsv-f2c.c $@ && rm -f tmp.BLAS_ztbsv-f2c.c*

$(SRC_PATH)/BLAS_dtbsv_s-f2c.c: ../f2c-bridge.m4 tbsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tbsv -DRARG="d, s" ../f2c-bridge.m4 > tmp.BLAS_dtbsv_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtbsv_s-f2c.c && \
    mv tmp.BLAS_dtbsv_s-f2c.c $@ && rm -f tmp.BLAS_dtbsv_s-f2c.c*

$(SRC_PATH)/BLAS_ztbsv_c-f2c.c: ../f2c-bridge.m4 tbsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tbsv -DRARG="z, c" ../f2c-bridge.m4 > tmp.BLAS_ztbsv_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztbsv_c-f2c.c && \
    mv tmp.BLAS_ztbsv_c-f2c.c $@ && rm -f tmp.BLAS_ztbsv_c-f2c.c*

$(SRC_PATH)/BLAS_ctbsv_s-f2c.c: ../f2c-bridge.m4 tbsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tbsv -DRARG="c, s" ../f2c-bridge.m4 > tmp.BLAS_ctbsv_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctbsv_s-f2c.c && \
    mv tmp.BLAS_ctbsv_s-f2c.c $@ && rm -f tmp.BLAS_ctbsv_s-f2c.c*

$(SRC_PATH)/BLAS_ztbsv_d-f2c.c: ../f2c-bridge.m4 tbsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tbsv -DRARG="z, d" ../f2c-bridge.m4 > tmp.BLAS_ztbsv_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztbsv_d-f2c.c && \
    mv tmp.BLAS_ztbsv_d-f2c.c $@ && rm -f tmp.BLAS_ztbsv_d-f2c.c*

$(SRC_PATH)/BLAS_stbsv_x-f2c.c: ../f2c-bridge.m4 tbsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tbsv -DRARG="s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_stbsv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_stbsv_x-f2c.c && \
    mv tmp.BLAS_stbsv_x-f2c.c $@ && rm -f tmp.BLAS_stbsv_x-f2c.c*

$(SRC_PATH)/BLAS_dtbsv_x-f2c.c: ../f2c-bridge.m4 tbsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tbsv -DRARG="d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dtbsv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtbsv_x-f2c.c && \
    mv tmp.BLAS_dtbsv_x-f2c.c $@ && rm -f tmp.BLAS_dtbsv_x-f2c.c*

$(SRC_PATH)/BLAS_dtbsv_s_x-f2c.c: ../f2c-bridge.m4 tbsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tbsv -DRARG="d, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dtbsv_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtbsv_s_x-f2c.c && \
    mv tmp.BLAS_dtbsv_s_x-f2c.c $@ && rm -f tmp.BLAS_dtbsv_s_x-f2c.c*

$(SRC_PATH)/BLAS_ctbsv_x-f2c.c: ../f2c-bridge.m4 tbsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tbsv -DRARG="c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_ctbsv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctbsv_x-f2c.c && \
    mv tmp.BLAS_ctbsv_x-f2c.c $@ && rm -f tmp.BLAS_ctbsv_x-f2c.c*

$(SRC_PATH)/BLAS_ztbsv_x-f2c.c: ../f2c-bridge.m4 tbsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tbsv -DRARG="z, z, _x" ../f2c-bridge.m4 > tmp.BLAS_ztbsv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztbsv_x-f2c.c && \
    mv tmp.BLAS_ztbsv_x-f2c.c $@ && rm -f tmp.BLAS_ztbsv_x-f2c.c*

$(SRC_PATH)/BLAS_ztbsv_c_x-f2c.c: ../f2c-bridge.m4 tbsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tbsv -DRARG="z, c, _x" ../f2c-bridge.m4 > tmp.BLAS_ztbsv_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztbsv_c_x-f2c.c && \
    mv tmp.BLAS_ztbsv_c_x-f2c.c $@ && rm -f tmp.BLAS_ztbsv_c_x-f2c.c*

$(SRC_PATH)/BLAS_ctbsv_s_x-f2c.c: ../f2c-bridge.m4 tbsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tbsv -DRARG="c, s, _x" ../f2c-bridge.m4 > tmp.BLAS_ctbsv_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctbsv_s_x-f2c.c && \
    mv tmp.BLAS_ctbsv_s_x-f2c.c $@ && rm -f tmp.BLAS_ctbsv_s_x-f2c.c*

$(SRC_PATH)/BLAS_ztbsv_d_x-f2c.c: ../f2c-bridge.m4 tbsv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tbsv -DRARG="z, d, _x" ../f2c-bridge.m4 > tmp.BLAS_ztbsv_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztbsv_d_x-f2c.c && \
    mv tmp.BLAS_ztbsv_d_x-f2c.c $@ && rm -f tmp.BLAS_ztbsv_d_x-f2c.c*


tbsv.h: ../cblas.m4 ../header.template tbsv-common.m4
	sed -e "s/routine/tbsv/g" -e "s/ROUTINE/TBSV/g" \
  ../header.template | $(M4) $(M4_OPTS) $(M4_PATH) - >$@

source-clean:
	rm -f $(TBSV_SRCS) tbsv.h

.PHONY: tbsv-source tbsv-header source-clean

