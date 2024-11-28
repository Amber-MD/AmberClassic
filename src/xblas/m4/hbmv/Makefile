include ../../make.conf
include ../../$(MAKEINC)


SRC_PATH = ../../src/hbmv
M4_PATH = -I. -I..

HBMV_SRCS = \
	$(SRC_PATH)/BLAS_chbmv.c\
	$(SRC_PATH)/BLAS_zhbmv.c\
	$(SRC_PATH)/BLAS_zhbmv_z_c.c\
	$(SRC_PATH)/BLAS_zhbmv_c_z.c\
	$(SRC_PATH)/BLAS_zhbmv_c_c.c\
	$(SRC_PATH)/BLAS_chbmv_c_s.c\
	$(SRC_PATH)/BLAS_zhbmv_z_d.c\
	$(SRC_PATH)/BLAS_chbmv_x.c\
	$(SRC_PATH)/BLAS_zhbmv_x.c\
	$(SRC_PATH)/BLAS_zhbmv_z_c_x.c\
	$(SRC_PATH)/BLAS_zhbmv_c_z_x.c\
	$(SRC_PATH)/BLAS_zhbmv_c_c_x.c\
	$(SRC_PATH)/BLAS_chbmv_c_s_x.c\
	$(SRC_PATH)/BLAS_zhbmv_z_d_x.c \
	\
	$(SRC_PATH)/BLAS_chbmv-f2c.c\
	$(SRC_PATH)/BLAS_zhbmv-f2c.c\
	$(SRC_PATH)/BLAS_zhbmv_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zhbmv_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zhbmv_c_c-f2c.c\
	$(SRC_PATH)/BLAS_chbmv_c_s-f2c.c\
	$(SRC_PATH)/BLAS_zhbmv_z_d-f2c.c\
	$(SRC_PATH)/BLAS_chbmv_x-f2c.c\
	$(SRC_PATH)/BLAS_zhbmv_x-f2c.c\
	$(SRC_PATH)/BLAS_zhbmv_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zhbmv_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zhbmv_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_chbmv_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zhbmv_z_d_x-f2c.c

hbmv-source: $(HBMV_SRCS)

hbmv-header: hbmv.h

$(SRC_PATH)/BLAS_chbmv.c: hbmv.m4 hbmv-common.m4 ../cblas.m4
	echo "HBMV(c, c, c)" | cat hbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_chbmv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chbmv.c && \
    mv tmp.BLAS_chbmv.c $@ && rm -f tmp.BLAS_chbmv.c*

$(SRC_PATH)/BLAS_zhbmv.c: hbmv.m4 hbmv-common.m4 ../cblas.m4
	echo "HBMV(z, z, z)" | cat hbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhbmv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhbmv.c && \
    mv tmp.BLAS_zhbmv.c $@ && rm -f tmp.BLAS_zhbmv.c*

$(SRC_PATH)/BLAS_zhbmv_z_c.c: hbmv.m4 hbmv-common.m4 ../cblas.m4
	echo "HBMV(z, z, c)" | cat hbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhbmv_z_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhbmv_z_c.c && \
    mv tmp.BLAS_zhbmv_z_c.c $@ && rm -f tmp.BLAS_zhbmv_z_c.c*

$(SRC_PATH)/BLAS_zhbmv_c_z.c: hbmv.m4 hbmv-common.m4 ../cblas.m4
	echo "HBMV(z, c, z)" | cat hbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhbmv_c_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhbmv_c_z.c && \
    mv tmp.BLAS_zhbmv_c_z.c $@ && rm -f tmp.BLAS_zhbmv_c_z.c*

$(SRC_PATH)/BLAS_zhbmv_c_c.c: hbmv.m4 hbmv-common.m4 ../cblas.m4
	echo "HBMV(z, c, c)" | cat hbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhbmv_c_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhbmv_c_c.c && \
    mv tmp.BLAS_zhbmv_c_c.c $@ && rm -f tmp.BLAS_zhbmv_c_c.c*

$(SRC_PATH)/BLAS_chbmv_c_s.c: hbmv.m4 hbmv-common.m4 ../cblas.m4
	echo "HBMV(c, c, s)" | cat hbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_chbmv_c_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chbmv_c_s.c && \
    mv tmp.BLAS_chbmv_c_s.c $@ && rm -f tmp.BLAS_chbmv_c_s.c*

$(SRC_PATH)/BLAS_zhbmv_z_d.c: hbmv.m4 hbmv-common.m4 ../cblas.m4
	echo "HBMV(z, z, d)" | cat hbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhbmv_z_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhbmv_z_d.c && \
    mv tmp.BLAS_zhbmv_z_d.c $@ && rm -f tmp.BLAS_zhbmv_z_d.c*

$(SRC_PATH)/BLAS_chbmv_x.c: hbmv.m4 hbmv-common.m4 ../cblas.m4
	echo "HBMV(c, c, c, _x)" | cat hbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_chbmv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chbmv_x.c && \
    mv tmp.BLAS_chbmv_x.c $@ && rm -f tmp.BLAS_chbmv_x.c*

$(SRC_PATH)/BLAS_zhbmv_x.c: hbmv.m4 hbmv-common.m4 ../cblas.m4
	echo "HBMV(z, z, z, _x)" | cat hbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhbmv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhbmv_x.c && \
    mv tmp.BLAS_zhbmv_x.c $@ && rm -f tmp.BLAS_zhbmv_x.c*

$(SRC_PATH)/BLAS_zhbmv_z_c_x.c: hbmv.m4 hbmv-common.m4 ../cblas.m4
	echo "HBMV(z, z, c, _x)" | cat hbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhbmv_z_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhbmv_z_c_x.c && \
    mv tmp.BLAS_zhbmv_z_c_x.c $@ && rm -f tmp.BLAS_zhbmv_z_c_x.c*

$(SRC_PATH)/BLAS_zhbmv_c_z_x.c: hbmv.m4 hbmv-common.m4 ../cblas.m4
	echo "HBMV(z, c, z, _x)" | cat hbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhbmv_c_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhbmv_c_z_x.c && \
    mv tmp.BLAS_zhbmv_c_z_x.c $@ && rm -f tmp.BLAS_zhbmv_c_z_x.c*

$(SRC_PATH)/BLAS_zhbmv_c_c_x.c: hbmv.m4 hbmv-common.m4 ../cblas.m4
	echo "HBMV(z, c, c, _x)" | cat hbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhbmv_c_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhbmv_c_c_x.c && \
    mv tmp.BLAS_zhbmv_c_c_x.c $@ && rm -f tmp.BLAS_zhbmv_c_c_x.c*

$(SRC_PATH)/BLAS_chbmv_c_s_x.c: hbmv.m4 hbmv-common.m4 ../cblas.m4
	echo "HBMV(c, c, s, _x)" | cat hbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_chbmv_c_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chbmv_c_s_x.c && \
    mv tmp.BLAS_chbmv_c_s_x.c $@ && rm -f tmp.BLAS_chbmv_c_s_x.c*

$(SRC_PATH)/BLAS_zhbmv_z_d_x.c: hbmv.m4 hbmv-common.m4 ../cblas.m4
	echo "HBMV(z, z, d, _x)" | cat hbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhbmv_z_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhbmv_z_d_x.c && \
    mv tmp.BLAS_zhbmv_z_d_x.c $@ && rm -f tmp.BLAS_zhbmv_z_d_x.c*



$(SRC_PATH)/BLAS_chbmv-f2c.c: ../f2c-bridge.m4 hbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hbmv -DRARG="c, c, c" ../f2c-bridge.m4 > tmp.BLAS_chbmv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chbmv-f2c.c && \
    mv tmp.BLAS_chbmv-f2c.c $@ && rm -f tmp.BLAS_chbmv-f2c.c*

$(SRC_PATH)/BLAS_zhbmv-f2c.c: ../f2c-bridge.m4 hbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hbmv -DRARG="z, z, z" ../f2c-bridge.m4 > tmp.BLAS_zhbmv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhbmv-f2c.c && \
    mv tmp.BLAS_zhbmv-f2c.c $@ && rm -f tmp.BLAS_zhbmv-f2c.c*

$(SRC_PATH)/BLAS_zhbmv_z_c-f2c.c: ../f2c-bridge.m4 hbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hbmv -DRARG="z, z, c" ../f2c-bridge.m4 > tmp.BLAS_zhbmv_z_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhbmv_z_c-f2c.c && \
    mv tmp.BLAS_zhbmv_z_c-f2c.c $@ && rm -f tmp.BLAS_zhbmv_z_c-f2c.c*

$(SRC_PATH)/BLAS_zhbmv_c_z-f2c.c: ../f2c-bridge.m4 hbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hbmv -DRARG="z, c, z" ../f2c-bridge.m4 > tmp.BLAS_zhbmv_c_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhbmv_c_z-f2c.c && \
    mv tmp.BLAS_zhbmv_c_z-f2c.c $@ && rm -f tmp.BLAS_zhbmv_c_z-f2c.c*

$(SRC_PATH)/BLAS_zhbmv_c_c-f2c.c: ../f2c-bridge.m4 hbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hbmv -DRARG="z, c, c" ../f2c-bridge.m4 > tmp.BLAS_zhbmv_c_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhbmv_c_c-f2c.c && \
    mv tmp.BLAS_zhbmv_c_c-f2c.c $@ && rm -f tmp.BLAS_zhbmv_c_c-f2c.c*

$(SRC_PATH)/BLAS_chbmv_c_s-f2c.c: ../f2c-bridge.m4 hbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hbmv -DRARG="c, c, s" ../f2c-bridge.m4 > tmp.BLAS_chbmv_c_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chbmv_c_s-f2c.c && \
    mv tmp.BLAS_chbmv_c_s-f2c.c $@ && rm -f tmp.BLAS_chbmv_c_s-f2c.c*

$(SRC_PATH)/BLAS_zhbmv_z_d-f2c.c: ../f2c-bridge.m4 hbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hbmv -DRARG="z, z, d" ../f2c-bridge.m4 > tmp.BLAS_zhbmv_z_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhbmv_z_d-f2c.c && \
    mv tmp.BLAS_zhbmv_z_d-f2c.c $@ && rm -f tmp.BLAS_zhbmv_z_d-f2c.c*

$(SRC_PATH)/BLAS_chbmv_x-f2c.c: ../f2c-bridge.m4 hbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hbmv -DRARG="c, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_chbmv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chbmv_x-f2c.c && \
    mv tmp.BLAS_chbmv_x-f2c.c $@ && rm -f tmp.BLAS_chbmv_x-f2c.c*

$(SRC_PATH)/BLAS_zhbmv_x-f2c.c: ../f2c-bridge.m4 hbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hbmv -DRARG="z, z, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zhbmv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhbmv_x-f2c.c && \
    mv tmp.BLAS_zhbmv_x-f2c.c $@ && rm -f tmp.BLAS_zhbmv_x-f2c.c*

$(SRC_PATH)/BLAS_zhbmv_z_c_x-f2c.c: ../f2c-bridge.m4 hbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hbmv -DRARG="z, z, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zhbmv_z_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhbmv_z_c_x-f2c.c && \
    mv tmp.BLAS_zhbmv_z_c_x-f2c.c $@ && rm -f tmp.BLAS_zhbmv_z_c_x-f2c.c*

$(SRC_PATH)/BLAS_zhbmv_c_z_x-f2c.c: ../f2c-bridge.m4 hbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hbmv -DRARG="z, c, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zhbmv_c_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhbmv_c_z_x-f2c.c && \
    mv tmp.BLAS_zhbmv_c_z_x-f2c.c $@ && rm -f tmp.BLAS_zhbmv_c_z_x-f2c.c*

$(SRC_PATH)/BLAS_zhbmv_c_c_x-f2c.c: ../f2c-bridge.m4 hbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hbmv -DRARG="z, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zhbmv_c_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhbmv_c_c_x-f2c.c && \
    mv tmp.BLAS_zhbmv_c_c_x-f2c.c $@ && rm -f tmp.BLAS_zhbmv_c_c_x-f2c.c*

$(SRC_PATH)/BLAS_chbmv_c_s_x-f2c.c: ../f2c-bridge.m4 hbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hbmv -DRARG="c, c, s, _x" ../f2c-bridge.m4 > tmp.BLAS_chbmv_c_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chbmv_c_s_x-f2c.c && \
    mv tmp.BLAS_chbmv_c_s_x-f2c.c $@ && rm -f tmp.BLAS_chbmv_c_s_x-f2c.c*

$(SRC_PATH)/BLAS_zhbmv_z_d_x-f2c.c: ../f2c-bridge.m4 hbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hbmv -DRARG="z, z, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zhbmv_z_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhbmv_z_d_x-f2c.c && \
    mv tmp.BLAS_zhbmv_z_d_x-f2c.c $@ && rm -f tmp.BLAS_zhbmv_z_d_x-f2c.c*


hbmv.h: ../cblas.m4 ../header.template hbmv-common.m4
	sed -e "s/routine/hbmv/g" -e "s/ROUTINE/HBMV/g" \
  ../header.template | $(M4) $(M4_OPTS) $(M4_PATH) - >$@

source-clean:
	rm -f $(HBMV_SRCS) hbmv.h

.PHONY: hbmv-source hbmv-header source-clean

