include ../../make.conf
include ../../$(MAKEINC)


SRC_PATH = ../../src/hpmv
M4_PATH = -I. -I..

HPMV_SRCS = \
	$(SRC_PATH)/BLAS_chpmv.c\
	$(SRC_PATH)/BLAS_zhpmv.c\
	$(SRC_PATH)/BLAS_zhpmv_z_c.c\
	$(SRC_PATH)/BLAS_zhpmv_c_z.c\
	$(SRC_PATH)/BLAS_zhpmv_c_c.c\
	$(SRC_PATH)/BLAS_chpmv_c_s.c\
	$(SRC_PATH)/BLAS_zhpmv_z_d.c\
	$(SRC_PATH)/BLAS_chpmv_x.c\
	$(SRC_PATH)/BLAS_zhpmv_x.c\
	$(SRC_PATH)/BLAS_zhpmv_z_c_x.c\
	$(SRC_PATH)/BLAS_zhpmv_c_z_x.c\
	$(SRC_PATH)/BLAS_zhpmv_c_c_x.c\
	$(SRC_PATH)/BLAS_chpmv_c_s_x.c\
	$(SRC_PATH)/BLAS_zhpmv_z_d_x.c \
	\
	$(SRC_PATH)/BLAS_chpmv-f2c.c\
	$(SRC_PATH)/BLAS_zhpmv-f2c.c\
	$(SRC_PATH)/BLAS_zhpmv_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zhpmv_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zhpmv_c_c-f2c.c\
	$(SRC_PATH)/BLAS_chpmv_c_s-f2c.c\
	$(SRC_PATH)/BLAS_zhpmv_z_d-f2c.c\
	$(SRC_PATH)/BLAS_chpmv_x-f2c.c\
	$(SRC_PATH)/BLAS_zhpmv_x-f2c.c\
	$(SRC_PATH)/BLAS_zhpmv_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zhpmv_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zhpmv_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_chpmv_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zhpmv_z_d_x-f2c.c

hpmv-source: $(HPMV_SRCS)

hpmv-header: hpmv.h

$(SRC_PATH)/BLAS_chpmv.c: hpmv.m4 hpmv-common.m4 ../cblas.m4
	echo "HPMV(c, c, c)" | cat hpmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_chpmv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chpmv.c && \
    mv tmp.BLAS_chpmv.c $@ && rm -f tmp.BLAS_chpmv.c*

$(SRC_PATH)/BLAS_zhpmv.c: hpmv.m4 hpmv-common.m4 ../cblas.m4
	echo "HPMV(z, z, z)" | cat hpmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhpmv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhpmv.c && \
    mv tmp.BLAS_zhpmv.c $@ && rm -f tmp.BLAS_zhpmv.c*

$(SRC_PATH)/BLAS_zhpmv_z_c.c: hpmv.m4 hpmv-common.m4 ../cblas.m4
	echo "HPMV(z, z, c)" | cat hpmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhpmv_z_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhpmv_z_c.c && \
    mv tmp.BLAS_zhpmv_z_c.c $@ && rm -f tmp.BLAS_zhpmv_z_c.c*

$(SRC_PATH)/BLAS_zhpmv_c_z.c: hpmv.m4 hpmv-common.m4 ../cblas.m4
	echo "HPMV(z, c, z)" | cat hpmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhpmv_c_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhpmv_c_z.c && \
    mv tmp.BLAS_zhpmv_c_z.c $@ && rm -f tmp.BLAS_zhpmv_c_z.c*

$(SRC_PATH)/BLAS_zhpmv_c_c.c: hpmv.m4 hpmv-common.m4 ../cblas.m4
	echo "HPMV(z, c, c)" | cat hpmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhpmv_c_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhpmv_c_c.c && \
    mv tmp.BLAS_zhpmv_c_c.c $@ && rm -f tmp.BLAS_zhpmv_c_c.c*

$(SRC_PATH)/BLAS_chpmv_c_s.c: hpmv.m4 hpmv-common.m4 ../cblas.m4
	echo "HPMV(c, c, s)" | cat hpmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_chpmv_c_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chpmv_c_s.c && \
    mv tmp.BLAS_chpmv_c_s.c $@ && rm -f tmp.BLAS_chpmv_c_s.c*

$(SRC_PATH)/BLAS_zhpmv_z_d.c: hpmv.m4 hpmv-common.m4 ../cblas.m4
	echo "HPMV(z, z, d)" | cat hpmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhpmv_z_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhpmv_z_d.c && \
    mv tmp.BLAS_zhpmv_z_d.c $@ && rm -f tmp.BLAS_zhpmv_z_d.c*

$(SRC_PATH)/BLAS_chpmv_x.c: hpmv.m4 hpmv-common.m4 ../cblas.m4
	echo "HPMV(c, c, c, _x)" | cat hpmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_chpmv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chpmv_x.c && \
    mv tmp.BLAS_chpmv_x.c $@ && rm -f tmp.BLAS_chpmv_x.c*

$(SRC_PATH)/BLAS_zhpmv_x.c: hpmv.m4 hpmv-common.m4 ../cblas.m4
	echo "HPMV(z, z, z, _x)" | cat hpmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhpmv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhpmv_x.c && \
    mv tmp.BLAS_zhpmv_x.c $@ && rm -f tmp.BLAS_zhpmv_x.c*

$(SRC_PATH)/BLAS_zhpmv_z_c_x.c: hpmv.m4 hpmv-common.m4 ../cblas.m4
	echo "HPMV(z, z, c, _x)" | cat hpmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhpmv_z_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhpmv_z_c_x.c && \
    mv tmp.BLAS_zhpmv_z_c_x.c $@ && rm -f tmp.BLAS_zhpmv_z_c_x.c*

$(SRC_PATH)/BLAS_zhpmv_c_z_x.c: hpmv.m4 hpmv-common.m4 ../cblas.m4
	echo "HPMV(z, c, z, _x)" | cat hpmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhpmv_c_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhpmv_c_z_x.c && \
    mv tmp.BLAS_zhpmv_c_z_x.c $@ && rm -f tmp.BLAS_zhpmv_c_z_x.c*

$(SRC_PATH)/BLAS_zhpmv_c_c_x.c: hpmv.m4 hpmv-common.m4 ../cblas.m4
	echo "HPMV(z, c, c, _x)" | cat hpmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhpmv_c_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhpmv_c_c_x.c && \
    mv tmp.BLAS_zhpmv_c_c_x.c $@ && rm -f tmp.BLAS_zhpmv_c_c_x.c*

$(SRC_PATH)/BLAS_chpmv_c_s_x.c: hpmv.m4 hpmv-common.m4 ../cblas.m4
	echo "HPMV(c, c, s, _x)" | cat hpmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_chpmv_c_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chpmv_c_s_x.c && \
    mv tmp.BLAS_chpmv_c_s_x.c $@ && rm -f tmp.BLAS_chpmv_c_s_x.c*

$(SRC_PATH)/BLAS_zhpmv_z_d_x.c: hpmv.m4 hpmv-common.m4 ../cblas.m4
	echo "HPMV(z, z, d, _x)" | cat hpmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhpmv_z_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhpmv_z_d_x.c && \
    mv tmp.BLAS_zhpmv_z_d_x.c $@ && rm -f tmp.BLAS_zhpmv_z_d_x.c*



$(SRC_PATH)/BLAS_chpmv-f2c.c: ../f2c-bridge.m4 hpmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hpmv -DRARG="c, c, c" ../f2c-bridge.m4 > tmp.BLAS_chpmv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chpmv-f2c.c && \
    mv tmp.BLAS_chpmv-f2c.c $@ && rm -f tmp.BLAS_chpmv-f2c.c*

$(SRC_PATH)/BLAS_zhpmv-f2c.c: ../f2c-bridge.m4 hpmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hpmv -DRARG="z, z, z" ../f2c-bridge.m4 > tmp.BLAS_zhpmv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhpmv-f2c.c && \
    mv tmp.BLAS_zhpmv-f2c.c $@ && rm -f tmp.BLAS_zhpmv-f2c.c*

$(SRC_PATH)/BLAS_zhpmv_z_c-f2c.c: ../f2c-bridge.m4 hpmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hpmv -DRARG="z, z, c" ../f2c-bridge.m4 > tmp.BLAS_zhpmv_z_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhpmv_z_c-f2c.c && \
    mv tmp.BLAS_zhpmv_z_c-f2c.c $@ && rm -f tmp.BLAS_zhpmv_z_c-f2c.c*

$(SRC_PATH)/BLAS_zhpmv_c_z-f2c.c: ../f2c-bridge.m4 hpmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hpmv -DRARG="z, c, z" ../f2c-bridge.m4 > tmp.BLAS_zhpmv_c_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhpmv_c_z-f2c.c && \
    mv tmp.BLAS_zhpmv_c_z-f2c.c $@ && rm -f tmp.BLAS_zhpmv_c_z-f2c.c*

$(SRC_PATH)/BLAS_zhpmv_c_c-f2c.c: ../f2c-bridge.m4 hpmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hpmv -DRARG="z, c, c" ../f2c-bridge.m4 > tmp.BLAS_zhpmv_c_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhpmv_c_c-f2c.c && \
    mv tmp.BLAS_zhpmv_c_c-f2c.c $@ && rm -f tmp.BLAS_zhpmv_c_c-f2c.c*

$(SRC_PATH)/BLAS_chpmv_c_s-f2c.c: ../f2c-bridge.m4 hpmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hpmv -DRARG="c, c, s" ../f2c-bridge.m4 > tmp.BLAS_chpmv_c_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chpmv_c_s-f2c.c && \
    mv tmp.BLAS_chpmv_c_s-f2c.c $@ && rm -f tmp.BLAS_chpmv_c_s-f2c.c*

$(SRC_PATH)/BLAS_zhpmv_z_d-f2c.c: ../f2c-bridge.m4 hpmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hpmv -DRARG="z, z, d" ../f2c-bridge.m4 > tmp.BLAS_zhpmv_z_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhpmv_z_d-f2c.c && \
    mv tmp.BLAS_zhpmv_z_d-f2c.c $@ && rm -f tmp.BLAS_zhpmv_z_d-f2c.c*

$(SRC_PATH)/BLAS_chpmv_x-f2c.c: ../f2c-bridge.m4 hpmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hpmv -DRARG="c, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_chpmv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chpmv_x-f2c.c && \
    mv tmp.BLAS_chpmv_x-f2c.c $@ && rm -f tmp.BLAS_chpmv_x-f2c.c*

$(SRC_PATH)/BLAS_zhpmv_x-f2c.c: ../f2c-bridge.m4 hpmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hpmv -DRARG="z, z, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zhpmv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhpmv_x-f2c.c && \
    mv tmp.BLAS_zhpmv_x-f2c.c $@ && rm -f tmp.BLAS_zhpmv_x-f2c.c*

$(SRC_PATH)/BLAS_zhpmv_z_c_x-f2c.c: ../f2c-bridge.m4 hpmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hpmv -DRARG="z, z, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zhpmv_z_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhpmv_z_c_x-f2c.c && \
    mv tmp.BLAS_zhpmv_z_c_x-f2c.c $@ && rm -f tmp.BLAS_zhpmv_z_c_x-f2c.c*

$(SRC_PATH)/BLAS_zhpmv_c_z_x-f2c.c: ../f2c-bridge.m4 hpmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hpmv -DRARG="z, c, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zhpmv_c_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhpmv_c_z_x-f2c.c && \
    mv tmp.BLAS_zhpmv_c_z_x-f2c.c $@ && rm -f tmp.BLAS_zhpmv_c_z_x-f2c.c*

$(SRC_PATH)/BLAS_zhpmv_c_c_x-f2c.c: ../f2c-bridge.m4 hpmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hpmv -DRARG="z, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zhpmv_c_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhpmv_c_c_x-f2c.c && \
    mv tmp.BLAS_zhpmv_c_c_x-f2c.c $@ && rm -f tmp.BLAS_zhpmv_c_c_x-f2c.c*

$(SRC_PATH)/BLAS_chpmv_c_s_x-f2c.c: ../f2c-bridge.m4 hpmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hpmv -DRARG="c, c, s, _x" ../f2c-bridge.m4 > tmp.BLAS_chpmv_c_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chpmv_c_s_x-f2c.c && \
    mv tmp.BLAS_chpmv_c_s_x-f2c.c $@ && rm -f tmp.BLAS_chpmv_c_s_x-f2c.c*

$(SRC_PATH)/BLAS_zhpmv_z_d_x-f2c.c: ../f2c-bridge.m4 hpmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hpmv -DRARG="z, z, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zhpmv_z_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhpmv_z_d_x-f2c.c && \
    mv tmp.BLAS_zhpmv_z_d_x-f2c.c $@ && rm -f tmp.BLAS_zhpmv_z_d_x-f2c.c*


hpmv.h: ../cblas.m4 ../header.template hpmv-common.m4
	sed -e "s/routine/hpmv/g" -e "s/ROUTINE/HPMV/g" \
  ../header.template | $(M4) $(M4_OPTS) $(M4_PATH) - >$@

source-clean:
	rm -f $(HPMV_SRCS) hpmv.h

.PHONY: hpmv-source hpmv-header source-clean

