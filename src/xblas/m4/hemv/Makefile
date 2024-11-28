include ../../make.conf
include ../../$(MAKEINC)


SRC_PATH = ../../src/hemv
M4_PATH = -I. -I..

HEMV_SRCS = \
	$(SRC_PATH)/BLAS_chemv.c\
	$(SRC_PATH)/BLAS_zhemv.c\
	$(SRC_PATH)/BLAS_zhemv_z_c.c\
	$(SRC_PATH)/BLAS_zhemv_c_z.c\
	$(SRC_PATH)/BLAS_zhemv_c_c.c\
	$(SRC_PATH)/BLAS_chemv_c_s.c\
	$(SRC_PATH)/BLAS_zhemv_z_d.c\
	$(SRC_PATH)/BLAS_chemv_x.c\
	$(SRC_PATH)/BLAS_zhemv_x.c\
	$(SRC_PATH)/BLAS_zhemv_z_c_x.c\
	$(SRC_PATH)/BLAS_zhemv_c_z_x.c\
	$(SRC_PATH)/BLAS_zhemv_c_c_x.c\
	$(SRC_PATH)/BLAS_chemv_c_s_x.c\
	$(SRC_PATH)/BLAS_zhemv_z_d_x.c \
	\
	$(SRC_PATH)/BLAS_chemv-f2c.c\
	$(SRC_PATH)/BLAS_zhemv-f2c.c\
	$(SRC_PATH)/BLAS_zhemv_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zhemv_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zhemv_c_c-f2c.c\
	$(SRC_PATH)/BLAS_chemv_c_s-f2c.c\
	$(SRC_PATH)/BLAS_zhemv_z_d-f2c.c\
	$(SRC_PATH)/BLAS_chemv_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemv_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemv_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemv_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemv_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_chemv_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemv_z_d_x-f2c.c

hemv-source: $(HEMV_SRCS)

hemv-header: hemv.h

$(SRC_PATH)/BLAS_chemv.c: hemv.m4 hemv-common.m4 ../cblas.m4
	echo "HEMV(c, c, c)" | cat hemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_chemv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chemv.c && \
    mv tmp.BLAS_chemv.c $@ && rm -f tmp.BLAS_chemv.c*

$(SRC_PATH)/BLAS_zhemv.c: hemv.m4 hemv-common.m4 ../cblas.m4
	echo "HEMV(z, z, z)" | cat hemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv.c && \
    mv tmp.BLAS_zhemv.c $@ && rm -f tmp.BLAS_zhemv.c*

$(SRC_PATH)/BLAS_zhemv_z_c.c: hemv.m4 hemv-common.m4 ../cblas.m4
	echo "HEMV(z, z, c)" | cat hemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemv_z_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv_z_c.c && \
    mv tmp.BLAS_zhemv_z_c.c $@ && rm -f tmp.BLAS_zhemv_z_c.c*

$(SRC_PATH)/BLAS_zhemv_c_z.c: hemv.m4 hemv-common.m4 ../cblas.m4
	echo "HEMV(z, c, z)" | cat hemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemv_c_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv_c_z.c && \
    mv tmp.BLAS_zhemv_c_z.c $@ && rm -f tmp.BLAS_zhemv_c_z.c*

$(SRC_PATH)/BLAS_zhemv_c_c.c: hemv.m4 hemv-common.m4 ../cblas.m4
	echo "HEMV(z, c, c)" | cat hemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemv_c_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv_c_c.c && \
    mv tmp.BLAS_zhemv_c_c.c $@ && rm -f tmp.BLAS_zhemv_c_c.c*

$(SRC_PATH)/BLAS_chemv_c_s.c: hemv.m4 hemv-common.m4 ../cblas.m4
	echo "HEMV(c, c, s)" | cat hemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_chemv_c_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chemv_c_s.c && \
    mv tmp.BLAS_chemv_c_s.c $@ && rm -f tmp.BLAS_chemv_c_s.c*

$(SRC_PATH)/BLAS_zhemv_z_d.c: hemv.m4 hemv-common.m4 ../cblas.m4
	echo "HEMV(z, z, d)" | cat hemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemv_z_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv_z_d.c && \
    mv tmp.BLAS_zhemv_z_d.c $@ && rm -f tmp.BLAS_zhemv_z_d.c*

$(SRC_PATH)/BLAS_chemv_x.c: hemv.m4 hemv-common.m4 ../cblas.m4
	echo "HEMV(c, c, c, _x)" | cat hemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_chemv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chemv_x.c && \
    mv tmp.BLAS_chemv_x.c $@ && rm -f tmp.BLAS_chemv_x.c*

$(SRC_PATH)/BLAS_zhemv_x.c: hemv.m4 hemv-common.m4 ../cblas.m4
	echo "HEMV(z, z, z, _x)" | cat hemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv_x.c && \
    mv tmp.BLAS_zhemv_x.c $@ && rm -f tmp.BLAS_zhemv_x.c*

$(SRC_PATH)/BLAS_zhemv_z_c_x.c: hemv.m4 hemv-common.m4 ../cblas.m4
	echo "HEMV(z, z, c, _x)" | cat hemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemv_z_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv_z_c_x.c && \
    mv tmp.BLAS_zhemv_z_c_x.c $@ && rm -f tmp.BLAS_zhemv_z_c_x.c*

$(SRC_PATH)/BLAS_zhemv_c_z_x.c: hemv.m4 hemv-common.m4 ../cblas.m4
	echo "HEMV(z, c, z, _x)" | cat hemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemv_c_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv_c_z_x.c && \
    mv tmp.BLAS_zhemv_c_z_x.c $@ && rm -f tmp.BLAS_zhemv_c_z_x.c*

$(SRC_PATH)/BLAS_zhemv_c_c_x.c: hemv.m4 hemv-common.m4 ../cblas.m4
	echo "HEMV(z, c, c, _x)" | cat hemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemv_c_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv_c_c_x.c && \
    mv tmp.BLAS_zhemv_c_c_x.c $@ && rm -f tmp.BLAS_zhemv_c_c_x.c*

$(SRC_PATH)/BLAS_chemv_c_s_x.c: hemv.m4 hemv-common.m4 ../cblas.m4
	echo "HEMV(c, c, s, _x)" | cat hemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_chemv_c_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chemv_c_s_x.c && \
    mv tmp.BLAS_chemv_c_s_x.c $@ && rm -f tmp.BLAS_chemv_c_s_x.c*

$(SRC_PATH)/BLAS_zhemv_z_d_x.c: hemv.m4 hemv-common.m4 ../cblas.m4
	echo "HEMV(z, z, d, _x)" | cat hemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemv_z_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv_z_d_x.c && \
    mv tmp.BLAS_zhemv_z_d_x.c $@ && rm -f tmp.BLAS_zhemv_z_d_x.c*



$(SRC_PATH)/BLAS_chemv-f2c.c: ../f2c-bridge.m4 hemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv -DRARG="c, c, c" ../f2c-bridge.m4 > tmp.BLAS_chemv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chemv-f2c.c && \
    mv tmp.BLAS_chemv-f2c.c $@ && rm -f tmp.BLAS_chemv-f2c.c*

$(SRC_PATH)/BLAS_zhemv-f2c.c: ../f2c-bridge.m4 hemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv -DRARG="z, z, z" ../f2c-bridge.m4 > tmp.BLAS_zhemv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv-f2c.c && \
    mv tmp.BLAS_zhemv-f2c.c $@ && rm -f tmp.BLAS_zhemv-f2c.c*

$(SRC_PATH)/BLAS_zhemv_z_c-f2c.c: ../f2c-bridge.m4 hemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv -DRARG="z, z, c" ../f2c-bridge.m4 > tmp.BLAS_zhemv_z_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv_z_c-f2c.c && \
    mv tmp.BLAS_zhemv_z_c-f2c.c $@ && rm -f tmp.BLAS_zhemv_z_c-f2c.c*

$(SRC_PATH)/BLAS_zhemv_c_z-f2c.c: ../f2c-bridge.m4 hemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv -DRARG="z, c, z" ../f2c-bridge.m4 > tmp.BLAS_zhemv_c_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv_c_z-f2c.c && \
    mv tmp.BLAS_zhemv_c_z-f2c.c $@ && rm -f tmp.BLAS_zhemv_c_z-f2c.c*

$(SRC_PATH)/BLAS_zhemv_c_c-f2c.c: ../f2c-bridge.m4 hemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv -DRARG="z, c, c" ../f2c-bridge.m4 > tmp.BLAS_zhemv_c_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv_c_c-f2c.c && \
    mv tmp.BLAS_zhemv_c_c-f2c.c $@ && rm -f tmp.BLAS_zhemv_c_c-f2c.c*

$(SRC_PATH)/BLAS_chemv_c_s-f2c.c: ../f2c-bridge.m4 hemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv -DRARG="c, c, s" ../f2c-bridge.m4 > tmp.BLAS_chemv_c_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chemv_c_s-f2c.c && \
    mv tmp.BLAS_chemv_c_s-f2c.c $@ && rm -f tmp.BLAS_chemv_c_s-f2c.c*

$(SRC_PATH)/BLAS_zhemv_z_d-f2c.c: ../f2c-bridge.m4 hemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv -DRARG="z, z, d" ../f2c-bridge.m4 > tmp.BLAS_zhemv_z_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv_z_d-f2c.c && \
    mv tmp.BLAS_zhemv_z_d-f2c.c $@ && rm -f tmp.BLAS_zhemv_z_d-f2c.c*

$(SRC_PATH)/BLAS_chemv_x-f2c.c: ../f2c-bridge.m4 hemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv -DRARG="c, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_chemv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chemv_x-f2c.c && \
    mv tmp.BLAS_chemv_x-f2c.c $@ && rm -f tmp.BLAS_chemv_x-f2c.c*

$(SRC_PATH)/BLAS_zhemv_x-f2c.c: ../f2c-bridge.m4 hemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv -DRARG="z, z, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zhemv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv_x-f2c.c && \
    mv tmp.BLAS_zhemv_x-f2c.c $@ && rm -f tmp.BLAS_zhemv_x-f2c.c*

$(SRC_PATH)/BLAS_zhemv_z_c_x-f2c.c: ../f2c-bridge.m4 hemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv -DRARG="z, z, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zhemv_z_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv_z_c_x-f2c.c && \
    mv tmp.BLAS_zhemv_z_c_x-f2c.c $@ && rm -f tmp.BLAS_zhemv_z_c_x-f2c.c*

$(SRC_PATH)/BLAS_zhemv_c_z_x-f2c.c: ../f2c-bridge.m4 hemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv -DRARG="z, c, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zhemv_c_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv_c_z_x-f2c.c && \
    mv tmp.BLAS_zhemv_c_z_x-f2c.c $@ && rm -f tmp.BLAS_zhemv_c_z_x-f2c.c*

$(SRC_PATH)/BLAS_zhemv_c_c_x-f2c.c: ../f2c-bridge.m4 hemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv -DRARG="z, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zhemv_c_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv_c_c_x-f2c.c && \
    mv tmp.BLAS_zhemv_c_c_x-f2c.c $@ && rm -f tmp.BLAS_zhemv_c_c_x-f2c.c*

$(SRC_PATH)/BLAS_chemv_c_s_x-f2c.c: ../f2c-bridge.m4 hemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv -DRARG="c, c, s, _x" ../f2c-bridge.m4 > tmp.BLAS_chemv_c_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chemv_c_s_x-f2c.c && \
    mv tmp.BLAS_chemv_c_s_x-f2c.c $@ && rm -f tmp.BLAS_chemv_c_s_x-f2c.c*

$(SRC_PATH)/BLAS_zhemv_z_d_x-f2c.c: ../f2c-bridge.m4 hemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv -DRARG="z, z, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zhemv_z_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv_z_d_x-f2c.c && \
    mv tmp.BLAS_zhemv_z_d_x-f2c.c $@ && rm -f tmp.BLAS_zhemv_z_d_x-f2c.c*


hemv.h: ../cblas.m4 ../header.template hemv-common.m4
	sed -e "s/routine/hemv/g" -e "s/ROUTINE/HEMV/g" \
  ../header.template | $(M4) $(M4_OPTS) $(M4_PATH) - >$@

source-clean:
	rm -f $(HEMV_SRCS) hemv.h

.PHONY: hemv-source hemv-header source-clean

