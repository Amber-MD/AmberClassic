include ../../make.conf
include ../../$(MAKEINC)


SRC_PATH = ../../src/hemv2
M4_PATH = -I. -I..

HEMV2_SRCS = \
	$(SRC_PATH)/BLAS_chemv2.c\
	$(SRC_PATH)/BLAS_zhemv2.c\
	$(SRC_PATH)/BLAS_zhemv2_z_c.c\
	$(SRC_PATH)/BLAS_zhemv2_c_z.c\
	$(SRC_PATH)/BLAS_zhemv2_c_c.c\
	$(SRC_PATH)/BLAS_chemv2_c_s.c\
	$(SRC_PATH)/BLAS_zhemv2_z_d.c\
	$(SRC_PATH)/BLAS_chemv2_x.c\
	$(SRC_PATH)/BLAS_zhemv2_x.c\
	$(SRC_PATH)/BLAS_zhemv2_z_c_x.c\
	$(SRC_PATH)/BLAS_zhemv2_c_z_x.c\
	$(SRC_PATH)/BLAS_zhemv2_c_c_x.c\
	$(SRC_PATH)/BLAS_chemv2_c_s_x.c\
	$(SRC_PATH)/BLAS_zhemv2_z_d_x.c \
	\
	$(SRC_PATH)/BLAS_chemv2-f2c.c\
	$(SRC_PATH)/BLAS_zhemv2-f2c.c\
	$(SRC_PATH)/BLAS_zhemv2_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zhemv2_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zhemv2_c_c-f2c.c\
	$(SRC_PATH)/BLAS_chemv2_c_s-f2c.c\
	$(SRC_PATH)/BLAS_zhemv2_z_d-f2c.c\
	$(SRC_PATH)/BLAS_chemv2_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemv2_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemv2_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemv2_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemv2_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_chemv2_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemv2_z_d_x-f2c.c

hemv2-source: $(HEMV2_SRCS)

hemv2-header: hemv2.h

$(SRC_PATH)/BLAS_chemv2.c: hemv2.m4 hemv2-common.m4 ../cblas.m4
	echo "HEMV2(c, c, c)" | cat hemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_chemv2.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chemv2.c && \
    mv tmp.BLAS_chemv2.c $@ && rm -f tmp.BLAS_chemv2.c*

$(SRC_PATH)/BLAS_zhemv2.c: hemv2.m4 hemv2-common.m4 ../cblas.m4
	echo "HEMV2(z, z, z)" | cat hemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemv2.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv2.c && \
    mv tmp.BLAS_zhemv2.c $@ && rm -f tmp.BLAS_zhemv2.c*

$(SRC_PATH)/BLAS_zhemv2_z_c.c: hemv2.m4 hemv2-common.m4 ../cblas.m4
	echo "HEMV2(z, z, c)" | cat hemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemv2_z_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv2_z_c.c && \
    mv tmp.BLAS_zhemv2_z_c.c $@ && rm -f tmp.BLAS_zhemv2_z_c.c*

$(SRC_PATH)/BLAS_zhemv2_c_z.c: hemv2.m4 hemv2-common.m4 ../cblas.m4
	echo "HEMV2(z, c, z)" | cat hemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemv2_c_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv2_c_z.c && \
    mv tmp.BLAS_zhemv2_c_z.c $@ && rm -f tmp.BLAS_zhemv2_c_z.c*

$(SRC_PATH)/BLAS_zhemv2_c_c.c: hemv2.m4 hemv2-common.m4 ../cblas.m4
	echo "HEMV2(z, c, c)" | cat hemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemv2_c_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv2_c_c.c && \
    mv tmp.BLAS_zhemv2_c_c.c $@ && rm -f tmp.BLAS_zhemv2_c_c.c*

$(SRC_PATH)/BLAS_chemv2_c_s.c: hemv2.m4 hemv2-common.m4 ../cblas.m4
	echo "HEMV2(c, c, s)" | cat hemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_chemv2_c_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chemv2_c_s.c && \
    mv tmp.BLAS_chemv2_c_s.c $@ && rm -f tmp.BLAS_chemv2_c_s.c*

$(SRC_PATH)/BLAS_zhemv2_z_d.c: hemv2.m4 hemv2-common.m4 ../cblas.m4
	echo "HEMV2(z, z, d)" | cat hemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemv2_z_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv2_z_d.c && \
    mv tmp.BLAS_zhemv2_z_d.c $@ && rm -f tmp.BLAS_zhemv2_z_d.c*

$(SRC_PATH)/BLAS_chemv2_x.c: hemv2.m4 hemv2-common.m4 ../cblas.m4
	echo "HEMV2(c, c, c, _x)" | cat hemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_chemv2_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chemv2_x.c && \
    mv tmp.BLAS_chemv2_x.c $@ && rm -f tmp.BLAS_chemv2_x.c*

$(SRC_PATH)/BLAS_zhemv2_x.c: hemv2.m4 hemv2-common.m4 ../cblas.m4
	echo "HEMV2(z, z, z, _x)" | cat hemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemv2_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv2_x.c && \
    mv tmp.BLAS_zhemv2_x.c $@ && rm -f tmp.BLAS_zhemv2_x.c*

$(SRC_PATH)/BLAS_zhemv2_z_c_x.c: hemv2.m4 hemv2-common.m4 ../cblas.m4
	echo "HEMV2(z, z, c, _x)" | cat hemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemv2_z_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv2_z_c_x.c && \
    mv tmp.BLAS_zhemv2_z_c_x.c $@ && rm -f tmp.BLAS_zhemv2_z_c_x.c*

$(SRC_PATH)/BLAS_zhemv2_c_z_x.c: hemv2.m4 hemv2-common.m4 ../cblas.m4
	echo "HEMV2(z, c, z, _x)" | cat hemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemv2_c_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv2_c_z_x.c && \
    mv tmp.BLAS_zhemv2_c_z_x.c $@ && rm -f tmp.BLAS_zhemv2_c_z_x.c*

$(SRC_PATH)/BLAS_zhemv2_c_c_x.c: hemv2.m4 hemv2-common.m4 ../cblas.m4
	echo "HEMV2(z, c, c, _x)" | cat hemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemv2_c_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv2_c_c_x.c && \
    mv tmp.BLAS_zhemv2_c_c_x.c $@ && rm -f tmp.BLAS_zhemv2_c_c_x.c*

$(SRC_PATH)/BLAS_chemv2_c_s_x.c: hemv2.m4 hemv2-common.m4 ../cblas.m4
	echo "HEMV2(c, c, s, _x)" | cat hemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_chemv2_c_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chemv2_c_s_x.c && \
    mv tmp.BLAS_chemv2_c_s_x.c $@ && rm -f tmp.BLAS_chemv2_c_s_x.c*

$(SRC_PATH)/BLAS_zhemv2_z_d_x.c: hemv2.m4 hemv2-common.m4 ../cblas.m4
	echo "HEMV2(z, z, d, _x)" | cat hemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemv2_z_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv2_z_d_x.c && \
    mv tmp.BLAS_zhemv2_z_d_x.c $@ && rm -f tmp.BLAS_zhemv2_z_d_x.c*



$(SRC_PATH)/BLAS_chemv2-f2c.c: ../f2c-bridge.m4 hemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv2 -DRARG="c, c, c" ../f2c-bridge.m4 > tmp.BLAS_chemv2-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chemv2-f2c.c && \
    mv tmp.BLAS_chemv2-f2c.c $@ && rm -f tmp.BLAS_chemv2-f2c.c*

$(SRC_PATH)/BLAS_zhemv2-f2c.c: ../f2c-bridge.m4 hemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv2 -DRARG="z, z, z" ../f2c-bridge.m4 > tmp.BLAS_zhemv2-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv2-f2c.c && \
    mv tmp.BLAS_zhemv2-f2c.c $@ && rm -f tmp.BLAS_zhemv2-f2c.c*

$(SRC_PATH)/BLAS_zhemv2_z_c-f2c.c: ../f2c-bridge.m4 hemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv2 -DRARG="z, z, c" ../f2c-bridge.m4 > tmp.BLAS_zhemv2_z_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv2_z_c-f2c.c && \
    mv tmp.BLAS_zhemv2_z_c-f2c.c $@ && rm -f tmp.BLAS_zhemv2_z_c-f2c.c*

$(SRC_PATH)/BLAS_zhemv2_c_z-f2c.c: ../f2c-bridge.m4 hemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv2 -DRARG="z, c, z" ../f2c-bridge.m4 > tmp.BLAS_zhemv2_c_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv2_c_z-f2c.c && \
    mv tmp.BLAS_zhemv2_c_z-f2c.c $@ && rm -f tmp.BLAS_zhemv2_c_z-f2c.c*

$(SRC_PATH)/BLAS_zhemv2_c_c-f2c.c: ../f2c-bridge.m4 hemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv2 -DRARG="z, c, c" ../f2c-bridge.m4 > tmp.BLAS_zhemv2_c_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv2_c_c-f2c.c && \
    mv tmp.BLAS_zhemv2_c_c-f2c.c $@ && rm -f tmp.BLAS_zhemv2_c_c-f2c.c*

$(SRC_PATH)/BLAS_chemv2_c_s-f2c.c: ../f2c-bridge.m4 hemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv2 -DRARG="c, c, s" ../f2c-bridge.m4 > tmp.BLAS_chemv2_c_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chemv2_c_s-f2c.c && \
    mv tmp.BLAS_chemv2_c_s-f2c.c $@ && rm -f tmp.BLAS_chemv2_c_s-f2c.c*

$(SRC_PATH)/BLAS_zhemv2_z_d-f2c.c: ../f2c-bridge.m4 hemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv2 -DRARG="z, z, d" ../f2c-bridge.m4 > tmp.BLAS_zhemv2_z_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv2_z_d-f2c.c && \
    mv tmp.BLAS_zhemv2_z_d-f2c.c $@ && rm -f tmp.BLAS_zhemv2_z_d-f2c.c*

$(SRC_PATH)/BLAS_chemv2_x-f2c.c: ../f2c-bridge.m4 hemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv2 -DRARG="c, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_chemv2_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chemv2_x-f2c.c && \
    mv tmp.BLAS_chemv2_x-f2c.c $@ && rm -f tmp.BLAS_chemv2_x-f2c.c*

$(SRC_PATH)/BLAS_zhemv2_x-f2c.c: ../f2c-bridge.m4 hemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv2 -DRARG="z, z, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zhemv2_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv2_x-f2c.c && \
    mv tmp.BLAS_zhemv2_x-f2c.c $@ && rm -f tmp.BLAS_zhemv2_x-f2c.c*

$(SRC_PATH)/BLAS_zhemv2_z_c_x-f2c.c: ../f2c-bridge.m4 hemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv2 -DRARG="z, z, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zhemv2_z_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv2_z_c_x-f2c.c && \
    mv tmp.BLAS_zhemv2_z_c_x-f2c.c $@ && rm -f tmp.BLAS_zhemv2_z_c_x-f2c.c*

$(SRC_PATH)/BLAS_zhemv2_c_z_x-f2c.c: ../f2c-bridge.m4 hemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv2 -DRARG="z, c, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zhemv2_c_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv2_c_z_x-f2c.c && \
    mv tmp.BLAS_zhemv2_c_z_x-f2c.c $@ && rm -f tmp.BLAS_zhemv2_c_z_x-f2c.c*

$(SRC_PATH)/BLAS_zhemv2_c_c_x-f2c.c: ../f2c-bridge.m4 hemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv2 -DRARG="z, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zhemv2_c_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv2_c_c_x-f2c.c && \
    mv tmp.BLAS_zhemv2_c_c_x-f2c.c $@ && rm -f tmp.BLAS_zhemv2_c_c_x-f2c.c*

$(SRC_PATH)/BLAS_chemv2_c_s_x-f2c.c: ../f2c-bridge.m4 hemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv2 -DRARG="c, c, s, _x" ../f2c-bridge.m4 > tmp.BLAS_chemv2_c_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chemv2_c_s_x-f2c.c && \
    mv tmp.BLAS_chemv2_c_s_x-f2c.c $@ && rm -f tmp.BLAS_chemv2_c_s_x-f2c.c*

$(SRC_PATH)/BLAS_zhemv2_z_d_x-f2c.c: ../f2c-bridge.m4 hemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemv2 -DRARG="z, z, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zhemv2_z_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemv2_z_d_x-f2c.c && \
    mv tmp.BLAS_zhemv2_z_d_x-f2c.c $@ && rm -f tmp.BLAS_zhemv2_z_d_x-f2c.c*


hemv2.h: ../cblas.m4 ../header.template hemv2-common.m4
	sed -e "s/routine/hemv2/g" -e "s/ROUTINE/HEMV2/g" \
  ../header.template | $(M4) $(M4_OPTS) $(M4_PATH) - >$@

source-clean:
	rm -f $(HEMV2_SRCS) hemv2.h

.PHONY: hemv2-source hemv2-header source-clean

