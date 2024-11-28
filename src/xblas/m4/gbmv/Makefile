include ../../make.conf
include ../../$(MAKEINC)


SRC_PATH = ../../src/gbmv
M4_PATH = -I. -I..

GBMV_SRCS = \
	$(SRC_PATH)/BLAS_dgbmv_d_s.c\
	$(SRC_PATH)/BLAS_dgbmv_s_d.c\
	$(SRC_PATH)/BLAS_dgbmv_s_s.c\
	$(SRC_PATH)/BLAS_zgbmv_z_c.c\
	$(SRC_PATH)/BLAS_zgbmv_c_z.c\
	$(SRC_PATH)/BLAS_zgbmv_c_c.c\
	$(SRC_PATH)/BLAS_cgbmv_c_s.c\
	$(SRC_PATH)/BLAS_cgbmv_s_c.c\
	$(SRC_PATH)/BLAS_cgbmv_s_s.c\
	$(SRC_PATH)/BLAS_zgbmv_z_d.c\
	$(SRC_PATH)/BLAS_zgbmv_d_z.c\
	$(SRC_PATH)/BLAS_zgbmv_d_d.c\
	$(SRC_PATH)/BLAS_sgbmv_x.c\
	$(SRC_PATH)/BLAS_dgbmv_x.c\
	$(SRC_PATH)/BLAS_cgbmv_x.c\
	$(SRC_PATH)/BLAS_zgbmv_x.c\
	$(SRC_PATH)/BLAS_dgbmv_d_s_x.c\
	$(SRC_PATH)/BLAS_dgbmv_s_d_x.c\
	$(SRC_PATH)/BLAS_dgbmv_s_s_x.c\
	$(SRC_PATH)/BLAS_zgbmv_z_c_x.c\
	$(SRC_PATH)/BLAS_zgbmv_c_z_x.c\
	$(SRC_PATH)/BLAS_zgbmv_c_c_x.c\
	$(SRC_PATH)/BLAS_cgbmv_c_s_x.c\
	$(SRC_PATH)/BLAS_cgbmv_s_c_x.c\
	$(SRC_PATH)/BLAS_cgbmv_s_s_x.c\
	$(SRC_PATH)/BLAS_zgbmv_z_d_x.c\
	$(SRC_PATH)/BLAS_zgbmv_d_z_x.c\
	$(SRC_PATH)/BLAS_zgbmv_d_d_x.c \
	\
	$(SRC_PATH)/BLAS_dgbmv_d_s-f2c.c\
	$(SRC_PATH)/BLAS_dgbmv_s_d-f2c.c\
	$(SRC_PATH)/BLAS_dgbmv_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv_c_c-f2c.c\
	$(SRC_PATH)/BLAS_cgbmv_c_s-f2c.c\
	$(SRC_PATH)/BLAS_cgbmv_s_c-f2c.c\
	$(SRC_PATH)/BLAS_cgbmv_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv_z_d-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv_d_z-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv_d_d-f2c.c\
	$(SRC_PATH)/BLAS_sgbmv_x-f2c.c\
	$(SRC_PATH)/BLAS_dgbmv_x-f2c.c\
	$(SRC_PATH)/BLAS_cgbmv_x-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv_x-f2c.c\
	$(SRC_PATH)/BLAS_dgbmv_d_s_x-f2c.c\
	$(SRC_PATH)/BLAS_dgbmv_s_d_x-f2c.c\
	$(SRC_PATH)/BLAS_dgbmv_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cgbmv_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_cgbmv_s_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cgbmv_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv_z_d_x-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv_d_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv_d_d_x-f2c.c

gbmv-source: $(GBMV_SRCS)

gbmv-header: gbmv.h

$(SRC_PATH)/BLAS_dgbmv_d_s.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(d, d, s)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgbmv_d_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv_d_s.c && \
    mv tmp.BLAS_dgbmv_d_s.c $@ && rm -f tmp.BLAS_dgbmv_d_s.c*

$(SRC_PATH)/BLAS_dgbmv_s_d.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(d, s, d)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgbmv_s_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv_s_d.c && \
    mv tmp.BLAS_dgbmv_s_d.c $@ && rm -f tmp.BLAS_dgbmv_s_d.c*

$(SRC_PATH)/BLAS_dgbmv_s_s.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(d, s, s)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgbmv_s_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv_s_s.c && \
    mv tmp.BLAS_dgbmv_s_s.c $@ && rm -f tmp.BLAS_dgbmv_s_s.c*

$(SRC_PATH)/BLAS_zgbmv_z_c.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(z, z, c)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgbmv_z_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv_z_c.c && \
    mv tmp.BLAS_zgbmv_z_c.c $@ && rm -f tmp.BLAS_zgbmv_z_c.c*

$(SRC_PATH)/BLAS_zgbmv_c_z.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(z, c, z)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgbmv_c_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv_c_z.c && \
    mv tmp.BLAS_zgbmv_c_z.c $@ && rm -f tmp.BLAS_zgbmv_c_z.c*

$(SRC_PATH)/BLAS_zgbmv_c_c.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(z, c, c)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgbmv_c_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv_c_c.c && \
    mv tmp.BLAS_zgbmv_c_c.c $@ && rm -f tmp.BLAS_zgbmv_c_c.c*

$(SRC_PATH)/BLAS_cgbmv_c_s.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(c, c, s)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgbmv_c_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv_c_s.c && \
    mv tmp.BLAS_cgbmv_c_s.c $@ && rm -f tmp.BLAS_cgbmv_c_s.c*

$(SRC_PATH)/BLAS_cgbmv_s_c.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(c, s, c)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgbmv_s_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv_s_c.c && \
    mv tmp.BLAS_cgbmv_s_c.c $@ && rm -f tmp.BLAS_cgbmv_s_c.c*

$(SRC_PATH)/BLAS_cgbmv_s_s.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(c, s, s)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgbmv_s_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv_s_s.c && \
    mv tmp.BLAS_cgbmv_s_s.c $@ && rm -f tmp.BLAS_cgbmv_s_s.c*

$(SRC_PATH)/BLAS_zgbmv_z_d.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(z, z, d)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgbmv_z_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv_z_d.c && \
    mv tmp.BLAS_zgbmv_z_d.c $@ && rm -f tmp.BLAS_zgbmv_z_d.c*

$(SRC_PATH)/BLAS_zgbmv_d_z.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(z, d, z)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgbmv_d_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv_d_z.c && \
    mv tmp.BLAS_zgbmv_d_z.c $@ && rm -f tmp.BLAS_zgbmv_d_z.c*

$(SRC_PATH)/BLAS_zgbmv_d_d.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(z, d, d)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgbmv_d_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv_d_d.c && \
    mv tmp.BLAS_zgbmv_d_d.c $@ && rm -f tmp.BLAS_zgbmv_d_d.c*

$(SRC_PATH)/BLAS_sgbmv_x.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(s, s, s, _x)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_sgbmv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sgbmv_x.c && \
    mv tmp.BLAS_sgbmv_x.c $@ && rm -f tmp.BLAS_sgbmv_x.c*

$(SRC_PATH)/BLAS_dgbmv_x.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(d, d, d, _x)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgbmv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv_x.c && \
    mv tmp.BLAS_dgbmv_x.c $@ && rm -f tmp.BLAS_dgbmv_x.c*

$(SRC_PATH)/BLAS_cgbmv_x.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(c, c, c, _x)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgbmv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv_x.c && \
    mv tmp.BLAS_cgbmv_x.c $@ && rm -f tmp.BLAS_cgbmv_x.c*

$(SRC_PATH)/BLAS_zgbmv_x.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(z, z, z, _x)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgbmv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv_x.c && \
    mv tmp.BLAS_zgbmv_x.c $@ && rm -f tmp.BLAS_zgbmv_x.c*

$(SRC_PATH)/BLAS_dgbmv_d_s_x.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(d, d, s, _x)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgbmv_d_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv_d_s_x.c && \
    mv tmp.BLAS_dgbmv_d_s_x.c $@ && rm -f tmp.BLAS_dgbmv_d_s_x.c*

$(SRC_PATH)/BLAS_dgbmv_s_d_x.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(d, s, d, _x)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgbmv_s_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv_s_d_x.c && \
    mv tmp.BLAS_dgbmv_s_d_x.c $@ && rm -f tmp.BLAS_dgbmv_s_d_x.c*

$(SRC_PATH)/BLAS_dgbmv_s_s_x.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(d, s, s, _x)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgbmv_s_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv_s_s_x.c && \
    mv tmp.BLAS_dgbmv_s_s_x.c $@ && rm -f tmp.BLAS_dgbmv_s_s_x.c*

$(SRC_PATH)/BLAS_zgbmv_z_c_x.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(z, z, c, _x)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgbmv_z_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv_z_c_x.c && \
    mv tmp.BLAS_zgbmv_z_c_x.c $@ && rm -f tmp.BLAS_zgbmv_z_c_x.c*

$(SRC_PATH)/BLAS_zgbmv_c_z_x.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(z, c, z, _x)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgbmv_c_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv_c_z_x.c && \
    mv tmp.BLAS_zgbmv_c_z_x.c $@ && rm -f tmp.BLAS_zgbmv_c_z_x.c*

$(SRC_PATH)/BLAS_zgbmv_c_c_x.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(z, c, c, _x)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgbmv_c_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv_c_c_x.c && \
    mv tmp.BLAS_zgbmv_c_c_x.c $@ && rm -f tmp.BLAS_zgbmv_c_c_x.c*

$(SRC_PATH)/BLAS_cgbmv_c_s_x.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(c, c, s, _x)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgbmv_c_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv_c_s_x.c && \
    mv tmp.BLAS_cgbmv_c_s_x.c $@ && rm -f tmp.BLAS_cgbmv_c_s_x.c*

$(SRC_PATH)/BLAS_cgbmv_s_c_x.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(c, s, c, _x)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgbmv_s_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv_s_c_x.c && \
    mv tmp.BLAS_cgbmv_s_c_x.c $@ && rm -f tmp.BLAS_cgbmv_s_c_x.c*

$(SRC_PATH)/BLAS_cgbmv_s_s_x.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(c, s, s, _x)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgbmv_s_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv_s_s_x.c && \
    mv tmp.BLAS_cgbmv_s_s_x.c $@ && rm -f tmp.BLAS_cgbmv_s_s_x.c*

$(SRC_PATH)/BLAS_zgbmv_z_d_x.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(z, z, d, _x)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgbmv_z_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv_z_d_x.c && \
    mv tmp.BLAS_zgbmv_z_d_x.c $@ && rm -f tmp.BLAS_zgbmv_z_d_x.c*

$(SRC_PATH)/BLAS_zgbmv_d_z_x.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(z, d, z, _x)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgbmv_d_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv_d_z_x.c && \
    mv tmp.BLAS_zgbmv_d_z_x.c $@ && rm -f tmp.BLAS_zgbmv_d_z_x.c*

$(SRC_PATH)/BLAS_zgbmv_d_d_x.c: gbmv.m4 gbmv-common.m4 ../cblas.m4
	echo "GBMV(z, d, d, _x)" | cat gbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgbmv_d_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv_d_d_x.c && \
    mv tmp.BLAS_zgbmv_d_d_x.c $@ && rm -f tmp.BLAS_zgbmv_d_d_x.c*



$(SRC_PATH)/BLAS_dgbmv_d_s-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="d, d, s" ../f2c-bridge.m4 > tmp.BLAS_dgbmv_d_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv_d_s-f2c.c && \
    mv tmp.BLAS_dgbmv_d_s-f2c.c $@ && rm -f tmp.BLAS_dgbmv_d_s-f2c.c*

$(SRC_PATH)/BLAS_dgbmv_s_d-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="d, s, d" ../f2c-bridge.m4 > tmp.BLAS_dgbmv_s_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv_s_d-f2c.c && \
    mv tmp.BLAS_dgbmv_s_d-f2c.c $@ && rm -f tmp.BLAS_dgbmv_s_d-f2c.c*

$(SRC_PATH)/BLAS_dgbmv_s_s-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="d, s, s" ../f2c-bridge.m4 > tmp.BLAS_dgbmv_s_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv_s_s-f2c.c && \
    mv tmp.BLAS_dgbmv_s_s-f2c.c $@ && rm -f tmp.BLAS_dgbmv_s_s-f2c.c*

$(SRC_PATH)/BLAS_zgbmv_z_c-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="z, z, c" ../f2c-bridge.m4 > tmp.BLAS_zgbmv_z_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv_z_c-f2c.c && \
    mv tmp.BLAS_zgbmv_z_c-f2c.c $@ && rm -f tmp.BLAS_zgbmv_z_c-f2c.c*

$(SRC_PATH)/BLAS_zgbmv_c_z-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="z, c, z" ../f2c-bridge.m4 > tmp.BLAS_zgbmv_c_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv_c_z-f2c.c && \
    mv tmp.BLAS_zgbmv_c_z-f2c.c $@ && rm -f tmp.BLAS_zgbmv_c_z-f2c.c*

$(SRC_PATH)/BLAS_zgbmv_c_c-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="z, c, c" ../f2c-bridge.m4 > tmp.BLAS_zgbmv_c_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv_c_c-f2c.c && \
    mv tmp.BLAS_zgbmv_c_c-f2c.c $@ && rm -f tmp.BLAS_zgbmv_c_c-f2c.c*

$(SRC_PATH)/BLAS_cgbmv_c_s-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="c, c, s" ../f2c-bridge.m4 > tmp.BLAS_cgbmv_c_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv_c_s-f2c.c && \
    mv tmp.BLAS_cgbmv_c_s-f2c.c $@ && rm -f tmp.BLAS_cgbmv_c_s-f2c.c*

$(SRC_PATH)/BLAS_cgbmv_s_c-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="c, s, c" ../f2c-bridge.m4 > tmp.BLAS_cgbmv_s_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv_s_c-f2c.c && \
    mv tmp.BLAS_cgbmv_s_c-f2c.c $@ && rm -f tmp.BLAS_cgbmv_s_c-f2c.c*

$(SRC_PATH)/BLAS_cgbmv_s_s-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="c, s, s" ../f2c-bridge.m4 > tmp.BLAS_cgbmv_s_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv_s_s-f2c.c && \
    mv tmp.BLAS_cgbmv_s_s-f2c.c $@ && rm -f tmp.BLAS_cgbmv_s_s-f2c.c*

$(SRC_PATH)/BLAS_zgbmv_z_d-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="z, z, d" ../f2c-bridge.m4 > tmp.BLAS_zgbmv_z_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv_z_d-f2c.c && \
    mv tmp.BLAS_zgbmv_z_d-f2c.c $@ && rm -f tmp.BLAS_zgbmv_z_d-f2c.c*

$(SRC_PATH)/BLAS_zgbmv_d_z-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="z, d, z" ../f2c-bridge.m4 > tmp.BLAS_zgbmv_d_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv_d_z-f2c.c && \
    mv tmp.BLAS_zgbmv_d_z-f2c.c $@ && rm -f tmp.BLAS_zgbmv_d_z-f2c.c*

$(SRC_PATH)/BLAS_zgbmv_d_d-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="z, d, d" ../f2c-bridge.m4 > tmp.BLAS_zgbmv_d_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv_d_d-f2c.c && \
    mv tmp.BLAS_zgbmv_d_d-f2c.c $@ && rm -f tmp.BLAS_zgbmv_d_d-f2c.c*

$(SRC_PATH)/BLAS_sgbmv_x-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="s, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_sgbmv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sgbmv_x-f2c.c && \
    mv tmp.BLAS_sgbmv_x-f2c.c $@ && rm -f tmp.BLAS_sgbmv_x-f2c.c*

$(SRC_PATH)/BLAS_dgbmv_x-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="d, d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dgbmv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv_x-f2c.c && \
    mv tmp.BLAS_dgbmv_x-f2c.c $@ && rm -f tmp.BLAS_dgbmv_x-f2c.c*

$(SRC_PATH)/BLAS_cgbmv_x-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="c, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_cgbmv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv_x-f2c.c && \
    mv tmp.BLAS_cgbmv_x-f2c.c $@ && rm -f tmp.BLAS_cgbmv_x-f2c.c*

$(SRC_PATH)/BLAS_zgbmv_x-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="z, z, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zgbmv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv_x-f2c.c && \
    mv tmp.BLAS_zgbmv_x-f2c.c $@ && rm -f tmp.BLAS_zgbmv_x-f2c.c*

$(SRC_PATH)/BLAS_dgbmv_d_s_x-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="d, d, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dgbmv_d_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv_d_s_x-f2c.c && \
    mv tmp.BLAS_dgbmv_d_s_x-f2c.c $@ && rm -f tmp.BLAS_dgbmv_d_s_x-f2c.c*

$(SRC_PATH)/BLAS_dgbmv_s_d_x-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="d, s, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dgbmv_s_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv_s_d_x-f2c.c && \
    mv tmp.BLAS_dgbmv_s_d_x-f2c.c $@ && rm -f tmp.BLAS_dgbmv_s_d_x-f2c.c*

$(SRC_PATH)/BLAS_dgbmv_s_s_x-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="d, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dgbmv_s_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv_s_s_x-f2c.c && \
    mv tmp.BLAS_dgbmv_s_s_x-f2c.c $@ && rm -f tmp.BLAS_dgbmv_s_s_x-f2c.c*

$(SRC_PATH)/BLAS_zgbmv_z_c_x-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="z, z, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zgbmv_z_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv_z_c_x-f2c.c && \
    mv tmp.BLAS_zgbmv_z_c_x-f2c.c $@ && rm -f tmp.BLAS_zgbmv_z_c_x-f2c.c*

$(SRC_PATH)/BLAS_zgbmv_c_z_x-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="z, c, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zgbmv_c_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv_c_z_x-f2c.c && \
    mv tmp.BLAS_zgbmv_c_z_x-f2c.c $@ && rm -f tmp.BLAS_zgbmv_c_z_x-f2c.c*

$(SRC_PATH)/BLAS_zgbmv_c_c_x-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="z, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zgbmv_c_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv_c_c_x-f2c.c && \
    mv tmp.BLAS_zgbmv_c_c_x-f2c.c $@ && rm -f tmp.BLAS_zgbmv_c_c_x-f2c.c*

$(SRC_PATH)/BLAS_cgbmv_c_s_x-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="c, c, s, _x" ../f2c-bridge.m4 > tmp.BLAS_cgbmv_c_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv_c_s_x-f2c.c && \
    mv tmp.BLAS_cgbmv_c_s_x-f2c.c $@ && rm -f tmp.BLAS_cgbmv_c_s_x-f2c.c*

$(SRC_PATH)/BLAS_cgbmv_s_c_x-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="c, s, c, _x" ../f2c-bridge.m4 > tmp.BLAS_cgbmv_s_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv_s_c_x-f2c.c && \
    mv tmp.BLAS_cgbmv_s_c_x-f2c.c $@ && rm -f tmp.BLAS_cgbmv_s_c_x-f2c.c*

$(SRC_PATH)/BLAS_cgbmv_s_s_x-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="c, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_cgbmv_s_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv_s_s_x-f2c.c && \
    mv tmp.BLAS_cgbmv_s_s_x-f2c.c $@ && rm -f tmp.BLAS_cgbmv_s_s_x-f2c.c*

$(SRC_PATH)/BLAS_zgbmv_z_d_x-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="z, z, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zgbmv_z_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv_z_d_x-f2c.c && \
    mv tmp.BLAS_zgbmv_z_d_x-f2c.c $@ && rm -f tmp.BLAS_zgbmv_z_d_x-f2c.c*

$(SRC_PATH)/BLAS_zgbmv_d_z_x-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="z, d, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zgbmv_d_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv_d_z_x-f2c.c && \
    mv tmp.BLAS_zgbmv_d_z_x-f2c.c $@ && rm -f tmp.BLAS_zgbmv_d_z_x-f2c.c*

$(SRC_PATH)/BLAS_zgbmv_d_d_x-f2c.c: ../f2c-bridge.m4 gbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv -DRARG="z, d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zgbmv_d_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv_d_d_x-f2c.c && \
    mv tmp.BLAS_zgbmv_d_d_x-f2c.c $@ && rm -f tmp.BLAS_zgbmv_d_d_x-f2c.c*


gbmv.h: ../cblas.m4 ../header.template gbmv-common.m4
	sed -e "s/routine/gbmv/g" -e "s/ROUTINE/GBMV/g" \
  ../header.template | $(M4) $(M4_OPTS) $(M4_PATH) - >$@

source-clean:
	rm -f $(GBMV_SRCS) gbmv.h

.PHONY: gbmv-source gbmv-header source-clean

