include ../../make.conf
include ../../$(MAKEINC)


SRC_PATH = ../../src/gbmv2
M4_PATH = -I. -I..

GBMV2_SRCS = \
	$(SRC_PATH)/BLAS_dgbmv2_d_s.c\
	$(SRC_PATH)/BLAS_dgbmv2_s_d.c\
	$(SRC_PATH)/BLAS_dgbmv2_s_s.c\
	$(SRC_PATH)/BLAS_zgbmv2_z_c.c\
	$(SRC_PATH)/BLAS_zgbmv2_c_z.c\
	$(SRC_PATH)/BLAS_zgbmv2_c_c.c\
	$(SRC_PATH)/BLAS_cgbmv2_c_s.c\
	$(SRC_PATH)/BLAS_cgbmv2_s_c.c\
	$(SRC_PATH)/BLAS_cgbmv2_s_s.c\
	$(SRC_PATH)/BLAS_zgbmv2_z_d.c\
	$(SRC_PATH)/BLAS_zgbmv2_d_z.c\
	$(SRC_PATH)/BLAS_zgbmv2_d_d.c\
	$(SRC_PATH)/BLAS_sgbmv2_x.c\
	$(SRC_PATH)/BLAS_dgbmv2_x.c\
	$(SRC_PATH)/BLAS_cgbmv2_x.c\
	$(SRC_PATH)/BLAS_zgbmv2_x.c\
	$(SRC_PATH)/BLAS_dgbmv2_d_s_x.c\
	$(SRC_PATH)/BLAS_dgbmv2_s_d_x.c\
	$(SRC_PATH)/BLAS_dgbmv2_s_s_x.c\
	$(SRC_PATH)/BLAS_zgbmv2_z_c_x.c\
	$(SRC_PATH)/BLAS_zgbmv2_c_z_x.c\
	$(SRC_PATH)/BLAS_zgbmv2_c_c_x.c\
	$(SRC_PATH)/BLAS_cgbmv2_c_s_x.c\
	$(SRC_PATH)/BLAS_cgbmv2_s_c_x.c\
	$(SRC_PATH)/BLAS_cgbmv2_s_s_x.c\
	$(SRC_PATH)/BLAS_zgbmv2_z_d_x.c\
	$(SRC_PATH)/BLAS_zgbmv2_d_z_x.c\
	$(SRC_PATH)/BLAS_zgbmv2_d_d_x.c \
	\
	$(SRC_PATH)/BLAS_dgbmv2_d_s-f2c.c\
	$(SRC_PATH)/BLAS_dgbmv2_s_d-f2c.c\
	$(SRC_PATH)/BLAS_dgbmv2_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv2_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv2_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv2_c_c-f2c.c\
	$(SRC_PATH)/BLAS_cgbmv2_c_s-f2c.c\
	$(SRC_PATH)/BLAS_cgbmv2_s_c-f2c.c\
	$(SRC_PATH)/BLAS_cgbmv2_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv2_z_d-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv2_d_z-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv2_d_d-f2c.c\
	$(SRC_PATH)/BLAS_sgbmv2_x-f2c.c\
	$(SRC_PATH)/BLAS_dgbmv2_x-f2c.c\
	$(SRC_PATH)/BLAS_cgbmv2_x-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv2_x-f2c.c\
	$(SRC_PATH)/BLAS_dgbmv2_d_s_x-f2c.c\
	$(SRC_PATH)/BLAS_dgbmv2_s_d_x-f2c.c\
	$(SRC_PATH)/BLAS_dgbmv2_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv2_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv2_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv2_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cgbmv2_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_cgbmv2_s_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cgbmv2_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv2_z_d_x-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv2_d_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zgbmv2_d_d_x-f2c.c

gbmv2-source: $(GBMV2_SRCS)

gbmv2-header: gbmv2.h

$(SRC_PATH)/BLAS_dgbmv2_d_s.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(d, d, s)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgbmv2_d_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv2_d_s.c && \
    mv tmp.BLAS_dgbmv2_d_s.c $@ && rm -f tmp.BLAS_dgbmv2_d_s.c*

$(SRC_PATH)/BLAS_dgbmv2_s_d.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(d, s, d)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgbmv2_s_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv2_s_d.c && \
    mv tmp.BLAS_dgbmv2_s_d.c $@ && rm -f tmp.BLAS_dgbmv2_s_d.c*

$(SRC_PATH)/BLAS_dgbmv2_s_s.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(d, s, s)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgbmv2_s_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv2_s_s.c && \
    mv tmp.BLAS_dgbmv2_s_s.c $@ && rm -f tmp.BLAS_dgbmv2_s_s.c*

$(SRC_PATH)/BLAS_zgbmv2_z_c.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(z, z, c)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgbmv2_z_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv2_z_c.c && \
    mv tmp.BLAS_zgbmv2_z_c.c $@ && rm -f tmp.BLAS_zgbmv2_z_c.c*

$(SRC_PATH)/BLAS_zgbmv2_c_z.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(z, c, z)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgbmv2_c_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv2_c_z.c && \
    mv tmp.BLAS_zgbmv2_c_z.c $@ && rm -f tmp.BLAS_zgbmv2_c_z.c*

$(SRC_PATH)/BLAS_zgbmv2_c_c.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(z, c, c)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgbmv2_c_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv2_c_c.c && \
    mv tmp.BLAS_zgbmv2_c_c.c $@ && rm -f tmp.BLAS_zgbmv2_c_c.c*

$(SRC_PATH)/BLAS_cgbmv2_c_s.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(c, c, s)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgbmv2_c_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv2_c_s.c && \
    mv tmp.BLAS_cgbmv2_c_s.c $@ && rm -f tmp.BLAS_cgbmv2_c_s.c*

$(SRC_PATH)/BLAS_cgbmv2_s_c.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(c, s, c)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgbmv2_s_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv2_s_c.c && \
    mv tmp.BLAS_cgbmv2_s_c.c $@ && rm -f tmp.BLAS_cgbmv2_s_c.c*

$(SRC_PATH)/BLAS_cgbmv2_s_s.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(c, s, s)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgbmv2_s_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv2_s_s.c && \
    mv tmp.BLAS_cgbmv2_s_s.c $@ && rm -f tmp.BLAS_cgbmv2_s_s.c*

$(SRC_PATH)/BLAS_zgbmv2_z_d.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(z, z, d)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgbmv2_z_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv2_z_d.c && \
    mv tmp.BLAS_zgbmv2_z_d.c $@ && rm -f tmp.BLAS_zgbmv2_z_d.c*

$(SRC_PATH)/BLAS_zgbmv2_d_z.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(z, d, z)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgbmv2_d_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv2_d_z.c && \
    mv tmp.BLAS_zgbmv2_d_z.c $@ && rm -f tmp.BLAS_zgbmv2_d_z.c*

$(SRC_PATH)/BLAS_zgbmv2_d_d.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(z, d, d)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgbmv2_d_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv2_d_d.c && \
    mv tmp.BLAS_zgbmv2_d_d.c $@ && rm -f tmp.BLAS_zgbmv2_d_d.c*

$(SRC_PATH)/BLAS_sgbmv2_x.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(s, s, s, _x)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_sgbmv2_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sgbmv2_x.c && \
    mv tmp.BLAS_sgbmv2_x.c $@ && rm -f tmp.BLAS_sgbmv2_x.c*

$(SRC_PATH)/BLAS_dgbmv2_x.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(d, d, d, _x)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgbmv2_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv2_x.c && \
    mv tmp.BLAS_dgbmv2_x.c $@ && rm -f tmp.BLAS_dgbmv2_x.c*

$(SRC_PATH)/BLAS_cgbmv2_x.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(c, c, c, _x)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgbmv2_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv2_x.c && \
    mv tmp.BLAS_cgbmv2_x.c $@ && rm -f tmp.BLAS_cgbmv2_x.c*

$(SRC_PATH)/BLAS_zgbmv2_x.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(z, z, z, _x)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgbmv2_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv2_x.c && \
    mv tmp.BLAS_zgbmv2_x.c $@ && rm -f tmp.BLAS_zgbmv2_x.c*

$(SRC_PATH)/BLAS_dgbmv2_d_s_x.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(d, d, s, _x)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgbmv2_d_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv2_d_s_x.c && \
    mv tmp.BLAS_dgbmv2_d_s_x.c $@ && rm -f tmp.BLAS_dgbmv2_d_s_x.c*

$(SRC_PATH)/BLAS_dgbmv2_s_d_x.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(d, s, d, _x)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgbmv2_s_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv2_s_d_x.c && \
    mv tmp.BLAS_dgbmv2_s_d_x.c $@ && rm -f tmp.BLAS_dgbmv2_s_d_x.c*

$(SRC_PATH)/BLAS_dgbmv2_s_s_x.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(d, s, s, _x)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgbmv2_s_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv2_s_s_x.c && \
    mv tmp.BLAS_dgbmv2_s_s_x.c $@ && rm -f tmp.BLAS_dgbmv2_s_s_x.c*

$(SRC_PATH)/BLAS_zgbmv2_z_c_x.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(z, z, c, _x)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgbmv2_z_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv2_z_c_x.c && \
    mv tmp.BLAS_zgbmv2_z_c_x.c $@ && rm -f tmp.BLAS_zgbmv2_z_c_x.c*

$(SRC_PATH)/BLAS_zgbmv2_c_z_x.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(z, c, z, _x)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgbmv2_c_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv2_c_z_x.c && \
    mv tmp.BLAS_zgbmv2_c_z_x.c $@ && rm -f tmp.BLAS_zgbmv2_c_z_x.c*

$(SRC_PATH)/BLAS_zgbmv2_c_c_x.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(z, c, c, _x)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgbmv2_c_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv2_c_c_x.c && \
    mv tmp.BLAS_zgbmv2_c_c_x.c $@ && rm -f tmp.BLAS_zgbmv2_c_c_x.c*

$(SRC_PATH)/BLAS_cgbmv2_c_s_x.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(c, c, s, _x)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgbmv2_c_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv2_c_s_x.c && \
    mv tmp.BLAS_cgbmv2_c_s_x.c $@ && rm -f tmp.BLAS_cgbmv2_c_s_x.c*

$(SRC_PATH)/BLAS_cgbmv2_s_c_x.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(c, s, c, _x)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgbmv2_s_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv2_s_c_x.c && \
    mv tmp.BLAS_cgbmv2_s_c_x.c $@ && rm -f tmp.BLAS_cgbmv2_s_c_x.c*

$(SRC_PATH)/BLAS_cgbmv2_s_s_x.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(c, s, s, _x)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgbmv2_s_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv2_s_s_x.c && \
    mv tmp.BLAS_cgbmv2_s_s_x.c $@ && rm -f tmp.BLAS_cgbmv2_s_s_x.c*

$(SRC_PATH)/BLAS_zgbmv2_z_d_x.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(z, z, d, _x)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgbmv2_z_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv2_z_d_x.c && \
    mv tmp.BLAS_zgbmv2_z_d_x.c $@ && rm -f tmp.BLAS_zgbmv2_z_d_x.c*

$(SRC_PATH)/BLAS_zgbmv2_d_z_x.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(z, d, z, _x)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgbmv2_d_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv2_d_z_x.c && \
    mv tmp.BLAS_zgbmv2_d_z_x.c $@ && rm -f tmp.BLAS_zgbmv2_d_z_x.c*

$(SRC_PATH)/BLAS_zgbmv2_d_d_x.c: gbmv2.m4 gbmv2-common.m4 ../cblas.m4
	echo "GBMV2(z, d, d, _x)" | cat gbmv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgbmv2_d_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv2_d_d_x.c && \
    mv tmp.BLAS_zgbmv2_d_d_x.c $@ && rm -f tmp.BLAS_zgbmv2_d_d_x.c*



$(SRC_PATH)/BLAS_dgbmv2_d_s-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="d, d, s" ../f2c-bridge.m4 > tmp.BLAS_dgbmv2_d_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv2_d_s-f2c.c && \
    mv tmp.BLAS_dgbmv2_d_s-f2c.c $@ && rm -f tmp.BLAS_dgbmv2_d_s-f2c.c*

$(SRC_PATH)/BLAS_dgbmv2_s_d-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="d, s, d" ../f2c-bridge.m4 > tmp.BLAS_dgbmv2_s_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv2_s_d-f2c.c && \
    mv tmp.BLAS_dgbmv2_s_d-f2c.c $@ && rm -f tmp.BLAS_dgbmv2_s_d-f2c.c*

$(SRC_PATH)/BLAS_dgbmv2_s_s-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="d, s, s" ../f2c-bridge.m4 > tmp.BLAS_dgbmv2_s_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv2_s_s-f2c.c && \
    mv tmp.BLAS_dgbmv2_s_s-f2c.c $@ && rm -f tmp.BLAS_dgbmv2_s_s-f2c.c*

$(SRC_PATH)/BLAS_zgbmv2_z_c-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="z, z, c" ../f2c-bridge.m4 > tmp.BLAS_zgbmv2_z_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv2_z_c-f2c.c && \
    mv tmp.BLAS_zgbmv2_z_c-f2c.c $@ && rm -f tmp.BLAS_zgbmv2_z_c-f2c.c*

$(SRC_PATH)/BLAS_zgbmv2_c_z-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="z, c, z" ../f2c-bridge.m4 > tmp.BLAS_zgbmv2_c_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv2_c_z-f2c.c && \
    mv tmp.BLAS_zgbmv2_c_z-f2c.c $@ && rm -f tmp.BLAS_zgbmv2_c_z-f2c.c*

$(SRC_PATH)/BLAS_zgbmv2_c_c-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="z, c, c" ../f2c-bridge.m4 > tmp.BLAS_zgbmv2_c_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv2_c_c-f2c.c && \
    mv tmp.BLAS_zgbmv2_c_c-f2c.c $@ && rm -f tmp.BLAS_zgbmv2_c_c-f2c.c*

$(SRC_PATH)/BLAS_cgbmv2_c_s-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="c, c, s" ../f2c-bridge.m4 > tmp.BLAS_cgbmv2_c_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv2_c_s-f2c.c && \
    mv tmp.BLAS_cgbmv2_c_s-f2c.c $@ && rm -f tmp.BLAS_cgbmv2_c_s-f2c.c*

$(SRC_PATH)/BLAS_cgbmv2_s_c-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="c, s, c" ../f2c-bridge.m4 > tmp.BLAS_cgbmv2_s_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv2_s_c-f2c.c && \
    mv tmp.BLAS_cgbmv2_s_c-f2c.c $@ && rm -f tmp.BLAS_cgbmv2_s_c-f2c.c*

$(SRC_PATH)/BLAS_cgbmv2_s_s-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="c, s, s" ../f2c-bridge.m4 > tmp.BLAS_cgbmv2_s_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv2_s_s-f2c.c && \
    mv tmp.BLAS_cgbmv2_s_s-f2c.c $@ && rm -f tmp.BLAS_cgbmv2_s_s-f2c.c*

$(SRC_PATH)/BLAS_zgbmv2_z_d-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="z, z, d" ../f2c-bridge.m4 > tmp.BLAS_zgbmv2_z_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv2_z_d-f2c.c && \
    mv tmp.BLAS_zgbmv2_z_d-f2c.c $@ && rm -f tmp.BLAS_zgbmv2_z_d-f2c.c*

$(SRC_PATH)/BLAS_zgbmv2_d_z-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="z, d, z" ../f2c-bridge.m4 > tmp.BLAS_zgbmv2_d_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv2_d_z-f2c.c && \
    mv tmp.BLAS_zgbmv2_d_z-f2c.c $@ && rm -f tmp.BLAS_zgbmv2_d_z-f2c.c*

$(SRC_PATH)/BLAS_zgbmv2_d_d-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="z, d, d" ../f2c-bridge.m4 > tmp.BLAS_zgbmv2_d_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv2_d_d-f2c.c && \
    mv tmp.BLAS_zgbmv2_d_d-f2c.c $@ && rm -f tmp.BLAS_zgbmv2_d_d-f2c.c*

$(SRC_PATH)/BLAS_sgbmv2_x-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="s, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_sgbmv2_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sgbmv2_x-f2c.c && \
    mv tmp.BLAS_sgbmv2_x-f2c.c $@ && rm -f tmp.BLAS_sgbmv2_x-f2c.c*

$(SRC_PATH)/BLAS_dgbmv2_x-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="d, d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dgbmv2_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv2_x-f2c.c && \
    mv tmp.BLAS_dgbmv2_x-f2c.c $@ && rm -f tmp.BLAS_dgbmv2_x-f2c.c*

$(SRC_PATH)/BLAS_cgbmv2_x-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="c, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_cgbmv2_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv2_x-f2c.c && \
    mv tmp.BLAS_cgbmv2_x-f2c.c $@ && rm -f tmp.BLAS_cgbmv2_x-f2c.c*

$(SRC_PATH)/BLAS_zgbmv2_x-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="z, z, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zgbmv2_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv2_x-f2c.c && \
    mv tmp.BLAS_zgbmv2_x-f2c.c $@ && rm -f tmp.BLAS_zgbmv2_x-f2c.c*

$(SRC_PATH)/BLAS_dgbmv2_d_s_x-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="d, d, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dgbmv2_d_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv2_d_s_x-f2c.c && \
    mv tmp.BLAS_dgbmv2_d_s_x-f2c.c $@ && rm -f tmp.BLAS_dgbmv2_d_s_x-f2c.c*

$(SRC_PATH)/BLAS_dgbmv2_s_d_x-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="d, s, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dgbmv2_s_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv2_s_d_x-f2c.c && \
    mv tmp.BLAS_dgbmv2_s_d_x-f2c.c $@ && rm -f tmp.BLAS_dgbmv2_s_d_x-f2c.c*

$(SRC_PATH)/BLAS_dgbmv2_s_s_x-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="d, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dgbmv2_s_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgbmv2_s_s_x-f2c.c && \
    mv tmp.BLAS_dgbmv2_s_s_x-f2c.c $@ && rm -f tmp.BLAS_dgbmv2_s_s_x-f2c.c*

$(SRC_PATH)/BLAS_zgbmv2_z_c_x-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="z, z, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zgbmv2_z_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv2_z_c_x-f2c.c && \
    mv tmp.BLAS_zgbmv2_z_c_x-f2c.c $@ && rm -f tmp.BLAS_zgbmv2_z_c_x-f2c.c*

$(SRC_PATH)/BLAS_zgbmv2_c_z_x-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="z, c, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zgbmv2_c_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv2_c_z_x-f2c.c && \
    mv tmp.BLAS_zgbmv2_c_z_x-f2c.c $@ && rm -f tmp.BLAS_zgbmv2_c_z_x-f2c.c*

$(SRC_PATH)/BLAS_zgbmv2_c_c_x-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="z, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zgbmv2_c_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv2_c_c_x-f2c.c && \
    mv tmp.BLAS_zgbmv2_c_c_x-f2c.c $@ && rm -f tmp.BLAS_zgbmv2_c_c_x-f2c.c*

$(SRC_PATH)/BLAS_cgbmv2_c_s_x-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="c, c, s, _x" ../f2c-bridge.m4 > tmp.BLAS_cgbmv2_c_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv2_c_s_x-f2c.c && \
    mv tmp.BLAS_cgbmv2_c_s_x-f2c.c $@ && rm -f tmp.BLAS_cgbmv2_c_s_x-f2c.c*

$(SRC_PATH)/BLAS_cgbmv2_s_c_x-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="c, s, c, _x" ../f2c-bridge.m4 > tmp.BLAS_cgbmv2_s_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv2_s_c_x-f2c.c && \
    mv tmp.BLAS_cgbmv2_s_c_x-f2c.c $@ && rm -f tmp.BLAS_cgbmv2_s_c_x-f2c.c*

$(SRC_PATH)/BLAS_cgbmv2_s_s_x-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="c, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_cgbmv2_s_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgbmv2_s_s_x-f2c.c && \
    mv tmp.BLAS_cgbmv2_s_s_x-f2c.c $@ && rm -f tmp.BLAS_cgbmv2_s_s_x-f2c.c*

$(SRC_PATH)/BLAS_zgbmv2_z_d_x-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="z, z, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zgbmv2_z_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv2_z_d_x-f2c.c && \
    mv tmp.BLAS_zgbmv2_z_d_x-f2c.c $@ && rm -f tmp.BLAS_zgbmv2_z_d_x-f2c.c*

$(SRC_PATH)/BLAS_zgbmv2_d_z_x-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="z, d, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zgbmv2_d_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv2_d_z_x-f2c.c && \
    mv tmp.BLAS_zgbmv2_d_z_x-f2c.c $@ && rm -f tmp.BLAS_zgbmv2_d_z_x-f2c.c*

$(SRC_PATH)/BLAS_zgbmv2_d_d_x-f2c.c: ../f2c-bridge.m4 gbmv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gbmv2 -DRARG="z, d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zgbmv2_d_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgbmv2_d_d_x-f2c.c && \
    mv tmp.BLAS_zgbmv2_d_d_x-f2c.c $@ && rm -f tmp.BLAS_zgbmv2_d_d_x-f2c.c*


gbmv2.h: ../cblas.m4 ../header.template gbmv2-common.m4
	sed -e "s/routine/gbmv2/g" -e "s/ROUTINE/GBMV2/g" \
  ../header.template | $(M4) $(M4_OPTS) $(M4_PATH) - >$@

source-clean:
	rm -f $(GBMV2_SRCS) gbmv2.h

.PHONY: gbmv2-source gbmv2-header source-clean

