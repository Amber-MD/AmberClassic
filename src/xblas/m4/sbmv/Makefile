include ../../make.conf
include ../../$(MAKEINC)


SRC_PATH = ../../src/sbmv
M4_PATH = -I. -I..

SBMV_SRCS = \
	$(SRC_PATH)/BLAS_ssbmv.c\
	$(SRC_PATH)/BLAS_dsbmv.c\
	$(SRC_PATH)/BLAS_csbmv.c\
	$(SRC_PATH)/BLAS_zsbmv.c\
	$(SRC_PATH)/BLAS_dsbmv_d_s.c\
	$(SRC_PATH)/BLAS_dsbmv_s_d.c\
	$(SRC_PATH)/BLAS_dsbmv_s_s.c\
	$(SRC_PATH)/BLAS_zsbmv_z_c.c\
	$(SRC_PATH)/BLAS_zsbmv_c_z.c\
	$(SRC_PATH)/BLAS_zsbmv_c_c.c\
	$(SRC_PATH)/BLAS_csbmv_c_s.c\
	$(SRC_PATH)/BLAS_csbmv_s_c.c\
	$(SRC_PATH)/BLAS_csbmv_s_s.c\
	$(SRC_PATH)/BLAS_zsbmv_z_d.c\
	$(SRC_PATH)/BLAS_zsbmv_d_z.c\
	$(SRC_PATH)/BLAS_zsbmv_d_d.c\
	$(SRC_PATH)/BLAS_ssbmv_x.c\
	$(SRC_PATH)/BLAS_dsbmv_x.c\
	$(SRC_PATH)/BLAS_csbmv_x.c\
	$(SRC_PATH)/BLAS_zsbmv_x.c\
	$(SRC_PATH)/BLAS_dsbmv_d_s_x.c\
	$(SRC_PATH)/BLAS_dsbmv_s_d_x.c\
	$(SRC_PATH)/BLAS_dsbmv_s_s_x.c\
	$(SRC_PATH)/BLAS_zsbmv_z_c_x.c\
	$(SRC_PATH)/BLAS_zsbmv_c_z_x.c\
	$(SRC_PATH)/BLAS_zsbmv_c_c_x.c\
	$(SRC_PATH)/BLAS_csbmv_c_s_x.c\
	$(SRC_PATH)/BLAS_csbmv_s_c_x.c\
	$(SRC_PATH)/BLAS_csbmv_s_s_x.c\
	$(SRC_PATH)/BLAS_zsbmv_z_d_x.c\
	$(SRC_PATH)/BLAS_zsbmv_d_z_x.c\
	$(SRC_PATH)/BLAS_zsbmv_d_d_x.c \
	\
	$(SRC_PATH)/BLAS_ssbmv-f2c.c\
	$(SRC_PATH)/BLAS_dsbmv-f2c.c\
	$(SRC_PATH)/BLAS_csbmv-f2c.c\
	$(SRC_PATH)/BLAS_zsbmv-f2c.c\
	$(SRC_PATH)/BLAS_dsbmv_d_s-f2c.c\
	$(SRC_PATH)/BLAS_dsbmv_s_d-f2c.c\
	$(SRC_PATH)/BLAS_dsbmv_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zsbmv_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zsbmv_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zsbmv_c_c-f2c.c\
	$(SRC_PATH)/BLAS_csbmv_c_s-f2c.c\
	$(SRC_PATH)/BLAS_csbmv_s_c-f2c.c\
	$(SRC_PATH)/BLAS_csbmv_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zsbmv_z_d-f2c.c\
	$(SRC_PATH)/BLAS_zsbmv_d_z-f2c.c\
	$(SRC_PATH)/BLAS_zsbmv_d_d-f2c.c\
	$(SRC_PATH)/BLAS_ssbmv_x-f2c.c\
	$(SRC_PATH)/BLAS_dsbmv_x-f2c.c\
	$(SRC_PATH)/BLAS_csbmv_x-f2c.c\
	$(SRC_PATH)/BLAS_zsbmv_x-f2c.c\
	$(SRC_PATH)/BLAS_dsbmv_d_s_x-f2c.c\
	$(SRC_PATH)/BLAS_dsbmv_s_d_x-f2c.c\
	$(SRC_PATH)/BLAS_dsbmv_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zsbmv_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zsbmv_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zsbmv_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_csbmv_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_csbmv_s_c_x-f2c.c\
	$(SRC_PATH)/BLAS_csbmv_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zsbmv_z_d_x-f2c.c\
	$(SRC_PATH)/BLAS_zsbmv_d_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zsbmv_d_d_x-f2c.c

sbmv-source: $(SBMV_SRCS)

sbmv-header: sbmv.h

$(SRC_PATH)/BLAS_ssbmv.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(s, s, s)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ssbmv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ssbmv.c && \
    mv tmp.BLAS_ssbmv.c $@ && rm -f tmp.BLAS_ssbmv.c*

$(SRC_PATH)/BLAS_dsbmv.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(d, d, d)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsbmv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsbmv.c && \
    mv tmp.BLAS_dsbmv.c $@ && rm -f tmp.BLAS_dsbmv.c*

$(SRC_PATH)/BLAS_csbmv.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(c, c, c)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csbmv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csbmv.c && \
    mv tmp.BLAS_csbmv.c $@ && rm -f tmp.BLAS_csbmv.c*

$(SRC_PATH)/BLAS_zsbmv.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(z, z, z)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsbmv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv.c && \
    mv tmp.BLAS_zsbmv.c $@ && rm -f tmp.BLAS_zsbmv.c*

$(SRC_PATH)/BLAS_dsbmv_d_s.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(d, d, s)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsbmv_d_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsbmv_d_s.c && \
    mv tmp.BLAS_dsbmv_d_s.c $@ && rm -f tmp.BLAS_dsbmv_d_s.c*

$(SRC_PATH)/BLAS_dsbmv_s_d.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(d, s, d)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsbmv_s_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsbmv_s_d.c && \
    mv tmp.BLAS_dsbmv_s_d.c $@ && rm -f tmp.BLAS_dsbmv_s_d.c*

$(SRC_PATH)/BLAS_dsbmv_s_s.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(d, s, s)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsbmv_s_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsbmv_s_s.c && \
    mv tmp.BLAS_dsbmv_s_s.c $@ && rm -f tmp.BLAS_dsbmv_s_s.c*

$(SRC_PATH)/BLAS_zsbmv_z_c.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(z, z, c)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsbmv_z_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv_z_c.c && \
    mv tmp.BLAS_zsbmv_z_c.c $@ && rm -f tmp.BLAS_zsbmv_z_c.c*

$(SRC_PATH)/BLAS_zsbmv_c_z.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(z, c, z)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsbmv_c_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv_c_z.c && \
    mv tmp.BLAS_zsbmv_c_z.c $@ && rm -f tmp.BLAS_zsbmv_c_z.c*

$(SRC_PATH)/BLAS_zsbmv_c_c.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(z, c, c)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsbmv_c_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv_c_c.c && \
    mv tmp.BLAS_zsbmv_c_c.c $@ && rm -f tmp.BLAS_zsbmv_c_c.c*

$(SRC_PATH)/BLAS_csbmv_c_s.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(c, c, s)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csbmv_c_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csbmv_c_s.c && \
    mv tmp.BLAS_csbmv_c_s.c $@ && rm -f tmp.BLAS_csbmv_c_s.c*

$(SRC_PATH)/BLAS_csbmv_s_c.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(c, s, c)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csbmv_s_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csbmv_s_c.c && \
    mv tmp.BLAS_csbmv_s_c.c $@ && rm -f tmp.BLAS_csbmv_s_c.c*

$(SRC_PATH)/BLAS_csbmv_s_s.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(c, s, s)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csbmv_s_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csbmv_s_s.c && \
    mv tmp.BLAS_csbmv_s_s.c $@ && rm -f tmp.BLAS_csbmv_s_s.c*

$(SRC_PATH)/BLAS_zsbmv_z_d.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(z, z, d)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsbmv_z_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv_z_d.c && \
    mv tmp.BLAS_zsbmv_z_d.c $@ && rm -f tmp.BLAS_zsbmv_z_d.c*

$(SRC_PATH)/BLAS_zsbmv_d_z.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(z, d, z)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsbmv_d_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv_d_z.c && \
    mv tmp.BLAS_zsbmv_d_z.c $@ && rm -f tmp.BLAS_zsbmv_d_z.c*

$(SRC_PATH)/BLAS_zsbmv_d_d.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(z, d, d)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsbmv_d_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv_d_d.c && \
    mv tmp.BLAS_zsbmv_d_d.c $@ && rm -f tmp.BLAS_zsbmv_d_d.c*

$(SRC_PATH)/BLAS_ssbmv_x.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(s, s, s, _x)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ssbmv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ssbmv_x.c && \
    mv tmp.BLAS_ssbmv_x.c $@ && rm -f tmp.BLAS_ssbmv_x.c*

$(SRC_PATH)/BLAS_dsbmv_x.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(d, d, d, _x)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsbmv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsbmv_x.c && \
    mv tmp.BLAS_dsbmv_x.c $@ && rm -f tmp.BLAS_dsbmv_x.c*

$(SRC_PATH)/BLAS_csbmv_x.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(c, c, c, _x)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csbmv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csbmv_x.c && \
    mv tmp.BLAS_csbmv_x.c $@ && rm -f tmp.BLAS_csbmv_x.c*

$(SRC_PATH)/BLAS_zsbmv_x.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(z, z, z, _x)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsbmv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv_x.c && \
    mv tmp.BLAS_zsbmv_x.c $@ && rm -f tmp.BLAS_zsbmv_x.c*

$(SRC_PATH)/BLAS_dsbmv_d_s_x.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(d, d, s, _x)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsbmv_d_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsbmv_d_s_x.c && \
    mv tmp.BLAS_dsbmv_d_s_x.c $@ && rm -f tmp.BLAS_dsbmv_d_s_x.c*

$(SRC_PATH)/BLAS_dsbmv_s_d_x.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(d, s, d, _x)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsbmv_s_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsbmv_s_d_x.c && \
    mv tmp.BLAS_dsbmv_s_d_x.c $@ && rm -f tmp.BLAS_dsbmv_s_d_x.c*

$(SRC_PATH)/BLAS_dsbmv_s_s_x.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(d, s, s, _x)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsbmv_s_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsbmv_s_s_x.c && \
    mv tmp.BLAS_dsbmv_s_s_x.c $@ && rm -f tmp.BLAS_dsbmv_s_s_x.c*

$(SRC_PATH)/BLAS_zsbmv_z_c_x.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(z, z, c, _x)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsbmv_z_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv_z_c_x.c && \
    mv tmp.BLAS_zsbmv_z_c_x.c $@ && rm -f tmp.BLAS_zsbmv_z_c_x.c*

$(SRC_PATH)/BLAS_zsbmv_c_z_x.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(z, c, z, _x)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsbmv_c_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv_c_z_x.c && \
    mv tmp.BLAS_zsbmv_c_z_x.c $@ && rm -f tmp.BLAS_zsbmv_c_z_x.c*

$(SRC_PATH)/BLAS_zsbmv_c_c_x.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(z, c, c, _x)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsbmv_c_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv_c_c_x.c && \
    mv tmp.BLAS_zsbmv_c_c_x.c $@ && rm -f tmp.BLAS_zsbmv_c_c_x.c*

$(SRC_PATH)/BLAS_csbmv_c_s_x.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(c, c, s, _x)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csbmv_c_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csbmv_c_s_x.c && \
    mv tmp.BLAS_csbmv_c_s_x.c $@ && rm -f tmp.BLAS_csbmv_c_s_x.c*

$(SRC_PATH)/BLAS_csbmv_s_c_x.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(c, s, c, _x)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csbmv_s_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csbmv_s_c_x.c && \
    mv tmp.BLAS_csbmv_s_c_x.c $@ && rm -f tmp.BLAS_csbmv_s_c_x.c*

$(SRC_PATH)/BLAS_csbmv_s_s_x.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(c, s, s, _x)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csbmv_s_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csbmv_s_s_x.c && \
    mv tmp.BLAS_csbmv_s_s_x.c $@ && rm -f tmp.BLAS_csbmv_s_s_x.c*

$(SRC_PATH)/BLAS_zsbmv_z_d_x.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(z, z, d, _x)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsbmv_z_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv_z_d_x.c && \
    mv tmp.BLAS_zsbmv_z_d_x.c $@ && rm -f tmp.BLAS_zsbmv_z_d_x.c*

$(SRC_PATH)/BLAS_zsbmv_d_z_x.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(z, d, z, _x)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsbmv_d_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv_d_z_x.c && \
    mv tmp.BLAS_zsbmv_d_z_x.c $@ && rm -f tmp.BLAS_zsbmv_d_z_x.c*

$(SRC_PATH)/BLAS_zsbmv_d_d_x.c: sbmv.m4 sbmv-common.m4 ../cblas.m4
	echo "SBMV(z, d, d, _x)" | cat sbmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsbmv_d_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv_d_d_x.c && \
    mv tmp.BLAS_zsbmv_d_d_x.c $@ && rm -f tmp.BLAS_zsbmv_d_d_x.c*



$(SRC_PATH)/BLAS_ssbmv-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="s, s, s" ../f2c-bridge.m4 > tmp.BLAS_ssbmv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ssbmv-f2c.c && \
    mv tmp.BLAS_ssbmv-f2c.c $@ && rm -f tmp.BLAS_ssbmv-f2c.c*

$(SRC_PATH)/BLAS_dsbmv-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="d, d, d" ../f2c-bridge.m4 > tmp.BLAS_dsbmv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsbmv-f2c.c && \
    mv tmp.BLAS_dsbmv-f2c.c $@ && rm -f tmp.BLAS_dsbmv-f2c.c*

$(SRC_PATH)/BLAS_csbmv-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="c, c, c" ../f2c-bridge.m4 > tmp.BLAS_csbmv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csbmv-f2c.c && \
    mv tmp.BLAS_csbmv-f2c.c $@ && rm -f tmp.BLAS_csbmv-f2c.c*

$(SRC_PATH)/BLAS_zsbmv-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="z, z, z" ../f2c-bridge.m4 > tmp.BLAS_zsbmv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv-f2c.c && \
    mv tmp.BLAS_zsbmv-f2c.c $@ && rm -f tmp.BLAS_zsbmv-f2c.c*

$(SRC_PATH)/BLAS_dsbmv_d_s-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="d, d, s" ../f2c-bridge.m4 > tmp.BLAS_dsbmv_d_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsbmv_d_s-f2c.c && \
    mv tmp.BLAS_dsbmv_d_s-f2c.c $@ && rm -f tmp.BLAS_dsbmv_d_s-f2c.c*

$(SRC_PATH)/BLAS_dsbmv_s_d-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="d, s, d" ../f2c-bridge.m4 > tmp.BLAS_dsbmv_s_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsbmv_s_d-f2c.c && \
    mv tmp.BLAS_dsbmv_s_d-f2c.c $@ && rm -f tmp.BLAS_dsbmv_s_d-f2c.c*

$(SRC_PATH)/BLAS_dsbmv_s_s-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="d, s, s" ../f2c-bridge.m4 > tmp.BLAS_dsbmv_s_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsbmv_s_s-f2c.c && \
    mv tmp.BLAS_dsbmv_s_s-f2c.c $@ && rm -f tmp.BLAS_dsbmv_s_s-f2c.c*

$(SRC_PATH)/BLAS_zsbmv_z_c-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="z, z, c" ../f2c-bridge.m4 > tmp.BLAS_zsbmv_z_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv_z_c-f2c.c && \
    mv tmp.BLAS_zsbmv_z_c-f2c.c $@ && rm -f tmp.BLAS_zsbmv_z_c-f2c.c*

$(SRC_PATH)/BLAS_zsbmv_c_z-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="z, c, z" ../f2c-bridge.m4 > tmp.BLAS_zsbmv_c_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv_c_z-f2c.c && \
    mv tmp.BLAS_zsbmv_c_z-f2c.c $@ && rm -f tmp.BLAS_zsbmv_c_z-f2c.c*

$(SRC_PATH)/BLAS_zsbmv_c_c-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="z, c, c" ../f2c-bridge.m4 > tmp.BLAS_zsbmv_c_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv_c_c-f2c.c && \
    mv tmp.BLAS_zsbmv_c_c-f2c.c $@ && rm -f tmp.BLAS_zsbmv_c_c-f2c.c*

$(SRC_PATH)/BLAS_csbmv_c_s-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="c, c, s" ../f2c-bridge.m4 > tmp.BLAS_csbmv_c_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csbmv_c_s-f2c.c && \
    mv tmp.BLAS_csbmv_c_s-f2c.c $@ && rm -f tmp.BLAS_csbmv_c_s-f2c.c*

$(SRC_PATH)/BLAS_csbmv_s_c-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="c, s, c" ../f2c-bridge.m4 > tmp.BLAS_csbmv_s_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csbmv_s_c-f2c.c && \
    mv tmp.BLAS_csbmv_s_c-f2c.c $@ && rm -f tmp.BLAS_csbmv_s_c-f2c.c*

$(SRC_PATH)/BLAS_csbmv_s_s-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="c, s, s" ../f2c-bridge.m4 > tmp.BLAS_csbmv_s_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csbmv_s_s-f2c.c && \
    mv tmp.BLAS_csbmv_s_s-f2c.c $@ && rm -f tmp.BLAS_csbmv_s_s-f2c.c*

$(SRC_PATH)/BLAS_zsbmv_z_d-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="z, z, d" ../f2c-bridge.m4 > tmp.BLAS_zsbmv_z_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv_z_d-f2c.c && \
    mv tmp.BLAS_zsbmv_z_d-f2c.c $@ && rm -f tmp.BLAS_zsbmv_z_d-f2c.c*

$(SRC_PATH)/BLAS_zsbmv_d_z-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="z, d, z" ../f2c-bridge.m4 > tmp.BLAS_zsbmv_d_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv_d_z-f2c.c && \
    mv tmp.BLAS_zsbmv_d_z-f2c.c $@ && rm -f tmp.BLAS_zsbmv_d_z-f2c.c*

$(SRC_PATH)/BLAS_zsbmv_d_d-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="z, d, d" ../f2c-bridge.m4 > tmp.BLAS_zsbmv_d_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv_d_d-f2c.c && \
    mv tmp.BLAS_zsbmv_d_d-f2c.c $@ && rm -f tmp.BLAS_zsbmv_d_d-f2c.c*

$(SRC_PATH)/BLAS_ssbmv_x-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="s, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_ssbmv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ssbmv_x-f2c.c && \
    mv tmp.BLAS_ssbmv_x-f2c.c $@ && rm -f tmp.BLAS_ssbmv_x-f2c.c*

$(SRC_PATH)/BLAS_dsbmv_x-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="d, d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dsbmv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsbmv_x-f2c.c && \
    mv tmp.BLAS_dsbmv_x-f2c.c $@ && rm -f tmp.BLAS_dsbmv_x-f2c.c*

$(SRC_PATH)/BLAS_csbmv_x-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="c, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_csbmv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csbmv_x-f2c.c && \
    mv tmp.BLAS_csbmv_x-f2c.c $@ && rm -f tmp.BLAS_csbmv_x-f2c.c*

$(SRC_PATH)/BLAS_zsbmv_x-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="z, z, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zsbmv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv_x-f2c.c && \
    mv tmp.BLAS_zsbmv_x-f2c.c $@ && rm -f tmp.BLAS_zsbmv_x-f2c.c*

$(SRC_PATH)/BLAS_dsbmv_d_s_x-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="d, d, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dsbmv_d_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsbmv_d_s_x-f2c.c && \
    mv tmp.BLAS_dsbmv_d_s_x-f2c.c $@ && rm -f tmp.BLAS_dsbmv_d_s_x-f2c.c*

$(SRC_PATH)/BLAS_dsbmv_s_d_x-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="d, s, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dsbmv_s_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsbmv_s_d_x-f2c.c && \
    mv tmp.BLAS_dsbmv_s_d_x-f2c.c $@ && rm -f tmp.BLAS_dsbmv_s_d_x-f2c.c*

$(SRC_PATH)/BLAS_dsbmv_s_s_x-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="d, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dsbmv_s_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsbmv_s_s_x-f2c.c && \
    mv tmp.BLAS_dsbmv_s_s_x-f2c.c $@ && rm -f tmp.BLAS_dsbmv_s_s_x-f2c.c*

$(SRC_PATH)/BLAS_zsbmv_z_c_x-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="z, z, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zsbmv_z_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv_z_c_x-f2c.c && \
    mv tmp.BLAS_zsbmv_z_c_x-f2c.c $@ && rm -f tmp.BLAS_zsbmv_z_c_x-f2c.c*

$(SRC_PATH)/BLAS_zsbmv_c_z_x-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="z, c, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zsbmv_c_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv_c_z_x-f2c.c && \
    mv tmp.BLAS_zsbmv_c_z_x-f2c.c $@ && rm -f tmp.BLAS_zsbmv_c_z_x-f2c.c*

$(SRC_PATH)/BLAS_zsbmv_c_c_x-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="z, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zsbmv_c_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv_c_c_x-f2c.c && \
    mv tmp.BLAS_zsbmv_c_c_x-f2c.c $@ && rm -f tmp.BLAS_zsbmv_c_c_x-f2c.c*

$(SRC_PATH)/BLAS_csbmv_c_s_x-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="c, c, s, _x" ../f2c-bridge.m4 > tmp.BLAS_csbmv_c_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csbmv_c_s_x-f2c.c && \
    mv tmp.BLAS_csbmv_c_s_x-f2c.c $@ && rm -f tmp.BLAS_csbmv_c_s_x-f2c.c*

$(SRC_PATH)/BLAS_csbmv_s_c_x-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="c, s, c, _x" ../f2c-bridge.m4 > tmp.BLAS_csbmv_s_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csbmv_s_c_x-f2c.c && \
    mv tmp.BLAS_csbmv_s_c_x-f2c.c $@ && rm -f tmp.BLAS_csbmv_s_c_x-f2c.c*

$(SRC_PATH)/BLAS_csbmv_s_s_x-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="c, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_csbmv_s_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csbmv_s_s_x-f2c.c && \
    mv tmp.BLAS_csbmv_s_s_x-f2c.c $@ && rm -f tmp.BLAS_csbmv_s_s_x-f2c.c*

$(SRC_PATH)/BLAS_zsbmv_z_d_x-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="z, z, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zsbmv_z_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv_z_d_x-f2c.c && \
    mv tmp.BLAS_zsbmv_z_d_x-f2c.c $@ && rm -f tmp.BLAS_zsbmv_z_d_x-f2c.c*

$(SRC_PATH)/BLAS_zsbmv_d_z_x-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="z, d, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zsbmv_d_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv_d_z_x-f2c.c && \
    mv tmp.BLAS_zsbmv_d_z_x-f2c.c $@ && rm -f tmp.BLAS_zsbmv_d_z_x-f2c.c*

$(SRC_PATH)/BLAS_zsbmv_d_d_x-f2c.c: ../f2c-bridge.m4 sbmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sbmv -DRARG="z, d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zsbmv_d_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsbmv_d_d_x-f2c.c && \
    mv tmp.BLAS_zsbmv_d_d_x-f2c.c $@ && rm -f tmp.BLAS_zsbmv_d_d_x-f2c.c*


sbmv.h: ../cblas.m4 ../header.template sbmv-common.m4
	sed -e "s/routine/sbmv/g" -e "s/ROUTINE/SBMV/g" \
  ../header.template | $(M4) $(M4_OPTS) $(M4_PATH) - >$@

source-clean:
	rm -f $(SBMV_SRCS) sbmv.h

.PHONY: sbmv-source sbmv-header source-clean

