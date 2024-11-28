include ../../make.conf
include ../../$(MAKEINC)


SRC_PATH = ../../src/gemm
M4_PATH = -I. -I..

GEMM_SRCS = \
	$(SRC_PATH)/BLAS_sgemm.c\
	$(SRC_PATH)/BLAS_dgemm.c\
	$(SRC_PATH)/BLAS_cgemm.c\
	$(SRC_PATH)/BLAS_zgemm.c\
	$(SRC_PATH)/BLAS_dgemm_d_s.c\
	$(SRC_PATH)/BLAS_dgemm_s_d.c\
	$(SRC_PATH)/BLAS_dgemm_s_s.c\
	$(SRC_PATH)/BLAS_zgemm_z_c.c\
	$(SRC_PATH)/BLAS_zgemm_c_z.c\
	$(SRC_PATH)/BLAS_zgemm_c_c.c\
	$(SRC_PATH)/BLAS_cgemm_c_s.c\
	$(SRC_PATH)/BLAS_cgemm_s_c.c\
	$(SRC_PATH)/BLAS_cgemm_s_s.c\
	$(SRC_PATH)/BLAS_zgemm_z_d.c\
	$(SRC_PATH)/BLAS_zgemm_d_z.c\
	$(SRC_PATH)/BLAS_zgemm_d_d.c\
	$(SRC_PATH)/BLAS_sgemm_x.c\
	$(SRC_PATH)/BLAS_dgemm_x.c\
	$(SRC_PATH)/BLAS_cgemm_x.c\
	$(SRC_PATH)/BLAS_zgemm_x.c\
	$(SRC_PATH)/BLAS_dgemm_d_s_x.c\
	$(SRC_PATH)/BLAS_dgemm_s_d_x.c\
	$(SRC_PATH)/BLAS_dgemm_s_s_x.c\
	$(SRC_PATH)/BLAS_zgemm_z_c_x.c\
	$(SRC_PATH)/BLAS_zgemm_c_z_x.c\
	$(SRC_PATH)/BLAS_zgemm_c_c_x.c\
	$(SRC_PATH)/BLAS_cgemm_c_s_x.c\
	$(SRC_PATH)/BLAS_cgemm_s_c_x.c\
	$(SRC_PATH)/BLAS_cgemm_s_s_x.c\
	$(SRC_PATH)/BLAS_zgemm_z_d_x.c\
	$(SRC_PATH)/BLAS_zgemm_d_z_x.c\
	$(SRC_PATH)/BLAS_zgemm_d_d_x.c \
	\
	$(SRC_PATH)/BLAS_sgemm-f2c.c\
	$(SRC_PATH)/BLAS_dgemm-f2c.c\
	$(SRC_PATH)/BLAS_cgemm-f2c.c\
	$(SRC_PATH)/BLAS_zgemm-f2c.c\
	$(SRC_PATH)/BLAS_dgemm_d_s-f2c.c\
	$(SRC_PATH)/BLAS_dgemm_s_d-f2c.c\
	$(SRC_PATH)/BLAS_dgemm_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zgemm_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zgemm_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zgemm_c_c-f2c.c\
	$(SRC_PATH)/BLAS_cgemm_c_s-f2c.c\
	$(SRC_PATH)/BLAS_cgemm_s_c-f2c.c\
	$(SRC_PATH)/BLAS_cgemm_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zgemm_z_d-f2c.c\
	$(SRC_PATH)/BLAS_zgemm_d_z-f2c.c\
	$(SRC_PATH)/BLAS_zgemm_d_d-f2c.c\
	$(SRC_PATH)/BLAS_sgemm_x-f2c.c\
	$(SRC_PATH)/BLAS_dgemm_x-f2c.c\
	$(SRC_PATH)/BLAS_cgemm_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemm_x-f2c.c\
	$(SRC_PATH)/BLAS_dgemm_d_s_x-f2c.c\
	$(SRC_PATH)/BLAS_dgemm_s_d_x-f2c.c\
	$(SRC_PATH)/BLAS_dgemm_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemm_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemm_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemm_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cgemm_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_cgemm_s_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cgemm_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemm_z_d_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemm_d_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemm_d_d_x-f2c.c

gemm-source: $(GEMM_SRCS)

gemm-header: gemm.h

$(SRC_PATH)/BLAS_sgemm.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(s, s, s)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_sgemm.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sgemm.c && \
    mv tmp.BLAS_sgemm.c $@ && rm -f tmp.BLAS_sgemm.c*

$(SRC_PATH)/BLAS_dgemm.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(d, d, d)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgemm.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemm.c && \
    mv tmp.BLAS_dgemm.c $@ && rm -f tmp.BLAS_dgemm.c*

$(SRC_PATH)/BLAS_cgemm.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(c, c, c)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgemm.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemm.c && \
    mv tmp.BLAS_cgemm.c $@ && rm -f tmp.BLAS_cgemm.c*

$(SRC_PATH)/BLAS_zgemm.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(z, z, z)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemm.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm.c && \
    mv tmp.BLAS_zgemm.c $@ && rm -f tmp.BLAS_zgemm.c*

$(SRC_PATH)/BLAS_dgemm_d_s.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(d, d, s)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgemm_d_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemm_d_s.c && \
    mv tmp.BLAS_dgemm_d_s.c $@ && rm -f tmp.BLAS_dgemm_d_s.c*

$(SRC_PATH)/BLAS_dgemm_s_d.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(d, s, d)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgemm_s_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemm_s_d.c && \
    mv tmp.BLAS_dgemm_s_d.c $@ && rm -f tmp.BLAS_dgemm_s_d.c*

$(SRC_PATH)/BLAS_dgemm_s_s.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(d, s, s)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgemm_s_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemm_s_s.c && \
    mv tmp.BLAS_dgemm_s_s.c $@ && rm -f tmp.BLAS_dgemm_s_s.c*

$(SRC_PATH)/BLAS_zgemm_z_c.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(z, z, c)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemm_z_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm_z_c.c && \
    mv tmp.BLAS_zgemm_z_c.c $@ && rm -f tmp.BLAS_zgemm_z_c.c*

$(SRC_PATH)/BLAS_zgemm_c_z.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(z, c, z)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemm_c_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm_c_z.c && \
    mv tmp.BLAS_zgemm_c_z.c $@ && rm -f tmp.BLAS_zgemm_c_z.c*

$(SRC_PATH)/BLAS_zgemm_c_c.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(z, c, c)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemm_c_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm_c_c.c && \
    mv tmp.BLAS_zgemm_c_c.c $@ && rm -f tmp.BLAS_zgemm_c_c.c*

$(SRC_PATH)/BLAS_cgemm_c_s.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(c, c, s)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgemm_c_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemm_c_s.c && \
    mv tmp.BLAS_cgemm_c_s.c $@ && rm -f tmp.BLAS_cgemm_c_s.c*

$(SRC_PATH)/BLAS_cgemm_s_c.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(c, s, c)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgemm_s_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemm_s_c.c && \
    mv tmp.BLAS_cgemm_s_c.c $@ && rm -f tmp.BLAS_cgemm_s_c.c*

$(SRC_PATH)/BLAS_cgemm_s_s.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(c, s, s)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgemm_s_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemm_s_s.c && \
    mv tmp.BLAS_cgemm_s_s.c $@ && rm -f tmp.BLAS_cgemm_s_s.c*

$(SRC_PATH)/BLAS_zgemm_z_d.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(z, z, d)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemm_z_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm_z_d.c && \
    mv tmp.BLAS_zgemm_z_d.c $@ && rm -f tmp.BLAS_zgemm_z_d.c*

$(SRC_PATH)/BLAS_zgemm_d_z.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(z, d, z)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemm_d_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm_d_z.c && \
    mv tmp.BLAS_zgemm_d_z.c $@ && rm -f tmp.BLAS_zgemm_d_z.c*

$(SRC_PATH)/BLAS_zgemm_d_d.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(z, d, d)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemm_d_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm_d_d.c && \
    mv tmp.BLAS_zgemm_d_d.c $@ && rm -f tmp.BLAS_zgemm_d_d.c*

$(SRC_PATH)/BLAS_sgemm_x.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(s, s, s, _x)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_sgemm_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sgemm_x.c && \
    mv tmp.BLAS_sgemm_x.c $@ && rm -f tmp.BLAS_sgemm_x.c*

$(SRC_PATH)/BLAS_dgemm_x.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(d, d, d, _x)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgemm_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemm_x.c && \
    mv tmp.BLAS_dgemm_x.c $@ && rm -f tmp.BLAS_dgemm_x.c*

$(SRC_PATH)/BLAS_cgemm_x.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(c, c, c, _x)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgemm_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemm_x.c && \
    mv tmp.BLAS_cgemm_x.c $@ && rm -f tmp.BLAS_cgemm_x.c*

$(SRC_PATH)/BLAS_zgemm_x.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(z, z, z, _x)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemm_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm_x.c && \
    mv tmp.BLAS_zgemm_x.c $@ && rm -f tmp.BLAS_zgemm_x.c*

$(SRC_PATH)/BLAS_dgemm_d_s_x.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(d, d, s, _x)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgemm_d_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemm_d_s_x.c && \
    mv tmp.BLAS_dgemm_d_s_x.c $@ && rm -f tmp.BLAS_dgemm_d_s_x.c*

$(SRC_PATH)/BLAS_dgemm_s_d_x.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(d, s, d, _x)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgemm_s_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemm_s_d_x.c && \
    mv tmp.BLAS_dgemm_s_d_x.c $@ && rm -f tmp.BLAS_dgemm_s_d_x.c*

$(SRC_PATH)/BLAS_dgemm_s_s_x.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(d, s, s, _x)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgemm_s_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemm_s_s_x.c && \
    mv tmp.BLAS_dgemm_s_s_x.c $@ && rm -f tmp.BLAS_dgemm_s_s_x.c*

$(SRC_PATH)/BLAS_zgemm_z_c_x.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(z, z, c, _x)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemm_z_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm_z_c_x.c && \
    mv tmp.BLAS_zgemm_z_c_x.c $@ && rm -f tmp.BLAS_zgemm_z_c_x.c*

$(SRC_PATH)/BLAS_zgemm_c_z_x.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(z, c, z, _x)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemm_c_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm_c_z_x.c && \
    mv tmp.BLAS_zgemm_c_z_x.c $@ && rm -f tmp.BLAS_zgemm_c_z_x.c*

$(SRC_PATH)/BLAS_zgemm_c_c_x.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(z, c, c, _x)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemm_c_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm_c_c_x.c && \
    mv tmp.BLAS_zgemm_c_c_x.c $@ && rm -f tmp.BLAS_zgemm_c_c_x.c*

$(SRC_PATH)/BLAS_cgemm_c_s_x.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(c, c, s, _x)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgemm_c_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemm_c_s_x.c && \
    mv tmp.BLAS_cgemm_c_s_x.c $@ && rm -f tmp.BLAS_cgemm_c_s_x.c*

$(SRC_PATH)/BLAS_cgemm_s_c_x.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(c, s, c, _x)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgemm_s_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemm_s_c_x.c && \
    mv tmp.BLAS_cgemm_s_c_x.c $@ && rm -f tmp.BLAS_cgemm_s_c_x.c*

$(SRC_PATH)/BLAS_cgemm_s_s_x.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(c, s, s, _x)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgemm_s_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemm_s_s_x.c && \
    mv tmp.BLAS_cgemm_s_s_x.c $@ && rm -f tmp.BLAS_cgemm_s_s_x.c*

$(SRC_PATH)/BLAS_zgemm_z_d_x.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(z, z, d, _x)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemm_z_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm_z_d_x.c && \
    mv tmp.BLAS_zgemm_z_d_x.c $@ && rm -f tmp.BLAS_zgemm_z_d_x.c*

$(SRC_PATH)/BLAS_zgemm_d_z_x.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(z, d, z, _x)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemm_d_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm_d_z_x.c && \
    mv tmp.BLAS_zgemm_d_z_x.c $@ && rm -f tmp.BLAS_zgemm_d_z_x.c*

$(SRC_PATH)/BLAS_zgemm_d_d_x.c: gemm.m4 gemm-common.m4 ../cblas.m4
	echo "GEMM(z, d, d, _x)" | cat gemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemm_d_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm_d_d_x.c && \
    mv tmp.BLAS_zgemm_d_d_x.c $@ && rm -f tmp.BLAS_zgemm_d_d_x.c*



$(SRC_PATH)/BLAS_sgemm-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="s, s, s" ../f2c-bridge.m4 > tmp.BLAS_sgemm-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sgemm-f2c.c && \
    mv tmp.BLAS_sgemm-f2c.c $@ && rm -f tmp.BLAS_sgemm-f2c.c*

$(SRC_PATH)/BLAS_dgemm-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="d, d, d" ../f2c-bridge.m4 > tmp.BLAS_dgemm-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemm-f2c.c && \
    mv tmp.BLAS_dgemm-f2c.c $@ && rm -f tmp.BLAS_dgemm-f2c.c*

$(SRC_PATH)/BLAS_cgemm-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="c, c, c" ../f2c-bridge.m4 > tmp.BLAS_cgemm-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemm-f2c.c && \
    mv tmp.BLAS_cgemm-f2c.c $@ && rm -f tmp.BLAS_cgemm-f2c.c*

$(SRC_PATH)/BLAS_zgemm-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="z, z, z" ../f2c-bridge.m4 > tmp.BLAS_zgemm-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm-f2c.c && \
    mv tmp.BLAS_zgemm-f2c.c $@ && rm -f tmp.BLAS_zgemm-f2c.c*

$(SRC_PATH)/BLAS_dgemm_d_s-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="d, d, s" ../f2c-bridge.m4 > tmp.BLAS_dgemm_d_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemm_d_s-f2c.c && \
    mv tmp.BLAS_dgemm_d_s-f2c.c $@ && rm -f tmp.BLAS_dgemm_d_s-f2c.c*

$(SRC_PATH)/BLAS_dgemm_s_d-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="d, s, d" ../f2c-bridge.m4 > tmp.BLAS_dgemm_s_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemm_s_d-f2c.c && \
    mv tmp.BLAS_dgemm_s_d-f2c.c $@ && rm -f tmp.BLAS_dgemm_s_d-f2c.c*

$(SRC_PATH)/BLAS_dgemm_s_s-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="d, s, s" ../f2c-bridge.m4 > tmp.BLAS_dgemm_s_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemm_s_s-f2c.c && \
    mv tmp.BLAS_dgemm_s_s-f2c.c $@ && rm -f tmp.BLAS_dgemm_s_s-f2c.c*

$(SRC_PATH)/BLAS_zgemm_z_c-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="z, z, c" ../f2c-bridge.m4 > tmp.BLAS_zgemm_z_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm_z_c-f2c.c && \
    mv tmp.BLAS_zgemm_z_c-f2c.c $@ && rm -f tmp.BLAS_zgemm_z_c-f2c.c*

$(SRC_PATH)/BLAS_zgemm_c_z-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="z, c, z" ../f2c-bridge.m4 > tmp.BLAS_zgemm_c_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm_c_z-f2c.c && \
    mv tmp.BLAS_zgemm_c_z-f2c.c $@ && rm -f tmp.BLAS_zgemm_c_z-f2c.c*

$(SRC_PATH)/BLAS_zgemm_c_c-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="z, c, c" ../f2c-bridge.m4 > tmp.BLAS_zgemm_c_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm_c_c-f2c.c && \
    mv tmp.BLAS_zgemm_c_c-f2c.c $@ && rm -f tmp.BLAS_zgemm_c_c-f2c.c*

$(SRC_PATH)/BLAS_cgemm_c_s-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="c, c, s" ../f2c-bridge.m4 > tmp.BLAS_cgemm_c_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemm_c_s-f2c.c && \
    mv tmp.BLAS_cgemm_c_s-f2c.c $@ && rm -f tmp.BLAS_cgemm_c_s-f2c.c*

$(SRC_PATH)/BLAS_cgemm_s_c-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="c, s, c" ../f2c-bridge.m4 > tmp.BLAS_cgemm_s_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemm_s_c-f2c.c && \
    mv tmp.BLAS_cgemm_s_c-f2c.c $@ && rm -f tmp.BLAS_cgemm_s_c-f2c.c*

$(SRC_PATH)/BLAS_cgemm_s_s-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="c, s, s" ../f2c-bridge.m4 > tmp.BLAS_cgemm_s_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemm_s_s-f2c.c && \
    mv tmp.BLAS_cgemm_s_s-f2c.c $@ && rm -f tmp.BLAS_cgemm_s_s-f2c.c*

$(SRC_PATH)/BLAS_zgemm_z_d-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="z, z, d" ../f2c-bridge.m4 > tmp.BLAS_zgemm_z_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm_z_d-f2c.c && \
    mv tmp.BLAS_zgemm_z_d-f2c.c $@ && rm -f tmp.BLAS_zgemm_z_d-f2c.c*

$(SRC_PATH)/BLAS_zgemm_d_z-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="z, d, z" ../f2c-bridge.m4 > tmp.BLAS_zgemm_d_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm_d_z-f2c.c && \
    mv tmp.BLAS_zgemm_d_z-f2c.c $@ && rm -f tmp.BLAS_zgemm_d_z-f2c.c*

$(SRC_PATH)/BLAS_zgemm_d_d-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="z, d, d" ../f2c-bridge.m4 > tmp.BLAS_zgemm_d_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm_d_d-f2c.c && \
    mv tmp.BLAS_zgemm_d_d-f2c.c $@ && rm -f tmp.BLAS_zgemm_d_d-f2c.c*

$(SRC_PATH)/BLAS_sgemm_x-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="s, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_sgemm_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sgemm_x-f2c.c && \
    mv tmp.BLAS_sgemm_x-f2c.c $@ && rm -f tmp.BLAS_sgemm_x-f2c.c*

$(SRC_PATH)/BLAS_dgemm_x-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="d, d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dgemm_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemm_x-f2c.c && \
    mv tmp.BLAS_dgemm_x-f2c.c $@ && rm -f tmp.BLAS_dgemm_x-f2c.c*

$(SRC_PATH)/BLAS_cgemm_x-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="c, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_cgemm_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemm_x-f2c.c && \
    mv tmp.BLAS_cgemm_x-f2c.c $@ && rm -f tmp.BLAS_cgemm_x-f2c.c*

$(SRC_PATH)/BLAS_zgemm_x-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="z, z, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zgemm_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm_x-f2c.c && \
    mv tmp.BLAS_zgemm_x-f2c.c $@ && rm -f tmp.BLAS_zgemm_x-f2c.c*

$(SRC_PATH)/BLAS_dgemm_d_s_x-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="d, d, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dgemm_d_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemm_d_s_x-f2c.c && \
    mv tmp.BLAS_dgemm_d_s_x-f2c.c $@ && rm -f tmp.BLAS_dgemm_d_s_x-f2c.c*

$(SRC_PATH)/BLAS_dgemm_s_d_x-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="d, s, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dgemm_s_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemm_s_d_x-f2c.c && \
    mv tmp.BLAS_dgemm_s_d_x-f2c.c $@ && rm -f tmp.BLAS_dgemm_s_d_x-f2c.c*

$(SRC_PATH)/BLAS_dgemm_s_s_x-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="d, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dgemm_s_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemm_s_s_x-f2c.c && \
    mv tmp.BLAS_dgemm_s_s_x-f2c.c $@ && rm -f tmp.BLAS_dgemm_s_s_x-f2c.c*

$(SRC_PATH)/BLAS_zgemm_z_c_x-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="z, z, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zgemm_z_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm_z_c_x-f2c.c && \
    mv tmp.BLAS_zgemm_z_c_x-f2c.c $@ && rm -f tmp.BLAS_zgemm_z_c_x-f2c.c*

$(SRC_PATH)/BLAS_zgemm_c_z_x-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="z, c, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zgemm_c_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm_c_z_x-f2c.c && \
    mv tmp.BLAS_zgemm_c_z_x-f2c.c $@ && rm -f tmp.BLAS_zgemm_c_z_x-f2c.c*

$(SRC_PATH)/BLAS_zgemm_c_c_x-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="z, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zgemm_c_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm_c_c_x-f2c.c && \
    mv tmp.BLAS_zgemm_c_c_x-f2c.c $@ && rm -f tmp.BLAS_zgemm_c_c_x-f2c.c*

$(SRC_PATH)/BLAS_cgemm_c_s_x-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="c, c, s, _x" ../f2c-bridge.m4 > tmp.BLAS_cgemm_c_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemm_c_s_x-f2c.c && \
    mv tmp.BLAS_cgemm_c_s_x-f2c.c $@ && rm -f tmp.BLAS_cgemm_c_s_x-f2c.c*

$(SRC_PATH)/BLAS_cgemm_s_c_x-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="c, s, c, _x" ../f2c-bridge.m4 > tmp.BLAS_cgemm_s_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemm_s_c_x-f2c.c && \
    mv tmp.BLAS_cgemm_s_c_x-f2c.c $@ && rm -f tmp.BLAS_cgemm_s_c_x-f2c.c*

$(SRC_PATH)/BLAS_cgemm_s_s_x-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="c, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_cgemm_s_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemm_s_s_x-f2c.c && \
    mv tmp.BLAS_cgemm_s_s_x-f2c.c $@ && rm -f tmp.BLAS_cgemm_s_s_x-f2c.c*

$(SRC_PATH)/BLAS_zgemm_z_d_x-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="z, z, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zgemm_z_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm_z_d_x-f2c.c && \
    mv tmp.BLAS_zgemm_z_d_x-f2c.c $@ && rm -f tmp.BLAS_zgemm_z_d_x-f2c.c*

$(SRC_PATH)/BLAS_zgemm_d_z_x-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="z, d, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zgemm_d_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm_d_z_x-f2c.c && \
    mv tmp.BLAS_zgemm_d_z_x-f2c.c $@ && rm -f tmp.BLAS_zgemm_d_z_x-f2c.c*

$(SRC_PATH)/BLAS_zgemm_d_d_x-f2c.c: ../f2c-bridge.m4 gemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemm -DRARG="z, d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zgemm_d_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemm_d_d_x-f2c.c && \
    mv tmp.BLAS_zgemm_d_d_x-f2c.c $@ && rm -f tmp.BLAS_zgemm_d_d_x-f2c.c*


gemm.h: ../cblas.m4 ../header.template gemm-common.m4
	sed -e "s/routine/gemm/g" -e "s/ROUTINE/GEMM/g" \
  ../header.template | $(M4) $(M4_OPTS) $(M4_PATH) - >$@

source-clean:
	rm -f $(GEMM_SRCS) gemm.h

.PHONY: gemm-source gemm-header source-clean

