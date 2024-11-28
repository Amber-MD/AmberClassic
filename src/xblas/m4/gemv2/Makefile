include ../../make.conf
include ../../$(MAKEINC)


SRC_PATH = ../../src/gemv2
M4_PATH = -I. -I..

GEMV2_SRCS = \
	$(SRC_PATH)/BLAS_sgemv2.c\
	$(SRC_PATH)/BLAS_dgemv2.c\
	$(SRC_PATH)/BLAS_cgemv2.c\
	$(SRC_PATH)/BLAS_zgemv2.c\
	$(SRC_PATH)/BLAS_dgemv2_d_s.c\
	$(SRC_PATH)/BLAS_dgemv2_s_d.c\
	$(SRC_PATH)/BLAS_dgemv2_s_s.c\
	$(SRC_PATH)/BLAS_zgemv2_z_c.c\
	$(SRC_PATH)/BLAS_zgemv2_c_z.c\
	$(SRC_PATH)/BLAS_zgemv2_c_c.c\
	$(SRC_PATH)/BLAS_cgemv2_c_s.c\
	$(SRC_PATH)/BLAS_cgemv2_s_c.c\
	$(SRC_PATH)/BLAS_cgemv2_s_s.c\
	$(SRC_PATH)/BLAS_zgemv2_z_d.c\
	$(SRC_PATH)/BLAS_zgemv2_d_z.c\
	$(SRC_PATH)/BLAS_zgemv2_d_d.c\
	$(SRC_PATH)/BLAS_sgemv2_x.c\
	$(SRC_PATH)/BLAS_dgemv2_x.c\
	$(SRC_PATH)/BLAS_cgemv2_x.c\
	$(SRC_PATH)/BLAS_zgemv2_x.c\
	$(SRC_PATH)/BLAS_dgemv2_d_s_x.c\
	$(SRC_PATH)/BLAS_dgemv2_s_d_x.c\
	$(SRC_PATH)/BLAS_dgemv2_s_s_x.c\
	$(SRC_PATH)/BLAS_zgemv2_z_c_x.c\
	$(SRC_PATH)/BLAS_zgemv2_c_z_x.c\
	$(SRC_PATH)/BLAS_zgemv2_c_c_x.c\
	$(SRC_PATH)/BLAS_cgemv2_c_s_x.c\
	$(SRC_PATH)/BLAS_cgemv2_s_c_x.c\
	$(SRC_PATH)/BLAS_cgemv2_s_s_x.c\
	$(SRC_PATH)/BLAS_zgemv2_z_d_x.c\
	$(SRC_PATH)/BLAS_zgemv2_d_z_x.c\
	$(SRC_PATH)/BLAS_zgemv2_d_d_x.c \
	\
	$(SRC_PATH)/BLAS_sgemv2-f2c.c\
	$(SRC_PATH)/BLAS_dgemv2-f2c.c\
	$(SRC_PATH)/BLAS_cgemv2-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2-f2c.c\
	$(SRC_PATH)/BLAS_dgemv2_d_s-f2c.c\
	$(SRC_PATH)/BLAS_dgemv2_s_d-f2c.c\
	$(SRC_PATH)/BLAS_dgemv2_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2_c_c-f2c.c\
	$(SRC_PATH)/BLAS_cgemv2_c_s-f2c.c\
	$(SRC_PATH)/BLAS_cgemv2_s_c-f2c.c\
	$(SRC_PATH)/BLAS_cgemv2_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2_z_d-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2_d_z-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2_d_d-f2c.c\
	$(SRC_PATH)/BLAS_sgemv2_x-f2c.c\
	$(SRC_PATH)/BLAS_dgemv2_x-f2c.c\
	$(SRC_PATH)/BLAS_cgemv2_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2_x-f2c.c\
	$(SRC_PATH)/BLAS_dgemv2_d_s_x-f2c.c\
	$(SRC_PATH)/BLAS_dgemv2_s_d_x-f2c.c\
	$(SRC_PATH)/BLAS_dgemv2_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cgemv2_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_cgemv2_s_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cgemv2_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2_z_d_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2_d_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2_d_d_x-f2c.c

gemv2-source: $(GEMV2_SRCS)

gemv2-header: gemv2.h

$(SRC_PATH)/BLAS_sgemv2.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(s, s, s)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_sgemv2.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sgemv2.c && \
    mv tmp.BLAS_sgemv2.c $@ && rm -f tmp.BLAS_sgemv2.c*

$(SRC_PATH)/BLAS_dgemv2.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(d, d, d)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgemv2.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv2.c && \
    mv tmp.BLAS_dgemv2.c $@ && rm -f tmp.BLAS_dgemv2.c*

$(SRC_PATH)/BLAS_cgemv2.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(c, c, c)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgemv2.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv2.c && \
    mv tmp.BLAS_cgemv2.c $@ && rm -f tmp.BLAS_cgemv2.c*

$(SRC_PATH)/BLAS_zgemv2.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(z, z, z)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv2.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2.c && \
    mv tmp.BLAS_zgemv2.c $@ && rm -f tmp.BLAS_zgemv2.c*

$(SRC_PATH)/BLAS_dgemv2_d_s.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(d, d, s)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgemv2_d_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv2_d_s.c && \
    mv tmp.BLAS_dgemv2_d_s.c $@ && rm -f tmp.BLAS_dgemv2_d_s.c*

$(SRC_PATH)/BLAS_dgemv2_s_d.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(d, s, d)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgemv2_s_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv2_s_d.c && \
    mv tmp.BLAS_dgemv2_s_d.c $@ && rm -f tmp.BLAS_dgemv2_s_d.c*

$(SRC_PATH)/BLAS_dgemv2_s_s.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(d, s, s)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgemv2_s_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv2_s_s.c && \
    mv tmp.BLAS_dgemv2_s_s.c $@ && rm -f tmp.BLAS_dgemv2_s_s.c*

$(SRC_PATH)/BLAS_zgemv2_z_c.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(z, z, c)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv2_z_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2_z_c.c && \
    mv tmp.BLAS_zgemv2_z_c.c $@ && rm -f tmp.BLAS_zgemv2_z_c.c*

$(SRC_PATH)/BLAS_zgemv2_c_z.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(z, c, z)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv2_c_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2_c_z.c && \
    mv tmp.BLAS_zgemv2_c_z.c $@ && rm -f tmp.BLAS_zgemv2_c_z.c*

$(SRC_PATH)/BLAS_zgemv2_c_c.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(z, c, c)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv2_c_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2_c_c.c && \
    mv tmp.BLAS_zgemv2_c_c.c $@ && rm -f tmp.BLAS_zgemv2_c_c.c*

$(SRC_PATH)/BLAS_cgemv2_c_s.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(c, c, s)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgemv2_c_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv2_c_s.c && \
    mv tmp.BLAS_cgemv2_c_s.c $@ && rm -f tmp.BLAS_cgemv2_c_s.c*

$(SRC_PATH)/BLAS_cgemv2_s_c.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(c, s, c)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgemv2_s_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv2_s_c.c && \
    mv tmp.BLAS_cgemv2_s_c.c $@ && rm -f tmp.BLAS_cgemv2_s_c.c*

$(SRC_PATH)/BLAS_cgemv2_s_s.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(c, s, s)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgemv2_s_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv2_s_s.c && \
    mv tmp.BLAS_cgemv2_s_s.c $@ && rm -f tmp.BLAS_cgemv2_s_s.c*

$(SRC_PATH)/BLAS_zgemv2_z_d.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(z, z, d)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv2_z_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2_z_d.c && \
    mv tmp.BLAS_zgemv2_z_d.c $@ && rm -f tmp.BLAS_zgemv2_z_d.c*

$(SRC_PATH)/BLAS_zgemv2_d_z.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(z, d, z)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv2_d_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2_d_z.c && \
    mv tmp.BLAS_zgemv2_d_z.c $@ && rm -f tmp.BLAS_zgemv2_d_z.c*

$(SRC_PATH)/BLAS_zgemv2_d_d.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(z, d, d)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv2_d_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2_d_d.c && \
    mv tmp.BLAS_zgemv2_d_d.c $@ && rm -f tmp.BLAS_zgemv2_d_d.c*

$(SRC_PATH)/BLAS_sgemv2_x.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(s, s, s, _x)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_sgemv2_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sgemv2_x.c && \
    mv tmp.BLAS_sgemv2_x.c $@ && rm -f tmp.BLAS_sgemv2_x.c*

$(SRC_PATH)/BLAS_dgemv2_x.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(d, d, d, _x)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgemv2_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv2_x.c && \
    mv tmp.BLAS_dgemv2_x.c $@ && rm -f tmp.BLAS_dgemv2_x.c*

$(SRC_PATH)/BLAS_cgemv2_x.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(c, c, c, _x)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgemv2_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv2_x.c && \
    mv tmp.BLAS_cgemv2_x.c $@ && rm -f tmp.BLAS_cgemv2_x.c*

$(SRC_PATH)/BLAS_zgemv2_x.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(z, z, z, _x)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv2_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2_x.c && \
    mv tmp.BLAS_zgemv2_x.c $@ && rm -f tmp.BLAS_zgemv2_x.c*

$(SRC_PATH)/BLAS_dgemv2_d_s_x.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(d, d, s, _x)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgemv2_d_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv2_d_s_x.c && \
    mv tmp.BLAS_dgemv2_d_s_x.c $@ && rm -f tmp.BLAS_dgemv2_d_s_x.c*

$(SRC_PATH)/BLAS_dgemv2_s_d_x.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(d, s, d, _x)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgemv2_s_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv2_s_d_x.c && \
    mv tmp.BLAS_dgemv2_s_d_x.c $@ && rm -f tmp.BLAS_dgemv2_s_d_x.c*

$(SRC_PATH)/BLAS_dgemv2_s_s_x.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(d, s, s, _x)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgemv2_s_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv2_s_s_x.c && \
    mv tmp.BLAS_dgemv2_s_s_x.c $@ && rm -f tmp.BLAS_dgemv2_s_s_x.c*

$(SRC_PATH)/BLAS_zgemv2_z_c_x.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(z, z, c, _x)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv2_z_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2_z_c_x.c && \
    mv tmp.BLAS_zgemv2_z_c_x.c $@ && rm -f tmp.BLAS_zgemv2_z_c_x.c*

$(SRC_PATH)/BLAS_zgemv2_c_z_x.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(z, c, z, _x)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv2_c_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2_c_z_x.c && \
    mv tmp.BLAS_zgemv2_c_z_x.c $@ && rm -f tmp.BLAS_zgemv2_c_z_x.c*

$(SRC_PATH)/BLAS_zgemv2_c_c_x.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(z, c, c, _x)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv2_c_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2_c_c_x.c && \
    mv tmp.BLAS_zgemv2_c_c_x.c $@ && rm -f tmp.BLAS_zgemv2_c_c_x.c*

$(SRC_PATH)/BLAS_cgemv2_c_s_x.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(c, c, s, _x)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgemv2_c_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv2_c_s_x.c && \
    mv tmp.BLAS_cgemv2_c_s_x.c $@ && rm -f tmp.BLAS_cgemv2_c_s_x.c*

$(SRC_PATH)/BLAS_cgemv2_s_c_x.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(c, s, c, _x)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgemv2_s_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv2_s_c_x.c && \
    mv tmp.BLAS_cgemv2_s_c_x.c $@ && rm -f tmp.BLAS_cgemv2_s_c_x.c*

$(SRC_PATH)/BLAS_cgemv2_s_s_x.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(c, s, s, _x)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgemv2_s_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv2_s_s_x.c && \
    mv tmp.BLAS_cgemv2_s_s_x.c $@ && rm -f tmp.BLAS_cgemv2_s_s_x.c*

$(SRC_PATH)/BLAS_zgemv2_z_d_x.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(z, z, d, _x)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv2_z_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2_z_d_x.c && \
    mv tmp.BLAS_zgemv2_z_d_x.c $@ && rm -f tmp.BLAS_zgemv2_z_d_x.c*

$(SRC_PATH)/BLAS_zgemv2_d_z_x.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(z, d, z, _x)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv2_d_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2_d_z_x.c && \
    mv tmp.BLAS_zgemv2_d_z_x.c $@ && rm -f tmp.BLAS_zgemv2_d_z_x.c*

$(SRC_PATH)/BLAS_zgemv2_d_d_x.c: gemv2.m4 gemv2-common.m4 ../cblas.m4
	echo "GEMV2(z, d, d, _x)" | cat gemv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv2_d_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2_d_d_x.c && \
    mv tmp.BLAS_zgemv2_d_d_x.c $@ && rm -f tmp.BLAS_zgemv2_d_d_x.c*



$(SRC_PATH)/BLAS_sgemv2-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="s, s, s" ../f2c-bridge.m4 > tmp.BLAS_sgemv2-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sgemv2-f2c.c && \
    mv tmp.BLAS_sgemv2-f2c.c $@ && rm -f tmp.BLAS_sgemv2-f2c.c*

$(SRC_PATH)/BLAS_dgemv2-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="d, d, d" ../f2c-bridge.m4 > tmp.BLAS_dgemv2-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv2-f2c.c && \
    mv tmp.BLAS_dgemv2-f2c.c $@ && rm -f tmp.BLAS_dgemv2-f2c.c*

$(SRC_PATH)/BLAS_cgemv2-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="c, c, c" ../f2c-bridge.m4 > tmp.BLAS_cgemv2-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv2-f2c.c && \
    mv tmp.BLAS_cgemv2-f2c.c $@ && rm -f tmp.BLAS_cgemv2-f2c.c*

$(SRC_PATH)/BLAS_zgemv2-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="z, z, z" ../f2c-bridge.m4 > tmp.BLAS_zgemv2-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2-f2c.c && \
    mv tmp.BLAS_zgemv2-f2c.c $@ && rm -f tmp.BLAS_zgemv2-f2c.c*

$(SRC_PATH)/BLAS_dgemv2_d_s-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="d, d, s" ../f2c-bridge.m4 > tmp.BLAS_dgemv2_d_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv2_d_s-f2c.c && \
    mv tmp.BLAS_dgemv2_d_s-f2c.c $@ && rm -f tmp.BLAS_dgemv2_d_s-f2c.c*

$(SRC_PATH)/BLAS_dgemv2_s_d-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="d, s, d" ../f2c-bridge.m4 > tmp.BLAS_dgemv2_s_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv2_s_d-f2c.c && \
    mv tmp.BLAS_dgemv2_s_d-f2c.c $@ && rm -f tmp.BLAS_dgemv2_s_d-f2c.c*

$(SRC_PATH)/BLAS_dgemv2_s_s-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="d, s, s" ../f2c-bridge.m4 > tmp.BLAS_dgemv2_s_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv2_s_s-f2c.c && \
    mv tmp.BLAS_dgemv2_s_s-f2c.c $@ && rm -f tmp.BLAS_dgemv2_s_s-f2c.c*

$(SRC_PATH)/BLAS_zgemv2_z_c-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="z, z, c" ../f2c-bridge.m4 > tmp.BLAS_zgemv2_z_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2_z_c-f2c.c && \
    mv tmp.BLAS_zgemv2_z_c-f2c.c $@ && rm -f tmp.BLAS_zgemv2_z_c-f2c.c*

$(SRC_PATH)/BLAS_zgemv2_c_z-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="z, c, z" ../f2c-bridge.m4 > tmp.BLAS_zgemv2_c_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2_c_z-f2c.c && \
    mv tmp.BLAS_zgemv2_c_z-f2c.c $@ && rm -f tmp.BLAS_zgemv2_c_z-f2c.c*

$(SRC_PATH)/BLAS_zgemv2_c_c-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="z, c, c" ../f2c-bridge.m4 > tmp.BLAS_zgemv2_c_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2_c_c-f2c.c && \
    mv tmp.BLAS_zgemv2_c_c-f2c.c $@ && rm -f tmp.BLAS_zgemv2_c_c-f2c.c*

$(SRC_PATH)/BLAS_cgemv2_c_s-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="c, c, s" ../f2c-bridge.m4 > tmp.BLAS_cgemv2_c_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv2_c_s-f2c.c && \
    mv tmp.BLAS_cgemv2_c_s-f2c.c $@ && rm -f tmp.BLAS_cgemv2_c_s-f2c.c*

$(SRC_PATH)/BLAS_cgemv2_s_c-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="c, s, c" ../f2c-bridge.m4 > tmp.BLAS_cgemv2_s_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv2_s_c-f2c.c && \
    mv tmp.BLAS_cgemv2_s_c-f2c.c $@ && rm -f tmp.BLAS_cgemv2_s_c-f2c.c*

$(SRC_PATH)/BLAS_cgemv2_s_s-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="c, s, s" ../f2c-bridge.m4 > tmp.BLAS_cgemv2_s_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv2_s_s-f2c.c && \
    mv tmp.BLAS_cgemv2_s_s-f2c.c $@ && rm -f tmp.BLAS_cgemv2_s_s-f2c.c*

$(SRC_PATH)/BLAS_zgemv2_z_d-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="z, z, d" ../f2c-bridge.m4 > tmp.BLAS_zgemv2_z_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2_z_d-f2c.c && \
    mv tmp.BLAS_zgemv2_z_d-f2c.c $@ && rm -f tmp.BLAS_zgemv2_z_d-f2c.c*

$(SRC_PATH)/BLAS_zgemv2_d_z-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="z, d, z" ../f2c-bridge.m4 > tmp.BLAS_zgemv2_d_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2_d_z-f2c.c && \
    mv tmp.BLAS_zgemv2_d_z-f2c.c $@ && rm -f tmp.BLAS_zgemv2_d_z-f2c.c*

$(SRC_PATH)/BLAS_zgemv2_d_d-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="z, d, d" ../f2c-bridge.m4 > tmp.BLAS_zgemv2_d_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2_d_d-f2c.c && \
    mv tmp.BLAS_zgemv2_d_d-f2c.c $@ && rm -f tmp.BLAS_zgemv2_d_d-f2c.c*

$(SRC_PATH)/BLAS_sgemv2_x-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="s, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_sgemv2_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sgemv2_x-f2c.c && \
    mv tmp.BLAS_sgemv2_x-f2c.c $@ && rm -f tmp.BLAS_sgemv2_x-f2c.c*

$(SRC_PATH)/BLAS_dgemv2_x-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="d, d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dgemv2_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv2_x-f2c.c && \
    mv tmp.BLAS_dgemv2_x-f2c.c $@ && rm -f tmp.BLAS_dgemv2_x-f2c.c*

$(SRC_PATH)/BLAS_cgemv2_x-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="c, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_cgemv2_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv2_x-f2c.c && \
    mv tmp.BLAS_cgemv2_x-f2c.c $@ && rm -f tmp.BLAS_cgemv2_x-f2c.c*

$(SRC_PATH)/BLAS_zgemv2_x-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="z, z, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zgemv2_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2_x-f2c.c && \
    mv tmp.BLAS_zgemv2_x-f2c.c $@ && rm -f tmp.BLAS_zgemv2_x-f2c.c*

$(SRC_PATH)/BLAS_dgemv2_d_s_x-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="d, d, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dgemv2_d_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv2_d_s_x-f2c.c && \
    mv tmp.BLAS_dgemv2_d_s_x-f2c.c $@ && rm -f tmp.BLAS_dgemv2_d_s_x-f2c.c*

$(SRC_PATH)/BLAS_dgemv2_s_d_x-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="d, s, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dgemv2_s_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv2_s_d_x-f2c.c && \
    mv tmp.BLAS_dgemv2_s_d_x-f2c.c $@ && rm -f tmp.BLAS_dgemv2_s_d_x-f2c.c*

$(SRC_PATH)/BLAS_dgemv2_s_s_x-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="d, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dgemv2_s_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv2_s_s_x-f2c.c && \
    mv tmp.BLAS_dgemv2_s_s_x-f2c.c $@ && rm -f tmp.BLAS_dgemv2_s_s_x-f2c.c*

$(SRC_PATH)/BLAS_zgemv2_z_c_x-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="z, z, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zgemv2_z_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2_z_c_x-f2c.c && \
    mv tmp.BLAS_zgemv2_z_c_x-f2c.c $@ && rm -f tmp.BLAS_zgemv2_z_c_x-f2c.c*

$(SRC_PATH)/BLAS_zgemv2_c_z_x-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="z, c, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zgemv2_c_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2_c_z_x-f2c.c && \
    mv tmp.BLAS_zgemv2_c_z_x-f2c.c $@ && rm -f tmp.BLAS_zgemv2_c_z_x-f2c.c*

$(SRC_PATH)/BLAS_zgemv2_c_c_x-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="z, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zgemv2_c_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2_c_c_x-f2c.c && \
    mv tmp.BLAS_zgemv2_c_c_x-f2c.c $@ && rm -f tmp.BLAS_zgemv2_c_c_x-f2c.c*

$(SRC_PATH)/BLAS_cgemv2_c_s_x-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="c, c, s, _x" ../f2c-bridge.m4 > tmp.BLAS_cgemv2_c_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv2_c_s_x-f2c.c && \
    mv tmp.BLAS_cgemv2_c_s_x-f2c.c $@ && rm -f tmp.BLAS_cgemv2_c_s_x-f2c.c*

$(SRC_PATH)/BLAS_cgemv2_s_c_x-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="c, s, c, _x" ../f2c-bridge.m4 > tmp.BLAS_cgemv2_s_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv2_s_c_x-f2c.c && \
    mv tmp.BLAS_cgemv2_s_c_x-f2c.c $@ && rm -f tmp.BLAS_cgemv2_s_c_x-f2c.c*

$(SRC_PATH)/BLAS_cgemv2_s_s_x-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="c, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_cgemv2_s_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv2_s_s_x-f2c.c && \
    mv tmp.BLAS_cgemv2_s_s_x-f2c.c $@ && rm -f tmp.BLAS_cgemv2_s_s_x-f2c.c*

$(SRC_PATH)/BLAS_zgemv2_z_d_x-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="z, z, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zgemv2_z_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2_z_d_x-f2c.c && \
    mv tmp.BLAS_zgemv2_z_d_x-f2c.c $@ && rm -f tmp.BLAS_zgemv2_z_d_x-f2c.c*

$(SRC_PATH)/BLAS_zgemv2_d_z_x-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="z, d, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zgemv2_d_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2_d_z_x-f2c.c && \
    mv tmp.BLAS_zgemv2_d_z_x-f2c.c $@ && rm -f tmp.BLAS_zgemv2_d_z_x-f2c.c*

$(SRC_PATH)/BLAS_zgemv2_d_d_x-f2c.c: ../f2c-bridge.m4 gemv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv2 -DRARG="z, d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zgemv2_d_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv2_d_d_x-f2c.c && \
    mv tmp.BLAS_zgemv2_d_d_x-f2c.c $@ && rm -f tmp.BLAS_zgemv2_d_d_x-f2c.c*


gemv2.h: ../cblas.m4 ../header.template gemv2-common.m4
	sed -e "s/routine/gemv2/g" -e "s/ROUTINE/GEMV2/g" \
  ../header.template | $(M4) $(M4_OPTS) $(M4_PATH) - >$@

source-clean:
	rm -f $(GEMV2_SRCS) gemv2.h

.PHONY: gemv2-source gemv2-header source-clean

