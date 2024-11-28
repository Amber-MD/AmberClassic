include ../../make.conf
include ../../$(MAKEINC)


SRC_PATH = ../../src/gemv
M4_PATH = -I. -I..

GEMV_SRCS = \
	$(SRC_PATH)/BLAS_sgemv.c\
	$(SRC_PATH)/BLAS_dgemv.c\
	$(SRC_PATH)/BLAS_cgemv.c\
	$(SRC_PATH)/BLAS_zgemv.c\
	$(SRC_PATH)/BLAS_dgemv_d_s.c\
	$(SRC_PATH)/BLAS_dgemv_s_d.c\
	$(SRC_PATH)/BLAS_dgemv_s_s.c\
	$(SRC_PATH)/BLAS_zgemv_z_c.c\
	$(SRC_PATH)/BLAS_zgemv_c_z.c\
	$(SRC_PATH)/BLAS_zgemv_c_c.c\
	$(SRC_PATH)/BLAS_cgemv_c_s.c\
	$(SRC_PATH)/BLAS_cgemv_s_c.c\
	$(SRC_PATH)/BLAS_cgemv_s_s.c\
	$(SRC_PATH)/BLAS_zgemv_z_d.c\
	$(SRC_PATH)/BLAS_zgemv_d_z.c\
	$(SRC_PATH)/BLAS_zgemv_d_d.c\
	$(SRC_PATH)/BLAS_sgemv_x.c\
	$(SRC_PATH)/BLAS_dgemv_x.c\
	$(SRC_PATH)/BLAS_cgemv_x.c\
	$(SRC_PATH)/BLAS_zgemv_x.c\
	$(SRC_PATH)/BLAS_dgemv_d_s_x.c\
	$(SRC_PATH)/BLAS_dgemv_s_d_x.c\
	$(SRC_PATH)/BLAS_dgemv_s_s_x.c\
	$(SRC_PATH)/BLAS_zgemv_z_c_x.c\
	$(SRC_PATH)/BLAS_zgemv_c_z_x.c\
	$(SRC_PATH)/BLAS_zgemv_c_c_x.c\
	$(SRC_PATH)/BLAS_cgemv_c_s_x.c\
	$(SRC_PATH)/BLAS_cgemv_s_c_x.c\
	$(SRC_PATH)/BLAS_cgemv_s_s_x.c\
	$(SRC_PATH)/BLAS_zgemv_z_d_x.c\
	$(SRC_PATH)/BLAS_zgemv_d_z_x.c\
	$(SRC_PATH)/BLAS_zgemv_d_d_x.c \
	\
	$(SRC_PATH)/BLAS_sgemv-f2c.c\
	$(SRC_PATH)/BLAS_dgemv-f2c.c\
	$(SRC_PATH)/BLAS_cgemv-f2c.c\
	$(SRC_PATH)/BLAS_zgemv-f2c.c\
	$(SRC_PATH)/BLAS_dgemv_d_s-f2c.c\
	$(SRC_PATH)/BLAS_dgemv_s_d-f2c.c\
	$(SRC_PATH)/BLAS_dgemv_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zgemv_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zgemv_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zgemv_c_c-f2c.c\
	$(SRC_PATH)/BLAS_cgemv_c_s-f2c.c\
	$(SRC_PATH)/BLAS_cgemv_s_c-f2c.c\
	$(SRC_PATH)/BLAS_cgemv_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zgemv_z_d-f2c.c\
	$(SRC_PATH)/BLAS_zgemv_d_z-f2c.c\
	$(SRC_PATH)/BLAS_zgemv_d_d-f2c.c\
	$(SRC_PATH)/BLAS_sgemv_x-f2c.c\
	$(SRC_PATH)/BLAS_dgemv_x-f2c.c\
	$(SRC_PATH)/BLAS_cgemv_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv_x-f2c.c\
	$(SRC_PATH)/BLAS_dgemv_d_s_x-f2c.c\
	$(SRC_PATH)/BLAS_dgemv_s_d_x-f2c.c\
	$(SRC_PATH)/BLAS_dgemv_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cgemv_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_cgemv_s_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cgemv_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv_z_d_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv_d_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv_d_d_x-f2c.c

gemv-source: $(GEMV_SRCS)

gemv-header: gemv.h

$(SRC_PATH)/BLAS_sgemv.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(s, s, s)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_sgemv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sgemv.c && \
    mv tmp.BLAS_sgemv.c $@ && rm -f tmp.BLAS_sgemv.c*

$(SRC_PATH)/BLAS_dgemv.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(d, d, d)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgemv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv.c && \
    mv tmp.BLAS_dgemv.c $@ && rm -f tmp.BLAS_dgemv.c*

$(SRC_PATH)/BLAS_cgemv.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(c, c, c)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgemv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv.c && \
    mv tmp.BLAS_cgemv.c $@ && rm -f tmp.BLAS_cgemv.c*

$(SRC_PATH)/BLAS_zgemv.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(z, z, z)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv.c && \
    mv tmp.BLAS_zgemv.c $@ && rm -f tmp.BLAS_zgemv.c*

$(SRC_PATH)/BLAS_dgemv_d_s.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(d, d, s)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgemv_d_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv_d_s.c && \
    mv tmp.BLAS_dgemv_d_s.c $@ && rm -f tmp.BLAS_dgemv_d_s.c*

$(SRC_PATH)/BLAS_dgemv_s_d.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(d, s, d)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgemv_s_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv_s_d.c && \
    mv tmp.BLAS_dgemv_s_d.c $@ && rm -f tmp.BLAS_dgemv_s_d.c*

$(SRC_PATH)/BLAS_dgemv_s_s.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(d, s, s)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgemv_s_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv_s_s.c && \
    mv tmp.BLAS_dgemv_s_s.c $@ && rm -f tmp.BLAS_dgemv_s_s.c*

$(SRC_PATH)/BLAS_zgemv_z_c.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(z, z, c)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv_z_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv_z_c.c && \
    mv tmp.BLAS_zgemv_z_c.c $@ && rm -f tmp.BLAS_zgemv_z_c.c*

$(SRC_PATH)/BLAS_zgemv_c_z.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(z, c, z)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv_c_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv_c_z.c && \
    mv tmp.BLAS_zgemv_c_z.c $@ && rm -f tmp.BLAS_zgemv_c_z.c*

$(SRC_PATH)/BLAS_zgemv_c_c.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(z, c, c)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv_c_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv_c_c.c && \
    mv tmp.BLAS_zgemv_c_c.c $@ && rm -f tmp.BLAS_zgemv_c_c.c*

$(SRC_PATH)/BLAS_cgemv_c_s.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(c, c, s)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgemv_c_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv_c_s.c && \
    mv tmp.BLAS_cgemv_c_s.c $@ && rm -f tmp.BLAS_cgemv_c_s.c*

$(SRC_PATH)/BLAS_cgemv_s_c.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(c, s, c)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgemv_s_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv_s_c.c && \
    mv tmp.BLAS_cgemv_s_c.c $@ && rm -f tmp.BLAS_cgemv_s_c.c*

$(SRC_PATH)/BLAS_cgemv_s_s.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(c, s, s)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgemv_s_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv_s_s.c && \
    mv tmp.BLAS_cgemv_s_s.c $@ && rm -f tmp.BLAS_cgemv_s_s.c*

$(SRC_PATH)/BLAS_zgemv_z_d.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(z, z, d)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv_z_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv_z_d.c && \
    mv tmp.BLAS_zgemv_z_d.c $@ && rm -f tmp.BLAS_zgemv_z_d.c*

$(SRC_PATH)/BLAS_zgemv_d_z.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(z, d, z)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv_d_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv_d_z.c && \
    mv tmp.BLAS_zgemv_d_z.c $@ && rm -f tmp.BLAS_zgemv_d_z.c*

$(SRC_PATH)/BLAS_zgemv_d_d.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(z, d, d)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv_d_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv_d_d.c && \
    mv tmp.BLAS_zgemv_d_d.c $@ && rm -f tmp.BLAS_zgemv_d_d.c*

$(SRC_PATH)/BLAS_sgemv_x.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(s, s, s, _x)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_sgemv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sgemv_x.c && \
    mv tmp.BLAS_sgemv_x.c $@ && rm -f tmp.BLAS_sgemv_x.c*

$(SRC_PATH)/BLAS_dgemv_x.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(d, d, d, _x)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgemv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv_x.c && \
    mv tmp.BLAS_dgemv_x.c $@ && rm -f tmp.BLAS_dgemv_x.c*

$(SRC_PATH)/BLAS_cgemv_x.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(c, c, c, _x)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgemv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv_x.c && \
    mv tmp.BLAS_cgemv_x.c $@ && rm -f tmp.BLAS_cgemv_x.c*

$(SRC_PATH)/BLAS_zgemv_x.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(z, z, z, _x)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv_x.c && \
    mv tmp.BLAS_zgemv_x.c $@ && rm -f tmp.BLAS_zgemv_x.c*

$(SRC_PATH)/BLAS_dgemv_d_s_x.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(d, d, s, _x)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgemv_d_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv_d_s_x.c && \
    mv tmp.BLAS_dgemv_d_s_x.c $@ && rm -f tmp.BLAS_dgemv_d_s_x.c*

$(SRC_PATH)/BLAS_dgemv_s_d_x.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(d, s, d, _x)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgemv_s_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv_s_d_x.c && \
    mv tmp.BLAS_dgemv_s_d_x.c $@ && rm -f tmp.BLAS_dgemv_s_d_x.c*

$(SRC_PATH)/BLAS_dgemv_s_s_x.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(d, s, s, _x)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dgemv_s_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv_s_s_x.c && \
    mv tmp.BLAS_dgemv_s_s_x.c $@ && rm -f tmp.BLAS_dgemv_s_s_x.c*

$(SRC_PATH)/BLAS_zgemv_z_c_x.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(z, z, c, _x)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv_z_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv_z_c_x.c && \
    mv tmp.BLAS_zgemv_z_c_x.c $@ && rm -f tmp.BLAS_zgemv_z_c_x.c*

$(SRC_PATH)/BLAS_zgemv_c_z_x.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(z, c, z, _x)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv_c_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv_c_z_x.c && \
    mv tmp.BLAS_zgemv_c_z_x.c $@ && rm -f tmp.BLAS_zgemv_c_z_x.c*

$(SRC_PATH)/BLAS_zgemv_c_c_x.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(z, c, c, _x)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv_c_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv_c_c_x.c && \
    mv tmp.BLAS_zgemv_c_c_x.c $@ && rm -f tmp.BLAS_zgemv_c_c_x.c*

$(SRC_PATH)/BLAS_cgemv_c_s_x.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(c, c, s, _x)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgemv_c_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv_c_s_x.c && \
    mv tmp.BLAS_cgemv_c_s_x.c $@ && rm -f tmp.BLAS_cgemv_c_s_x.c*

$(SRC_PATH)/BLAS_cgemv_s_c_x.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(c, s, c, _x)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgemv_s_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv_s_c_x.c && \
    mv tmp.BLAS_cgemv_s_c_x.c $@ && rm -f tmp.BLAS_cgemv_s_c_x.c*

$(SRC_PATH)/BLAS_cgemv_s_s_x.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(c, s, s, _x)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cgemv_s_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv_s_s_x.c && \
    mv tmp.BLAS_cgemv_s_s_x.c $@ && rm -f tmp.BLAS_cgemv_s_s_x.c*

$(SRC_PATH)/BLAS_zgemv_z_d_x.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(z, z, d, _x)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv_z_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv_z_d_x.c && \
    mv tmp.BLAS_zgemv_z_d_x.c $@ && rm -f tmp.BLAS_zgemv_z_d_x.c*

$(SRC_PATH)/BLAS_zgemv_d_z_x.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(z, d, z, _x)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv_d_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv_d_z_x.c && \
    mv tmp.BLAS_zgemv_d_z_x.c $@ && rm -f tmp.BLAS_zgemv_d_z_x.c*

$(SRC_PATH)/BLAS_zgemv_d_d_x.c: gemv.m4 gemv-common.m4 ../cblas.m4
	echo "GEMV(z, d, d, _x)" | cat gemv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zgemv_d_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv_d_d_x.c && \
    mv tmp.BLAS_zgemv_d_d_x.c $@ && rm -f tmp.BLAS_zgemv_d_d_x.c*



$(SRC_PATH)/BLAS_sgemv-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="s, s, s" ../f2c-bridge.m4 > tmp.BLAS_sgemv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sgemv-f2c.c && \
    mv tmp.BLAS_sgemv-f2c.c $@ && rm -f tmp.BLAS_sgemv-f2c.c*

$(SRC_PATH)/BLAS_dgemv-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="d, d, d" ../f2c-bridge.m4 > tmp.BLAS_dgemv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv-f2c.c && \
    mv tmp.BLAS_dgemv-f2c.c $@ && rm -f tmp.BLAS_dgemv-f2c.c*

$(SRC_PATH)/BLAS_cgemv-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="c, c, c" ../f2c-bridge.m4 > tmp.BLAS_cgemv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv-f2c.c && \
    mv tmp.BLAS_cgemv-f2c.c $@ && rm -f tmp.BLAS_cgemv-f2c.c*

$(SRC_PATH)/BLAS_zgemv-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="z, z, z" ../f2c-bridge.m4 > tmp.BLAS_zgemv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv-f2c.c && \
    mv tmp.BLAS_zgemv-f2c.c $@ && rm -f tmp.BLAS_zgemv-f2c.c*

$(SRC_PATH)/BLAS_dgemv_d_s-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="d, d, s" ../f2c-bridge.m4 > tmp.BLAS_dgemv_d_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv_d_s-f2c.c && \
    mv tmp.BLAS_dgemv_d_s-f2c.c $@ && rm -f tmp.BLAS_dgemv_d_s-f2c.c*

$(SRC_PATH)/BLAS_dgemv_s_d-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="d, s, d" ../f2c-bridge.m4 > tmp.BLAS_dgemv_s_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv_s_d-f2c.c && \
    mv tmp.BLAS_dgemv_s_d-f2c.c $@ && rm -f tmp.BLAS_dgemv_s_d-f2c.c*

$(SRC_PATH)/BLAS_dgemv_s_s-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="d, s, s" ../f2c-bridge.m4 > tmp.BLAS_dgemv_s_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv_s_s-f2c.c && \
    mv tmp.BLAS_dgemv_s_s-f2c.c $@ && rm -f tmp.BLAS_dgemv_s_s-f2c.c*

$(SRC_PATH)/BLAS_zgemv_z_c-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="z, z, c" ../f2c-bridge.m4 > tmp.BLAS_zgemv_z_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv_z_c-f2c.c && \
    mv tmp.BLAS_zgemv_z_c-f2c.c $@ && rm -f tmp.BLAS_zgemv_z_c-f2c.c*

$(SRC_PATH)/BLAS_zgemv_c_z-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="z, c, z" ../f2c-bridge.m4 > tmp.BLAS_zgemv_c_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv_c_z-f2c.c && \
    mv tmp.BLAS_zgemv_c_z-f2c.c $@ && rm -f tmp.BLAS_zgemv_c_z-f2c.c*

$(SRC_PATH)/BLAS_zgemv_c_c-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="z, c, c" ../f2c-bridge.m4 > tmp.BLAS_zgemv_c_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv_c_c-f2c.c && \
    mv tmp.BLAS_zgemv_c_c-f2c.c $@ && rm -f tmp.BLAS_zgemv_c_c-f2c.c*

$(SRC_PATH)/BLAS_cgemv_c_s-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="c, c, s" ../f2c-bridge.m4 > tmp.BLAS_cgemv_c_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv_c_s-f2c.c && \
    mv tmp.BLAS_cgemv_c_s-f2c.c $@ && rm -f tmp.BLAS_cgemv_c_s-f2c.c*

$(SRC_PATH)/BLAS_cgemv_s_c-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="c, s, c" ../f2c-bridge.m4 > tmp.BLAS_cgemv_s_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv_s_c-f2c.c && \
    mv tmp.BLAS_cgemv_s_c-f2c.c $@ && rm -f tmp.BLAS_cgemv_s_c-f2c.c*

$(SRC_PATH)/BLAS_cgemv_s_s-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="c, s, s" ../f2c-bridge.m4 > tmp.BLAS_cgemv_s_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv_s_s-f2c.c && \
    mv tmp.BLAS_cgemv_s_s-f2c.c $@ && rm -f tmp.BLAS_cgemv_s_s-f2c.c*

$(SRC_PATH)/BLAS_zgemv_z_d-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="z, z, d" ../f2c-bridge.m4 > tmp.BLAS_zgemv_z_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv_z_d-f2c.c && \
    mv tmp.BLAS_zgemv_z_d-f2c.c $@ && rm -f tmp.BLAS_zgemv_z_d-f2c.c*

$(SRC_PATH)/BLAS_zgemv_d_z-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="z, d, z" ../f2c-bridge.m4 > tmp.BLAS_zgemv_d_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv_d_z-f2c.c && \
    mv tmp.BLAS_zgemv_d_z-f2c.c $@ && rm -f tmp.BLAS_zgemv_d_z-f2c.c*

$(SRC_PATH)/BLAS_zgemv_d_d-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="z, d, d" ../f2c-bridge.m4 > tmp.BLAS_zgemv_d_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv_d_d-f2c.c && \
    mv tmp.BLAS_zgemv_d_d-f2c.c $@ && rm -f tmp.BLAS_zgemv_d_d-f2c.c*

$(SRC_PATH)/BLAS_sgemv_x-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="s, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_sgemv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sgemv_x-f2c.c && \
    mv tmp.BLAS_sgemv_x-f2c.c $@ && rm -f tmp.BLAS_sgemv_x-f2c.c*

$(SRC_PATH)/BLAS_dgemv_x-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="d, d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dgemv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv_x-f2c.c && \
    mv tmp.BLAS_dgemv_x-f2c.c $@ && rm -f tmp.BLAS_dgemv_x-f2c.c*

$(SRC_PATH)/BLAS_cgemv_x-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="c, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_cgemv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv_x-f2c.c && \
    mv tmp.BLAS_cgemv_x-f2c.c $@ && rm -f tmp.BLAS_cgemv_x-f2c.c*

$(SRC_PATH)/BLAS_zgemv_x-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="z, z, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zgemv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv_x-f2c.c && \
    mv tmp.BLAS_zgemv_x-f2c.c $@ && rm -f tmp.BLAS_zgemv_x-f2c.c*

$(SRC_PATH)/BLAS_dgemv_d_s_x-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="d, d, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dgemv_d_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv_d_s_x-f2c.c && \
    mv tmp.BLAS_dgemv_d_s_x-f2c.c $@ && rm -f tmp.BLAS_dgemv_d_s_x-f2c.c*

$(SRC_PATH)/BLAS_dgemv_s_d_x-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="d, s, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dgemv_s_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv_s_d_x-f2c.c && \
    mv tmp.BLAS_dgemv_s_d_x-f2c.c $@ && rm -f tmp.BLAS_dgemv_s_d_x-f2c.c*

$(SRC_PATH)/BLAS_dgemv_s_s_x-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="d, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dgemv_s_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dgemv_s_s_x-f2c.c && \
    mv tmp.BLAS_dgemv_s_s_x-f2c.c $@ && rm -f tmp.BLAS_dgemv_s_s_x-f2c.c*

$(SRC_PATH)/BLAS_zgemv_z_c_x-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="z, z, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zgemv_z_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv_z_c_x-f2c.c && \
    mv tmp.BLAS_zgemv_z_c_x-f2c.c $@ && rm -f tmp.BLAS_zgemv_z_c_x-f2c.c*

$(SRC_PATH)/BLAS_zgemv_c_z_x-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="z, c, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zgemv_c_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv_c_z_x-f2c.c && \
    mv tmp.BLAS_zgemv_c_z_x-f2c.c $@ && rm -f tmp.BLAS_zgemv_c_z_x-f2c.c*

$(SRC_PATH)/BLAS_zgemv_c_c_x-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="z, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zgemv_c_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv_c_c_x-f2c.c && \
    mv tmp.BLAS_zgemv_c_c_x-f2c.c $@ && rm -f tmp.BLAS_zgemv_c_c_x-f2c.c*

$(SRC_PATH)/BLAS_cgemv_c_s_x-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="c, c, s, _x" ../f2c-bridge.m4 > tmp.BLAS_cgemv_c_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv_c_s_x-f2c.c && \
    mv tmp.BLAS_cgemv_c_s_x-f2c.c $@ && rm -f tmp.BLAS_cgemv_c_s_x-f2c.c*

$(SRC_PATH)/BLAS_cgemv_s_c_x-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="c, s, c, _x" ../f2c-bridge.m4 > tmp.BLAS_cgemv_s_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv_s_c_x-f2c.c && \
    mv tmp.BLAS_cgemv_s_c_x-f2c.c $@ && rm -f tmp.BLAS_cgemv_s_c_x-f2c.c*

$(SRC_PATH)/BLAS_cgemv_s_s_x-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="c, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_cgemv_s_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cgemv_s_s_x-f2c.c && \
    mv tmp.BLAS_cgemv_s_s_x-f2c.c $@ && rm -f tmp.BLAS_cgemv_s_s_x-f2c.c*

$(SRC_PATH)/BLAS_zgemv_z_d_x-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="z, z, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zgemv_z_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv_z_d_x-f2c.c && \
    mv tmp.BLAS_zgemv_z_d_x-f2c.c $@ && rm -f tmp.BLAS_zgemv_z_d_x-f2c.c*

$(SRC_PATH)/BLAS_zgemv_d_z_x-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="z, d, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zgemv_d_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv_d_z_x-f2c.c && \
    mv tmp.BLAS_zgemv_d_z_x-f2c.c $@ && rm -f tmp.BLAS_zgemv_d_z_x-f2c.c*

$(SRC_PATH)/BLAS_zgemv_d_d_x-f2c.c: ../f2c-bridge.m4 gemv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=gemv -DRARG="z, d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zgemv_d_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zgemv_d_d_x-f2c.c && \
    mv tmp.BLAS_zgemv_d_d_x-f2c.c $@ && rm -f tmp.BLAS_zgemv_d_d_x-f2c.c*


gemv.h: ../cblas.m4 ../header.template gemv-common.m4
	sed -e "s/routine/gemv/g" -e "s/ROUTINE/GEMV/g" \
  ../header.template | $(M4) $(M4_OPTS) $(M4_PATH) - >$@

source-clean:
	rm -f $(GEMV_SRCS) gemv.h

.PHONY: gemv-source gemv-header source-clean

