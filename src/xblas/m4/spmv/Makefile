include ../../make.conf
include ../../$(MAKEINC)


SRC_PATH = ../../src/spmv
M4_PATH = -I. -I..

SPMV_SRCS = \
	$(SRC_PATH)/BLAS_sspmv.c\
	$(SRC_PATH)/BLAS_dspmv.c\
	$(SRC_PATH)/BLAS_cspmv.c\
	$(SRC_PATH)/BLAS_zspmv.c\
	$(SRC_PATH)/BLAS_dspmv_d_s.c\
	$(SRC_PATH)/BLAS_dspmv_s_d.c\
	$(SRC_PATH)/BLAS_dspmv_s_s.c\
	$(SRC_PATH)/BLAS_zspmv_z_c.c\
	$(SRC_PATH)/BLAS_zspmv_c_z.c\
	$(SRC_PATH)/BLAS_zspmv_c_c.c\
	$(SRC_PATH)/BLAS_cspmv_c_s.c\
	$(SRC_PATH)/BLAS_cspmv_s_c.c\
	$(SRC_PATH)/BLAS_cspmv_s_s.c\
	$(SRC_PATH)/BLAS_zspmv_z_d.c\
	$(SRC_PATH)/BLAS_zspmv_d_z.c\
	$(SRC_PATH)/BLAS_zspmv_d_d.c\
	$(SRC_PATH)/BLAS_sspmv_x.c\
	$(SRC_PATH)/BLAS_dspmv_x.c\
	$(SRC_PATH)/BLAS_cspmv_x.c\
	$(SRC_PATH)/BLAS_zspmv_x.c\
	$(SRC_PATH)/BLAS_dspmv_d_s_x.c\
	$(SRC_PATH)/BLAS_dspmv_s_d_x.c\
	$(SRC_PATH)/BLAS_dspmv_s_s_x.c\
	$(SRC_PATH)/BLAS_zspmv_z_c_x.c\
	$(SRC_PATH)/BLAS_zspmv_c_z_x.c\
	$(SRC_PATH)/BLAS_zspmv_c_c_x.c\
	$(SRC_PATH)/BLAS_cspmv_c_s_x.c\
	$(SRC_PATH)/BLAS_cspmv_s_c_x.c\
	$(SRC_PATH)/BLAS_cspmv_s_s_x.c\
	$(SRC_PATH)/BLAS_zspmv_z_d_x.c\
	$(SRC_PATH)/BLAS_zspmv_d_z_x.c\
	$(SRC_PATH)/BLAS_zspmv_d_d_x.c \
	\
	$(SRC_PATH)/BLAS_sspmv-f2c.c\
	$(SRC_PATH)/BLAS_dspmv-f2c.c\
	$(SRC_PATH)/BLAS_cspmv-f2c.c\
	$(SRC_PATH)/BLAS_zspmv-f2c.c\
	$(SRC_PATH)/BLAS_dspmv_d_s-f2c.c\
	$(SRC_PATH)/BLAS_dspmv_s_d-f2c.c\
	$(SRC_PATH)/BLAS_dspmv_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zspmv_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zspmv_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zspmv_c_c-f2c.c\
	$(SRC_PATH)/BLAS_cspmv_c_s-f2c.c\
	$(SRC_PATH)/BLAS_cspmv_s_c-f2c.c\
	$(SRC_PATH)/BLAS_cspmv_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zspmv_z_d-f2c.c\
	$(SRC_PATH)/BLAS_zspmv_d_z-f2c.c\
	$(SRC_PATH)/BLAS_zspmv_d_d-f2c.c\
	$(SRC_PATH)/BLAS_sspmv_x-f2c.c\
	$(SRC_PATH)/BLAS_dspmv_x-f2c.c\
	$(SRC_PATH)/BLAS_cspmv_x-f2c.c\
	$(SRC_PATH)/BLAS_zspmv_x-f2c.c\
	$(SRC_PATH)/BLAS_dspmv_d_s_x-f2c.c\
	$(SRC_PATH)/BLAS_dspmv_s_d_x-f2c.c\
	$(SRC_PATH)/BLAS_dspmv_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zspmv_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zspmv_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zspmv_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cspmv_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_cspmv_s_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cspmv_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zspmv_z_d_x-f2c.c\
	$(SRC_PATH)/BLAS_zspmv_d_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zspmv_d_d_x-f2c.c

spmv-source: $(SPMV_SRCS)

spmv-header: spmv.h

$(SRC_PATH)/BLAS_sspmv.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(s, s, s)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_sspmv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sspmv.c && \
    mv tmp.BLAS_sspmv.c $@ && rm -f tmp.BLAS_sspmv.c*

$(SRC_PATH)/BLAS_dspmv.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(d, d, d)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dspmv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dspmv.c && \
    mv tmp.BLAS_dspmv.c $@ && rm -f tmp.BLAS_dspmv.c*

$(SRC_PATH)/BLAS_cspmv.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(c, c, c)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cspmv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cspmv.c && \
    mv tmp.BLAS_cspmv.c $@ && rm -f tmp.BLAS_cspmv.c*

$(SRC_PATH)/BLAS_zspmv.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(z, z, z)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zspmv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv.c && \
    mv tmp.BLAS_zspmv.c $@ && rm -f tmp.BLAS_zspmv.c*

$(SRC_PATH)/BLAS_dspmv_d_s.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(d, d, s)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dspmv_d_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dspmv_d_s.c && \
    mv tmp.BLAS_dspmv_d_s.c $@ && rm -f tmp.BLAS_dspmv_d_s.c*

$(SRC_PATH)/BLAS_dspmv_s_d.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(d, s, d)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dspmv_s_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dspmv_s_d.c && \
    mv tmp.BLAS_dspmv_s_d.c $@ && rm -f tmp.BLAS_dspmv_s_d.c*

$(SRC_PATH)/BLAS_dspmv_s_s.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(d, s, s)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dspmv_s_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dspmv_s_s.c && \
    mv tmp.BLAS_dspmv_s_s.c $@ && rm -f tmp.BLAS_dspmv_s_s.c*

$(SRC_PATH)/BLAS_zspmv_z_c.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(z, z, c)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zspmv_z_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv_z_c.c && \
    mv tmp.BLAS_zspmv_z_c.c $@ && rm -f tmp.BLAS_zspmv_z_c.c*

$(SRC_PATH)/BLAS_zspmv_c_z.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(z, c, z)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zspmv_c_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv_c_z.c && \
    mv tmp.BLAS_zspmv_c_z.c $@ && rm -f tmp.BLAS_zspmv_c_z.c*

$(SRC_PATH)/BLAS_zspmv_c_c.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(z, c, c)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zspmv_c_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv_c_c.c && \
    mv tmp.BLAS_zspmv_c_c.c $@ && rm -f tmp.BLAS_zspmv_c_c.c*

$(SRC_PATH)/BLAS_cspmv_c_s.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(c, c, s)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cspmv_c_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cspmv_c_s.c && \
    mv tmp.BLAS_cspmv_c_s.c $@ && rm -f tmp.BLAS_cspmv_c_s.c*

$(SRC_PATH)/BLAS_cspmv_s_c.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(c, s, c)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cspmv_s_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cspmv_s_c.c && \
    mv tmp.BLAS_cspmv_s_c.c $@ && rm -f tmp.BLAS_cspmv_s_c.c*

$(SRC_PATH)/BLAS_cspmv_s_s.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(c, s, s)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cspmv_s_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cspmv_s_s.c && \
    mv tmp.BLAS_cspmv_s_s.c $@ && rm -f tmp.BLAS_cspmv_s_s.c*

$(SRC_PATH)/BLAS_zspmv_z_d.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(z, z, d)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zspmv_z_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv_z_d.c && \
    mv tmp.BLAS_zspmv_z_d.c $@ && rm -f tmp.BLAS_zspmv_z_d.c*

$(SRC_PATH)/BLAS_zspmv_d_z.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(z, d, z)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zspmv_d_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv_d_z.c && \
    mv tmp.BLAS_zspmv_d_z.c $@ && rm -f tmp.BLAS_zspmv_d_z.c*

$(SRC_PATH)/BLAS_zspmv_d_d.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(z, d, d)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zspmv_d_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv_d_d.c && \
    mv tmp.BLAS_zspmv_d_d.c $@ && rm -f tmp.BLAS_zspmv_d_d.c*

$(SRC_PATH)/BLAS_sspmv_x.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(s, s, s, _x)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_sspmv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sspmv_x.c && \
    mv tmp.BLAS_sspmv_x.c $@ && rm -f tmp.BLAS_sspmv_x.c*

$(SRC_PATH)/BLAS_dspmv_x.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(d, d, d, _x)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dspmv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dspmv_x.c && \
    mv tmp.BLAS_dspmv_x.c $@ && rm -f tmp.BLAS_dspmv_x.c*

$(SRC_PATH)/BLAS_cspmv_x.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(c, c, c, _x)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cspmv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cspmv_x.c && \
    mv tmp.BLAS_cspmv_x.c $@ && rm -f tmp.BLAS_cspmv_x.c*

$(SRC_PATH)/BLAS_zspmv_x.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(z, z, z, _x)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zspmv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv_x.c && \
    mv tmp.BLAS_zspmv_x.c $@ && rm -f tmp.BLAS_zspmv_x.c*

$(SRC_PATH)/BLAS_dspmv_d_s_x.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(d, d, s, _x)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dspmv_d_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dspmv_d_s_x.c && \
    mv tmp.BLAS_dspmv_d_s_x.c $@ && rm -f tmp.BLAS_dspmv_d_s_x.c*

$(SRC_PATH)/BLAS_dspmv_s_d_x.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(d, s, d, _x)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dspmv_s_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dspmv_s_d_x.c && \
    mv tmp.BLAS_dspmv_s_d_x.c $@ && rm -f tmp.BLAS_dspmv_s_d_x.c*

$(SRC_PATH)/BLAS_dspmv_s_s_x.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(d, s, s, _x)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dspmv_s_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dspmv_s_s_x.c && \
    mv tmp.BLAS_dspmv_s_s_x.c $@ && rm -f tmp.BLAS_dspmv_s_s_x.c*

$(SRC_PATH)/BLAS_zspmv_z_c_x.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(z, z, c, _x)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zspmv_z_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv_z_c_x.c && \
    mv tmp.BLAS_zspmv_z_c_x.c $@ && rm -f tmp.BLAS_zspmv_z_c_x.c*

$(SRC_PATH)/BLAS_zspmv_c_z_x.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(z, c, z, _x)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zspmv_c_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv_c_z_x.c && \
    mv tmp.BLAS_zspmv_c_z_x.c $@ && rm -f tmp.BLAS_zspmv_c_z_x.c*

$(SRC_PATH)/BLAS_zspmv_c_c_x.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(z, c, c, _x)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zspmv_c_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv_c_c_x.c && \
    mv tmp.BLAS_zspmv_c_c_x.c $@ && rm -f tmp.BLAS_zspmv_c_c_x.c*

$(SRC_PATH)/BLAS_cspmv_c_s_x.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(c, c, s, _x)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cspmv_c_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cspmv_c_s_x.c && \
    mv tmp.BLAS_cspmv_c_s_x.c $@ && rm -f tmp.BLAS_cspmv_c_s_x.c*

$(SRC_PATH)/BLAS_cspmv_s_c_x.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(c, s, c, _x)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cspmv_s_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cspmv_s_c_x.c && \
    mv tmp.BLAS_cspmv_s_c_x.c $@ && rm -f tmp.BLAS_cspmv_s_c_x.c*

$(SRC_PATH)/BLAS_cspmv_s_s_x.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(c, s, s, _x)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cspmv_s_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cspmv_s_s_x.c && \
    mv tmp.BLAS_cspmv_s_s_x.c $@ && rm -f tmp.BLAS_cspmv_s_s_x.c*

$(SRC_PATH)/BLAS_zspmv_z_d_x.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(z, z, d, _x)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zspmv_z_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv_z_d_x.c && \
    mv tmp.BLAS_zspmv_z_d_x.c $@ && rm -f tmp.BLAS_zspmv_z_d_x.c*

$(SRC_PATH)/BLAS_zspmv_d_z_x.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(z, d, z, _x)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zspmv_d_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv_d_z_x.c && \
    mv tmp.BLAS_zspmv_d_z_x.c $@ && rm -f tmp.BLAS_zspmv_d_z_x.c*

$(SRC_PATH)/BLAS_zspmv_d_d_x.c: spmv.m4 spmv-common.m4 ../cblas.m4
	echo "SPMV(z, d, d, _x)" | cat spmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zspmv_d_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv_d_d_x.c && \
    mv tmp.BLAS_zspmv_d_d_x.c $@ && rm -f tmp.BLAS_zspmv_d_d_x.c*



$(SRC_PATH)/BLAS_sspmv-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="s, s, s" ../f2c-bridge.m4 > tmp.BLAS_sspmv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sspmv-f2c.c && \
    mv tmp.BLAS_sspmv-f2c.c $@ && rm -f tmp.BLAS_sspmv-f2c.c*

$(SRC_PATH)/BLAS_dspmv-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="d, d, d" ../f2c-bridge.m4 > tmp.BLAS_dspmv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dspmv-f2c.c && \
    mv tmp.BLAS_dspmv-f2c.c $@ && rm -f tmp.BLAS_dspmv-f2c.c*

$(SRC_PATH)/BLAS_cspmv-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="c, c, c" ../f2c-bridge.m4 > tmp.BLAS_cspmv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cspmv-f2c.c && \
    mv tmp.BLAS_cspmv-f2c.c $@ && rm -f tmp.BLAS_cspmv-f2c.c*

$(SRC_PATH)/BLAS_zspmv-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="z, z, z" ../f2c-bridge.m4 > tmp.BLAS_zspmv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv-f2c.c && \
    mv tmp.BLAS_zspmv-f2c.c $@ && rm -f tmp.BLAS_zspmv-f2c.c*

$(SRC_PATH)/BLAS_dspmv_d_s-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="d, d, s" ../f2c-bridge.m4 > tmp.BLAS_dspmv_d_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dspmv_d_s-f2c.c && \
    mv tmp.BLAS_dspmv_d_s-f2c.c $@ && rm -f tmp.BLAS_dspmv_d_s-f2c.c*

$(SRC_PATH)/BLAS_dspmv_s_d-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="d, s, d" ../f2c-bridge.m4 > tmp.BLAS_dspmv_s_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dspmv_s_d-f2c.c && \
    mv tmp.BLAS_dspmv_s_d-f2c.c $@ && rm -f tmp.BLAS_dspmv_s_d-f2c.c*

$(SRC_PATH)/BLAS_dspmv_s_s-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="d, s, s" ../f2c-bridge.m4 > tmp.BLAS_dspmv_s_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dspmv_s_s-f2c.c && \
    mv tmp.BLAS_dspmv_s_s-f2c.c $@ && rm -f tmp.BLAS_dspmv_s_s-f2c.c*

$(SRC_PATH)/BLAS_zspmv_z_c-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="z, z, c" ../f2c-bridge.m4 > tmp.BLAS_zspmv_z_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv_z_c-f2c.c && \
    mv tmp.BLAS_zspmv_z_c-f2c.c $@ && rm -f tmp.BLAS_zspmv_z_c-f2c.c*

$(SRC_PATH)/BLAS_zspmv_c_z-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="z, c, z" ../f2c-bridge.m4 > tmp.BLAS_zspmv_c_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv_c_z-f2c.c && \
    mv tmp.BLAS_zspmv_c_z-f2c.c $@ && rm -f tmp.BLAS_zspmv_c_z-f2c.c*

$(SRC_PATH)/BLAS_zspmv_c_c-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="z, c, c" ../f2c-bridge.m4 > tmp.BLAS_zspmv_c_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv_c_c-f2c.c && \
    mv tmp.BLAS_zspmv_c_c-f2c.c $@ && rm -f tmp.BLAS_zspmv_c_c-f2c.c*

$(SRC_PATH)/BLAS_cspmv_c_s-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="c, c, s" ../f2c-bridge.m4 > tmp.BLAS_cspmv_c_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cspmv_c_s-f2c.c && \
    mv tmp.BLAS_cspmv_c_s-f2c.c $@ && rm -f tmp.BLAS_cspmv_c_s-f2c.c*

$(SRC_PATH)/BLAS_cspmv_s_c-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="c, s, c" ../f2c-bridge.m4 > tmp.BLAS_cspmv_s_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cspmv_s_c-f2c.c && \
    mv tmp.BLAS_cspmv_s_c-f2c.c $@ && rm -f tmp.BLAS_cspmv_s_c-f2c.c*

$(SRC_PATH)/BLAS_cspmv_s_s-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="c, s, s" ../f2c-bridge.m4 > tmp.BLAS_cspmv_s_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cspmv_s_s-f2c.c && \
    mv tmp.BLAS_cspmv_s_s-f2c.c $@ && rm -f tmp.BLAS_cspmv_s_s-f2c.c*

$(SRC_PATH)/BLAS_zspmv_z_d-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="z, z, d" ../f2c-bridge.m4 > tmp.BLAS_zspmv_z_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv_z_d-f2c.c && \
    mv tmp.BLAS_zspmv_z_d-f2c.c $@ && rm -f tmp.BLAS_zspmv_z_d-f2c.c*

$(SRC_PATH)/BLAS_zspmv_d_z-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="z, d, z" ../f2c-bridge.m4 > tmp.BLAS_zspmv_d_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv_d_z-f2c.c && \
    mv tmp.BLAS_zspmv_d_z-f2c.c $@ && rm -f tmp.BLAS_zspmv_d_z-f2c.c*

$(SRC_PATH)/BLAS_zspmv_d_d-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="z, d, d" ../f2c-bridge.m4 > tmp.BLAS_zspmv_d_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv_d_d-f2c.c && \
    mv tmp.BLAS_zspmv_d_d-f2c.c $@ && rm -f tmp.BLAS_zspmv_d_d-f2c.c*

$(SRC_PATH)/BLAS_sspmv_x-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="s, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_sspmv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sspmv_x-f2c.c && \
    mv tmp.BLAS_sspmv_x-f2c.c $@ && rm -f tmp.BLAS_sspmv_x-f2c.c*

$(SRC_PATH)/BLAS_dspmv_x-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="d, d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dspmv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dspmv_x-f2c.c && \
    mv tmp.BLAS_dspmv_x-f2c.c $@ && rm -f tmp.BLAS_dspmv_x-f2c.c*

$(SRC_PATH)/BLAS_cspmv_x-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="c, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_cspmv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cspmv_x-f2c.c && \
    mv tmp.BLAS_cspmv_x-f2c.c $@ && rm -f tmp.BLAS_cspmv_x-f2c.c*

$(SRC_PATH)/BLAS_zspmv_x-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="z, z, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zspmv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv_x-f2c.c && \
    mv tmp.BLAS_zspmv_x-f2c.c $@ && rm -f tmp.BLAS_zspmv_x-f2c.c*

$(SRC_PATH)/BLAS_dspmv_d_s_x-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="d, d, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dspmv_d_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dspmv_d_s_x-f2c.c && \
    mv tmp.BLAS_dspmv_d_s_x-f2c.c $@ && rm -f tmp.BLAS_dspmv_d_s_x-f2c.c*

$(SRC_PATH)/BLAS_dspmv_s_d_x-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="d, s, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dspmv_s_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dspmv_s_d_x-f2c.c && \
    mv tmp.BLAS_dspmv_s_d_x-f2c.c $@ && rm -f tmp.BLAS_dspmv_s_d_x-f2c.c*

$(SRC_PATH)/BLAS_dspmv_s_s_x-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="d, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dspmv_s_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dspmv_s_s_x-f2c.c && \
    mv tmp.BLAS_dspmv_s_s_x-f2c.c $@ && rm -f tmp.BLAS_dspmv_s_s_x-f2c.c*

$(SRC_PATH)/BLAS_zspmv_z_c_x-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="z, z, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zspmv_z_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv_z_c_x-f2c.c && \
    mv tmp.BLAS_zspmv_z_c_x-f2c.c $@ && rm -f tmp.BLAS_zspmv_z_c_x-f2c.c*

$(SRC_PATH)/BLAS_zspmv_c_z_x-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="z, c, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zspmv_c_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv_c_z_x-f2c.c && \
    mv tmp.BLAS_zspmv_c_z_x-f2c.c $@ && rm -f tmp.BLAS_zspmv_c_z_x-f2c.c*

$(SRC_PATH)/BLAS_zspmv_c_c_x-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="z, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zspmv_c_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv_c_c_x-f2c.c && \
    mv tmp.BLAS_zspmv_c_c_x-f2c.c $@ && rm -f tmp.BLAS_zspmv_c_c_x-f2c.c*

$(SRC_PATH)/BLAS_cspmv_c_s_x-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="c, c, s, _x" ../f2c-bridge.m4 > tmp.BLAS_cspmv_c_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cspmv_c_s_x-f2c.c && \
    mv tmp.BLAS_cspmv_c_s_x-f2c.c $@ && rm -f tmp.BLAS_cspmv_c_s_x-f2c.c*

$(SRC_PATH)/BLAS_cspmv_s_c_x-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="c, s, c, _x" ../f2c-bridge.m4 > tmp.BLAS_cspmv_s_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cspmv_s_c_x-f2c.c && \
    mv tmp.BLAS_cspmv_s_c_x-f2c.c $@ && rm -f tmp.BLAS_cspmv_s_c_x-f2c.c*

$(SRC_PATH)/BLAS_cspmv_s_s_x-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="c, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_cspmv_s_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cspmv_s_s_x-f2c.c && \
    mv tmp.BLAS_cspmv_s_s_x-f2c.c $@ && rm -f tmp.BLAS_cspmv_s_s_x-f2c.c*

$(SRC_PATH)/BLAS_zspmv_z_d_x-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="z, z, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zspmv_z_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv_z_d_x-f2c.c && \
    mv tmp.BLAS_zspmv_z_d_x-f2c.c $@ && rm -f tmp.BLAS_zspmv_z_d_x-f2c.c*

$(SRC_PATH)/BLAS_zspmv_d_z_x-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="z, d, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zspmv_d_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv_d_z_x-f2c.c && \
    mv tmp.BLAS_zspmv_d_z_x-f2c.c $@ && rm -f tmp.BLAS_zspmv_d_z_x-f2c.c*

$(SRC_PATH)/BLAS_zspmv_d_d_x-f2c.c: ../f2c-bridge.m4 spmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=spmv -DRARG="z, d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zspmv_d_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zspmv_d_d_x-f2c.c && \
    mv tmp.BLAS_zspmv_d_d_x-f2c.c $@ && rm -f tmp.BLAS_zspmv_d_d_x-f2c.c*


spmv.h: ../cblas.m4 ../header.template spmv-common.m4
	sed -e "s/routine/spmv/g" -e "s/ROUTINE/SPMV/g" \
  ../header.template | $(M4) $(M4_OPTS) $(M4_PATH) - >$@

source-clean:
	rm -f $(SPMV_SRCS) spmv.h

.PHONY: spmv-source spmv-header source-clean

