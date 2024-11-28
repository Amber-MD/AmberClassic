include ../../make.conf
include ../../$(MAKEINC)


SRC_PATH = ../../src/symv
M4_PATH = -I. -I..

SYMV_SRCS = \
	$(SRC_PATH)/BLAS_ssymv.c\
	$(SRC_PATH)/BLAS_dsymv.c\
	$(SRC_PATH)/BLAS_csymv.c\
	$(SRC_PATH)/BLAS_zsymv.c\
	$(SRC_PATH)/BLAS_dsymv_d_s.c\
	$(SRC_PATH)/BLAS_dsymv_s_d.c\
	$(SRC_PATH)/BLAS_dsymv_s_s.c\
	$(SRC_PATH)/BLAS_zsymv_z_c.c\
	$(SRC_PATH)/BLAS_zsymv_c_z.c\
	$(SRC_PATH)/BLAS_zsymv_c_c.c\
	$(SRC_PATH)/BLAS_csymv_c_s.c\
	$(SRC_PATH)/BLAS_csymv_s_c.c\
	$(SRC_PATH)/BLAS_csymv_s_s.c\
	$(SRC_PATH)/BLAS_zsymv_z_d.c\
	$(SRC_PATH)/BLAS_zsymv_d_z.c\
	$(SRC_PATH)/BLAS_zsymv_d_d.c\
	$(SRC_PATH)/BLAS_ssymv_x.c\
	$(SRC_PATH)/BLAS_dsymv_x.c\
	$(SRC_PATH)/BLAS_csymv_x.c\
	$(SRC_PATH)/BLAS_zsymv_x.c\
	$(SRC_PATH)/BLAS_dsymv_d_s_x.c\
	$(SRC_PATH)/BLAS_dsymv_s_d_x.c\
	$(SRC_PATH)/BLAS_dsymv_s_s_x.c\
	$(SRC_PATH)/BLAS_zsymv_z_c_x.c\
	$(SRC_PATH)/BLAS_zsymv_c_z_x.c\
	$(SRC_PATH)/BLAS_zsymv_c_c_x.c\
	$(SRC_PATH)/BLAS_csymv_c_s_x.c\
	$(SRC_PATH)/BLAS_csymv_s_c_x.c\
	$(SRC_PATH)/BLAS_csymv_s_s_x.c\
	$(SRC_PATH)/BLAS_zsymv_z_d_x.c\
	$(SRC_PATH)/BLAS_zsymv_d_z_x.c\
	$(SRC_PATH)/BLAS_zsymv_d_d_x.c \
	\
	$(SRC_PATH)/BLAS_ssymv-f2c.c\
	$(SRC_PATH)/BLAS_dsymv-f2c.c\
	$(SRC_PATH)/BLAS_csymv-f2c.c\
	$(SRC_PATH)/BLAS_zsymv-f2c.c\
	$(SRC_PATH)/BLAS_dsymv_d_s-f2c.c\
	$(SRC_PATH)/BLAS_dsymv_s_d-f2c.c\
	$(SRC_PATH)/BLAS_dsymv_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zsymv_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zsymv_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zsymv_c_c-f2c.c\
	$(SRC_PATH)/BLAS_csymv_c_s-f2c.c\
	$(SRC_PATH)/BLAS_csymv_s_c-f2c.c\
	$(SRC_PATH)/BLAS_csymv_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zsymv_z_d-f2c.c\
	$(SRC_PATH)/BLAS_zsymv_d_z-f2c.c\
	$(SRC_PATH)/BLAS_zsymv_d_d-f2c.c\
	$(SRC_PATH)/BLAS_ssymv_x-f2c.c\
	$(SRC_PATH)/BLAS_dsymv_x-f2c.c\
	$(SRC_PATH)/BLAS_csymv_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv_x-f2c.c\
	$(SRC_PATH)/BLAS_dsymv_d_s_x-f2c.c\
	$(SRC_PATH)/BLAS_dsymv_s_d_x-f2c.c\
	$(SRC_PATH)/BLAS_dsymv_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_csymv_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_csymv_s_c_x-f2c.c\
	$(SRC_PATH)/BLAS_csymv_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv_z_d_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv_d_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv_d_d_x-f2c.c

symv-source: $(SYMV_SRCS)

symv-header: symv.h

$(SRC_PATH)/BLAS_ssymv.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(s, s, s)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ssymv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ssymv.c && \
    mv tmp.BLAS_ssymv.c $@ && rm -f tmp.BLAS_ssymv.c*

$(SRC_PATH)/BLAS_dsymv.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(d, d, d)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsymv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv.c && \
    mv tmp.BLAS_dsymv.c $@ && rm -f tmp.BLAS_dsymv.c*

$(SRC_PATH)/BLAS_csymv.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(c, c, c)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csymv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv.c && \
    mv tmp.BLAS_csymv.c $@ && rm -f tmp.BLAS_csymv.c*

$(SRC_PATH)/BLAS_zsymv.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(z, z, z)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv.c && \
    mv tmp.BLAS_zsymv.c $@ && rm -f tmp.BLAS_zsymv.c*

$(SRC_PATH)/BLAS_dsymv_d_s.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(d, d, s)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsymv_d_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv_d_s.c && \
    mv tmp.BLAS_dsymv_d_s.c $@ && rm -f tmp.BLAS_dsymv_d_s.c*

$(SRC_PATH)/BLAS_dsymv_s_d.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(d, s, d)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsymv_s_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv_s_d.c && \
    mv tmp.BLAS_dsymv_s_d.c $@ && rm -f tmp.BLAS_dsymv_s_d.c*

$(SRC_PATH)/BLAS_dsymv_s_s.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(d, s, s)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsymv_s_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv_s_s.c && \
    mv tmp.BLAS_dsymv_s_s.c $@ && rm -f tmp.BLAS_dsymv_s_s.c*

$(SRC_PATH)/BLAS_zsymv_z_c.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(z, z, c)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv_z_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv_z_c.c && \
    mv tmp.BLAS_zsymv_z_c.c $@ && rm -f tmp.BLAS_zsymv_z_c.c*

$(SRC_PATH)/BLAS_zsymv_c_z.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(z, c, z)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv_c_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv_c_z.c && \
    mv tmp.BLAS_zsymv_c_z.c $@ && rm -f tmp.BLAS_zsymv_c_z.c*

$(SRC_PATH)/BLAS_zsymv_c_c.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(z, c, c)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv_c_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv_c_c.c && \
    mv tmp.BLAS_zsymv_c_c.c $@ && rm -f tmp.BLAS_zsymv_c_c.c*

$(SRC_PATH)/BLAS_csymv_c_s.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(c, c, s)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csymv_c_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv_c_s.c && \
    mv tmp.BLAS_csymv_c_s.c $@ && rm -f tmp.BLAS_csymv_c_s.c*

$(SRC_PATH)/BLAS_csymv_s_c.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(c, s, c)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csymv_s_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv_s_c.c && \
    mv tmp.BLAS_csymv_s_c.c $@ && rm -f tmp.BLAS_csymv_s_c.c*

$(SRC_PATH)/BLAS_csymv_s_s.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(c, s, s)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csymv_s_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv_s_s.c && \
    mv tmp.BLAS_csymv_s_s.c $@ && rm -f tmp.BLAS_csymv_s_s.c*

$(SRC_PATH)/BLAS_zsymv_z_d.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(z, z, d)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv_z_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv_z_d.c && \
    mv tmp.BLAS_zsymv_z_d.c $@ && rm -f tmp.BLAS_zsymv_z_d.c*

$(SRC_PATH)/BLAS_zsymv_d_z.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(z, d, z)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv_d_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv_d_z.c && \
    mv tmp.BLAS_zsymv_d_z.c $@ && rm -f tmp.BLAS_zsymv_d_z.c*

$(SRC_PATH)/BLAS_zsymv_d_d.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(z, d, d)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv_d_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv_d_d.c && \
    mv tmp.BLAS_zsymv_d_d.c $@ && rm -f tmp.BLAS_zsymv_d_d.c*

$(SRC_PATH)/BLAS_ssymv_x.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(s, s, s, _x)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ssymv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ssymv_x.c && \
    mv tmp.BLAS_ssymv_x.c $@ && rm -f tmp.BLAS_ssymv_x.c*

$(SRC_PATH)/BLAS_dsymv_x.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(d, d, d, _x)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsymv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv_x.c && \
    mv tmp.BLAS_dsymv_x.c $@ && rm -f tmp.BLAS_dsymv_x.c*

$(SRC_PATH)/BLAS_csymv_x.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(c, c, c, _x)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csymv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv_x.c && \
    mv tmp.BLAS_csymv_x.c $@ && rm -f tmp.BLAS_csymv_x.c*

$(SRC_PATH)/BLAS_zsymv_x.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(z, z, z, _x)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv_x.c && \
    mv tmp.BLAS_zsymv_x.c $@ && rm -f tmp.BLAS_zsymv_x.c*

$(SRC_PATH)/BLAS_dsymv_d_s_x.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(d, d, s, _x)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsymv_d_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv_d_s_x.c && \
    mv tmp.BLAS_dsymv_d_s_x.c $@ && rm -f tmp.BLAS_dsymv_d_s_x.c*

$(SRC_PATH)/BLAS_dsymv_s_d_x.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(d, s, d, _x)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsymv_s_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv_s_d_x.c && \
    mv tmp.BLAS_dsymv_s_d_x.c $@ && rm -f tmp.BLAS_dsymv_s_d_x.c*

$(SRC_PATH)/BLAS_dsymv_s_s_x.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(d, s, s, _x)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsymv_s_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv_s_s_x.c && \
    mv tmp.BLAS_dsymv_s_s_x.c $@ && rm -f tmp.BLAS_dsymv_s_s_x.c*

$(SRC_PATH)/BLAS_zsymv_z_c_x.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(z, z, c, _x)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv_z_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv_z_c_x.c && \
    mv tmp.BLAS_zsymv_z_c_x.c $@ && rm -f tmp.BLAS_zsymv_z_c_x.c*

$(SRC_PATH)/BLAS_zsymv_c_z_x.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(z, c, z, _x)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv_c_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv_c_z_x.c && \
    mv tmp.BLAS_zsymv_c_z_x.c $@ && rm -f tmp.BLAS_zsymv_c_z_x.c*

$(SRC_PATH)/BLAS_zsymv_c_c_x.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(z, c, c, _x)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv_c_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv_c_c_x.c && \
    mv tmp.BLAS_zsymv_c_c_x.c $@ && rm -f tmp.BLAS_zsymv_c_c_x.c*

$(SRC_PATH)/BLAS_csymv_c_s_x.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(c, c, s, _x)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csymv_c_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv_c_s_x.c && \
    mv tmp.BLAS_csymv_c_s_x.c $@ && rm -f tmp.BLAS_csymv_c_s_x.c*

$(SRC_PATH)/BLAS_csymv_s_c_x.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(c, s, c, _x)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csymv_s_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv_s_c_x.c && \
    mv tmp.BLAS_csymv_s_c_x.c $@ && rm -f tmp.BLAS_csymv_s_c_x.c*

$(SRC_PATH)/BLAS_csymv_s_s_x.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(c, s, s, _x)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csymv_s_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv_s_s_x.c && \
    mv tmp.BLAS_csymv_s_s_x.c $@ && rm -f tmp.BLAS_csymv_s_s_x.c*

$(SRC_PATH)/BLAS_zsymv_z_d_x.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(z, z, d, _x)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv_z_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv_z_d_x.c && \
    mv tmp.BLAS_zsymv_z_d_x.c $@ && rm -f tmp.BLAS_zsymv_z_d_x.c*

$(SRC_PATH)/BLAS_zsymv_d_z_x.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(z, d, z, _x)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv_d_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv_d_z_x.c && \
    mv tmp.BLAS_zsymv_d_z_x.c $@ && rm -f tmp.BLAS_zsymv_d_z_x.c*

$(SRC_PATH)/BLAS_zsymv_d_d_x.c: symv.m4 symv-common.m4 ../cblas.m4
	echo "SYMV(z, d, d, _x)" | cat symv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv_d_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv_d_d_x.c && \
    mv tmp.BLAS_zsymv_d_d_x.c $@ && rm -f tmp.BLAS_zsymv_d_d_x.c*



$(SRC_PATH)/BLAS_ssymv-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="s, s, s" ../f2c-bridge.m4 > tmp.BLAS_ssymv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ssymv-f2c.c && \
    mv tmp.BLAS_ssymv-f2c.c $@ && rm -f tmp.BLAS_ssymv-f2c.c*

$(SRC_PATH)/BLAS_dsymv-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="d, d, d" ../f2c-bridge.m4 > tmp.BLAS_dsymv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv-f2c.c && \
    mv tmp.BLAS_dsymv-f2c.c $@ && rm -f tmp.BLAS_dsymv-f2c.c*

$(SRC_PATH)/BLAS_csymv-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="c, c, c" ../f2c-bridge.m4 > tmp.BLAS_csymv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv-f2c.c && \
    mv tmp.BLAS_csymv-f2c.c $@ && rm -f tmp.BLAS_csymv-f2c.c*

$(SRC_PATH)/BLAS_zsymv-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="z, z, z" ../f2c-bridge.m4 > tmp.BLAS_zsymv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv-f2c.c && \
    mv tmp.BLAS_zsymv-f2c.c $@ && rm -f tmp.BLAS_zsymv-f2c.c*

$(SRC_PATH)/BLAS_dsymv_d_s-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="d, d, s" ../f2c-bridge.m4 > tmp.BLAS_dsymv_d_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv_d_s-f2c.c && \
    mv tmp.BLAS_dsymv_d_s-f2c.c $@ && rm -f tmp.BLAS_dsymv_d_s-f2c.c*

$(SRC_PATH)/BLAS_dsymv_s_d-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="d, s, d" ../f2c-bridge.m4 > tmp.BLAS_dsymv_s_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv_s_d-f2c.c && \
    mv tmp.BLAS_dsymv_s_d-f2c.c $@ && rm -f tmp.BLAS_dsymv_s_d-f2c.c*

$(SRC_PATH)/BLAS_dsymv_s_s-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="d, s, s" ../f2c-bridge.m4 > tmp.BLAS_dsymv_s_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv_s_s-f2c.c && \
    mv tmp.BLAS_dsymv_s_s-f2c.c $@ && rm -f tmp.BLAS_dsymv_s_s-f2c.c*

$(SRC_PATH)/BLAS_zsymv_z_c-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="z, z, c" ../f2c-bridge.m4 > tmp.BLAS_zsymv_z_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv_z_c-f2c.c && \
    mv tmp.BLAS_zsymv_z_c-f2c.c $@ && rm -f tmp.BLAS_zsymv_z_c-f2c.c*

$(SRC_PATH)/BLAS_zsymv_c_z-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="z, c, z" ../f2c-bridge.m4 > tmp.BLAS_zsymv_c_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv_c_z-f2c.c && \
    mv tmp.BLAS_zsymv_c_z-f2c.c $@ && rm -f tmp.BLAS_zsymv_c_z-f2c.c*

$(SRC_PATH)/BLAS_zsymv_c_c-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="z, c, c" ../f2c-bridge.m4 > tmp.BLAS_zsymv_c_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv_c_c-f2c.c && \
    mv tmp.BLAS_zsymv_c_c-f2c.c $@ && rm -f tmp.BLAS_zsymv_c_c-f2c.c*

$(SRC_PATH)/BLAS_csymv_c_s-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="c, c, s" ../f2c-bridge.m4 > tmp.BLAS_csymv_c_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv_c_s-f2c.c && \
    mv tmp.BLAS_csymv_c_s-f2c.c $@ && rm -f tmp.BLAS_csymv_c_s-f2c.c*

$(SRC_PATH)/BLAS_csymv_s_c-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="c, s, c" ../f2c-bridge.m4 > tmp.BLAS_csymv_s_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv_s_c-f2c.c && \
    mv tmp.BLAS_csymv_s_c-f2c.c $@ && rm -f tmp.BLAS_csymv_s_c-f2c.c*

$(SRC_PATH)/BLAS_csymv_s_s-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="c, s, s" ../f2c-bridge.m4 > tmp.BLAS_csymv_s_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv_s_s-f2c.c && \
    mv tmp.BLAS_csymv_s_s-f2c.c $@ && rm -f tmp.BLAS_csymv_s_s-f2c.c*

$(SRC_PATH)/BLAS_zsymv_z_d-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="z, z, d" ../f2c-bridge.m4 > tmp.BLAS_zsymv_z_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv_z_d-f2c.c && \
    mv tmp.BLAS_zsymv_z_d-f2c.c $@ && rm -f tmp.BLAS_zsymv_z_d-f2c.c*

$(SRC_PATH)/BLAS_zsymv_d_z-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="z, d, z" ../f2c-bridge.m4 > tmp.BLAS_zsymv_d_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv_d_z-f2c.c && \
    mv tmp.BLAS_zsymv_d_z-f2c.c $@ && rm -f tmp.BLAS_zsymv_d_z-f2c.c*

$(SRC_PATH)/BLAS_zsymv_d_d-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="z, d, d" ../f2c-bridge.m4 > tmp.BLAS_zsymv_d_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv_d_d-f2c.c && \
    mv tmp.BLAS_zsymv_d_d-f2c.c $@ && rm -f tmp.BLAS_zsymv_d_d-f2c.c*

$(SRC_PATH)/BLAS_ssymv_x-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="s, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_ssymv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ssymv_x-f2c.c && \
    mv tmp.BLAS_ssymv_x-f2c.c $@ && rm -f tmp.BLAS_ssymv_x-f2c.c*

$(SRC_PATH)/BLAS_dsymv_x-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="d, d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dsymv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv_x-f2c.c && \
    mv tmp.BLAS_dsymv_x-f2c.c $@ && rm -f tmp.BLAS_dsymv_x-f2c.c*

$(SRC_PATH)/BLAS_csymv_x-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="c, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_csymv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv_x-f2c.c && \
    mv tmp.BLAS_csymv_x-f2c.c $@ && rm -f tmp.BLAS_csymv_x-f2c.c*

$(SRC_PATH)/BLAS_zsymv_x-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="z, z, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zsymv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv_x-f2c.c && \
    mv tmp.BLAS_zsymv_x-f2c.c $@ && rm -f tmp.BLAS_zsymv_x-f2c.c*

$(SRC_PATH)/BLAS_dsymv_d_s_x-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="d, d, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dsymv_d_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv_d_s_x-f2c.c && \
    mv tmp.BLAS_dsymv_d_s_x-f2c.c $@ && rm -f tmp.BLAS_dsymv_d_s_x-f2c.c*

$(SRC_PATH)/BLAS_dsymv_s_d_x-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="d, s, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dsymv_s_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv_s_d_x-f2c.c && \
    mv tmp.BLAS_dsymv_s_d_x-f2c.c $@ && rm -f tmp.BLAS_dsymv_s_d_x-f2c.c*

$(SRC_PATH)/BLAS_dsymv_s_s_x-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="d, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dsymv_s_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv_s_s_x-f2c.c && \
    mv tmp.BLAS_dsymv_s_s_x-f2c.c $@ && rm -f tmp.BLAS_dsymv_s_s_x-f2c.c*

$(SRC_PATH)/BLAS_zsymv_z_c_x-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="z, z, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zsymv_z_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv_z_c_x-f2c.c && \
    mv tmp.BLAS_zsymv_z_c_x-f2c.c $@ && rm -f tmp.BLAS_zsymv_z_c_x-f2c.c*

$(SRC_PATH)/BLAS_zsymv_c_z_x-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="z, c, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zsymv_c_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv_c_z_x-f2c.c && \
    mv tmp.BLAS_zsymv_c_z_x-f2c.c $@ && rm -f tmp.BLAS_zsymv_c_z_x-f2c.c*

$(SRC_PATH)/BLAS_zsymv_c_c_x-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="z, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zsymv_c_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv_c_c_x-f2c.c && \
    mv tmp.BLAS_zsymv_c_c_x-f2c.c $@ && rm -f tmp.BLAS_zsymv_c_c_x-f2c.c*

$(SRC_PATH)/BLAS_csymv_c_s_x-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="c, c, s, _x" ../f2c-bridge.m4 > tmp.BLAS_csymv_c_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv_c_s_x-f2c.c && \
    mv tmp.BLAS_csymv_c_s_x-f2c.c $@ && rm -f tmp.BLAS_csymv_c_s_x-f2c.c*

$(SRC_PATH)/BLAS_csymv_s_c_x-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="c, s, c, _x" ../f2c-bridge.m4 > tmp.BLAS_csymv_s_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv_s_c_x-f2c.c && \
    mv tmp.BLAS_csymv_s_c_x-f2c.c $@ && rm -f tmp.BLAS_csymv_s_c_x-f2c.c*

$(SRC_PATH)/BLAS_csymv_s_s_x-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="c, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_csymv_s_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv_s_s_x-f2c.c && \
    mv tmp.BLAS_csymv_s_s_x-f2c.c $@ && rm -f tmp.BLAS_csymv_s_s_x-f2c.c*

$(SRC_PATH)/BLAS_zsymv_z_d_x-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="z, z, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zsymv_z_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv_z_d_x-f2c.c && \
    mv tmp.BLAS_zsymv_z_d_x-f2c.c $@ && rm -f tmp.BLAS_zsymv_z_d_x-f2c.c*

$(SRC_PATH)/BLAS_zsymv_d_z_x-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="z, d, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zsymv_d_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv_d_z_x-f2c.c && \
    mv tmp.BLAS_zsymv_d_z_x-f2c.c $@ && rm -f tmp.BLAS_zsymv_d_z_x-f2c.c*

$(SRC_PATH)/BLAS_zsymv_d_d_x-f2c.c: ../f2c-bridge.m4 symv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv -DRARG="z, d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zsymv_d_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv_d_d_x-f2c.c && \
    mv tmp.BLAS_zsymv_d_d_x-f2c.c $@ && rm -f tmp.BLAS_zsymv_d_d_x-f2c.c*


symv.h: ../cblas.m4 ../header.template symv-common.m4
	sed -e "s/routine/symv/g" -e "s/ROUTINE/SYMV/g" \
  ../header.template | $(M4) $(M4_OPTS) $(M4_PATH) - >$@

source-clean:
	rm -f $(SYMV_SRCS) symv.h

.PHONY: symv-source symv-header source-clean

