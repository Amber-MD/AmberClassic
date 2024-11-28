include ../../make.conf
include ../../$(MAKEINC)


SRC_PATH = ../../src/symv2
M4_PATH = -I. -I..

SYMV2_SRCS = \
	$(SRC_PATH)/BLAS_ssymv2.c\
	$(SRC_PATH)/BLAS_dsymv2.c\
	$(SRC_PATH)/BLAS_csymv2.c\
	$(SRC_PATH)/BLAS_zsymv2.c\
	$(SRC_PATH)/BLAS_dsymv2_d_s.c\
	$(SRC_PATH)/BLAS_dsymv2_s_d.c\
	$(SRC_PATH)/BLAS_dsymv2_s_s.c\
	$(SRC_PATH)/BLAS_zsymv2_z_c.c\
	$(SRC_PATH)/BLAS_zsymv2_c_z.c\
	$(SRC_PATH)/BLAS_zsymv2_c_c.c\
	$(SRC_PATH)/BLAS_csymv2_c_s.c\
	$(SRC_PATH)/BLAS_csymv2_s_c.c\
	$(SRC_PATH)/BLAS_csymv2_s_s.c\
	$(SRC_PATH)/BLAS_zsymv2_z_d.c\
	$(SRC_PATH)/BLAS_zsymv2_d_z.c\
	$(SRC_PATH)/BLAS_zsymv2_d_d.c\
	$(SRC_PATH)/BLAS_ssymv2_x.c\
	$(SRC_PATH)/BLAS_dsymv2_x.c\
	$(SRC_PATH)/BLAS_csymv2_x.c\
	$(SRC_PATH)/BLAS_zsymv2_x.c\
	$(SRC_PATH)/BLAS_dsymv2_d_s_x.c\
	$(SRC_PATH)/BLAS_dsymv2_s_d_x.c\
	$(SRC_PATH)/BLAS_dsymv2_s_s_x.c\
	$(SRC_PATH)/BLAS_zsymv2_z_c_x.c\
	$(SRC_PATH)/BLAS_zsymv2_c_z_x.c\
	$(SRC_PATH)/BLAS_zsymv2_c_c_x.c\
	$(SRC_PATH)/BLAS_csymv2_c_s_x.c\
	$(SRC_PATH)/BLAS_csymv2_s_c_x.c\
	$(SRC_PATH)/BLAS_csymv2_s_s_x.c\
	$(SRC_PATH)/BLAS_zsymv2_z_d_x.c\
	$(SRC_PATH)/BLAS_zsymv2_d_z_x.c\
	$(SRC_PATH)/BLAS_zsymv2_d_d_x.c \
	\
	$(SRC_PATH)/BLAS_ssymv2-f2c.c\
	$(SRC_PATH)/BLAS_dsymv2-f2c.c\
	$(SRC_PATH)/BLAS_csymv2-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2-f2c.c\
	$(SRC_PATH)/BLAS_dsymv2_d_s-f2c.c\
	$(SRC_PATH)/BLAS_dsymv2_s_d-f2c.c\
	$(SRC_PATH)/BLAS_dsymv2_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2_c_c-f2c.c\
	$(SRC_PATH)/BLAS_csymv2_c_s-f2c.c\
	$(SRC_PATH)/BLAS_csymv2_s_c-f2c.c\
	$(SRC_PATH)/BLAS_csymv2_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2_z_d-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2_d_z-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2_d_d-f2c.c\
	$(SRC_PATH)/BLAS_ssymv2_x-f2c.c\
	$(SRC_PATH)/BLAS_dsymv2_x-f2c.c\
	$(SRC_PATH)/BLAS_csymv2_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2_x-f2c.c\
	$(SRC_PATH)/BLAS_dsymv2_d_s_x-f2c.c\
	$(SRC_PATH)/BLAS_dsymv2_s_d_x-f2c.c\
	$(SRC_PATH)/BLAS_dsymv2_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_csymv2_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_csymv2_s_c_x-f2c.c\
	$(SRC_PATH)/BLAS_csymv2_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2_z_d_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2_d_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2_d_d_x-f2c.c

symv2-source: $(SYMV2_SRCS)

symv2-header: symv2.h

$(SRC_PATH)/BLAS_ssymv2.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(s, s, s)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ssymv2.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ssymv2.c && \
    mv tmp.BLAS_ssymv2.c $@ && rm -f tmp.BLAS_ssymv2.c*

$(SRC_PATH)/BLAS_dsymv2.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(d, d, d)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsymv2.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv2.c && \
    mv tmp.BLAS_dsymv2.c $@ && rm -f tmp.BLAS_dsymv2.c*

$(SRC_PATH)/BLAS_csymv2.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(c, c, c)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csymv2.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv2.c && \
    mv tmp.BLAS_csymv2.c $@ && rm -f tmp.BLAS_csymv2.c*

$(SRC_PATH)/BLAS_zsymv2.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(z, z, z)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv2.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2.c && \
    mv tmp.BLAS_zsymv2.c $@ && rm -f tmp.BLAS_zsymv2.c*

$(SRC_PATH)/BLAS_dsymv2_d_s.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(d, d, s)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsymv2_d_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv2_d_s.c && \
    mv tmp.BLAS_dsymv2_d_s.c $@ && rm -f tmp.BLAS_dsymv2_d_s.c*

$(SRC_PATH)/BLAS_dsymv2_s_d.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(d, s, d)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsymv2_s_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv2_s_d.c && \
    mv tmp.BLAS_dsymv2_s_d.c $@ && rm -f tmp.BLAS_dsymv2_s_d.c*

$(SRC_PATH)/BLAS_dsymv2_s_s.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(d, s, s)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsymv2_s_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv2_s_s.c && \
    mv tmp.BLAS_dsymv2_s_s.c $@ && rm -f tmp.BLAS_dsymv2_s_s.c*

$(SRC_PATH)/BLAS_zsymv2_z_c.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(z, z, c)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv2_z_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2_z_c.c && \
    mv tmp.BLAS_zsymv2_z_c.c $@ && rm -f tmp.BLAS_zsymv2_z_c.c*

$(SRC_PATH)/BLAS_zsymv2_c_z.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(z, c, z)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv2_c_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2_c_z.c && \
    mv tmp.BLAS_zsymv2_c_z.c $@ && rm -f tmp.BLAS_zsymv2_c_z.c*

$(SRC_PATH)/BLAS_zsymv2_c_c.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(z, c, c)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv2_c_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2_c_c.c && \
    mv tmp.BLAS_zsymv2_c_c.c $@ && rm -f tmp.BLAS_zsymv2_c_c.c*

$(SRC_PATH)/BLAS_csymv2_c_s.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(c, c, s)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csymv2_c_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv2_c_s.c && \
    mv tmp.BLAS_csymv2_c_s.c $@ && rm -f tmp.BLAS_csymv2_c_s.c*

$(SRC_PATH)/BLAS_csymv2_s_c.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(c, s, c)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csymv2_s_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv2_s_c.c && \
    mv tmp.BLAS_csymv2_s_c.c $@ && rm -f tmp.BLAS_csymv2_s_c.c*

$(SRC_PATH)/BLAS_csymv2_s_s.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(c, s, s)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csymv2_s_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv2_s_s.c && \
    mv tmp.BLAS_csymv2_s_s.c $@ && rm -f tmp.BLAS_csymv2_s_s.c*

$(SRC_PATH)/BLAS_zsymv2_z_d.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(z, z, d)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv2_z_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2_z_d.c && \
    mv tmp.BLAS_zsymv2_z_d.c $@ && rm -f tmp.BLAS_zsymv2_z_d.c*

$(SRC_PATH)/BLAS_zsymv2_d_z.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(z, d, z)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv2_d_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2_d_z.c && \
    mv tmp.BLAS_zsymv2_d_z.c $@ && rm -f tmp.BLAS_zsymv2_d_z.c*

$(SRC_PATH)/BLAS_zsymv2_d_d.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(z, d, d)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv2_d_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2_d_d.c && \
    mv tmp.BLAS_zsymv2_d_d.c $@ && rm -f tmp.BLAS_zsymv2_d_d.c*

$(SRC_PATH)/BLAS_ssymv2_x.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(s, s, s, _x)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ssymv2_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ssymv2_x.c && \
    mv tmp.BLAS_ssymv2_x.c $@ && rm -f tmp.BLAS_ssymv2_x.c*

$(SRC_PATH)/BLAS_dsymv2_x.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(d, d, d, _x)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsymv2_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv2_x.c && \
    mv tmp.BLAS_dsymv2_x.c $@ && rm -f tmp.BLAS_dsymv2_x.c*

$(SRC_PATH)/BLAS_csymv2_x.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(c, c, c, _x)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csymv2_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv2_x.c && \
    mv tmp.BLAS_csymv2_x.c $@ && rm -f tmp.BLAS_csymv2_x.c*

$(SRC_PATH)/BLAS_zsymv2_x.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(z, z, z, _x)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv2_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2_x.c && \
    mv tmp.BLAS_zsymv2_x.c $@ && rm -f tmp.BLAS_zsymv2_x.c*

$(SRC_PATH)/BLAS_dsymv2_d_s_x.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(d, d, s, _x)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsymv2_d_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv2_d_s_x.c && \
    mv tmp.BLAS_dsymv2_d_s_x.c $@ && rm -f tmp.BLAS_dsymv2_d_s_x.c*

$(SRC_PATH)/BLAS_dsymv2_s_d_x.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(d, s, d, _x)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsymv2_s_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv2_s_d_x.c && \
    mv tmp.BLAS_dsymv2_s_d_x.c $@ && rm -f tmp.BLAS_dsymv2_s_d_x.c*

$(SRC_PATH)/BLAS_dsymv2_s_s_x.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(d, s, s, _x)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsymv2_s_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv2_s_s_x.c && \
    mv tmp.BLAS_dsymv2_s_s_x.c $@ && rm -f tmp.BLAS_dsymv2_s_s_x.c*

$(SRC_PATH)/BLAS_zsymv2_z_c_x.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(z, z, c, _x)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv2_z_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2_z_c_x.c && \
    mv tmp.BLAS_zsymv2_z_c_x.c $@ && rm -f tmp.BLAS_zsymv2_z_c_x.c*

$(SRC_PATH)/BLAS_zsymv2_c_z_x.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(z, c, z, _x)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv2_c_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2_c_z_x.c && \
    mv tmp.BLAS_zsymv2_c_z_x.c $@ && rm -f tmp.BLAS_zsymv2_c_z_x.c*

$(SRC_PATH)/BLAS_zsymv2_c_c_x.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(z, c, c, _x)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv2_c_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2_c_c_x.c && \
    mv tmp.BLAS_zsymv2_c_c_x.c $@ && rm -f tmp.BLAS_zsymv2_c_c_x.c*

$(SRC_PATH)/BLAS_csymv2_c_s_x.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(c, c, s, _x)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csymv2_c_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv2_c_s_x.c && \
    mv tmp.BLAS_csymv2_c_s_x.c $@ && rm -f tmp.BLAS_csymv2_c_s_x.c*

$(SRC_PATH)/BLAS_csymv2_s_c_x.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(c, s, c, _x)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csymv2_s_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv2_s_c_x.c && \
    mv tmp.BLAS_csymv2_s_c_x.c $@ && rm -f tmp.BLAS_csymv2_s_c_x.c*

$(SRC_PATH)/BLAS_csymv2_s_s_x.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(c, s, s, _x)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csymv2_s_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv2_s_s_x.c && \
    mv tmp.BLAS_csymv2_s_s_x.c $@ && rm -f tmp.BLAS_csymv2_s_s_x.c*

$(SRC_PATH)/BLAS_zsymv2_z_d_x.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(z, z, d, _x)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv2_z_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2_z_d_x.c && \
    mv tmp.BLAS_zsymv2_z_d_x.c $@ && rm -f tmp.BLAS_zsymv2_z_d_x.c*

$(SRC_PATH)/BLAS_zsymv2_d_z_x.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(z, d, z, _x)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv2_d_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2_d_z_x.c && \
    mv tmp.BLAS_zsymv2_d_z_x.c $@ && rm -f tmp.BLAS_zsymv2_d_z_x.c*

$(SRC_PATH)/BLAS_zsymv2_d_d_x.c: symv2.m4 symv2-common.m4 ../cblas.m4
	echo "SYMV2(z, d, d, _x)" | cat symv2.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymv2_d_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2_d_d_x.c && \
    mv tmp.BLAS_zsymv2_d_d_x.c $@ && rm -f tmp.BLAS_zsymv2_d_d_x.c*



$(SRC_PATH)/BLAS_ssymv2-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="s, s, s" ../f2c-bridge.m4 > tmp.BLAS_ssymv2-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ssymv2-f2c.c && \
    mv tmp.BLAS_ssymv2-f2c.c $@ && rm -f tmp.BLAS_ssymv2-f2c.c*

$(SRC_PATH)/BLAS_dsymv2-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="d, d, d" ../f2c-bridge.m4 > tmp.BLAS_dsymv2-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv2-f2c.c && \
    mv tmp.BLAS_dsymv2-f2c.c $@ && rm -f tmp.BLAS_dsymv2-f2c.c*

$(SRC_PATH)/BLAS_csymv2-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="c, c, c" ../f2c-bridge.m4 > tmp.BLAS_csymv2-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv2-f2c.c && \
    mv tmp.BLAS_csymv2-f2c.c $@ && rm -f tmp.BLAS_csymv2-f2c.c*

$(SRC_PATH)/BLAS_zsymv2-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="z, z, z" ../f2c-bridge.m4 > tmp.BLAS_zsymv2-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2-f2c.c && \
    mv tmp.BLAS_zsymv2-f2c.c $@ && rm -f tmp.BLAS_zsymv2-f2c.c*

$(SRC_PATH)/BLAS_dsymv2_d_s-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="d, d, s" ../f2c-bridge.m4 > tmp.BLAS_dsymv2_d_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv2_d_s-f2c.c && \
    mv tmp.BLAS_dsymv2_d_s-f2c.c $@ && rm -f tmp.BLAS_dsymv2_d_s-f2c.c*

$(SRC_PATH)/BLAS_dsymv2_s_d-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="d, s, d" ../f2c-bridge.m4 > tmp.BLAS_dsymv2_s_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv2_s_d-f2c.c && \
    mv tmp.BLAS_dsymv2_s_d-f2c.c $@ && rm -f tmp.BLAS_dsymv2_s_d-f2c.c*

$(SRC_PATH)/BLAS_dsymv2_s_s-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="d, s, s" ../f2c-bridge.m4 > tmp.BLAS_dsymv2_s_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv2_s_s-f2c.c && \
    mv tmp.BLAS_dsymv2_s_s-f2c.c $@ && rm -f tmp.BLAS_dsymv2_s_s-f2c.c*

$(SRC_PATH)/BLAS_zsymv2_z_c-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="z, z, c" ../f2c-bridge.m4 > tmp.BLAS_zsymv2_z_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2_z_c-f2c.c && \
    mv tmp.BLAS_zsymv2_z_c-f2c.c $@ && rm -f tmp.BLAS_zsymv2_z_c-f2c.c*

$(SRC_PATH)/BLAS_zsymv2_c_z-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="z, c, z" ../f2c-bridge.m4 > tmp.BLAS_zsymv2_c_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2_c_z-f2c.c && \
    mv tmp.BLAS_zsymv2_c_z-f2c.c $@ && rm -f tmp.BLAS_zsymv2_c_z-f2c.c*

$(SRC_PATH)/BLAS_zsymv2_c_c-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="z, c, c" ../f2c-bridge.m4 > tmp.BLAS_zsymv2_c_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2_c_c-f2c.c && \
    mv tmp.BLAS_zsymv2_c_c-f2c.c $@ && rm -f tmp.BLAS_zsymv2_c_c-f2c.c*

$(SRC_PATH)/BLAS_csymv2_c_s-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="c, c, s" ../f2c-bridge.m4 > tmp.BLAS_csymv2_c_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv2_c_s-f2c.c && \
    mv tmp.BLAS_csymv2_c_s-f2c.c $@ && rm -f tmp.BLAS_csymv2_c_s-f2c.c*

$(SRC_PATH)/BLAS_csymv2_s_c-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="c, s, c" ../f2c-bridge.m4 > tmp.BLAS_csymv2_s_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv2_s_c-f2c.c && \
    mv tmp.BLAS_csymv2_s_c-f2c.c $@ && rm -f tmp.BLAS_csymv2_s_c-f2c.c*

$(SRC_PATH)/BLAS_csymv2_s_s-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="c, s, s" ../f2c-bridge.m4 > tmp.BLAS_csymv2_s_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv2_s_s-f2c.c && \
    mv tmp.BLAS_csymv2_s_s-f2c.c $@ && rm -f tmp.BLAS_csymv2_s_s-f2c.c*

$(SRC_PATH)/BLAS_zsymv2_z_d-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="z, z, d" ../f2c-bridge.m4 > tmp.BLAS_zsymv2_z_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2_z_d-f2c.c && \
    mv tmp.BLAS_zsymv2_z_d-f2c.c $@ && rm -f tmp.BLAS_zsymv2_z_d-f2c.c*

$(SRC_PATH)/BLAS_zsymv2_d_z-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="z, d, z" ../f2c-bridge.m4 > tmp.BLAS_zsymv2_d_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2_d_z-f2c.c && \
    mv tmp.BLAS_zsymv2_d_z-f2c.c $@ && rm -f tmp.BLAS_zsymv2_d_z-f2c.c*

$(SRC_PATH)/BLAS_zsymv2_d_d-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="z, d, d" ../f2c-bridge.m4 > tmp.BLAS_zsymv2_d_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2_d_d-f2c.c && \
    mv tmp.BLAS_zsymv2_d_d-f2c.c $@ && rm -f tmp.BLAS_zsymv2_d_d-f2c.c*

$(SRC_PATH)/BLAS_ssymv2_x-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="s, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_ssymv2_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ssymv2_x-f2c.c && \
    mv tmp.BLAS_ssymv2_x-f2c.c $@ && rm -f tmp.BLAS_ssymv2_x-f2c.c*

$(SRC_PATH)/BLAS_dsymv2_x-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="d, d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dsymv2_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv2_x-f2c.c && \
    mv tmp.BLAS_dsymv2_x-f2c.c $@ && rm -f tmp.BLAS_dsymv2_x-f2c.c*

$(SRC_PATH)/BLAS_csymv2_x-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="c, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_csymv2_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv2_x-f2c.c && \
    mv tmp.BLAS_csymv2_x-f2c.c $@ && rm -f tmp.BLAS_csymv2_x-f2c.c*

$(SRC_PATH)/BLAS_zsymv2_x-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="z, z, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zsymv2_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2_x-f2c.c && \
    mv tmp.BLAS_zsymv2_x-f2c.c $@ && rm -f tmp.BLAS_zsymv2_x-f2c.c*

$(SRC_PATH)/BLAS_dsymv2_d_s_x-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="d, d, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dsymv2_d_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv2_d_s_x-f2c.c && \
    mv tmp.BLAS_dsymv2_d_s_x-f2c.c $@ && rm -f tmp.BLAS_dsymv2_d_s_x-f2c.c*

$(SRC_PATH)/BLAS_dsymv2_s_d_x-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="d, s, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dsymv2_s_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv2_s_d_x-f2c.c && \
    mv tmp.BLAS_dsymv2_s_d_x-f2c.c $@ && rm -f tmp.BLAS_dsymv2_s_d_x-f2c.c*

$(SRC_PATH)/BLAS_dsymv2_s_s_x-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="d, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dsymv2_s_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymv2_s_s_x-f2c.c && \
    mv tmp.BLAS_dsymv2_s_s_x-f2c.c $@ && rm -f tmp.BLAS_dsymv2_s_s_x-f2c.c*

$(SRC_PATH)/BLAS_zsymv2_z_c_x-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="z, z, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zsymv2_z_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2_z_c_x-f2c.c && \
    mv tmp.BLAS_zsymv2_z_c_x-f2c.c $@ && rm -f tmp.BLAS_zsymv2_z_c_x-f2c.c*

$(SRC_PATH)/BLAS_zsymv2_c_z_x-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="z, c, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zsymv2_c_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2_c_z_x-f2c.c && \
    mv tmp.BLAS_zsymv2_c_z_x-f2c.c $@ && rm -f tmp.BLAS_zsymv2_c_z_x-f2c.c*

$(SRC_PATH)/BLAS_zsymv2_c_c_x-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="z, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zsymv2_c_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2_c_c_x-f2c.c && \
    mv tmp.BLAS_zsymv2_c_c_x-f2c.c $@ && rm -f tmp.BLAS_zsymv2_c_c_x-f2c.c*

$(SRC_PATH)/BLAS_csymv2_c_s_x-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="c, c, s, _x" ../f2c-bridge.m4 > tmp.BLAS_csymv2_c_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv2_c_s_x-f2c.c && \
    mv tmp.BLAS_csymv2_c_s_x-f2c.c $@ && rm -f tmp.BLAS_csymv2_c_s_x-f2c.c*

$(SRC_PATH)/BLAS_csymv2_s_c_x-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="c, s, c, _x" ../f2c-bridge.m4 > tmp.BLAS_csymv2_s_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv2_s_c_x-f2c.c && \
    mv tmp.BLAS_csymv2_s_c_x-f2c.c $@ && rm -f tmp.BLAS_csymv2_s_c_x-f2c.c*

$(SRC_PATH)/BLAS_csymv2_s_s_x-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="c, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_csymv2_s_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymv2_s_s_x-f2c.c && \
    mv tmp.BLAS_csymv2_s_s_x-f2c.c $@ && rm -f tmp.BLAS_csymv2_s_s_x-f2c.c*

$(SRC_PATH)/BLAS_zsymv2_z_d_x-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="z, z, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zsymv2_z_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2_z_d_x-f2c.c && \
    mv tmp.BLAS_zsymv2_z_d_x-f2c.c $@ && rm -f tmp.BLAS_zsymv2_z_d_x-f2c.c*

$(SRC_PATH)/BLAS_zsymv2_d_z_x-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="z, d, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zsymv2_d_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2_d_z_x-f2c.c && \
    mv tmp.BLAS_zsymv2_d_z_x-f2c.c $@ && rm -f tmp.BLAS_zsymv2_d_z_x-f2c.c*

$(SRC_PATH)/BLAS_zsymv2_d_d_x-f2c.c: ../f2c-bridge.m4 symv2-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symv2 -DRARG="z, d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zsymv2_d_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymv2_d_d_x-f2c.c && \
    mv tmp.BLAS_zsymv2_d_d_x-f2c.c $@ && rm -f tmp.BLAS_zsymv2_d_d_x-f2c.c*


symv2.h: ../cblas.m4 ../header.template symv2-common.m4
	sed -e "s/routine/symv2/g" -e "s/ROUTINE/SYMV2/g" \
  ../header.template | $(M4) $(M4_OPTS) $(M4_PATH) - >$@

source-clean:
	rm -f $(SYMV2_SRCS) symv2.h

.PHONY: symv2-source symv2-header source-clean

