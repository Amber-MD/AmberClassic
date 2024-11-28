include ../../make.conf
include ../../$(MAKEINC)


SRC_PATH = ../../src/symm
M4_PATH = -I. -I..

SYMM_SRCS = \
	$(SRC_PATH)/BLAS_ssymm.c\
	$(SRC_PATH)/BLAS_dsymm.c\
	$(SRC_PATH)/BLAS_csymm.c\
	$(SRC_PATH)/BLAS_zsymm.c\
	$(SRC_PATH)/BLAS_dsymm_d_s.c\
	$(SRC_PATH)/BLAS_dsymm_s_d.c\
	$(SRC_PATH)/BLAS_dsymm_s_s.c\
	$(SRC_PATH)/BLAS_zsymm_z_c.c\
	$(SRC_PATH)/BLAS_zsymm_c_z.c\
	$(SRC_PATH)/BLAS_zsymm_c_c.c\
	$(SRC_PATH)/BLAS_csymm_c_s.c\
	$(SRC_PATH)/BLAS_csymm_s_c.c\
	$(SRC_PATH)/BLAS_csymm_s_s.c\
	$(SRC_PATH)/BLAS_zsymm_z_d.c\
	$(SRC_PATH)/BLAS_zsymm_d_z.c\
	$(SRC_PATH)/BLAS_zsymm_d_d.c\
	$(SRC_PATH)/BLAS_ssymm_x.c\
	$(SRC_PATH)/BLAS_dsymm_x.c\
	$(SRC_PATH)/BLAS_csymm_x.c\
	$(SRC_PATH)/BLAS_zsymm_x.c\
	$(SRC_PATH)/BLAS_dsymm_d_s_x.c\
	$(SRC_PATH)/BLAS_dsymm_s_d_x.c\
	$(SRC_PATH)/BLAS_dsymm_s_s_x.c\
	$(SRC_PATH)/BLAS_zsymm_z_c_x.c\
	$(SRC_PATH)/BLAS_zsymm_c_z_x.c\
	$(SRC_PATH)/BLAS_zsymm_c_c_x.c\
	$(SRC_PATH)/BLAS_csymm_c_s_x.c\
	$(SRC_PATH)/BLAS_csymm_s_c_x.c\
	$(SRC_PATH)/BLAS_csymm_s_s_x.c\
	$(SRC_PATH)/BLAS_zsymm_z_d_x.c\
	$(SRC_PATH)/BLAS_zsymm_d_z_x.c\
	$(SRC_PATH)/BLAS_zsymm_d_d_x.c \
	\
	$(SRC_PATH)/BLAS_ssymm-f2c.c\
	$(SRC_PATH)/BLAS_dsymm-f2c.c\
	$(SRC_PATH)/BLAS_csymm-f2c.c\
	$(SRC_PATH)/BLAS_zsymm-f2c.c\
	$(SRC_PATH)/BLAS_dsymm_d_s-f2c.c\
	$(SRC_PATH)/BLAS_dsymm_s_d-f2c.c\
	$(SRC_PATH)/BLAS_dsymm_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zsymm_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zsymm_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zsymm_c_c-f2c.c\
	$(SRC_PATH)/BLAS_csymm_c_s-f2c.c\
	$(SRC_PATH)/BLAS_csymm_s_c-f2c.c\
	$(SRC_PATH)/BLAS_csymm_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zsymm_z_d-f2c.c\
	$(SRC_PATH)/BLAS_zsymm_d_z-f2c.c\
	$(SRC_PATH)/BLAS_zsymm_d_d-f2c.c\
	$(SRC_PATH)/BLAS_ssymm_x-f2c.c\
	$(SRC_PATH)/BLAS_dsymm_x-f2c.c\
	$(SRC_PATH)/BLAS_csymm_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymm_x-f2c.c\
	$(SRC_PATH)/BLAS_dsymm_d_s_x-f2c.c\
	$(SRC_PATH)/BLAS_dsymm_s_d_x-f2c.c\
	$(SRC_PATH)/BLAS_dsymm_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymm_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymm_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymm_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_csymm_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_csymm_s_c_x-f2c.c\
	$(SRC_PATH)/BLAS_csymm_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymm_z_d_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymm_d_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymm_d_d_x-f2c.c

symm-source: $(SYMM_SRCS)

symm-header: symm.h

$(SRC_PATH)/BLAS_ssymm.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(s, s, s)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ssymm.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ssymm.c && \
    mv tmp.BLAS_ssymm.c $@ && rm -f tmp.BLAS_ssymm.c*

$(SRC_PATH)/BLAS_dsymm.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(d, d, d)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsymm.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymm.c && \
    mv tmp.BLAS_dsymm.c $@ && rm -f tmp.BLAS_dsymm.c*

$(SRC_PATH)/BLAS_csymm.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(c, c, c)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csymm.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymm.c && \
    mv tmp.BLAS_csymm.c $@ && rm -f tmp.BLAS_csymm.c*

$(SRC_PATH)/BLAS_zsymm.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(z, z, z)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymm.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm.c && \
    mv tmp.BLAS_zsymm.c $@ && rm -f tmp.BLAS_zsymm.c*

$(SRC_PATH)/BLAS_dsymm_d_s.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(d, d, s)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsymm_d_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymm_d_s.c && \
    mv tmp.BLAS_dsymm_d_s.c $@ && rm -f tmp.BLAS_dsymm_d_s.c*

$(SRC_PATH)/BLAS_dsymm_s_d.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(d, s, d)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsymm_s_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymm_s_d.c && \
    mv tmp.BLAS_dsymm_s_d.c $@ && rm -f tmp.BLAS_dsymm_s_d.c*

$(SRC_PATH)/BLAS_dsymm_s_s.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(d, s, s)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsymm_s_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymm_s_s.c && \
    mv tmp.BLAS_dsymm_s_s.c $@ && rm -f tmp.BLAS_dsymm_s_s.c*

$(SRC_PATH)/BLAS_zsymm_z_c.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(z, z, c)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymm_z_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm_z_c.c && \
    mv tmp.BLAS_zsymm_z_c.c $@ && rm -f tmp.BLAS_zsymm_z_c.c*

$(SRC_PATH)/BLAS_zsymm_c_z.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(z, c, z)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymm_c_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm_c_z.c && \
    mv tmp.BLAS_zsymm_c_z.c $@ && rm -f tmp.BLAS_zsymm_c_z.c*

$(SRC_PATH)/BLAS_zsymm_c_c.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(z, c, c)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymm_c_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm_c_c.c && \
    mv tmp.BLAS_zsymm_c_c.c $@ && rm -f tmp.BLAS_zsymm_c_c.c*

$(SRC_PATH)/BLAS_csymm_c_s.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(c, c, s)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csymm_c_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymm_c_s.c && \
    mv tmp.BLAS_csymm_c_s.c $@ && rm -f tmp.BLAS_csymm_c_s.c*

$(SRC_PATH)/BLAS_csymm_s_c.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(c, s, c)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csymm_s_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymm_s_c.c && \
    mv tmp.BLAS_csymm_s_c.c $@ && rm -f tmp.BLAS_csymm_s_c.c*

$(SRC_PATH)/BLAS_csymm_s_s.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(c, s, s)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csymm_s_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymm_s_s.c && \
    mv tmp.BLAS_csymm_s_s.c $@ && rm -f tmp.BLAS_csymm_s_s.c*

$(SRC_PATH)/BLAS_zsymm_z_d.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(z, z, d)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymm_z_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm_z_d.c && \
    mv tmp.BLAS_zsymm_z_d.c $@ && rm -f tmp.BLAS_zsymm_z_d.c*

$(SRC_PATH)/BLAS_zsymm_d_z.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(z, d, z)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymm_d_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm_d_z.c && \
    mv tmp.BLAS_zsymm_d_z.c $@ && rm -f tmp.BLAS_zsymm_d_z.c*

$(SRC_PATH)/BLAS_zsymm_d_d.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(z, d, d)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymm_d_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm_d_d.c && \
    mv tmp.BLAS_zsymm_d_d.c $@ && rm -f tmp.BLAS_zsymm_d_d.c*

$(SRC_PATH)/BLAS_ssymm_x.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(s, s, s, _x)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ssymm_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ssymm_x.c && \
    mv tmp.BLAS_ssymm_x.c $@ && rm -f tmp.BLAS_ssymm_x.c*

$(SRC_PATH)/BLAS_dsymm_x.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(d, d, d, _x)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsymm_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymm_x.c && \
    mv tmp.BLAS_dsymm_x.c $@ && rm -f tmp.BLAS_dsymm_x.c*

$(SRC_PATH)/BLAS_csymm_x.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(c, c, c, _x)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csymm_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymm_x.c && \
    mv tmp.BLAS_csymm_x.c $@ && rm -f tmp.BLAS_csymm_x.c*

$(SRC_PATH)/BLAS_zsymm_x.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(z, z, z, _x)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymm_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm_x.c && \
    mv tmp.BLAS_zsymm_x.c $@ && rm -f tmp.BLAS_zsymm_x.c*

$(SRC_PATH)/BLAS_dsymm_d_s_x.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(d, d, s, _x)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsymm_d_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymm_d_s_x.c && \
    mv tmp.BLAS_dsymm_d_s_x.c $@ && rm -f tmp.BLAS_dsymm_d_s_x.c*

$(SRC_PATH)/BLAS_dsymm_s_d_x.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(d, s, d, _x)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsymm_s_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymm_s_d_x.c && \
    mv tmp.BLAS_dsymm_s_d_x.c $@ && rm -f tmp.BLAS_dsymm_s_d_x.c*

$(SRC_PATH)/BLAS_dsymm_s_s_x.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(d, s, s, _x)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsymm_s_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymm_s_s_x.c && \
    mv tmp.BLAS_dsymm_s_s_x.c $@ && rm -f tmp.BLAS_dsymm_s_s_x.c*

$(SRC_PATH)/BLAS_zsymm_z_c_x.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(z, z, c, _x)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymm_z_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm_z_c_x.c && \
    mv tmp.BLAS_zsymm_z_c_x.c $@ && rm -f tmp.BLAS_zsymm_z_c_x.c*

$(SRC_PATH)/BLAS_zsymm_c_z_x.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(z, c, z, _x)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymm_c_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm_c_z_x.c && \
    mv tmp.BLAS_zsymm_c_z_x.c $@ && rm -f tmp.BLAS_zsymm_c_z_x.c*

$(SRC_PATH)/BLAS_zsymm_c_c_x.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(z, c, c, _x)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymm_c_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm_c_c_x.c && \
    mv tmp.BLAS_zsymm_c_c_x.c $@ && rm -f tmp.BLAS_zsymm_c_c_x.c*

$(SRC_PATH)/BLAS_csymm_c_s_x.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(c, c, s, _x)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csymm_c_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymm_c_s_x.c && \
    mv tmp.BLAS_csymm_c_s_x.c $@ && rm -f tmp.BLAS_csymm_c_s_x.c*

$(SRC_PATH)/BLAS_csymm_s_c_x.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(c, s, c, _x)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csymm_s_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymm_s_c_x.c && \
    mv tmp.BLAS_csymm_s_c_x.c $@ && rm -f tmp.BLAS_csymm_s_c_x.c*

$(SRC_PATH)/BLAS_csymm_s_s_x.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(c, s, s, _x)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csymm_s_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymm_s_s_x.c && \
    mv tmp.BLAS_csymm_s_s_x.c $@ && rm -f tmp.BLAS_csymm_s_s_x.c*

$(SRC_PATH)/BLAS_zsymm_z_d_x.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(z, z, d, _x)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymm_z_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm_z_d_x.c && \
    mv tmp.BLAS_zsymm_z_d_x.c $@ && rm -f tmp.BLAS_zsymm_z_d_x.c*

$(SRC_PATH)/BLAS_zsymm_d_z_x.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(z, d, z, _x)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymm_d_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm_d_z_x.c && \
    mv tmp.BLAS_zsymm_d_z_x.c $@ && rm -f tmp.BLAS_zsymm_d_z_x.c*

$(SRC_PATH)/BLAS_zsymm_d_d_x.c: symm.m4 symm-common.m4 ../cblas.m4
	echo "SYMM(z, d, d, _x)" | cat symm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsymm_d_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm_d_d_x.c && \
    mv tmp.BLAS_zsymm_d_d_x.c $@ && rm -f tmp.BLAS_zsymm_d_d_x.c*



$(SRC_PATH)/BLAS_ssymm-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="s, s, s" ../f2c-bridge.m4 > tmp.BLAS_ssymm-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ssymm-f2c.c && \
    mv tmp.BLAS_ssymm-f2c.c $@ && rm -f tmp.BLAS_ssymm-f2c.c*

$(SRC_PATH)/BLAS_dsymm-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="d, d, d" ../f2c-bridge.m4 > tmp.BLAS_dsymm-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymm-f2c.c && \
    mv tmp.BLAS_dsymm-f2c.c $@ && rm -f tmp.BLAS_dsymm-f2c.c*

$(SRC_PATH)/BLAS_csymm-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="c, c, c" ../f2c-bridge.m4 > tmp.BLAS_csymm-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymm-f2c.c && \
    mv tmp.BLAS_csymm-f2c.c $@ && rm -f tmp.BLAS_csymm-f2c.c*

$(SRC_PATH)/BLAS_zsymm-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="z, z, z" ../f2c-bridge.m4 > tmp.BLAS_zsymm-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm-f2c.c && \
    mv tmp.BLAS_zsymm-f2c.c $@ && rm -f tmp.BLAS_zsymm-f2c.c*

$(SRC_PATH)/BLAS_dsymm_d_s-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="d, d, s" ../f2c-bridge.m4 > tmp.BLAS_dsymm_d_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymm_d_s-f2c.c && \
    mv tmp.BLAS_dsymm_d_s-f2c.c $@ && rm -f tmp.BLAS_dsymm_d_s-f2c.c*

$(SRC_PATH)/BLAS_dsymm_s_d-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="d, s, d" ../f2c-bridge.m4 > tmp.BLAS_dsymm_s_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymm_s_d-f2c.c && \
    mv tmp.BLAS_dsymm_s_d-f2c.c $@ && rm -f tmp.BLAS_dsymm_s_d-f2c.c*

$(SRC_PATH)/BLAS_dsymm_s_s-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="d, s, s" ../f2c-bridge.m4 > tmp.BLAS_dsymm_s_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymm_s_s-f2c.c && \
    mv tmp.BLAS_dsymm_s_s-f2c.c $@ && rm -f tmp.BLAS_dsymm_s_s-f2c.c*

$(SRC_PATH)/BLAS_zsymm_z_c-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="z, z, c" ../f2c-bridge.m4 > tmp.BLAS_zsymm_z_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm_z_c-f2c.c && \
    mv tmp.BLAS_zsymm_z_c-f2c.c $@ && rm -f tmp.BLAS_zsymm_z_c-f2c.c*

$(SRC_PATH)/BLAS_zsymm_c_z-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="z, c, z" ../f2c-bridge.m4 > tmp.BLAS_zsymm_c_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm_c_z-f2c.c && \
    mv tmp.BLAS_zsymm_c_z-f2c.c $@ && rm -f tmp.BLAS_zsymm_c_z-f2c.c*

$(SRC_PATH)/BLAS_zsymm_c_c-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="z, c, c" ../f2c-bridge.m4 > tmp.BLAS_zsymm_c_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm_c_c-f2c.c && \
    mv tmp.BLAS_zsymm_c_c-f2c.c $@ && rm -f tmp.BLAS_zsymm_c_c-f2c.c*

$(SRC_PATH)/BLAS_csymm_c_s-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="c, c, s" ../f2c-bridge.m4 > tmp.BLAS_csymm_c_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymm_c_s-f2c.c && \
    mv tmp.BLAS_csymm_c_s-f2c.c $@ && rm -f tmp.BLAS_csymm_c_s-f2c.c*

$(SRC_PATH)/BLAS_csymm_s_c-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="c, s, c" ../f2c-bridge.m4 > tmp.BLAS_csymm_s_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymm_s_c-f2c.c && \
    mv tmp.BLAS_csymm_s_c-f2c.c $@ && rm -f tmp.BLAS_csymm_s_c-f2c.c*

$(SRC_PATH)/BLAS_csymm_s_s-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="c, s, s" ../f2c-bridge.m4 > tmp.BLAS_csymm_s_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymm_s_s-f2c.c && \
    mv tmp.BLAS_csymm_s_s-f2c.c $@ && rm -f tmp.BLAS_csymm_s_s-f2c.c*

$(SRC_PATH)/BLAS_zsymm_z_d-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="z, z, d" ../f2c-bridge.m4 > tmp.BLAS_zsymm_z_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm_z_d-f2c.c && \
    mv tmp.BLAS_zsymm_z_d-f2c.c $@ && rm -f tmp.BLAS_zsymm_z_d-f2c.c*

$(SRC_PATH)/BLAS_zsymm_d_z-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="z, d, z" ../f2c-bridge.m4 > tmp.BLAS_zsymm_d_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm_d_z-f2c.c && \
    mv tmp.BLAS_zsymm_d_z-f2c.c $@ && rm -f tmp.BLAS_zsymm_d_z-f2c.c*

$(SRC_PATH)/BLAS_zsymm_d_d-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="z, d, d" ../f2c-bridge.m4 > tmp.BLAS_zsymm_d_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm_d_d-f2c.c && \
    mv tmp.BLAS_zsymm_d_d-f2c.c $@ && rm -f tmp.BLAS_zsymm_d_d-f2c.c*

$(SRC_PATH)/BLAS_ssymm_x-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="s, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_ssymm_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ssymm_x-f2c.c && \
    mv tmp.BLAS_ssymm_x-f2c.c $@ && rm -f tmp.BLAS_ssymm_x-f2c.c*

$(SRC_PATH)/BLAS_dsymm_x-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="d, d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dsymm_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymm_x-f2c.c && \
    mv tmp.BLAS_dsymm_x-f2c.c $@ && rm -f tmp.BLAS_dsymm_x-f2c.c*

$(SRC_PATH)/BLAS_csymm_x-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="c, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_csymm_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymm_x-f2c.c && \
    mv tmp.BLAS_csymm_x-f2c.c $@ && rm -f tmp.BLAS_csymm_x-f2c.c*

$(SRC_PATH)/BLAS_zsymm_x-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="z, z, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zsymm_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm_x-f2c.c && \
    mv tmp.BLAS_zsymm_x-f2c.c $@ && rm -f tmp.BLAS_zsymm_x-f2c.c*

$(SRC_PATH)/BLAS_dsymm_d_s_x-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="d, d, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dsymm_d_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymm_d_s_x-f2c.c && \
    mv tmp.BLAS_dsymm_d_s_x-f2c.c $@ && rm -f tmp.BLAS_dsymm_d_s_x-f2c.c*

$(SRC_PATH)/BLAS_dsymm_s_d_x-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="d, s, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dsymm_s_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymm_s_d_x-f2c.c && \
    mv tmp.BLAS_dsymm_s_d_x-f2c.c $@ && rm -f tmp.BLAS_dsymm_s_d_x-f2c.c*

$(SRC_PATH)/BLAS_dsymm_s_s_x-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="d, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dsymm_s_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsymm_s_s_x-f2c.c && \
    mv tmp.BLAS_dsymm_s_s_x-f2c.c $@ && rm -f tmp.BLAS_dsymm_s_s_x-f2c.c*

$(SRC_PATH)/BLAS_zsymm_z_c_x-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="z, z, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zsymm_z_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm_z_c_x-f2c.c && \
    mv tmp.BLAS_zsymm_z_c_x-f2c.c $@ && rm -f tmp.BLAS_zsymm_z_c_x-f2c.c*

$(SRC_PATH)/BLAS_zsymm_c_z_x-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="z, c, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zsymm_c_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm_c_z_x-f2c.c && \
    mv tmp.BLAS_zsymm_c_z_x-f2c.c $@ && rm -f tmp.BLAS_zsymm_c_z_x-f2c.c*

$(SRC_PATH)/BLAS_zsymm_c_c_x-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="z, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zsymm_c_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm_c_c_x-f2c.c && \
    mv tmp.BLAS_zsymm_c_c_x-f2c.c $@ && rm -f tmp.BLAS_zsymm_c_c_x-f2c.c*

$(SRC_PATH)/BLAS_csymm_c_s_x-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="c, c, s, _x" ../f2c-bridge.m4 > tmp.BLAS_csymm_c_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymm_c_s_x-f2c.c && \
    mv tmp.BLAS_csymm_c_s_x-f2c.c $@ && rm -f tmp.BLAS_csymm_c_s_x-f2c.c*

$(SRC_PATH)/BLAS_csymm_s_c_x-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="c, s, c, _x" ../f2c-bridge.m4 > tmp.BLAS_csymm_s_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymm_s_c_x-f2c.c && \
    mv tmp.BLAS_csymm_s_c_x-f2c.c $@ && rm -f tmp.BLAS_csymm_s_c_x-f2c.c*

$(SRC_PATH)/BLAS_csymm_s_s_x-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="c, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_csymm_s_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csymm_s_s_x-f2c.c && \
    mv tmp.BLAS_csymm_s_s_x-f2c.c $@ && rm -f tmp.BLAS_csymm_s_s_x-f2c.c*

$(SRC_PATH)/BLAS_zsymm_z_d_x-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="z, z, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zsymm_z_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm_z_d_x-f2c.c && \
    mv tmp.BLAS_zsymm_z_d_x-f2c.c $@ && rm -f tmp.BLAS_zsymm_z_d_x-f2c.c*

$(SRC_PATH)/BLAS_zsymm_d_z_x-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="z, d, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zsymm_d_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm_d_z_x-f2c.c && \
    mv tmp.BLAS_zsymm_d_z_x-f2c.c $@ && rm -f tmp.BLAS_zsymm_d_z_x-f2c.c*

$(SRC_PATH)/BLAS_zsymm_d_d_x-f2c.c: ../f2c-bridge.m4 symm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=symm -DRARG="z, d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zsymm_d_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsymm_d_d_x-f2c.c && \
    mv tmp.BLAS_zsymm_d_d_x-f2c.c $@ && rm -f tmp.BLAS_zsymm_d_d_x-f2c.c*


symm.h: ../cblas.m4 ../header.template symm-common.m4
	sed -e "s/routine/symm/g" -e "s/ROUTINE/SYMM/g" \
  ../header.template | $(M4) $(M4_OPTS) $(M4_PATH) - >$@

source-clean:
	rm -f $(SYMM_SRCS) symm.h

.PHONY: symm-source symm-header source-clean

