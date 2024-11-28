include ../../make.conf
include ../../$(MAKEINC)


SRC_PATH = ../../src/dot
M4_PATH = -I. -I..

DOT_SRCS = \
	$(SRC_PATH)/BLAS_sdot.c\
	$(SRC_PATH)/BLAS_ddot.c\
	$(SRC_PATH)/BLAS_cdot.c\
	$(SRC_PATH)/BLAS_zdot.c\
	$(SRC_PATH)/BLAS_ddot_d_s.c\
	$(SRC_PATH)/BLAS_ddot_s_d.c\
	$(SRC_PATH)/BLAS_ddot_s_s.c\
	$(SRC_PATH)/BLAS_zdot_z_c.c\
	$(SRC_PATH)/BLAS_zdot_c_z.c\
	$(SRC_PATH)/BLAS_zdot_c_c.c\
	$(SRC_PATH)/BLAS_cdot_c_s.c\
	$(SRC_PATH)/BLAS_cdot_s_c.c\
	$(SRC_PATH)/BLAS_cdot_s_s.c\
	$(SRC_PATH)/BLAS_zdot_z_d.c\
	$(SRC_PATH)/BLAS_zdot_d_z.c\
	$(SRC_PATH)/BLAS_zdot_d_d.c\
	$(SRC_PATH)/BLAS_sdot_x.c\
	$(SRC_PATH)/BLAS_ddot_x.c\
	$(SRC_PATH)/BLAS_cdot_x.c\
	$(SRC_PATH)/BLAS_zdot_x.c\
	$(SRC_PATH)/BLAS_ddot_d_s_x.c\
	$(SRC_PATH)/BLAS_ddot_s_d_x.c\
	$(SRC_PATH)/BLAS_ddot_s_s_x.c\
	$(SRC_PATH)/BLAS_zdot_z_c_x.c\
	$(SRC_PATH)/BLAS_zdot_c_z_x.c\
	$(SRC_PATH)/BLAS_zdot_c_c_x.c\
	$(SRC_PATH)/BLAS_cdot_c_s_x.c\
	$(SRC_PATH)/BLAS_cdot_s_c_x.c\
	$(SRC_PATH)/BLAS_cdot_s_s_x.c\
	$(SRC_PATH)/BLAS_zdot_z_d_x.c\
	$(SRC_PATH)/BLAS_zdot_d_z_x.c\
	$(SRC_PATH)/BLAS_zdot_d_d_x.c \
	\
	$(SRC_PATH)/BLAS_sdot-f2c.c\
	$(SRC_PATH)/BLAS_ddot-f2c.c\
	$(SRC_PATH)/BLAS_cdot-f2c.c\
	$(SRC_PATH)/BLAS_zdot-f2c.c\
	$(SRC_PATH)/BLAS_ddot_d_s-f2c.c\
	$(SRC_PATH)/BLAS_ddot_s_d-f2c.c\
	$(SRC_PATH)/BLAS_ddot_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zdot_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zdot_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zdot_c_c-f2c.c\
	$(SRC_PATH)/BLAS_cdot_c_s-f2c.c\
	$(SRC_PATH)/BLAS_cdot_s_c-f2c.c\
	$(SRC_PATH)/BLAS_cdot_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zdot_z_d-f2c.c\
	$(SRC_PATH)/BLAS_zdot_d_z-f2c.c\
	$(SRC_PATH)/BLAS_zdot_d_d-f2c.c\
	$(SRC_PATH)/BLAS_sdot_x-f2c.c\
	$(SRC_PATH)/BLAS_ddot_x-f2c.c\
	$(SRC_PATH)/BLAS_cdot_x-f2c.c\
	$(SRC_PATH)/BLAS_zdot_x-f2c.c\
	$(SRC_PATH)/BLAS_ddot_d_s_x-f2c.c\
	$(SRC_PATH)/BLAS_ddot_s_d_x-f2c.c\
	$(SRC_PATH)/BLAS_ddot_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zdot_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zdot_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zdot_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cdot_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_cdot_s_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cdot_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zdot_z_d_x-f2c.c\
	$(SRC_PATH)/BLAS_zdot_d_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zdot_d_d_x-f2c.c

dot-source: $(DOT_SRCS)

dot-header: dot.h

$(SRC_PATH)/BLAS_sdot.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(s, s, s)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_sdot.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sdot.c && \
    mv tmp.BLAS_sdot.c $@ && rm -f tmp.BLAS_sdot.c*

$(SRC_PATH)/BLAS_ddot.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(d, d, d)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ddot.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ddot.c && \
    mv tmp.BLAS_ddot.c $@ && rm -f tmp.BLAS_ddot.c*

$(SRC_PATH)/BLAS_cdot.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(c, c, c)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cdot.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cdot.c && \
    mv tmp.BLAS_cdot.c $@ && rm -f tmp.BLAS_cdot.c*

$(SRC_PATH)/BLAS_zdot.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(z, z, z)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zdot.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot.c && \
    mv tmp.BLAS_zdot.c $@ && rm -f tmp.BLAS_zdot.c*

$(SRC_PATH)/BLAS_ddot_d_s.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(d, d, s)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ddot_d_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ddot_d_s.c && \
    mv tmp.BLAS_ddot_d_s.c $@ && rm -f tmp.BLAS_ddot_d_s.c*

$(SRC_PATH)/BLAS_ddot_s_d.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(d, s, d)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ddot_s_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ddot_s_d.c && \
    mv tmp.BLAS_ddot_s_d.c $@ && rm -f tmp.BLAS_ddot_s_d.c*

$(SRC_PATH)/BLAS_ddot_s_s.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(d, s, s)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ddot_s_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ddot_s_s.c && \
    mv tmp.BLAS_ddot_s_s.c $@ && rm -f tmp.BLAS_ddot_s_s.c*

$(SRC_PATH)/BLAS_zdot_z_c.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(z, z, c)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zdot_z_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot_z_c.c && \
    mv tmp.BLAS_zdot_z_c.c $@ && rm -f tmp.BLAS_zdot_z_c.c*

$(SRC_PATH)/BLAS_zdot_c_z.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(z, c, z)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zdot_c_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot_c_z.c && \
    mv tmp.BLAS_zdot_c_z.c $@ && rm -f tmp.BLAS_zdot_c_z.c*

$(SRC_PATH)/BLAS_zdot_c_c.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(z, c, c)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zdot_c_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot_c_c.c && \
    mv tmp.BLAS_zdot_c_c.c $@ && rm -f tmp.BLAS_zdot_c_c.c*

$(SRC_PATH)/BLAS_cdot_c_s.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(c, c, s)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cdot_c_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cdot_c_s.c && \
    mv tmp.BLAS_cdot_c_s.c $@ && rm -f tmp.BLAS_cdot_c_s.c*

$(SRC_PATH)/BLAS_cdot_s_c.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(c, s, c)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cdot_s_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cdot_s_c.c && \
    mv tmp.BLAS_cdot_s_c.c $@ && rm -f tmp.BLAS_cdot_s_c.c*

$(SRC_PATH)/BLAS_cdot_s_s.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(c, s, s)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cdot_s_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cdot_s_s.c && \
    mv tmp.BLAS_cdot_s_s.c $@ && rm -f tmp.BLAS_cdot_s_s.c*

$(SRC_PATH)/BLAS_zdot_z_d.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(z, z, d)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zdot_z_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot_z_d.c && \
    mv tmp.BLAS_zdot_z_d.c $@ && rm -f tmp.BLAS_zdot_z_d.c*

$(SRC_PATH)/BLAS_zdot_d_z.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(z, d, z)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zdot_d_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot_d_z.c && \
    mv tmp.BLAS_zdot_d_z.c $@ && rm -f tmp.BLAS_zdot_d_z.c*

$(SRC_PATH)/BLAS_zdot_d_d.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(z, d, d)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zdot_d_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot_d_d.c && \
    mv tmp.BLAS_zdot_d_d.c $@ && rm -f tmp.BLAS_zdot_d_d.c*

$(SRC_PATH)/BLAS_sdot_x.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(s, s, s, _x)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_sdot_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sdot_x.c && \
    mv tmp.BLAS_sdot_x.c $@ && rm -f tmp.BLAS_sdot_x.c*

$(SRC_PATH)/BLAS_ddot_x.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(d, d, d, _x)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ddot_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ddot_x.c && \
    mv tmp.BLAS_ddot_x.c $@ && rm -f tmp.BLAS_ddot_x.c*

$(SRC_PATH)/BLAS_cdot_x.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(c, c, c, _x)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cdot_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cdot_x.c && \
    mv tmp.BLAS_cdot_x.c $@ && rm -f tmp.BLAS_cdot_x.c*

$(SRC_PATH)/BLAS_zdot_x.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(z, z, z, _x)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zdot_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot_x.c && \
    mv tmp.BLAS_zdot_x.c $@ && rm -f tmp.BLAS_zdot_x.c*

$(SRC_PATH)/BLAS_ddot_d_s_x.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(d, d, s, _x)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ddot_d_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ddot_d_s_x.c && \
    mv tmp.BLAS_ddot_d_s_x.c $@ && rm -f tmp.BLAS_ddot_d_s_x.c*

$(SRC_PATH)/BLAS_ddot_s_d_x.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(d, s, d, _x)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ddot_s_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ddot_s_d_x.c && \
    mv tmp.BLAS_ddot_s_d_x.c $@ && rm -f tmp.BLAS_ddot_s_d_x.c*

$(SRC_PATH)/BLAS_ddot_s_s_x.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(d, s, s, _x)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ddot_s_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ddot_s_s_x.c && \
    mv tmp.BLAS_ddot_s_s_x.c $@ && rm -f tmp.BLAS_ddot_s_s_x.c*

$(SRC_PATH)/BLAS_zdot_z_c_x.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(z, z, c, _x)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zdot_z_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot_z_c_x.c && \
    mv tmp.BLAS_zdot_z_c_x.c $@ && rm -f tmp.BLAS_zdot_z_c_x.c*

$(SRC_PATH)/BLAS_zdot_c_z_x.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(z, c, z, _x)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zdot_c_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot_c_z_x.c && \
    mv tmp.BLAS_zdot_c_z_x.c $@ && rm -f tmp.BLAS_zdot_c_z_x.c*

$(SRC_PATH)/BLAS_zdot_c_c_x.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(z, c, c, _x)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zdot_c_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot_c_c_x.c && \
    mv tmp.BLAS_zdot_c_c_x.c $@ && rm -f tmp.BLAS_zdot_c_c_x.c*

$(SRC_PATH)/BLAS_cdot_c_s_x.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(c, c, s, _x)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cdot_c_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cdot_c_s_x.c && \
    mv tmp.BLAS_cdot_c_s_x.c $@ && rm -f tmp.BLAS_cdot_c_s_x.c*

$(SRC_PATH)/BLAS_cdot_s_c_x.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(c, s, c, _x)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cdot_s_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cdot_s_c_x.c && \
    mv tmp.BLAS_cdot_s_c_x.c $@ && rm -f tmp.BLAS_cdot_s_c_x.c*

$(SRC_PATH)/BLAS_cdot_s_s_x.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(c, s, s, _x)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cdot_s_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cdot_s_s_x.c && \
    mv tmp.BLAS_cdot_s_s_x.c $@ && rm -f tmp.BLAS_cdot_s_s_x.c*

$(SRC_PATH)/BLAS_zdot_z_d_x.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(z, z, d, _x)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zdot_z_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot_z_d_x.c && \
    mv tmp.BLAS_zdot_z_d_x.c $@ && rm -f tmp.BLAS_zdot_z_d_x.c*

$(SRC_PATH)/BLAS_zdot_d_z_x.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(z, d, z, _x)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zdot_d_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot_d_z_x.c && \
    mv tmp.BLAS_zdot_d_z_x.c $@ && rm -f tmp.BLAS_zdot_d_z_x.c*

$(SRC_PATH)/BLAS_zdot_d_d_x.c: dot.m4 dot-common.m4 ../cblas.m4
	echo "DOT(z, d, d, _x)" | cat dot.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zdot_d_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot_d_d_x.c && \
    mv tmp.BLAS_zdot_d_d_x.c $@ && rm -f tmp.BLAS_zdot_d_d_x.c*



$(SRC_PATH)/BLAS_sdot-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="s, s, s" ../f2c-bridge.m4 > tmp.BLAS_sdot-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sdot-f2c.c && \
    mv tmp.BLAS_sdot-f2c.c $@ && rm -f tmp.BLAS_sdot-f2c.c*

$(SRC_PATH)/BLAS_ddot-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="d, d, d" ../f2c-bridge.m4 > tmp.BLAS_ddot-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ddot-f2c.c && \
    mv tmp.BLAS_ddot-f2c.c $@ && rm -f tmp.BLAS_ddot-f2c.c*

$(SRC_PATH)/BLAS_cdot-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="c, c, c" ../f2c-bridge.m4 > tmp.BLAS_cdot-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cdot-f2c.c && \
    mv tmp.BLAS_cdot-f2c.c $@ && rm -f tmp.BLAS_cdot-f2c.c*

$(SRC_PATH)/BLAS_zdot-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="z, z, z" ../f2c-bridge.m4 > tmp.BLAS_zdot-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot-f2c.c && \
    mv tmp.BLAS_zdot-f2c.c $@ && rm -f tmp.BLAS_zdot-f2c.c*

$(SRC_PATH)/BLAS_ddot_d_s-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="d, d, s" ../f2c-bridge.m4 > tmp.BLAS_ddot_d_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ddot_d_s-f2c.c && \
    mv tmp.BLAS_ddot_d_s-f2c.c $@ && rm -f tmp.BLAS_ddot_d_s-f2c.c*

$(SRC_PATH)/BLAS_ddot_s_d-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="d, s, d" ../f2c-bridge.m4 > tmp.BLAS_ddot_s_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ddot_s_d-f2c.c && \
    mv tmp.BLAS_ddot_s_d-f2c.c $@ && rm -f tmp.BLAS_ddot_s_d-f2c.c*

$(SRC_PATH)/BLAS_ddot_s_s-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="d, s, s" ../f2c-bridge.m4 > tmp.BLAS_ddot_s_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ddot_s_s-f2c.c && \
    mv tmp.BLAS_ddot_s_s-f2c.c $@ && rm -f tmp.BLAS_ddot_s_s-f2c.c*

$(SRC_PATH)/BLAS_zdot_z_c-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="z, z, c" ../f2c-bridge.m4 > tmp.BLAS_zdot_z_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot_z_c-f2c.c && \
    mv tmp.BLAS_zdot_z_c-f2c.c $@ && rm -f tmp.BLAS_zdot_z_c-f2c.c*

$(SRC_PATH)/BLAS_zdot_c_z-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="z, c, z" ../f2c-bridge.m4 > tmp.BLAS_zdot_c_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot_c_z-f2c.c && \
    mv tmp.BLAS_zdot_c_z-f2c.c $@ && rm -f tmp.BLAS_zdot_c_z-f2c.c*

$(SRC_PATH)/BLAS_zdot_c_c-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="z, c, c" ../f2c-bridge.m4 > tmp.BLAS_zdot_c_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot_c_c-f2c.c && \
    mv tmp.BLAS_zdot_c_c-f2c.c $@ && rm -f tmp.BLAS_zdot_c_c-f2c.c*

$(SRC_PATH)/BLAS_cdot_c_s-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="c, c, s" ../f2c-bridge.m4 > tmp.BLAS_cdot_c_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cdot_c_s-f2c.c && \
    mv tmp.BLAS_cdot_c_s-f2c.c $@ && rm -f tmp.BLAS_cdot_c_s-f2c.c*

$(SRC_PATH)/BLAS_cdot_s_c-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="c, s, c" ../f2c-bridge.m4 > tmp.BLAS_cdot_s_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cdot_s_c-f2c.c && \
    mv tmp.BLAS_cdot_s_c-f2c.c $@ && rm -f tmp.BLAS_cdot_s_c-f2c.c*

$(SRC_PATH)/BLAS_cdot_s_s-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="c, s, s" ../f2c-bridge.m4 > tmp.BLAS_cdot_s_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cdot_s_s-f2c.c && \
    mv tmp.BLAS_cdot_s_s-f2c.c $@ && rm -f tmp.BLAS_cdot_s_s-f2c.c*

$(SRC_PATH)/BLAS_zdot_z_d-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="z, z, d" ../f2c-bridge.m4 > tmp.BLAS_zdot_z_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot_z_d-f2c.c && \
    mv tmp.BLAS_zdot_z_d-f2c.c $@ && rm -f tmp.BLAS_zdot_z_d-f2c.c*

$(SRC_PATH)/BLAS_zdot_d_z-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="z, d, z" ../f2c-bridge.m4 > tmp.BLAS_zdot_d_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot_d_z-f2c.c && \
    mv tmp.BLAS_zdot_d_z-f2c.c $@ && rm -f tmp.BLAS_zdot_d_z-f2c.c*

$(SRC_PATH)/BLAS_zdot_d_d-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="z, d, d" ../f2c-bridge.m4 > tmp.BLAS_zdot_d_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot_d_d-f2c.c && \
    mv tmp.BLAS_zdot_d_d-f2c.c $@ && rm -f tmp.BLAS_zdot_d_d-f2c.c*

$(SRC_PATH)/BLAS_sdot_x-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="s, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_sdot_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sdot_x-f2c.c && \
    mv tmp.BLAS_sdot_x-f2c.c $@ && rm -f tmp.BLAS_sdot_x-f2c.c*

$(SRC_PATH)/BLAS_ddot_x-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="d, d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_ddot_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ddot_x-f2c.c && \
    mv tmp.BLAS_ddot_x-f2c.c $@ && rm -f tmp.BLAS_ddot_x-f2c.c*

$(SRC_PATH)/BLAS_cdot_x-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="c, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_cdot_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cdot_x-f2c.c && \
    mv tmp.BLAS_cdot_x-f2c.c $@ && rm -f tmp.BLAS_cdot_x-f2c.c*

$(SRC_PATH)/BLAS_zdot_x-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="z, z, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zdot_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot_x-f2c.c && \
    mv tmp.BLAS_zdot_x-f2c.c $@ && rm -f tmp.BLAS_zdot_x-f2c.c*

$(SRC_PATH)/BLAS_ddot_d_s_x-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="d, d, s, _x" ../f2c-bridge.m4 > tmp.BLAS_ddot_d_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ddot_d_s_x-f2c.c && \
    mv tmp.BLAS_ddot_d_s_x-f2c.c $@ && rm -f tmp.BLAS_ddot_d_s_x-f2c.c*

$(SRC_PATH)/BLAS_ddot_s_d_x-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="d, s, d, _x" ../f2c-bridge.m4 > tmp.BLAS_ddot_s_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ddot_s_d_x-f2c.c && \
    mv tmp.BLAS_ddot_s_d_x-f2c.c $@ && rm -f tmp.BLAS_ddot_s_d_x-f2c.c*

$(SRC_PATH)/BLAS_ddot_s_s_x-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="d, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_ddot_s_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ddot_s_s_x-f2c.c && \
    mv tmp.BLAS_ddot_s_s_x-f2c.c $@ && rm -f tmp.BLAS_ddot_s_s_x-f2c.c*

$(SRC_PATH)/BLAS_zdot_z_c_x-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="z, z, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zdot_z_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot_z_c_x-f2c.c && \
    mv tmp.BLAS_zdot_z_c_x-f2c.c $@ && rm -f tmp.BLAS_zdot_z_c_x-f2c.c*

$(SRC_PATH)/BLAS_zdot_c_z_x-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="z, c, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zdot_c_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot_c_z_x-f2c.c && \
    mv tmp.BLAS_zdot_c_z_x-f2c.c $@ && rm -f tmp.BLAS_zdot_c_z_x-f2c.c*

$(SRC_PATH)/BLAS_zdot_c_c_x-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="z, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zdot_c_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot_c_c_x-f2c.c && \
    mv tmp.BLAS_zdot_c_c_x-f2c.c $@ && rm -f tmp.BLAS_zdot_c_c_x-f2c.c*

$(SRC_PATH)/BLAS_cdot_c_s_x-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="c, c, s, _x" ../f2c-bridge.m4 > tmp.BLAS_cdot_c_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cdot_c_s_x-f2c.c && \
    mv tmp.BLAS_cdot_c_s_x-f2c.c $@ && rm -f tmp.BLAS_cdot_c_s_x-f2c.c*

$(SRC_PATH)/BLAS_cdot_s_c_x-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="c, s, c, _x" ../f2c-bridge.m4 > tmp.BLAS_cdot_s_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cdot_s_c_x-f2c.c && \
    mv tmp.BLAS_cdot_s_c_x-f2c.c $@ && rm -f tmp.BLAS_cdot_s_c_x-f2c.c*

$(SRC_PATH)/BLAS_cdot_s_s_x-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="c, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_cdot_s_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cdot_s_s_x-f2c.c && \
    mv tmp.BLAS_cdot_s_s_x-f2c.c $@ && rm -f tmp.BLAS_cdot_s_s_x-f2c.c*

$(SRC_PATH)/BLAS_zdot_z_d_x-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="z, z, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zdot_z_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot_z_d_x-f2c.c && \
    mv tmp.BLAS_zdot_z_d_x-f2c.c $@ && rm -f tmp.BLAS_zdot_z_d_x-f2c.c*

$(SRC_PATH)/BLAS_zdot_d_z_x-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="z, d, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zdot_d_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot_d_z_x-f2c.c && \
    mv tmp.BLAS_zdot_d_z_x-f2c.c $@ && rm -f tmp.BLAS_zdot_d_z_x-f2c.c*

$(SRC_PATH)/BLAS_zdot_d_d_x-f2c.c: ../f2c-bridge.m4 dot-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=dot -DRARG="z, d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zdot_d_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zdot_d_d_x-f2c.c && \
    mv tmp.BLAS_zdot_d_d_x-f2c.c $@ && rm -f tmp.BLAS_zdot_d_d_x-f2c.c*


dot.h: ../cblas.m4 ../header.template dot-common.m4
	sed -e "s/routine/dot/g" -e "s/ROUTINE/DOT/g" \
  ../header.template | $(M4) $(M4_OPTS) $(M4_PATH) - >$@

source-clean:
	rm -f $(DOT_SRCS) dot.h

.PHONY: dot-source dot-header source-clean

