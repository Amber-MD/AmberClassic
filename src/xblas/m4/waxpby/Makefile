include ../../make.conf
include ../../$(MAKEINC)


SRC_PATH = ../../src/waxpby
M4_PATH = -I. -I..

WAXPBY_SRCS = \
	$(SRC_PATH)/BLAS_swaxpby.c\
	$(SRC_PATH)/BLAS_dwaxpby.c\
	$(SRC_PATH)/BLAS_cwaxpby.c\
	$(SRC_PATH)/BLAS_zwaxpby.c\
	$(SRC_PATH)/BLAS_dwaxpby_d_s.c\
	$(SRC_PATH)/BLAS_dwaxpby_s_d.c\
	$(SRC_PATH)/BLAS_dwaxpby_s_s.c\
	$(SRC_PATH)/BLAS_zwaxpby_z_c.c\
	$(SRC_PATH)/BLAS_zwaxpby_c_z.c\
	$(SRC_PATH)/BLAS_zwaxpby_c_c.c\
	$(SRC_PATH)/BLAS_cwaxpby_c_s.c\
	$(SRC_PATH)/BLAS_cwaxpby_s_c.c\
	$(SRC_PATH)/BLAS_cwaxpby_s_s.c\
	$(SRC_PATH)/BLAS_zwaxpby_z_d.c\
	$(SRC_PATH)/BLAS_zwaxpby_d_z.c\
	$(SRC_PATH)/BLAS_zwaxpby_d_d.c\
	$(SRC_PATH)/BLAS_swaxpby_x.c\
	$(SRC_PATH)/BLAS_dwaxpby_x.c\
	$(SRC_PATH)/BLAS_cwaxpby_x.c\
	$(SRC_PATH)/BLAS_zwaxpby_x.c\
	$(SRC_PATH)/BLAS_dwaxpby_d_s_x.c\
	$(SRC_PATH)/BLAS_dwaxpby_s_d_x.c\
	$(SRC_PATH)/BLAS_dwaxpby_s_s_x.c\
	$(SRC_PATH)/BLAS_zwaxpby_z_c_x.c\
	$(SRC_PATH)/BLAS_zwaxpby_c_z_x.c\
	$(SRC_PATH)/BLAS_zwaxpby_c_c_x.c\
	$(SRC_PATH)/BLAS_cwaxpby_c_s_x.c\
	$(SRC_PATH)/BLAS_cwaxpby_s_c_x.c\
	$(SRC_PATH)/BLAS_cwaxpby_s_s_x.c\
	$(SRC_PATH)/BLAS_zwaxpby_z_d_x.c\
	$(SRC_PATH)/BLAS_zwaxpby_d_z_x.c\
	$(SRC_PATH)/BLAS_zwaxpby_d_d_x.c \
	\
	$(SRC_PATH)/BLAS_swaxpby-f2c.c\
	$(SRC_PATH)/BLAS_dwaxpby-f2c.c\
	$(SRC_PATH)/BLAS_cwaxpby-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby-f2c.c\
	$(SRC_PATH)/BLAS_dwaxpby_d_s-f2c.c\
	$(SRC_PATH)/BLAS_dwaxpby_s_d-f2c.c\
	$(SRC_PATH)/BLAS_dwaxpby_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby_c_c-f2c.c\
	$(SRC_PATH)/BLAS_cwaxpby_c_s-f2c.c\
	$(SRC_PATH)/BLAS_cwaxpby_s_c-f2c.c\
	$(SRC_PATH)/BLAS_cwaxpby_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby_z_d-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby_d_z-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby_d_d-f2c.c\
	$(SRC_PATH)/BLAS_swaxpby_x-f2c.c\
	$(SRC_PATH)/BLAS_dwaxpby_x-f2c.c\
	$(SRC_PATH)/BLAS_cwaxpby_x-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby_x-f2c.c\
	$(SRC_PATH)/BLAS_dwaxpby_d_s_x-f2c.c\
	$(SRC_PATH)/BLAS_dwaxpby_s_d_x-f2c.c\
	$(SRC_PATH)/BLAS_dwaxpby_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cwaxpby_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_cwaxpby_s_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cwaxpby_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby_z_d_x-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby_d_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zwaxpby_d_d_x-f2c.c

waxpby-source: $(WAXPBY_SRCS)

waxpby-header: waxpby.h

$(SRC_PATH)/BLAS_swaxpby.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(s, s, s)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_swaxpby.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_swaxpby.c && \
    mv tmp.BLAS_swaxpby.c $@ && rm -f tmp.BLAS_swaxpby.c*

$(SRC_PATH)/BLAS_dwaxpby.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(d, d, d)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dwaxpby.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dwaxpby.c && \
    mv tmp.BLAS_dwaxpby.c $@ && rm -f tmp.BLAS_dwaxpby.c*

$(SRC_PATH)/BLAS_cwaxpby.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(c, c, c)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cwaxpby.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cwaxpby.c && \
    mv tmp.BLAS_cwaxpby.c $@ && rm -f tmp.BLAS_cwaxpby.c*

$(SRC_PATH)/BLAS_zwaxpby.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(z, z, z)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zwaxpby.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby.c && \
    mv tmp.BLAS_zwaxpby.c $@ && rm -f tmp.BLAS_zwaxpby.c*

$(SRC_PATH)/BLAS_dwaxpby_d_s.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(d, d, s)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dwaxpby_d_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dwaxpby_d_s.c && \
    mv tmp.BLAS_dwaxpby_d_s.c $@ && rm -f tmp.BLAS_dwaxpby_d_s.c*

$(SRC_PATH)/BLAS_dwaxpby_s_d.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(d, s, d)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dwaxpby_s_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dwaxpby_s_d.c && \
    mv tmp.BLAS_dwaxpby_s_d.c $@ && rm -f tmp.BLAS_dwaxpby_s_d.c*

$(SRC_PATH)/BLAS_dwaxpby_s_s.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(d, s, s)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dwaxpby_s_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dwaxpby_s_s.c && \
    mv tmp.BLAS_dwaxpby_s_s.c $@ && rm -f tmp.BLAS_dwaxpby_s_s.c*

$(SRC_PATH)/BLAS_zwaxpby_z_c.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(z, z, c)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zwaxpby_z_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby_z_c.c && \
    mv tmp.BLAS_zwaxpby_z_c.c $@ && rm -f tmp.BLAS_zwaxpby_z_c.c*

$(SRC_PATH)/BLAS_zwaxpby_c_z.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(z, c, z)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zwaxpby_c_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby_c_z.c && \
    mv tmp.BLAS_zwaxpby_c_z.c $@ && rm -f tmp.BLAS_zwaxpby_c_z.c*

$(SRC_PATH)/BLAS_zwaxpby_c_c.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(z, c, c)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zwaxpby_c_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby_c_c.c && \
    mv tmp.BLAS_zwaxpby_c_c.c $@ && rm -f tmp.BLAS_zwaxpby_c_c.c*

$(SRC_PATH)/BLAS_cwaxpby_c_s.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(c, c, s)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cwaxpby_c_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cwaxpby_c_s.c && \
    mv tmp.BLAS_cwaxpby_c_s.c $@ && rm -f tmp.BLAS_cwaxpby_c_s.c*

$(SRC_PATH)/BLAS_cwaxpby_s_c.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(c, s, c)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cwaxpby_s_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cwaxpby_s_c.c && \
    mv tmp.BLAS_cwaxpby_s_c.c $@ && rm -f tmp.BLAS_cwaxpby_s_c.c*

$(SRC_PATH)/BLAS_cwaxpby_s_s.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(c, s, s)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cwaxpby_s_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cwaxpby_s_s.c && \
    mv tmp.BLAS_cwaxpby_s_s.c $@ && rm -f tmp.BLAS_cwaxpby_s_s.c*

$(SRC_PATH)/BLAS_zwaxpby_z_d.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(z, z, d)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zwaxpby_z_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby_z_d.c && \
    mv tmp.BLAS_zwaxpby_z_d.c $@ && rm -f tmp.BLAS_zwaxpby_z_d.c*

$(SRC_PATH)/BLAS_zwaxpby_d_z.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(z, d, z)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zwaxpby_d_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby_d_z.c && \
    mv tmp.BLAS_zwaxpby_d_z.c $@ && rm -f tmp.BLAS_zwaxpby_d_z.c*

$(SRC_PATH)/BLAS_zwaxpby_d_d.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(z, d, d)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zwaxpby_d_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby_d_d.c && \
    mv tmp.BLAS_zwaxpby_d_d.c $@ && rm -f tmp.BLAS_zwaxpby_d_d.c*

$(SRC_PATH)/BLAS_swaxpby_x.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(s, s, s, _x)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_swaxpby_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_swaxpby_x.c && \
    mv tmp.BLAS_swaxpby_x.c $@ && rm -f tmp.BLAS_swaxpby_x.c*

$(SRC_PATH)/BLAS_dwaxpby_x.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(d, d, d, _x)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dwaxpby_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dwaxpby_x.c && \
    mv tmp.BLAS_dwaxpby_x.c $@ && rm -f tmp.BLAS_dwaxpby_x.c*

$(SRC_PATH)/BLAS_cwaxpby_x.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(c, c, c, _x)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cwaxpby_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cwaxpby_x.c && \
    mv tmp.BLAS_cwaxpby_x.c $@ && rm -f tmp.BLAS_cwaxpby_x.c*

$(SRC_PATH)/BLAS_zwaxpby_x.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(z, z, z, _x)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zwaxpby_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby_x.c && \
    mv tmp.BLAS_zwaxpby_x.c $@ && rm -f tmp.BLAS_zwaxpby_x.c*

$(SRC_PATH)/BLAS_dwaxpby_d_s_x.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(d, d, s, _x)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dwaxpby_d_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dwaxpby_d_s_x.c && \
    mv tmp.BLAS_dwaxpby_d_s_x.c $@ && rm -f tmp.BLAS_dwaxpby_d_s_x.c*

$(SRC_PATH)/BLAS_dwaxpby_s_d_x.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(d, s, d, _x)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dwaxpby_s_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dwaxpby_s_d_x.c && \
    mv tmp.BLAS_dwaxpby_s_d_x.c $@ && rm -f tmp.BLAS_dwaxpby_s_d_x.c*

$(SRC_PATH)/BLAS_dwaxpby_s_s_x.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(d, s, s, _x)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dwaxpby_s_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dwaxpby_s_s_x.c && \
    mv tmp.BLAS_dwaxpby_s_s_x.c $@ && rm -f tmp.BLAS_dwaxpby_s_s_x.c*

$(SRC_PATH)/BLAS_zwaxpby_z_c_x.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(z, z, c, _x)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zwaxpby_z_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby_z_c_x.c && \
    mv tmp.BLAS_zwaxpby_z_c_x.c $@ && rm -f tmp.BLAS_zwaxpby_z_c_x.c*

$(SRC_PATH)/BLAS_zwaxpby_c_z_x.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(z, c, z, _x)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zwaxpby_c_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby_c_z_x.c && \
    mv tmp.BLAS_zwaxpby_c_z_x.c $@ && rm -f tmp.BLAS_zwaxpby_c_z_x.c*

$(SRC_PATH)/BLAS_zwaxpby_c_c_x.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(z, c, c, _x)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zwaxpby_c_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby_c_c_x.c && \
    mv tmp.BLAS_zwaxpby_c_c_x.c $@ && rm -f tmp.BLAS_zwaxpby_c_c_x.c*

$(SRC_PATH)/BLAS_cwaxpby_c_s_x.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(c, c, s, _x)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cwaxpby_c_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cwaxpby_c_s_x.c && \
    mv tmp.BLAS_cwaxpby_c_s_x.c $@ && rm -f tmp.BLAS_cwaxpby_c_s_x.c*

$(SRC_PATH)/BLAS_cwaxpby_s_c_x.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(c, s, c, _x)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cwaxpby_s_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cwaxpby_s_c_x.c && \
    mv tmp.BLAS_cwaxpby_s_c_x.c $@ && rm -f tmp.BLAS_cwaxpby_s_c_x.c*

$(SRC_PATH)/BLAS_cwaxpby_s_s_x.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(c, s, s, _x)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cwaxpby_s_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cwaxpby_s_s_x.c && \
    mv tmp.BLAS_cwaxpby_s_s_x.c $@ && rm -f tmp.BLAS_cwaxpby_s_s_x.c*

$(SRC_PATH)/BLAS_zwaxpby_z_d_x.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(z, z, d, _x)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zwaxpby_z_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby_z_d_x.c && \
    mv tmp.BLAS_zwaxpby_z_d_x.c $@ && rm -f tmp.BLAS_zwaxpby_z_d_x.c*

$(SRC_PATH)/BLAS_zwaxpby_d_z_x.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(z, d, z, _x)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zwaxpby_d_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby_d_z_x.c && \
    mv tmp.BLAS_zwaxpby_d_z_x.c $@ && rm -f tmp.BLAS_zwaxpby_d_z_x.c*

$(SRC_PATH)/BLAS_zwaxpby_d_d_x.c: waxpby.m4 waxpby-common.m4 ../cblas.m4
	echo "WAXPBY(z, d, d, _x)" | cat waxpby.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zwaxpby_d_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby_d_d_x.c && \
    mv tmp.BLAS_zwaxpby_d_d_x.c $@ && rm -f tmp.BLAS_zwaxpby_d_d_x.c*



$(SRC_PATH)/BLAS_swaxpby-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="s, s, s" ../f2c-bridge.m4 > tmp.BLAS_swaxpby-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_swaxpby-f2c.c && \
    mv tmp.BLAS_swaxpby-f2c.c $@ && rm -f tmp.BLAS_swaxpby-f2c.c*

$(SRC_PATH)/BLAS_dwaxpby-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="d, d, d" ../f2c-bridge.m4 > tmp.BLAS_dwaxpby-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dwaxpby-f2c.c && \
    mv tmp.BLAS_dwaxpby-f2c.c $@ && rm -f tmp.BLAS_dwaxpby-f2c.c*

$(SRC_PATH)/BLAS_cwaxpby-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="c, c, c" ../f2c-bridge.m4 > tmp.BLAS_cwaxpby-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cwaxpby-f2c.c && \
    mv tmp.BLAS_cwaxpby-f2c.c $@ && rm -f tmp.BLAS_cwaxpby-f2c.c*

$(SRC_PATH)/BLAS_zwaxpby-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="z, z, z" ../f2c-bridge.m4 > tmp.BLAS_zwaxpby-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby-f2c.c && \
    mv tmp.BLAS_zwaxpby-f2c.c $@ && rm -f tmp.BLAS_zwaxpby-f2c.c*

$(SRC_PATH)/BLAS_dwaxpby_d_s-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="d, d, s" ../f2c-bridge.m4 > tmp.BLAS_dwaxpby_d_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dwaxpby_d_s-f2c.c && \
    mv tmp.BLAS_dwaxpby_d_s-f2c.c $@ && rm -f tmp.BLAS_dwaxpby_d_s-f2c.c*

$(SRC_PATH)/BLAS_dwaxpby_s_d-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="d, s, d" ../f2c-bridge.m4 > tmp.BLAS_dwaxpby_s_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dwaxpby_s_d-f2c.c && \
    mv tmp.BLAS_dwaxpby_s_d-f2c.c $@ && rm -f tmp.BLAS_dwaxpby_s_d-f2c.c*

$(SRC_PATH)/BLAS_dwaxpby_s_s-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="d, s, s" ../f2c-bridge.m4 > tmp.BLAS_dwaxpby_s_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dwaxpby_s_s-f2c.c && \
    mv tmp.BLAS_dwaxpby_s_s-f2c.c $@ && rm -f tmp.BLAS_dwaxpby_s_s-f2c.c*

$(SRC_PATH)/BLAS_zwaxpby_z_c-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="z, z, c" ../f2c-bridge.m4 > tmp.BLAS_zwaxpby_z_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby_z_c-f2c.c && \
    mv tmp.BLAS_zwaxpby_z_c-f2c.c $@ && rm -f tmp.BLAS_zwaxpby_z_c-f2c.c*

$(SRC_PATH)/BLAS_zwaxpby_c_z-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="z, c, z" ../f2c-bridge.m4 > tmp.BLAS_zwaxpby_c_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby_c_z-f2c.c && \
    mv tmp.BLAS_zwaxpby_c_z-f2c.c $@ && rm -f tmp.BLAS_zwaxpby_c_z-f2c.c*

$(SRC_PATH)/BLAS_zwaxpby_c_c-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="z, c, c" ../f2c-bridge.m4 > tmp.BLAS_zwaxpby_c_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby_c_c-f2c.c && \
    mv tmp.BLAS_zwaxpby_c_c-f2c.c $@ && rm -f tmp.BLAS_zwaxpby_c_c-f2c.c*

$(SRC_PATH)/BLAS_cwaxpby_c_s-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="c, c, s" ../f2c-bridge.m4 > tmp.BLAS_cwaxpby_c_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cwaxpby_c_s-f2c.c && \
    mv tmp.BLAS_cwaxpby_c_s-f2c.c $@ && rm -f tmp.BLAS_cwaxpby_c_s-f2c.c*

$(SRC_PATH)/BLAS_cwaxpby_s_c-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="c, s, c" ../f2c-bridge.m4 > tmp.BLAS_cwaxpby_s_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cwaxpby_s_c-f2c.c && \
    mv tmp.BLAS_cwaxpby_s_c-f2c.c $@ && rm -f tmp.BLAS_cwaxpby_s_c-f2c.c*

$(SRC_PATH)/BLAS_cwaxpby_s_s-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="c, s, s" ../f2c-bridge.m4 > tmp.BLAS_cwaxpby_s_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cwaxpby_s_s-f2c.c && \
    mv tmp.BLAS_cwaxpby_s_s-f2c.c $@ && rm -f tmp.BLAS_cwaxpby_s_s-f2c.c*

$(SRC_PATH)/BLAS_zwaxpby_z_d-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="z, z, d" ../f2c-bridge.m4 > tmp.BLAS_zwaxpby_z_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby_z_d-f2c.c && \
    mv tmp.BLAS_zwaxpby_z_d-f2c.c $@ && rm -f tmp.BLAS_zwaxpby_z_d-f2c.c*

$(SRC_PATH)/BLAS_zwaxpby_d_z-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="z, d, z" ../f2c-bridge.m4 > tmp.BLAS_zwaxpby_d_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby_d_z-f2c.c && \
    mv tmp.BLAS_zwaxpby_d_z-f2c.c $@ && rm -f tmp.BLAS_zwaxpby_d_z-f2c.c*

$(SRC_PATH)/BLAS_zwaxpby_d_d-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="z, d, d" ../f2c-bridge.m4 > tmp.BLAS_zwaxpby_d_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby_d_d-f2c.c && \
    mv tmp.BLAS_zwaxpby_d_d-f2c.c $@ && rm -f tmp.BLAS_zwaxpby_d_d-f2c.c*

$(SRC_PATH)/BLAS_swaxpby_x-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="s, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_swaxpby_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_swaxpby_x-f2c.c && \
    mv tmp.BLAS_swaxpby_x-f2c.c $@ && rm -f tmp.BLAS_swaxpby_x-f2c.c*

$(SRC_PATH)/BLAS_dwaxpby_x-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="d, d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dwaxpby_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dwaxpby_x-f2c.c && \
    mv tmp.BLAS_dwaxpby_x-f2c.c $@ && rm -f tmp.BLAS_dwaxpby_x-f2c.c*

$(SRC_PATH)/BLAS_cwaxpby_x-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="c, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_cwaxpby_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cwaxpby_x-f2c.c && \
    mv tmp.BLAS_cwaxpby_x-f2c.c $@ && rm -f tmp.BLAS_cwaxpby_x-f2c.c*

$(SRC_PATH)/BLAS_zwaxpby_x-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="z, z, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zwaxpby_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby_x-f2c.c && \
    mv tmp.BLAS_zwaxpby_x-f2c.c $@ && rm -f tmp.BLAS_zwaxpby_x-f2c.c*

$(SRC_PATH)/BLAS_dwaxpby_d_s_x-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="d, d, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dwaxpby_d_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dwaxpby_d_s_x-f2c.c && \
    mv tmp.BLAS_dwaxpby_d_s_x-f2c.c $@ && rm -f tmp.BLAS_dwaxpby_d_s_x-f2c.c*

$(SRC_PATH)/BLAS_dwaxpby_s_d_x-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="d, s, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dwaxpby_s_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dwaxpby_s_d_x-f2c.c && \
    mv tmp.BLAS_dwaxpby_s_d_x-f2c.c $@ && rm -f tmp.BLAS_dwaxpby_s_d_x-f2c.c*

$(SRC_PATH)/BLAS_dwaxpby_s_s_x-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="d, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dwaxpby_s_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dwaxpby_s_s_x-f2c.c && \
    mv tmp.BLAS_dwaxpby_s_s_x-f2c.c $@ && rm -f tmp.BLAS_dwaxpby_s_s_x-f2c.c*

$(SRC_PATH)/BLAS_zwaxpby_z_c_x-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="z, z, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zwaxpby_z_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby_z_c_x-f2c.c && \
    mv tmp.BLAS_zwaxpby_z_c_x-f2c.c $@ && rm -f tmp.BLAS_zwaxpby_z_c_x-f2c.c*

$(SRC_PATH)/BLAS_zwaxpby_c_z_x-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="z, c, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zwaxpby_c_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby_c_z_x-f2c.c && \
    mv tmp.BLAS_zwaxpby_c_z_x-f2c.c $@ && rm -f tmp.BLAS_zwaxpby_c_z_x-f2c.c*

$(SRC_PATH)/BLAS_zwaxpby_c_c_x-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="z, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zwaxpby_c_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby_c_c_x-f2c.c && \
    mv tmp.BLAS_zwaxpby_c_c_x-f2c.c $@ && rm -f tmp.BLAS_zwaxpby_c_c_x-f2c.c*

$(SRC_PATH)/BLAS_cwaxpby_c_s_x-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="c, c, s, _x" ../f2c-bridge.m4 > tmp.BLAS_cwaxpby_c_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cwaxpby_c_s_x-f2c.c && \
    mv tmp.BLAS_cwaxpby_c_s_x-f2c.c $@ && rm -f tmp.BLAS_cwaxpby_c_s_x-f2c.c*

$(SRC_PATH)/BLAS_cwaxpby_s_c_x-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="c, s, c, _x" ../f2c-bridge.m4 > tmp.BLAS_cwaxpby_s_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cwaxpby_s_c_x-f2c.c && \
    mv tmp.BLAS_cwaxpby_s_c_x-f2c.c $@ && rm -f tmp.BLAS_cwaxpby_s_c_x-f2c.c*

$(SRC_PATH)/BLAS_cwaxpby_s_s_x-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="c, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_cwaxpby_s_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cwaxpby_s_s_x-f2c.c && \
    mv tmp.BLAS_cwaxpby_s_s_x-f2c.c $@ && rm -f tmp.BLAS_cwaxpby_s_s_x-f2c.c*

$(SRC_PATH)/BLAS_zwaxpby_z_d_x-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="z, z, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zwaxpby_z_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby_z_d_x-f2c.c && \
    mv tmp.BLAS_zwaxpby_z_d_x-f2c.c $@ && rm -f tmp.BLAS_zwaxpby_z_d_x-f2c.c*

$(SRC_PATH)/BLAS_zwaxpby_d_z_x-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="z, d, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zwaxpby_d_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby_d_z_x-f2c.c && \
    mv tmp.BLAS_zwaxpby_d_z_x-f2c.c $@ && rm -f tmp.BLAS_zwaxpby_d_z_x-f2c.c*

$(SRC_PATH)/BLAS_zwaxpby_d_d_x-f2c.c: ../f2c-bridge.m4 waxpby-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=waxpby -DRARG="z, d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zwaxpby_d_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zwaxpby_d_d_x-f2c.c && \
    mv tmp.BLAS_zwaxpby_d_d_x-f2c.c $@ && rm -f tmp.BLAS_zwaxpby_d_d_x-f2c.c*


waxpby.h: ../cblas.m4 ../header.template waxpby-common.m4
	sed -e "s/routine/waxpby/g" -e "s/ROUTINE/WAXPBY/g" \
  ../header.template | $(M4) $(M4_OPTS) $(M4_PATH) - >$@

source-clean:
	rm -f $(WAXPBY_SRCS) waxpby.h

.PHONY: waxpby-source waxpby-header source-clean

