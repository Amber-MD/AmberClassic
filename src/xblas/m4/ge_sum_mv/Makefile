include ../../make.conf
include ../../$(MAKEINC)


SRC_PATH = ../../src/ge_sum_mv
M4_PATH = -I. -I..

GE_SUM_MV_SRCS = \
	$(SRC_PATH)/BLAS_sge_sum_mv.c\
	$(SRC_PATH)/BLAS_dge_sum_mv.c\
	$(SRC_PATH)/BLAS_cge_sum_mv.c\
	$(SRC_PATH)/BLAS_zge_sum_mv.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_d_s.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_s_d.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_s_s.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_z_c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_c_z.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_c_c.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_c_s.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_s_c.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_s_s.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_z_d.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_d_z.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_d_d.c\
	$(SRC_PATH)/BLAS_sge_sum_mv_x.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_x.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_x.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_x.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_d_s_x.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_s_d_x.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_s_s_x.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_z_c_x.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_c_z_x.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_c_c_x.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_c_s_x.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_s_c_x.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_s_s_x.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_z_d_x.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_d_z_x.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_d_d_x.c \
	\
	$(SRC_PATH)/BLAS_sge_sum_mv-f2c.c\
	$(SRC_PATH)/BLAS_dge_sum_mv-f2c.c\
	$(SRC_PATH)/BLAS_cge_sum_mv-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv-f2c.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_d_s-f2c.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_s_d-f2c.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_c_c-f2c.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_c_s-f2c.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_s_c-f2c.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_z_d-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_d_z-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_d_d-f2c.c\
	$(SRC_PATH)/BLAS_sge_sum_mv_x-f2c.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_x-f2c.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_x-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_x-f2c.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_d_s_x-f2c.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_s_d_x-f2c.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_s_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_z_d_x-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_d_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_d_d_x-f2c.c

ge_sum_mv-source: $(GE_SUM_MV_SRCS)

ge_sum_mv-header: ge_sum_mv.h

$(SRC_PATH)/BLAS_sge_sum_mv.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(s, s, s)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_sge_sum_mv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sge_sum_mv.c && \
    mv tmp.BLAS_sge_sum_mv.c $@ && rm -f tmp.BLAS_sge_sum_mv.c*

$(SRC_PATH)/BLAS_dge_sum_mv.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(d, d, d)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dge_sum_mv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dge_sum_mv.c && \
    mv tmp.BLAS_dge_sum_mv.c $@ && rm -f tmp.BLAS_dge_sum_mv.c*

$(SRC_PATH)/BLAS_cge_sum_mv.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(c, c, c)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cge_sum_mv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cge_sum_mv.c && \
    mv tmp.BLAS_cge_sum_mv.c $@ && rm -f tmp.BLAS_cge_sum_mv.c*

$(SRC_PATH)/BLAS_zge_sum_mv.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(z, z, z)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zge_sum_mv.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv.c && \
    mv tmp.BLAS_zge_sum_mv.c $@ && rm -f tmp.BLAS_zge_sum_mv.c*

$(SRC_PATH)/BLAS_dge_sum_mv_d_s.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(d, d, s)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dge_sum_mv_d_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dge_sum_mv_d_s.c && \
    mv tmp.BLAS_dge_sum_mv_d_s.c $@ && rm -f tmp.BLAS_dge_sum_mv_d_s.c*

$(SRC_PATH)/BLAS_dge_sum_mv_s_d.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(d, s, d)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dge_sum_mv_s_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dge_sum_mv_s_d.c && \
    mv tmp.BLAS_dge_sum_mv_s_d.c $@ && rm -f tmp.BLAS_dge_sum_mv_s_d.c*

$(SRC_PATH)/BLAS_dge_sum_mv_s_s.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(d, s, s)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dge_sum_mv_s_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dge_sum_mv_s_s.c && \
    mv tmp.BLAS_dge_sum_mv_s_s.c $@ && rm -f tmp.BLAS_dge_sum_mv_s_s.c*

$(SRC_PATH)/BLAS_zge_sum_mv_z_c.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(z, z, c)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zge_sum_mv_z_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv_z_c.c && \
    mv tmp.BLAS_zge_sum_mv_z_c.c $@ && rm -f tmp.BLAS_zge_sum_mv_z_c.c*

$(SRC_PATH)/BLAS_zge_sum_mv_c_z.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(z, c, z)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zge_sum_mv_c_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv_c_z.c && \
    mv tmp.BLAS_zge_sum_mv_c_z.c $@ && rm -f tmp.BLAS_zge_sum_mv_c_z.c*

$(SRC_PATH)/BLAS_zge_sum_mv_c_c.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(z, c, c)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zge_sum_mv_c_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv_c_c.c && \
    mv tmp.BLAS_zge_sum_mv_c_c.c $@ && rm -f tmp.BLAS_zge_sum_mv_c_c.c*

$(SRC_PATH)/BLAS_cge_sum_mv_c_s.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(c, c, s)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cge_sum_mv_c_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cge_sum_mv_c_s.c && \
    mv tmp.BLAS_cge_sum_mv_c_s.c $@ && rm -f tmp.BLAS_cge_sum_mv_c_s.c*

$(SRC_PATH)/BLAS_cge_sum_mv_s_c.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(c, s, c)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cge_sum_mv_s_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cge_sum_mv_s_c.c && \
    mv tmp.BLAS_cge_sum_mv_s_c.c $@ && rm -f tmp.BLAS_cge_sum_mv_s_c.c*

$(SRC_PATH)/BLAS_cge_sum_mv_s_s.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(c, s, s)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cge_sum_mv_s_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cge_sum_mv_s_s.c && \
    mv tmp.BLAS_cge_sum_mv_s_s.c $@ && rm -f tmp.BLAS_cge_sum_mv_s_s.c*

$(SRC_PATH)/BLAS_zge_sum_mv_z_d.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(z, z, d)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zge_sum_mv_z_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv_z_d.c && \
    mv tmp.BLAS_zge_sum_mv_z_d.c $@ && rm -f tmp.BLAS_zge_sum_mv_z_d.c*

$(SRC_PATH)/BLAS_zge_sum_mv_d_z.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(z, d, z)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zge_sum_mv_d_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv_d_z.c && \
    mv tmp.BLAS_zge_sum_mv_d_z.c $@ && rm -f tmp.BLAS_zge_sum_mv_d_z.c*

$(SRC_PATH)/BLAS_zge_sum_mv_d_d.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(z, d, d)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zge_sum_mv_d_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv_d_d.c && \
    mv tmp.BLAS_zge_sum_mv_d_d.c $@ && rm -f tmp.BLAS_zge_sum_mv_d_d.c*

$(SRC_PATH)/BLAS_sge_sum_mv_x.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(s, s, s, _x)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_sge_sum_mv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sge_sum_mv_x.c && \
    mv tmp.BLAS_sge_sum_mv_x.c $@ && rm -f tmp.BLAS_sge_sum_mv_x.c*

$(SRC_PATH)/BLAS_dge_sum_mv_x.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(d, d, d, _x)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dge_sum_mv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dge_sum_mv_x.c && \
    mv tmp.BLAS_dge_sum_mv_x.c $@ && rm -f tmp.BLAS_dge_sum_mv_x.c*

$(SRC_PATH)/BLAS_cge_sum_mv_x.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(c, c, c, _x)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cge_sum_mv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cge_sum_mv_x.c && \
    mv tmp.BLAS_cge_sum_mv_x.c $@ && rm -f tmp.BLAS_cge_sum_mv_x.c*

$(SRC_PATH)/BLAS_zge_sum_mv_x.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(z, z, z, _x)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zge_sum_mv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv_x.c && \
    mv tmp.BLAS_zge_sum_mv_x.c $@ && rm -f tmp.BLAS_zge_sum_mv_x.c*

$(SRC_PATH)/BLAS_dge_sum_mv_d_s_x.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(d, d, s, _x)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dge_sum_mv_d_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dge_sum_mv_d_s_x.c && \
    mv tmp.BLAS_dge_sum_mv_d_s_x.c $@ && rm -f tmp.BLAS_dge_sum_mv_d_s_x.c*

$(SRC_PATH)/BLAS_dge_sum_mv_s_d_x.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(d, s, d, _x)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dge_sum_mv_s_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dge_sum_mv_s_d_x.c && \
    mv tmp.BLAS_dge_sum_mv_s_d_x.c $@ && rm -f tmp.BLAS_dge_sum_mv_s_d_x.c*

$(SRC_PATH)/BLAS_dge_sum_mv_s_s_x.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(d, s, s, _x)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dge_sum_mv_s_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dge_sum_mv_s_s_x.c && \
    mv tmp.BLAS_dge_sum_mv_s_s_x.c $@ && rm -f tmp.BLAS_dge_sum_mv_s_s_x.c*

$(SRC_PATH)/BLAS_zge_sum_mv_z_c_x.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(z, z, c, _x)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zge_sum_mv_z_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv_z_c_x.c && \
    mv tmp.BLAS_zge_sum_mv_z_c_x.c $@ && rm -f tmp.BLAS_zge_sum_mv_z_c_x.c*

$(SRC_PATH)/BLAS_zge_sum_mv_c_z_x.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(z, c, z, _x)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zge_sum_mv_c_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv_c_z_x.c && \
    mv tmp.BLAS_zge_sum_mv_c_z_x.c $@ && rm -f tmp.BLAS_zge_sum_mv_c_z_x.c*

$(SRC_PATH)/BLAS_zge_sum_mv_c_c_x.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(z, c, c, _x)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zge_sum_mv_c_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv_c_c_x.c && \
    mv tmp.BLAS_zge_sum_mv_c_c_x.c $@ && rm -f tmp.BLAS_zge_sum_mv_c_c_x.c*

$(SRC_PATH)/BLAS_cge_sum_mv_c_s_x.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(c, c, s, _x)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cge_sum_mv_c_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cge_sum_mv_c_s_x.c && \
    mv tmp.BLAS_cge_sum_mv_c_s_x.c $@ && rm -f tmp.BLAS_cge_sum_mv_c_s_x.c*

$(SRC_PATH)/BLAS_cge_sum_mv_s_c_x.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(c, s, c, _x)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cge_sum_mv_s_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cge_sum_mv_s_c_x.c && \
    mv tmp.BLAS_cge_sum_mv_s_c_x.c $@ && rm -f tmp.BLAS_cge_sum_mv_s_c_x.c*

$(SRC_PATH)/BLAS_cge_sum_mv_s_s_x.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(c, s, s, _x)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_cge_sum_mv_s_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cge_sum_mv_s_s_x.c && \
    mv tmp.BLAS_cge_sum_mv_s_s_x.c $@ && rm -f tmp.BLAS_cge_sum_mv_s_s_x.c*

$(SRC_PATH)/BLAS_zge_sum_mv_z_d_x.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(z, z, d, _x)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zge_sum_mv_z_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv_z_d_x.c && \
    mv tmp.BLAS_zge_sum_mv_z_d_x.c $@ && rm -f tmp.BLAS_zge_sum_mv_z_d_x.c*

$(SRC_PATH)/BLAS_zge_sum_mv_d_z_x.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(z, d, z, _x)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zge_sum_mv_d_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv_d_z_x.c && \
    mv tmp.BLAS_zge_sum_mv_d_z_x.c $@ && rm -f tmp.BLAS_zge_sum_mv_d_z_x.c*

$(SRC_PATH)/BLAS_zge_sum_mv_d_d_x.c: ge_sum_mv.m4 ge_sum_mv-common.m4 ../cblas.m4
	echo "GE_SUM_MV(z, d, d, _x)" | cat ge_sum_mv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zge_sum_mv_d_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv_d_d_x.c && \
    mv tmp.BLAS_zge_sum_mv_d_d_x.c $@ && rm -f tmp.BLAS_zge_sum_mv_d_d_x.c*



$(SRC_PATH)/BLAS_sge_sum_mv-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="s, s, s" ../f2c-bridge.m4 > tmp.BLAS_sge_sum_mv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sge_sum_mv-f2c.c && \
    mv tmp.BLAS_sge_sum_mv-f2c.c $@ && rm -f tmp.BLAS_sge_sum_mv-f2c.c*

$(SRC_PATH)/BLAS_dge_sum_mv-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="d, d, d" ../f2c-bridge.m4 > tmp.BLAS_dge_sum_mv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dge_sum_mv-f2c.c && \
    mv tmp.BLAS_dge_sum_mv-f2c.c $@ && rm -f tmp.BLAS_dge_sum_mv-f2c.c*

$(SRC_PATH)/BLAS_cge_sum_mv-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="c, c, c" ../f2c-bridge.m4 > tmp.BLAS_cge_sum_mv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cge_sum_mv-f2c.c && \
    mv tmp.BLAS_cge_sum_mv-f2c.c $@ && rm -f tmp.BLAS_cge_sum_mv-f2c.c*

$(SRC_PATH)/BLAS_zge_sum_mv-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="z, z, z" ../f2c-bridge.m4 > tmp.BLAS_zge_sum_mv-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv-f2c.c && \
    mv tmp.BLAS_zge_sum_mv-f2c.c $@ && rm -f tmp.BLAS_zge_sum_mv-f2c.c*

$(SRC_PATH)/BLAS_dge_sum_mv_d_s-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="d, d, s" ../f2c-bridge.m4 > tmp.BLAS_dge_sum_mv_d_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dge_sum_mv_d_s-f2c.c && \
    mv tmp.BLAS_dge_sum_mv_d_s-f2c.c $@ && rm -f tmp.BLAS_dge_sum_mv_d_s-f2c.c*

$(SRC_PATH)/BLAS_dge_sum_mv_s_d-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="d, s, d" ../f2c-bridge.m4 > tmp.BLAS_dge_sum_mv_s_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dge_sum_mv_s_d-f2c.c && \
    mv tmp.BLAS_dge_sum_mv_s_d-f2c.c $@ && rm -f tmp.BLAS_dge_sum_mv_s_d-f2c.c*

$(SRC_PATH)/BLAS_dge_sum_mv_s_s-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="d, s, s" ../f2c-bridge.m4 > tmp.BLAS_dge_sum_mv_s_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dge_sum_mv_s_s-f2c.c && \
    mv tmp.BLAS_dge_sum_mv_s_s-f2c.c $@ && rm -f tmp.BLAS_dge_sum_mv_s_s-f2c.c*

$(SRC_PATH)/BLAS_zge_sum_mv_z_c-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="z, z, c" ../f2c-bridge.m4 > tmp.BLAS_zge_sum_mv_z_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv_z_c-f2c.c && \
    mv tmp.BLAS_zge_sum_mv_z_c-f2c.c $@ && rm -f tmp.BLAS_zge_sum_mv_z_c-f2c.c*

$(SRC_PATH)/BLAS_zge_sum_mv_c_z-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="z, c, z" ../f2c-bridge.m4 > tmp.BLAS_zge_sum_mv_c_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv_c_z-f2c.c && \
    mv tmp.BLAS_zge_sum_mv_c_z-f2c.c $@ && rm -f tmp.BLAS_zge_sum_mv_c_z-f2c.c*

$(SRC_PATH)/BLAS_zge_sum_mv_c_c-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="z, c, c" ../f2c-bridge.m4 > tmp.BLAS_zge_sum_mv_c_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv_c_c-f2c.c && \
    mv tmp.BLAS_zge_sum_mv_c_c-f2c.c $@ && rm -f tmp.BLAS_zge_sum_mv_c_c-f2c.c*

$(SRC_PATH)/BLAS_cge_sum_mv_c_s-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="c, c, s" ../f2c-bridge.m4 > tmp.BLAS_cge_sum_mv_c_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cge_sum_mv_c_s-f2c.c && \
    mv tmp.BLAS_cge_sum_mv_c_s-f2c.c $@ && rm -f tmp.BLAS_cge_sum_mv_c_s-f2c.c*

$(SRC_PATH)/BLAS_cge_sum_mv_s_c-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="c, s, c" ../f2c-bridge.m4 > tmp.BLAS_cge_sum_mv_s_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cge_sum_mv_s_c-f2c.c && \
    mv tmp.BLAS_cge_sum_mv_s_c-f2c.c $@ && rm -f tmp.BLAS_cge_sum_mv_s_c-f2c.c*

$(SRC_PATH)/BLAS_cge_sum_mv_s_s-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="c, s, s" ../f2c-bridge.m4 > tmp.BLAS_cge_sum_mv_s_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cge_sum_mv_s_s-f2c.c && \
    mv tmp.BLAS_cge_sum_mv_s_s-f2c.c $@ && rm -f tmp.BLAS_cge_sum_mv_s_s-f2c.c*

$(SRC_PATH)/BLAS_zge_sum_mv_z_d-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="z, z, d" ../f2c-bridge.m4 > tmp.BLAS_zge_sum_mv_z_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv_z_d-f2c.c && \
    mv tmp.BLAS_zge_sum_mv_z_d-f2c.c $@ && rm -f tmp.BLAS_zge_sum_mv_z_d-f2c.c*

$(SRC_PATH)/BLAS_zge_sum_mv_d_z-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="z, d, z" ../f2c-bridge.m4 > tmp.BLAS_zge_sum_mv_d_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv_d_z-f2c.c && \
    mv tmp.BLAS_zge_sum_mv_d_z-f2c.c $@ && rm -f tmp.BLAS_zge_sum_mv_d_z-f2c.c*

$(SRC_PATH)/BLAS_zge_sum_mv_d_d-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="z, d, d" ../f2c-bridge.m4 > tmp.BLAS_zge_sum_mv_d_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv_d_d-f2c.c && \
    mv tmp.BLAS_zge_sum_mv_d_d-f2c.c $@ && rm -f tmp.BLAS_zge_sum_mv_d_d-f2c.c*

$(SRC_PATH)/BLAS_sge_sum_mv_x-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="s, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_sge_sum_mv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_sge_sum_mv_x-f2c.c && \
    mv tmp.BLAS_sge_sum_mv_x-f2c.c $@ && rm -f tmp.BLAS_sge_sum_mv_x-f2c.c*

$(SRC_PATH)/BLAS_dge_sum_mv_x-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="d, d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dge_sum_mv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dge_sum_mv_x-f2c.c && \
    mv tmp.BLAS_dge_sum_mv_x-f2c.c $@ && rm -f tmp.BLAS_dge_sum_mv_x-f2c.c*

$(SRC_PATH)/BLAS_cge_sum_mv_x-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="c, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_cge_sum_mv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cge_sum_mv_x-f2c.c && \
    mv tmp.BLAS_cge_sum_mv_x-f2c.c $@ && rm -f tmp.BLAS_cge_sum_mv_x-f2c.c*

$(SRC_PATH)/BLAS_zge_sum_mv_x-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="z, z, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zge_sum_mv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv_x-f2c.c && \
    mv tmp.BLAS_zge_sum_mv_x-f2c.c $@ && rm -f tmp.BLAS_zge_sum_mv_x-f2c.c*

$(SRC_PATH)/BLAS_dge_sum_mv_d_s_x-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="d, d, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dge_sum_mv_d_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dge_sum_mv_d_s_x-f2c.c && \
    mv tmp.BLAS_dge_sum_mv_d_s_x-f2c.c $@ && rm -f tmp.BLAS_dge_sum_mv_d_s_x-f2c.c*

$(SRC_PATH)/BLAS_dge_sum_mv_s_d_x-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="d, s, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dge_sum_mv_s_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dge_sum_mv_s_d_x-f2c.c && \
    mv tmp.BLAS_dge_sum_mv_s_d_x-f2c.c $@ && rm -f tmp.BLAS_dge_sum_mv_s_d_x-f2c.c*

$(SRC_PATH)/BLAS_dge_sum_mv_s_s_x-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="d, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dge_sum_mv_s_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dge_sum_mv_s_s_x-f2c.c && \
    mv tmp.BLAS_dge_sum_mv_s_s_x-f2c.c $@ && rm -f tmp.BLAS_dge_sum_mv_s_s_x-f2c.c*

$(SRC_PATH)/BLAS_zge_sum_mv_z_c_x-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="z, z, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zge_sum_mv_z_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv_z_c_x-f2c.c && \
    mv tmp.BLAS_zge_sum_mv_z_c_x-f2c.c $@ && rm -f tmp.BLAS_zge_sum_mv_z_c_x-f2c.c*

$(SRC_PATH)/BLAS_zge_sum_mv_c_z_x-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="z, c, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zge_sum_mv_c_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv_c_z_x-f2c.c && \
    mv tmp.BLAS_zge_sum_mv_c_z_x-f2c.c $@ && rm -f tmp.BLAS_zge_sum_mv_c_z_x-f2c.c*

$(SRC_PATH)/BLAS_zge_sum_mv_c_c_x-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="z, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zge_sum_mv_c_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv_c_c_x-f2c.c && \
    mv tmp.BLAS_zge_sum_mv_c_c_x-f2c.c $@ && rm -f tmp.BLAS_zge_sum_mv_c_c_x-f2c.c*

$(SRC_PATH)/BLAS_cge_sum_mv_c_s_x-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="c, c, s, _x" ../f2c-bridge.m4 > tmp.BLAS_cge_sum_mv_c_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cge_sum_mv_c_s_x-f2c.c && \
    mv tmp.BLAS_cge_sum_mv_c_s_x-f2c.c $@ && rm -f tmp.BLAS_cge_sum_mv_c_s_x-f2c.c*

$(SRC_PATH)/BLAS_cge_sum_mv_s_c_x-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="c, s, c, _x" ../f2c-bridge.m4 > tmp.BLAS_cge_sum_mv_s_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cge_sum_mv_s_c_x-f2c.c && \
    mv tmp.BLAS_cge_sum_mv_s_c_x-f2c.c $@ && rm -f tmp.BLAS_cge_sum_mv_s_c_x-f2c.c*

$(SRC_PATH)/BLAS_cge_sum_mv_s_s_x-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="c, s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_cge_sum_mv_s_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_cge_sum_mv_s_s_x-f2c.c && \
    mv tmp.BLAS_cge_sum_mv_s_s_x-f2c.c $@ && rm -f tmp.BLAS_cge_sum_mv_s_s_x-f2c.c*

$(SRC_PATH)/BLAS_zge_sum_mv_z_d_x-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="z, z, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zge_sum_mv_z_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv_z_d_x-f2c.c && \
    mv tmp.BLAS_zge_sum_mv_z_d_x-f2c.c $@ && rm -f tmp.BLAS_zge_sum_mv_z_d_x-f2c.c*

$(SRC_PATH)/BLAS_zge_sum_mv_d_z_x-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="z, d, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zge_sum_mv_d_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv_d_z_x-f2c.c && \
    mv tmp.BLAS_zge_sum_mv_d_z_x-f2c.c $@ && rm -f tmp.BLAS_zge_sum_mv_d_z_x-f2c.c*

$(SRC_PATH)/BLAS_zge_sum_mv_d_d_x-f2c.c: ../f2c-bridge.m4 ge_sum_mv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=ge_sum_mv -DRARG="z, d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zge_sum_mv_d_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zge_sum_mv_d_d_x-f2c.c && \
    mv tmp.BLAS_zge_sum_mv_d_d_x-f2c.c $@ && rm -f tmp.BLAS_zge_sum_mv_d_d_x-f2c.c*


ge_sum_mv.h: ../cblas.m4 ../header.template ge_sum_mv-common.m4
	sed -e "s/routine/ge_sum_mv/g" -e "s/ROUTINE/GE_SUM_MV/g" \
  ../header.template | $(M4) $(M4_OPTS) $(M4_PATH) - >$@

source-clean:
	rm -f $(GE_SUM_MV_SRCS) ge_sum_mv.h

.PHONY: ge_sum_mv-source ge_sum_mv-header source-clean

