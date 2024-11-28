include ../../make.conf
include ../../$(MAKEINC)


SRC_PATH = ../../src/hemm
M4_PATH = -I. -I..

HEMM_SRCS = \
	$(SRC_PATH)/BLAS_chemm.c\
	$(SRC_PATH)/BLAS_zhemm.c\
	$(SRC_PATH)/BLAS_zhemm_z_c.c\
	$(SRC_PATH)/BLAS_zhemm_c_z.c\
	$(SRC_PATH)/BLAS_zhemm_c_c.c\
	$(SRC_PATH)/BLAS_chemm_c_s.c\
	$(SRC_PATH)/BLAS_zhemm_z_d.c\
	$(SRC_PATH)/BLAS_chemm_x.c\
	$(SRC_PATH)/BLAS_zhemm_x.c\
	$(SRC_PATH)/BLAS_zhemm_z_c_x.c\
	$(SRC_PATH)/BLAS_zhemm_c_z_x.c\
	$(SRC_PATH)/BLAS_zhemm_c_c_x.c\
	$(SRC_PATH)/BLAS_chemm_c_s_x.c\
	$(SRC_PATH)/BLAS_zhemm_z_d_x.c \
	\
	$(SRC_PATH)/BLAS_chemm-f2c.c\
	$(SRC_PATH)/BLAS_zhemm-f2c.c\
	$(SRC_PATH)/BLAS_zhemm_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zhemm_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zhemm_c_c-f2c.c\
	$(SRC_PATH)/BLAS_chemm_c_s-f2c.c\
	$(SRC_PATH)/BLAS_zhemm_z_d-f2c.c\
	$(SRC_PATH)/BLAS_chemm_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemm_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemm_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemm_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemm_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_chemm_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemm_z_d_x-f2c.c

hemm-source: $(HEMM_SRCS)

hemm-header: hemm.h

$(SRC_PATH)/BLAS_chemm.c: hemm.m4 hemm-common.m4 ../cblas.m4
	echo "HEMM(c, c, c)" | cat hemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_chemm.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chemm.c && \
    mv tmp.BLAS_chemm.c $@ && rm -f tmp.BLAS_chemm.c*

$(SRC_PATH)/BLAS_zhemm.c: hemm.m4 hemm-common.m4 ../cblas.m4
	echo "HEMM(z, z, z)" | cat hemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemm.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemm.c && \
    mv tmp.BLAS_zhemm.c $@ && rm -f tmp.BLAS_zhemm.c*

$(SRC_PATH)/BLAS_zhemm_z_c.c: hemm.m4 hemm-common.m4 ../cblas.m4
	echo "HEMM(z, z, c)" | cat hemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemm_z_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemm_z_c.c && \
    mv tmp.BLAS_zhemm_z_c.c $@ && rm -f tmp.BLAS_zhemm_z_c.c*

$(SRC_PATH)/BLAS_zhemm_c_z.c: hemm.m4 hemm-common.m4 ../cblas.m4
	echo "HEMM(z, c, z)" | cat hemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemm_c_z.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemm_c_z.c && \
    mv tmp.BLAS_zhemm_c_z.c $@ && rm -f tmp.BLAS_zhemm_c_z.c*

$(SRC_PATH)/BLAS_zhemm_c_c.c: hemm.m4 hemm-common.m4 ../cblas.m4
	echo "HEMM(z, c, c)" | cat hemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemm_c_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemm_c_c.c && \
    mv tmp.BLAS_zhemm_c_c.c $@ && rm -f tmp.BLAS_zhemm_c_c.c*

$(SRC_PATH)/BLAS_chemm_c_s.c: hemm.m4 hemm-common.m4 ../cblas.m4
	echo "HEMM(c, c, s)" | cat hemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_chemm_c_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chemm_c_s.c && \
    mv tmp.BLAS_chemm_c_s.c $@ && rm -f tmp.BLAS_chemm_c_s.c*

$(SRC_PATH)/BLAS_zhemm_z_d.c: hemm.m4 hemm-common.m4 ../cblas.m4
	echo "HEMM(z, z, d)" | cat hemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemm_z_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemm_z_d.c && \
    mv tmp.BLAS_zhemm_z_d.c $@ && rm -f tmp.BLAS_zhemm_z_d.c*

$(SRC_PATH)/BLAS_chemm_x.c: hemm.m4 hemm-common.m4 ../cblas.m4
	echo "HEMM(c, c, c, _x)" | cat hemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_chemm_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chemm_x.c && \
    mv tmp.BLAS_chemm_x.c $@ && rm -f tmp.BLAS_chemm_x.c*

$(SRC_PATH)/BLAS_zhemm_x.c: hemm.m4 hemm-common.m4 ../cblas.m4
	echo "HEMM(z, z, z, _x)" | cat hemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemm_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemm_x.c && \
    mv tmp.BLAS_zhemm_x.c $@ && rm -f tmp.BLAS_zhemm_x.c*

$(SRC_PATH)/BLAS_zhemm_z_c_x.c: hemm.m4 hemm-common.m4 ../cblas.m4
	echo "HEMM(z, z, c, _x)" | cat hemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemm_z_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemm_z_c_x.c && \
    mv tmp.BLAS_zhemm_z_c_x.c $@ && rm -f tmp.BLAS_zhemm_z_c_x.c*

$(SRC_PATH)/BLAS_zhemm_c_z_x.c: hemm.m4 hemm-common.m4 ../cblas.m4
	echo "HEMM(z, c, z, _x)" | cat hemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemm_c_z_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemm_c_z_x.c && \
    mv tmp.BLAS_zhemm_c_z_x.c $@ && rm -f tmp.BLAS_zhemm_c_z_x.c*

$(SRC_PATH)/BLAS_zhemm_c_c_x.c: hemm.m4 hemm-common.m4 ../cblas.m4
	echo "HEMM(z, c, c, _x)" | cat hemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemm_c_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemm_c_c_x.c && \
    mv tmp.BLAS_zhemm_c_c_x.c $@ && rm -f tmp.BLAS_zhemm_c_c_x.c*

$(SRC_PATH)/BLAS_chemm_c_s_x.c: hemm.m4 hemm-common.m4 ../cblas.m4
	echo "HEMM(c, c, s, _x)" | cat hemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_chemm_c_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chemm_c_s_x.c && \
    mv tmp.BLAS_chemm_c_s_x.c $@ && rm -f tmp.BLAS_chemm_c_s_x.c*

$(SRC_PATH)/BLAS_zhemm_z_d_x.c: hemm.m4 hemm-common.m4 ../cblas.m4
	echo "HEMM(z, z, d, _x)" | cat hemm.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zhemm_z_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemm_z_d_x.c && \
    mv tmp.BLAS_zhemm_z_d_x.c $@ && rm -f tmp.BLAS_zhemm_z_d_x.c*



$(SRC_PATH)/BLAS_chemm-f2c.c: ../f2c-bridge.m4 hemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemm -DRARG="c, c, c" ../f2c-bridge.m4 > tmp.BLAS_chemm-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chemm-f2c.c && \
    mv tmp.BLAS_chemm-f2c.c $@ && rm -f tmp.BLAS_chemm-f2c.c*

$(SRC_PATH)/BLAS_zhemm-f2c.c: ../f2c-bridge.m4 hemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemm -DRARG="z, z, z" ../f2c-bridge.m4 > tmp.BLAS_zhemm-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemm-f2c.c && \
    mv tmp.BLAS_zhemm-f2c.c $@ && rm -f tmp.BLAS_zhemm-f2c.c*

$(SRC_PATH)/BLAS_zhemm_z_c-f2c.c: ../f2c-bridge.m4 hemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemm -DRARG="z, z, c" ../f2c-bridge.m4 > tmp.BLAS_zhemm_z_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemm_z_c-f2c.c && \
    mv tmp.BLAS_zhemm_z_c-f2c.c $@ && rm -f tmp.BLAS_zhemm_z_c-f2c.c*

$(SRC_PATH)/BLAS_zhemm_c_z-f2c.c: ../f2c-bridge.m4 hemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemm -DRARG="z, c, z" ../f2c-bridge.m4 > tmp.BLAS_zhemm_c_z-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemm_c_z-f2c.c && \
    mv tmp.BLAS_zhemm_c_z-f2c.c $@ && rm -f tmp.BLAS_zhemm_c_z-f2c.c*

$(SRC_PATH)/BLAS_zhemm_c_c-f2c.c: ../f2c-bridge.m4 hemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemm -DRARG="z, c, c" ../f2c-bridge.m4 > tmp.BLAS_zhemm_c_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemm_c_c-f2c.c && \
    mv tmp.BLAS_zhemm_c_c-f2c.c $@ && rm -f tmp.BLAS_zhemm_c_c-f2c.c*

$(SRC_PATH)/BLAS_chemm_c_s-f2c.c: ../f2c-bridge.m4 hemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemm -DRARG="c, c, s" ../f2c-bridge.m4 > tmp.BLAS_chemm_c_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chemm_c_s-f2c.c && \
    mv tmp.BLAS_chemm_c_s-f2c.c $@ && rm -f tmp.BLAS_chemm_c_s-f2c.c*

$(SRC_PATH)/BLAS_zhemm_z_d-f2c.c: ../f2c-bridge.m4 hemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemm -DRARG="z, z, d" ../f2c-bridge.m4 > tmp.BLAS_zhemm_z_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemm_z_d-f2c.c && \
    mv tmp.BLAS_zhemm_z_d-f2c.c $@ && rm -f tmp.BLAS_zhemm_z_d-f2c.c*

$(SRC_PATH)/BLAS_chemm_x-f2c.c: ../f2c-bridge.m4 hemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemm -DRARG="c, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_chemm_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chemm_x-f2c.c && \
    mv tmp.BLAS_chemm_x-f2c.c $@ && rm -f tmp.BLAS_chemm_x-f2c.c*

$(SRC_PATH)/BLAS_zhemm_x-f2c.c: ../f2c-bridge.m4 hemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemm -DRARG="z, z, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zhemm_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemm_x-f2c.c && \
    mv tmp.BLAS_zhemm_x-f2c.c $@ && rm -f tmp.BLAS_zhemm_x-f2c.c*

$(SRC_PATH)/BLAS_zhemm_z_c_x-f2c.c: ../f2c-bridge.m4 hemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemm -DRARG="z, z, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zhemm_z_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemm_z_c_x-f2c.c && \
    mv tmp.BLAS_zhemm_z_c_x-f2c.c $@ && rm -f tmp.BLAS_zhemm_z_c_x-f2c.c*

$(SRC_PATH)/BLAS_zhemm_c_z_x-f2c.c: ../f2c-bridge.m4 hemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemm -DRARG="z, c, z, _x" ../f2c-bridge.m4 > tmp.BLAS_zhemm_c_z_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemm_c_z_x-f2c.c && \
    mv tmp.BLAS_zhemm_c_z_x-f2c.c $@ && rm -f tmp.BLAS_zhemm_c_z_x-f2c.c*

$(SRC_PATH)/BLAS_zhemm_c_c_x-f2c.c: ../f2c-bridge.m4 hemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemm -DRARG="z, c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_zhemm_c_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemm_c_c_x-f2c.c && \
    mv tmp.BLAS_zhemm_c_c_x-f2c.c $@ && rm -f tmp.BLAS_zhemm_c_c_x-f2c.c*

$(SRC_PATH)/BLAS_chemm_c_s_x-f2c.c: ../f2c-bridge.m4 hemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemm -DRARG="c, c, s, _x" ../f2c-bridge.m4 > tmp.BLAS_chemm_c_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_chemm_c_s_x-f2c.c && \
    mv tmp.BLAS_chemm_c_s_x-f2c.c $@ && rm -f tmp.BLAS_chemm_c_s_x-f2c.c*

$(SRC_PATH)/BLAS_zhemm_z_d_x-f2c.c: ../f2c-bridge.m4 hemm-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=hemm -DRARG="z, z, d, _x" ../f2c-bridge.m4 > tmp.BLAS_zhemm_z_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zhemm_z_d_x-f2c.c && \
    mv tmp.BLAS_zhemm_z_d_x-f2c.c $@ && rm -f tmp.BLAS_zhemm_z_d_x-f2c.c*


hemm.h: ../cblas.m4 ../header.template hemm-common.m4
	sed -e "s/routine/hemm/g" -e "s/ROUTINE/HEMM/g" \
  ../header.template | $(M4) $(M4_OPTS) $(M4_PATH) - >$@

source-clean:
	rm -f $(HEMM_SRCS) hemm.h

.PHONY: hemm-source hemm-header source-clean

