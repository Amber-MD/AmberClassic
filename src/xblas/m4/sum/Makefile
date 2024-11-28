include ../../make.conf
include ../../$(MAKEINC)


SRC_PATH = ../../src/sum
M4_PATH = -I. -I..

SUM_SRCS = \
	$(SRC_PATH)/BLAS_ssum.c\
	$(SRC_PATH)/BLAS_dsum.c\
	$(SRC_PATH)/BLAS_csum.c\
	$(SRC_PATH)/BLAS_zsum.c\
	$(SRC_PATH)/BLAS_ssum_x.c\
	$(SRC_PATH)/BLAS_dsum_x.c\
	$(SRC_PATH)/BLAS_csum_x.c\
	$(SRC_PATH)/BLAS_zsum_x.c \
	\
	$(SRC_PATH)/BLAS_ssum-f2c.c\
	$(SRC_PATH)/BLAS_dsum-f2c.c\
	$(SRC_PATH)/BLAS_csum-f2c.c\
	$(SRC_PATH)/BLAS_zsum-f2c.c\
	$(SRC_PATH)/BLAS_ssum_x-f2c.c\
	$(SRC_PATH)/BLAS_dsum_x-f2c.c\
	$(SRC_PATH)/BLAS_csum_x-f2c.c\
	$(SRC_PATH)/BLAS_zsum_x-f2c.c

sum-source: $(SUM_SRCS)

sum-header: sum.h

$(SRC_PATH)/BLAS_ssum.c: sum.m4 sum-common.m4 ../cblas.m4
	echo "SUM(s)" | cat sum.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ssum.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ssum.c && \
    mv tmp.BLAS_ssum.c $@ && rm -f tmp.BLAS_ssum.c*

$(SRC_PATH)/BLAS_dsum.c: sum.m4 sum-common.m4 ../cblas.m4
	echo "SUM(d)" | cat sum.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsum.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsum.c && \
    mv tmp.BLAS_dsum.c $@ && rm -f tmp.BLAS_dsum.c*

$(SRC_PATH)/BLAS_csum.c: sum.m4 sum-common.m4 ../cblas.m4
	echo "SUM(c)" | cat sum.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csum.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csum.c && \
    mv tmp.BLAS_csum.c $@ && rm -f tmp.BLAS_csum.c*

$(SRC_PATH)/BLAS_zsum.c: sum.m4 sum-common.m4 ../cblas.m4
	echo "SUM(z)" | cat sum.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsum.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsum.c && \
    mv tmp.BLAS_zsum.c $@ && rm -f tmp.BLAS_zsum.c*

$(SRC_PATH)/BLAS_ssum_x.c: sum.m4 sum-common.m4 ../cblas.m4
	echo "SUM(s, _x)" | cat sum.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ssum_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ssum_x.c && \
    mv tmp.BLAS_ssum_x.c $@ && rm -f tmp.BLAS_ssum_x.c*

$(SRC_PATH)/BLAS_dsum_x.c: sum.m4 sum-common.m4 ../cblas.m4
	echo "SUM(d, _x)" | cat sum.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dsum_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsum_x.c && \
    mv tmp.BLAS_dsum_x.c $@ && rm -f tmp.BLAS_dsum_x.c*

$(SRC_PATH)/BLAS_csum_x.c: sum.m4 sum-common.m4 ../cblas.m4
	echo "SUM(c, _x)" | cat sum.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_csum_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csum_x.c && \
    mv tmp.BLAS_csum_x.c $@ && rm -f tmp.BLAS_csum_x.c*

$(SRC_PATH)/BLAS_zsum_x.c: sum.m4 sum-common.m4 ../cblas.m4
	echo "SUM(z, _x)" | cat sum.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_zsum_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsum_x.c && \
    mv tmp.BLAS_zsum_x.c $@ && rm -f tmp.BLAS_zsum_x.c*



$(SRC_PATH)/BLAS_ssum-f2c.c: ../f2c-bridge.m4 sum-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sum -DRARG="s" ../f2c-bridge.m4 > tmp.BLAS_ssum-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ssum-f2c.c && \
    mv tmp.BLAS_ssum-f2c.c $@ && rm -f tmp.BLAS_ssum-f2c.c*

$(SRC_PATH)/BLAS_dsum-f2c.c: ../f2c-bridge.m4 sum-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sum -DRARG="d" ../f2c-bridge.m4 > tmp.BLAS_dsum-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsum-f2c.c && \
    mv tmp.BLAS_dsum-f2c.c $@ && rm -f tmp.BLAS_dsum-f2c.c*

$(SRC_PATH)/BLAS_csum-f2c.c: ../f2c-bridge.m4 sum-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sum -DRARG="c" ../f2c-bridge.m4 > tmp.BLAS_csum-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csum-f2c.c && \
    mv tmp.BLAS_csum-f2c.c $@ && rm -f tmp.BLAS_csum-f2c.c*

$(SRC_PATH)/BLAS_zsum-f2c.c: ../f2c-bridge.m4 sum-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sum -DRARG="z" ../f2c-bridge.m4 > tmp.BLAS_zsum-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsum-f2c.c && \
    mv tmp.BLAS_zsum-f2c.c $@ && rm -f tmp.BLAS_zsum-f2c.c*

$(SRC_PATH)/BLAS_ssum_x-f2c.c: ../f2c-bridge.m4 sum-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sum -DRARG="s, _x" ../f2c-bridge.m4 > tmp.BLAS_ssum_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ssum_x-f2c.c && \
    mv tmp.BLAS_ssum_x-f2c.c $@ && rm -f tmp.BLAS_ssum_x-f2c.c*

$(SRC_PATH)/BLAS_dsum_x-f2c.c: ../f2c-bridge.m4 sum-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sum -DRARG="d, _x" ../f2c-bridge.m4 > tmp.BLAS_dsum_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dsum_x-f2c.c && \
    mv tmp.BLAS_dsum_x-f2c.c $@ && rm -f tmp.BLAS_dsum_x-f2c.c*

$(SRC_PATH)/BLAS_csum_x-f2c.c: ../f2c-bridge.m4 sum-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sum -DRARG="c, _x" ../f2c-bridge.m4 > tmp.BLAS_csum_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_csum_x-f2c.c && \
    mv tmp.BLAS_csum_x-f2c.c $@ && rm -f tmp.BLAS_csum_x-f2c.c*

$(SRC_PATH)/BLAS_zsum_x-f2c.c: ../f2c-bridge.m4 sum-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=sum -DRARG="z, _x" ../f2c-bridge.m4 > tmp.BLAS_zsum_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_zsum_x-f2c.c && \
    mv tmp.BLAS_zsum_x-f2c.c $@ && rm -f tmp.BLAS_zsum_x-f2c.c*


sum.h: ../cblas.m4 ../header.template sum-common.m4
	sed -e "s/routine/sum/g" -e "s/ROUTINE/SUM/g" \
  ../header.template | $(M4) $(M4_OPTS) $(M4_PATH) - >$@

source-clean:
	rm -f $(SUM_SRCS) sum.h

.PHONY: sum-source sum-header source-clean

