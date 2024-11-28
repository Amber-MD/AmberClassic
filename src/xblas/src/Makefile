include ../make.conf
include ../$(MAKEINC)

MKFILES = Makefile dot/Makefile sum/Makefile axpby/Makefile waxpby/Makefile gemv/Makefile ge_sum_mv/Makefile gbmv/Makefile symv/Makefile spmv/Makefile sbmv/Makefile hemv/Makefile hpmv/Makefile hbmv/Makefile trmv/Makefile tpmv/Makefile trsv/Makefile tbsv/Makefile gemm/Makefile symm/Makefile hemm/Makefile gemv2/Makefile symv2/Makefile hemv2/Makefile gbmv2/Makefile common/Makefile

makefiles: $(MKFILES)

Makefile: Makefile.m4 ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) Makefile.m4 >$@


dot:
	mkdir dot

dot/Makefile: dot Makefile.template ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=dot Makefile.template >$@

sum:
	mkdir sum

sum/Makefile: sum Makefile.template ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=sum Makefile.template >$@

axpby:
	mkdir axpby

axpby/Makefile: axpby Makefile.template ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=axpby Makefile.template >$@

waxpby:
	mkdir waxpby

waxpby/Makefile: waxpby Makefile.template ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=waxpby Makefile.template >$@

gemv:
	mkdir gemv

gemv/Makefile: gemv Makefile.template ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=gemv Makefile.template >$@

ge_sum_mv:
	mkdir ge_sum_mv

ge_sum_mv/Makefile: ge_sum_mv Makefile.template ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=ge_sum_mv Makefile.template >$@

gbmv:
	mkdir gbmv

gbmv/Makefile: gbmv Makefile.template ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=gbmv Makefile.template >$@

symv:
	mkdir symv

symv/Makefile: symv Makefile.template ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=symv Makefile.template >$@

spmv:
	mkdir spmv

spmv/Makefile: spmv Makefile.template ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=spmv Makefile.template >$@

sbmv:
	mkdir sbmv

sbmv/Makefile: sbmv Makefile.template ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=sbmv Makefile.template >$@

hemv:
	mkdir hemv

hemv/Makefile: hemv Makefile.template ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=hemv Makefile.template >$@

hpmv:
	mkdir hpmv

hpmv/Makefile: hpmv Makefile.template ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=hpmv Makefile.template >$@

hbmv:
	mkdir hbmv

hbmv/Makefile: hbmv Makefile.template ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=hbmv Makefile.template >$@

trmv:
	mkdir trmv

trmv/Makefile: trmv Makefile.template ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=trmv Makefile.template >$@

tpmv:
	mkdir tpmv

tpmv/Makefile: tpmv Makefile.template ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=tpmv Makefile.template >$@

trsv:
	mkdir trsv

trsv/Makefile: trsv Makefile.template ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=trsv Makefile.template >$@

tbsv:
	mkdir tbsv

tbsv/Makefile: tbsv Makefile.template ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=tbsv Makefile.template >$@

gemm:
	mkdir gemm

gemm/Makefile: gemm Makefile.template ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=gemm Makefile.template >$@

symm:
	mkdir symm

symm/Makefile: symm Makefile.template ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=symm Makefile.template >$@

hemm:
	mkdir hemm

hemm/Makefile: hemm Makefile.template ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=hemm Makefile.template >$@

gemv2:
	mkdir gemv2

gemv2/Makefile: gemv2 Makefile.template ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=gemv2 Makefile.template >$@

symv2:
	mkdir symv2

symv2/Makefile: symv2 Makefile.template ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=symv2 Makefile.template >$@

hemv2:
	mkdir hemv2

hemv2/Makefile: hemv2 Makefile.template ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=hemv2 Makefile.template >$@

gbmv2:
	mkdir gbmv2

gbmv2/Makefile: gbmv2 Makefile.template ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=gbmv2 Makefile.template >$@


common/Makefile: common/Makefile.m4 ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) common/Makefile.m4 >$@

clean:
	rm -f temp *~

maintainer-clean: clean
	rm -f  dot/Makefile sum/Makefile axpby/Makefile waxpby/Makefile gemv/Makefile ge_sum_mv/Makefile gbmv/Makefile symv/Makefile spmv/Makefile sbmv/Makefile hemv/Makefile hpmv/Makefile hbmv/Makefile trmv/Makefile tpmv/Makefile trsv/Makefile tbsv/Makefile gemm/Makefile symm/Makefile hemm/Makefile gemv2/Makefile symv2/Makefile hemv2/Makefile gbmv2/Makefile common/Makefile

