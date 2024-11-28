dnl M4 file for generating src/Makefile.
dnl
include(../m4/cblas.m4)dnl
dnl
dnl
`include' ../make.conf
`include' ../$(MAKEINC)

MKFILES = Makefile`'FOREACH(`LIBNAMES', ` arg/Makefile') common/Makefile

dnl Can't use 'MAKEFILES' since it is used by GNU make.
makefiles: $(MKFILES)

Makefile: Makefile.m4 ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) Makefile.m4 >$`'@

FOREACH(`LIBNAMES', `
arg:
	mkdir arg

arg/Makefile: arg Makefile.template ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) -Droutine=arg Makefile.template >$`'@
')

common/Makefile: common/Makefile.m4 ../m4/cblas.m4 ../$(MAKEINC)
	$(M4) $(M4_OPTS) common/Makefile.m4 >$`'@

clean:
	rm -f temp *~

maintainer-clean: clean
	rm -f FOREACH(`LIBNAMES', ` arg/Makefile') common/Makefile

