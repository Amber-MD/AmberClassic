include ../../config.h

CUSTOMBUILDFLAGS = -I../sff

.SUFFIXES:
.SUFFIXES:  .c .o 

.c.o:
	@echo "[NABC]  CC $<"
	$(VB)$(CC) -c $(COPTFLAGS) $(CFLAGS) $(NETCDFINC) -o $@ $<

NABCOBJS=	\
	axis2frame.o	\
	bonded_atoms.o \
	deform_ener.o \
    fd_helix.o \
	frameutil.o	\
	getdist.o	\
	jacobi.o	\
	match.o		\
	matop.o		\
	molio.o		\
	molsurf.o   \
    molutil.o   \
	pointutil.o	\
	prm.o \
	regexp.o \
	select_atoms.o	\
    stringutil.o \
	traceback.o \
    wc_basepair.o \
    wc_complement.o \
    wc_helix.o

#===========================================================================

nabc: $(LIBDIR)/libnabc.a

$(LIBDIR)/libnabc.a: $(NABCOBJS)
	@echo "[NABC]  AR $@"
	$(VB)/bin/rm -f $@
	$(VB)$(AR) $@ $(NABCOBJS) > /dev/null
	$(VB)$(RANLIB) $@
	$(VB)cp nabc.h $(INCDIR)

lex.dg_options.c: dg_options.l
	@echo "[NABC]  $(LEX) $<"
	$(VB)$(LEX) -t dg_options.l | sed -e 's/{stdout}/{DONTchangeSTDOUThere}/' -e 's/yy/dgo/g' -e 's/stdout/nabout/' -e 's/DONTchangeSTDOUThere/stdout/' > lex.dg_options.c

#===========================================================================

clean:
	-rm -f *.o

uninstall: clean
	-rm -f  $(LIBDIR)/libnabc.a  $(INCDIR)/nabc.h
	-rm -f *.o

#===========================================================================

depend: $(NABCOBJS:.o=.c)
	@echo "[NABC] make depend"
	$(VB)$(CC) -MM $(CFLAGS) $^ > $@

include depend
