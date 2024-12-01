include ../../config.h

.SUFFIXES:
.SUFFIXES:  .nab .o .c

%.o: %.nab $(BINDIR)/$(MPI)nab$(SFX) $(BINDIR)/$(MPI)nab2c$(SFX)
	@echo "[NAB]  NAB $<"
	$(VB) $(BINDIR)/$(MPI)nab$(SFX) -c $< -o $@
	@rm -f $*.c

NAB2COBJS=	\
	nab2c.o	\
	cgen.o	\
	checkexpr.o	\
	dumpnode.o\
	errormsg.o	\
	fixexpr.o	\
	node.o	\
	symbol.o	\
	nabgrm.tab.o

LIBNABOBJS=	\
	axis2frame.o	\
	bonded_atoms.o \
	boundsutil.o	\
	chirvol.o \
	database.o  \
	deform_ener.o \
	dump.o		\
    embed.o     \
	fd_helix.o  \
	frameutil.o	\
	getdist.o	\
	gsub.o	\
	hashutil.o	\
	jacobi.o	\
	match.o		\
	matop.o		\
	molio.o		\
	molsurf.o   \
	molutil.o	\
	pointutil.o	\
	prm.o \
	regexp.o	\
	reslib.o	\
	select_atoms.o	\
	traceback.o \
	mme_init.o \
	stringutil.o	\
	temputil.o	\
    time.o \
	tm_malloc.o	\
	bdna.o	\
	dg_helix.o \
	getchivol.o \
	getseq_from_pdb.o \
	getxyz_from_pdb.o \
	getpdb_prm.o \
	linkprot.o	\
	link_na.o	\
	na_anal.o 	\
	setboundsfromdb.o \
	setchiplane.o	\
	wc_basepair.o	\
	wc_complement.o	\
	wc_helix.o

#===========================================================================

install_nab: nabhome.h $(BINDIR)/$(MPI)nab$(SFX) \
         $(BINDIR)/$(MPI)nab2c$(SFX) $(LIBDIR)/libnab.a
	@cp -p  nabcode.h nabtypes.h defreal.h nab.h nab_netcdf.h \
               xmin_opt.h lmod_opt.h $(INCDIR)
#
#   Following for nab.MPI, not yet supported here:
#	@mv $(BINDIR)/$(MPI)nab$(SFX) $(BINDIR)/to_be_dispatched/
#	@cat nab_template | sed "s/REPLACE_ME/$(MPI)nab$(SFX)/g" > $(BINDIR)/$(MPI)nab$(SFX)
#	@chmod +x $(BINDIR)/$(MPI)nab$(SFX)
#	cp lmodprmtop $(BINDIR)

#  Force rebuilding of nab, since there is no easy to to see if FLIBS has
#  changed.
$(BINDIR)/$(MPI)nab$(SFX)::
	@echo "[NAB]  CC $@"
	$(VB)$(CC) -DCC='"$(CC)"' -DCPP='"$(CPP)"' -DFLIBS='"$(FLIBS)"' \
		$(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) \
		$(LDFLAGS) $(AMBERLDFLAGS) -o $(BINDIR)/$(MPI)nab$(SFX) nab.c

$(BINDIR)/$(MPI)nab2c$(SFX):	$(NAB2COBJS)
	@echo "[NAB]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
		-o $(BINDIR)/$(MPI)nab2c$(SFX) $(NAB2COBJS)

#===========================================================================

$(LIBDIR)/libnab.a: $(LIBNABOBJS) $(BINDIR)/$(MPI)nab$(SFX) $(BINDIR)/$(MPI)nab2c$(SFX)
	@echo "[NAB]  CC $@"
	$(VB)-rm -f $@
	$(VB)$(AR) $@ $(LIBNABOBJS)  > /dev/null
	$(VB)-$(RANLIB) $@

nabhome.h::
	@echo "#define NABHOME \"$(BASEDIR)\"" > nabhome.h

nabgrm.tab.h: nabgrm.tab.c

nabgrm.tab.c:  nabgrm.y lex.yy.c
	@echo "[NAB]  BISON $<"
	$(VB)bison -v -d nabgrm.y

nabgrm.tab.o: nabgrm.tab.c
	@echo "[NAB]  CC $<"
	$(VB)$(CC) -c -D$(LEX) $(CFLAGS) $(AMBERCFLAGS) $<

lex.yy.c: nablex.l
	@echo "[NAB]  FLEX $<"
	$(VB)flex nablex.l

lex.dgo.c: dg_options.l
	@echo "[NAB]  FLEX $<"
	$(VB)$(LEX) dg_options.l

embed.o:  embed.c
	@echo "[NAB]  CC $<"
	$(VB)$(CC) -c -D$(LEX) $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) embed.c

checkexpr.c: nabgrm.tab.o
	(cd semantics; $(MAKE))
	mv semantics/checkexpr.c .

.c.o:
	@echo "[NAB]  CC $<"
	$(VB)$(CC) -c $(CFLAGS) $(AMBERCFLAGS) -o $@ $<

#    Compiles with optimized code:

molio.o:  molio.c
	@echo "[NAB]  CC $<"
	$(CC) -c $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) molio.c

mme_init.o:  mme_init.c
	@echo "[NAB]  CC $<"
	$(VB)$(CC) -c -D$(LEX) $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) mme_init.c

boundsutil.o:  boundsutil.c
	@echo "[NAB]  CC $<"
	$(VB)$(CC) -c $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) boundsutil.c

bonded_atoms.o:  bonded_atoms.c
	@echo "[NAB]  CC $<"
	$(VB)$(CC) -c $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) bonded_atoms.c

#===========================================================================

clean:
	-rm -f *.o nab$(SFX) nab2c$(SFX) libnab.a checkexpr.c nabhome.h
	-rm -f nabgrm.tab.c nabgrm.tab.h lex.dgo.c lex.yy.c 
	(cd semantics && $(MAKE) clean )

uninstall: clean
	-rm -f \
	$(BINDIR)/to_be_dispatched/nab$(SFX) \
	$(BINDIR)/to_be_dispatched/mpinab$(SFX) \
	$(BINDIR)/nab$(SFX) \
        $(BINDIR)/nab2c$(SFX) \
	$(BINDIR)/mpinab$(SFX) \
        $(BINDIR)/mpinab2c$(SFX) \
	$(LIBDIR)/libnab.a \
	$(INCDIR)/nabcode.h \
	$(INCDIR)/nabtypes.h \
	$(INCDIR)/defreal.h \
        $(INCDIR)/nab.h \
	$(INCDIR)/xmin_opt.h \
	$(INCDIR)/lmod_opt.h \
        $(BINDIR)/lmodprmtop

examples: serial 
	@for i in $(EXAMPLES); do $(BINDIR)/nab $(NABFLAGS) -o `ls $$i | sed 's/.nab/$(SFX)/'` $$i; done
	(cd ../examples; $(MAKE) examples )

#===========================================================================

axis2frame.o:	\
	nab.h
binpos.o:  \
	memutil.h
bonded_atoms.o: \
	nab.h
boundsutil.o:	\
	nab.h  \
	memutil.h \
    traceback.h \
	errormsg.h
cgen.o:	\
	nab.h	\
	errormsg.h	\
	nabgrm.tab.h
chirvol.o:  \
	nab.h
database.o:  \
	database.h
dump.o:		\
	nab.h
dumpnode.o:	\
	nab.h	\
	nabgrm.tab.h
embed.o: \
	nab.h  \
	errormsg.h  \
	memutil.h  \
	traceback.h  \
	chirvol.h  \
	lex.dgo.c 
errormsg.o:	\
	errormsg.h
fixexpr.o:	\
	nab.h	\
	errormsg.h	\
	nabgrm.tab.h
frameutil.o:	\
	nab.h	\
	memutil.h
genavscode.o:	\
	nab.h	\
	nabgrm.tab.h
gsub.o:    \
	nab.h
hashutil.o:	\
	nab.h
jacobi.o:	\
	memutil.h
matop.o: \
	nabcode.h
match.o:	\
	nab.h
molio.o:	\
	nab.h	\
	errormsg.h
molsurf.o:  \
	nab.h   \
	molsurf.h
molutil.o:	\
	nab.h	\
	errormsg.h	\
	memutil.h
nab2c.o:	\
	nab.h
node.o:	\
	nab.h	\
	errormsg.h	\
	nabgrm.tab.h
pointutil.o: \
	nab.h
prm.o: \
	nab.h
reslib.o: \
	nab.h \
	errormsg.h \
	memutil.h  \
	database.h
select_atoms.o:	\
	nab.h
mme_init.o: \
	nab.h 
stringutil.o: \
	nab.h
symbol.o:	\
	nab.h	\
	errormsg.h
temputil.o: \
	nab.h
traceback.o: \
	errormsg.h
nabgrm.tab.o:	\
	nab.h	\
	errormsg.h \
	lex.yy.c

