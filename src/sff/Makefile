include ../../config.h

%.o: %.c
	@echo "[SFF]  CC $<"
	$(VB)$(CC) -c $(COPTFLAGS) $(CFLAGS) $(NETCDFINC) -o $@ $<

OBJS = conjgrad.o lmodC.o memutil.o nblist.o nmode.o \
	   prm.o rand2.o sasad.o sff.o xminC.o AmberNetcdf.o ptrajmask.o \
       wallclock.o

# install: libsff hcp_getpdb$(SFX)

libsff: $(LIBDIR)/libsff.a

$(LIBDIR)/libsff.a: $(OBJS)
	@echo "[SFF]  AR $@"
	$(VB)/bin/rm -f $(LIBDIR)/libsff.a
	$(VB)$(AR) $(LIBDIR)/libsff.a $(OBJS)	> /dev/null
	$(VB)$(RANLIB) $(LIBDIR)/libsff.a
	$(VB)cp sff.h AmberNetcdf.h $(INCDIR)

libsff_mpi: $(OBJS)
	@echo "[SFF]  AR $@"
	$(VB)/bin/rm -f $(LIBDIR)/libsff_mpi.a
	$(VB)$(AR) $(LIBDIR)/libsff_mpi.a $(OBJS)	> /dev/null
	$(VB)$(RANLIB) $(LIBDIR)/libsff_mpi.a
	$(VB)cp sff.h AmberNetcdf.h $(INCDIR)

uninstall:
	/bin/rm -f $(BINDIR)/hcp_getpdb$(SFX)
	/bin/rm -f $(LIBDIR)/libsff.a $(LIBDIR)/libsff_mpi.a
	/bin/rm -f $(INCDIR)/sff.h $(INCDIR)/AmberNetcdf.h

clean:
	/bin/rm -f $(OBJS) lex.mmo.c

lex.mmo.c: mm_options.l
	@echo "[SFF]  FLEX $<"
	$(VB)$(LEX) mm_options.l

hcp_getpdb$(SFX): hcp_getpdb.c
	@echo "[SFF]  CC $<"
	$(VB)$(CC) -o $(BINDIR)/hcp_getpdb$(SFX) hcp_getpdb.c

sff.c: lex.mmo.c

#depend: $(OBJS:.o=.c)
#	@echo "[SFF] make depend"
#	$(VB)$(CC) -MM $^ > $@

include depend

