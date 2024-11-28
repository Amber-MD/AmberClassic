include ../../config.h

.c.o:
	@echo "[SFF]  CC $<"
	$(VB)$(CC) -c $(COPTFLAGS) $(CFLAGS) -o $@ $<

OBJS = conjgrad.o lmodC.o memutil.o nblist.o newton.o nmode.o wallclock.o \
	   prm.o rand2.o sasad.o sff.o xminC.o ptrajmask.o

sff: $(LIBDIR)/libsff.a
	$(VB)cp sff.h $(INCDIR)

$(LIBDIR)/libsff.a: $(OBJS)
	@echo "[SFF]  AR $@"
	$(VB)-rm -f $(LIBDIR)/libsff.a
	$(VB)$(AR) $(LIBDIR)/libsff.a $(OBJS)	> /dev/null
	$(VB)$(RANLIB) $(LIBDIR)/libsff.a

uninstall:
	-rm -f $(LIBDIR)/libsff.a
	-rm -f $(LIBDIR)/libsff_mpi.a
	-rm -f $(INCDIR)/sff.h

clean:
	-/bin/rm -f *.o lex.mmo.c 

lex.mmo.c: mm_options.l
	$(LEX) mm_options.l 

depend: $(OBJS:.o=.c) lex.mmo.c 
	@echo "[SFF]  make depend"
	$(VB)$(CC) -MM $^ > $@

include depend

