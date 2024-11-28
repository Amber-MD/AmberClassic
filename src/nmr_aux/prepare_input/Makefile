include ../../../config.h

$(BINDIR)/makeDIST_RST:	makeDIST_RST.o
	@echo "[NMR_AUX]  CC $@"
	$(VB)$(CC) $(CFLAGS) -o $@ makeDIST_RST.o \
		$(LDFLAGS) $(AMBERLDFLAGS) $(LM)

$(BINDIR)/nef_to_RST:	nef_to_RST.o $(LIBDIR)/libcifparse.a
	@echo "[NMR_AUX]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) -o $@ nef_to_RST.o \
		$(LDFLAGS) $(AMBERLDFLAGS) $(LM) -L$(LIBDIR) -lcifparse

$(BINDIR)/makeANG_RST:	makeANG_RST.o
	@echo "[NMR_AUX]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) -o $@ makeANG_RST.o \
		$(LDFLAGS) $(AMBERLDFLAGS) $(LM)

nmr_aux: 	$(BINDIR)/makeDIST_RST $(BINDIR)/makeANG_RST $(BINDIR)/nef_to_RST
	$(VB)cp makeCHIR_RST $(BINDIR)
	$(VB)cp map.NEF-AMBER map.DG-AMBER $(DATDIR)

uninstall:
	/bin/rm -f $(BINDIR)/makeDIST_RST \
	           $(BINDIR)/makeANG_RST $(BINDIR)/nef_to_RST
	/bin/rm -f $(BINDIR)/makeCHIR_RST $(BINDIR)/makeDIP_RST.cyana
	/bin/rm -f $(DATDIR)/map.NEF-AMBER $(DATDIR)/map.DG-AMBER

clean:
	/bin/rm -f *.o 

makeDIST_RST.o:  global.h
nef_to_RST.o:  global.h

#  default rule for C compilation:

.c.o:
	@echo "[NMR_AUX]  CC $<"
	$(VB)$(CC) -c $(CPPFLAGS) $(CNOOPTFLAGS) $(CFLAGS) -o $@ $<

$(LIBDIR)/libcifparse.a:
	cd ../../cifparse && $(MAKE) install
