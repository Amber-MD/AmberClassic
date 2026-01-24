include ../../config.h

PROGS=	$(BINDIR)/pol_h $(BINDIR)/protonate $(BINDIR)/gwh

install_protonate:	$(PROGS)

uninstall:  clean
		-rm -f $(BINDIR)/pol_h
		-rm -f $(BINDIR)/protonate
		-rm -f $(BINDIR)/gwh

LIBOBJ= mexit.o amopen.o

SRCPH=	pol_h.F stm_sel.F sortf_as.F pdbatm.F
OBJPH=	pol_h.o stm_sel.o sortf_as.o pdbatm.o

$(BINDIR)/pol_h: 		$(OBJPH) $(LIBOBJ)
	@echo "[PROTONATE]  FC $@"
	$(VB)$(FC) $(FFLAGS) -o $@ $(OBJPH) $(LIBOBJ)
        

$(BINDIR)/protonate:  	protonate.o $(LIBOBJ)
	@echo "[PROTONATE]  FC $@"
	$(VB)$(FC) $(FFLAGS) -o $@ protonate.o $(LIBOBJ)

OBJGW=	gwh.F pdbatm.F pdbhet2.F top_3.F stm_sel.F sortf_as.F
OBJGW=	gwh.o pdbatm.o pdbhet2.o top_3.o stm_sel.o sortf_as.o

$(BINDIR)/gwh:		$(OBJGW) $(LIBOBJ)
	@echo "[PROTONATE]  FC $@"
	$(VB)$(FC) $(FFLAGS) -o $@ $(OBJGW) $(LIBOBJ)

clean:
		/bin/rm -f *.o 

protonate.o:  protonate.h

.F.o:   $<
	@echo "[PROTONATE]  FC $<"
	$(VB)$(FC) -c $(FOPTFLAGS) $(FFLAGS) -o $@ $<

