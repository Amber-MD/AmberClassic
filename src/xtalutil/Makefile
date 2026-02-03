include ../../config.h

.SUFFIXES : .c .F90 .o

.c.o:
	@echo "[XTALUTIL]  CC $<"
	@$(CC) -c $(CFLAGS) -o $@ $< 

.F90.o:
	@echo "[XTALUTIL]  FC $<"
	@$(FC) -c $(FFLAGS) -o $@ $< 

OBJS = crdmanip.o grid.o matrix.o pdbRead.o rand2.o topRead.o vector.o

AOBJ=   ambpdb.o forFIRSTtype.o forFIRSTbond.o forFIRSThbond.o \
        forFIRSTteth.o forFIRSToutp.o forFIRSTsugar.o nxtsec.o \
        amopen.o mexit.o

xtalutil: $(BINDIR)/ChBox $(BINDIR)/UnitCell $(BINDIR)/CheckMD \
     $(BINDIR)/AddToBox $(BINDIR)/PropPDB $(BINDIR)/ambpdb
	cp ensemble_scripts/find_alts ensemble_scripts/make_addles.in $(BINDIR)
	cp phenix_scripts/XrayPrep phenix_scripts/XrayPrep_redo $(BINDIR)

$(BINDIR)/ambpdb: $(AOBJ)
	@echo "[XTALUTIL]  FC $@"
	@$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(AOBJ)

$(BINDIR)/AddToBox: AddToBox.o $(OBJS)
	@echo "[XTALUTIL]  CC $@"
	@$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lm

$(BINDIR)/ChBox: ChBox.o $(OBJS)
	@echo "[XTALUTIL]  CC $@"
	@$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lm

$(BINDIR)/CheckMD: CheckMD.o $(OBJS)
	@echo "[XTALUTIL]  CC $@"
	@$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lm

$(BINDIR)/PropPDB: PropPDB.o $(OBJS)
	@echo "[XTALUTIL]  CC $@"
	@$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lm

$(BINDIR)/UnitCell: UnitCell.o $(OBJS)
	@echo "[XTALUTIL]  CC $@"
	@$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lm

$(BINDIR)/pdb2rdf: pdb2rdf.o $(OBJS)
	@echo "[XTALUTIL]  CC $@"
	@$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lm

uninstall: clean
	/bin/rm -f $(BINDIR)/AddToBox $(BINDIR)/ChBox \
       $(BINDIR)/PropPDB $(BINDIR)/UnitCell

clean:
	/bin/rm -f *.o

