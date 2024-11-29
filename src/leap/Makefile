include ../../config.h
LOCALFLAGS=-Wno-unused-value

.c.o:
	@echo "[TLEAP] CC $<"
	$(VB)$(CC) -c -I$(XHOME)/include $(NETCDFINC) $(CFLAGS) $(AMBERCFLAGS) -o $@ $<

# 
#  Absolutely necessary routines
# 

OA1= basics.o sysdepend.o stringExtra.o varArray.o getline.o avl.o
OA = $(OA1)

# 
#  General tools
# 

OG1= pdb_format.o pdb_read.o pdb_sprntf.o 
OG2= pdb_sscanf.o pdb_write.o
OG3= vector.o zMatrix.o sort.o bag.o hash.o 
OG4= dictionary.o database.o nVector.o ring.o matrix.o
OG5= fortran.o displayer.o
OG = $(OG1) $(OG2) $(OG3) $(OG4) $(OG5)

# 
#  Object hierarchy objekts
# 

OO1= assoc.o atom.o byteArray.o 
OO2= collection.o container.o internal.o
OO3= list.o loop.o molecule.o 
OO4= oDouble.o oInteger.o oString.o
OO5= objekt.o parmSet.o residue.o unit.o unitio.o tripos.o
OO6= graphUtil.o select.o
OO = $(OO1) $(OO2) $(OO3) $(OO4) $(OO5) $(OO6)

# 
#  xaLeap and leap specific objects and tools
# 

OL1= amber.o build.o elements.o library.o chirality.o
OL2= minimizer.o model.o parmLib.o pdbFile.o mol2File.o mol3File.o tools.o
OL3= variables.o parser.o help.o helptext.o octree.o
OL4= commands.o mathop.o block.o restraint.o hybrid.o
OL = $(OL1) $(OL2) $(OL3) $(OL4)


# 
#  Terminal interface specifics
# 

OC1= tLeap.o
OC  = $(OC1)


#  utilities specifics
# 
#        OXUMF:  xaUtilMessageFilter
#        OUCWB:  utilConvertWaterBox
#        OUC2O:  utilCrd2Off
#        OUDT:   utilDescribeTopology
#        OUMH:   utilMakeHelp
#        OUT:	mol2_to_off
# 

OUCWB= utilConvertWaterBox.o
OUDT= utilDescribeTopology.o
OUT= mol2_to_off.o
OUC2O= utilCrd2Off.o
OUMH= utilMakeHelp.o
OUL2P= utilLib2Pdb.o

# 
#  Other dependancy files
# 

OBJ = $(OA) $(OG) $(OO) $(OL)

TLEAP_OBJ = $(OA) $(OG) $(OO) $(OL) $(OC)

UMF_OBJ = $(OA) $(OXUMF)
UMH_OBJ = $(OUMH)
UED_OBJ = $(OA) $(OG) $(OUED)
UCWB_OBJ = $(OA) $(OG) $(OO) $(OL) $(OUCWB)
UDT_OBJ = $(OA) $(OG) $(OO) $(OL) $(OUDT)
UT_OBJ = $(OA) $(OG) $(OO) $(OL) $(OUT)
UL2P_OBJ = $(OA) $(OG) $(OO) $(OL) $(OUL2P)
UC2O_OBJ = $(OA) $(OG) $(OO) $(OL) $(OUC2O)
UWR_OBJ = $(OA) $(OUWR)

T_OBJ = $(OA) $(OG) $(OO) $(OL) test_parmSet.o

UTILPROGS = utilMakeHelp$(SFX) utilDescribeTopology$(SFX) utilCrd2Off$(SFX) mol2_to_off$(SFX) utilConvertWaterBox$(SFX) utilLib2Pdb$(SFX)

LEAPPROGS = teLeap$(SFX)
PROGS = $(LEAPPROGS) $(UTILPROGS)

# ##############################################################
#
#  "make all" makes EVERYTHING, Wcl libraries, teLeap, xaLeap, util progs.
#  "make t" makes teLeap only
#  "make x" makes x windows stuff only
#  "make util" makes utility programs
#
#  "make all" must be used when the system is first ported to a
#           machine to make Wcl stuff.
#

util:   $(UTILPROGS)

leap: $(BINDIR)/teLeap$(SFX) $(BINDIR)/tleap

uninstall: clean
	-rm -f $(BINDIR)/teLeap$(SFX) $(BINDIR)/tleap

$(BINDIR)/teLeap$(SFX): $(TLEAP_OBJ)
	@echo "[TLEAP] CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(TLEAP_OBJ) $(LM) $(NETCDFLIB)

$(BINDIR)/tleap: tleap
	cp tleap $(BINDIR)

utilMakeHelp$(SFX): $(UMH_OBJ)
	@echo "[TLEAP] CC $<"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(UMH_OBJ) $(LM)

utilConvertWaterBox$(SFX): $(UCWB_OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(UCWB_OBJ) $(LM)

utilDescribeTopology$(SFX): $(UDT_OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(UDT_OBJ) $(LM)

utilCrd2Off$(SFX): $(UC2O_OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(UC2O_OBJ) $(LM)

mol2_to_off$(SFX): $(UT_OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(UT_OBJ) $(LM)

utilLib2Pdb$(SFX): $(UL2P_OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(UL2P_OBJ) $(LM)

clean::
	-rm -f *.o 
	-rm -f helptext.c parser.c

parser.c:       parser.y parser.h
	@echo "[TLEAP] BISON $@"
	$(VB)$(BISON) parser.y
#	@echo "(The warning 'conflicts: 6 shift/reduce' is normal)"
	$(VB)/bin/rm -f parser.c
	@mv parser.tab.c parser.c

helptext.c: helptext.text utilMakeHelp$(SFX)
	$(VB)./utilMakeHelp$(SFX) helptext.text helptext.c

depend: $(TLEAP_OBJ:.o=.c)
	@echo "[TLEAP] make depend"
	$(VB)$(CC) -MM $^ > $@

include depend
