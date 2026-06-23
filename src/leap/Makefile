include ../../config.h

CIFPARSEINC = -I../cifparse
CIFPARSE = $(LIBDIR)/libcifparse.a

.c.o:
	@echo "[TLEAP]  CC $<"
	$(VB)$(CC) -c -I$(XHOME)/include $(NETCDFINC) $(CIFPARSEINC) $(CFLAGS) -o $@ $<

# 
#  Absolutely necessary routines
# 

OA1= basics.o sysdepend.o stringExtra.o varArray.o getline.o avl.o
OA = $(OA1)

# 
#  General tools
# 

OG1= pdb_format.o pdb_read.o pdb_sprntf.o hybrid36.o neighbors.o symmetry.o
OG2= pdb_sscanf.o pdb_write.o cifFile.o
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
OO5= objekt.o parmSet.o residue.o unit.o unitio.o unitio_prmtop.o tripos.o
OO6= graphUtil.o select.o
OO = $(OO1) $(OO2) $(OO3) $(OO4) $(OO5) $(OO6)

# 
#  leap specific objects and tools
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

OXUMF= mess.o function.o
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

XALEAP_OBJ = $(OA) $(OG) $(OO) $(OL) $(OXT) $(OXAW)
XALEAP_DEP = $(XBITMAPS) $(ADDDEPS)

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

# SRCS= StrToPmap.c amber.c assoc.c atom.c bag.c basics.c block.c build.c byteArray.c chirality.c collection.c commands.c container.c database.c dictionary.c displayer.c elements.c fortran.c function.c getline.c avl.c graphUtil.c hash.c help.c hybrid.c internal.c library.c list.c loop.c mathop.c matrix.c mess.c minimizer.c model.c molecule.c nVector.c oDouble.c oInteger.c oString.c objekt.c octree.c parmLib.c parmSet.c pdbFile.c mol2File.c mol3File.c pdb_format.c pdb_read.c pdb_sprntf.c pdb_sscanf.c pdb_write.c residue.c restraint.c ring.c select.c sort.c stringExtra.c sysdepend.c tLeap.c test_dataBase.c test_parmSet.c tools.c tripos.c unit.c unitio.c utilConvertWaterBox.c utilCrd2Off.c utilDescribeTopology.c utilFixOff.c utilLib2Pdb.c utilMakeHelp.c mol2_to_off.c utilWrapRes.c varArray.c variables.c vector.c x3d.c xAction.c xBasics.c xTank.c xaAngleParmTable.c xaAtomParmTable.c xaAtomTable.c xaBondParmTable.c xaCommand.c xaHBondParmTable.c xaImproperParmTable.c xaLeapc.c xaParmEditor.c xaTable.c xaTools.c xaTorsionParmTable.c xaUnitEditor.c xaUtilMessageFilter.c XrawRegistr.c zMatrix.c

# HEADERS = amber.h assoc.h atom.h bag.h basics.h block.h build.h byteArray.h chirality.h classes.h collection.h commands.h config.h container.h database.h defaults.h dictionary.h displayer.h elements.h fortran.h function.h getline.h avl.h avl2.h graphUtil.h hash.h help.h hybrid.h internal.h leap.h library.h list.h loop.h machine.h mathop.h matrix.h mess.h minimizer.h model.h molecule.h nVector.h oDouble.h oInteger.h oString.h objekt.h octree.h parmLib.h parmSet.h parser.h pdb.h pdbFile.h mol2File.h mol3File.h pdb_int.h residue.h restraint.h ring.h select.h sort.h stringExtra.h sysdepend.h threed.h tools.h tripos.h unit.h unitio.h varArray.h vector.h x3d.h xAction.h xDisplayList.h xTank.h xTankP.h xaAngleParmTable.h xaAtomParmTable.h xaAtomTable.h xaBondParmTable.h xaCommand.h xaHBondParmTable.h xaLeap.h xaResources.h xaTable.h xaTools.h xaTorsionParmTable.h xaUnitEditor.h XrawRegistr.h zMatrix.h

# ##############################################################
#
util:   $(UTILPROGS)

obj:	$(OBJ)

leap: $(BINDIR)/teLeap$(SFX) $(BINDIR)/tleap

uninstall: clean
	-rm -f $(BINDIR)/teLeap$(SFX) $(BINDIR)/tleap

$(BINDIR)/teLeap$(SFX): $(TLEAP_OBJ)
	@echo "[TLEAP]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
	      -o $@ $(TLEAP_OBJ) $(LM) $(NETCDFLIB) $(CIFPARSE)

$(BINDIR)/tleap: tleap
	cp tleap $(BINDIR)

utilMakeHelp$(SFX): $(UMH_OBJ)
	@echo "[TLEAP]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(UMH_OBJ) $(LM)

utilConvertWaterBox$(SFX): $(UCWB_OBJ)
	@echo "[TEAP]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(UCWB_OBJ) $(LM)

utilDescribeTopology$(SFX): $(UDT_OBJ)
	@echo "[TLEAP]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(UDT_OBJ) $(LM)

utilCrd2Off$(SFX): $(UC2O_OBJ)
	@echo "[TLEAP]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(UC2O_OBJ) $(LM)

mol2_to_off$(SFX): $(UT_OBJ)
	@echo "[TLEAP]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(UT_OBJ) $(LM)

utilLib2Pdb$(SFX): $(UL2P_OBJ)
	@echo "[TLEAP]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(UL2P_OBJ) $(LM)

clean::
	-rm -f *.o $(PROGS) *\~
	-rm -f helptext.c parser.c

parser.o:       parser.y parser.h
	$(YACC) parser.y
	@echo "(The warning 'conflicts: 6 shift/reduce' is normal)"
	@/bin/rm -f parser.c
	@mv y.tab.c parser.c
	@echo "[TLEAP]  CC parser.c"
	$(VB)$(CC) -c -I$(XHOME)/include $(NETCDFINC) -DYYERROR_VERBOSE \
	      $(CFLAGS) $(LDFLAGS) -o $@ parser.c

helptext.c: helptext.text utilMakeHelp$(SFX)
	$(VB)./utilMakeHelp$(SFX) helptext.text helptext.c

# depend::
#	./makedepend > depend

depend: $(TLEAP_OBJ:.o=.c)
	@echo "[TLEAP] make depend"
	$(VB)$(CC) -MM $(CIFPARSEINC) $^ > $@

include depend
