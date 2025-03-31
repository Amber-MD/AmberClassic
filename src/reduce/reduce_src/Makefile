include ../../../config.h

.SUFFIXES: .c .cpp .o

.c.o:
	@echo "[REDUCE] CXX $<"
	$(VB)$(CC) $(CFLAGS) -c $< -o $@

.cpp.o:
	@echo "[REDUCE] CXX $<"
	$(VB)$(CXX) $(CFLAGS) -c $< -o $@

CFLAGS		= -O -DAROMATICS_ACCEPT_HBONDS -DCHARFUNCMACROS \
          -DBOOLPREDEFINED  -DBRACKETOPERPARMS -DLEFT_JUSTIFY_NUC_RES_OK \
          -DLISTFRIENDFIX -DINCTEMPLATEDEFNS \
          -DHET_DICTIONARY='"$(DATDIR)/reduce_wwPDB_het_dict.txt"' \
          -DHET_DICTOLD='"$(DATDIR)/reduce_het_dict.txt"' \
          -I$(TOOLHOME) -I$(PDBHOME)

OBJS = main.o GraphToHoldScores.o reduce.o CTab.o ElementInfo.o StdResH.o \
       ResBlk.o AtomConn.o AtomPositions.o DotSph.o Mover.o Rot.o \
       Rot3Fold.o RotAromMethyl.o RotDonor.o FlipMemo.o CliqueList.o \
       AtomDescr.o PDBrec.o MoveableNode.o hybrid_36_c.o

TOOLHOME = ../toolclasses
PDBHOME = ../libpdb

install_reduce: $(BINDIR)/reduce
	@gunzip -c ../reduce_wwPDB_het_dict.txt.gz > $(DATDIR)/reduce_wwPDB_het_dict.txt

$(BINDIR)/reduce:     $(OBJS)
	@echo "[REDUCE] CXX $@"
	$(VB)$(CXX) -o $@ $(CFLAGS) $(OBJS) $(LIBDIR)/toolclasses.a \
          $(LIBDIR)/pdb++.a $(LM)

clean:;		@rm -f $(OBJS)

# DO NOT DELETE THIS LINE -- make depend uses it
AtomConn.o: $(TOOLHOME)/Point3d.h $(TOOLHOME)/utility.cpp \
         $(TOOLHOME)/utility.h AtomConn.cpp \
	 AtomConn.h ElementInfo.h
AtomDescr.o: AtomDescr.cpp AtomDescr.h
AtomPositions.o: $(PDBHOME)/pdb++.h $(TOOLHOME)/DisjointSets.h \
	 $(TOOLHOME)/Point3d.h $(TOOLHOME)/utility.cpp \
	 $(TOOLHOME)/utility.h AtomConn.h AtomPositions.cpp \
	 AtomPositions.h BumperPoint.h CliqueList.h DotSph.h \
	 ElementInfo.h FlipMemo.h LocBlk.h Mover.h PDBrec.h ResBlk.h \
	 Rot.h RotDonor.h Rot3Fold.h StdResH.h neighbors.cpp neighbors.h \
	 AtomDescr.h MoveableNode.h GraphToHoldScores.h
CTab.o: $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.cpp $(TOOLHOME)/utility.h AtomConn.h \
	 CTab.cpp CTab.h ElementInfo.h
CliqueList.o: $(PDBHOME)/pdb++.h \
	 $(TOOLHOME)/Point3d.h $(TOOLHOME)/utility.cpp \
	 $(TOOLHOME)/utility.h BumperPoint.h CliqueList.cpp \
	 CliqueList.h DotSph.h ElementInfo.h LocBlk.h Mover.h PDBrec.h \
	 neighbors.cpp neighbors.h
DotSph.o: $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.cpp $(TOOLHOME)/utility.h DotSph.cpp \
	 DotSph.h
ElementInfo.o: $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.cpp $(TOOLHOME)/utility.h AtomConn.h \
	 ElementInfo.cpp ElementInfo.h StdResH.h
FlipMemo.o: $(PDBHOME)/pdb++.h $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.cpp $(TOOLHOME)/utility.h AtomConn.h \
	 AtomPositions.h BumperPoint.h CliqueList.h DotSph.h \
	 ElementInfo.h FlipMemo.cpp FlipMemo.h LocBlk.h Mover.h PDBrec.h \
	 ResBlk.h neighbors.cpp neighbors.h
GraphToHoldScores.o: GraphToHoldScores.h GraphToHoldScores.cpp AtomPositions.h \
    AtomDescr.h
hybrid_36_c.o: hybrid_36_c.c hybrid_36_c.h
Mover.o: $(PDBHOME)/pdb++.h $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.cpp $(TOOLHOME)/utility.h \
	 BumperPoint.h DotSph.h ElementInfo.h LocBlk.h Mover.cpp Mover.h \
	 PDBrec.h ResBlk.h neighbors.cpp neighbors.h
MoveableNode.o: MoveableNode.cpp MoveableNode.h AtomDescr.h \
    GraphToHoldScores.h AtomPositions.h
PDBrec.o: $(PDBHOME)/pdb++.h \
	 $(TOOLHOME)/Point3d.h $(TOOLHOME)/utility.cpp \
	 $(TOOLHOME)/utility.h ElementInfo.h AtomDescr.h PDBrec.cpp PDBrec.h hybrid_36_c.h
ResBlk.o: $(PDBHOME)/pdb++.h $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.cpp $(TOOLHOME)/utility.h \
	 ElementInfo.h PDBrec.h ResBlk.cpp ResBlk.h
Rot.o: $(PDBHOME)/pdb++.h $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.cpp $(TOOLHOME)/utility.h \
	 AtomPositions.h BumperPoint.h CliqueList.h DotSph.h \
	 ElementInfo.h FlipMemo.h LocBlk.h Mover.h PDBrec.h ResBlk.h \
	 Rot.cpp Rot.h neighbors.cpp neighbors.h
RotDonor.o: $(PDBHOME)/pdb++.h $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.cpp $(TOOLHOME)/utility.h \
	 AtomPositions.h BumperPoint.h CliqueList.h DotSph.h \
	 ElementInfo.h FlipMemo.h LocBlk.h Mover.h PDBrec.h ResBlk.h \
	 RotDonor.cpp RotDonor.h neighbors.cpp neighbors.h
Rot3Fold.o: $(PDBHOME)/pdb++.h $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.cpp $(TOOLHOME)/utility.h \
	 AtomPositions.h BumperPoint.h CliqueList.h DotSph.h \
	 ElementInfo.h LocBlk.h Mover.h PDBrec.h ResBlk.h Rot3Fold.cpp \
	 Rot3Fold.h neighbors.cpp neighbors.h
RotAromMethyl.o: $(PDBHOME)/pdb++.h $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.cpp $(TOOLHOME)/utility.h \
	 AtomPositions.h BumperPoint.h CliqueList.h DotSph.h \
	 ElementInfo.h LocBlk.h Mover.h PDBrec.h ResBlk.h RotAromMethyl.cpp \
	 RotAromMethyl.h neighbors.cpp neighbors.h
StdResH.o: $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.cpp $(TOOLHOME)/utility.h AtomConn.h \
	 ElementInfo.h StdResH.cpp StdResH.h
neighbors.o: $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.cpp $(TOOLHOME)/utility.h LocBlk.h \
	 neighbors.cpp neighbors.h
reduce.o: $(PDBHOME)/pdb++.h $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.cpp $(TOOLHOME)/utility.h AtomConn.h \
	 AtomPositions.h BumperPoint.h CTab.h CliqueList.h DotSph.h \
	 ElementInfo.h LocBlk.h Mover.h PDBrec.h ResBlk.h StdResH.h \
	 reduce.cpp
main.o: $(PDBHOME)/pdb++.h $(TOOLHOME)/Point3d.h \
	 $(TOOLHOME)/utility.cpp $(TOOLHOME)/utility.h AtomConn.h \
	 AtomPositions.h BumperPoint.h CTab.h CliqueList.h DotSph.h \
	 ElementInfo.h LocBlk.h Mover.h PDBrec.h ResBlk.h StdResH.h \
	 main.cpp
# DO NOT DELETE THIS 2nd LINE -- make depend uses it
