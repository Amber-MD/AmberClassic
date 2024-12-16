include ../../../config.h

.SUFFIXES:	.cpp .o

.cpp.o:
	@echo "[REDUCE] CXX $<"
	$(VB)$(CXX) $(CFLAGS) -c $< -o $@

CFLAGS		= -O -DBOOLPREDEFINED -DHASSSCANFEXTERN

OBJS		= pdb_read.o pdb_sprntf.o pdb_sscanf.o pdb_chars.o \
		  pdb_type.o pdb++.o pdbinput.o

libpdb: $(LIBDIR)/pdb++.a

$(LIBDIR)/pdb++.a:     $(OBJS)
		@echo "[REDUCE] AR $@"
		$(VB)$(AR) $@ $(OBJS)
		$(VB)$(RANLIB) $@

clean:;		@rm -f $(OBJS)

# DO NOT DELETE THIS LINE -- make depend uses it
pdb++.o: pdb++.cpp pdb++.h
pdb_chars.o: pdb++.h pdb_chars.cpp pdbrun5_write.i pdbrun6_write.i \
	 write_format.i
pdb_read.o: pdb++.h pdb_read.cpp pdbrun5_read.i pdbrun6_read.i \
	 read_format.i
pdb_sprntf.o: pdb++.h pdb_sprntf.cpp
pdb_sscanf.o: pdb++.h pdb_sscanf.cpp
pdb_type.o: pdb++.h pdb_type.cpp
pdbinput.o: pdb++.h pdbinput.cpp
# DO NOT DELETE THIS 2nd LINE -- make depend uses it
