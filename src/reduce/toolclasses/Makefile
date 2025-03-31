include ../../../config.h

.SUFFIXES:	.cpp .o

.cpp.o:
	@echo "[REDUCE] CXX $<"
	$(VB)$(CXX) $(CFLAGS) -c $< -o $@

CFLAGS		= -O -DCHARFUNCMACROS -DINCTEMPLATEDEFNS \
		   -DBOOLPREDEFINED -DBRACKETOPERPARMS -DLISTFRIENDFIX

LIBARCH		= lib$(LIBRARY).a

OBJS		= UseCount.o Point3d.o utility.o \
		  column_sprintf.o column_sscanf.o DisjointSets.o

toolclasses: $(LIBDIR)/toolclasses.a

$(LIBDIR)/toolclasses.a:     $(OBJS)
		@echo "[REDUCE] AR $@"
		$(VB)$(AR) $@ $(OBJS)
		$(VB)$(RANLIB) $@

clean:;		@rm -f $(OBJS)

# DO NOT DELETE THIS LINE -- make depend uses it
DisjointSets.o: DisjointSets.cpp DisjointSets.h
Point3d.o: Point3d.cpp Point3d.h
UseCount.o: UseCount.cpp UseCount.h
column_sprintf.o: column_sprintf.cpp
column_sscanf.o: column_sscanf.cpp
utility.o: utility.cpp
# DO NOT DELETE THIS 2nd LINE -- make depend uses it
