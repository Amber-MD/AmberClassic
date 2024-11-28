include ../../config.h

saxs: $(BINDIR)/saxs_rism $(BINDIR)/saxs_md

lib: $(LIBDIR)/libsaxs_md.so

$(BINDIR)/saxs_rism: saxs_rism.o lebedev.o atm_f.o dx.o f_calc.o IO.o pdb.o vol_decomp.o
	@echo "[SAXS]  CXX $@"
	$(VB)$(CXX) $(COPTFLAGS) $(CFLAGS) $(LDFLAGS) \
	   -o $@ saxs_rism.o lebedev.o atm_f.o dx.o f_calc.o IO.o pdb.o vol_decomp.o

$(BINDIR)/saxs_md: saxs_md.o lebedev.o atm_f.o dx.o f_calc.o IO.o pdb.o
	@echo "[SAXS]  CXX $@"
	$(VB)$(CXX) $(COPTFLAGS) $(CFLAGS) $(LDFLAGS) \
	   -o $@ saxs_md.o lebedev.o atm_f.o dx.o f_calc.o IO.o pdb.o

libsaxs_md.so: saxs_md.o lebedev.o atm_f.o dx.o f_calc.o IO.o pdb.o
	@echo "[SAXS]  CXX $@"
	$(VB)$(CXX) -shared $(COPTFLAGS) $(CFLAGS) $(LDFLAGS) \
		-o $@ saxs_md.o lebedev.o atm_f.o dx.o f_calc.o IO.o pdb.o

libsaxs_rism.so: saxs_rism.o lebedev.o atm_f.o dx.o f_calc.o IO.o pdb.o vol_decomp.o
	@echo "[SAXS]  CXX $@"
	$(VB)$(CXX) -shared $(COPTFLAGS) $(CFLAGS) $(LDFLAGS) \
		-o $@ saxs_rism.o lebedev.o atm_f.o dx.o f_calc.o IO.o pdb.o vol_decomp.o

.cpp.o:
	@echo "[SAXS]  CXX $<"
	$(VB)$(CXX) -c $(CXXOPTFLAGS) $(CXXFLAGS) $(AMBERCXXFLAGS) -I$(INCDIR) -o $@ $<
	
#lebedev.cpp:  lebedev.hpp

clean:
	/bin/rm -f *.o

uninstall:
	/bin/rm -f $(BINDIR)/saxs_rism $(BINDIR)/saxs_md
	/bin/rm -f $(LIBDIR)/libsaxs_md.so
	/bin/rm -f $(LIBDIR)/libsaxs_rism.so
