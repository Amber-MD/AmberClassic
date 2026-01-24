include ../../config.h

.SUFFIXES:
.SUFFIXES: .F90 .o

SRC=\
	gbnsr6.F90 pb_read.F90  gb_read.F90 pb_write.F90 getcoor.F90 \
	runmd.F90 runmin.F90 force.F90 \
	pb_init.F90 sa_driver.F90 np_force.F90 variable_module.F90 \
	pb_exmol.F90 NSR6routines.F90  pb_list.F90 timer.F90 \
	egb.F90 ene.F90 locmem.F90 myopen.F90 rdparm.F90 \
	decomp.F90 rgroup.F90 rfree.F90 debug.F90 \
        svdcmp.F90 svbksb.F90 pythag.F90  memory_module.F90 \
        gen_dx_file.F90 parms.F90

OBJ=\
	gbnsr6.o pb_read.o  gb_read.o  pb_write.o getcoor.o \
	runmd.o runmin.o force.o\
	pb_init.o sa_driver.o np_force.o variable_module.o\
	 pb_exmol.o NSR6routines.o pb_list.o timer.o \
	egb.o ene.o locmem.o myopen.o rdparm.o \
	decomp.o rgroup.o rfree.o debug.o \
        svdcmp.o svbksb.o \
        parms.o pythag.o memory_module.o \
        gen_dx_file.o nxtsec.o

install_gbnsr6: $(BINDIR)/gbnsr6

$(BINDIR)/gbnsr6: $(OBJ)
	@echo "[GBNSR6]  FC $@"
	$(VB)$(FC) $(FFLAGS) -o $@ $(OBJ) -L$(LIBDIR) $(FLIBSF) 

.F90.o:   $<
	@echo "[GBNSR6]  FC $<"
	$(VB)$(FC) -c $(FOPTFLAGS) $(FFLAGS) -o $@ $<

nxtsec.o:  ../msander/nxtsec.F90
	@echo "[XRAY]  FC $<"
	$(VB)$(FC) -c $(FNOOPTFLAGS) $(FFLAGS) -o $@ $<

depend::
	./makedepend > depend

clean:
	/bin/rm -f *.o *.mod

uninstall: 
	-/bin/rm -f $(BINDIR)/gbnsr6 

distclean: clean uninstall



#----------HEADER DEPENDENCIES (currently from makedepend) 
include depend
# DO NOT DELETE
