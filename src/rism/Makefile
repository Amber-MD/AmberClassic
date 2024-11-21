SHELL=/bin/sh
include ../../config.h

# Add search paths for .mod files
LOCALFLAGS = -I../msander

# ------- rism3d information: ----------------------------------------------

LIBOBJ = array_util.o rism_util.o \
	rism3d_grid_c.o rism3d_closure_c.o \
	rism3d_kh_c.o rism3d_hnc_c.o rism3d_psen_c.o \
	rism3d_c.o rism3d_potential_c.o rism3d_restart.o \
	rism3d_xyzv.o rism3d_opendx.o rism3d_mrc.o \
	rism3d_solvent_c.o rism3d_solute_c.o \
	pubfft.o rism3d_fft.o rism_parm.o rism_io.o \
	mdiis_orig_c.o mdiis_blas_c.o mdiis_blas2_c.o mdiis_c.o \
	erfcfun.o safemem.o constants_rism.o \
	bspline.o  fftw3.o rism_report_c.o amber_rism_interface.o \
	rism1d_c.o solvmdl_c.o rism1d_potential_c.o\
	rism1d_closure_c.o rism1d_kh_c.o rism1d_hnc_c.o rism1d_py_c.o \
	rism1d_mv0_c.o rism1d_psen_c.o rism1d_nub_c.o rism1d_polyt_c.o \
	rism1d_devb_c.o rism_timer_c.o constants.o quaternion.o

librism: $(LIBDIR)/librism.a

rism1d: $(BINDIR)/rism1d

$(BINDIR)/rism1d: rism1d_driver.o $(LIBDIR)/librism.a $(LIBDIR)/libxblas-amb.a
	@echo "[RISM]  FC $@"
	$(VB)$(FC) $(FFLAGS) -o $@ rism1d_driver.o \
	        -L$(LIBDIR) -lxblas-amb $(FLIBSF)

$(LIBDIR)/librism.a: $(LIBOBJ)
	@echo "[RISM]  AR $@"
	$(VB)$(AR) $@ $(LIBOBJ) > /dev/null
	$(VB)$(RANLIB) $@

$(LIBDIR)/libxblas-amb.a:
	cd ../xblas && $(MAKE) lib-amb

uninstall:
	/bin/rm -f $(LIBDIR)/librism.a $(BINDIR)/rism1d

clean:
	-/bin/rm -f *.o *.mod

depend::
	./makedepend > depend

skip::
	@test -e amber_rism_interface.o || touch amber_rism_interface.o 

include depend

# DO NOT DELETE
