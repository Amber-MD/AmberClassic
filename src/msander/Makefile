
include ../../config.h

SHELL=/bin/sh
VPATH=

# Add search paths for .mod files
LOCALFLAGS = $(CUDAFFLAGS) -I../rism -I../xray 

NFE_OBJECTS = \
    nfe-abmd-hooks.o nfe-bbmd-ctxt.o nfe-bbmd-hooks.o  \
    nfe-colvar.o nfe-colvar-math.o nfe-colvar-type.o nfe-colvar-utils.o \
    nfe-constants.o nfe-cv-ANGLE.o nfe-cv-COM_ANGLE.o \
    nfe-cv-COM_DISTANCE.o nfe-cv-COM_TORSION.o nfe-cv-COS_OF_DIHEDRAL.o \
    nfe-cv-DISTANCE.o nfe-cv-HANDEDNESS.o nfe-cv-LCOD.o \
    nfe-cv-MULTI_RMSD.o nfe-cv-N_OF_BONDS.o nfe-cv-N_OF_STRUCTURES.o \
    nfe-cv-R_OF_GYRATION.o nfe-cv-TORSION.o  \
    nfe-pmd-hooks.o nfe-rmsd.o nfe-sander-hooks.o nfe-sander-proxy.o \
    nfe-stsm-hooks.o nfe-cv-PCA.o nfe-read-pca.o nfe-cv-SIN_OF_DIHEDRAL.o \
    nfe-cv-PAIR_DIHEDRAL.o nfe-cv-PATTERN_DIHEDRAL.o nfe-cv-DF_COM_DISTANCE.o \
    nfe-smd-hooks.o nfe-umbrella.o nfe-utils.o nfe-value.o mt19937.o \
    nfe-cv-ORIENTATION_ANGLE.o nfe-cv-ORIENTATION_PROJ.o \
	nfe-cv-QUATERNION0.o nfe-cv-QUATERNION1.o nfe-cv-QUATERNION2.o \
	nfe-cv-QUATERNION3.o nfe-cv-SPINANGLE.o nfe-cv-TILT.o

MMOBJ = \
    align.o AmberNetcdf.o amopen.o assert.o \
    barostats.o binrestart.o bintraj.o charmm.o \
    constants.o csa.o cshf.o \
    debug.o decnvh.o degcnt.o dssp.o \
    dynlib.o egb.o emap.o ene.o erfcfun.o \
    ew_box.o ew_bspline.o ew_fft.o \
    ew_force.o ew_recip.o ew_setup.o extra_pts.o \
    fastwt.o file_io_dat.o findmask.o force.o getcor.o \
    icosasurf.o iwrap2.o les.o linear_response.o \
    lmodC.o lmod.o locmem.o matinv.o mbar.o \
    mdfil.o mdread.o md_scheme.o mdwrit.o memory_module.o \
    mexit.o minrit.o molecule.o mtmdcall.o multisander.o \
    multitmd.o music_module.o new_time.o nmlsrc.o nmrcal.o \
    nmr.o nonbond_list.o parallel.o parms.o pcshift.o \
    pearsn.o printe.o quench.o rand2.o \
    random.o rattle.o rdparm.o relax_mat.o \
    remd.o rfree.o rgroup.o rmsgrd.o runmd.o \
    sander_lib.o sander.o set.o get_atomic_number.o utilitiesModule.o \
    sglds.o shake.o $(SHORT_ENE) softcore.o \
    stack.o state.o thermo_int.o trajene.o \
    wallclock.o xminC.o veclib.o dsarpack.o nxtsec.o sinr_t.o

LESOBJ= \
    align.o AmberNetcdf.o amopen.o \
    assert.o barostats.o binrestart.LES.o bintraj.LES.o charmm.o \
    constants.o csa.o \
    cshf.o debug.o decnvh.o degcnt.LES.o \
    dssp.o dynlib.LES.o egb.LES.o emap.o ene.LES.o \
    erfcfun.o ew_box.LES.o ew_bspline.o \
    ew_fft.o ew_force.LES.o ew_recip.LES.o ew_setup.LES.o \
    extra_pts.LES.o fastwt.o file_io_dat.o findmask.o force.LES.o \
    get_atomic_number.o utilitiesModule.o getcor.LES.o icosasurf.o iwrap2.o \
    les.o linear_response.o lmodC.o lmod.LES.o locmem.LES.o \
    matinv.o mbar.o mdfil.LES.o mdread.LES.o md_scheme.LES.o \
    mdwrit.LES.o memory_module.o mexit.o minrit.LES.o molecule.o \
    mtmdcall.o multisander.LES.o multitmd.o music_module.o new_time.o \
    nmlsrc.o nmrcal.o nmr.LES.o nonbond_list.LES.o parallel.LES.o \
    parms.o pcshift.o pearsn.o printe.o \
    quench.o rand2.o random.o \
    rattle.o rdparm.LES.o relax_mat.o remd.LES.o \
    rfree.o rgroup.o rmsgrd.o runmd.LES.o \
    sander.LES.o  sander_lib.o set.LES.o sglds.o \
    shake.o short_ene.LES.o softcore.o \
    stack.o state.o thermo_int.o trajene.o wallclock.o \
    xminC.o veclib.o dsarpack.o nxtsec.o sinr_t.o

XRAY_OBJS = \
   xray_globals.o xray_interface.o xray_interface_impl_cpu.o \
   xray_interface_pre_init_data.o

.PHONY: dummy

msander$(SFX): $(BINDIR)/msander$(SFX)

msander.LES$(SFX): $(BINDIR)/msander.LES$(SFX)

msander.cuda$(SFX): $(BINDIR)/msander.cuda$(SFX)

dummy:

#---------------------------------------------------------------------------
$(BINDIR)/msander$(SFX): libs $(MMOBJ) \
         $(NFE_OBJECTS) $(XRAY_OBJS) $(PLUMED_DEPENDENCIES) \
         dummy
	@echo "[MSANDER]  FC $@"
	$(VB)$(FC) $(FFLAGS) -o $@  $(MMOBJ) \
	     $(NFE_OBJECTS) $(XRAY_OBJS) \
	     -L$(LIBDIR) $(NETCDFLIBF) $(FLIBSF) $(LDFLAGS) $(PLUMED_LOAD)

#---------------------------------------------------------------------------
$(BINDIR)/msander.cuda$(SFX): libs lib_cuda $(MMOBJ) \
         $(NFE_OBJECTS) $(XRAY_OBJS) xray_interface_impl_gpu.o \
         $(PLUMED_DEPENDENCIES) dummy
	@echo "[MSANDER]  FC $@"
	$(VB)$(FC) -o $@  $(MMOBJ) \
	     $(NFE_OBJECTS) $(XRAY_OBJS) xray_interface_impl_gpu.o  \
	     -L$(LIBDIR) $(NETCDFLIBF) $(FLIBSF) $(LDFLAGS) $(PLUMED_LOAD) \
         -lxray_cuda -L$(CUDA_HOME)/lib64 -L$(CUDA_HOME)/lib \
         -lcurand -lcufft -lcudart -lcudadevrt -lgfortran -lstdc++

#---------------------------------------------------------------------------
$(BINDIR)/msander.LES$(SFX): libs $(LESOBJ) \
         $(NFE_OBJECTS) $(XRAY_OBJS) \
         $(PLUMED_DEPENDENCIES) dummy
	@echo "[MSANDER_LES]  FC $@"
	$(VB)$(FC) $(FFLAGS) -o $@ $(LESOBJ) \
	     $(XRAY_OBJS) $(NFE_OBJECTS) -L$(LIBDIR) $(NETCDFLIBF) \
	     $(FLIBSF) $(LDFLAGS) $(PLUMED_LOAD)

#---------------------------------------------------------------------------

$(INCDIR)/sander.h: sander.h
	$(VB)/bin/cp sander.h $(INCDIR)

#---------------------------------------------------------------------------
libs:
	$(VB)cd ../blas && $(MAKE) $(BLAS)
	$(VB)cd ../lapack && $(MAKE) $(LAPACK)
	$(VB)cd ../arpack && $(MAKE) $(ARPACK)
	$(VB)$(MAKE) binrestart.o    # needed by rism_io 
	$(VB)cd ../rism && $(MAKE) librism
	$(VB)cd ../xray && $(MAKE) libxray
    
lib_cuda:
	$(VB)cd ../xray/cuda && $(MAKE) libxray_cuda
#---------------------------------------------------------------------------
Plumed.o: Plumed.c
	@echo "[MSANDER]  CC $<"
	$(VB)$(CC) -c $(CNOOPTFLAGS) $(CFLAGS) \
	        -o Plumed.o Plumed.c

rand2.o:  rand2.c
	@echo "[MSANDER]  CC $<"
	$(VB)$(CC) -c $(CNOOPTFLAGS) $(CFLAGS) -o $@ $<

lmodC.o:  lmodC.c
	@echo "[MSANDER]  CC $<"
	$(VB)$(CC) -c -DSQM $(COPTFLAGS) $(CFLAGS) -o $@ $<

xminC.o:  xminC.c
	@echo "[MSANDER]  CC $<"
	$(VB)$(CC) -c -DSQM $(COPTFLAGS) $(CFLAGS) -o $@ $<

#---------------------------------------------------------------------------
clean:
	/bin/rm -f *.o
	/bin/rm -f *.mod
	/bin/rm -f *__genmod.f90  # intel compiler-generated interface modules

#---------------------------------------------------------------------------
depend::
	./makedepend > depend

#----------HEADER DEPENDENCIES (constructed from ./makedepend > depend)
include depend
# DO NOT DELETE
