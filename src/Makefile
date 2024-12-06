#   top-level Makefile for AmberClassic

include ../config.h

install: $(INSTALLTYPE)

serial:
	@cd msander && $(MAKE) msander
	@cd msander && $(MAKE) msander.LES
	@cd xray && $(MAKE) xray
	@cd addles && $(MAKE) addles
	@cd leap && $(MAKE) leap
	@cd sff && $(MAKE) libsff
	@cd nabc && $(MAKE) nabc
	@cd antechamber && $(MAKE) install_ac
	@cd sqm && $(MAKE) install_sqm
	@cd nmr_aux && $(MAKE) nmr_aux
	@cd xtalutil && $(MAKE) xtalutil
	@cd paramfit && $(MAKE) paramfit
	@cd moft && $(MAKE) $(MOFT)
	@cd saxs && $(MAKE) saxs
	@cd rism && $(MAKE) rism1d
	@cd gbnsr6 && $(MAKE) install_gbnsr6
	@cd protonate && $(MAKE) install_protonate
	@cd ucpp-1.3 && make ucpp
	@cd nab && make install_nab
#	@cd nss & make install_nss  # compiles, but not tested yet

parallel:
	@cd msander && $(MAKE) msander.MPI
#	@cd msander && $(MAKE) msander.LES.MPI

cuda:
	@cd msander && $(MAKE) msander.cuda

mpicuda:
	@cd msander && $(MAKE) msander.cuda.MPI

conda:
	@cd msander && $(MAKE) msander
	@cd msander && $(MAKE) msander.LES
	@cd xray && $(MAKE) xray
	@cd addles && $(MAKE) addles
	@cd leap && $(MAKE) leap
	@cd sff && $(MAKE) libsff
	@cd antechamber && $(MAKE) install_ac
	@cd sqm && $(MAKE) install_sqm
	@cd nmr_aux && $(MAKE) nmr_aux
	@cd xtalutil && $(MAKE) xtalutil
	@cd paramfit && $(MAKE) paramfit
	@cd moft && $(MAKE) $(MOFT)
	@cd saxs && $(MAKE) saxs
	@cd rism && $(MAKE) rism1d
	@cd gbnsr6 && $(MAKE) install_gbnsr6
	@cd protonate && $(MAKE) install_protonate

clean:
	-(cd blas && $(MAKE) clean )
	-(cd lapack && $(MAKE) clean )
	-(cd arpack && $(MAKE) clean )
	-(cd rism && $(MAKE) clean )
	-(cd msander && $(MAKE) clean )
	-(cd xray && $(MAKE) clean )
	-(cd addles && $(MAKE) clean )
	-(cd xray/cuda && $(MAKE) clean )
	-(cd leap && $(MAKE) clean )
	-(cd sff && $(MAKE) clean )
	-(cd nabc && $(MAKE) clean )
	-(cd antechamber && $(MAKE) clean )
	-(cd sqm && $(MAKE) clean )
	-(cd nmr_aux && $(MAKE) clean )
	-(cd xtalutil && $(MAKE) clean )
	-(cd paramfit && $(MAKE) clean )
	-(cd moft && $(MAKE) clean )
	-(cd saxs && $(MAKE) clean )
	-(cd xblas && $(MAKE) clean )
	-(cd gbnsr6 && $(MAKE) clean )
	-(cd protonate && $(MAKE) clean )
	-(cd ucpp-1.3 && $(MAKE) clean )
	-(cd nab && $(MAKE) clean )
#	-(cd nss && $(MAKE) clean )
	-(cd netcdf-c-4.7.2 && make clean)
	-(cd netcdf-fortran-4.5.2 && make clean)
	-(cd fftw-3.3 && make clean)
	-(cd boost_1_86_0 && make clean)

uninstall: clean
	-(cd $(BINDIR) && /bin/rm -rf *)
	-(cd $(LIBDIR) && /bin/rm -f liblapack.a libblas.a \
        libarpack.a librism.a libsff.a )
	/bin/rm -f $(INCDIR)/* ../config.h

distclean: uninstall
	-(/bin/rm -f netcdf*.log )
	-(cd netcdf-c-4.7.2 && make uninstall)
	-(cd netcdf-fortran-4.5.2 && make uninstall)
	-(cd fftw-3.3 && make uninstall)
	-(cd .. && /bin/rm -rf bin lib inc share logs config.h)

