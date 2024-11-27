#   top-level Makefile for msander

include ../config.h

install: $(INSTALLTYPE)

serial:
	@cd msander && $(MAKE) msander
	@cd msander && $(MAKE) msander.LES
	@cd xray && $(MAKE) xray
	@cd addles && $(MAKE) addles
	@cd leap && $(MAKE) leap
	@cd sff && $(MAKE) sff
	@cd antechamber && $(MAKE) install_ac
	@cd sqm && $(MAKE) install_sqm
	@cd nmr_aux && $(MAKE) nmr_aux
	@cd xtalutil && $(MAKE) xtalutil
	@cd paramfit && $(MAKE) paramfit
	@cd moft && $(MAKE) moft
	@cd saxs && $(MAKE) saxs
	@cd rism && $(MAKE) rism1d
	@cd gbnsr6 && $(MAKE) install_gbnsr6

parallel:
	@cd msander && $(MAKE) msander.MPI
#	@cd msander && $(MAKE) msander.LES.MPI

cuda:
	@cd msander && $(MAKE) msander.cuda

mpicuda:
	@cd msander && $(MAKE) msander.cuda.MPI

clean:
	-(cd blas && $(MAKE) clean )
	-(cd lapack && $(MAKE) clean )
	-(cd arpack && $(MAKE) clean )
	-(cd rism && $(MAKE) clean )
	-(cd msander && $(MAKE) clean )
	-(cd xray && $(MAKE) clean )
	-(cd xray/cuda && $(MAKE) clean )
	-(cd leap && $(MAKE) clean )
	-(cd sff && $(MAKE) clean )
	-(cd antechamber && $(MAKE) clean )
	-(cd nmr_aux && $(MAKE) clean )
	-(cd xtalutil && $(MAKE) clean )
	-(cd paramfit && $(MAKE) clean )
	-(cd moft && $(MAKE) clean )
	-(cd saxs && $(MAKE) clean )

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

