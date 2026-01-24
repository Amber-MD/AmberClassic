include ../../config.h

LOCALFLAGS = $(CUDAFFLAGS) 

SRC = add_xray.F90 add_pdb.F90
OBJ = add_xray.o add_pdb.o

LIBOBJ = xray_atomic_scatter_factor.o xray_atomic_scatter_factor_impl_cpu.o \
xray_atomic_scatter_factor_impl_gpu.o xray_bulk_mask_data.o \
xray_bulk_mask.o xray_bulk_mask_impl_cpu.o xray_bulk_mask_impl_gpu.o \
xray_bulk_model_afonine_2013.o xray_bulk_model.o xray_bulk_model_none.o \
xray_bulk_model_simple.o xray_contracts.o \
xray_dpartial_data.o xray_dpartial.o xray_dpartial_impl_cpu.o \
xray_dpartial_impl_gpu.o xray_interface2_data.o xray_interface2.o \
xray_non_bulk_data.o xray_non_bulk.o xray_non_bulk_impl_cpu.o \
xray_non_bulk_impl_gpu.o xray_pure_utils.o xray_scaling_data.o \
xray_scaling.o xray_scaling_impl_cpu.o xray_scaling_impl_gpu.o \
xray_target.o xray_target_least_squares.o xray_target_max_likelihood_data.o \
xray_target_max_likelihood.o xray_target_max_likelihood_impl_cpu.o \
xray_target_vector_least_squares_data.o xray_target_vector_least_squares.o \
xray_target_vector_least_squares_impl_cpu.o \
xray_unit_cell.o constants_xray.o \
xray_fft_interface_impl_fftw.o xray_fft_interface_impl_mkl.o \
xray_fft_interface_impl_none.o xray_fft_interface.o


libxray: $(LIBDIR)/libxray.a

$(LIBDIR)/libxray.a: $(LIBOBJ)
	@echo "[XRAY]  AR $@"
	$(VB)$(AR) $@ $(LIBOBJ) > /dev/null
	$(VB)$(RANLIB) $@

$(BINDIR)/add_xray: add_xray.o nxtsec.o
	@echo "[XRAY]  FC $@"
	$(VB)$(FC) $(FFLAGS) -o $@ $(LDFLAGS) add_xray.o nxtsec.o

$(BINDIR)/add_pdb: add_pdb.o nxtsec.o
	@echo "[XRAY]  FC $@"
	$(VB)$(FC) $(FFLAGS) -o $@ $(LDFLAGS) add_pdb.o nxtsec.o

xray: $(BINDIR)/add_xray $(BINDIR)/add_pdb libxray

clean:
	/bin/rm -f *.o *.mod

nxtsec.o:  ../msander/nxtsec.F90
	@echo "[XRAY]  FC $<"
	$(VB)$(FC) -c $(FNOOPTFLAGS) $(FFLAGS) -o $@ $<

depend::
	./makedepend > depend

uninstall: clean
	-rm -f $(BINDIR)/add_xray
	-rm -f $(BINDIR)/add_pdb

include depend
