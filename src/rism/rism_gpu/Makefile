SHELL=/bin/sh
include ../../../config.h

.SUFFIXES: .cpp .cu .o

#CUDA Specific build flags  -- use nvcc in the PATH for now; 
#                              may want to set with configure later
NVCC=nvcc \
   -Wno-deprecated-declarations -use_fast_math -O3 \
   -I/usr/local/cuda/include -IB40C $(CUDAFLAGS) $(NVCC_CFLAGS) --std c++14

.cu.o:
	@echo "[RISM] nvcc $<"
	$(VB)$(NVCC) -c -o $@ $<

.cpp.o:
	@echo "[RISM] nvcc $<"
	$(VB)$(NVCC) -c -o $@ $<

LIBOBJ_CUDA = mdiis.o k_rism3d.o \
    rism3d_grid.o rism3d_closure.o rism3d_closure_kh.o \
    k_kh.o rism3d_closure_hnc.o rism3d_closure_psen.o \
    k_grid.o rism3d_fft.o rism_util.o k_rism_util.o \
    rism3d.o rism3d_solute_class.o rism3d_solvent_class.o \
    rism3d_potential.o k_potential.o 

librism_cuda: $(LIBDIR)/librism_cuda.a

$(LIBDIR)/librism_cuda.a: $(LIBOBJ_CUDA)
	@echo "[RISM]  AR $@"
	$(VB)$(AR) $@ $(LIBOBJ_CUDA) > /dev/null
	$(VB)$(RANLIB) $@

uninstall:
	/bin/rm -f $(LIBDIR)/librism_cuda.a

clean:
	-/bin/rm -f *.o

