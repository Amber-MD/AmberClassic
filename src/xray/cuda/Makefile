include ../../../config.h

.SUFFIXES: .cu .o

#CUDA Specific build flags  -- use nvcc in the PATH for now; 
#                              may want to set with configure later
NVCC=nvcc \
   -Wno-deprecated-declarations -use_fast_math -O3 \
   -I/usr/local/cuda/include -IB40C $(CUDAFLAGS) $(NVCC_CFLAGS) --std c++14

CUOBJS = BulkMask.o BulkMaskCPU.o BulkMaskGPU.o \
         DPartial.o DPartialCPU.o DPartialGPU.o \
         NonBulk.o NonBulkCPU.o NonBulkGPU.o \
         UnitCell.o xray_bulk_mask.o \
         xray_dpartial.o xray_non_bulk.o

.cu.o:
	@echo "[XRAY] nvcc $<"
	$(VB)$(NVCC) -c -I$(INCDIR) -o $@ $<

xray_cuda.o: $(CUOBJS)
	@echo "[XRAY] nvcc -dlink $@"
	$(VB)$(NVCC) -dlink -o $@ $(CUOBJS)

libxray_cuda: $(LIBDIR)/libxray_cuda.a

$(LIBDIR)/libxray_cuda.a: xray_cuda.o
	@echo "[XRAY]  AR $@"
	$(VB)ar rvs $@ $< $(CUOBJS) $(CPPOBJS) > /dev/null

clean:
	/bin/rm -f *.o
