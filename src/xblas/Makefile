include ../../config.h

OBJS = BLAS_error.o blas_malloc.o BLAS_fpinfo_x.o BLAS_fpinfo_x-f2c.o \
	BLAS_ddot.o BLAS_ddot_x.o  BLAS_ddot-f2c.o BLAS_ddot_x-f2c.o  \
	BLAS_daxpby.o BLAS_daxpby_x.o  BLAS_daxpby-f2c.o BLAS_daxpby_x-f2c.o  \
	BLAS_dwaxpby.o BLAS_dwaxpby_x.o  BLAS_dwaxpby-f2c.o BLAS_dwaxpby_x-f2c.o \
	BLAS_dgemm.o BLAS_dgemm_x.o  BLAS_dgemm-f2c.o BLAS_dgemm_x-f2c.o  \
	BLAS_dgemv.o BLAS_dgemv_x.o BLAS_dgemv-f2c.o BLAS_dgemv_x-f2c.o  \
	BLAS_dgemv2.o BLAS_dgemv2_x.o  BLAS_dgemv2-f2c.o BLAS_dgemv2_x-f2c.o  \

lib-amb: $(LIBDIR)/libxblas-amb.a

$(LIBDIR)/libxblas-amb.a: $(OBJS)
	@echo "[XBLAS] AR $@"
	$(VB)$(AR) $@ $(OBJS) > /dev/null
	$(VB)$(RANLIB) $@

.c.o:
	@echo "[XBLAS] CC $<"
	$(VB)$(CC) $(CFLAGS)  -c -o $@ $<

clean:
	rm -f *.o *~ core


