include ../../../../config.h

COMM_OBJS = BLAS_error.o blas_malloc.o BLAS_fpinfo_x.o BLAS_fpinfo_x-f2c.o

lib-amb: $(LIBDIR)/libxblas-amb.a

$(LIBDIR)/libxblas-amb.a: $(COMM_OBJS)
	@echo "[XBLAS] AR $@"
	$(VB)$(AR) $@ $(COMM_OBJS)
	$(VB)$(RANLIB) $@

.c.o:
	@echo "[XBLAS] CC $<"
	$(VB)$(CC) $(CFLAGS) -I.. -c -o $@ $<

clean:
	rm -f *.o *~ core


