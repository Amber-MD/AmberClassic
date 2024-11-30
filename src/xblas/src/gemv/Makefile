include ../../../../config.h

GEMV_AMB_SRCS =\
	BLAS_dgemv.c\
	BLAS_dgemv_x.c \
	BLAS_dgemv-f2c.c\
	BLAS_dgemv_x-f2c.c 

GEMV_AMB_OBJS = $(GEMV_AMB_SRCS:.c=.o)

lib-amb: $(LIBDIR)/libxblas-amb.a

$(LIBDIR)/libxblas-amb.a: $(GEMV_AMB_OBJS)
	@echo "[XBLAS] AR $@"
	$(VB)$(AR) $@ $(GEMV_AMB_OBJS)
	$(VB)$(RANLIB) $@

.c.o:
	@echo "[XBLAS] CC $<"
	$(VB)$(CC) $(CFLAGS) -I.. -c -o $@ $<

clean:
	rm -f *.o *~ core


