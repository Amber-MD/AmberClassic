include ../../../../config.h

GEMV2_AMB_SRCS =\
	BLAS_dgemv2.c\
	BLAS_dgemv2_x.c \
	BLAS_dgemv2-f2c.c\
	BLAS_dgemv2_x-f2c.c 

GEMV2_AMB_OBJS = $(GEMV2_AMB_SRCS:.c=.o)

lib-amb: $(LIBDIR)/libxblas-amb.a

$(LIBDIR)/libxblas-amb.a: $(GEMV2_AMB_OBJS)
	@echo "[XBLAS] AR $@"
	$(VB)$(AR) $@ $(GEMV2_AMB_OBJS)
	$(VB)$(RANLIB) $@

.c.o:
	@echo "[XBLAS] CC $<"
	$(VB)$(CC) $(CFLAGS) -I.. -c -o $@ $<

clean:
	rm -f *.o *~ core


