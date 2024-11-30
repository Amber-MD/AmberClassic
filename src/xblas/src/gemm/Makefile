include ../../../../config.h

GEMM_AMB_SRCS =\
	BLAS_dgemm.c\
	BLAS_dgemm_x.c \
	BLAS_dgemm-f2c.c\
	BLAS_dgemm_x-f2c.c 

GEMM_AMB_OBJS = $(GEMM_AMB_SRCS:.c=.o)

lib-amb: $(LIBDIR)/libxblas-amb.a

$(LIBDIR)/libxblas-amb.a: $(GEMM_AMB_OBJS)
	$(VB)$(AR) $@ $(GEMM_AMB_OBJS)
	$(VB)$(RANLIB) $@

.c.o:
	@echo "[XBLAS] CC $<"
	$(VB)$(CC) $(CFLAGS) -I.. -c -o $@ $<

clean:
	rm -f *.o *~ core


