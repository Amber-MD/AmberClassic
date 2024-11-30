include ../../../../config.h

AXPBY_AMB_SRCS =\
	BLAS_daxpby.c\
	BLAS_daxpby_x.c \
	BLAS_daxpby-f2c.c\
	BLAS_daxpby_x-f2c.c 

AXPBY_AMB_OBJS = $(AXPBY_AMB_SRCS:.c=.o)

lib-amb: $(LIBDIR)/libxblas-amb.a

$(LIBDIR)/libxblas-amb.a: $(AXPBY_AMB_OBJS)
	$(VB)$(AR) $@ $(AXPBY_AMB_OBJS)
	$(VB)$(RANLIB) $@

.c.o:
	@echo "[XBLAS] CC $<"
	$(VB)$(CC) $(CFLAGS) -I.. -c -o $@ $<

clean:
	rm -f *.o *~ core


