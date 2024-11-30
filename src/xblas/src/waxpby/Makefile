include ../../../../config.h

WAXPBY_AMB_SRCS =\
	BLAS_dwaxpby.c\
	BLAS_dwaxpby_x.c \
	BLAS_dwaxpby-f2c.c\
	BLAS_dwaxpby_x-f2c.c 

WAXPBY_AMB_OBJS = $(WAXPBY_AMB_SRCS:.c=.o)

lib-amb: $(LIBDIR)/libxblas-amb.a

$(LIBDIR)/libxblas-amb.a: $(WAXPBY_AMB_OBJS)
	@echo "[XBLAS] AR $@"
	$(VB)$(AR) $@ $(WAXPBY_AMB_OBJS)
	$(VB)$(RANLIB) $@

.c.o:
	@echo "[XBLAS] CC $<"
	$(VB)$(CC) $(CFLAGS) -I.. -c -o $@ $<

clean:
	rm -f *.o *~ core


