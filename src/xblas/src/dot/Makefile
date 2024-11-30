include ../../../../config.h

DOT_AMB_SRCS =\
	BLAS_ddot.c\
	BLAS_ddot_x.c \
	BLAS_ddot-f2c.c\
	BLAS_ddot_x-f2c.c 

DOT_AMB_OBJS = $(DOT_AMB_SRCS:.c=.o)

lib-amb: $(LIBDIR)/libxblas-amb.a

$(LIBDIR)/libxblas-amb.a: $(DOT_AMB_OBJS)
	@echo "[XBLAS] AR $@"
	$(VB)$(AR) $@ $(DOT_AMB_OBJS)
	$(VB)$(RANLIB) $@

.c.o:
	@echo "[XBLAS] CC $<"
	$(VB)$(CC) $(CFLAGS) -I.. -c -o $@ $<

clean:
	rm -f *.o *~ core


