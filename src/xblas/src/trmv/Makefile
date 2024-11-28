include ../../make.conf
include ../../$(MAKEINC)

LIB_PATH = ../../$(OUTPUT_DIR)
SRC_PATH = .
HEADER_PATH = ..

TRMV_SRCS =\
	$(SRC_PATH)/BLAS_strmv.c\
	$(SRC_PATH)/BLAS_dtrmv.c\
	$(SRC_PATH)/BLAS_ctrmv.c\
	$(SRC_PATH)/BLAS_ztrmv.c\
	$(SRC_PATH)/BLAS_dtrmv_s.c\
	$(SRC_PATH)/BLAS_ztrmv_c.c\
	$(SRC_PATH)/BLAS_ctrmv_s.c\
	$(SRC_PATH)/BLAS_ztrmv_d.c\
	$(SRC_PATH)/BLAS_strmv_x.c\
	$(SRC_PATH)/BLAS_dtrmv_x.c\
	$(SRC_PATH)/BLAS_ctrmv_x.c\
	$(SRC_PATH)/BLAS_ztrmv_x.c\
	$(SRC_PATH)/BLAS_dtrmv_s_x.c\
	$(SRC_PATH)/BLAS_ztrmv_c_x.c\
	$(SRC_PATH)/BLAS_ctrmv_s_x.c\
	$(SRC_PATH)/BLAS_ztrmv_d_x.c \
	\
	$(SRC_PATH)/BLAS_strmv-f2c.c\
	$(SRC_PATH)/BLAS_dtrmv-f2c.c\
	$(SRC_PATH)/BLAS_ctrmv-f2c.c\
	$(SRC_PATH)/BLAS_ztrmv-f2c.c\
	$(SRC_PATH)/BLAS_dtrmv_s-f2c.c\
	$(SRC_PATH)/BLAS_ztrmv_c-f2c.c\
	$(SRC_PATH)/BLAS_ctrmv_s-f2c.c\
	$(SRC_PATH)/BLAS_ztrmv_d-f2c.c\
	$(SRC_PATH)/BLAS_strmv_x-f2c.c\
	$(SRC_PATH)/BLAS_dtrmv_x-f2c.c\
	$(SRC_PATH)/BLAS_ctrmv_x-f2c.c\
	$(SRC_PATH)/BLAS_ztrmv_x-f2c.c\
	$(SRC_PATH)/BLAS_dtrmv_s_x-f2c.c\
	$(SRC_PATH)/BLAS_ztrmv_c_x-f2c.c\
	$(SRC_PATH)/BLAS_ctrmv_s_x-f2c.c\
	$(SRC_PATH)/BLAS_ztrmv_d_x-f2c.c 

TRMV_AMB_SRCS =\
	$(SRC_PATH)/BLAS_TRMV_AMB_ARGStrmv_.c \
	\
	$(SRC_PATH)/BLAS_TRMV_AMB_ARGStrmv_-f2c.c 

TRMV_OBJS = $(TRMV_SRCS:.c=.o)

TRMV_AMB_OBJS = $(TRMV_AMB_SRCS:.c=.o)

all: $(TRMV_OBJS)

all-amb: $(TRMV_AMB_OBJS)

lib: $(TRMV_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(TRMV_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

lib-amb: $(TRMV_AMB_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB_AMB) $(TRMV_AMB_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB_AMB)

.c.o: 
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK

