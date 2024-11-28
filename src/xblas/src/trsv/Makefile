include ../../make.conf
include ../../$(MAKEINC)

LIB_PATH = ../../$(OUTPUT_DIR)
SRC_PATH = .
HEADER_PATH = ..

TRSV_SRCS =\
	$(SRC_PATH)/BLAS_strsv.c\
	$(SRC_PATH)/BLAS_dtrsv.c\
	$(SRC_PATH)/BLAS_ctrsv.c\
	$(SRC_PATH)/BLAS_ztrsv.c\
	$(SRC_PATH)/BLAS_dtrsv_s.c\
	$(SRC_PATH)/BLAS_ztrsv_c.c\
	$(SRC_PATH)/BLAS_ctrsv_s.c\
	$(SRC_PATH)/BLAS_ztrsv_d.c\
	$(SRC_PATH)/BLAS_strsv_x.c\
	$(SRC_PATH)/BLAS_dtrsv_x.c\
	$(SRC_PATH)/BLAS_dtrsv_s_x.c\
	$(SRC_PATH)/BLAS_ctrsv_x.c\
	$(SRC_PATH)/BLAS_ztrsv_x.c\
	$(SRC_PATH)/BLAS_ztrsv_c_x.c\
	$(SRC_PATH)/BLAS_ctrsv_s_x.c\
	$(SRC_PATH)/BLAS_ztrsv_d_x.c \
	\
	$(SRC_PATH)/BLAS_strsv-f2c.c\
	$(SRC_PATH)/BLAS_dtrsv-f2c.c\
	$(SRC_PATH)/BLAS_ctrsv-f2c.c\
	$(SRC_PATH)/BLAS_ztrsv-f2c.c\
	$(SRC_PATH)/BLAS_dtrsv_s-f2c.c\
	$(SRC_PATH)/BLAS_ztrsv_c-f2c.c\
	$(SRC_PATH)/BLAS_ctrsv_s-f2c.c\
	$(SRC_PATH)/BLAS_ztrsv_d-f2c.c\
	$(SRC_PATH)/BLAS_strsv_x-f2c.c\
	$(SRC_PATH)/BLAS_dtrsv_x-f2c.c\
	$(SRC_PATH)/BLAS_dtrsv_s_x-f2c.c\
	$(SRC_PATH)/BLAS_ctrsv_x-f2c.c\
	$(SRC_PATH)/BLAS_ztrsv_x-f2c.c\
	$(SRC_PATH)/BLAS_ztrsv_c_x-f2c.c\
	$(SRC_PATH)/BLAS_ctrsv_s_x-f2c.c\
	$(SRC_PATH)/BLAS_ztrsv_d_x-f2c.c 

TRSV_AMB_SRCS =\
	$(SRC_PATH)/BLAS_TRSV_AMB_ARGStrsv_.c \
	\
	$(SRC_PATH)/BLAS_TRSV_AMB_ARGStrsv_-f2c.c 

TRSV_OBJS = $(TRSV_SRCS:.c=.o)

TRSV_AMB_OBJS = $(TRSV_AMB_SRCS:.c=.o)

all: $(TRSV_OBJS)

all-amb: $(TRSV_AMB_OBJS)

lib: $(TRSV_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(TRSV_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

lib-amb: $(TRSV_AMB_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB_AMB) $(TRSV_AMB_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB_AMB)

.c.o: 
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK

