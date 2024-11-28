include ../../make.conf
include ../../$(MAKEINC)

LIB_PATH = ../../$(OUTPUT_DIR)
SRC_PATH = .
HEADER_PATH = ..

GEMV_SRCS =\
	$(SRC_PATH)/BLAS_sgemv.c\
	$(SRC_PATH)/BLAS_dgemv.c\
	$(SRC_PATH)/BLAS_cgemv.c\
	$(SRC_PATH)/BLAS_zgemv.c\
	$(SRC_PATH)/BLAS_dgemv_d_s.c\
	$(SRC_PATH)/BLAS_dgemv_s_d.c\
	$(SRC_PATH)/BLAS_dgemv_s_s.c\
	$(SRC_PATH)/BLAS_zgemv_z_c.c\
	$(SRC_PATH)/BLAS_zgemv_c_z.c\
	$(SRC_PATH)/BLAS_zgemv_c_c.c\
	$(SRC_PATH)/BLAS_cgemv_c_s.c\
	$(SRC_PATH)/BLAS_cgemv_s_c.c\
	$(SRC_PATH)/BLAS_cgemv_s_s.c\
	$(SRC_PATH)/BLAS_zgemv_z_d.c\
	$(SRC_PATH)/BLAS_zgemv_d_z.c\
	$(SRC_PATH)/BLAS_zgemv_d_d.c\
	$(SRC_PATH)/BLAS_sgemv_x.c\
	$(SRC_PATH)/BLAS_dgemv_x.c\
	$(SRC_PATH)/BLAS_cgemv_x.c\
	$(SRC_PATH)/BLAS_zgemv_x.c\
	$(SRC_PATH)/BLAS_dgemv_d_s_x.c\
	$(SRC_PATH)/BLAS_dgemv_s_d_x.c\
	$(SRC_PATH)/BLAS_dgemv_s_s_x.c\
	$(SRC_PATH)/BLAS_zgemv_z_c_x.c\
	$(SRC_PATH)/BLAS_zgemv_c_z_x.c\
	$(SRC_PATH)/BLAS_zgemv_c_c_x.c\
	$(SRC_PATH)/BLAS_cgemv_c_s_x.c\
	$(SRC_PATH)/BLAS_cgemv_s_c_x.c\
	$(SRC_PATH)/BLAS_cgemv_s_s_x.c\
	$(SRC_PATH)/BLAS_zgemv_z_d_x.c\
	$(SRC_PATH)/BLAS_zgemv_d_z_x.c\
	$(SRC_PATH)/BLAS_zgemv_d_d_x.c \
	\
	$(SRC_PATH)/BLAS_sgemv-f2c.c\
	$(SRC_PATH)/BLAS_dgemv-f2c.c\
	$(SRC_PATH)/BLAS_cgemv-f2c.c\
	$(SRC_PATH)/BLAS_zgemv-f2c.c\
	$(SRC_PATH)/BLAS_dgemv_d_s-f2c.c\
	$(SRC_PATH)/BLAS_dgemv_s_d-f2c.c\
	$(SRC_PATH)/BLAS_dgemv_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zgemv_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zgemv_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zgemv_c_c-f2c.c\
	$(SRC_PATH)/BLAS_cgemv_c_s-f2c.c\
	$(SRC_PATH)/BLAS_cgemv_s_c-f2c.c\
	$(SRC_PATH)/BLAS_cgemv_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zgemv_z_d-f2c.c\
	$(SRC_PATH)/BLAS_zgemv_d_z-f2c.c\
	$(SRC_PATH)/BLAS_zgemv_d_d-f2c.c\
	$(SRC_PATH)/BLAS_sgemv_x-f2c.c\
	$(SRC_PATH)/BLAS_dgemv_x-f2c.c\
	$(SRC_PATH)/BLAS_cgemv_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv_x-f2c.c\
	$(SRC_PATH)/BLAS_dgemv_d_s_x-f2c.c\
	$(SRC_PATH)/BLAS_dgemv_s_d_x-f2c.c\
	$(SRC_PATH)/BLAS_dgemv_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cgemv_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_cgemv_s_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cgemv_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv_z_d_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv_d_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv_d_d_x-f2c.c 

GEMV_AMB_SRCS =\
	$(SRC_PATH)/BLAS_dgemv.c\
	$(SRC_PATH)/BLAS_dgemv_x.c \
	\
	$(SRC_PATH)/BLAS_dgemv-f2c.c\
	$(SRC_PATH)/BLAS_dgemv_x-f2c.c 

GEMV_OBJS = $(GEMV_SRCS:.c=.o)

GEMV_AMB_OBJS = $(GEMV_AMB_SRCS:.c=.o)

all: $(GEMV_OBJS)

all-amb: $(GEMV_AMB_OBJS)

lib: $(GEMV_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(GEMV_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

lib-amb: $(GEMV_AMB_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB_AMB) $(GEMV_AMB_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB_AMB)

.c.o: 
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK

