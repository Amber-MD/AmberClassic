include ../../make.conf
include ../../$(MAKEINC)

LIB_PATH = ../../$(OUTPUT_DIR)
SRC_PATH = .
HEADER_PATH = ..

GEMV2_SRCS =\
	$(SRC_PATH)/BLAS_sgemv2.c\
	$(SRC_PATH)/BLAS_dgemv2.c\
	$(SRC_PATH)/BLAS_cgemv2.c\
	$(SRC_PATH)/BLAS_zgemv2.c\
	$(SRC_PATH)/BLAS_dgemv2_d_s.c\
	$(SRC_PATH)/BLAS_dgemv2_s_d.c\
	$(SRC_PATH)/BLAS_dgemv2_s_s.c\
	$(SRC_PATH)/BLAS_zgemv2_z_c.c\
	$(SRC_PATH)/BLAS_zgemv2_c_z.c\
	$(SRC_PATH)/BLAS_zgemv2_c_c.c\
	$(SRC_PATH)/BLAS_cgemv2_c_s.c\
	$(SRC_PATH)/BLAS_cgemv2_s_c.c\
	$(SRC_PATH)/BLAS_cgemv2_s_s.c\
	$(SRC_PATH)/BLAS_zgemv2_z_d.c\
	$(SRC_PATH)/BLAS_zgemv2_d_z.c\
	$(SRC_PATH)/BLAS_zgemv2_d_d.c\
	$(SRC_PATH)/BLAS_sgemv2_x.c\
	$(SRC_PATH)/BLAS_dgemv2_x.c\
	$(SRC_PATH)/BLAS_cgemv2_x.c\
	$(SRC_PATH)/BLAS_zgemv2_x.c\
	$(SRC_PATH)/BLAS_dgemv2_d_s_x.c\
	$(SRC_PATH)/BLAS_dgemv2_s_d_x.c\
	$(SRC_PATH)/BLAS_dgemv2_s_s_x.c\
	$(SRC_PATH)/BLAS_zgemv2_z_c_x.c\
	$(SRC_PATH)/BLAS_zgemv2_c_z_x.c\
	$(SRC_PATH)/BLAS_zgemv2_c_c_x.c\
	$(SRC_PATH)/BLAS_cgemv2_c_s_x.c\
	$(SRC_PATH)/BLAS_cgemv2_s_c_x.c\
	$(SRC_PATH)/BLAS_cgemv2_s_s_x.c\
	$(SRC_PATH)/BLAS_zgemv2_z_d_x.c\
	$(SRC_PATH)/BLAS_zgemv2_d_z_x.c\
	$(SRC_PATH)/BLAS_zgemv2_d_d_x.c \
	\
	$(SRC_PATH)/BLAS_sgemv2-f2c.c\
	$(SRC_PATH)/BLAS_dgemv2-f2c.c\
	$(SRC_PATH)/BLAS_cgemv2-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2-f2c.c\
	$(SRC_PATH)/BLAS_dgemv2_d_s-f2c.c\
	$(SRC_PATH)/BLAS_dgemv2_s_d-f2c.c\
	$(SRC_PATH)/BLAS_dgemv2_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2_c_c-f2c.c\
	$(SRC_PATH)/BLAS_cgemv2_c_s-f2c.c\
	$(SRC_PATH)/BLAS_cgemv2_s_c-f2c.c\
	$(SRC_PATH)/BLAS_cgemv2_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2_z_d-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2_d_z-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2_d_d-f2c.c\
	$(SRC_PATH)/BLAS_sgemv2_x-f2c.c\
	$(SRC_PATH)/BLAS_dgemv2_x-f2c.c\
	$(SRC_PATH)/BLAS_cgemv2_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2_x-f2c.c\
	$(SRC_PATH)/BLAS_dgemv2_d_s_x-f2c.c\
	$(SRC_PATH)/BLAS_dgemv2_s_d_x-f2c.c\
	$(SRC_PATH)/BLAS_dgemv2_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cgemv2_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_cgemv2_s_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cgemv2_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2_z_d_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2_d_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zgemv2_d_d_x-f2c.c 

GEMV2_AMB_SRCS =\
	$(SRC_PATH)/BLAS_dgemv2.c\
	$(SRC_PATH)/BLAS_dgemv2_x.c \
	\
	$(SRC_PATH)/BLAS_dgemv2-f2c.c\
	$(SRC_PATH)/BLAS_dgemv2_x-f2c.c 

GEMV2_OBJS = $(GEMV2_SRCS:.c=.o)

GEMV2_AMB_OBJS = $(GEMV2_AMB_SRCS:.c=.o)

all: $(GEMV2_OBJS)

all-amb: $(GEMV2_AMB_OBJS)

lib: $(GEMV2_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(GEMV2_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

lib-amb: $(GEMV2_AMB_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB_AMB) $(GEMV2_AMB_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB_AMB)

.c.o: 
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK

