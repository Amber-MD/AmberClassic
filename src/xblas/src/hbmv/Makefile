include ../../make.conf
include ../../$(MAKEINC)

LIB_PATH = ../../$(OUTPUT_DIR)
SRC_PATH = .
HEADER_PATH = ..

HBMV_SRCS =\
	$(SRC_PATH)/BLAS_chbmv.c\
	$(SRC_PATH)/BLAS_zhbmv.c\
	$(SRC_PATH)/BLAS_zhbmv_z_c.c\
	$(SRC_PATH)/BLAS_zhbmv_c_z.c\
	$(SRC_PATH)/BLAS_zhbmv_c_c.c\
	$(SRC_PATH)/BLAS_chbmv_c_s.c\
	$(SRC_PATH)/BLAS_zhbmv_z_d.c\
	$(SRC_PATH)/BLAS_chbmv_x.c\
	$(SRC_PATH)/BLAS_zhbmv_x.c\
	$(SRC_PATH)/BLAS_zhbmv_z_c_x.c\
	$(SRC_PATH)/BLAS_zhbmv_c_z_x.c\
	$(SRC_PATH)/BLAS_zhbmv_c_c_x.c\
	$(SRC_PATH)/BLAS_chbmv_c_s_x.c\
	$(SRC_PATH)/BLAS_zhbmv_z_d_x.c \
	\
	$(SRC_PATH)/BLAS_chbmv-f2c.c\
	$(SRC_PATH)/BLAS_zhbmv-f2c.c\
	$(SRC_PATH)/BLAS_zhbmv_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zhbmv_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zhbmv_c_c-f2c.c\
	$(SRC_PATH)/BLAS_chbmv_c_s-f2c.c\
	$(SRC_PATH)/BLAS_zhbmv_z_d-f2c.c\
	$(SRC_PATH)/BLAS_chbmv_x-f2c.c\
	$(SRC_PATH)/BLAS_zhbmv_x-f2c.c\
	$(SRC_PATH)/BLAS_zhbmv_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zhbmv_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zhbmv_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_chbmv_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zhbmv_z_d_x-f2c.c 

HBMV_AMB_SRCS =\
	$(SRC_PATH)/BLAS_HBMV_AMB_ARGShbmv__.c \
	\
	$(SRC_PATH)/BLAS_HBMV_AMB_ARGShbmv__-f2c.c 

HBMV_OBJS = $(HBMV_SRCS:.c=.o)

HBMV_AMB_OBJS = $(HBMV_AMB_SRCS:.c=.o)

all: $(HBMV_OBJS)

all-amb: $(HBMV_AMB_OBJS)

lib: $(HBMV_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(HBMV_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

lib-amb: $(HBMV_AMB_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB_AMB) $(HBMV_AMB_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB_AMB)

.c.o: 
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK

