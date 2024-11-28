include ../../make.conf
include ../../$(MAKEINC)

LIB_PATH = ../../$(OUTPUT_DIR)
SRC_PATH = .
HEADER_PATH = ..

HEMV_SRCS =\
	$(SRC_PATH)/BLAS_chemv.c\
	$(SRC_PATH)/BLAS_zhemv.c\
	$(SRC_PATH)/BLAS_zhemv_z_c.c\
	$(SRC_PATH)/BLAS_zhemv_c_z.c\
	$(SRC_PATH)/BLAS_zhemv_c_c.c\
	$(SRC_PATH)/BLAS_chemv_c_s.c\
	$(SRC_PATH)/BLAS_zhemv_z_d.c\
	$(SRC_PATH)/BLAS_chemv_x.c\
	$(SRC_PATH)/BLAS_zhemv_x.c\
	$(SRC_PATH)/BLAS_zhemv_z_c_x.c\
	$(SRC_PATH)/BLAS_zhemv_c_z_x.c\
	$(SRC_PATH)/BLAS_zhemv_c_c_x.c\
	$(SRC_PATH)/BLAS_chemv_c_s_x.c\
	$(SRC_PATH)/BLAS_zhemv_z_d_x.c \
	\
	$(SRC_PATH)/BLAS_chemv-f2c.c\
	$(SRC_PATH)/BLAS_zhemv-f2c.c\
	$(SRC_PATH)/BLAS_zhemv_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zhemv_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zhemv_c_c-f2c.c\
	$(SRC_PATH)/BLAS_chemv_c_s-f2c.c\
	$(SRC_PATH)/BLAS_zhemv_z_d-f2c.c\
	$(SRC_PATH)/BLAS_chemv_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemv_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemv_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemv_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemv_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_chemv_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemv_z_d_x-f2c.c 

HEMV_AMB_SRCS =\
	$(SRC_PATH)/BLAS_HEMV_AMB_ARGShemv__.c \
	\
	$(SRC_PATH)/BLAS_HEMV_AMB_ARGShemv__-f2c.c 

HEMV_OBJS = $(HEMV_SRCS:.c=.o)

HEMV_AMB_OBJS = $(HEMV_AMB_SRCS:.c=.o)

all: $(HEMV_OBJS)

all-amb: $(HEMV_AMB_OBJS)

lib: $(HEMV_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(HEMV_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

lib-amb: $(HEMV_AMB_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB_AMB) $(HEMV_AMB_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB_AMB)

.c.o: 
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK

