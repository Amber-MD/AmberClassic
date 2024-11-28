include ../../make.conf
include ../../$(MAKEINC)

LIB_PATH = ../../$(OUTPUT_DIR)
SRC_PATH = .
HEADER_PATH = ..

HEMV2_SRCS =\
	$(SRC_PATH)/BLAS_chemv2.c\
	$(SRC_PATH)/BLAS_zhemv2.c\
	$(SRC_PATH)/BLAS_zhemv2_z_c.c\
	$(SRC_PATH)/BLAS_zhemv2_c_z.c\
	$(SRC_PATH)/BLAS_zhemv2_c_c.c\
	$(SRC_PATH)/BLAS_chemv2_c_s.c\
	$(SRC_PATH)/BLAS_zhemv2_z_d.c\
	$(SRC_PATH)/BLAS_chemv2_x.c\
	$(SRC_PATH)/BLAS_zhemv2_x.c\
	$(SRC_PATH)/BLAS_zhemv2_z_c_x.c\
	$(SRC_PATH)/BLAS_zhemv2_c_z_x.c\
	$(SRC_PATH)/BLAS_zhemv2_c_c_x.c\
	$(SRC_PATH)/BLAS_chemv2_c_s_x.c\
	$(SRC_PATH)/BLAS_zhemv2_z_d_x.c \
	\
	$(SRC_PATH)/BLAS_chemv2-f2c.c\
	$(SRC_PATH)/BLAS_zhemv2-f2c.c\
	$(SRC_PATH)/BLAS_zhemv2_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zhemv2_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zhemv2_c_c-f2c.c\
	$(SRC_PATH)/BLAS_chemv2_c_s-f2c.c\
	$(SRC_PATH)/BLAS_zhemv2_z_d-f2c.c\
	$(SRC_PATH)/BLAS_chemv2_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemv2_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemv2_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemv2_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemv2_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_chemv2_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemv2_z_d_x-f2c.c 

HEMV2_AMB_SRCS =\
	$(SRC_PATH)/BLAS_HEMV2_AMB_ARGShemv2__.c \
	\
	$(SRC_PATH)/BLAS_HEMV2_AMB_ARGShemv2__-f2c.c 

HEMV2_OBJS = $(HEMV2_SRCS:.c=.o)

HEMV2_AMB_OBJS = $(HEMV2_AMB_SRCS:.c=.o)

all: $(HEMV2_OBJS)

all-amb: $(HEMV2_AMB_OBJS)

lib: $(HEMV2_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(HEMV2_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

lib-amb: $(HEMV2_AMB_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB_AMB) $(HEMV2_AMB_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB_AMB)

.c.o: 
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK

