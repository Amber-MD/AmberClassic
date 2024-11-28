include ../../make.conf
include ../../$(MAKEINC)

LIB_PATH = ../../$(OUTPUT_DIR)
SRC_PATH = .
HEADER_PATH = ..

HEMM_SRCS =\
	$(SRC_PATH)/BLAS_chemm.c\
	$(SRC_PATH)/BLAS_zhemm.c\
	$(SRC_PATH)/BLAS_zhemm_z_c.c\
	$(SRC_PATH)/BLAS_zhemm_c_z.c\
	$(SRC_PATH)/BLAS_zhemm_c_c.c\
	$(SRC_PATH)/BLAS_chemm_c_s.c\
	$(SRC_PATH)/BLAS_zhemm_z_d.c\
	$(SRC_PATH)/BLAS_chemm_x.c\
	$(SRC_PATH)/BLAS_zhemm_x.c\
	$(SRC_PATH)/BLAS_zhemm_z_c_x.c\
	$(SRC_PATH)/BLAS_zhemm_c_z_x.c\
	$(SRC_PATH)/BLAS_zhemm_c_c_x.c\
	$(SRC_PATH)/BLAS_chemm_c_s_x.c\
	$(SRC_PATH)/BLAS_zhemm_z_d_x.c \
	\
	$(SRC_PATH)/BLAS_chemm-f2c.c\
	$(SRC_PATH)/BLAS_zhemm-f2c.c\
	$(SRC_PATH)/BLAS_zhemm_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zhemm_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zhemm_c_c-f2c.c\
	$(SRC_PATH)/BLAS_chemm_c_s-f2c.c\
	$(SRC_PATH)/BLAS_zhemm_z_d-f2c.c\
	$(SRC_PATH)/BLAS_chemm_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemm_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemm_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemm_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemm_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_chemm_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zhemm_z_d_x-f2c.c 

HEMM_AMB_SRCS =\
	$(SRC_PATH)/BLAS_HEMM_AMB_ARGShemm__.c \
	\
	$(SRC_PATH)/BLAS_HEMM_AMB_ARGShemm__-f2c.c 

HEMM_OBJS = $(HEMM_SRCS:.c=.o)

HEMM_AMB_OBJS = $(HEMM_AMB_SRCS:.c=.o)

all: $(HEMM_OBJS)

all-amb: $(HEMM_AMB_OBJS)

lib: $(HEMM_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(HEMM_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

lib-amb: $(HEMM_AMB_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB_AMB) $(HEMM_AMB_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB_AMB)

.c.o: 
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK

