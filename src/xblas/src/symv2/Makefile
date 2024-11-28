include ../../make.conf
include ../../$(MAKEINC)

LIB_PATH = ../../$(OUTPUT_DIR)
SRC_PATH = .
HEADER_PATH = ..

SYMV2_SRCS =\
	$(SRC_PATH)/BLAS_ssymv2.c\
	$(SRC_PATH)/BLAS_dsymv2.c\
	$(SRC_PATH)/BLAS_csymv2.c\
	$(SRC_PATH)/BLAS_zsymv2.c\
	$(SRC_PATH)/BLAS_dsymv2_d_s.c\
	$(SRC_PATH)/BLAS_dsymv2_s_d.c\
	$(SRC_PATH)/BLAS_dsymv2_s_s.c\
	$(SRC_PATH)/BLAS_zsymv2_z_c.c\
	$(SRC_PATH)/BLAS_zsymv2_c_z.c\
	$(SRC_PATH)/BLAS_zsymv2_c_c.c\
	$(SRC_PATH)/BLAS_csymv2_c_s.c\
	$(SRC_PATH)/BLAS_csymv2_s_c.c\
	$(SRC_PATH)/BLAS_csymv2_s_s.c\
	$(SRC_PATH)/BLAS_zsymv2_z_d.c\
	$(SRC_PATH)/BLAS_zsymv2_d_z.c\
	$(SRC_PATH)/BLAS_zsymv2_d_d.c\
	$(SRC_PATH)/BLAS_ssymv2_x.c\
	$(SRC_PATH)/BLAS_dsymv2_x.c\
	$(SRC_PATH)/BLAS_csymv2_x.c\
	$(SRC_PATH)/BLAS_zsymv2_x.c\
	$(SRC_PATH)/BLAS_dsymv2_d_s_x.c\
	$(SRC_PATH)/BLAS_dsymv2_s_d_x.c\
	$(SRC_PATH)/BLAS_dsymv2_s_s_x.c\
	$(SRC_PATH)/BLAS_zsymv2_z_c_x.c\
	$(SRC_PATH)/BLAS_zsymv2_c_z_x.c\
	$(SRC_PATH)/BLAS_zsymv2_c_c_x.c\
	$(SRC_PATH)/BLAS_csymv2_c_s_x.c\
	$(SRC_PATH)/BLAS_csymv2_s_c_x.c\
	$(SRC_PATH)/BLAS_csymv2_s_s_x.c\
	$(SRC_PATH)/BLAS_zsymv2_z_d_x.c\
	$(SRC_PATH)/BLAS_zsymv2_d_z_x.c\
	$(SRC_PATH)/BLAS_zsymv2_d_d_x.c \
	\
	$(SRC_PATH)/BLAS_ssymv2-f2c.c\
	$(SRC_PATH)/BLAS_dsymv2-f2c.c\
	$(SRC_PATH)/BLAS_csymv2-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2-f2c.c\
	$(SRC_PATH)/BLAS_dsymv2_d_s-f2c.c\
	$(SRC_PATH)/BLAS_dsymv2_s_d-f2c.c\
	$(SRC_PATH)/BLAS_dsymv2_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2_c_c-f2c.c\
	$(SRC_PATH)/BLAS_csymv2_c_s-f2c.c\
	$(SRC_PATH)/BLAS_csymv2_s_c-f2c.c\
	$(SRC_PATH)/BLAS_csymv2_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2_z_d-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2_d_z-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2_d_d-f2c.c\
	$(SRC_PATH)/BLAS_ssymv2_x-f2c.c\
	$(SRC_PATH)/BLAS_dsymv2_x-f2c.c\
	$(SRC_PATH)/BLAS_csymv2_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2_x-f2c.c\
	$(SRC_PATH)/BLAS_dsymv2_d_s_x-f2c.c\
	$(SRC_PATH)/BLAS_dsymv2_s_d_x-f2c.c\
	$(SRC_PATH)/BLAS_dsymv2_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_csymv2_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_csymv2_s_c_x-f2c.c\
	$(SRC_PATH)/BLAS_csymv2_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2_z_d_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2_d_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymv2_d_d_x-f2c.c 

SYMV2_AMB_SRCS =\
	$(SRC_PATH)/BLAS_SYMV2_AMB_ARGSsymv2__.c \
	\
	$(SRC_PATH)/BLAS_SYMV2_AMB_ARGSsymv2__-f2c.c 

SYMV2_OBJS = $(SYMV2_SRCS:.c=.o)

SYMV2_AMB_OBJS = $(SYMV2_AMB_SRCS:.c=.o)

all: $(SYMV2_OBJS)

all-amb: $(SYMV2_AMB_OBJS)

lib: $(SYMV2_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(SYMV2_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

lib-amb: $(SYMV2_AMB_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB_AMB) $(SYMV2_AMB_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB_AMB)

.c.o: 
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK

