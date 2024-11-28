include ../../make.conf
include ../../$(MAKEINC)

LIB_PATH = ../../$(OUTPUT_DIR)
SRC_PATH = .
HEADER_PATH = ..

SYMM_SRCS =\
	$(SRC_PATH)/BLAS_ssymm.c\
	$(SRC_PATH)/BLAS_dsymm.c\
	$(SRC_PATH)/BLAS_csymm.c\
	$(SRC_PATH)/BLAS_zsymm.c\
	$(SRC_PATH)/BLAS_dsymm_d_s.c\
	$(SRC_PATH)/BLAS_dsymm_s_d.c\
	$(SRC_PATH)/BLAS_dsymm_s_s.c\
	$(SRC_PATH)/BLAS_zsymm_z_c.c\
	$(SRC_PATH)/BLAS_zsymm_c_z.c\
	$(SRC_PATH)/BLAS_zsymm_c_c.c\
	$(SRC_PATH)/BLAS_csymm_c_s.c\
	$(SRC_PATH)/BLAS_csymm_s_c.c\
	$(SRC_PATH)/BLAS_csymm_s_s.c\
	$(SRC_PATH)/BLAS_zsymm_z_d.c\
	$(SRC_PATH)/BLAS_zsymm_d_z.c\
	$(SRC_PATH)/BLAS_zsymm_d_d.c\
	$(SRC_PATH)/BLAS_ssymm_x.c\
	$(SRC_PATH)/BLAS_dsymm_x.c\
	$(SRC_PATH)/BLAS_csymm_x.c\
	$(SRC_PATH)/BLAS_zsymm_x.c\
	$(SRC_PATH)/BLAS_dsymm_d_s_x.c\
	$(SRC_PATH)/BLAS_dsymm_s_d_x.c\
	$(SRC_PATH)/BLAS_dsymm_s_s_x.c\
	$(SRC_PATH)/BLAS_zsymm_z_c_x.c\
	$(SRC_PATH)/BLAS_zsymm_c_z_x.c\
	$(SRC_PATH)/BLAS_zsymm_c_c_x.c\
	$(SRC_PATH)/BLAS_csymm_c_s_x.c\
	$(SRC_PATH)/BLAS_csymm_s_c_x.c\
	$(SRC_PATH)/BLAS_csymm_s_s_x.c\
	$(SRC_PATH)/BLAS_zsymm_z_d_x.c\
	$(SRC_PATH)/BLAS_zsymm_d_z_x.c\
	$(SRC_PATH)/BLAS_zsymm_d_d_x.c \
	\
	$(SRC_PATH)/BLAS_ssymm-f2c.c\
	$(SRC_PATH)/BLAS_dsymm-f2c.c\
	$(SRC_PATH)/BLAS_csymm-f2c.c\
	$(SRC_PATH)/BLAS_zsymm-f2c.c\
	$(SRC_PATH)/BLAS_dsymm_d_s-f2c.c\
	$(SRC_PATH)/BLAS_dsymm_s_d-f2c.c\
	$(SRC_PATH)/BLAS_dsymm_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zsymm_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zsymm_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zsymm_c_c-f2c.c\
	$(SRC_PATH)/BLAS_csymm_c_s-f2c.c\
	$(SRC_PATH)/BLAS_csymm_s_c-f2c.c\
	$(SRC_PATH)/BLAS_csymm_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zsymm_z_d-f2c.c\
	$(SRC_PATH)/BLAS_zsymm_d_z-f2c.c\
	$(SRC_PATH)/BLAS_zsymm_d_d-f2c.c\
	$(SRC_PATH)/BLAS_ssymm_x-f2c.c\
	$(SRC_PATH)/BLAS_dsymm_x-f2c.c\
	$(SRC_PATH)/BLAS_csymm_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymm_x-f2c.c\
	$(SRC_PATH)/BLAS_dsymm_d_s_x-f2c.c\
	$(SRC_PATH)/BLAS_dsymm_s_d_x-f2c.c\
	$(SRC_PATH)/BLAS_dsymm_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymm_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymm_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymm_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_csymm_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_csymm_s_c_x-f2c.c\
	$(SRC_PATH)/BLAS_csymm_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymm_z_d_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymm_d_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zsymm_d_d_x-f2c.c 

SYMM_AMB_SRCS =\
	$(SRC_PATH)/BLAS_SYMM_AMB_ARGSsymm__.c \
	\
	$(SRC_PATH)/BLAS_SYMM_AMB_ARGSsymm__-f2c.c 

SYMM_OBJS = $(SYMM_SRCS:.c=.o)

SYMM_AMB_OBJS = $(SYMM_AMB_SRCS:.c=.o)

all: $(SYMM_OBJS)

all-amb: $(SYMM_AMB_OBJS)

lib: $(SYMM_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(SYMM_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

lib-amb: $(SYMM_AMB_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB_AMB) $(SYMM_AMB_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB_AMB)

.c.o: 
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK

