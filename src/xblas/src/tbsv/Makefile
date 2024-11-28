include ../../make.conf
include ../../$(MAKEINC)

LIB_PATH = ../../$(OUTPUT_DIR)
SRC_PATH = .
HEADER_PATH = ..

TBSV_SRCS =\
	$(SRC_PATH)/BLAS_stbsv.c\
	$(SRC_PATH)/BLAS_dtbsv.c\
	$(SRC_PATH)/BLAS_ctbsv.c\
	$(SRC_PATH)/BLAS_ztbsv.c\
	$(SRC_PATH)/BLAS_dtbsv_s.c\
	$(SRC_PATH)/BLAS_ztbsv_c.c\
	$(SRC_PATH)/BLAS_ctbsv_s.c\
	$(SRC_PATH)/BLAS_ztbsv_d.c\
	$(SRC_PATH)/BLAS_stbsv_x.c\
	$(SRC_PATH)/BLAS_dtbsv_x.c\
	$(SRC_PATH)/BLAS_dtbsv_s_x.c\
	$(SRC_PATH)/BLAS_ctbsv_x.c\
	$(SRC_PATH)/BLAS_ztbsv_x.c\
	$(SRC_PATH)/BLAS_ztbsv_c_x.c\
	$(SRC_PATH)/BLAS_ctbsv_s_x.c\
	$(SRC_PATH)/BLAS_ztbsv_d_x.c \
	\
	$(SRC_PATH)/BLAS_stbsv-f2c.c\
	$(SRC_PATH)/BLAS_dtbsv-f2c.c\
	$(SRC_PATH)/BLAS_ctbsv-f2c.c\
	$(SRC_PATH)/BLAS_ztbsv-f2c.c\
	$(SRC_PATH)/BLAS_dtbsv_s-f2c.c\
	$(SRC_PATH)/BLAS_ztbsv_c-f2c.c\
	$(SRC_PATH)/BLAS_ctbsv_s-f2c.c\
	$(SRC_PATH)/BLAS_ztbsv_d-f2c.c\
	$(SRC_PATH)/BLAS_stbsv_x-f2c.c\
	$(SRC_PATH)/BLAS_dtbsv_x-f2c.c\
	$(SRC_PATH)/BLAS_dtbsv_s_x-f2c.c\
	$(SRC_PATH)/BLAS_ctbsv_x-f2c.c\
	$(SRC_PATH)/BLAS_ztbsv_x-f2c.c\
	$(SRC_PATH)/BLAS_ztbsv_c_x-f2c.c\
	$(SRC_PATH)/BLAS_ctbsv_s_x-f2c.c\
	$(SRC_PATH)/BLAS_ztbsv_d_x-f2c.c 

TBSV_AMB_SRCS =\
	$(SRC_PATH)/BLAS_TBSV_AMB_ARGStbsv_.c \
	\
	$(SRC_PATH)/BLAS_TBSV_AMB_ARGStbsv_-f2c.c 

TBSV_OBJS = $(TBSV_SRCS:.c=.o)

TBSV_AMB_OBJS = $(TBSV_AMB_SRCS:.c=.o)

all: $(TBSV_OBJS)

all-amb: $(TBSV_AMB_OBJS)

lib: $(TBSV_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(TBSV_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

lib-amb: $(TBSV_AMB_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB_AMB) $(TBSV_AMB_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB_AMB)

.c.o: 
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK

