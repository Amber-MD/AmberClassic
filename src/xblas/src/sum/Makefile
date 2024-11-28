include ../../make.conf
include ../../$(MAKEINC)

LIB_PATH = ../../$(OUTPUT_DIR)
SRC_PATH = .
HEADER_PATH = ..

SUM_SRCS =\
	$(SRC_PATH)/BLAS_ssum.c\
	$(SRC_PATH)/BLAS_dsum.c\
	$(SRC_PATH)/BLAS_csum.c\
	$(SRC_PATH)/BLAS_zsum.c\
	$(SRC_PATH)/BLAS_ssum_x.c\
	$(SRC_PATH)/BLAS_dsum_x.c\
	$(SRC_PATH)/BLAS_csum_x.c\
	$(SRC_PATH)/BLAS_zsum_x.c \
	\
	$(SRC_PATH)/BLAS_ssum-f2c.c\
	$(SRC_PATH)/BLAS_dsum-f2c.c\
	$(SRC_PATH)/BLAS_csum-f2c.c\
	$(SRC_PATH)/BLAS_zsum-f2c.c\
	$(SRC_PATH)/BLAS_ssum_x-f2c.c\
	$(SRC_PATH)/BLAS_dsum_x-f2c.c\
	$(SRC_PATH)/BLAS_csum_x-f2c.c\
	$(SRC_PATH)/BLAS_zsum_x-f2c.c 

SUM_AMB_SRCS =\
	$(SRC_PATH)/BLAS_SUM_AMB_ARGSsum.c \
	\
	$(SRC_PATH)/BLAS_SUM_AMB_ARGSsum-f2c.c 

SUM_OBJS = $(SUM_SRCS:.c=.o)

SUM_AMB_OBJS = $(SUM_AMB_SRCS:.c=.o)

all: $(SUM_OBJS)

all-amb: $(SUM_AMB_OBJS)

lib: $(SUM_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(SUM_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

lib-amb: $(SUM_AMB_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB_AMB) $(SUM_AMB_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB_AMB)

.c.o: 
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK

