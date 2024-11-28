include ../../make.conf
include ../../$(MAKEINC)

LIB_PATH = ../../$(OUTPUT_DIR)
SRC_PATH = .
HEADER_PATH = ..

GE_SUM_MV_SRCS =\
	$(SRC_PATH)/BLAS_sge_sum_mv.c\
	$(SRC_PATH)/BLAS_dge_sum_mv.c\
	$(SRC_PATH)/BLAS_cge_sum_mv.c\
	$(SRC_PATH)/BLAS_zge_sum_mv.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_d_s.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_s_d.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_s_s.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_z_c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_c_z.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_c_c.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_c_s.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_s_c.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_s_s.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_z_d.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_d_z.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_d_d.c\
	$(SRC_PATH)/BLAS_sge_sum_mv_x.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_x.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_x.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_x.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_d_s_x.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_s_d_x.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_s_s_x.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_z_c_x.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_c_z_x.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_c_c_x.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_c_s_x.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_s_c_x.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_s_s_x.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_z_d_x.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_d_z_x.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_d_d_x.c \
	\
	$(SRC_PATH)/BLAS_sge_sum_mv-f2c.c\
	$(SRC_PATH)/BLAS_dge_sum_mv-f2c.c\
	$(SRC_PATH)/BLAS_cge_sum_mv-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv-f2c.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_d_s-f2c.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_s_d-f2c.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_z_c-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_c_z-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_c_c-f2c.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_c_s-f2c.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_s_c-f2c.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_s_s-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_z_d-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_d_z-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_d_d-f2c.c\
	$(SRC_PATH)/BLAS_sge_sum_mv_x-f2c.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_x-f2c.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_x-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_x-f2c.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_d_s_x-f2c.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_s_d_x-f2c.c\
	$(SRC_PATH)/BLAS_dge_sum_mv_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_z_c_x-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_c_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_c_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_c_s_x-f2c.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_s_c_x-f2c.c\
	$(SRC_PATH)/BLAS_cge_sum_mv_s_s_x-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_z_d_x-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_d_z_x-f2c.c\
	$(SRC_PATH)/BLAS_zge_sum_mv_d_d_x-f2c.c 

GE_SUM_MV_AMB_SRCS =\
	$(SRC_PATH)/BLAS_GE_SUM_MV_AMB_ARGSge_sum_mv__.c \
	\
	$(SRC_PATH)/BLAS_GE_SUM_MV_AMB_ARGSge_sum_mv__-f2c.c 

GE_SUM_MV_OBJS = $(GE_SUM_MV_SRCS:.c=.o)

GE_SUM_MV_AMB_OBJS = $(GE_SUM_MV_AMB_SRCS:.c=.o)

all: $(GE_SUM_MV_OBJS)

all-amb: $(GE_SUM_MV_AMB_OBJS)

lib: $(GE_SUM_MV_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(GE_SUM_MV_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

lib-amb: $(GE_SUM_MV_AMB_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB_AMB) $(GE_SUM_MV_AMB_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB_AMB)

.c.o: 
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK

