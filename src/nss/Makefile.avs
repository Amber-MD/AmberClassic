###
#
# Makefile for AVS symmetry server modules
#
###

AVS_LIBS=	$(AVS_PATH)/lib
AVS_INC=	$(AVS_PATH)/include
NAB_INC=	$(NABHOME)/include
NAB_LIBS=	$(NABHOME)/lib/$(ARCH)
# SYM_LIBS=	-L$(SYM_HOME) -lsym -L$(NAB_LIBS) -lnab
SYM_LIBS=	-L$(NAB_LIBS) -lsym -lnab
BASELIBS=	-lgeom -lutil -lm $(LASTLIBS)
CFLAGS=		-I$(AVS_INC) -I$(NAB_INC) $(AOPTCFLAGS) $(LOCAL_CFLAGS) $(G)
CFLOWLIBS=	-L$(AVS_LIBS) -lflow_c $(BASELIBS)

avs:	\
	ma		\
	mg_2_fold	\
	mg_3_fold	\
	mg_4_fold	\
	mg_5_fold	\
	mg_6_fold	\
	mg_N_fold	\
	mg_cube		\
	mg_dihedral	\
	mg_helix	\
	mg_ico		\
	mg_octa		\
	mg_orient	\
	mg_orient2	\
	mg_rand		\
	mg_rotate	\
	mg_tetra	\
	mg_translate	\
	mm		\
	mr		\
	mw		\
	mx		\
	ron

ma:	ma.o
	$(CC) $(CFLAGS) -o ma ma.o $(SYM_LIBS) $(CFLOWLIBS)

mg_2_fold:	mg_2_fold.o
	$(CC) $(CFLAGS) -o mg_2_fold mg_2_fold.o $(SYM_LIBS) $(CFLOWLIBS)

mg_3_fold:	mg_3_fold.o
	$(CC) $(CFLAGS) -o mg_3_fold mg_3_fold.o $(SYM_LIBS) $(CFLOWLIBS)

mg_4_fold:	mg_4_fold.o
	$(CC) $(CFLAGS) -o mg_4_fold mg_4_fold.o $(SYM_LIBS) $(CFLOWLIBS)

mg_5_fold:	mg_5_fold.o
	$(CC) $(CFLAGS) -o mg_5_fold mg_5_fold.o $(SYM_LIBS) $(CFLOWLIBS)

mg_6_fold:	mg_6_fold.o
	$(CC) $(CFLAGS) -o mg_6_fold mg_6_fold.o $(SYM_LIBS) $(CFLOWLIBS)

mg_N_fold:	mg_N_fold.o
	$(CC) $(CFLAGS) -o mg_N_fold mg_N_fold.o $(SYM_LIBS) $(CFLOWLIBS)

mg_cube:	mg_cube.o
	$(CC) $(CFLAGS) -o mg_cube mg_cube.o $(SYM_LIBS) $(CFLOWLIBS)

mg_dihedral:	mg_dihedral.o
	$(CC) $(CFLAGS) -o mg_dihedral mg_dihedral.o $(SYM_LIBS) $(CFLOWLIBS)

mg_helix:	mg_helix.o
	$(CC) $(CFLAGS) -o mg_helix mg_helix.o $(SYM_LIBS) $(CFLOWLIBS)

mg_ico:	mg_ico.o
	$(CC) $(CFLAGS) -o mg_ico mg_ico.o $(SYM_LIBS) $(CFLOWLIBS)

mg_octa:	mg_octa.o
	$(CC) $(CFLAGS) -o mg_octa mg_octa.o $(SYM_LIBS) $(CFLOWLIBS)

mg_orient:	mg_orient.o
	$(CC) $(CFLAGS) -o mg_orient mg_orient.o $(SYM_LIBS) $(CFLOWLIBS)

mg_orient2:	mg_orient2.o
	$(CC) $(CFLAGS) -o mg_orient2 mg_orient2.o $(SYM_LIBS) $(CFLOWLIBS)

mg_rand:	mg_rand.o
	$(CC) $(CFLAGS) -o mg_rand mg_rand.o $(SYM_LIBS) $(CFLOWLIBS)

mg_rotate:	mg_rotate.o
	$(CC) $(CFLAGS) -o mg_rotate mg_rotate.o $(SYM_LIBS) $(CFLOWLIBS)

mg_tetra:	mg_tetra.o
	$(CC) $(CFLAGS) -o mg_tetra mg_tetra.o $(SYM_LIBS) $(CFLOWLIBS)

mg_translate:	mg_translate.o
	$(CC) $(CFLAGS) -o mg_translate mg_translate.o $(SYM_LIBS) $(CFLOWLIBS)

mm:	mm.o
	$(CC) $(CFLAGS) -o mm mm.o $(SYM_LIBS) $(CFLOWLIBS)

mr:	mr.o
	$(CC) $(CFLAGS) -o mr mr.o $(SYM_LIBS) $(CFLOWLIBS)

mw:	mw.o
	$(CC) $(CFLAGS) -o mw mw.o $(SYM_LIBS) $(CFLOWLIBS)

mx:	mx.o
	$(CC) $(CFLAGS) -o mx mx.o $(SYM_LIBS) $(CFLOWLIBS)

ron:	ron.o
	$(CC) $(CFLAGS) -o ron ron.o $(SYM_LIBS) $(CFLOWLIBS)

clean:
	rm -f		\
	core 		\
	*.o		\
	ma		\
	mg_2_fold	\
	mg_3_fold	\
	mg_4_fold	\
	mg_5_fold	\
	mg_6_fold	\
	mg_N_fold	\
	mg_cube		\
	mg_dihedral	\
	mg_helix	\
	mg_ico		\
	mg_octa		\
	mg_orient	\
	mg_orient2	\
	mg_rand		\
	mg_rotate	\
	mg_tetra	\
	mg_translate	\
	mm		\
	mr		\
	mw		\
	mx		\
	ron
