include ../../config.h

####
#
#       Make stuff:
#
NAB=$(BINDIR)/nab

.SUFFIXES:
.SUFFIXES:	.nab .c .o

.nab.o:
	$(NAB) -c $<
	@rm $*.c

.c.o:
	$(NAB) -c $<
#
#       END of Make stuff
#
####

SYMOBJS=     		\
	symop.o		\
	tss.o

all:	libsym	\
	matextract$(SFX)	\
	matgen$(SFX)  		\
	matmerge$(SFX)		\
	matmul$(SFX)		\
	transform$(SFX)		\
	tss_init$(SFX)		\
	tss_main$(SFX)		\
	tss_next$(SFX)

install:  all
	mv matextract$(SFX) matgen$(SFX) matmerge$(SFX) matmul$(SFX) \
		transform$(SFX) tss_init$(SFX) tss_main$(SFX) \
		tss_next$(SFX) $(BINDIR)

uninstall: clean
	/bin/rm -f $(BINDIR)/matextract$(SFX) $(BINDIR)/matgen$(SFX) \
		$(BINDIR)/matmerge$(SFX) \
		$(BINDIR)/matmul$(SFX) $(BINDIR)/transform$(SFX) \
		$(BINDIR)/tss_init$(SFX) \
		$(BINDIR)/tss_main$(SFX) $(BINDIR)/tss_next$(SFX)

libsym:	$(SYMOBJS) 
	$(AR) $(LIBDIR)/libnab.a $(SYMOBJS)
	-ranlib $(LIBDIR)/libnab.a

matextract$(SFX): matextract.o
	$(NAB) -o matextract$(SFX) matextract.o

matgen$(SFX):	matgen.o
	$(NAB) -o matgen$(SFX) matgen.o

matmerge$(SFX):	matmerge.o
	$(NAB) -o matmerge$(SFX) matmerge.o

matmul$(SFX):	matmul.o
	$(NAB) -o matmul$(SFX) matmul.o

transform$(SFX):	transform.o
	$(NAB) -o transform$(SFX) transform.o

tss_init$(SFX):	tss_init.o
	$(NAB) -o tss_init$(SFX) tss_init.o

tss_main$(SFX):	tss_main.o
	$(NAB) -o tss_main$(SFX) tss_main.o

tss_next$(SFX):	tss_next.o
	$(NAB) -o tss_next$(SFX) tss_next.o

clean:
	/bin/rm -f	\
	*.o 		\
	matextract$(SFX)	\
	matgen$(SFX)		\
	matmerge$(SFX)	\
	matmul$(SFX)		\
	transform$(SFX)	\
	tss_init$(SFX)	\
	tss_main$(SFX)	\
	tss_next$(SFX)	\
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

test: matgen transform
	@echo =====================================================
	@echo   "Running test create an icosahderal virus"
	@echo
	@./matgen -create ../test/ico.sym |\
		 ./transform -pdb ../test/tbv.com.pdb > tbv.60.com.pdb
	@-diff ../test/tbv.60.com.pdb.check tbv.60.com.pdb \
		&& echo "	PASSED" \
		|| echo "	FAILED (possibly: see if diffs above look OK)"
	@/bin/rm -f tbv.60.com.pdb

matextract.o: libsym
matgen.o: libsym
matmerge.o: libsym
matmul.o: libsym
transform.o: libsym
tss_init.o: libsym
tss_main.o: libsym
tss_next.o: libsym
