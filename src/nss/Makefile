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
	matextract	\
	matgen  		\
	matmerge		\
	matmul		\
	transform		\
	tss_init		\
	tss_main		\
	tss_next

install_nss:  all
	mv matextract matgen matmerge matmul \
		transform tss_init tss_main \
		tss_next $(BINDIR)

uninstall: clean
	/bin/rm -f $(BINDIR)/matextract $(BINDIR)/matgen \
		$(BINDIR)/matmerge \
		$(BINDIR)/matmul $(BINDIR)/transform \
		$(BINDIR)/tss_init \
		$(BINDIR)/tss_main $(BINDIR)/tss_next

libsym:	$(SYMOBJS) 
	$(AR) $(LIBDIR)/libnab.a $(SYMOBJS)
	-ranlib $(LIBDIR)/libnab.a

matextract: matextract.o
	$(NAB) -o matextract matextract.o

matgen:	matgen.o
	$(NAB) -o matgen matgen.o

matmerge:	matmerge.o
	$(NAB) -o matmerge matmerge.o

matmul:	matmul.o
	$(NAB) -o matmul matmul.o

transform:	transform.o
	$(NAB) -o transform transform.o

tss_init:	tss_init.o
	$(NAB) -o tss_init tss_init.o

tss_main:	tss_main.o
	$(NAB) -o tss_main tss_main.o

tss_next:	tss_next.o
	$(NAB) -o tss_next tss_next.o

clean:
	/bin/rm -f	\
	*.o 		\
	matextract	\
	matgen		\
	matmerge	\
	matmul		\
	transform	\
	tss_init	\
	tss_main	\
	tss_next
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
