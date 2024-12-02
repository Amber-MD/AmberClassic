include ../../config.h

SUFFIXES:
SUFFIXES: .c .cu .o

test:  test_mask test_xmin test_md test_nmode

test_xmin:  txmin
	./txmin xfin.parm7 xfin.rst7 foo.rst7 > txmin.out
	../dacdif -r 1.e-5 txmin.out.check txmin.out
	@/bin/rm -f foo.rst7 txmin txmin.o

test_lmod:  tlmod
	./tlmod trpcage.parm7 trpcage.rst7 > tlmod.out
	../dacdif tlmod.out.check tlmod.out
	@/bin/rm -f conflib.dat tlmod tlmod.o

test_md:  tmd
	./tmd xfin.parm7 xfin.md1.rst7 > tmd.out
	../dacdif tmd.out.check tmd.out
	@/bin/rm -f xfin.md2.x tmd tmd.o

test_md.MPI:  tmd.MPI
	time mpirun -np 4 ./tmd.MPI xfin.parm7 xfin.md1.rst7

test_nmode:  tnmode
	./tnmode crambin.parm7 crambin.rst7 > tnmode.out
	../dacdif tnmode.out.check tnmode.out
	@/bin/rm -f vecs tnmode tnmode.o

test_mask:  tmask
	./tmask xfin.parm7 xfin.md1.x ":1-3" > tmask.out
	../dacdif tmask.out.check tmask.out
	@/bin/rm -f tmask tmask.o

test_hcp1:  thcp1
	./thcp1

test_hcp2:  thcp2
	./thcp2

test_nabc:  tnabc
	./tnabc 2.pdb 2.rst7
	@/bin/rm -f tnabc tnabc.o

test_netcdf:  tnetcdf
	./tnetcdf

test_cif:  tcif
	./tcif
	/bin/rm tcif tcif.o

LIBS = -L$(LIBDIR) $(FLIBS) $(NETCDFLIBF)  $(LM) 

#---------------------------------------------------------------------------
program_2:  program_2.o
	@echo "[NABC] CC $@"
	$(VB)$(CC) -o $@ $< $(LIBS)

txmin:  txmin.o
	@echo "[NABC] CC $@"
	$(VB)$(CC) -o $@ $< $(LIBS) $(FLIBS_ARCH)

tlmod:  tlmod.o
	@echo "[NABC] CC $@"
	$(VB)$(CC) -o $@ $< $(LIBS) $(FLIBS_ARCH)
                
tmd:  tmd.o
	@echo "[NABC] CC $@"
	$(VB)$(CC) -o $@ $<  $(LIBS) $(FLIBS_ARCH)

tmd.MPI:  tmd.MPI.o
	@echo "[NABC] CC $@"
	$(VB)mpicc -o $@ $<  $(LIBS) $(FLIBS_ARCH)

tnmode:  tnmode.o
	@echo "[NABC] CC $@"
	$(VB)$(CC) -o $@ $<  $(LIBS) $(FLIBS_ARCH)

tmask:  tmask.o
	@echo "[NABC] CC $@"
	$(VB)$(CC) -o $@ $<  $(LIBS) $(FLIBS_ARCH)

thcp1:  thcp1.o
	@echo "[NABC] CC $@"
	$(VB)$(CC) -o $@ $<  $(LIBS) $(FLIBS_ARCH)

thcp2:  thcp2.o
	@echo "[NABC] CC $@"
	$(VB)$(CC) -o $@ $<  $(LIBS) $(FLIBS_ARCH)

tnabc:  tnabc.o
	@echo "[NABC] CC $@"
	$(VB)$(CC) -o $@ $<  $(LIBS) $(FLIBS_ARCH)

tcif:  tcif.o
	@echo "[NABC] CC $@"
	$(VB)$(CC) -o $@ $<  -L$(LIBDIR)  -lnabc -lcifparse $(LIBS) $(FLIBS_ARCH)

tnetcdf:  tnetcdf.o
	@echo "[NABC] CC $@"
	$(VB)$(CC) -o $@ $<  $(LIBS) $(FLIBS_ARCH)

#---------------------------------------------------------------------------
#  non-tests:

fd:  fd.c
	@echo "[NABC] CC $<"
	$(VB)$(CC) $(CFLAGS) -I$(INCDIR) -o $@ $< -L$(LIBDIR) \
        -lnabc -lcifparse  $(LIBS)

wch:  wch.c
	@echo "[NABC] CC $<"
	$(VB)$(CC) $(CFLAGS) -I$(INCDIR) -o $@ $< -L$(LIBDIR) \
        -lnabc -lcifparse $(LIBS)

circle:  circle.c
	@echo "[NABC] CC $<"
	$(VB)$(CC) $(CFLAGS) -I$(INCDIR) -o $@ $< -L$(LIBDIR) \
        -lnabc -lcifparse $(LIBS)

loop:  loop.c
	@echo "[NABC] CC $<"
	$(VB)$(CC) $(CFLAGS) -I$(INCDIR) -o $@ $< -L$(LIBDIR) \
        -lnabc -lcifparse $(LIBS)

#---------------------------------------------------------------------------
%.o: %.c
	@echo "[NABC] CC $<"
	$(VB)$(CC) $(COPTFLAGS) $(CFLAGS) -I$(INCDIR) -c $< -o $@ 

%.o: %.cu
	@echo "[NABC] nvcc $<"
	nvcc $(NVFLAGS) -I$(INCDIR) -c $< -o $@

tmd.MPI.o: tmd.c
	@echo "[NABC] CC $<"
	$(VB)mpicc -DMPI -c $(COPTFLAGS) $(CFLAGS) -I$(INCDIR) -o $@ $<

clean:
	/bin/rm -f *.o txmin tmd tmd.MPI tnabc tlmod tnmode tmask thcp1 thcp2 \
       trpcageTest.nc
