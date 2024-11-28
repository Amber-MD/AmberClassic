include ../../config.h

.SUFFIXES:
.SUFFIXES:  .c .o

LOCALFLAGS=-Wno-format

.c.o:
	@echo "[CIFPARSE]  CC $<"
	$(VB)$(CC) -c $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) $<

CIFOBS=	\
	lex.cifp.o \
	cifp.tab.o \
	cifparse.o

$(LIBDIR)/libcifparse.a: $(CIFOBS)
	@echo "[CIFPARSE]  AR $@"
	$(VB)/bin/rm -f $@
	$(VB)$(AR) $@ $(CIFOBS) > /dev/null
	$(VB)-$(RANLIB) $@

cifp.tab.h: cifp.tab.c

cifp.tab.c:  cifparse.y lex.cifp.c
	@echo "[CIFPARSE]  BISON $<"
	$(VB)bison -v -d -pcifp -o $@ cifparse.y

cifp.tab.o: cifp.tab.c
	@echo "[CIFPARSE]  CC $<"
	$(VB)$(CC) -c -D$(LEX) $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) cifp.tab.c

lex.cifp.c: cifparse.l 
	@echo "[CIFPARSE]  FLEX $<"
	$(VB)flex -Pcifp cifparse.l

install:  $(LIBDIR)/libcifparse.a

clean:
	/bin/rm -f *.o  cifp.output
	/bin/rm -f lex.cifp.c cifp.tab.c cifp.tab.h

uninstall: clean
	/bin/rm -f $(LIBDIR)/libcifparse.a

cif2cif:  cif2cif.c $(LIBDIR)/libcifparse.a
	$(CC) -o cif2cif cif2cif.c $(LIBDIR)/libcifparse.a

nef2nef:  nef2nef.c $(LIBDIR)/libcifparse.a
	$(CC) -o nef2nef nef2nef.c $(LIBDIR)/libcifparse.a

lex.cifp.o:	\
	cifp.tab.h
