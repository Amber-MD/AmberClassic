include ../../config.h

install:  randvec rmscorr tensorfit
	/bin/mv randvec rmscorr tensorfit $(BINDIR)

randvec:  randvec.o
	$(FC) $(FPPFLAGS) $(LDFLAGS) -o randvec randvec.o

rmscorr:  rmscorr.o
	$(FC) $(FPPFLAGS) $(LDFLAGS) -o rmscorr rmscorr.o

tensorfit:  tensorfit.o
	$(FC) $(FPPFLAGS) $(LDFLAGS) -o tensorfit tensorfit.o

clean:
	/bin/rm -f *.o randvec rmscorr rotdifsvd tensorfit 

uninstall:
	-rm -f $(BINDIR)/randvec
	-rm -f $(BINDIR)/rmscorr
	-rm -f $(BINDIR)/tensorfit

randvec.o:	randvec.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

rmscorr.o:	rmscorr.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) -o $@ $<

tensorfit.o:	tensorfit.F90
	$(FC) $(FPPFLAGS) -c $(FOPTFLAGS) $(FREEFORMAT_FLAG) -o $@ $<

