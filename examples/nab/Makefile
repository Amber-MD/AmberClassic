include ../../../config.h

NAB=$(BINDIR)/nab

##############################################################################
######### no configuration should be needed beyond this line #################
##############################################################################

all: p1 p2 p3 p4 p7 p8 p9 p10 p11 p12 p14

p1:
	$(NAB) $(NABFLAGS) -o program_1$(SFX) program_1.nab

p2:
	$(NAB) $(NABFLAGS) -o program_2$(SFX) program_2.nab

p3:
	$(NAB) $(NABFLAGS) -o program_3$(SFX) program_3.nab

p4:
	$(NAB) $(NABFLAGS) -o program_4$(SFX) program_4.nab

p5:
	$(NAB) $(NABFLAGS) -o program_5$(SFX) program_5.nab

p6:
	$(NAB) $(NABFLAGS) -o program_6$(SFX) program_6.nab

p7:
	$(NAB) $(NABFLAGS) -o program_7$(SFX) program_7.nab

p8:
	$(NAB) $(NABFLAGS) -o program_8$(SFX) program_8.nab

p9: 
	$(NAB) $(NABFLAGS) -o program_9$(SFX) program_9.nab getbase.nab

p10: 
	$(NAB) $(NABFLAGS) -o program_10$(SFX) program_10.nab getbase.nab

p11: 
	$(NAB) $(NABFLAGS) -o program_11$(SFX) program_11.nab spline.nab

p12: 
	$(NAB) $(NABFLAGS) -o program_12$(SFX) program_12.nab getbase.nab program_13.nab

p14:
	$(NAB) $(NABFLAGS) -o program_14$(SFX) program_14.nab

test.p1: p1
	@echo =====================================================
	@echo  "Testing example: program_1.nab"
	@NABHOME=`pwd`/.. ; export NABHOME;
	@-./program_1$(SFX) && echo "     PASSED" || echo "**** FAILED ****"
#	@-rm -f program_1$(SFX)

test.p2: p2
	@echo =====================================================
	@echo  "Testing example: program_2.nab"
	@./program_2$(SFX) > program_2.test;
	@-diff program_2.test program_2.ref && echo "     PASSED" || echo "**** FAILED ****"
#	@-rm -f program_2$(SFX) program_2.test

test.p3: p3
	@echo =====================================================
	@echo  "Testing example: program_3.nab"
	@-./program_3$(SFX) && echo "     PASSED" || echo "**** FAILED ****"
#	@-rm -f program_3$(SFX)

test.p4: p4
	@echo =====================================================
	@echo  "Testing example: program_4.nab"
	@-./program_4$(SFX) 2.0 > program_4.test;
	@-diff program_4.ref program_4.test && echo "     PASSED" || echo "**** FAILED ****"
#	@-rm -f program_4$(SFX) program_4.test

test.p5: p5
	@echo =====================================================
	@echo  "Testing example: program_5.nab"
	@-./program_5$(SFX) < program_5.tat.in > program_5.test;
	@-diff taDT.energy.ref taDT.energy.dat && echo "     PASSED" || echo "**** FAILED ****"
#	@-rm -f program_5$(SFX) taDT.energy.dat program_5.test

clean::
	rm -f *.o *.c program_1$(SFX) program_2$(SFX) program_3$(SFX) \
		program_4$(SFX) program_5$(SFX) program_6$(SFX) program_7$(SFX) \
		program_8$(SFX) program_9$(SFX) program_10$(SFX) program_11$(SFX) \
		program_12$(SFX) program_14$(SFX) *.test test.sup.pdb *.std.pdb

.SUFFIXES:  .nab

#  default rules for compilation:

.nab.o:
	$(NAB) -c $(NABFLAGS) -I../include $<

