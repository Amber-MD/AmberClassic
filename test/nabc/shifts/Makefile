test: 6pti_test 132d_test 1acf_test 1mbc_test 2lef_test clean

6pti_test::
	@echo "====================================================="
	@echo "Running test for proton shifts in BPTI"
	@echo
	../../bin/shifts -readobs -sander -swap '::H*' 6pti > messages 
	@./dacdif -a 1.e-4 6pti.emp.save 6pti.emp
	@echo "Checking sander input file:"
	@./dacdif 6pti.shf.save 6pti.shf
	@-/bin/rm 6pti.bmrb

1mbc_test::
	@echo "====================================================="
	@echo "Running test for proton shifts on myoglobin"
	@echo
	../../bin/shifts ':64-68:H*' 1mbc > messages 
	@./dacdif 1mbc.emp.save 1mbc.emp
	@-/bin/rm 1mbc.bmrb

132d_test::
	@echo =====================================================
	@echo "Running test for proton shifts in a DNA duplex"
	@echo
	../../bin/shifts -reslib -readobs '::H*' 132d  > messages 2> 132d.errs
	@./dacdif -a 0.015 132d.emp.save 132d.emp
#	../bin/shifts -readobs '::H*' 132d  > messages 2> 132d.errs
#	@./dacdif -a 0.015 132d.noreslib.emp.save 132d.emp
#	@-/bin/rm 132d.bmrb

2lef_test::
	@echo =====================================================
	@echo "Running test for shifts in a protein-DNA duplex"
	@echo
	../../bin/shifts -reslib -noreduce '::H*' 2lef  > messages 2> 2lef.errs
	@./dacdif -a 0.015 2lef.emp.save 2lef.emp
	@-/bin/rm -f 2lef.bmrb messages 2lef.emp.rdb tleap.out 2lef.errs

1ehb_test::
	@echo =====================================================
	@echo "Running test for HN neural network shifts"
	@echo
	../../bin/shifts -HN 1ehb 
	@./dacdif 1ehb.HN.save 1ehb.HN
	@-/bin/rm 1ehb.bmrb 1ehb.emp

1acf_test::
	@echo =====================================================
	@echo "Running test for carbon/nitrogen shifts in a protein"
	@echo
	../../bin/shifts -qdb -refine 1acf  
	@sed 's/-0.00/ 0.00/g' 1acf.qdb  > 1acf.out1b
#	@sed 's/-0.00/ 0.00/g' 1acf.qdb2 > 1acf.out2b
	@./dacdif ../../src/shifts/examples/1acf.out1 1acf.out1b
	@-/bin/rm 1acf.qdb 1acf.qdb2

benzene_test::
	@echo =====================================================
	@echo "Running test for dummy points around benzene:"
	@echo
	../../bin/shifts -nocoil -noelec -csa :DUM: contact.2.5  | grep total > contact.2.5.emp
	@./dacdif contact.2.5.emp.save contact.2.5.emp

susc_fit_test::
	@echo =====================================================
	@echo "Running test to fit fluorobenzene susceptibility"
	@echo
	TEST_SHIFTS=1; export TEST_SHIFTS; ../bin/susc_fit 2.0 fbenzene.pdb fbG
	@./dacdif fbG.lsq.out.save fbG.lsq.out

clean:
	@/bin/rm -f *.par? *.out1 *.out2 *.errs TEST_FAILURES.diff messages \
        reduce.info pdb4amber.messages mdinfo *.rdb
