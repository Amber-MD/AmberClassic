include ../config_testing.h

# -------- Targets from ../Makefile ----------------------------------------

test:  test.$(INSTALLTYPE)

test.serial:
	@./run_serial_tests.sh

test.conda:
	@./run_conda_tests.sh

test.parallel:
	@./run_parallel_tests.sh

test.cuda:
	@./run_cuda_tests.sh

test.mpicuda:
	@./run_mpicuda_tests.sh

# -------- Targets for run_* scripts, making the logfiles --------------------

test.allserial: test.general test.rism test.xray test.cryoem test.nmr \
          test.leap test.antechamber test.paramfit test.rism1d test.xtalutil \
          test.gbnsr6 test.nabc test.nab

test.allparallel: test.general test.rism test.xray test.cryoem test.nmr 

test.allcuda: test.xray

test.allconda: test.general test.rism test.xray test.cryoem test.nmr \
          test.leap test.antechamber test.paramfit test.rism1d test.xtalutil \
          test.gbnsr6

# -------- Lower-level test targets: ---------------------------------------

test.general:
	-(cd dhfr && ./Run.dhfr)
	-(cd dhfr && ./Run.dhfr.xmin)
	-(cd dhfr && ./Run.dhfr_synctraj)
	-(cd gb8_trx && ./Run.trxox)
	-(cd gb8_trx && ./Run.trxox_md)
	-(cd tip4p && ./Run.tip4p)
	-(cd tip4p && ./Run.tip4p_mcbar)
	-(cd middle-scheme/MD_Constr_Langevin_water && ./Run.LGV_LF)
	-(cd middle-scheme/MD_Unconstr_Langevin_water && ./Run.LGV_LF)
#	-(cd softcore && ./Run_sc)   # parallel only!
	-(cd gb_rna && ./Run.gbrna)
	-(cd gb_rna && ./Run.gbrna.xmin)
	-(cd gb_rna && ./Run.gbrna.sgld)
	-(cd gb_rna && ./Run.gbrna.sgldg)
#	-(cd LES && ./Run.PME_LES)

test.rism:
	-(cd rism3d/4lzt && ./Run.4lzt.pme)
	-(cd rism3d/2igd && ./Run.2igd.pme)
	-(cd xray/480d && ./Run.480d_pse2)

test.nmr:
	-(cd rdc && Run.nmr)
	-(cd rdc && Run.dip)
	-(cd rdc && Run.csa_min)
	-(cd rdc && Run.dip2)

test.xray:
#	-(cd xray/2wcu && ./Run.2wcu)
	-(cd xray/3onk && ./Run.3onk_ml)
	-(cd xray/3onk && ./Run.3onk_mlopt)
#	-(cd xray/480d && ./Run.480d_xray)
	-(cd xray/6o2h && ./Run.6o2h_min)
	-(cd xray/8dz7 && ./Run.8dz7)
	-(cd xray/4lzt && ./Run.4lzt)
	-(cd xray/2ol9 && ./Run.2ol9)
	-(cd xray/2ol9 && ./Run.2ol9.min)
#	-(cd xray/2msi && ./Run.2msi)

test.cryoem:
	-(cd cryoem/3j9d && ./Run.3j9d)
	-(cd cryoem/emap && ./Run.emap)

test.leap:
	-(cd leap && ./Run.cmds)
	-(cd leap && ./Run.mbondi3)

test.antechamber:
	-(cd antechamber/abcg2 && Run.abcg2)
	-(cd antechamber/fpph && Run.fpph)
	-(cd antechamber/parmchk2 && Run.parmchk2)
	-(cd antechamber/sustiva && Run.sustiva_g2)
	-(cd antechamber/tp && Run.tp2)

test.paramfit:
	-(cd paramfit/write_qminput_nma; ./run_test.x)
	-(cd paramfit/simplex_perfect_fit_NMA; ./run_test.x)
	-(cd paramfit/ga_perfect_fit_NMA; ./run_test.x)
	-(cd paramfit/dihedral_least_squares_NMA; ./run_test.x)

test.rism1d:
	-(cd rism1d/spc-kh && ./Run.spc)
	-(cd rism1d/spc-psen && ./Run.spc-nacl-3)
	-(cd rism1d/spc-hnc && ./Run.spc-lj)

test.xtalutil:
	-(cd xtalutil/AddToBox && ./Run.AddToBox)
#	-(cd xtalutil/CheckMD && ./Run.CheckMD)
	-(cd xtalutil/PropPDB && ./Run.PropPDB)
	-(cd xtalutil/UnitCell && ./Run.UnitCell)

test.gbnsr6:
	-(cd gbnsr6 && ./Run.gbnsr6.min)

test.nabc:
	-(cd nabc && $(MAKE) test)

test.nab:
	-(cd nab && $(MAKE) test)

clean::
	/bin/rm -f TEST_FAILURES.diff
