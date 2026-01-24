#ifndef NonbondedHeadings
#define NonbondedHeadings

#include "pmeDirectDS.h"
#include "TopologyDS.h"
#include "CompFrcDS.h"
#include "CellManipDS.h"
#include "TrajectoryDS.h"

int TestBondedExclusion(int a1, int a2, prmtop *tp);

void CompSameVgtECR(cell *C, atomc *atm1, int natm1, FrcTab *EFrc,
		    dircon *dcinp, prmtop *tp);

void CompIntrVgtECR(cell *C, atomc *atm1, atomc *atm2, int *ordr1, int *ordr2,
		    int natm1, int natm2, FrcTab *EFrc, dircon *dcinp,
		    prmtop *tp);

void CompIntrVgtE(cell *C, atomc *atm1, atomc *atm2, int *ordr1, int *ordr2,
		  int natm1, int natm2, FrcTab *EFrc, dircon *dcinp,
		  prmtop *tp);

void CompSameVeqECR(cell *C, atomc *atm1, int natm1, FrcTab *EFrc,
		    dircon *dcinp, prmtop *tp);

void CompIntrVeqECR(cell *C, atomc *atm1, atomc *atm2, int *ordr1, int *ordr2,
		    int natm1, int natm2, FrcTab *EFrc, dircon *dcinp,
		    prmtop *tp);

void CompIntrVeqE(cell *C, atomc *atm1, atomc *atm2, int *ordr1, int *ordr2,
		  int natm1, int natm2, FrcTab *EFrc, dircon *dcinp,
		  prmtop *tp);

void CompSameVgtECRnrg(cell *C, atomc *atm1, int natm1, FrcTab *EFrc,
		       dircon *dcinp, prmtop *tp, Energy *sysU);

void CompIntrVgtECRnrg(cell *C, atomc *atm1, atomc *atm2, int *ordr1,
		       int *ordr2, int natm1, int natm2, FrcTab *EFrc,
		       dircon *dcinp, prmtop *tp, Energy *sysU);

void CompIntrVgtEnrg(cell *C, atomc *atm1, atomc *atm2, int *ordr1, int *ordr2,
		     int natm1, int natm2, FrcTab *EFrc, dircon *dcinp,
		     prmtop *tp, Energy *sysU);

void CompSameVeqECRnrg(cell *C, atomc *atm1, int natm1, FrcTab *EFrc,
		       dircon *dcinp, prmtop *tp, Energy *sysU);

void CompIntrVeqECRnrg(cell *C, atomc *atm1, atomc *atm2, int *ordr1,
		       int *ordr2, int natm1, int natm2, FrcTab *EFrc,
		       dircon *dcinp, prmtop *tp, Energy *sysU);

void CompIntrVeqEnrg(cell *C, atomc *atm1, atomc *atm2, int *ordr1, int *ordr2,
		     int natm1, int natm2, FrcTab *EFrc, dircon *dcinp,
		     prmtop *tp, Energy *sysU);

void CompSameVgtECRnrgvir(cell *C, atomc *atm1, int natm1, FrcTab *EFrc,
			  dircon *dcinp, prmtop *tp, Energy *sysU);

void CompIntrVgtECRnrgvir(cell *C, atomc *atm1, atomc *atm2, int *ordr1,
			  int *ordr2, int natm1, int natm2, FrcTab *EFrc,
			  dircon *dcinp, prmtop *tp, Energy *sysU);

void CompIntrVgtEnrgvir(cell *C, atomc *atm1, atomc *atm2, int *ordr1,
			int *ordr2, int natm1, int natm2, FrcTab *EFrc,
			dircon *dcinp,prmtop *tp, Energy *sysU);

void CompSameVeqECRnrgvir(cell *C, atomc *atm1, int natm1, FrcTab *EFrc,
			  dircon *dcinp, prmtop *tp, Energy *sysU);

void CompIntrVeqECRnrgvir(cell *C, atomc *atm1, atomc *atm2, int *ordr1,
			  int *ordr2, int natm1, int natm2, FrcTab *EFrc,
			  dircon *dcinp, prmtop *tp, Energy *sysU);

void CompIntrVeqEnrgvir(cell *C, atomc *atm1, atomc *atm2, int *ordr1,
			int *ordr2, int natm1, int natm2, FrcTab *EFrc,
			dircon *dcinp, prmtop *tp, Energy *sysU);

void CompSameVgtECRnrgx(cell *C, atomc *atm1, int natm1, FrcTab *EFrc,
			dircon *dcinp, prmtop *tp, Energy *sysU);

void CompIntrVgtECRnrgx(cell *C, atomc *atm1, atomc *atm2,int *ordr1,
			int *ordr2, int natm1, int natm2, FrcTab *EFrc,
			dircon *dcinp, prmtop *tp, Energy *sysU);

void CompIntrVgtEnrgx(cell *C, atomc *atm1, atomc *atm2, int *ordr1,
		      int *ordr2, int natm1, int natm2, FrcTab *EFrc,
		      dircon *dcinp, prmtop *tp, Energy *sysU);

void CompSameVeqECRnrgx(cell *C, atomc *atm1, int natm1, FrcTab *EFrc,
			dircon *dcinp, prmtop *tp, Energy *sysU);

void CompIntrVeqECRnrgx(cell *C, atomc *atm1, atomc *atm2, int *ordr1,
			int *ordr2, int natm1, int natm2, FrcTab *EFrc,
			dircon *dcinp, prmtop *tp, Energy *sysU);

void CompIntrVeqEnrgx(cell *C, atomc *atm1, atomc *atm2, int *ordr1,
		      int *ordr2, int natm1, int natm2, FrcTab *EFrc,
		      dircon *dcinp, prmtop *tp, Energy *sysU);

#endif
