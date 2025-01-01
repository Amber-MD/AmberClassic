#ifndef BondedHeadings
#define BondedHeadings

#include "CrdManipDS.h"
#include "TopologyDS.h"
#include "CompFrcDS.h"
#include "TrajectoryDS.h"
#include "CellManipDS.h"

void CellBondedIntr(cell *C, cellgrid *CG, coord *crd, prmtop *tp,
		    FrcTab *Cfrc, FrcTab *Hfrc, Energy *sysUV, int qform);

void AttenuatePairVir(prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc, int atmA,
		      int atmB, double dx, double dy, double dz, double fmag,
		      double* afrc, double* bfrc, Energy *sysUV,
		      double elec14fac, double lj14fac, int qform);

void AttenuatePairFrcNrg(prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc, int atmA,
                         int atmB, double dx, double dy, double dz,
			 double fmag, double* afrc, double* bfrc,
			 Energy *sysUV, double elec14fac, double lj14fac,
			 int qform);

void AttenuatePairFrc(prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc, int atmA,
		      int atmB, double dx, double dy, double dz, double fmag,
		      double* afrc, double* bfrc, double elec14fac,
		      double lj14fac, int qform);

void AttenuatePairNrg(prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc, int atmA,
		      int atmB, double dx, double dy, double dz,
		      Energy *sysUV, double elec14fac, double lj14fac,
		      int qform);

void BondVir(double* aptr, double *bptr, double* afrc, double* bfrc,
             bondcomm *bcom, prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc,
	     Energy *sysUV, int qform);

void BondFrcNrg(double* aptr, double *bptr, double* afrc, double* bfrc,
                bondcomm *bcom, prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc,
		Energy *sysUV, int qform);

void BondFrc(double* aptr, double *bptr, double* afrc, double* bfrc,
             bondcomm *bcom, prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc,
	     int qform);

void BondNrg(double* aptr, double *bptr, bondcomm *bcom, prmtop *tp,
             FrcTab *Cfrc, FrcTab *Hfrc, Energy *sysUV, int qform);

void AnglVir(double *aptr, double *bptr, double *ctpr, double *afrc,
             double *bfrc, double *cfrc, anglcomm *acom, prmtop *tp,
	     FrcTab *Cfrc, FrcTab *Hfrc, Energy *sysUV, int qform);

void AnglFrcNrg(double *aptr, double *bptr, double *ctpr, double *afrc,
		double *bfrc, double *cfrc, anglcomm *acom, prmtop *tp,
		FrcTab *Cfrc, FrcTab *Hfrc, Energy *sysUV, int qform);

void AnglFrc(double *aptr, double *bptr, double *ctpr, double *afrc,
	     double *bfrc, double *cfrc, anglcomm *acom, prmtop *tp,
	     FrcTab *Cfrc, FrcTab *Hfrc, int qform);

void AnglNrg(double *aptr, double *bptr, double *ctpr, anglcomm *acom,
	     prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc, Energy *sysUV, int qform);

void DiheVir(double *aptr, double *bptr, double *cptr, double *dptr,
             double *afrc, double *bfrc, double *cfrc, double *dfrc,
             dihecomm *hcom, prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc,
	     Energy *sysUV, int qform);

void DiheFrcNrg(double *aptr, double *bptr, double *cptr, double *dptr,
		double *afrc, double *bfrc, double *cfrc, double *dfrc,
		dihecomm *hcom, prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc,
		Energy *sysUV, int qform);

void DiheFrc(double *aptr, double *bptr, double *cptr, double *dptr,
	     double *afrc, double *bfrc, double *cfrc, double *dfrc,
	     dihecomm *hcom, prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc,
	     int qform);

void DiheNrg(double *aptr, double *bptr, double *cptr, double *dptr,
	     dihecomm *hcom, prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc,
	     Energy *sysUV, int qform);

void CmapVir(double *aptr, double *bptr, double *cptr, double *dptr,
             double *eptr, double *afrc, double *bfrc, double *cfrc,
	     double *dfrc, double *efrc, cmapcomm *mcom, prmtop *tp,
	     Energy *sysUV, int qform);

void CmapFrcNrg(double *aptr, double *bptr, double *cptr, double *dptr,
		double *eptr, double *afrc, double *bfrc, double *cfrc,
		double *dfrc, double *efrc, cmapcomm *mcom, prmtop *tp,
		Energy *sysUV, int qform);

void CmapFrc(double *aptr, double *bptr, double *cptr, double *dptr,
	     double *eptr, double *afrc, double *bfrc, double *cfrc,
	     double *dfrc, double *efrc, cmapcomm *mcom, prmtop *tp,
	     int qform);

void CmapNrg(double *aptr, double *bptr, double *cptr, double *dptr,
	     double *eptr, cmapcomm *mcom, prmtop *tp, Energy *sysUV,
	     int qform);
#endif
