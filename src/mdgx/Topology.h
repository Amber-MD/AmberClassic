#ifndef TopologyHeadings
#define TopologyHeadings

#include "Constants.h"
#include "TopologyDS.h"
#include "CrdManipDS.h"
#include "pmeDirectDS.h"
#include "TrajectoryDS.h"

cmat Ascii2Mem(char* fname, int maxw, int stride, char* errmsg);

void OrderLJParameters(prmtop *tp);

void ComputeLinearCmapDerivatives(dmat *esrf, dmat *dphi, dmat *dpsi,
                                  dmat *d2pp);

dmat Compute1DStencil(int npts, double* tvec);

void SetBasicTopology(prmtop *tp, char* dsource, char* epsource, int dorattle,
                      int dosettle);

void GetPrmTop(prmtop *tp, trajcon *tj, int adjbnd);

prmtop CopyTopology(prmtop *tp);

prmtop InterpolateTopology(prmtop *tpA, prmtop *tpB, double lambda);

void FreeTopology(prmtop *tp);

int FindAtom(prmtop *tp, int il, int ih, char* aname);

void PutPrmTop(prmtop *tp, char* fname, char* title);

int ProteinChiralityCheck(prmtop *tp, coord *tc, FILE *outp);

int FindDisulfides(prmtop *tp, coord *tc, FILE *outp);

int CheckLJRules(prmtop *tp, FILE *outp);

int LocateResID(prmtop *tp, int ai, int llim, int hlim);

void LongRangeVDW(prmtop *tp, dircon *dcinp);

#endif
