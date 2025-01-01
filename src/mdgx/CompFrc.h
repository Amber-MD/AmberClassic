#ifndef CompFrcHeadings
#define CompFrcHeadings

#include "CompFrcDS.h"

CSpln* SplineTab(double* u, double* du, int nbin, double dr);

CSpln* SplineTab2(double* u, double* du, int nbin, double dr);

FrcTab DirectSpaceR2(double range, double spc, double ewcoeff, int contderiv,
		     int doerf);

void FreeFrcTab(FrcTab *F);
#endif
