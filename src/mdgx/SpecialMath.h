#ifndef SpecialMathHeadings
#define SpecialMathHeadings

#include "CompFrcDS.h"

double perfc(double x);

double PerfcFast(double x);

double SinFast(double x, CSpln* stab);

double CosFast(double x, CSpln* ctab);

double AsinFast(double x, CSpln* astab);

double AcosFast(double x, CSpln* actab);

double SafeSqrt(double x);

#endif
