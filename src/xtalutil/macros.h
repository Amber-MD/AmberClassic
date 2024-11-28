#ifndef SIGN
#define SIGN(x)            ((x >= 0.0) ? 1.0 : -1.0)
#endif

#ifndef SIGN2
#define SIGN2(a,b)         ((b) >= 0.0 ? fabs(a) : -fabs(a))
#endif

#ifndef MAX
#define MAX(a, b)          ((a) > (b) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a, b)          ((a) > (b) ? (b) : (a))
#endif

#ifndef SWAP
#define SWAP(x, y, tmp)    (tmp) = (x); (x) = (y); (y) = (tmp)
#endif

#ifndef DIST
#define DIST(x, y, z)      (sqrt(((x)*(x) + (y)*(y) + (z)*(z))))
#endif

#ifndef SQ_DIST
#define SQ_DIST(x, y, z)   ((x)*(x) + (y)*(y) + (z)*(z))
#endif

#ifndef NOTEQ
#define NOTEQ(x, y)        (fabs((x) - (y)) > MACH_PREC_ZERO)
#endif

#ifndef DNEQ
#define DNEQ(x, y)         ((fabs(x) > 1.0 && fabs(y) > 1.0 && fabs((x) - (y))/fabs(0.5*((x) + (y))) > 1.0e-12) || fabs((x) - (y)) > 1.0e-12)
#endif
