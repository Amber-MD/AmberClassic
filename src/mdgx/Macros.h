#ifndef MacroDefinitions
#define MacroDefinitions

//-----------------------------------------------------------------------------
// Useful functions
//-----------------------------------------------------------------------------
#define SIGN(x)            ((x >= 0.0) ? 1.0 : -1.0)

#define SIGN2(a,b)         ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define MAX(a, b)          ((a) > (b) ? (a) : (b))

#define MIN(a, b)          ((a) > (b) ? (b) : (a))

#define SWAP(x, y, tmp)    (tmp) = (x); (x) = (y); (y) = (tmp)

#endif
