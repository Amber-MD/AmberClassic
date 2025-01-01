#ifndef CompFrcStructs
#define CompFrcStructs

//-----------------------------------------------------------------------------
// CubicSpline: this structure holds coefficients for a piecewise cubic spline
//              approximation of some real-valued function of one variable.
//              An array of these structs forms the complete approximation
//              over some range of the variable:          
//                                                                      
//              f(x) = Ax^3 + Bx^2 + Cx + D                             
//-----------------------------------------------------------------------------
struct CubicSpline {
  double A;
  double B;
  double C;
  double D;
};
typedef struct CubicSpline CSpln;

//-----------------------------------------------------------------------------
// ForceTable: this structure holds a force (and energy) piecewise cubic
//             spline table.                                            
//-----------------------------------------------------------------------------
struct ForceTable {

  // Attributes needed for force and energy evaluation
  int nbin;
  double rmax;
  double dr;
  double ivdr;
  CSpln* SD;
  CSpln* dSD;

  // Error analysis
  int FitType;            // The type of spline interpolation
                          //   0: best fit at intervals 0, 1/3, 2/3, 1, ...
                          //   1: fit with continuous derivatives at 0, 1, ...
  double fmaxerr;         // Maximum absolute error in the derivative (force)
  double fmaxerrloc;      // Location of the maximum error in the derivative
  double fmaxrelerr;      // Maximum relative error in the derivative
  double fmaxrelerrloc;   // Location of maximum relative error...
  double umaxerr;         // Maximum absolute error in the function (potential)
  double umaxerrloc;      // Location of the maximum error in the function
  double umaxrelerr;      // Maximum relative error in the function
  double umaxrelerrloc;   // Location of maximum relative error in the function
};
typedef struct ForceTable FrcTab;

#endif
