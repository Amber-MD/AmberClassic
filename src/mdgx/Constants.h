#ifndef MdgxConstants
#define MdgxConstants

// Pi and some useful fractions
#define PI        3.1415926535898
#define TWOPI     6.2831853071796
#define HALFPI    1.5707963267949
#define INVPI     3.1830988618379e-1
#define INV_TWOPI 1.5915494309189e-1
#define ONE_THIRD 0.3333333333333
#define ONE_SIXTH 0.1666666666667
#define ELE_SIXTH 1.8333333333333
#define TN_BASE_E 2.3025850929940

// The tetrahedral angle, in radians
#define TETANG    1.9106332356470

// Base of the natural logarithm
#define LNBASE    2.7182818284590

// Discretizations of lookup tables
#define SIN_TAB_DSC  6.2831853071796e-2
#define INV_SIN_DSC  1.5915494309189e+1
#define ASIN_TAB_DSC 0.02
#define ACOS_TAB_BUF 0.76
#define INV_ASIN_DSC 50.0
#define COS_OFFSET   25

// Coulomb's constant in units of kcal/mol-e^2
#define BIOQ 332.0522

// The universal gas constant in kcal/mol-K
#define GASCNST 1.987216e-3

// Avogadro's number
#define AVOGADRO 6.0221418e+23

// Factor for converting pressure in kcal/mol-A^2 to bar
#define PCONVFAC 6.94770014e4

// Conversion of Hartrees to kcal/mol
#define H2KCAL 627.509649

// Conversion of e-A to Debye
#define EA2DEBYE 4.80245625

// Conversion of Bohr to Angstrom
#define B2ANG 0.529177249

// The nonbonded minimum distance
#define MINNB  2.0
#define MINNB2 4.0

// Factor for internal velocities
#define VEL_INT_SCALING 20.4548282809

// Double-to-integer conversion factors
#define LOC2IFAC   8388608.0
#define I2LOCFAC   1.1920928955078e-07
#define I2LOCBITS  23
#define Q2IFAC     268435456.0
#define I2QFAC     3.7252902984619e-09
#define I2QBITS    28

// String length variables
#define MAXTITL 128
#define MAXNAME 512
#define MAXLINE 1024

// Maximum number of systems for replica exchange
#define MAXSYS 512

// Detect debugging flags
#ifdef DEBUG
#define IDEBUG    1
#else
#define IDEBUG    0
#endif

// Maximum valgrind-acceptable message size
#ifdef MPI
#define VALGR_MAX_BUFF   64000
#endif

#endif
