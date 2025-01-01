#ifndef RestraintStructs
#define RestraintStructs

#ifndef PREP_API
#include "MatrixDS.h"
#include "GridDS.h"
#endif

struct ElectronDensityGaussian {
  double amp;       // Amplitude of the Gaussian
  double sig;       // Gaussian standard deviation (width)
  double xtns;      // The extensivity of the Gaussian (in theory it should
                    //   extend forever, but sometimes these 
};
typedef struct ElectronDensityGaussian edgauss;

struct SymmetryOperation {
  int nt;           // The number of symmetry operations encoded in this set
  dmat* rmat;       // Rotation matrices for all symmetry operations
  double* tvec;     // Translation vectors for all symmetry operations
};
typedef struct SymmetryOperation symop;

struct ClingRule {
  int ngss;         // The number of Gaussians that describe this atom
  edgauss* gss;     // gaussian electron density masks
};
typedef struct ClingRule cling;

struct RestraintControls {
  int active;          // Flag to indicate that any restraints are active
  int usegrid;         // Flag to indicate a grid-based restraint is in use
  int usebelly;        // Flag to activate a belly mask and freeze all other
                       //   atoms
  int XpandGrid;       // Number to indicate what grid promotion should be used
                       //   in order to get smoother lookup tables (default of
                       //   1 means no expansion, expand by a factor of 1)
  double GridScale;    // Scaling factor to make the grid on file uniformly
                       //   stronger or weaker at run time
  char* GridFile;      // Name of the grid restraint file
  char* GridDefsFile;  // Name of file describing how atoms of the simulation
                       //   interact with the grid
  char* BellyMask;     // The belly mask string (ambmask format)
  char* FrozenMask;    // The belly mask string (ambmask format)
  fbook Rgrd;          // The restraint grid
};
typedef struct RestraintControls rstrcon;

#endif
