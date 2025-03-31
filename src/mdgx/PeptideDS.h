#ifndef PeptideDataStructures
#define PeptideDataStructures

#include "MatrixDS.h"
#include "GpuDS.h"
#include "AmberNetcdfDS.h"

struct SharedMemoryHostMockup {
  int*   ixfrc;     // Fixed-precision integer representations of forces in X,
  int*   iyfrc;     // Y, and Z dimensions.  These representations will break
  int*   izfrc;     // at +/- 4096 kcal/mol-A, but that should be enough.
  int*   ixbuff;    // Buffers for perturbations that will be contributed
  int*   iybuff;    //   to ixfrc, iyfrc, and izfrc at the end of each batch,
  int*   izbuff;    //   replicating the process on the GPU
  int*   active;    // Array to show which atoms are affected by any set of 32
                    //   constraint operations
  int*   ljid;      // Lennard-Jones type IDs for all atoms
  int*   neckid;    // Indices into the neck GB parameters array
  float* xcrd;      //
  float* ycrd;      // Atom X, Y, and Z coordinates
  float* zcrd;      //
  float* xprvcrd;   //
  float* yprvcrd;   // Previous X, Y, and Z coordinates of each atom
  float* zprvcrd;   //
  float* qval;      // Atomic partial charges
  float* rborn;     // Baseline Born radii
  float* reff;      // Effective Born radii
  float* sumdeijda; // Generalized Born derivative accumulators
  float* gbalpha;   // Generalized Born alpha parameters for each atom
  float* gbbeta;    // Generalized Born beta parameters for each atom
  float* gbgamma;   // Generalized Born gamma parameters for each atom
  float* fs;        // Screening parameters for all atoms in Generalized Born
  int*   psi;       // Pre-computations for computing effective Born radii
};
typedef struct SharedMemoryHostMockup slmem;

struct PeptideControl {
  int nsys;           // The number of simulations to manage simultaneously
  int remd;           // Flag to indicate that remd is to be performed
  int ntexch;         // The number of steps between exchanges (must be a
                      //   factor of ntwx, if not mdgx will set it to the
                      //   nearest such factor, rounding up)
  int nsegment;       // The number of dynamics segments to perform: each
                      //   segment is a number of MD steps.  The number of
                      //   steps per segment must be a multiple of ntpr and
                      //   equal to ntwx from the &cntrl namelist.
  int nsgmdout;       // The number of times diagnostics (mdout) will be
                      //   printed per segment
  int notraj;         // Flag to have no trajectory output at all (triggered
                      //   by ntwx = 0, as a nonzero ntwx must be present for
                      //   managing the dynamics segments)
  int nBondSteps;     // The number of steps taken to update bonded
                      //   interactions in between the major (non-bonded and
                      //   dihedral) force updates.  Angles are updated at
                      //   half this frequency.  
  int prngRefresh;    // The pseudo-random number refresh rate (GPU only)
  int totalPRNG;      // The total number of pseudo-random numbers to generate
                      //   each time the heap is refreshed
  int rattle;         // Flag to engage RATTLE for all systems
  int igb;            // The type of Generalized Born calculations to include
  int blockDim;       // The block dimension to use (this will determine the
                      //   choice of kernel, one size does not fit all
  int blocks;         // The number of blocks to launch (this will vary
                      //   depending on the architecture)
  double dielectric;  // Solvent dielectric constant
  double gboffset;    // Genralized Born radius offset parameter
  int* Treplicas;     // The numbers of temperature replicas to make for each
                      //   system, at temperatures on a uniform grid between
                      //   values in Tranges.x and Tranges.y
  int* Preplicas;     // The numbers of topology replicas to make for each
                      //   system, interpolating between the -pi and -pf
                      //   topologies
  int* Sreplicas;     // The numbers of standard replicas to make for each
                      //   system.  If specified, this will kill REMD and just
                      //   make lots of copies running equilibrium MD at the
                      //   one temperature Tmin or temp0 from &cntrl.
  double* T;          // Temperature settings for each instance
  double* Tmix;       // Mixing factors for each instance, between temperatures
  double* Pmix;       // Mixing factors for each instance, between topologies
  double* starttime;  // Starting times for each replica (read from restart
                      //   files if available)
  double2* Tranges;   // Temperature ranges for simulations on each system
                      //   (x field = low value, y field = high value)
  cdftrj* CDFHandles; // Pointers to NetCDF structs (trajectory output)
  char blockreq[64];  // User-specified block size requirement, options are
                      //   "SMALL", "MEDIUM", and "LARGE." This directive will
                      //   only take effect if the selection is viable for the
                      //   systems at hand.
  cmat tpinames;      // Paths to the initial topologies (these will be
                      //   interpolated to the "final" topologies if many
                      //   simulations rely on the same pair of initial and
                      //   final topologies)
  cmat tpfnames;      // Paths to the final topologies (these will be
                      //   interpolated to the "final" topologies if many
                      //   simulations rely on the same pair of initial and
                      //   final topologies)
  cmat crdnames;      // Paths to the initial coordinates for each simulation
  cmat outbases;      // Base names for output files from each simulation
                      //   (a suffix for these output files from the &files
                      //   namelist will be applied)
  cmat trjbases;      // Base names for trajectory output files (a suffix for
                      //   these output files from the &files namelist will be
                      //   applied)
  cmat rstrtbases;    // Base names for restart (checkpoint) files for each run
};
typedef struct PeptideControl pepcon;

#endif
