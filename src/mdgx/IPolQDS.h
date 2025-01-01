#ifndef IPOLQ_STRUCTS
#define IPOLQ_STRUCTS

#ifndef PREP_API
#include "MatrixDS.h"
#include "CrdManipDS.h"
#include "BSplineDS.h"
#endif

struct IPolQControlData {
  int ntqs;              // The frequency of simulation sampling
  int nqframe;           // Number of frames of simulation used to construct
                         //   the charge density
  int nqmod;             // The number of charge modifications to implement
  int nQshell;           // The number of charge shells used to correct the
                         //   electrostatics of explicit charges
  int nVshell;           // The number of charge shells used to sample the
                         //   electrostatics at and around solute atom sites
  int neqstep;           // The number of steps to be used in equilibration,
                         //   before any data collection
  int nQphpt;            // The number of equidistant points on the sphere
                         //   surface used to place charges for reproducing the
                         //   the solvent reaction field potential
  int nVphpt;            // The number of equidistant points on the sphere
                         //   surface used to sample the solvent reaction field
                         //   potential near the solute atom sites
  int nblock;            // Counter to direct accumulation of block averages
  int nQcloud;           // The number of explicit charges in the charge cloud
                         //   generating the solvent reaction field
  int stotq;             // Total solute charge (convenient just to store this)
  int excitation;        // Excitation state to target in the QM calculations
                         //   (default zero, may not work the same in all QM
                         //    packages)
  int verbose;           // Output level (default 0, no special messages)
  int retqminp;          // Flag to retain QM input
  int retqmchk;          // Flag to retain QM checkpoint file
  int retqmout;          // Flag to retain QM output
  int retptfi;           // Flag to retain charge cloud
  int checkex;           // Flag to activate existence checks on QM executables
                         //   (default 1, ON)
  int CenterGrid;        // Indicates the type of grid centering to perform
                         //   The default behavior depends on the specific
                         //   quantum program, but user input can override the
                         //   default.
  int MaxCore;           // The maximum memory to be allocated for arrays in
                         //   QM calculations, in units of MB (default 512)
  int nQMThreads;        // The number of threads to be used by the QM program
  int DoSinglePoint;     // Flag to submit the conformation for single point
                         //   energy and gradient evaluation.  This will return
                         //   the average of the energy of the polarized wave
                         //   function in a vacuum plus the energy of the
                         //   unpolarized wave function in a vacuum.  The
                         //   gradients are computed likewise.
  int gdim[3];           // Electrostatic potential grid dimensions
  int* qnbrs;            // Array storing atoms that participate in the
                         //   explicit charge pool
  double econv;          // The convergence criterion for the electrostatic
                         //   potential
  double minqfac;        // Factor for tethering fitted shell charges to zero
  double Qshell[4];      // Up to three charge shells will be placed at these
                         //   distances from the solute.  Charges sitting at
                         //   less than qshell[0] from the solute are included
                         //   explicitly in the calculation
  double Vshell[4];      // Up to three shells will be placed at these
                         //   distances from the solute atom sites.
  double gspc[3];        // Electrostatic potential grid spacings
  double gorig[3];       // Origin of the electrostatic potential grid
  double* Vfrc;          // Vector storing the running averages of forces at
                         //   points of Vsurf
  double* QModVal;       // Vector storing temporary charge modifications
  char* SoluteMask;      // A description of the solute, which will be kept
                         //   immobile; the solvent is defined reciprocally
  char* qmprog;          // The quantum mechanics program called by mdgx
  char* qmpath;          // Path to the QM executable
  char* uvpath;          // Path to the potential evaluator executable
  char* inpfile;         // Quantum mechanics input file
  char* outfile;         // Quantum mechanics input file
  char* ptqfile;         // Point charge file for solvated quantum calculation
  char* finfile;         // File written at the end of quantum calculations
  char* grdfile;         // Grid file (Gaussian format) of electrostatic U
  char* scrdir;          // Scratch directory in which to do the quantum work
  char* qmmeth;          // Quantum mechanics theory level (default mp2)
  char* basis;           // Quantum mechanics basis set
  char* fmpath;          // Checkpoint file formatting program (for Gaussian)
  coord Solute;          // Coordinates of the solute molecule
  coord Qsurf;           // A set of equidistant points around the solute where
                         //   charges are placed to reproduce electrostatics
  coord Vsurf;           // A set of equidistant points around the solute where
                         //   the electrostatic potential and gradient are to
                         //   be calculated
  dmat Qcloud;           // The explicit charges from the simulation that will
                         //   be included in quantum calculations
  dmat SAfrc;            // Matrix storing forces due to the solvent reaction
                         //   field at solute atom sites from all frames
                         //   sampled
  cmat prepcalls;        // System calls to be made in preparation for QM
                         //   calculations
  cmat postcalls;        // System calls to be made after QM calculations
  cmat QModMask;         // AmbMask strings for charge modification
};
typedef struct IPolQControlData ipqcon;

#endif
