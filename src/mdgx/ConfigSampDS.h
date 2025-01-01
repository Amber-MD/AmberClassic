#ifndef ConfigStructs
#define ConfigStructs

#include "MatrixDS.h"

//-----------------------------------------------------------------------------
// CoordinateOperation: an operation to perform on coordinates in order to
//                      generate a new configuration.  Each operation is
//                      essentially the application of a harmonic restraint to
//                      coax the relative positions of two, three, or four
//                      atoms into a particular orientation.  Many operations
//                      carries with their own range, however: the idea is to
//                      apply this operation to many copies of a molecule, and
//                      in each particular case to pick a value within the
//                      given range and apply a custom harmonic restraint to
//                      the given atoms that moves the one copy of the system
//                      towards that value.
//-----------------------------------------------------------------------------
struct CoordinateOperation {
  int order;        // Order of the operation: 2 is bond, 3 angle, 4 dihedral.
  int absrange;     // Flag to have values selected from the range [min, max]
                    //   irrespective of what geometry the coordinates
                    //   initially have for the atoms of interest, if set to 1.
                    //   Otherwise values will be selected from the range
                    //   [ X-min, X+max ], where X is the initial value found
                    //   in the coordinates.
  int scatter;      // Flag to perform random sampling on a uniform
                    //   distribution within the specified range, if set to 1.
                    //   Otherwise uniform sampling will be performed in the
                    //   same range.
  int anchors[4];   // The anchor atoms for this operation (up to four may be
                    //   specified, but the number is the same as the order)
  int pegtype;      // Flag to have the one atom of this coordinate operation
                    //   iteratively pegged to points in space that scan a box
                    //   volume defined by min(x,y,z) and max(x,y,z)
                    //   (pegtype == 0) or in a line between the min and max
                    //   coordinate (pegtype == 1) 
  double Krst;      // Restraint force--units vary by the operation's order:
                    //   kcal/mol-A^2 for bonds, kcal/mol-rad^2 otherwise.
  double quadwin;   // The length of the quadratic potential window (implied
                    //   by Krst and Ftop--either value can be set by the user)
  double fbhw;      // Half Width of the Flat Bottom part of the restraint
                    //   (set this to zero and Ftop to something very high
                    //   for purely harmonic restraints)
  double distance;  // The distance at which to hold an atom from its target
                    //   position in the context of a positional restraint
                    //   (default 0.0)
  double minval;    // Minimum (or the sole value) of the sampling range
  double maxval;    // Maximum value of the sampling range  
  double minvx;     //
  double minvy;     // Minimum and maximum values of a coordinate range, used
  double minvz;     //   to specify points in space to which a particular atom
  double maxvx;     //   shall be bound.
  double maxvy;     //
  double maxvz;     //
  char label[32];   // A label to give this operation (optional, but lets it
                    //   be tracked later)
  cmat atommasks;   // Ambmask strings to denote atoms in this operation
  double* targets;  // Target equilibria for this coordinate operation (these
                    //   may be computed after getting a look at the actual
                    //   structures)
  dmat refcrd;      // The reference coordinates to use in the case of a
                    //   positional restraint
};
typedef struct CoordinateOperation cfigop;

//-----------------------------------------------------------------------------
// CoordinateTrack: a track restraint is sort of a generalized point position
//                  restraint.
//-----------------------------------------------------------------------------
struct CoordinateTrack {
  int nknot;          // The number of knots in the path
  int ntp;            // The number of track particles
  int atmA;           // First atom to join along the path
  int atmB;           // Second atom to join along the path
  int* monotonic;     // Array of flags to indicate whether the path is
  double* crd;        // Coordinates of each knot in the path
  double* xcoef;      // Coeffficients of quadratic splies connecting the knots
  double* ycoef;      //   in x, y, and z, in the order A1, B1, C1, A2, B2, C2,
  double* zcoef;      //   ... for splines between knots 1 and 2, 2 and 3, ...
  double p1;          // Location of the first particle along the path
  double p2;          // Location of the second particle along the path
  double fbhw;        // Flat bottom half width of the piecewise harmonic
                      //   restraint pulling the track particles together along
                      //   the path
  double l0;          // Equilibrium length of the piecewise harmonic restraint
                      //   pulling the two particles together along the path
  double quadwin;     // Distance for which the retraint between p1 and p2 is
                      //   harmonic (beyond this it is a linear potential)
};
typedef struct CoordinateTrack track;

//-----------------------------------------------------------------------------
// QuantumSettings: a structure to hold the quantum program inputs for
//                  configuration outputs.
//-----------------------------------------------------------------------------
struct QuantumSettings {
  int spin;            // Spin multiplicity of the molecule
  int MaxCore;         // The maximum memory (per core) to be allocated by the
                       //   quantum program
  int ncpu;            // The number of CPUs to run the QM program on
                       //   (parallel MPI setting)
  int opt;             // Flag to have the quantum calculations also set up for
                       //   geometry optimization (default 0, false)
  int xferRst;         // Flag to have restraints stipulated in the mdgx input
                       //   transferred to the quantum program for optimization
                       //   (default 1, true, and will thus be invoked if opt
                       //   is set to true)
  char theory[64];     // Level of QM theory (e.g. MP2)
  char opttheory[64];  // Level of QM theory for optimization 
  char basis[64];      // Basis set to use (e.g. 6-31++g)
  char optbasis[64];   // Basis set to use in optimization
  char* checkpoint;    // Checkpoint file for QM calculations
};
typedef struct QuantumSettings qmset;

//-----------------------------------------------------------------------------
// AtomExclusions: a structure for storing the exclusions of an atom in a
//                 manner that is highly efficient in all-to-all calculations.
//-----------------------------------------------------------------------------
struct AtomExclusions {
  int idx;          // The index of this atom in the master topology
  int llim;         // The low and high limits of excluded atom indices.  These
  int hlim;         //   values give a range over which a more thorough check
                    //   will be necessary, assuming that in most cases the
                    //   system's bonds extend between atoms that are nearby in
                    //   the order.
  int nexcl;        // The number of exclusions for this atom       
  int* nblist;      // The list of non-bonded exclusions for atom idx
  double* qqval;    // Strengths of the electrostatic and Lennard-Jones
  double* ljval;    //   non-bonded interactions (these will never be 1.0, as
                    //   only interactions with some degree of attenuation will
                    //   be counted among the list, but they may be greater
                    //   than 0.0)
};
typedef struct AtomExclusions excltab;

//-----------------------------------------------------------------------------
// ConvergenceStatistic: a record of convergence information on a particular
//                       configuration's energy minimization.  An array of
//                       these beats having to collect all types of the
//                       information into separate vectors.  This will also
//                       hold information on whether the configuration passes
//                       various sanity checks.
//-----------------------------------------------------------------------------
struct ConvergenceStatistics {
  int nstep;        // Number of steps of line minimization needed to
                    //   converge--a value of -1 indicates this did not
                    //   converge
  int maxFatom;     // Atom upon which the maximum force for the configuration
                    //   is acting
  int pass;         // Flag to indicate whether this configuration passes the
                    //   sanity checks after energy minimization
  double lastE;     // The most recently recorded (non-restraint) energy for
                    //   this configuration
  double lastER;    // The most recently recorded restraint energy for
                    //   this configuration
  double lastdE;    // Change in configuration (non-restraint) energy for the
                    //   last step
  double lastdER;   // Change in configuration restraint energy for the last
                    //   step
  double maxF;      // Largest force on any atom of this configuration in the
                    //   final state (this is a magnitude, not a vector)
  int bstrnloc[3];  // Index of the atom (first element) and index of the bond
                    //   that atom controls (second element) which identify the
                    //   bond most strained in this configuration.  The third
                    //   element is a fixed-precision, approximate value of the
                    //   offending bond.
  int astrnloc[3];  // Index of the atom (first element) and index of the angle
                    //   that atom controls (second element) which identify the
                    //   angle most strained in this configuration.  The third
                    //   element is a fixed-precision, approximate value of the
                    //   offending angle.
  int maxrloc;      // Index of the highest-energy restraint found in this
                    //   configuration
  double maxbstrn;  // The maximum bond strain found in this configuration
  double maxastrn;  // The maximum angle strain found in this configuration
  double maxrst;    // The maximum restraint energy found in this configuration
  double dEshuffle; // Cumulative changes in energy due to shuffling and
                    //   reoptimizing
};
typedef struct ConvergenceStatistics ConvStat;

//-----------------------------------------------------------------------------
// OperationCoupler: a struct for detailing that two or three operations shall
//                   be combined for grid-based sampling.
//-----------------------------------------------------------------------------
struct OperationCoupler {
  int npcs;           // The number of pieces (two or three)
  int pcs[3];         // The pieces
  char labels[3][32]; // A two-dimensional array of strings for the labels
                      //   of operations that this will couple
};
typedef struct OperationCoupler coupler;

//-----------------------------------------------------------------------------
// NonbondedBox: a struct for storing the local information needed to compute
//               the non-bonded energy of movable atoms in the presence of the
//               molecule's rigid components.
//-----------------------------------------------------------------------------
struct NonbondedBox {
  int ngrdx;
  int ngrdy;
  int ngrdz;
  double gspc;
  double origx;
  double origy;
  double origz;
  int natom;
  int* atmidx;
  double* crd;
  double* elecU;
  double* vdwU;
};
typedef struct NonbondedBox nbbox;

//-----------------------------------------------------------------------------
// NonbondedStack: a struct for storing an intricate grid of potentials and
//                 atom information needed to approximate (to a very high
//                 precision) the nonbonded energy of a molecule with a
//                 significant rigid component.
//-----------------------------------------------------------------------------
struct NonbondedStack {
  int nbinx;
  int nbiny;
  int nbinz;
  double gspc;
};
typedef struct NonbondedStack nbstack;

//-----------------------------------------------------------------------------
// ConfigurationSearch: a data structure to store directives for searching
//                      different configurations / conformations of a molecule.
//-----------------------------------------------------------------------------
struct ConfigurationSearch {
  int count;         // The number of configurations ultimately to generate
  int nops;          // The length of ops; the number of operations to perform
  int ncombo;        // The number of combinations binding operations together
  int ncyc;          // Number of steepest descent steps for energy
                     //   minimization
  int maxcyc;        // Maximum number of steps (cycles) for energy
                     //   minimization
  int nbelim;        // Number of configurations eliminated (not printed)
                     //   because they failed to meet the bond strain criterion
  int naelim;        // Number of configurations eliminated (not printed)
                     //   because they failed to meet the angle strain
                     //   criterion
  int nrelim;        // Number of configurations eliminated (not printed)
                     //   because they failed to meet the overall restraint
                     //   energy criterion
  int npass;         // The number of configurations passing all basic sanity
                     //   checks
  int verbose;       // Verbosity level (0, default, is silent, 1 will generate
                     //   information on the progress of the run
  int freezeH;       // Flag to freeze the lengths of bonds to hydrogen atoms
                     //   at their equilibrium values
  int showorigins;   // Flag to have mdgx print a list of origins for each
                     //   configuration it generates
  int nbelly;        // The number of bellymasks provided (if zero, the default
                     //   is to have all atoms be mobile)
  int allmove;       // Flag to indicate that all atoms are indeed mobile
  int atomlimit;     // The maximum number of atoms that the system may have
                     //   before mdgx starts working off of an exclusions list
                     //   rather than a scaling matrix for non-bonded
                     //   interactions (default 512)
  int rattle;        // Flag to have bonds to hydrogen constrained by RATTLE
  int MaxRattleIter; // Maximum number of iterations for RATTLE (exceeding this
                     //   will trigger a statistic, not an error, as this is
                     //   about minimization)
  int FreezePolarH;  // Flag to freeze rotation about bonds connecting to a
                     //   polar hydrogen
  int PreventHBonds; // Flag to have hydrogen bonding generally prevented via
                     //   automated restraint generation.  If activated (it is
                     //   set to 0, OFF by default), it will identify all
                     //   polar hydrogens in the system (using atomic numbers)
                     //   and lone pair donors, then establish a set of
                     //   distance restraints between all members of each set.
  int SpinHydroxyls; // Flag to have polar hydrogens randomly rotated if they
                     //    are part of hydroxyl groups.  This will reorient
                     //    dipoles if sampling electrostatics is important.
  int* movable;      // Flags for all atoms in the system to tell whether they
                     //   are moveable (set to 1) or frozen (set to 0)
  double fconv;      // Convergence criterion for energy minimization: when the
                     //   net force on each atom (including restraints) is less
                     //   than this number the minimization for a configuration
                     //   will be considered converged.
  double stepconv;   // Convergence criterion for energy minimization: the step
                     //   length (defined as the absolute magnitude of the
                     //   movement of all atoms in the configuration) goes
                     //   below this value and we declare the energy
                     //   minimization to be complete.
  double step0;      // Initial step size (default 0.01)
  double rattletol;  // Convergence goal for RATTLE algorithm
  double strainlim;  // The limit of restraint energy that is permitted in each
                     //   configuration; structures which violate this maximum
                     //   will not be printed.
  double HBMargin;   // The margin for preventing hydrogen bonding
  double KSpinOH;    // The stiffness of restraint to impose when spinning
                     //   hydroxyl groups
  double fbhwSpinOH; // Flat bottom half-width to apply when spinning hydroxyl
                     //   groups
  double maxbstrn;   // The maximum bond strain that will be permitted for any
                     //   particular bond in any configuration
  double maxastrn;   // The maximum angle strain that will be permitted for any
                     //   particular angle in any configuration
  double simEtol;    // The minimum mutual energy differences that structures
                     //   in the final list must have, unless they have already
                     //   passed the mutual rmsd criterion
  double rmsdtol;    // The minimum mutual rmsd that structures in the final
                     //   list must have
  cmat ostyle;       // Styles of the output (available options are enumerated
                     //   in the Command library when input is read).  Multiple
                     //   styles may be specified to have the configurations in
                     //   many forms.
  cmat outbase;      // Base names of the output files (configurations, either
                     //   as coordinates or quantum calculation input files).
                     //   One base name may be specified for every output
                     //   style, or multiple base names may be given for each
                     //   style.
  cmat outsuff;      // Suffix of the output file names.  As with the base
                     //   name, one suffix may be given, or many may be
                     //   specified.
  cmat origins;      // A list of the origins of each configuration.  When
                     //   processing a list of files in a directory, the
                     //   opendir function tends to jumble things up as a
                     //   human would see it.  This makes it possible to trace
                     //   things back to where they began.
  cmat belly;        // Bellymask strings to hold atoms fixed during energy
                     //   optimization
  cfigop* ops;       // A series of operations to perform on each configuration
                     //   to differentiate the various copies of the
                     //   coordinates or alter them all relative to their
                     //   original values
  coupler* combo;    // Combinations of operations to perform multi-dimensional
                     //   sampling
  qmset QMsettings;  // Inputs to pass on to the quantum program

  // Attributes to control how configurations can be reshuffled to try and
  // achieve lower energy states
  int reshuffle;            // The number of reshuffling rounds to perform
  int nreopt;               // The number of reoptimizations performed as a
                            //   result of reshuffling the initial states
  double proximity;         // Proximity required of replacement candidates (in
                            //   terms of restraint energy)
  double Ereplace;          // Threshold of total energy improvement (including
                            //   restraints) at which a reoptimized
                            //   configuration will be admitted (note that
  char shfdir[32];          // Direction in which to take the shuffling (can be
                            //   "up" or "down", to amke shuffling and
                            //   reoptimization take the energies of solutions
                            //   for each set of restraints up or down)
  char shuffletype[32];     // The type of shuffling to perform (available
                            //   options include "jackknife", "bootstrap", and
                            //   "proximity."
};
typedef struct ConfigurationSearch configs;

#endif
