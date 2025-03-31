#ifndef SinglePointEvalDataStructures
#define SinglePointEvalDataStructures

struct QuantumResults {
  int pass;          // Flag to indicate that this is a valid QM result
  int natom;         // The number of atoms involved
  int nbfunc;        // The number of basis functions involved
  int haserror;      // Flag to indicate that errors in the output were found
  int ccremoved;     // Flag to indicate that this point was disqualified by
                     //   clustering analysis, despite no obvious problems in
                     //   its source files
  double scftol;     // The self-consistent field energy convergence criterion 
  double spE;        // The final energy (after post-HF methods, if applicable)
  double nucrepE;    // The nuclear repulsion energy
  double mp2corrE;   // The MP2 correlation energy
  double* crd;       // The coordinates found in the quantum output (converted
                     //   to units of A)
  char prog[16];     // The program used to do the calculation (this is not
                     //   filled by the user directly, so it can be short)
  char method[32];   // The quantum method used (this is also assigned by mdgx,
                     //   not the user directly)
  char* source;      // Source information for this data point (needed because
                     //   the data point could be found as part of a regular
                     //   expression search string, or part of a trajectory
                     //   among many other data points)
  cmat errors;       // Matrix to hold messages on any errors encountered  
};
typedef struct QuantumResults qmresult;

struct EvaluationItem {
  int takeqmcrd;     // Flag to have coordinates taken from the QM files
  int npts;          // The number of data points in this set.  This number is
                     //   determined after having read and processed all of the
                     //   items, and stands for the number of valid
                     //   conformations of the molecule for which single-point
                     //   energies are known.
  int npass;         // The number of points that passed sanity checks 
  char* crdsrc;      // Source file for coordinates (this must be a regular
                     //   file, though it can be a trajectory if and only if
                     //   the energy source file is a list of energies
                     //   containing the proper number of real values.
  char* esrc;        // The energy source file (could be a QM output file,
                     //   could be an mdgx-readable list of energies)
  qmresult* points;  // The data points collected from all the various items
};
typedef struct EvaluationItem sptask;

struct SinglePointData {
  int nitems;        // The number of items to read (this is not necessarily
                     //   the number of data points
  int enfCluster;    // Flag to have energy clustering enforced (this will cull
                     //   data points which cannot fit within a cluster of a
                     //   specified size) 
  int maxclst;       // Size of the largest cluster found within the data 
  int clstCull;      // The number of data points removed by clustering (these
                     //   points could not be included in a cluster of the
                     //   specified width)
  int verbose;       // Flag to control the level of output
  int PrintItems;    // Flag to have all items in the input file reprinted to
                     //   the output, regardless of input file length
  double clstWidth;  // The width of the cluster to enforce
  sptask* items;     // The single point data to hunt down (each item could be
                     //   a number of things--a coordinates file containing a
                     //   single conformation of the molecule paired with a QM
                     //   output file, a lone QM output file containing
                     //   coordinates, or a regular expression or directory
                     //   full of QM output files (containing coordinates) that
                     //   can be related to the topology on hand.
};
typedef struct SinglePointData spdata;

#endif
