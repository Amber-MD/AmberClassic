#ifndef ParamFitDataStructures
#define ParamFitDataStructures

#ifndef PREP_API
#include "Constants.h"

#include "ChargeFitDS.h"
#include "CrdManipDS.h"
#endif

struct EquivAtomGroup {
  int natom;
  char* types;
};
typedef struct EquivAtomGroup eagrp;

struct ExtendedAtomDef {
  int inreport;   // Flag to indicate that this atom type is needed in the
                  //   final report
  int dup;        // Flag to indicate that this atom type is a duplicate,
                  //   branched and renamed from one of the original types
  double mass;    // The atom mass
  double apol;    // The atomic polarizability
  double ljsig;   // The Lennard-Jones Sigma parameter
  double ljeps;   // The Lennard-Jones Sigma parameter
  char atype[8];  // The atom type
  char* comment;  // Comment from the parameter file
};
typedef struct ExtendedAtomDef xatomdef;

struct ExtendedBondDef {
  int fitcolX;    // Columns of the fitting matrix to which angles of this sort
  int fitcolY;    //   are mapped: when fitting only the stiffness constant,
  int fitcolZ;    //   the X column is used exclusively.  When fitting the
  int fitcolW;    //   stiffness and equilibrium constants, the Y column is
                  //   also used.  The Z ans W columns are used for fitting
                  //   pulling and compression augmentations, respectively.
  int dup;        // Flag to indicate that this bond type is a duplicate,
                  //   branched and renamed from one of the original types
  int samprow;    // Instances of the bond populate this row of the sampling
                  //   tables (coarse-grained and fine-grained)
  int ninst;      // The number of instances of this bond in the fit
  int inreport;   // Flag to indicate that this bond type is needed in the
                  //   final report
  int restKeq;    // Row of the restraint equation applied to the stiffness
                  //   in the main fitting matrix
  int restLeq;    // Row of the restraint equation applied to the equilibrium
                  //   length of the bond in the main fitting matrix
  int isAug;      // Flag to indicate that the bond has augmentations.  The
                  //   stiffness constants of these augmentations, but not
                  //   their equilibria, will be fitted if these terms are
                  //   present.  The input force field need not contain any
                  //   augmentation force constants if these parameters are to
                  //   be fitted: rather, the &param input can specify that
                  //   all bonds are to receive augmentation.
  double K;       // Spring constant in the original (input) force field
  double l0;      // Equilibrium length in the original (input) force field
  double Kpull;   // Stiffness and equilibrium length of the bond pulling
  double lpull0;  //   augmentation in the original (input) force field
  double Kpress;  // Stiffness and equilibrium length of the bond compression
  double lpress0; //    augmentation in the original (input) force field
  double Kfit;    // The fitted value of the stiffness constant
  double l0fit;   // The fitted value of the equilibrium length
  double lbasisX; // Equilibrium length for the X basis function
  double lbasisY; // Equilibrium length for the Y basis function
  double rstwK;   // Weight of the restraint equation to apply for the
                  //   stiffness constant
  double rstwl0;  // Weight of the restraint equation to apply for the
                  //   equilibrium constant
  double targK;   // Target of the stiffness restraint equation
  double targl0;  // Target of the equilibrium restraint equation
  char atype[8];  // The type of atom A
  char btype[8];  // The type of atom B
  char* comment;  // Comment from the parameter file
};
typedef struct ExtendedBondDef xbonddef;

struct ExtendedHBondDef {
  int fitcol;     // Column of the fitting matrix to which this maps
  int dup;        // Flag to indicate that this H-bond type is a duplicate,
                  //   branched and renamed from one of the original types
  int ninst;      // Instances of this H-bond
  int inreport;   // Flag to indicate that this is needed in the final report
  double Aterm;   // r^12, r^10 related constants 
  double Bterm;   //
  char atype[8];  // The type of atom A
  char btype[8];  // The type of atom B
  char* comment;  // Comment from the parameter file
};
typedef struct ExtendedHBondDef xhb1012def;

struct ExtendedAngleDef {
  int fitcolX;    // Columns of the fitting matrix to which angles of this sort
  int fitcolY;    //   are mapped: when fitting only the stiffness constant,
                  //   the X column is used exclusively.  When fitting the
                  //   stiffness and equilibrium constants, the Y column is
                  //   also used.
  int dup;        // Flag to indicate that this angle type is a duplicate,
                  //   branched and renamed from one of the original types
  int samprow;    // Instances of the angle populate this row of the sampling
                  //   tables (coarse-grained and fine-grained)
  int ninst;      // The number of instances of this angle in the fit
  int inreport;   // Flag to indicate that this bond type is needed in the
                  //   final report
  int restKeq;    // Row of the restraint equation applied to the stiffness
                  //   in the main fitting matrix
  int restTh0;    // Row of the restraint equation applied to the equilibrium
                  //   value of the angle in the main fitting matrix
  double K;       // Spring constant in the original (input) force field
  double th0;     // Equilibrium length in the original (input) force field
  double Kfit;    // The fitted value of the stiffness constant
  double th0fit;  // The fitted value of the equilibrium angle
  double tbasisX; // Equilibrium angle for the X basis function
  double tbasisY; // Equilibrium angle for the Y basis function
  double rstwK;   // Weight of the restraint equation to apply for the
                  //   stiffness constant
  double rstwTh0; // Weight of the restraint equation to apply for the
                  //   equilibrium constant
  double targK;   // Target of the stiffness restraint equation
  double targTh0; // Target of the equilibrium restraint equation
  char atype[8];  // The type of atom A
  char btype[8];  // The type of atom B
  char ctype[8];  // The type of atom C
  char* comment;  // Comment from the parameter file
};
typedef struct ExtendedAngleDef xangldef;

struct TorsionTerm {
  int fitcol;     // Column of the fitting matrix to which torsions of this
                  //   sort are mapped
  int dup;        // Flag to indicate that this torsion type is a duplicate,
                  //   branched and renamed from one of the original types
  int samprow;    // Instances of the torsion populate this row of the sampling
                  //   tables (coarse-grained and fine-grained)
  int ninst;      // The number of instances of this dihedral in the fit
  int inreport;   // Flag to indicate that this torsion type is needed in the
                  //   final report
  int impr;       // Flag to indicate that this dihedral is an improper
  int restKval;   // Row of the restraint equation applied to this torsion in
                  //   the main fitting matrix
  double K;       // The basic torsional potential for an angle phi is defined:
  double phase;   //   U = (pk/idivf) * (1 + cos(pn*phi - phase)), K = pk/idivf
  double pn;      //   All of these values apply to the input force field.
  double Kfit;    // The fitted amplitude of this term
  double singlet; // Indicates whether the dihedral is singlet or not; +1.0 if
                  //   there is only one term, -1.0 if there are more terms 
  double rstw;    // Weight of the restraint equation to apply
  double target;  // Target value of the restraint to apply
  char atype[8];  // 
  char btype[8];  // The A, B, C, and D atom types (i.e. CT, HC, N3) in a
  char ctype[8];  //   torsion arrangement  A--B--C--D or improper A(D)--B--C
  char dtype[8];  // 
  char* comment;  // Comment from the parameter file
};
typedef struct TorsionTerm torterm;

struct ExtendedCmapDef {
  int nres;        // The number of unique central residues to which this CMAP
                   //   applies (any residues labeled 0 in the relResID list
                   //   below is the "central" residue)
  int fitcoli;     // Initial fitting column for the CMAP in the stiffness
                   //   matrix
  int fitcolf;     // Final fitting column, the CMAP occupies columns
                   //   [ fitcoli, fitcolf ] inclusive
  int ninst;       // Number of instances of this CMAP in all structures from
                   //   the fitting data
  int inreport;    // Flag to indicate that this CMAP type is to be included
                   //   in the final report
  double rstw;     // Restraint holding points on the CMAP towards the target
  dmat target;     // Restraint target (this will probably be a mtrix of all
                   //   zeros or the original force field's CMAP)
  dmat stencil;    // Stencil defining the relationship between derivatives and
                   //   energy values in the map (this is a bridge to a linear
                   //   least squares fitting)
  dmat esrf;       // The energy surface of the input force field
  dmat dphi;       // The partial derivative with respect to phi (input FF)
  dmat dpsi;       // The partial derivative with respect to phi (input FF)
  dmat d2pp;       // The mixed partial derivative (input FF)
  dmat esrfFit;    // The fitted potential surface (this comprises all of the
                   //   fitted parameters)
  dmat dphiFit;    // Phi partial derivatives of the fitted potential surface
  dmat dpsiFit;    // Psi partial derivatives of the fitted potential surface
  dmat d2ppFit;    // Psi partial derivatives of the fitted potential surface
  cmat resnames;   // Names of the central residue (residue 0) to which this
                   //   CMAP applies.  This does NOT cover names of leading or
                   //   trailing residues.
  int relResID[5]; // Integers denoting whether atoms in the CMAP are to be
                   //   located on the central residue (0) or some leading or
                   //   trailin gresidue (i.e. leading -1, trailing +1)
  char atoms[24];  // Names of the atoms in this CMAP (four characters each),
                   //   defaults to C, N, CA, C, N
  char* title;     // Title of the CMAP
};
typedef struct ExtendedCmapDef xcmapdef;

struct NMROptimizationTerm {
  int fitcol;     // Column of the fitting matrix to which torsions of this
                  //   sort are mapped
  int samprow;    // Instances of the torsion populate this row of the sampling
                  //   tables (coarse-grained and fine-grained)
  int ninst;      // The number of instances of this NMR restraint in the fit
  int inreport;   // Flag to indicate that this torsion type is needed in the
                  //   final report
  int restKval;   // Row of the amplitude restraint equation applied to this 
                  //   NMR restraint the main fitting matrix
  int order;      // The order of this NMR restraint (or, the number of atoms
                  //   it applies to, as identified by atom names)
  int stylecode;  // The integer code for the style of this restraint (values
                  //   include 0 (halfcup, default) or 1 (sfunc))
  int nlayer;     // The number of separate NMR restraints at work in this one
                  //   operation
  int nsibling;   // The number of copies of this restraint that apply to other
                  //   contexts, likely equivalent atoms in the same molecule,
                  //   but could also be different combinations of atom types
  int usetypes;   // Flag to have mdgx seek out atom TYPES rather than atom
                  //   names for application of a restraint.  In either case
                  //   the atoms in question must be connected by a bond, bond
                  //   angle, or dihedral term in order for the NMR operation
                  //   to kick in.  Default 0 (use names), set to 1 to make use
                  //   of types instead.
  double rstw;    // The weight of the restraining potential used in fitting--
                  //   this is an NMR RESTRAINT, but there's also a restraining
                  //   potential applied to the amplitude during fitting to
                  //   keep things from going crazy.  The rstw attribute is
                  //   about the fitting restaint.
  double target;  // Desired stiffness constant for this restraint, analogous
                  //   to the eponymous attribute in structs such as
                  //   TorsionTerm and ExtendedAngleDef.
  double r1[6];   // Critical points along the FIRST NMR restraint potential,
  double r2[6];   //   matching r1, r2, r3, and r4 in the standard NMR
  double r3[6];   //   restraint namelist jargon
  double r4[6];   //
  double rk2[6];  // NMR restraint potential stiffness between r1 and r2, or
  double rk3[6];  //   between r3 and r4, as given in the INPUT force field
  cmat amask;     // The A, B, C, and D atom NAMEs (i.e. HA, CA, N) and their
  cmat bmask;     //   respective residue names in the NMR restraint A -> B ->
  cmat cmask;     //   C -> D.  The atoms must be either on the same residue
  cmat dmask;     //   or on consecutive residues in the order those residues
                  //   appear in the topology file.
  char* style;    // The style of restraint to use (acceptable values include
                  //   "halfcup" and "sfunc")
  char* label;    // The label for this NMR optimization term (essential if
                  //   the user wants to specify particular operations to
                  //   optimize)
};
typedef struct NMROptimizationTerm nmroper;

struct BondRestraint {
  int has_wK;     // Boolean integer (0 = False, 1 = True) indicating whether
                  //   there is information on the stiffness restraint weight
  int has_wl0;    // Boolean integer on the equilibrium restraint weight
  int has_tK;     // Boolean integer on the stiffness target
  int has_tl0;    // Boolean integer on the equilibrium target
  char atype[8];  // The A and B atom types (i.e. CT, HC, N3) in the bond
  char btype[8];  //
  double rstwK;   // Weight that the user wishes to assign to the restraint of
                  //   the bond stiffness for these atom types
  double rstwl0;  // Weight that the user wishes to assign to the restraint of
                  //   the bond equilibrium for these atom types
  double targK;   // Target value for the bond stiffness
  double targl0;  // Target value for the bond equilibrium
};
typedef struct BondRestraint bondrst;

struct AngleRestraint {
  int has_wK;     //
  int has_wTh0;   // Boolean integers similar to those in the BondRestraint 
  int has_tK;     //   struct indicating what restraint information is present
  int has_tTh0;   //
  char atype[8];  //
  char btype[8];  // The A, B, and C atom types (i.e. CT, HC, N3) in the angle
  char ctype[8];  //
  double rstwK;   // Weight that the user wishes to assign to the restraint of
                  //   the angle stiffness for these atom types
  double rstwTh0; // Weight that the user wishes to assign to the restraint of
                  //   the angle equilibrium for these atom types
  double targK;   // Target value for the angle stiffness
  double targTh0; // Target value for the angle equilibrium
};
typedef struct AngleRestraint anglrst;

struct TorsionRestraint {
  char atype[8];  //
  char btype[8];  // The A, B, C, and D atom types (i.e. CT, HC, N3) in the
  char ctype[8];  //   torsion
  char dtype[8];  //
  double rstw;    // Weight that the user wishes to assign to the restraint of
                  //   the torsion for these atom types
  double pn;      // Periodicity to which this restraint applies (values < 0
                  //   imply that it pertains to all frequencies and any target
                  //   value should be ignored)
  double target;  // Target value for the restraint
};
typedef struct TorsionRestraint torrst;

struct GeometryRestraint {
  int nvar;       // Number of angles being collectively restrained
  double target;  // The target sum for the angle equilibria in this restraint
  double srange;  // The search range for the constrained optimization
  int* fitcol;    // Fitting columns for each angle variable
  int* rstrow;    // Restraint rows for each angle variable
  double* basis;  // Basis function minima for each angle variable
  cmat atmtype;   // The atom types of angles in this restraint, listed in the
                  //   order (typeA, typeB, typeC), (typeA, typeB, typeC), ...
};
typedef struct GeometryRestraint geomrst;

struct BondIndex {
  int a;          // The A and B atoms of the bond, indexed according to the
  int b;          //   system's own topology
  int key;        // The index into the master list of bonds, stored in the
                  //   parameter set and spanning all systems
};
typedef struct BondIndex bidx;

struct AnglIndex {
  int a;          // The A, B, and C atoms of the angle, indexed according to
  int b;          //   the system's own topology
  int c;          //
  int key;        // The index into the master list of angles, stored in the
                  //   parameter set and spanning all systems
};
typedef struct AnglIndex aidx;

struct TorsionIndex {
  int a;          //
  int b;          // The A, B, C, and D atoms of the torsion term, indexed
  int c;          //   according to the system's own topology
  int d;          //
  int key;        // The index into the master list of torsions, stored in the
                  //   parameter set and spanning all systems
};
typedef struct TorsionIndex hidx;

struct CmapIndex {
  int a;          //
  int b;          // The A, B, C, D, and E atoms of the CMAP term, indexed
  int c;          //   according to the system's own topology
  int d;          //
  int e;          //
  int key;        // The index into the master list of cmaps, stored in the
                  //   parameter set and spanning all systems
};
typedef struct CmapIndex midx;

struct NMROperationIndex {
  int order;      // The number of atoms in the NMR restraint
  int a;          //
  int b;          // The A, B, C, and D atoms of the torsion term, indexed
  int c;          //   according to the system's own topology
  int d;          //
  int key;        // The index into the master list of NMR operations, stored
                  //   in the parameter set and spanning all systems
};
typedef struct NMROperationIndex ridx;

struct BondMap {
  int nbond;        // The number of bonds in this system
  bidx* id;         // Numbers of atoms participating in each bond, and indices
                    //   into the master list of bonds
  double* val;      // Values of the underlying bond length in the fitting set
  double* UkernelX; // The contributions to the fitting matrix; X and Y columns
  double* UkernelY; //   may be invoked depending on the type of bond fitting,
  double* UkernelZ; //   and Z and W columns describe bond augmentations which
  double* UkernelW; //   may or may not be up for fitting
  double* Ucontrib; // The kernel x stiffness constant = contribution to energy
};
typedef struct BondMap bondmap;

struct AngleMap {
  int nangl;        // The number of angles in this system
  aidx* id;         // Indices into the master list of angles
  double* val;      // Values of the underlying angle in the fitting set
  double* UkernelX; // The contributions to the fitting matrix; X and Y columns
  double* UkernelY; //   may be invoked depending on the type of angle fitting
  double* Ucontrib; // The kernel x stiffness constant = contribution to energy
};
typedef struct AngleMap anglmap;

struct TorsionMap {
  int ntterm;       // The number of torsional terms in this system
  hidx* id;         // Indices into the master list of torsion terms
  double* val;      // Values of the underlying torsion angle in the
                    //   fitting set
  double* Ukernel;  // The contributions to the fitting matrix
  double* Ucontrib; // The kernel x stiffness constant = contribution to energy
};
typedef struct TorsionMap tormap;

struct SplineMap {
  int ncmap;        // The number of CMAP terms in this system (conformation)
  midx* id;         // Indices into the master list of CMAP surfaces
  double* val;      // Values of the underlying pairs of dihedrals.  This array
                    //   has 2 * ncmap elements, and coordinates in the CMAP
                    //   plane are listed in (phi, psi) pairs.
  double* Ucontrib; // The kernel times contribution stencil = one scalar for
                    //   the contribution to the energy
};
typedef struct SplineMap splnmap;

struct NMROperationMap {
  int nops;         // The number of NMR operations in this system
  ridx* id;         // Indices into the master list of operations
  double* val;      // Values of the underlying coordinates in the fitting set
                    //   (the nature of these quantities depends on the number
                    //   of atoms in the NMR restraint)
  double* Ukernel;  // The contributions to the fitting matrix
  double* Ucontrib; // The contributions to the molecular mechanics energy
};
typedef struct NMROperationMap nmrmap;

struct EnergyContributor {
  int fitcol;       // The fitting column that this packet pertains to
  double eave;      // The average energy
  double estd;      // Standard deviation in the energy
};
typedef struct EnergyContributor epacket;

struct AtomTypeSwitch {
  char orig[8];     // The original atom type name
  char pnew[8];     // The new atom type name
};
typedef struct AtomTypeSwitch typeswitch;

struct AtomTypeBranch {
  char instances[MAXLINE];  // The instances in which the atom type of name
                            //   orig is to be recast to new
  char orig[8];             // The original atom type name
  char pnew[8];             // The new atom type name
};
typedef struct AtomTypeBranch typebranch;

struct MMSystem {
  int GroupNum;    // Number of the topology group, computed in GroupSystems 
                   //   to correlate systems with similar topologies
  int PassedEtol;  // Flag to indicate whether this conformation passed the
                   //   check on conformational energy according to the
                   //   original MM force field (0 if not, 1 if the
                   //   conformation passed without problems, 2 if the
                   //   conformation passed after rearrangement)
  prmtop *tp;      // System topology pointer (points to tpencyc)
  prmtop *vactp;   // Pointer to a secondary topology containing vacuum phase
                   //   charge parameters.  This topology is obtained from a
                   //   file name in the esrc field, if the file in the crdsrc
                   //   field points to an IPolQ energy + coordinates bundle
                   //   or a regular expression / directory of such things.
  coord crd;       // System coordinates
  bondmap bmap;    // Bond terms mapping to a list spanning the entire fit
  anglmap amap;    // Angle terms mapping to a list spanning the entire fit
  tormap hmap;     // Torsion terms mapping to a list spanning the entire fit
  splnmap mmap;    // CMAP terms mapping to a list spanning the entire fit
  nmrmap rmap;     // NMR operations mapping to a list spannign the entire fit
  dmat excl;       // Matrix of exclusions; 0.0 for total exclusion (1:1
                   //   virtual site anchoring, 1:2 bonded, and 1:3 angle
                   //   interactions fall under this category), 1.0 for no
                   //   exclusion (nonbonded interactions 1:5 and more distal
                   //   have no exclusions), and any other value for a partial
                   //   exclusion (1:4 interactions)
  dmat nbnrg;      // The nonbonded energy matrix, electrostatics above the
                   //   diagonal and van-der Waals interactions below it
  double EEkernel; // The kernel for 1-4 adjustable electrostatic interactions 
  double LJkernel; // The kernel for 1-4 adjustable Lennard-Jones interactions 
  double EEnonfit; // The unfitted, unscaled non-bonded electrostatic energy
  double LJnonfit; // The unfitted, unscaled non-bonded Lennard-Jones energy
  double etrg;     // The target energy of the conformation, derived from
                   //   quantum calculations most likely
  double enorm;    // The "normalized" energy of this conformation--obtained by
                   //   subtracting from this system's etrg the average of etrg
                   //   from all conformations sharing the same topology 
  double eorig;    // The final energy of this conformation according to the
                   //   input Hamiltonian (topology file)
  double efin;     // The final energy of this conformation according to the
                   //   fitted parameters
  double nonfitmm; // The energy of non-adjustable molecular mechanics terms
  double wt;       // The weight that this conformation will get in the fit
  char* crdsrc;    // Source file for the coordinates
  char* tpsrc;     // The source file for the topology (prmtop format)
  char* esrc;      // The source file or number-as-string for the energy
};
typedef struct MMSystem mmsys;

struct InstanceTracker {
  int sysid;       // The system ID number
  int sysno;       // The system number in the master list of conformations
  int order;       // The order of this instance (2, 3, or 4 for bonds, angles,
                   //   and dihedrals)
  int tnum;        // The number of the term within system structs for orders
                   //   2, 3, or 4
  char res[32];    // The residue names of atoms involved in this instance
  char atom[32];   // The atom names of atoms involved in this instance
};
typedef struct InstanceTracker itrack;

struct SpectrumRequest {
  int level;       // The level at which parameters matching this request shall
                   //   be included in reoptimization: 1 for inclusion, 2 for
                   //   explicit spectral resampling.
  int order;       // The order of the spectrum request (2 = bonds, 3 = angles,
                   //   4 = dihedrals, and the default is 4)
  double minval;   // The minimum value of the spectrum to which parameters
                   //   shall be constrained
  double maxval;   // The maximum value of the spectrum to which parameters
                   //   shall be constrained
  double gspc;     // The grid spacing for the spectrum
  cmat atmtypes;   // The A, B, C, and D atom types (i.e. CT, HC, N3) in this 
                   //   request.  At least one, but up to four atom types may
                   //   be specified, and all types must be present in order
                   //   for a dihedral to be subject to spectral evaluation.
};
typedef struct SpectrumRequest specreq;

struct SpectrumControl {
  int level;       // The level at which to include   
  int order;       // The order of the parameter (2 = bonds, 3 = angles, 4 =
                   //   torsions, 0 = other)
  int colXmain;    // The (first) column of the fitting matrix in which this
                   //   parameter may be found in the full fitting matrix
  int colYmain;    // The (second) column of the fitting matrix in which this
                   //   parameter may be found in the full fitting matrix
  int colXred;     // Like colXmain but for the reduced matrix
  int colYred;     // Like colYmain but for the reduced matrix
  int rstXmain;    // Index of the restraint equation that tethers the X
                   //   column variable, or the sum of X and Y, in the main
                   //   (unreduced) matrix
  int rstYmain;    // Index of the restraint equation that tethers the ratio
                   //   of X and Y column variables in the main (unreduced)
                   //   matrix
  int rstXred;     // Like rstXmain but for the reduced matrix
  int rstYred;     // Like rstYmain but for the reduced matrix
  double llimX;    // Lower limit of sampling the value of X for this variable
  double hlimX;    // Upper limit of sampling the value of X for this variable
  double llimY;    // Lower limit of sampling for the ratio of X and Y
  double hlimY;    // Upper limit of sampling for the ratio of X and Y
  double gspcX;    // Sampling discretization for the value of X
  double gspcY;    // Sampling discretization for the ratio of X and Y
};
typedef struct SpectrumControl speccon;

struct AlternateSolution {
  int ndim;        // Number of dimensions in the solution (this information
                   //   will likely be available by simply referencing the
                   //   column count of the main fitting matrix, but is copied
                   //   for convenience)
  double* x;       // The alternative solution
  double* b;       // The result of the alternative solution when applied to
                   //   all conformations in the training set
  imat rstpos;     // Indices of restraint information, given in the order
                   //   (row, column1, column2, ..., columnN), all pertaining
                   //   to the main fitting matrix
  dmat rstval;     // Values of restraint equations, given in the order
                   //   (column1, column2, ..., columnN, target).  Each row of
                   //   rstval corresponds to a row of rstpos; the first index
                   //   of the row in rstpos provides the row of the fitting
                   //   matrix in which the restraint shall be applied.
};
typedef struct AlternateSolution altsol;

struct SpectrumResampler {
  int naltsoln;      // The number of alternative solutions found by spectral
                     //   resampling (NOT the number actually requested, which
                     //   is found in the ParameterFit struct below and
                     //   transferred to this data struct in the form of the
                     //   row or column count of ssd1/2).
  dmat *A;           // The original coefficient matrix in the linear least-
                     //   squares problem (LLSP).
  double* x;         // The original solution to the LLSP.  As in, Ax = b.
  double* b;         // The original target in the LLSP.
  speccon* specvar;  // Spectrum control data array (this will be saved from a
                     //   free so that it can be useful later in the output) 
  altsol* xalt;      // Array of alternate solutions to Ax = b, with critical
                     //   information on what to do to the matrix A to achieve
                     //   a similar result.  These solutions are not exactly
                     //   what will come out of Ax = b because they were made
                     //   with a reduced coefficient matrix, but they will be
                     //   close.  The first element of this list is the global
                     //   minimum.
  dmat ssd1;         // Matrix of RMSDs encoding the differences in predicted
                     //   energies, across the entire training set, among the
                     //   alternative solutions selected by mutual
                     //   dissimilarity.  These RMSDs are computed from
                     //   alternative solutions emerging from the reduced
                     //   matrix equation.
  dmat ssd2;         // Like ssd1 but for alternative solutions refined after
                     //   solving the full matrix equation.  Compare to ssd1 to
                     //   verify that the predictions did not change much and
                     //   that the matrix reduction was sound.
};
typedef struct SpectrumResampler sresamp;

struct ParameterWarning {
  int alarm;         // The level of alarm to raise
  int nminima;       // The number of minima detected
  double* minloc;    // The locations of the minima
  char* output;      // Output string to describe this situation
};
typedef struct ParameterWarning parmwarn;

struct ParameterFit {
  int nconf;          // The number of conformations in this fitting set
  int natom;          //
  int nbond;          // The number of unique atom, bond, hydrogen bond, angle,
  int nhb1012;        //   and torsion terms in the parameter file and (if 
  int nangl;          //   supplied) in the frcmod file
  int ntor;           // 
  int ndihe;          // The number of proper dihedrals
  int nimpr;          // The number of impropers (ndihe + nimpr = ntor)
  int ncmaps;         // The number of CMAP surfaces
  int cmapFitResln;   // The resolution for fitted CMAPs (may differ from their
                      //   resolution in the input force field)
  int nops;           // The number of NMR operations (these will be entered in
                      //   a file separate from the standard parm and frcmod
                      //   parameter files, but in spirit it is like another
                      //   parameter file and the number of NMR operations is
                      //   analogous to the number of bonds, angles, etc.)
  int nrestraint;     // The total number of restraints
  int nparm;          // The number of adjustable parameters, which must be no
                      //   greater than nbond + nangl + ntor + nscee + nscnb
  int nunisys;        // The number of unique systems (each topology file with
                      //   a unique name identifies a unique system)
  int nchng;          // The number of changes that have been made to atom
                      //   types in specific cases across all topologies
  int ljbuck;         // This parameter serves as a placeholder for the prmtop
                      //   struct attribute of the same name.  The value found
                      //   in this field will be passed down into every other
                      //   topology file in the conf array
  int ipolq;          // Flag to indicate that the fitting data contains
                      //   IPolQ charges, and in fact pairs of topologies with
                      //   hyperpolarized and vacuum phase charges.
  int FitAllBonds;    //
  int FitAllAngles;   // Flags to activate fitting for all identified bonds,
  int FitAllTorsions; //   bond angles, torsions, CMAPs, and NMR operations
  int FitAllCMAPs;    // 
  int FitAllNMROps;   //
  int FitBondEq;      // Flags to activate bond and angle equilibrium fitting
  int FitAnglEq;      //
  int FitAllAugs;     // Flag to fit augmentations on all bonds involved in the
                      //   fitting
  int FitCurrentAugs; // Flag to fit augmentations on bonds with augmentations
                      //   in the input parameter set
  int reportall;      // Flag to indicate that all parameters encountered
                      //   should be reported
  int zeroNonfit;     // Flag to indicate that only fitted energy terms should
                      //   contribute to the molecular mechanics energy
  int nbadj;          //
  int naadj;          // The numbers of adjustable bonds, angles, torsions,
  int nhadj;          //   CMAPs, and NMR operations specified by the user
  int nmadj;          //
  int nradj;          //
  int nbvar;          //
  int navar;          // The number of adjustable bonds, angles, torsions,
  int nhvar;          //   and NMR operations actually found
  int nmvar;          //
  int nrvar;          //
  int nuserrstB;      // Number of user-specified bond restraint stiffnesses
  int nuserrstA;      // Number of user-specified angle restraint stiffnesses
  int nuserrstH;      // Number of user-specified torsion restraint stiffnesses
  int nuserrstR;      // Number of user-specified restraints for fitting NMR
                      //   operation scaling factors
  int nanglsum;       // Number of geometry restraints
  int nzerocol;       // The number of zero-values columns in the fitting
                      //   matrix; this is generally not a problem as all
                      //   fitted parameters are restrained
  int ncorrcol;       // The number of highly correlated columns in the fitting
                      //   matrix
  int nrecast;        // The number of atom types to recast
  int ncleave;        // The number of atom types to branch (cleave)
  int fitscnb;        // Flag to activate fitting of the scnb (Lennard-Jones
                      //   1:4 scaling) term
  int fitscee;        // Flag to activate fitting of the scee (electrostatic
                      //   1:4 scaling) term
  int neqgroups;      // Number of equivalent atom groups
  int verbose;        // Display progress for the user (default 1, yes)
  int RemoveOutliers; // Flag to activate removal of conformations whose
                      //   energies are outliers (default 0, no)
  int spectrum;       // Solve the matrix equation multiple times given a
                      //   spectrum of values for particular variables
  int nspectralH;     // The number of dihedrals for which spectrum results
                      //   are desired
  int nspectralx;     // The number of alternative solutions that spectrum
                      //   resampling is to produce
  int PrintFitPoints; // Directive to have all of the names of the original
                      //   fitting data reprinted to the output file
  int RestScaleMeth;  // Method for determining the restraint scaling.  Default
                      //   (0) behavior is to scale by the summing the number
                      //   of instances of each optimizable parameter in any
                      //   one system times the number of conformations of that
                      //   system, for all systems in the fitting problem.
                      //   Legacy behavior, which looks at the fitting matrix
                      //   itself and counts the number of non-zero entries in
                      //   that parameter's column, can be enabled by setting
                      //   this to 1.
  int* zerocol;       // List of zero-valued columns
  int* corrcol;       // List of pairs of highly correlated columns
  int* GroupCount;    // The number of each system, as differentiated by 
                      //   topologies
  int* FirstConf;     // The first conformations of each topology / system
  int* LastConf;      // The last conformations of each topology / system
  double lj14fac;     // The Lennard-Jones 1:4 scaling default value
  double elec14fac;   // The electrostatic 1:4 scaling default value
  double grstB;       // The general constraints on bond, angle, torsion, CMAP,
  double grstA;       //   and NMR operation restraint stiffnesses, to keep
  double grstH;       //   these values tethered to zero and thus small in the
  double grstM;       //   result
  double grstR;       //
  double grstBcpl;    // The coupling factor between restraints on the bond
                      //   spring constants and equilibria (default penalizes
                      //   a 0.01A change in the equilibrium by the same amount
                      //   as a 50 kcal/mol change in spring constant)
  double grstAcpl;    // The coupling factor between restraints on the angle
                      //   spring constants and equilibria (default penalizes
                      //   a 1-degree change in equilibrium by the same amount
                      //   as a 2 kcal/mol change in spring constant)
  double grst14;      // The 1:4 scaling term constraint factor (if 1:4 scaling
                      //   terms are being fitted)
  double mmtol;       // The molecular mechanics (bonded terms) energy
                      //   tolerance; configurations with energies higher than
                      //   this will be flagged for re-arrangement and reported
                      //   if the energy cannot be reduced
  double esigtol;     // Tolerance for the deviation of the total energy of any
                      //   particular conformation from the mean for all
                      //   conformations of that system
  double fsigtol;     // Tolerance for the deviation of the fitted energy of
                      //   any particular conformation from the target value
  double fdevfloor;   // Minimum value at which a fitted energy can be flagged
                      //   as an outlier
  double wtfloor;     // The minimum weight that any conformation may have;
                      //   default 0.5.  Set to smaller values to increase the
                      //   slant of the training set towards data points with
                      //   favorable energies.  Larger values make the training
                      //   set treat all conformations with equanimity.
  double lpost;       // The distance from the original equilibrium value, plus
                      //   and minus, to separate the equilibrium bond lengths
                      //   in each of two basis functions for fitting a new
                      //   bond.  All basis functions for every angle are
                      //   spaced by the same value.  
  double thpost;      // The distance from the original equilibrium value, plus
                      //   and minus, to separate the equilibrium angles in
                      //   each of two basis functions for fitting a new angle.
                      //   All basis functions for every angle are spaced by
                      //   the same value.  
  double  spvtol;     // The tolerance for new solutions to the parameter set
                      //   when conducting spectrum resampling, a relative 
                      //   factor based on the global optimum RMSE for any
                      //   given system
  double* corrval;    // Correlations among highly aligned columns
  mmsys* conf;        // The system conformations, each with its own 
                      //   coordinates and topology
  xbonddef* badj;     // The array of adjustable bonds
  xangldef* aadj;     // The array of adjustable angles
  torterm* hadj;      // The array of adjustable torsions
  xcmapdef* madj;     // The array of adjustable CMAP surfaces
  nmroper* radj;      // Adjustable NMR restraints
                      //   (the arrays of adjustable bonds, bond angles,
                      //    torsions, and NMR operations will store the final
                      //    fitted parameters at the end of the fit)
  xatomdef* atoms;    // Atom definitions (masses and type names are included)
  xbonddef* bonds;    // Bond definitions (stiffnesses and lengths can be
                      //   fitted parameters in a linear least-squares fit)
  xhb1012def* hb1012; // Hydrogen bondind 10-12 parameters (currently only kept
                      //   for the purposes of printing them out in the final
                      //   parameter file)
  xangldef* angls;    // Angle definitions (stiffnesses and equilibria can be
                      //   fitted by linear least-squares optimization)
  torterm* torsions;  // Torsion terms (amplitudes are fitted by linear
                      //   least-squares optimization)
  xcmapdef* cmaps;    // CMAP terms (derivatives at any one grid point are
                      //   derived from energies at other grid points)
  nmroper* nmrops;    // NMR restraints (amplitudes are fitted by linear least-
                      //   least-squares optimization)
  bondrst* userrstB;  // An array of user-specified restraint constants for
                      //   fitting particular bonds
  anglrst* userrstA;  // An array of user-specified restraint constants for
                      //   fitting particular angles
  torrst* userrstH;   // An array of user-specified restraint constants for
                      //   fitting particular torsions
  nmroper* userrstR;  // An array of user-specified restraints for fitting
                      //   the scaling constants for NMR operations
  geomrst* anglsum;   // Angle sum restraints (geomrst structs), used to keep
                      //   the sum of multiple angles constrained to particular
                      //   values (i.e. 360 degrees)
  specreq* spectralH; // An array of definitions for dihedrals that are to be
                      //   subjected to spectrum restrained optimization.
  dmat LJAcoef;       // Lennard-Jones A and B coefficient matrices
  dmat LJBcoef;       //
  cmat Hydrophilics;  // Hydrophilic atom types
  cmat ChangeLog;     // Log of changes to atom types
  eagrp* eqgroups;    // Equivalent atom groups
  prmtop* tpencyc;    // The encyclopedia of topologies in this fitting run
  prmtop* vactpencyc; // Encyclopedia of topologies with vacuum phase charges,
                      //   ordered the same as tpencyc
  typebranch* cleave; // Array of types that are to be branched
  typeswitch* recast; // Array of types whose names are to be recast
  char* ititl;        // Title for output parameter file
  char* icomm;        // Comment for fitted parameters
  char* NMROpsFile;   // NMR operations file
  char* NMRParmFile;  // NMR parameter file (follows the same format as the ops
                      //   file, but is written after fitting and may have more
                      //   detail in terms of the stiffness coefficients)
  char NrgUnits[64];  // Units of the energies for conformations supplied in
                      //   the input file
  char WaterName[8];  // The name of water molecules (for implementing SETTLE),
                      //   serves as a placeholder for the prmtop attribute of
                      //   the same name
  char ep[MAXNAME];   // The name of the EP rules file, a placeholder for the
                      //   prmtop eprulesource attribute; only one file may be
                      //   specified, as it will be passed down to all
                      //   topologies in the conf array
  char sr[MAXNAME];   // The series report file; if specified every torsion
                      //   series solved in the fit will be printed
  char ao[MAXNAME];   // The accuracy output file; if specified a MatLab
                      //   script will be printed to display the results,
                      //   system by system

  // Repositories of all the systems' parameter and term contributions.  These
  // matrices are allocated for the sake of keeping memory as contiguous as
  // possible without wasting space.
  bidx** BondMapIDBucket;       // Bond atom and term type indices
  aidx** AnglMapIDBucket;       // Angle atom and term type indices
  hidx** DiheMapIDBucket;       // Dihedral atom and term type indices
  ridx** NOpsMapIDBucket;       // NMR Pperation atom and term type indices
  midx** SplnMapIDBucket;       // CMAP atom and surface indices
  dmat* BondMapValBucket;       //
  dmat* BondMapUxBucket;        //
  dmat* BondMapUyBucket;        // Bond data
  dmat* BondMapUzBucket;        //
  dmat* BondMapUwBucket;        //
  dmat* BondMapCntrbBucket;     //
  dmat* AnglMapValBucket;       //
  dmat* AnglMapUxBucket;        // Angle data
  dmat* AnglMapUyBucket;        //
  dmat* AnglMapCntrbBucket;     //
  dmat* DiheMapValBucket;       //
  dmat* DiheMapUBucket;         // Dihedral data
  dmat* DiheMapCntrbBucket;     //
  dmat* NOpsMapValBucket;       //
  dmat* NOpsMapKernelBucket;    // NMR Operations data
  dmat* NOpsMapCntrbBucket;     //
  dmat* SplnMapValBucket;       // CMAP data
  dmat* SplnMapCntrbBucket;     //
};
typedef struct ParameterFit prmset;

#endif
