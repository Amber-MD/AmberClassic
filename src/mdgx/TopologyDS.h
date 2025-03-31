#ifndef TopologyStructs
#define TopologyStructs

#ifndef PREP_API
#include "Constants.h"
#include "MatrixDS.h"
#endif

//-----------------------------------------------------------------------------
// BNDANGDIHE: a group of three integers, bonds, angles, and dihedrals. 
//-----------------------------------------------------------------------------
struct bndangdihe {
  int nbond;
  int nangl;
  int ndihe;
};
typedef struct bndangdihe bah;

//-----------------------------------------------------------------------------
// BONDIDX: a bond index structure matching the AMBER topology file format.
//-----------------------------------------------------------------------------
struct bondidx {
  int a;          // The "A" atom in the A-B bond, as shown below:
                  //
                  //         A--B
                  //
                  // The "A" atom serves as the "anchor" when determining
                  // which cell or group of cells is charged with computing
                  // the bond force and energy
  int b;          // The "B" atom
  int idx;        // An index into the bond parameter arrays
};
typedef struct bondidx bond;

//-----------------------------------------------------------------------------
// BONDCOMMAND: the bond command structure, containing information needed
//              during program execution for evaluating a particular bond.
//-----------------------------------------------------------------------------
struct bondcommand {
  int a;     // The "A" atom in the bond
  int b;     // The "B" atom in the bond
  int t;     // An index into the bond parameter arrays
  char H;    // Flag to tell whether a bond has hydrogen in it
};
typedef struct bondcommand bondcomm;

//-----------------------------------------------------------------------------
// ATOMBONDLIST: a list of bond commands for a particular atom.         
//-----------------------------------------------------------------------------
struct atombondlist {
  int nbond;     // The number of bonds this atom controls
  bondcomm* BC;  // The bond commands
};
typedef struct atombondlist bondlist;

//-----------------------------------------------------------------------------
// BONDDEF: a structure for defining a bond in terms that are more convenient
//          to reference during program execution.           
//-----------------------------------------------------------------------------
struct BondDef {
  double K;        // Spring constant
  double l0;       // Equilibrium length
  double Kpull;    // Spring constant and equilibrium length that kick in
  double lpull0;   //   when the bond is stretched a lot
  double Kpress;   // Spring constant and equilibrium length that kick in
  double lpress0;  //   when the bond is compressed a lot
};
typedef struct BondDef bonddef;

//-----------------------------------------------------------------------------
// ANGLIDX: an angle index structure matching the AMBER topology file format.
//-----------------------------------------------------------------------------
struct anglidx {
  int a;          // The "A" atom in the A-B-C angle, as shown below:
                  //
                  //              C
                  //             /
                  //         A--B
                  //
  int b;          // The "B" atom (this atom serves as the "anchor" when
                  // determining which cell or group of cells is charged with
                  // computing the angle force and energy)
  int c;          // The "C" atom
  int idx;        // An index into the angle parameter arrays
};
typedef struct anglidx angle;

//-----------------------------------------------------------------------------
// ANGLECOMMAND: the angle command structure, containing information needed
//               during program execution for evaluating a particular angle.
//-----------------------------------------------------------------------------
struct anglecommand {
  int a;     // The "A" atom in the angle
  int c;     // The "C" atom in the angle
  int t;     // An index into the angle parameter arrays
  int excl;  // Flag to active exclusion of the 1:3 interaction
};
typedef struct anglecommand anglcomm;

//-----------------------------------------------------------------------------
// ATOMANGLELIST: a list of bond commands for a particular atom.        
//-----------------------------------------------------------------------------
struct atomanglelist {
  int nangl;     // The number of angles this atom controls
  anglcomm* AC;  // The angle commands
};
typedef struct atomanglelist angllist;

//-----------------------------------------------------------------------------
// ANGLEDEF: a structure for defining an angle in terms that are more   
//           convenient to reference during program execution.          
//-----------------------------------------------------------------------------
struct AngleDef {
  double K;    // Spring constant
  double th0;  // Equilibrium angle
};
typedef struct AngleDef angldef;

//-----------------------------------------------------------------------------
// DIHEIDX: a dihedral index structure matching the AMBER topology file format.
//-----------------------------------------------------------------------------
struct diheidx {
  int a;          // The "A" atom in the A-B-C-D dihedral, as shown below:
                  //
                  //              D
                  //             /
                  //         B--C
                  //        /
                  //       A
                  //
  int b;          // The "B" atom (this atom serves as the "anchor," as
                  // discussed below)
  int c;          // The "C" atom
  int d;          // The "D" atom
  int idx;        // An index into the dihedral parameter arrays
};
typedef struct diheidx dihedral;

//-----------------------------------------------------------------------------
// DIHEDRALCOMMAND: this structure lists multiple Fourier terms that may be
//                  associated with a dihedral, and also whether it is an 
//                  improper.  This structure is what mdgx relies on when
//                  actually evaluating dihedrals.              
//-----------------------------------------------------------------------------
struct dihedralcommand {
  int a;          // The "A" atom (note that the b atom is implicit in the
                  // ordered list in which this dihedralcomm structure appears)
  int c;          // The "C" atom
  int d;          // The "D" atom
  int impr;       // Flag to tell us whether this is an improper
  int eval14;     // Flag to tell us whether to evaluate the 1:4 interactions 
  double scnb;    // The 1-4 van-der Waals scaling factor
  double scee;    // The 1-4 electrostatic scaling factor
  int nt;         // The number of Fourier terms
  int* t;         // A list of indices into the dihedral definitions of this
                  //   topology (see the diehdraldefs structure)
};
typedef struct dihedralcommand dihecomm;

//-----------------------------------------------------------------------------
// DIHEDRALLIST: a list of dihedrals associated with a particular atom, more
//               suitable for my dihedral calculation schemes.     
//-----------------------------------------------------------------------------
struct dihedrallist {
  int ndihe;          // The number of dihedral angles this atom controls
  dihecomm* HC;       // The list of dihedral commands
};
typedef struct dihedrallist dihelist;

//-----------------------------------------------------------------------------
// DIHEDRALDEF: a structure for defining a single Fourier series term in a
//              dihedral interaction.  This structure is used, rather than the
//              data read directly from the topology file, during program
//              execution.                               
//-----------------------------------------------------------------------------
struct DihedralDef {
  double K;     // The amplitude
  double N;     // The periodicity
  double Phi;   // The phase angle
};
typedef struct DihedralDef dihedef;

//-----------------------------------------------------------------------------
// CMAPIDX: a structure for storing a CMAP, taken directly from the topology.
//          This, like the dihedral indexing, will be converted into CMAP
//          commands.
//-----------------------------------------------------------------------------
struct cmapidx {
  int a;     //
  int b;     // The CMAP spans two dihedral interactions, A-B-C-D and B-C-D-E.
  int c;     //   These are the atom indices for those interactions within the
  int d;     //   topology.
  int e;     //
  int idx;   // Index into the stored array of CMAP surfaces
};
typedef struct cmapidx cmapterm;

//-----------------------------------------------------------------------------
// CMAPCOMMAND: a structure to store CMAP interactions and surface references
//-----------------------------------------------------------------------------
struct CMAPcommand {
  int a;          // The "A" atom of the CMAP, which comprises two coupled
                  //   dihedral interactions A-B-C-D and B-C-D-E.  The "C"
                  //   atom (not listed here) controls the CMAP and is known
                  //   by the same mechanism that the "B" atom controlling a
                  //   dihedral is known: mdgx loops over all atoms and stops
                  //   at each one to compute the bonds / angles / dihedrals /
                  //   CMAPs that the atom controls.
  int b;          //
  int d;          // The "B", "D", and "E" atoms of the CMAP
  int e;          //
  int idx;        // Index of the CMAP in the stored array of surfaces
  int eval14i;    // Flags to evaluate 1-4 interactions for the A-B-C-D or
  int eval14ii;   //   B-C-D-E dihedrals, respectively
  double scnbi;   // 
  double scnbii;  //
  double sceei;   //
  double sceeii;  //
};
typedef struct CMAPcommand cmapcomm;

//-----------------------------------------------------------------------------
// CORRECTIONMAPLIST: a list of CMAP commands, similar to the dihedral list
//                    and dihedral commands.
//-----------------------------------------------------------------------------
struct CorrectionMapList {
  int ncmap;      // The number of CMAP terms this atom controls
  cmapcomm* MC;   // The list of CMAP terms
};
typedef struct CorrectionMapList cmaplist;

//-----------------------------------------------------------------------------
// CMAPDEFINITION: a structure for holding a CMAP surface.
//-----------------------------------------------------------------------------
struct CMAPDefinition {
  char* note;     // A note taken from the topology file
  dmat esrf;      // The energy surface itself
  dmat dphi;      // The partial derivative with respect to phi
  dmat dpsi;      // The partial derivative with respect to phi
  dmat d2pp;      // The mixed partial derivative with respect to phi and psi
  dmat Ldphi;     // Equivalent of dphi obtained by linear, exact solution
  dmat Ldpsi;     // Equivalent of dpsi obtained by linear, exact solution
  dmat Ld2pp;     // Equivalent of d2pp obtained by linear, exact solution
  dmat stencil;   // Stencil for showing how the points of esrf construct dphi
                  //   or dpsi.  The outer product of the stencil with itself
                  //   shows how to construct the dphi/dpsi mixed partial
                  //   derivative.
  double muePhi;  // Error analysis between linear (exact) and iterative
  double muePsi;  //   (approximate) solutions
  double mue2PP;  //
  double mueFPhi; // Error analysis between linear (exact) and formula
  double mueFPsi; //   (eact, but for numerical precision) solutions
  double mueF2PP; //
};
typedef struct CMAPDefinition cmapdef;

//-----------------------------------------------------------------------------
// CONSTRAINTCOMMAND: command to implement constraints, whether they be 
//                    to use SETTLE for a rigid water molecule, to use  
//                    RATTLE to constrain a bond, or to use P-RATTLE to 
//                    constrain a complex rigid or semi-rigid structure.
//-----------------------------------------------------------------------------
struct ConstraintCommand {
  int exe;     // Flag to execute SETTLE or other constraints with this atom 
               //   as the base
  int* blist;  // Atom identification number of the other atoms in the group
};
typedef struct ConstraintCommand cnstcomm;

//-----------------------------------------------------------------------------
// SETTLEPARAMETERS: three parameters that can be conveniently pre-computed
//                   for SETTLE rigid water constraints.       
//-----------------------------------------------------------------------------
struct SettleParameters {
  double ra;
  double rb;
  double rc;
  double mO;
  double mH;
};
typedef struct SettleParameters settleparm;

//-----------------------------------------------------------------------------
// EXTRAPOINT: this structure identifies a massless, "extra point" particle
//             to distinguish it from the rest of the system.  Frame atoms,
//             frame style, and other features that show what needs to be done.
//-----------------------------------------------------------------------------
struct ExtraPoint {
  int atomid;          // The atom number in the regular topology
  int frstyle;         // The frame style:
                       //   0: extra point on TIP4P-style water
                       //   1: extra point on a line between two frame atoms
                       //  (Other frame types have yet to be included)
  int nfratm;          // The number of atoms in the frame
  int fr1;             // The first frame atom
  int fr2;             // The second frame atom
  int fr3;             // The (optional) third frame atom
  int fr4;             // The (optional) fourth frame atom 
  int atm11;           // Identifier to indicate whether this extra point lies
                       //   1:1 to fr1 (atm11 = 1), 1:1 to fr2 (atm11 += 2),
                       //   1:1 to fr3 (atm11 += 4), and / or 1:1 to fr4
                       //   (atm11 += 8)
  double d1;           // The parameters d1, d2, d3, and d4 mean different
  double d2;           //   things depending on the frame type.  Not all of
  double d3;           //   these parameters are used in all frame types.
                       //   In frame type 4, for instance, d2 refers to an
                       //   angle.
};
typedef struct ExtraPoint expt;

//-----------------------------------------------------------------------------
// EXTRAPOINTRULE: stores a rule for dealing with certain classes of extra
//                 points; rules specified by the user override built-in rules
//                 for objects such as four-point waters.
//-----------------------------------------------------------------------------
struct ExtraPointRule {
  int nfratm;          // The number of frame atoms
  int frstyle;         // The frame style (choose from the built-in options)
  int LJIdx;           // The Lennard-Jones type
  int Join;            // To make the rule mesh with the Join prmtop array
                       //   (see the ambprmtop struct, below)
  int Rotat;           // To make the rule mesh with the Rotat prmtop array
                       //   (see the ambprmtop struct, below)
  int excl2;           // Flags to indicate that frame atoms 2, 3, or 4 should 
  int excl3;           //   also be treated as 1:1 to the extra point when 
  int excl4;           //   computing its exclusions
  int modElec;         // Flag to modify electrostatic properties (default 0,
                       //   do not modify electrostatic properties) 
  int modLJ;           // Flag to modify Lennard-Jones properties (default 0,
                       //   do not modify Lennard-Jones properties)
  char epname[5];      // The name of the atom that is an extra point
  char fr1[5];         // The name of the first frame atom
  char fr2[5];         // The name of the second frame atom
  char fr3[5];         // The name of the (optional) third frame atom
  char fr4[5];         // The name of the (optional) fourth frame atom
  char resname[5];     // The residue name to which this rule applies
  double d1;           // The parameters d1, d2, and d3 mean different
  double d2;           //   things depending on the frame type.  Not all of
  double d3;           //   these parameters are used in all frame types.
  double Charge;       // The charge on the extra point
  double sig;          // The Lennard-Jones sigma value of the "extra point."
                       //   This field (as well as eps) is only used if the
                       //   frame type is zero, meaning that the rule modifies
                       //   the nonbonded properties of an atom
  double eps;          // The Lennard-Jones epsilon value of the "extra point."
};
typedef struct ExtraPointRule eprule;

//-----------------------------------------------------------------------------
// CONNECTORMATRIX: stores lists of atoms that are 1:1, 1:2, 1:3, and 1:4
//                  bonded to an atom.                              
//-----------------------------------------------------------------------------
struct ConnectorMatrix {
  int n11;   // The current number of 1:1 interactions
  int n12;   // The current number of 1:2 interactions
  int n13;   // The current number of 1:3 interactions
  int n14;   // The current number of 1:4 interactions
  int mx11;  // The maximum storable number of 1:1 interactions
  int mx12;  // The maximum storable number of 1:2 interactions
  int mx13;  // The maximum storable number of 1:3 interactions
  int mx14;  // The maximum storable number of 1:4 interactions
  int* L11;  // The list of 1:1 interactions
  int* L12;  // The list of 1:2 interactions
  int* L13;  // The list of 1:3 interactions
  int* L14;  // The list of 1:4 interactions
};
typedef struct ConnectorMatrix map1234;

//-----------------------------------------------------------------------------
// LISTOFGROUP: a list of grouped atoms, such as atoms linked by bonds. 
//-----------------------------------------------------------------------------
struct ListOfGroup {
  int natom;        // The number of atoms in this group
  int* atoms;       // The topology numbers of the atoms in this group
};
typedef struct ListOfGroup lgrp;

//-----------------------------------------------------------------------------
// GROUPMAP: a map of a group of atoms, detailing names and bonds.  This struct
//           cannot really stand alone; examples of it are only created and
//           destroyed for the purposes of comparing lgrp structs.
//-----------------------------------------------------------------------------
struct GroupMap {
  int natom;
  int* BondID;      // Internal ID numbers of atoms participating in all bonds
                    // that hold the atom group together
  char* AtomNames;  // Names of all atoms in the group
  char* ResNames;   // Residues names of all atoms in the group
};
typedef struct GroupMap grpmap;

//-----------------------------------------------------------------------------
// NB14PAIR: a struct for storing a pair of atom ID numbers which share a 1:1,
//           1:2, 1:3, or 1:4 interaction, and some indication of the nature of
//           the interaction.  If extra points are added to a topology at
//           run-time by mdgx, the extra points will imply new interactions,
//           but some of these will need to be nullified or scaled down.
//-----------------------------------------------------------------------------
struct NB14Pair {
  int atmX;       // The first atom in the pair, set to the negative of the
                  //   atom ID number if a vdW interaction must be evaluated
  int atmY;       // The second atom in the pair, set to the negative of the
                  //   atom ID number if electrostatics must be evaluated
};
typedef struct NB14Pair nixpr;

//-----------------------------------------------------------------------------
// AUXILIARYELIMINATION: a list of auxiliary 1:1, 1:2, 1:3, and 1:4
//                       interactions that must be eliminated if they are not
//                       skipped altogether, but are not implied by bonds,
//                       angles, or dihedrals.        
//-----------------------------------------------------------------------------
struct AuxiliaryElimination {
  int n11;
  nixpr* list11;
  int n12;
  nixpr* list12;
  int n13;
  nixpr* list13;
  int n14;
  nixpr* list14;
};
typedef struct AuxiliaryElimination auxelim;

//-----------------------------------------------------------------------------
// CONSTRAINTGRAPH: struct for making graphs of constraint groups; the struct
//                  stores ordered lists of all constrained bonds as well as
//                  mappings between the lists and reservations to indicate
//                  whether certain bonds are already in use as part of this
//                  or another constraint group.
//-----------------------------------------------------------------------------
struct ConstraintGraph {
  int nratl;          // The number of bonds in the system subject to RATTLE
  int grpsize;        // The number of bonds in this particular constraint
                      //   group subject to RATTLE
  int* openend;       // Flags 0 (false) or 1 (true) indicating whether the
                      //   RATTLE'd bonds have been checked for possible
                      //   connections to other RATTLE'd bonds
  int* nodeorder;     // The order of all nodes in the constraint group; a list
                      //   of edges (bonds between atoms) is first compiled,
                      //   then all unique nodes (atoms) are found and the
                      //   number of edges leading to each nodes is counted
  int* rsrv;          // Reservations 0 (false) or 1 (true) indicating whether
                      //   edges have already been assigned to this group or
                      //   some other one
  int* map12;         // Map indicating that edge A as found in list ratlbond
                      //   is the same as edge B in list ratlbond2, where
                      //   map12[A] = B
  int* map21;         // Map indicating that edge A as found in list ratlbond2
                      //   is the same as edge B in list ratlbond, where
                      //   map21[A] = B
  nixpr* ratlbond;    // 
  nixpr* ratlbond2;
  nixpr* ratlgrp;
};
typedef struct ConstraintGraph cnstgrp;

//-----------------------------------------------------------------------------
// AMBPRMTOP: a structure for AMBER 7 topology file data.               
//-----------------------------------------------------------------------------
struct ambprmtop {

  // Integers
  int natom;    // The number of atoms
  int nres;     // The number of residues
  int ntypes;   // Total number of distinct atom types
  int nhparm;   // Currently not used
  int nparm;    // Currently not used
  int tnexcl;   // The number of excluded atoms
  int natyp;    // The number of atom types in the parameter file (see Solty)
  int nphb;     // The number of distinct 10-12 interaction types
  int ifpert;   // Set to 1 if perturbation information is to be read
  int ifbox;    // Set to 1 if standard periodic box, 2 if truncated octahedron
  int nmxrs;    // Number of atoms in the largest residue
  int ifcap;    // Set to 1 if the CAP option from edit was specified
  int iptres;   // The final residues that is considered part of the solute
  int iptatom;  // The final atom that is considered part of the solute
  int nspm;     // The total number of molecules
  int nspsol;   // The first solvent molecule
  int natcap;   // The last atom before the start of cap waters placed by edit
  int blank;    // Not sure what this is
  int rattle;   // Flag to activate RATTLE
  int settle;   // Flag to activate SETTLE
  int nwat;     // The number of rigid three point water molecules
  int ncnst;    // The number of constraints in the system
  int ndf;      // The number of degrees of freedom in the system
  int nprtcl;   // The number of particles in the system (each constrained
                //   group of atoms counts as one particle, and each atom not
                //   in a constrained group, counts as one particle)
  int ljbuck;   // Flag to signal whether a Lennard-Jones 12-6 or Buckingham
                //   exp-6 potential is to be used
  int qshape;   // Flag to signal whether point charges or spherical Gaussian-
                //   distributed charges are to be used
  int neprule;  // The number of extra point rules, obtained from an auxiliary
                //   file (the name of the file is stored in the eprulesource
                //   field)
  int nxtrapt;  // The number of extra points, as counted by the extra point
                //   rule implementation routines
  int maxtrapt; // The maximum number of extra points that can be stored; kept
                //   to indicate when the extra points array must be expanded
  int numextra; // The number of extra points, as stated in the preamble to the
                //   topology file
  int ncopy;    // The number of Path Integral Molecular Dynamics slices /
                //   number of beads
  int ngrp;     // The number of bonded atom groups
  int ipol;     // Flag to indicate that polarizability is present

  int EPInserted;    // Flag to indicate that extra points have been added to
                     //   this topology
  int norigatom;     // The original number of atoms, relevant if new extra
                     //   points are inserted
  int nCustomVS;     // The number of custom virtual site frame definitions in
                     //   the topology
  int nclingrule;    // The number of cling rules defining how atoms react to
                     //   a grid-based restraint potential
  int RattleGrpMax;  // The maximum number of atoms that could be in any of the
                     //   RATTLE-constrained groups in this topology
  int AtomsNoElec;   // Flags to indicate that this topology contains atoms
                     //   which have no electrostatic or Lennard-Jones
  int AtomsNoLJ;     //   properties, respectively
  int ExclMarked;    // Flag to indicate that exclusions must be eliminated
                     //   from non-bonded interactions apart from what bond,
                     //   angle, and dihedral routines would already exclude
  int ncmap;         // The number of CMAP terms present.  Default 0, again
                     //   indicating that no CMAP is present.  This cannot be
                     //   part of the bah struct as that defines the number of
                     //   unique bond / angle / dihedral parameters, whereas
                     //   this describes the number of atom groups actually
                     //   participating in CMAP interactions.  It is the length
                     //   of the array Mterms.
  int ncmapsurf;     // The number of unique CMAP surfaces.  This could be part
                     //   of the bah struct.

  // Reals
  double cutcap;    // The distance from the center of the cap to the outside
  double xcap;      // X coordinate for the center of the cap
  double ycap;      // Y coordinate for the center of the cap
  double zcap;      // Z coordinate for the center of the cap
  double lj14fac;   // The van-der Waals 1:4 interaction adjustment factor
  double elec14fac; // The electrostatic 1:4 interaction adjustment factor
  double TotalMass; // The total mass of the system (atomic mass units) 
  double initq;     // The initial total charge on the system (for output
                    //   record keeping--system charges are typically rounded
                    //   to the nearest integer, often zero, by adding a small
                    //   amount of charge to every atom)

  // In what follows, (b, a, h) means "(bonds, angles, and dihedrals)"
  bah withH;    // Number of (b, a, h) with hydrogen
  bah woH;      // Number of (b, a, h) without hydrogen
  bah woHC;     // Number of (b, a, h) without hydrogen plus constraints
  bah nBAH;     // Number of unique (b, a, h) types
  bah pert;     // Number of (b, a, h) to be perturbed
  bah wpert;    // Number of (b, a, h) with all atoms in the perturbed group

  // Arrays of elemental types
  int* ZNumber;        // The atomic Z-number for each atom (set to zero if not
                       //   present in the topology file)
  int* LJIdx;          // The atom type index for each atom
  int* NExcl;          // The number of excluded atoms for each atom
  int* ConExcl;        // Start positions in the excluded atom list
  int* ExclList;       // The excluded atom list
  int* NBParmIdx;      // Nonbonded parameter index, ntypes*ntypes elements
  int* ResLims;        // The residue limits
  int* Join;           // Tree joining information, used by ancient programs
  int* Rotat;          // No longer in use
  int* Nsp;            // The total number of atoms in each molecule
  int* OldAtomNum;     // Old atom number sequence (initially not allocated,
                       //   only allocated if EPInserted is 1); this array
                       //   stores the number of the atom in the topology as
                       //   initially read from a file, before any extra points
                       //   have been inserted, or -1 for extra points that
                       //   were not part of the original topology
  int* MobileAtoms;    // Identifies mobile atoms in the system
  int* customVSi;      // Custom virtual site frame indexing
  double* Charges;     // The atomic charges
  double* OrigCharges; // Copy of the original atomic charges (stored for
                       //   reference when re-writing the topology)
  double* Masses;      // The atomic masses (atomic mass units)
  double* InvMasses;   // The square root inverse atomic masses
  double* BondK;       // The bond force constants (kcal/mol-A^2)
  double* BondEq;      // The equilibrium bond lengths (Angstroms)
  double* BondPullK;   // Force constants (kcal/mol-A^2) and equilibrium bond
  double* BondPullEq;  //   lengths (Angstroms) that kick in when bonds are
                       //   pulled farther than usual
  double* BondPressK;  // Force constants (kcal/mol-A^2) and equilibrium bond
  double* BondPressEq; //   lengths (Angstroms) that kick in when bonds are
                       //   compressed farther than usual
  double* AnglK;       // The angle force constants (kcal/mol-rad^2)
  double* AnglEq;      // The equilibrium angle values (radians)
  double* DiheK;       // The dihedral force constants (kcal/mol)
  double* DiheN;       // The dihedral periodicities
  double* DihePhi;     // The dihedral phase
  double* scee;        // The 1-4 electrostatic adjustments for all dihedrals
  double* scnb;        // The 1-4 van-der Waals adjustments for all dihedrals
  double* solty;       // The SOLTY field is reserved for future use
  double* LJA;         // The Lennard-Jones A parameters
  double* LJB;         // The Lennard-Jones B parameters
  double* LJC;         // The "Lennard-Jones" C parameters (if these are given,
                       //   it means that the van-der Waals model is a modified
                       //   Buckingham potential, and thus the Lennard-Jones A
                       //   and B parameters are really Buckingham A and B
                       //   parameters)
  double* lVDWc;       // Long-Range van-der Waals corrections, stored by type,
                       //   to be divided by the instantaneous unit cell volume
                       //   when implemented
  double* SolA;        // The HBond A parameters
  double* SolB;        // The HBond B parameters
  double* HBCut;       // No longer in use
  double* Radii;       // Atomic PB radii
  double* Screen;      // Atomic screening terms for Generalized Born
  double* customVSd;   // Custom virtual site frame details
  char RadSet[128];    // Radius set name
  char vstamp[128];    // Unknown         
  char WaterName[8];   // The name of water molecules (for implementing SETTLE)
  char* AtomNames;     // The atom names
  char* ResNames;      // The residue names
  char* AtomTypes;     // The atom type names
  char* TreeSymbols;   // The atom tree chain symbols
  char* rattlemask;    // Mask to specially request RATTLE bond constraints
  char* norattlemask;  // Mask to specially omit RATTLE bond constraints
  double gdim[6];      // The box dimensions

  // Arrays of structures
  bond* BIncH;       // Bonds including Hydrogen
  bond* BNoH;        // Bonds not including Hydrogen
  angle* AIncH;      // Angles including Hydrogen
  angle* ANoH;       // Angles not including Hydrogen
  dihedral* HIncH;   // Dihedrals including Hydrogen
  dihedral* HNoH;    // Dihedrals not including Hydrogen
  cmapterm* Mterms;  // All CMAP terms
  bondlist* BLC;     // The list of bonds controlled by each atom
  angllist* ALC;     // The list of angles controlled by each atom
  dihelist* HLC;     // The list of dihedral angles controlled by each atom
  cmaplist* MLC;     // The list of CMAP terms controlled by each atom
  bonddef* BParam;   // The re-organized, unified bond parameter array
  angldef* AParam;   // The re-organized, unified angle parameter array
  dihedef* HParam;   // The re-organized, unified dihedral parameter array
  cmapdef* MParam;   // The re-organized, unified list of CMAPs
  cnstcomm* SHL;     // Constraint commands for each individual atom (SETTLE is
                     //   performed starting with the oxygens)
  eprule* eprules;   // List of rules for extra points
  expt* xtrapts;     // List of extra points (field nxtrapts holds the size)
  lgrp* lgrps;       // List of atom groups in this topology, defining
                     //   separate molecules within the simulation
  map1234* nb1234;   // A connectivity array detailing 1:1, 1:2, 1:3, and 1:4
                     //   interactions which is reflexive and symmetric
  auxelim* ElimPair; // A list of 1:1, 1:2, and 1:3 bonded or 1:4 nonbonded
                     //   interactions not involving dihedrals that must be
                     //   calculated based on added extra points
  lgrp* FR1Idx;      // Index into extra point array xtrapts (each atom has 
                     //   is own element of FR1Idx, but it points to no 
                     //   elements of the xtrapts array unless the atom is
                     //   frame atom 1 of one or more extra points).  When
                     //   atoms go into cells, extra points are ultimately
                     //   assigned after their frame atoms, and frame atom 1,
                     //   to which the extra point is attached, is the first
                     //   thing to look for.

  // Lennard-Jones pre-computed parameter tables
  dmat LJftab;                 // Forces
  dmat LJutab;                 // Energies

  // SETTLE pre-computed parameters
  settleparm FWtab;            

  // Source information
  char source[MAXNAME];        // Topology source file
  char eprulesource[MAXNAME];  // Extra points rules source file
};
typedef struct ambprmtop prmtop;

//-----------------------------------------------------------------------------
// UniqueBondList: this struct, similar in form to the atom group struct but
//                 given a new name to avoid confusion in the field names,
//                 stores (for a particular atom) a total count and the ID
//                 numbers of bonds which an controls.  The ID number of a
//                 bond is the its number in the list of bonds controlled by
//                 an atom--see the BLC field of the ambprmtop struct above.
//-----------------------------------------------------------------------------
struct UniqueBondList {
  int nitem;
  int* items;
};
typedef struct UniqueBondList ublist;

//-----------------------------------------------------------------------------
// PrmtopCorrespondence: stores the correspondence between two topology 
//                       structs, in terms of the atom-to-atom lineup and   
//                       global matching style.                     
//-----------------------------------------------------------------------------
struct PrmtopCorrespondence {
  int relate;           // How these topologies correspond generally...
                        //   0: atom N of topology A corresponds to atom N of
                        //      topology B
                        //   1: if atoms N and N+1 of topology A correspond to
                        //      anything in topology B, they will correspond to
                        //      atoms K and K+P of topology B, where P > 0
                        //   2: atoms of topology A may or may not correspond
                        //      atoms of topology B, but there is no guarantee
                        //      of ordering
  int vdw;              // Flag to indicate whether any van-der Waals
                        //   interactions change between the corresponding
                        //   atoms of topologies A and B
  int vdw14;            // Flag to indicate whether any 1:4 van-der Waals
                        //   interactions change between the corresponding
                        //   atoms of topologies A and B
  int elec;             // Flag to indicate whether any charges change between
                        //   the corresponding atoms of topologies A and B
  int elec14;           // Flag to indicate whether any 1:4 electrostatics
                        //   change between corresponding atoms of topologies
                        //   A and B
  int bond;             // Flag to indicate whether any bond terms change
                        //   between corresponding atoms of topologies A and B
                        //   (also set to 1 if RATTLE is toggled)
  int angl;             // Flag to indicate whether any angle terms change
                        //   between corresponding atoms of topologies A and B
  int dihe;             // Flag to indicate whether any dihedral terms change
                        //   between corresponding atoms of topologies A and B
  int appear;           // Flag to indicate that topology B has atoms that
                        //   have appeared amidst the contents of topology A
  int disappear;        // Flag to indicate that topology B has atoms that
                        //   have appeared from the contents of topology A
  int nxdf;             // The number of degrees of freedom in the extended
                        //   system which is the union of particles in A and B
  int nxprtcl;          // The number of particles in the extended system  
  int uniA;             // Number of unique atoms in system A
  int uniB;             // Number of unique atoms in system B
  int comAB;            // Number of common atoms between systems A and B
  int* matchA;          // The correspondences of atoms in topology A to B
  int* matchB;          // The correspondences of atoms in topology B to A
  int* corrA;           // The correspondences of atoms in topology A to B
  int* corrB;           // The correspondences of atoms in topology B to A
  double* dQ;
  ublist* SubBondA;     // Bonds that must be subtracted from topology A to
                        //   get topology B (bonds unique to topology A)
  ublist* AddBondB;     // Bonds unique to topology B
  ublist* SubAnglA;     // Angles that must be subtracted from topology A to
                        //   get topology B (angles unique to topology A)
  ublist* AddAnglB;     // Angles unique to topology B
  ublist* SubDiheA;     // Dihedrals that must be subtracted from topology A to
                        //   get topology B (dihedrals unique to topology A)
  ublist* AddDiheB;     // Dihedrals unique to topology B
  prmtop *tpA;          // Pointers to topologies A and B are stored for
  prmtop *tpB;          //   reference, so things like the number of atoms in
                        //   each topology is known
};
typedef struct PrmtopCorrespondence prmcorr;
#endif
