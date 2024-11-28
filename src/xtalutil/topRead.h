#include "stringDefs.h"

#ifndef TOPREAD_STRUCTS
#define TOPREAD_STRUCTS

/***=======================================================================***/
/*** BOXDIMENSIONS: the dimensions of the periodic box.  This information  ***/
/***                is NOT updated during MD runs--updates to the X, Y,    ***/
/***                and Z dimensios are stored elsewhere during the        ***/
/***                simulations.                                           ***/
/***=======================================================================***/
struct boxdimensions {
  double x;       // Length in X
  double y;       // Length in Y
  double z;       // Length in Z
  double alpha;   // Angle between XY and XZ planes
  double beta;    // Angle betweer XY and YZ planes
  double gamma;   // Angle between XZ and YZ planes
};
typedef struct boxdimensions boxdim;

/***=======================================================================***/
/*** BNDANGDIHE: a group of three integers, bonds, angles, and dihedrals.  ***/
/***=======================================================================***/
struct bndangdihe {
  int nbond;
  int nangl;
  int ndihe;
};
typedef struct bndangdihe bah;

/***=======================================================================***/
/*** BONDIDX: a bond index structure.                                      ***/
/***=======================================================================***/
struct bondidx {
  int a;
  int b;
  int idx;
};
typedef struct bondidx bond;

/***=======================================================================***/
/*** ANGLIDX: an angle index structure.                                    ***/
/***=======================================================================***/
struct anglidx {
  int a;
  int b;
  int c;
  int idx;
};
typedef struct anglidx angle;

/***=======================================================================***/
/*** DIHEIDX: a dihedral index structure.                                  ***/
/***=======================================================================***/
struct diheidx {
  int a;
  int b;
  int c;
  int d;
  int idx;
};
typedef struct diheidx dihedral;

/***=======================================================================***/
/*** AMBPRMTOP: a structure for AMBER 7 topology file data.                ***/
/***=======================================================================***/
struct ambprmtop {

  /*** Integers ***/
  int natom;   // The number of atoms
  int nres;    // The number of residues
  int ntypes;  // Total number of distinct atom types
  int nhparm;  // Currently not used
  int nparm;   // Currently not used
  int tnexcl;  // The number of excluded atoms
  int natyp;   // The number of atom types in the parameter file (see Solty)
  int nphb;    // The number of distinct 10-12 interaction types
  int ifpert;  // Set to 1 if perturbation information is to be read
  int ifbox;   // Set to 1 if standard periodic box, 2 if truncated octahedron
  int nmxrs;   // Number of atoms in the largest residue
  int ifcap;   // Set to 1 if the CAP option from edit was specified
  int iptres;  // The final residues that is considered part of the solute
  int iptatom; // The final atom that is considered part of the solute
  int nspm;    // The total number of molecules
  int nspsol;  // The first solvent molecule
  int natcap;  // The last atom before the start of cap waters placed by edit
  int blank;   // Not sure what this is

  /*** Reals ***/
  double cutcap; // The distance from the center of the cap to the outside
  double xcap;   // X coordinate for the center of the cap
  double ycap;   // Y coordinate for the center of the cap
  double zcap;   // Z coordinate for the center of the cap

  // In what follows, (b, a, h) means "(bonds, angles, and dihedrals)"
  bah withH;    // Number of (b, a, h) with hydrogen
  bah woH;      // Number of (b, a, h) without hydrogen
  bah woHC;     // Number of (b, a, h) without hydrogen plus constraints
  bah nBAH;     // Number of unique (b, a, h) types
  bah pert;     // Number of (b, a, h) to be perturbed
  bah wpert;    // Number of (b, a, h) with all atoms in the perturbed group

  // Arrays of elemental types
  int* LJIdx;        // The atom type index
  int* NExcl;        // The number of excluded atoms for each atom
  int* ConExcl;      // Start positions in the excluded atom list
  int* ExclList;     // The excluded atom list
  int* NBParmIdx;    // Nonbonded parameter index, ntypes*ntypes elements
  int* ResLims;      // The residue limits
  int* Solty;        // Currently unused (reserved for future use)
  int* Join;         // Tree joining information, used by ancient programs
  int* Rotat;        // No longer in use
  int* Nsp;          // The total number of atoms in each molecule
  int* IAPert;       // 1 if the atom is perturbed, 0 if not
  int* LJidxPert;    // The atom type index at lambda = 0
  double* Charges;   // The atomic charges
  double* Charges0;  // The atomic charges at lambda = 0
  double* Masses;    // The atomic masses
  double* BondK;     // The bond force constants (kcal/mol-A^2)
  double* BondEq;    // The equilibrium bond lengths (Angstroms)
  double* AnglK;     // The angle force constants (kcal/mol-rad^2)
  double* AnglEq;    // The equilibrium angle values (radians)
  double* DiheK;     // The dihedral force constants (kcal/mol)
  double* DiheN;     // The dihedral periodicities
  double* DihePhi;   // The dihedral phase
  double* LJA;       // The Lennard-Jones A parameters
  double* LJB;       // The Lennard-Jones B parameters
  double* SolA;      // The HBond A parameters
  double* SolB;      // The HBond B parameters
  double* HBCut;     // No longer in use
  double* AtomPol;   // Atomic polarizabilities
  double* AtomPol0;  // Atomic polarizabilities at lambda = 0
  double* Radii;     // Atomic radii
  double* Screen;    // Atomic screening terms (?)
  char RadSet[128];   // Radius set name
  char vstamp[128];   // Radius set name
  char* AtomNames;   // The atom names
  char* ResNames;    // The residue names
  char* AtomTypes;   // The atom type names
  char* AtomNamesP0; // The atom names
  char* ResNamesP0;  // The residue names
  char* AtomTypesP0; // The atom type names
  char* TreeSymbols; // The atom tree chain symbols

  /*** Arrays of structures ***/
  bond* BIncH;       // Bonds including Hydrogen
  bond* BNoH;        // Bonds not including Hydrogen
  bond* BPert1;      // Perturbed bonds with lambda = 1
  bond* BPert0;      // Perturbed bonds with lambda = 0
  angle* AIncH;      // Angles including Hydrogen
  angle* ANoH;       // Angles not including Hydrogen
  angle* APert1;     // Perturbed angles with lambda = 1
  angle* APert0;     // Perturbed angles with lambda = 0
  dihedral* HIncH;   // Dihedrals including Hydrogen
  dihedral* HNoH;    // Dihedrals not including Hydrogen
  dihedral* HPert1;  // Perturbed dihedrals with lambda = 1
  dihedral* HPert0;  // Perturbed dihedrals with lambda = 0
  boxdim smbx;       // The periodic box information

  /*** Source information ***/
  char source[MAXNAME];
};
typedef struct ambprmtop prmtop;

/***=======================================================================***/
/*** COORDINATES: a set of coordinates from an AMBER restart or .crd file. ***/
/***=======================================================================***/
struct coordinates {
  int natom;
  int rst;
  double t;
  double* loc;
  double* vel;
  double box[6];
  char source[MAXNAME];
};
typedef struct coordinates coord;

#endif

#ifndef TOPREAD_FUNCS
#define TOPREAD_FUNCS

void GetPrmTop(prmtop *tp, int verbosity, int adjbnd);

void GetRst(coord *tc, prmtop *tp);

void ReadToFlag(char* fname, FILE *inp);

char* LoadChar(FILE *inp, int N, int P, int Q);

int* LoadInteger(FILE *inp, int N);

double* LoadDouble(FILE *inp, int N);

void MapResidues(prmtop *tp, coord *tc, char* outname);

void AdjustBondArray(int* A, int N, int P);

void PrintBondInfo(prmtop *tp, int* B, int nB, int P, int finst, double* K,
                   double* Eq, double* Pd, FILE *outp, int* cpylist,
		   coord *tc);

void BondScan(prmtop *tp, int* cpylist, int finst, int A, int B, int C,
              int D, coord *tc, double *mval, double *rmsval, int bondorder);

void FindDisulfides(prmtop *tp, coord *tc);

void VacuumBubble(prmtop *tp, coord *tc, double gspc, double bubspc);

int FindAtom(prmtop *tp, int il, int ih, char* aname);

int Chirality(prmtop *tp, coord *tc, int CA, int A, int B, int C, int D);

void ProteinChiralityCheck(prmtop *tp, coord *tc);

void PutPrmTop(prmtop *tp, char* fname, char* title);

#endif
