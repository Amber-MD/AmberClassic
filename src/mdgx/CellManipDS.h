#ifndef CellManipStructs
#define CellManipStructs

#ifndef PREP_API
#include "MatrixDS.h"
#include "BSplineDS.h"
#include "MPIMapDS.h"
#include "GridDS.h"
#endif

struct AtomInCell {
  int id;
  int lj;
  double q;
  double loc[3];
  double frc[3];
};
typedef struct AtomInCell atomc;

struct AtomBuffer {
  int id;           // The atom ID number, its index in the topology
  int dreg;         // The destination region of this atom
  double loc[3];    // The location of the atom in the cell
};
typedef struct AtomBuffer atomb;

struct AtomBufferPlusVelocity {
  int id;           // The atom ID number, its index in the topology
  int dreg;         // The destination region of this atom
  double loc[3];    // The location of the atom in the cell
  double vel[3];    // The velocity of the atom (information in this field
                    //   will go straight to the associated coord struct)
};
typedef struct AtomBufferPlusVelocity atombv;

struct AtomBufferPlusAllInfo {
  int id;           // The atom ID number, its index in the topology
  int dreg;         // The destination region of this atom
  double loc[3];    // The location of the atom in the cell
  double vel[3];    // The velocity of the atom (information in this field
                    //   will go straight to the associated coord struct)
  double sysloc[3]; // The location of the atom in the system
  double prvloc[3]; // The previous location of the atom
  double prvvel[3]; // The previous velocity of the atom
  double prvfrc[3]; // The previous force on the atom
};
typedef struct AtomBufferPlusAllInfo atombx;

struct RangeBuffer {
  double dx;      // The displacement in x
  double dy;      // The displacement in y
  double dz;      // The displacement in z
  double r2;      // The squared distance
};
typedef struct RangeBuffer rngbuff;

struct CellBlock {
  int nexp;            // The number of atoms the cell is exporting
  int nimp;            // The number of atoms the cell has imported
  int maxatom;         // The maximum number of atoms that this cell can hold
  int CGRank;          // The rank of the process responsible for this cell
                       //   within the cell grid private communicator
#ifdef MPI
  int nbond;           // The number of bonds
  int nangl;           // The number of angles
  int ndihe;           // The number of dihedrals
  int ncmap;           // The number of CMAPs
  int nelim;           // The number of eliminations
#endif
  int pmordr[3];       // Reciprocal space particle -> mesh interpolation order
                       //   (stored in the cell for convenience; otherwise it
                       //   becomes a headache to keep referencing that
                       //   information from the reccon struct)
  int gbin[4];         // The location of the cell in the cell grid (X, Y, and
                       //   Z position followed by absolute number in the
                       //   range 0... nX*nY*nZ)
  int nFscr;           // The number of scoring grids maintained by this cell
                       //   (initialized to zero by default, and only set to
                       //   nonzero values by restraint module functions)
  int* nr;             // Total number of atoms in each sector
  int* nsr;            // Total number of atoms in each sub-sector list
  int* ljIDbuff;       // Buffer array for CELL ID numbers of atoms that have
                       //   passed all the necessary cutoff tests; corresponds
                       //   to r2 values stored in the ljr2buff array
  int* qIDbuff;        // Buffer array for CELL ID numbers of atoms that have
                       //   passed the electrostatic cutoff test; corresponds
                       //   to r2 values stored in the qr2buff array
  int* GPSptr;         // Pointer into parent cell grid's AtmGPS 2D array
  rngbuff* ljr2buff;   // Buffer array of r2 values for LJ interactions
  rngbuff* qr2buff;    // Buffer array of r2 values for Q-Q interactions
  imat ordr;           // Ordering array for direct space interactions
  imat supordr;        // "Super-ordering" array for direct space interactions,
                       //   holds additional sorting after ordr is computed
  double orig[3];      // Coordinate origin of the cell's primary sector
  double midp[3];      // Coordinate origin of the cell's midpoint

  /*** Auxiliary arrays to store data that is not always needed ***/
  bcof* xcof;           // xcof, ycof, and zcof store B-spline coefficients for
  bcof* ycof;           //   mapping atoms in this cell to a mesh.
  bcof* zcof;           //
  atomc* atmscr;        // Scratch space for merge sort

  /*** Atom information used during nonbonded loop.  The elements of the  ***/
  /*** map array point to regions of data, similar to a matrix structure. ***/
  atomc* data;
  atomc** map;

  /*** Restraint grid information; a series of  ***/
  fbook* Fscr;

  /*** Buffer arrays ***/
  atomb* import;
  atomb* pexport;
  atombv* Vimport;
  atombv* Vexport;
  atombx* Ximport;
  atombx* Xexport;
};
typedef struct CellBlock cell;

struct CellBlockGrid {
  int nsend;            // The maximum number of sends or receives for which
  int nrecv;            //   this cell grid is prepared
  int MyCellCount;      // Number of cells controlled by a particular process
  int sysID;            // The number of the system to which this grid pertains
  int tid;              // The rank of this particular process in the cell
                        //   grid's private communicator  
  int nthreads;         // The number of (MPI) threads assigned to this cell
                        //   grid, or 1 if no MPI
  int MasterHalfLoad;   // Half the master process's cell load
  int* nexp;            // The number of atoms that this cell grid holds in its
                        //   its export buffers
  int* maxexp;          // The maximum number of atoms that may be sent or
  int* maximp;          //   received in any of the expected messages
  int* MyCellDomain;    // Cells controlled by this particular process
  int* CrdPoolSize;     // 
  atomb** import;       // Pooled import and export buffers for all cells,
  atomb** pexport;      //   process-to-process communication involving atomb,
  atombv** Vimport;     //   atombv, and atombx types, similar to their
  atombv** Vexport;     //   counterparts in the individual cell structs
  atombx** Ximport;     //
  atombx** Xexport;     //
  atombx** CrdPool;     // Buffer for pooling coordinates to the master process
  dcplan DirCommPlan;   // Direct space communication plan
  int maxatom;          // The maximum number of atoms that any one cell may
                        //   hold
  int ncell;            // The number of cells in this grid
  int ng[3];            // The dimensions of the cell block grid
  double dbng[3];       // Dimensions of the cell block grid recast as doubles
  double celldim[7];    // celldim contains the lengths, in Angstroms, of the
                        //   edges of each cell, then the length of the
                        //   maximum direct space cutoff (Mcut copied over from
                        //   a dircon struct found in pmeDirectDS.h), and
                        //   finally (in the last three elements) the location
                        //   of the origin of the cell's central region in
                        //   units of the cell lengths
  cell* data;           // The linear array of cells in this grid
  cell*** map;          // Map to cells by x / y / z indices
#ifdef MPI
  MPI_Comm dspcomm;     // Private communicator for this cell grid
#endif
  gsplc* MeshCommPlan;  // Mesh communication plans (allocates nthreads
                        //   on the master process, 1 on every other
                        //   process)
  imat AtmGPS;          // Lists of cell data indices where atoms ordered
                        //   according to the master topology file can be found
};
typedef struct CellBlockGrid cellgrid;

#endif
