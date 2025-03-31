#ifndef CrdManipStructs
#define CrdManipStructs

#ifndef PREP_API
#include "MatrixDS.h"
#endif

struct Coordinates {
  int natom;        // The number of atoms in the system
  int isortho;      // Flag to indicate whether the unit cell is orthorhombic
  int* atmid;       // The ID numbers of atoms in the system; starts out as a
                    //   series of numbers from 0 to natom-1, but may be useful
                    //   for tracking atoms if coord structs are copied and
                    //   rearranged
  double* loc;      // The (current) locations of the atoms
  double* prvloc;   // The previous locations of the atoms
  double* scrloc;   // Scratch space for locations of the atoms
  double* vel;      // The (current) velocities of the atoms
  double* prvvel;   // The previous velocities of the atoms
  double* frc;      // The (current) forces on atoms
  double* prvfrc;   // The previous forces on atoms
  double* scrfrc;   // Scratch array of forces on atoms
  double gdim[6];   // The simulation cell dimensions
  double hgdim[3];  // The simulation cell vector half lengths, stored for
                    //   convenience but crucial to update
  dmat U;           // The transformation matrix for taking real-space
                    //   coordinates into fractional coordinate space
  dmat invU;        // The inverse transformation matrix, invU = U^-1
  dmat fcnorm;      // Normals to each of the pair interaction interfaces; rows
                    //   of this matrix contain the unit normal vectors for
                    //   easy memory access
};
typedef struct Coordinates coord;

#endif
