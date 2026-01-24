#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "Macros.h"

#include "TrajectoryDS.h"
#include "pmeDirectDS.h"
#include "CellManipDS.h"

#ifdef MPI
//-----------------------------------------------------------------------------
// DefineMPITypes: this function defines a set of MPI types for data structs
//                 that must be passed around.  It allows mdgx structs to be
//                 passed by MPI without deconstruction.  
//
// Arguments:                                                           
//   tj:      trajectory control data, where the integer pointers of the
//            custom-defined MPI types will be stored                   
//-----------------------------------------------------------------------------
void DefineMPITypes(trajcon *tj)
{
  int count, iextent, dextent, cextent;
  int lengths[64];
  MPI_Aint offsets[64];
  MPI_Datatype mtypes[64];

  // Sizes of integers and doubles
  iextent = sizeof(int);
  dextent = sizeof(double);
  cextent = sizeof(char);

  // Define the pmeDirectControlData type
  count = 2;
  lengths[0] = 1;
  lengths[1] = 10;
  offsets[0] = 0;
  offsets[1] = MAX(iextent, dextent);
  mtypes[0] = MPI_INT;
  mtypes[1] = MPI_DOUBLE;
  MPI_Type_create_struct(count, lengths, offsets, mtypes, &tj->MPI_DIRCON);
  MPI_Type_commit(&tj->MPI_DIRCON);
  if ((int)offsets[1] + 10*dextent != (int)sizeof(dircon)) {
    printf("DefineMPITypes >> Error.  sizeof(dircon) = %d, MPI type has size "
	   "%d.\n", (int)sizeof(dircon), (int)offsets[1] + 10*dextent);
    exit(1);
  }

  // Define the AtomBuffer type
  lengths[0] = 2;
  lengths[1] = 3;
  offsets[1] = MAX(2*iextent, dextent);
  MPI_Type_create_struct(count, lengths, offsets, mtypes, &tj->MPI_ATOMB);
  MPI_Type_commit(&tj->MPI_ATOMB);
  if ((int)offsets[1] + 3*dextent != (int)sizeof(atomb)) {
    printf("DefineMPITypes >> Error.  sizeof(atomb) = %d, MPI type has size "
	   "%d.\n", (int)sizeof(atomb), (int)offsets[1] + 3*dextent);
    exit(1);
  }

  // Define the AtomBufferPlusVelocity type
  lengths[1] = 6;
  MPI_Type_create_struct(count, lengths, offsets, mtypes, &tj->MPI_ATOMV);
  MPI_Type_commit(&tj->MPI_ATOMV);
  if ((int)offsets[1] + 6*dextent != (int)sizeof(atombv)) {
    printf("DefineMPITypes >> Error.  sizeof(atombv) = %d, MPI type has size "
           "%d.\n", (int)sizeof(atombv), (int)offsets[1] + 6*dextent);
    exit(1);
  }

  // Define the AtomBufferPlusAllInfo type
  lengths[1] = 18;
  MPI_Type_create_struct(count, lengths, offsets, mtypes, &tj->MPI_ATOMX);
  MPI_Type_commit(&tj->MPI_ATOMX);
  if ((int)offsets[1] + 18*dextent != (int)sizeof(atombx)) {
    printf("DefineMPITypes >> Error.  sizeof(atombx) = %d, MPI type has size "
           "%d.\n", (int)sizeof(atombx), (int)offsets[1] + 18*dextent);
    exit(1);
  }

  // Define the ExtraPointRule type
  count = 3;
  lengths[0] = 10;
  lengths[1] = 30;
  lengths[2] = 7;
  offsets[1] = MAX(10*iextent, 5*dextent);
  offsets[2] = offsets[1] + MAX(30*cextent, 4*dextent);
  mtypes[1] = MPI_CHAR;
  mtypes[2] = MPI_DOUBLE;
  MPI_Type_create_struct(count, lengths, offsets, mtypes, &tj->MPI_EPRULE);
  MPI_Type_commit(&tj->MPI_EPRULE);
  if ((int)offsets[2] + 7*dextent != (int)sizeof(eprule)) {
    printf("DefineMPITypes >> Error.  sizeof(eprule) = %d, MPI type has size "
           "%d.\n", (int)sizeof(eprule), (int)offsets[2] + 7*dextent);
    exit(1);
  }
}

//-----------------------------------------------------------------------------// FreeMPITypes: destroys custom-defined MPI types to prevent memory leaks.
//
// Arguments:                                                           
//   tj:      trajectory control data, where the integer pointers of the
//            custom-defined MPI types will be stored                   
//-----------------------------------------------------------------------------
void FreeMPITypes(trajcon *tj)
{
  MPI_Type_free(&tj->MPI_DIRCON);
  MPI_Type_free(&tj->MPI_ATOMB);
  MPI_Type_free(&tj->MPI_ATOMV);
  MPI_Type_free(&tj->MPI_ATOMX);
  MPI_Type_free(&tj->MPI_EPRULE);
}
#endif
