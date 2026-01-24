#ifndef BroadcastCommandHeadings
#define BroadcastCommandHeadings

#include "pmeDirectDS.h"
#include "pmeRecipDS.h"
#include "TrajectoryDS.h"
#include "ChargeFitDS.h"
#include "IPolQDS.h"

#ifdef MPI
void BroadcastInputData(dircon *dcinp, reccon *rcinp, trajcon *tj, prmtop* tp,
			ipqcon *ipqinp);

void BroadcastCoordinates(coord *crd, trajcon *tj, int sID);

void BroadcastEPInfo(prmtop *tp, trajcon *tj);

int MPI_Valgrind_bcast(void *buffer, int count, MPI_Datatype datatype,
		       int root, MPI_Comm comm, int idebug);
#endif

#endif
