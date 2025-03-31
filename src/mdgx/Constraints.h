#ifndef ConstraintHeadings
#define ConstraintHeadings

#include "TopologyDS.h"
#include "CrdManipDS.h"
#include "CellManipDS.h"
#include "TrajectoryDS.h"

//-----------------------------------------------------------------------------
// Diagram for SETTLE:                                                  
//                                                                      
// Oxygen at a0, hydrogens at b0, c0, center of mass at the origin.     
//                                                                      
//                   |                                                  
//                   |                                                  
//                   |                                                  
//                a0 x -----                                            
//                   |   |                                              
//                   |   ra                                             
//                   |   |                                              
// ----------------Orig.------------------                              
//           |       |        |                                         
//           rb      |---rc---|                                         
//           |       |        |                                         
//    b0 x-----------|        c0 x                                      
//                                                                      
//-----------------------------------------------------------------------------

void CellPositionRscl(cellgrid *CG, int cellid, coord *crd, prmtop *tp,
                      double* chi);

void UnpackProcessRsclUpdate(cell *D, coord *crd);

void C2CProcessCnstUpdate(cellgrid *CG, int cellid, int imove, int nmsg);

void ApplyGridCnst(cellgrid *CG, coord *crd, prmtop *tp, trajcon *tj,
                   double Mcut, int app, double* chi);

void VVConstraintVirial(coord *crd, prmtop *tp, trajcon *tj, Energy *sysUV);

#endif
