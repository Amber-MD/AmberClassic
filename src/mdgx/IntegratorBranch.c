//-----------------------------------------------------------------------------
// RemoveMomentum: this function removes the net momentum of the system given
//                 all atoms within a cell framework.  The TI variant applies
//                 to the special case of thermodynamic integration involving
//                 two coupled systems: the net momentum of the union of both
//                 systems is computed and removed from both systems.
//
// Arguments:                                                           
//   CG:    the cell grid (directs which atoms accumulate momentum)     
//   crd:   the coordinates (for velocities)                            
//   tp:    the topology array (for masses)                             
//-----------------------------------------------------------------------------
#if NEED_TI == 0
static void RemoveMomentum(cellgrid *CG, coord *crd, prmtop *tp)
#else
static void RemoveMomentumTI(cellgrid* CG, coord* crd, prmtop* tp,
                             trajcon *tj)
#endif
{
#if NEED_TI == 1
  int h;
#endif
  int i, j, g3con;
  double mvx, mvy, mvz, invn, jmass;
  double *vtmp, *mtmp;
  cell *C;
#if NEED_TI == 1
  prmcorr *prc;
#endif

  // In the case of TI, crd, tp, and CG are actually array 
  // variables, but we can treat them as pointers to single
  // variables for now as we only need the first element of
  // the array in the case of TI.                          
  vtmp = crd->vel;
  mtmp = tp->Masses;
  mvx = 0.0;
  mvy = 0.0;
  mvz = 0.0;
  for (i = 0; i < CG->MyCellCount; i++) {
    C = &CG->data[CG->MyCellDomain[i]];
    for (j = 0; j < C->nr[0]; j++) {
      jmass = mtmp[C->data[j].id];
      g3con = 3*C->data[j].id;
      mvx += vtmp[g3con]*jmass;
      mvy += vtmp[g3con+1]*jmass;
      mvz += vtmp[g3con+2]*jmass;
    }
  }

#if NEED_TI == 0
  invn = 1.0/tp[0].TotalMass;
#else
  // If either system has unique atoms, the
  // total mass must be computed with care.
  prc = &tj->prc;
  if (prc->relate > 0) {
    invn = 0.0;
    vtmp = crd[1].vel;
    mtmp = tp[1].Masses;
    for (i = 0; i < CG[1].MyCellCount; i++) {
      C = &CG[1].data[CG[1].MyCellDomain[i]];
      for (j = 0; j < C->nr[0]; j++) {
        if (prc->matchB[C->data[j].id] < 0) {
          jmass = mtmp[C->data[j].id];
          g3con = 3*C->data[j].id;
          mvx += vtmp[g3con]*jmass;
          mvy += vtmp[g3con+1]*jmass;
          mvz += vtmp[g3con+2]*jmass;
          invn += jmass;
        }
      }
    }
    invn += tp[0].TotalMass;
    invn = 1.0/invn;
  }
  else {
    invn = 1.0/tp[0].TotalMass;
  }
#endif
#ifdef MPI
  double mxport[3], mbuff[3];
  mxport[0] = mvx;
  mxport[1] = mvy;
  mxport[2] = mvz;
  MPI_Allreduce(mxport, mbuff, 4, MPI_DOUBLE, MPI_SUM, CG->dspcomm);
  mvx = mbuff[0];
  mvy = mbuff[1];
  mvz = mbuff[2];
#endif

  // Compute velocity of the center of mass
  mvx *= invn;
  mvy *= invn;
  mvz *= invn;
#if NEED_TI == 1
  for (h = 0; h < 2; h++) {
    vtmp = crd[h].vel;
#endif
#if NEED_TI == 0
    for (i = 0; i < CG->MyCellCount; i++) {
      C = &CG->data[CG->MyCellDomain[i]];
#else
    for (i = 0; i < CG[h].MyCellCount; i++) {
      C = &CG[h].data[CG[h].MyCellDomain[i]];
#endif
      for (j = 0; j < C->nr[0]; j++) {
        g3con = 3*C->data[j].id;
        vtmp[g3con] -= mvx;
        vtmp[g3con+1] -= mvy;
        vtmp[g3con+2] -= mvz;
      }
    }
#if NEED_TI == 1
  }
#endif
}

