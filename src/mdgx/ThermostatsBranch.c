//-----------------------------------------------------------------------------
// KineticEnergy: compute the kinetic energy of the system from a basis of
//                cells.  The TI variant comptues the kinetic energy of the
//                union of particles between both systems.       
//                                                                      
// Arguments:                                                           
//   CG:     the cell grid                                              
//   crd:    the coordinates (cells reference the ordered list of atom  
//           velocities in this struct)                                 
//   tp:     the topology (for masses)                                  
//   tj:     trajectory control data (flag for Langevin integrator)     
//-----------------------------------------------------------------------------
#if NEED_TI == 0
double KineticEnergy(cellgrid *CG, coord *crd, prmtop *tp, trajcon *tj)
#else
double KineticEnergyTI(cellgrid* CG, coord* crd, prmtop* tp, trajcon *tj)
#endif
{
  int i, j, g3con;
  double result, vx, vy, vz;
  double *mass, *vel;
  cell *C;

  mass = tp->Masses;
  vel = crd->vel;
  result = 0.0;
#ifdef MPI
  for (i = 0; i < CG->MyCellCount; i++) {
    C = &CG->data[CG->MyCellDomain[i]];
#else
  const int ncell = CG->ng[0]*CG->ng[1]*CG->ng[2];
  for (i = 0; i < ncell; i++) {
    C = &CG->data[i];
#endif
    for (j = 0; j < C->nr[0]; j++) {
      if (tp->MobileAtoms[C->data[j].id] == 1) {
	g3con = 3*C->data[j].id;
	vx = vel[g3con];
	vy = vel[g3con+1];
	vz = vel[g3con+2];
	result += mass[C->data[j].id]*(vx*vx + vy*vy + vz*vz);
      }
    }
  }
#if NEED_TI == 1
  prmcorr *prc;
  prc = &tj->prc;
  mass = tp[1].Masses;
  vel = crd[1].vel;
#ifdef MPI
  for (i = 0; i < CG->MyCellCount; i++) {
    C = &CG->data[CG->MyCellDomain[i]];
#else
  for (i = 0; i < ncell; i++) {
    C = &CG[1].data[i];
#endif
    for (j = 0; j < C->nr[0]; j++) {
      if (prc->matchB[C->data[j].id] >= 0 ||
	  tp->MobileAtoms[C->data[j].id] != 1) {
	continue;
      }
      g3con = 3*C->data[j].id;
      vx = vel[g3con];
      vy = vel[g3con+1];
      vz = vel[g3con+2];
      result += mass[C->data[j].id]*(vx*vx + vy*vy + vz*vz);
    }
  }
#endif

#ifdef MPI
  double accresult;
  MPI_Allreduce(&result, &accresult, 1, MPI_DOUBLE, MPI_SUM, CG->dspcomm);
  result = accresult;
#endif

  // Correct the temperature if a Langevin integrator is in use
  if (tj->ntt == 3) {
    result *= 1.0 + 0.5 * tj->lnth.gamma_ln * tj->dt;
  }

  return 0.5 * result;
}

//-----------------------------------------------------------------------------
// SystemTemperature: computes the instantaneous temperature of the system
//                    using the particle velocities.  Assumes the system has no
//                    net velocity.  The TI variant computes the temperature of
//                    the union of particles in both systems.  In this manner,
//                    both systems in a TI calculation have identical
//                    temperatures; the unique atoms in either system are
//                    merely treated as dummy particles, which do not interact
//                    with their surroundings, in the other system.
//                                                                      
// Arguments:                                                           
//   CG:        the cell grid (for velocity information)                
//   tp:        the topology (for masses)                               
//   sysUV:     the system energy (and virial) record                   
//   updateKE:  flag to recompute kinetic energy (if it may be stale)   
//-----------------------------------------------------------------------------
#if NEED_TI == 0
double SystemTemperature(cellgrid *CG, coord *crd, prmtop *tp, Energy *sysUV,
                         trajcon *tj, int updateKE)
#else
double SystemTemperatureTI(cellgrid* CG, coord* crd, prmtop* tp, Energy* sysUV,
			   trajcon *tj, int updateKE)
#endif
{
  if (updateKE == 1) {
#if NEED_TI == 0
    sysUV->kine = KineticEnergy(CG, crd, tp, tj);
#else
    sysUV[0].kine = KineticEnergyTI(CG, crd, tp, tj);
    sysUV[1].kine = sysUV[0].kine;
#endif
  }

#if NEED_TI == 0
  return 2.0*sysUV->kine/(tp->ndf * GASCNST);
#else
  return 2.0*sysUV[0].kine/(tj->prc.nxdf * GASCNST);
#endif
}

//-----------------------------------------------------------------------------
// BerendsenThermostat: the dirt-simple Berendsen thermostat working in the
//                      context of cells.  Beware, this does not sample the
//                      canonical ensemble!  Clamped rescaling values were
//                      introduced to prevent catastrophic rescaling.
//                                                                      
// Reference: Allen and Tildesley, Eq. 7.59, pg 231                     
//                                                                      
// Arguments:                                                           
//   CG:      the cell grid                                             
//   crd:     the coordinates (cells reference coordinates for velocity 
//            information)                                              
//   tj:      trajectory control information                            
//   sysUV:   energy and virial information                             
//-----------------------------------------------------------------------------
#if NEED_TI == 0
void BerendsenThermostat(cellgrid *CG, coord *crd, trajcon *tj, Energy *sysUV)
#else
void BerendsenThermostatTI(cellgrid* CG, coord* crd, trajcon *tj,
			   Energy* sysUV)
#endif
{
  int i, j, g3con;
  double tscale, chi, T;
  double *vtmp;
  cell *C;

  // Temperature was computed after the force computation;
  // velocities have been updated since the temperature   
  // computation but that is part of the algorithm.       
  T = sysUV->T;
  if (T < 0.5 * tj->Ttarget) {
    tscale = 1.0;
  }

  // Too hot!  Let's call it twice the target temperature
  else if (T > 2.0 * tj->Ttarget) {
    tscale = -0.5;
  }

  // Ahh, just right...
  else {
    tscale = (tj->Ttarget / T) - 1.0;
  }
  chi = sqrt(1.0 + tscale * tj->dt / tj->BerendsenTCoupl);
  vtmp = crd->vel;
  for (i = 0; i < CG->MyCellCount; i++) {

    // Velocity scaling is not restricted to atoms with mass,
    // as atoms without mass (virtual sites) should have zero
    // velocity already.                                     
    C = &CG->data[CG->MyCellDomain[i]];
    for (j = 0; j < C->nr[0]; j++) {
      g3con = 3*C->data[j].id;
      vtmp[g3con] *= chi;
      vtmp[g3con+1] *= chi;
      vtmp[g3con+2] *= chi;
    }
  }
#if NEED_TI == 1
  vtmp = crd[1].vel;
  for (i = 0; i < CG[1].MyCellCount; i++) {

    // Velocity scaling is not restricted to atoms with mass,
    // as atoms without mass (virtual sites) should have zero
    // velocity already.                                     
    C = &CG[1].data[CG->MyCellDomain[i]];
    for (j = 0; j < C->nr[0]; j++) {
      g3con = 3*C->data[j].id;
      vtmp[g3con] *= chi;
      vtmp[g3con+1] *= chi;
      vtmp[g3con+2] *= chi;
    }
  }
#endif
}
