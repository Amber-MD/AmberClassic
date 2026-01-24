//-----------------------------------------------------------------------------
// ????UpdateExport: when atom positions and velocities are affected by 
//                   constraints, a decision must be made as to whether the
//                   information must be exported.  This routine makes that
//                   decision.  It will update atom locations in the cell and
//                   system and update velocities in the system if the atom is
//                   (or was, before application of the constraint) owned by
//                   the same cell that applied the constraint.  If the atom is
//                   not owned by the process that applied the constraint it
//                   will load the atom for export.  When ???? is "Cnst" the 
//                   function deals with position constraints; when ???? is
//                   Rscl, however, the function deals with updates for volume
//                   rescaling involving constrained groups. 
//
// Arguments:
//   CG:      the cell grid                                             
//   cloc:    the location of the cell in the cell grid                 
//   atmreg:  the region of the cell in which the atom myatm is found   
//   myatm:   the atom whose properties must be exported or updated     
//   xpos:    the updated (current) position of the atom                
//   xposc:   the original position of the atom                         
//   sysloc:  pointer to system coordinates for this atom               
//   sysvel:  pointer to system velocities for this atom (Cnst only)    
//   invdt:   the inverse time step (Cnst only)                         
//-----------------------------------------------------------------------------
#if CNSTTYPE == 0
static void CnstUpdateExport(cellgrid *CG, int* cloc, int atmreg, atomc *myatm,
			     double *xpos, double *xposc, double *sysloc,
			     double *sysvel, double invdt)
#elif CNSTTYPE == 1
static void RsclUpdateExport(cellgrid *CG, int* cloc, int atmreg, atomc *myatm,
			     double *xpos, double *xposc, double *sysloc)
#endif
{
  int i, nexp;
  double dx[3];
  cell *C;

  // The cell of interest
  C = &CG->map[cloc[0]][cloc[1]][cloc[2]];

  // This atom will need to be exported if 
  // it is not in the primary sector,      
  // regardless of whether the same process
  // controls the cell it exports to.      
  dx[0] = xpos[0] - xposc[0];
  dx[1] = xpos[1] - xposc[1];
  dx[2] = xpos[2] - xposc[2];
  if (atmreg != 0) {
    nexp = C->nexp;
    C->Vexport[nexp].id = myatm->id;
    C->Vexport[nexp].dreg = atmreg;
    for (i = 0; i < 3; i++) {
      C->Vexport[nexp].loc[i] = dx[i];
    }
    C->nexp += 1;
  }

  // This atom is controlled by the cell
  else {
    for (i = 0; i < 3; i++) {
      myatm->loc[i] = xpos[i];
      sysloc[i] += dx[i];
#if CNSTTYPE == 0
      sysvel[i] += dx[i]*invdt;
#endif
    }
  }
}

//-----------------------------------------------------------------------------
// RattleGroup: position constraints stage of the RATTLE algorithm for rigid
//              bonds in the context of Velocity-Verlet integration.
//
// This function is adapted from code in the tinker software package and
// the AMBER sff program.                                               
//                                                                      
// Arguments:                                                           
//   CG:       the cell grid                                            
//   crd:      the coordinates (this is where velocity information      
//             resides)                                                
//   tj:       trajectory control information                           
//   tp:       the topology                                             
//   cloc:     location of the cell of interest                         
//   myCC:     vector containing codified information about the atoms   
//             and connectivity in this RATTLE group                    
//   cellidx:  pre-allocated array to hold the indices of all atoms in  
//             this constraint group, indices into the cell's data array
//   xref:     pre-allocated array holding previous coordinates of atoms
//   xpos:     pre-allocated array holding coordinates of atoms         
//   xposc:    pre-allocated array holding copy of xpos to compare after
//             constraints are applied                                  
//   hsctr:    the home sector of the RATTLE group's controlling atom   
//-----------------------------------------------------------------------------
#if CNSTTYPE == 0
static void RattleGroupCnst(cellgrid *CG, coord *crd, trajcon *tj, prmtop *tp,
			    cell *C, int* cloc, int* myCC, int* cellidx,
			    double* xref, double* xpos, double* xposc)
#elif CNSTTYPE == 1
static void RattleGroupRscl(cellgrid *CG, coord *crd, prmtop *tp, cell *C,
			    int* cloc, int* myCC, int* cellidx, double* xref,
			    double* xpos, double* xposc, double* chi,
			    double* invchi)
#endif
{
  int i, j, i3, ididx, atma, nnode;
#if CNSTTYPE == 0
  int done, niter, atmb, atma3, atmb3;
  double delta, l0, rx, ry, rz, rrefx, rrefy, rrefz, dot;
  double r2, rma, rmb, term, xterm, yterm, zterm;
#endif
  atomc *Cdidx;

  // First, this function must find all of
  // the atoms in the constraint group.   
  ididx = 3*myCC[0]+2;
  nnode = myCC[ididx-1];
  for (i = 0; i < nnode; i++) {
    atma = myCC[ididx+i];
    cellidx[i] = C->GPSptr[atma];
    i3 = 3*i;
    Cdidx = &C->data[cellidx[i]];
    for (j = 0; j < 3; j++) {
#if CNSTTYPE == 0
      xref[i3+j] = crd->prvloc[3*atma+j];
#endif
      xpos[i3+j] = Cdidx->loc[j];
      xposc[i3+j] = xpos[i3+j];
    }
  }

#if CNSTTYPE == 0
  // Apply RATTLE constraints to the group
  const double invdt = 1.0/(sqrt(418.4)*tj->dt);
  const double rtol = tj->rattletol;
  niter = 0;
  done = 0;
  const int maxiter = tj->MaxRattleIter;
  while (done == 0 && niter < maxiter) {
    done = 1;
    for (i = 0; i < myCC[0]; i++) {
      atma = myCC[3*i+1];
      atmb = myCC[3*i+2];
      atma3 = 3*atma;
      atmb3 = 3*atmb;
      rx = xpos[atmb3] - xpos[atma3];
      ry = xpos[atmb3+1] - xpos[atma3+1];
      rz = xpos[atmb3+2] - xpos[atma3+2];
      r2 = rx*rx + ry*ry + rz*rz;
      l0 = myCC[3*i+3]/1.0e8;
      l0 *= l0;
      delta = l0 - r2;
      if (fabs(delta) > rtol) {
	done = 0;
	rrefx = xref[atmb3] - xref[atma3];
	rrefy = xref[atmb3+1] - xref[atma3+1];
	rrefz = xref[atmb3+2] - xref[atma3+2];
	dot = rx*rrefx + ry*rrefy + rz*rrefz;
	rma = tp->InvMasses[myCC[ididx + atma]];
	rmb = tp->InvMasses[myCC[ididx + atmb]];
	term = 1.2 * delta / (2.0*dot*(rma+rmb));
	xterm = rrefx*term;
	yterm = rrefy*term;
	zterm = rrefz*term;
	xpos[atma3] -= xterm*rma;
	xpos[atma3+1] -= yterm*rma;
	xpos[atma3+2] -= zterm*rma;
	xpos[atmb3] += xterm*rmb;
	xpos[atmb3+1] += yterm*rmb;
	xpos[atmb3+2] += zterm*rmb;
      }
    }
    niter++;
  }

  // Test to make sure that the maximum  
  // number of iterations was not reached
  if (niter == maxiter && done == 0) {
    printf("RattleGroupCnst >> Error.  Maximum RATTLE iterations exceeded.\n\n"
	   "Group = [\n");
    for (i = 0; i < myCC[0]; i++) {

      // The use of the variables atma and atmb changes here.
      // The program is about to exit in error, so there is  
      // no risk of confusing the variables later on.        
      atma = myCC[ididx + myCC[3*i+1]];
      atmb = myCC[ididx + myCC[3*i+2]];
      printf("    %4d %.4s %.4s -> %4d %.4s %.4s\n", atma,
	     &tp->ResNames[4*LocateResID(tp, atma, 0, tp->nres)],
	     &tp->AtomNames[4*atma], atmb,
	     &tp->ResNames[4*LocateResID(tp, atmb, 0, tp->nres)],
	     &tp->AtomNames[4*atmb]);
    }
    printf("];\n");
    exit(1);
  }
#elif CNSTTYPE == 1

  // When coordinate rescaling is being performed, xref
  // holds the atomic masses of the constraint group   
  double tmass, scofm;
  tmass = 0.0;
  for (i = 0; i < nnode; i++) {
    xref[i] = tp->Masses[myCC[ididx+i]];
    tmass += xref[i];
  }
  tmass = 1.0/tmass;
  for (i = 0; i < nnode; i++) {
    xref[i] *= tmass;
  }
  for (i = 0; i < 3; i++) {

    // For rescaling, the first three indices of xref  
    // are used to hold the shift of the center of mass
    scofm = 0.0;
    for (j = 0; j < nnode; j++) {
      scofm += xpos[3*j+i]*xref[j];
    }
    scofm *= (chi[i] - 1.0);

    // Coordinates are shifted appropriately and
    // then divided by chi so that a subseqeuent
    // multiplication by chi will result in the 
    // correct final position.                  
    for (j = 0; j < nnode; j++) {
      xpos[3*j+i] = (xpos[3*j+i] + scofm)*invchi[i];
    }
  }
#endif

  // Prepare information about the constrained atoms for export
  for (i = 0; i < nnode; i++) {
#if CNSTTYPE == 0
    CnstUpdateExport(CG, cloc, cellidx[i]/C->maxatom, &C->data[cellidx[i]],
		     &xpos[3*i], &xposc[3*i], &crd->loc[3*myCC[ididx+i]],
		     &crd->vel[3*myCC[ididx+i]], invdt);
#elif CNSTTYPE == 1
    RsclUpdateExport(CG, cloc, cellidx[i]/C->maxatom, &C->data[cellidx[i]],
                     &xpos[3*i], &xposc[3*i], &crd->loc[3*myCC[ididx+i]]);
#endif
  }
}

//-----------------------------------------------------------------------------
// CellPosition????: this is a wrapper for SETTLE, RATTLE, and other constraint
//                   functions.  Like the CompBondedIntr function in the Bonded
//                   library, this function loops over all atoms and performs
//                   whatever constraints are required in each case.  More
//                   precisely, each constraint is "owned" by a particular
//                   atom, and whatever cell owns the atom on a particular step
//                   is responsible for computing its constraints.  In the
//                   ???? = "Cnst" form, this will apply position constraints;
//                   in the ???? = "Rscl" form, this will move constrained
//                   groups to prepare for coordinate rescaling.
//
// Arguments:                                                           
//   CG:      the cell grid (this is where position information and     
//            assignments as to which cell will perform the constraint  
//            comes from)                                               
//   cellid:  the location of cell C within the cell grid               
//   crd:     the coordinates                                           
//   tp:      the topology                                              
//   tj:      in the case of cell position constraints, this provides   
//            trajectory control data for things like the time step     
//   chi:     in the case of constrained group center-of-mass rescaling,
//            this provides the scaling factors for all box dimensions  
//-----------------------------------------------------------------------------
#if CNSTTYPE == 0
static void CellPositionCnst(cellgrid *CG, int cellid, coord *crd,
                             prmtop *tp, trajcon *tj)
#elif CNSTTYPE == 1
void CellPositionRscl(cellgrid *CG, int cellid, coord *crd, prmtop *tp,
		      double* chi)
#endif
{
  int h, i, j, atmid, atmid3, nh1, nh1id3, nh1reg;
  int nh2, nh2id3, nh2reg;
  int* cellidx;
  double *ccen, *dbng;
  double* xref;
  double* xpos;
  double* xposc;
  atomc *atmhi, *atmnh1, *atmnh2;
  cnstcomm *tSHL;
  cell *C;

#if CNSTTYPE == 0
  double invdt;
  invdt = 1.0/(sqrt(418.4)*tj->dt);
#elif CNSTTYPE == 1
  double invchi[3];
  for (i = 0; i < 3; i++) {
    invchi[i] = 1.0/chi[i];
  }
#endif

  // Pointers
  C = &CG->data[cellid];
  tSHL = tp->SHL;
  ccen = &CG->celldim[4];
  dbng = CG->dbng;

  // Allocate scratch arrays for temporary locations
  xref = (double*)malloc(tp->RattleGrpMax*sizeof(double));
  xpos = (double*)malloc(tp->RattleGrpMax*sizeof(double));
  xposc = (double*)malloc(tp->RattleGrpMax*sizeof(double));
  cellidx = (int*)malloc(tp->RattleGrpMax/3*sizeof(int));

  // Loop over all atoms in all sectors, apply constraints 
  // that this cell is responsible for, record some results
  // directly and accumulate an export table.              
  C->nexp = 0;
  for (h = 0; h < 8; h++) {
    for (i = 0; i < C->nr[h]; i++) {

      // Test this atom to see if it
      // owns any constraint groups 
      atmhi = &C->map[h][i];
      if (tSHL[atmhi->id].exe == 1 &&
	  IsCentralAtom(C->orig, ccen, atmhi->loc, crd->U.data,
			dbng, crd->isortho) >= 0) {
        atmid = atmhi->id;
        nh1 = C->GPSptr[tSHL[atmid].blist[0]];
        nh1reg = nh1 / C->maxatom;
        atmnh1 = &C->data[nh1];
        nh2 = C->GPSptr[tSHL[atmid].blist[1]];
        nh2reg = nh2 / C->maxatom;
        atmnh2 = &C->data[nh2];
        atmid3 = 3*atmid;
        nh1id3 = 3*atmnh1->id;
        nh2id3 = 3*atmnh2->id;

#if CNSTTYPE == 1
	// When coordinate rescaling is being performed, xref
	// holds the atomic masses of the constraint group   
	xref[0] = tp->Masses[atmid];
	xref[1] = tp->Masses[atmnh1->id];
	xref[2] = tp->Masses[atmnh2->id];
	xref[8] = 1.0/(xref[0] + xref[1] + xref[2]);
	xref[0] *= xref[8];
	xref[1] *= xref[8];
	xref[2] *= xref[8];
#endif
        for (j = 0; j < 3; j++) {
#if CNSTTYPE == 0
          xref[j] = crd->prvloc[atmid3+j];
          xref[3+j] = crd->prvloc[nh1id3+j];
          xref[6+j] = crd->prvloc[nh2id3+j];
#endif
          xpos[j] = atmhi->loc[j];
          xpos[3+j] = atmnh1->loc[j];
          xpos[6+j] = atmnh2->loc[j];
        }
        for (j = 0; j < 9; j++) {
          xposc[j] = xpos[j];
        }

#if CNSTTYPE == 0
        settlec(xref, xpos, &tp->FWtab);

        // The positions and velocities of three atoms
        // should now change.  However, only if the   
        // atoms were in this cell's primary sector   
        // before applying constraints are the system 
        // locations, particle velocities, and other  
        // information authoritative.  Therefore, if  
        // the atoms started in the primary sector (h,
        // nh1reg, or nh2reg is zero) the necessary   
        // information can be updated right away, but 
        // otherwise particular bits of information   
        // such as CHANGE in particle velocity must be
        // transmitted to other cells.                
        CnstUpdateExport(CG, C->gbin, h, atmhi, xpos, xposc, &crd->loc[atmid3],
                         &crd->vel[atmid3], invdt);
        CnstUpdateExport(CG, C->gbin, nh1reg, atmnh1, &xpos[3], &xposc[3],
                         &crd->loc[nh1id3], &crd->vel[nh1id3], invdt);
        CnstUpdateExport(CG, C->gbin, nh2reg, atmnh2, &xpos[6], &xposc[6],
                         &crd->loc[nh2id3], &crd->vel[nh2id3], invdt);
#elif CNSTTYPE == 1
	double scofm;
	for (j = 0; j < 3; j++) {

	  // For rescaling, the first three indices of xref  
	  // are used to hold the shift of the center of mass
	  scofm = (chi[j] - 1.0) *
	    (xpos[j]*xref[0] + xpos[3+j]*xref[1] + xpos[6+j]*xref[2]);

	  // Coordinates are shifted appropriately and
	  // then divided by chi so that a subseqeuent
	  // multiplication by chi will result in the 
	  // correct final position.                  
	  xpos[j] = (xpos[j] + scofm)*invchi[j];
	  xpos[3+j] = (xpos[3+j] + scofm)*invchi[j];
	  xpos[6+j] = (xpos[6+j] + scofm)*invchi[j];
	}
        RsclUpdateExport(CG, C->gbin, h, atmhi, xpos, xposc,
			 &crd->loc[atmid3]);
        RsclUpdateExport(CG, C->gbin, nh1reg, atmnh1, &xpos[3], &xposc[3],
                         &crd->loc[nh1id3]);
        RsclUpdateExport(CG, C->gbin, nh2reg, atmnh2, &xpos[6], &xposc[6],
                         &crd->loc[nh2id3]);
#endif
      }
      else if (tSHL[atmhi->id].exe == 2 &&
	       IsCentralAtom(C->orig, ccen, atmhi->loc,
			     crd->U.data, dbng, crd->isortho) >= 0) {
#if CNSTTYPE == 0
	RattleGroupCnst(CG, crd, tj, tp, C, C->gbin, tSHL[atmhi->id].blist,
			cellidx, xref, xpos, xposc);
#elif CNSTTYPE == 1
	RattleGroupRscl(CG, crd, tp, C, C->gbin, tSHL[atmhi->id].blist,
			cellidx, xref, xpos, xposc, chi, invchi);
#endif
      }
    }
  }

  // Free allocated memory
  free(xref);
  free(xpos);
  free(xposc);
  free(cellidx);
}

//-----------------------------------------------------------------------------
// UnpackProcess???Update: unpacks data packaged (and possibly          
//                         transferred) by C2CProcessCnstUpdate.        
//
// Arguments:
//   D:        the cell to unpack, data which came in from a different process
//   crd:      the coordinates of the system, as known by this process
//   invdt:    the inverse time step
//-----------------------------------------------------------------------------
#if CNSTTYPE == 0
static void UnpackProcessCnstUpdate(cell *D, coord *crd, double invdt)
#elif CNSTTYPE == 1
void UnpackProcessRsclUpdate(cell *D, coord *crd)
#endif
{
  int i, j, atmid, nexp;
  double *ltmp;
#if CNSTTYPE == 0
  double *vtmp;
#endif
  atomc *datm;
  atomb *dRimp;

  // Unpack in cell D
  ltmp = crd->loc;
#if CNSTTYPE == 0
  vtmp = crd->vel;
#endif
  for (i = 1; i < D->import[0].id; i++) {

    // If this atom is to be added to the primary sector,
    // make the official update.  Otherwise, it's just   
    // passing through so update the export array.       
    dRimp = &D->import[i];
    if (dRimp->dreg == 0) {
      atmid = dRimp->id;
      datm = &D->data[D->GPSptr[atmid]];

      // Update cell location, system location,
      // and velocity if required              
      for (j = 0; j < 3; j++) {
        datm->loc[j] += dRimp->loc[j];
        ltmp[3*atmid+j] += dRimp->loc[j];
#if CNSTTYPE == 0
        vtmp[3*atmid+j] += dRimp->loc[j]*invdt;
#endif
      }
    }
    else {

      // Here, nexp refers not to the number of atomb structs 
      // that D is about to export but instead to the number  
      // of atombv structs that D has accumulated and WILL    
      // filter for export on the next C2CProcessTrasfer().   
      // D->nexp was initialized by CellPositionCnst() above, 
      // and holds the number of atoms in constraint groups   
      // whose positions have changed by applying constraints.
      // D->nexp is NOT padded by 1 as it is not exported.    
      nexp = D->nexp;
      D->Vexport[nexp].id = dRimp->id;
      D->Vexport[nexp].dreg = dRimp->dreg;
      D->Vexport[nexp].loc[0] = dRimp->loc[0];
      D->Vexport[nexp].loc[1] = dRimp->loc[1];
      D->Vexport[nexp].loc[2] = dRimp->loc[2];
      D->nexp += 1;
    }
  }
}
