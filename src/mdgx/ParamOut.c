#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "mdgxVector.h"
#include "Matrix.h"
#include "ParamFit.h"
#include "Topology.h"
#include "CrdManip.h"
#include "Parse.h"
#include "ChargeFit.h"
#include "VirtualSites.h"
#include "Manual.h"
#include "Nonbonded.h"
#include "Trajectory.h"
#include "Macros.h"
#include "ParamRead.h"
#include "ParamFit.h"

#include "CrdManipDS.h"

//-----------------------------------------------------------------------------
// SetParmOutputFlags: make decisions regarding whether to print the parameters
//                     to output.  By default, all parameters which were
//                     present in the original parm.dat and frcmod files are
//                     printed to output, but if there are cloned parameters
//                     then they will be printed only if they are found in one
//                     of the systems under consideration.  By turning up the
//                     report level, all cloned parameters can appear in the   
//                     output, and by turning down the report level only
//                     parameters needed for the systems involved in the
//                     fitting exercise will make it into the fit.      
//                                                                      
// Arguments:                                                           
//   mp:       fitting control data (contains atom types and examples of
//             each atom in various topologies)                         
//-----------------------------------------------------------------------------
static void SetParmOutputFlags(prmset *mp)
{
  int i, j, k, atmid, replvl;
  prmtop *tp;

  replvl = mp->reportall;

  // Atoms
  for (i = 0; i < mp->natom; i++) {
    if (replvl == 2 || (replvl == 1 && mp->atoms[i].dup == 0)) {
      mp->atoms[i].inreport = 1;
    }
    else {
      mp->atoms[i].inreport = 0;
    }
  }

  // Determine whether these atoms are present in ANY system used in this
  // fit.  If so, they get a pass under the standard level of reporting.
  if (replvl == 1) {
    for (i = 0; i < mp->nunisys; i++) {
      for (j = 0; j < mp->nconf; j++) {
        if (mp->conf[j].GroupNum == i) {
          tp = mp->conf[j].tp;
          for (k = 0; k < tp->natom; k++) {
            atmid = CrossRefAtomType(mp, &tp->AtomTypes[4*k]);
            mp->atoms[atmid].inreport = 1;
          }
          break;
        }
      }
    }
  }

  // Bonds, angles, dihedrals, and CMAPs
  for (i = 0; i < mp->nbond; i++) {
    if (replvl == 2 || (replvl == 1 && mp->bonds[i].dup == 0)) {
      mp->bonds[i].inreport = 1;
    }
    else {
      mp->bonds[i].inreport = 0;
    }
  }
  for (i = 0; i < mp->nangl; i++) {
    if (replvl == 2 || (replvl == 1 && mp->angls[i].dup == 0)) {
      mp->angls[i].inreport = 1;
    }
    else {
      mp->angls[i].inreport = 0;
    }
  }
  for (i = 0; i < mp->ntor; i++) {
    if (replvl == 2 || (replvl == 1 && mp->torsions[i].dup == 0)) {
      mp->torsions[i].inreport = 1;
    }
    else {
      mp->torsions[i].inreport = 0;
    }
  }
  for (i = 0; i < mp->nhb1012; i++) {
    if (replvl == 2 || (replvl == 1 && mp->hb1012[i].dup == 0)) {
      mp->hb1012[i].inreport = 1;
    }
    else {
      mp->hb1012[i].inreport = 0;
    }
  }
  for (i = 0; i < mp->ncmaps; i++) {
    if (replvl == 0) {
      mp->cmaps[i].inreport = 1;
    }
    else {
      mp->cmaps[i].inreport = 0;
    }
  }
  
  // Determine whether these bonds, angles, or dihedrals are
  // present in ANY system used in this fit.  If so, they get
  // a pass and are included in standard parameter output.  This
  // is what really allows parameters with branched atom types
  // to be printed.
  if (replvl == 1) {
    for (i = 0; i < mp->nunisys; i++) {
      for (j = 0; j < mp->nconf; j++) {
        if (mp->conf[j].GroupNum == i) {
          for (k = 0; k < mp->conf[j].bmap.nbond; k++) {
            mp->bonds[mp->conf[j].bmap.id[k].key].inreport = 1;
          }
          for (k = 0; k < mp->conf[j].amap.nangl; k++) {
            mp->angls[mp->conf[j].amap.id[k].key].inreport = 1;
          }
          for (k = 0; k < mp->conf[j].hmap.ntterm; k++) {
            mp->torsions[mp->conf[j].hmap.id[k].key].inreport = 1;
          }
          break;
        }
      }
    }
  }

  // For printing frcmod files, find the fitted parameters.  CMAPs will
  // always get printed, as they can only appear in frcmod type files.
  if (replvl == 0) {
    for (i = 0; i < mp->nbond; i++) {
      if (mp->bonds[i].fitcolX >= 0) {
        atmid = CrossRefAtomType(mp, mp->bonds[i].atype);
        if (atmid >= 0) {
          mp->atoms[atmid].inreport = mp->atoms[atmid].dup;
        }
        atmid = CrossRefAtomType(mp, mp->bonds[i].btype);
        if (atmid >= 0) {
          mp->atoms[atmid].inreport = mp->atoms[atmid].dup;
        }
        mp->bonds[i].inreport = 1;
      }
    }
    for (i = 0; i < mp->nangl; i++) {
      if (mp->angls[i].fitcolX >= 0) {
        atmid = CrossRefAtomType(mp, mp->angls[i].atype);
        if (atmid >= 0) {
          mp->atoms[atmid].inreport = mp->atoms[atmid].dup;
        }
        atmid = CrossRefAtomType(mp, mp->angls[i].btype);
        if (atmid >= 0) {
          mp->atoms[atmid].inreport = mp->atoms[atmid].dup;
        }
        atmid = CrossRefAtomType(mp, mp->angls[i].ctype);
        if (atmid >= 0) {
          mp->atoms[atmid].inreport = mp->atoms[atmid].dup;
        }
        mp->angls[i].inreport = 1;
      }
    }
    for (i = 0; i < mp->ntor; i++) {
      if (mp->torsions[i].fitcol >= 0) {
        atmid = CrossRefAtomType(mp, mp->torsions[i].atype);
        if (atmid >= 0) {
          mp->atoms[atmid].inreport = mp->atoms[atmid].dup;
        }
        atmid = CrossRefAtomType(mp, mp->torsions[i].btype);
        if (atmid >= 0) {
          mp->atoms[atmid].inreport = mp->atoms[atmid].dup;
        }
        atmid = CrossRefAtomType(mp, mp->torsions[i].ctype);
        if (atmid >= 0) {
          mp->atoms[atmid].inreport = mp->atoms[atmid].dup;
        }
        atmid = CrossRefAtomType(mp, mp->torsions[i].dtype);
        if (atmid >= 0) {
          mp->atoms[atmid].inreport = mp->atoms[atmid].dup;
        }
        mp->torsions[i].inreport = 1;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// PrintParmAtoms: print the masses involved in a frcmod file for this  
//                 parameter fit.                                     
//                                                                      
// Arguments:                                                           
//   mp:       fitting control data (contains atom types and examples of
//             each atom in various topologies)                         
//   outp:     the output frcmod file (pointer)                         
//-----------------------------------------------------------------------------
static void PrintParmAtoms(prmset *mp, FILE *outp)
{
  int i;

  for (i = 0; i < mp->natom; i++) {
    if (mp->atoms[i].inreport == 1) {
      fprintf(outp, "%.2s               %10.4lf  %10.4lf  %s\n",
              mp->atoms[i].atype, mp->atoms[i].mass, mp->atoms[i].apol,
              mp->atoms[i].comment);
    }
  }
  fprintf(outp, "\n");
}

//-----------------------------------------------------------------------------
// PrintParmBond: this procedure is much like PrintParmAtoms above.     
//                                                                      
// Arguments:                                                           
//   mp:       fitting control data (contains atom types and examples of
//             each atom in various topologies)                         
//   order:    the order of bonded terms to seek out (2 = bonds,        
//             3 = angles, 4 = dihedrals)                               
//   outp:     the output file                                          
//   x:        the vector of fitted parameters                          
//-----------------------------------------------------------------------------
static void PrintParmBond(prmset *mp, int order, double* x, FILE *outp)
{
  int i, ilim, slen;
  double stiffness, equilibrium;
  char tmpline[1024];

  // Print bonds, angles, and dihedrals to the new parameter file
  if (order == 2) ilim = mp->nbond;
  else if (order == 3) ilim = mp->nangl;
  else if (order == 4) ilim = mp->ntor;
  else if (order == 5) ilim = mp->ntor;
  else if (order == 9) ilim = mp->nhb1012;
  for (i = 0; i < ilim; i++) {
    if (order == 2 && mp->bonds[i].inreport == 1) {
      if (mp->bonds[i].fitcolX >= 0) {
        stiffness = x[mp->bonds[i].fitcolX];
        if (mp->FitBondEq == 1) {
          stiffness += x[mp->bonds[i].fitcolY];
          equilibrium = (x[mp->bonds[i].fitcolX] * mp->bonds[i].lbasisX + 
                         x[mp->bonds[i].fitcolY] * mp->bonds[i].lbasisY) /
            stiffness;
        }
        else {
          equilibrium = mp->bonds[i].l0;
        }
        if (strncmp(mp->bonds[i].comment, "Branched from ", 14) == 0) {
          slen = strlen(mp->bonds[i].comment);
	  if (slen > 0) {
            mp->bonds[i].comment[slen] = ',';
            mp->bonds[i].comment[slen+1] = ' ';
            mp->bonds[i].comment[slen+2] = '\0';
            slen += 2;
	  }
        }
        else {
          slen = 0;
        }
        strcpy(mp->bonds[i].comment + slen, mp->icomm);
      }
      else {
        stiffness = mp->bonds[i].K;
        equilibrium = mp->bonds[i].l0;
      }
      sprintf(tmpline, "%.2s-%.2s            %10.4lf  %10.4lf  %s\n",
              mp->bonds[i].atype, mp->bonds[i].btype, stiffness,
              equilibrium, mp->bonds[i].comment);
      fprintf(outp, "%s", tmpline);
    }
    if (order == 3 && mp->angls[i].inreport == 1) {
      if (mp->angls[i].fitcolX >= 0) {
        stiffness = x[mp->angls[i].fitcolX];
        if (mp->FitAnglEq == 1) {
          stiffness += x[mp->angls[i].fitcolY];
          equilibrium = (x[mp->angls[i].fitcolX] * mp->angls[i].tbasisX + 
                         x[mp->angls[i].fitcolY] * mp->angls[i].tbasisY) /
            stiffness;
        }
        else {
          equilibrium = mp->angls[i].th0;
        }
        if (strncmp(mp->angls[i].comment, "Branched from ", 14) == 0) {
          slen = strlen(mp->angls[i].comment);
	  if (slen > 0) {
            mp->angls[i].comment[slen] = ',';
            mp->angls[i].comment[slen+1] = ' ';
            mp->angls[i].comment[slen+2] = '\0';
            slen += 2;
	  }
        }
        else {
          slen = 0;
        }
        strcpy(mp->angls[i].comment + slen, mp->icomm);
      }
      else {
        stiffness = mp->angls[i].K;
        equilibrium = mp->angls[i].th0;
      }
      sprintf(tmpline, "%.2s-%.2s-%.2s         %10.4lf  %10.2lf  %s\n",
              mp->angls[i].atype, mp->angls[i].btype, mp->angls[i].ctype,
              stiffness, (180.0/PI)*equilibrium, mp->angls[i].comment);
      fprintf(outp, "%s", tmpline);
    }
    if (order == 4 && mp->torsions[i].impr == 0 &&
        mp->torsions[i].inreport == 1) {
      if (mp->torsions[i].fitcol >= 0) {
        stiffness = x[mp->torsions[i].fitcol];
        if (strendswith(mp->torsions[i].comment, mp->icomm) == 0) {
          slen = strlen(mp->torsions[i].comment);
          if (strendswith(mp->torsions[i].comment, " ") == 0 && slen > 0) {
            if (strendswith(mp->torsions[i].comment, ",") == 0) {
              sprintf(&mp->torsions[i].comment[slen], ",");
              slen += 1;
            }
            sprintf(&mp->torsions[i].comment[slen], " ");
            slen += 1;
          }
          sprintf(&mp->torsions[i].comment[slen], "%s", mp->icomm);
        }
      }
      else {
        stiffness = mp->torsions[i].K;
      }
      sprintf(tmpline, "%.2s-%.2s-%.2s-%.2s   1  %10.5lf %6.1lf %4.1lf  "
              "%s\n", mp->torsions[i].atype, mp->torsions[i].btype,
              mp->torsions[i].ctype, mp->torsions[i].dtype, stiffness,
              mp->torsions[i].phase*180.0/PI,
              mp->torsions[i].singlet * mp->torsions[i].pn,
              mp->torsions[i].comment);
      fprintf(outp, "%s", tmpline);
    }
    if (order == 5 && mp->torsions[i].impr == 1 &&
        mp->torsions[i].inreport == 1) {
      if (mp->torsions[i].fitcol >= 0) {
        stiffness = x[mp->torsions[i].fitcol];
        strcpy(mp->torsions[i].comment, mp->icomm);
      }
      else {
        stiffness = mp->torsions[i].K;
      }
      sprintf(tmpline, "%.2s-%.2s-%.2s-%.2s      %10.5lf %6.1lf %4.1lf  "
              "%s\n", mp->torsions[i].atype, mp->torsions[i].btype,
              mp->torsions[i].ctype, mp->torsions[i].dtype, stiffness,
              mp->torsions[i].phase*180.0/PI,
              mp->torsions[i].singlet * mp->torsions[i].pn,
              mp->torsions[i].comment);
      fprintf(outp, "%s", tmpline);
    }
    if (order == 9 && mp->hb1012[i].inreport == 1) {

      // Currently there is no support for
      // fitting H-bond 10-12 potentials. 
      sprintf(tmpline, "  %.2s  %.2s         %10.4lf  %10.4lf  %s\n",
              mp->hb1012[i].atype, mp->hb1012[i].btype, mp->hb1012[i].Aterm,
              mp->hb1012[i].Bterm, mp->hb1012[i].comment);
      fprintf(outp, "%s", tmpline);
    }
  }
  fprintf(outp, "\n");
}

//-----------------------------------------------------------------------------
// PrintParmCMAP: print a fitted CMAP to a frcmod file.
//
// Arguments:
//   mp:       fitting control data (contains atom types and examples of
//             each atom in various topologies)
//   x:        the vector of fitted parameters
//   outp:     the output file
//-----------------------------------------------------------------------------
static void PrintParmCMAP(prmset *mp, double* x, FILE *outp)
{
  int i, j, k, m, npr, nbin, nval;
  double val1, val2, evalue;
  
  npr = 0;
  for (i = 0; i < mp->ncmaps; i++) {

    // Skip if this is not flagged for inclusion in the parameter output
    if (mp->cmaps[i].inreport == 0) {
      continue;
    }
    npr++;
    fprintf(outp, "%%FLAG CMAP_COUNT %3d\n", npr);
    fprintf(outp, "%%FLAG CMAP_TITLE\n");
    fprintf(outp, "%s\n", mp->cmaps[i].title);
    fprintf(outp, "%%FLAG CMAP_RESLIST %2d\n", mp->cmaps[i].resnames.row);
    for (j = 0; j < mp->cmaps[i].resnames.row; j++) {
      fprintf(outp, "%4.4s", mp->cmaps[i].resnames.map[j]);
    }
    fprintf(outp, "\n%%FLAG CMAP_RESOLUTION %2d\n", mp->cmaps[i].esrf.row);
    fprintf(outp, "%%FLAG CMAP_PARAMETER\n");

    // Extract the fitted CMAP from the solution vector, or if it was
    // not fitted then just copy the input surface.
    if (mp->cmaps[i].fitcoli >= 0) {
      nbin = mp->cmaps[i].esrfFit.row;
      m = 0;
      for (j = 0; j < nbin; j++) {
        for (k = 0; k < nbin; k++) {
          mp->cmaps[i].esrfFit.map[j][k] = x[mp->cmaps[i].fitcoli + m];
          m++;
        }
      }
    }
    else {
      DestroyDmat(&mp->cmaps[i].esrfFit);
      CopyDmat(&mp->cmaps[i].esrfFit, &mp->cmaps[i].esrf, 0);
    }

    // If the CMAP's fitted resolution does not match its original (and thus,
    // printed resolution), interpolate the result from the fitted surface.
    nbin = mp->cmaps[i].esrf.row;
    nval = 0;
    if (mp->cmaps[i].esrfFit.row != mp->cmaps[i].esrf.row) {
      ComputeLinearCmapDerivatives(&mp->cmaps[i].esrfFit,
                                   &mp->cmaps[i].dphiFit,
                                   &mp->cmaps[i].dpsiFit,
                                   &mp->cmaps[i].d2ppFit);
      for (j = 0; j < nbin; j++) {
        val1 = (double)j / (double)nbin * (2.0 * PI);
        for (k = 0; k < nbin; k++) {
          val2 = (double)k / (double)nbin * (2.0 * PI);
          evalue = EvaluateCMAPEnergy(val1, val2, &mp->cmaps[i].esrfFit,
                                      &mp->cmaps[i].dphiFit,
                                      &mp->cmaps[i].dpsiFit,
                                      &mp->cmaps[i].d2ppFit);
          fprintf(outp, "%9.5lf", evalue);
          nval++;
          if (nval == 8) {
            fprintf(outp, "\n");
	    nval = 0;
          }
        }
      }
    }
    else {
      for (j = 0; j < nbin; j++) {
        for (k = 0; k < nbin; k++) {
          fprintf(outp, "%9.5lf", mp->cmaps[i].esrfFit.map[j][k]);
          nval++;
          if (nval == 8) {
            fprintf(outp, "\n");
	    nval = 0;
          }
        }
      }
    }
    if (nval > 0) {
      fprintf(outp, "\n");
    }
  }
}


//-----------------------------------------------------------------------------
// MapSamplingIndices: while the fitcol indices in each bond, angle, or 
//                     torsion map each adjustable term into the fitting
//                     matrix, the samprow indices are needed to map    
//                     each term into the distinct sampling matrices.   
//                                                                      
// Arguments:                                                           
//   mp:       fitting control data (contains atom types and examples of
//             each atom in various topologies)                         
//   order:    the order of bonded terms to seek out (2 = bonds,        
//             3 = angles, 4 = dihedrals)                               
//-----------------------------------------------------------------------------
static void MapSamplingIndices(prmset *mp, int order)
{
  int i, j, ilim;

  // Map the indices of the sampling matrices
  j = 0;
  ilim = (order == 2) ? mp->nbond : (order == 3) ? mp->nangl : mp->ntor;
  for (i = 0; i < ilim; i++) {
    if (order == 2) {
      if (mp->bonds[i].fitcolX >= 0) {
        mp->bonds[i].samprow = j;
        j++;
      }
      else {
        mp->bonds[i].samprow = -1;
      }
    }
    else if (order == 3) {
      if (mp->angls[i].fitcolX >= 0) {
        mp->angls[i].samprow = j;
        j++;
      }
      else {
        mp->angls[i].samprow = -1;
      }
    }
    else if (order == 4) {
      if (mp->torsions[i].fitcol >= 0) {
        mp->torsions[i].samprow = j;
        j++;
      }
      else {
        mp->torsions[i].samprow = -1;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// AccumulateSamplingTable: accumulate the sampling histogram based on all
//                          conformations.  This routine is encapsulated to
//                          allow it to be used in multiple cases.
//                                                                      
// Arguments:                                                           
//   mp:       fitting control data (contains atom types and examples of
//             each atom in various topologies)                         
//   S:        the short-form sampling histogram (set to 18 bins for    
//             the purpose of output file reporting)                    
//   Sv:       the long-form sampling (set to 720 bins for the purpose  
//             of minimum checking)                                     
//   order:    the order of bonded terms to seek out (2 = bonds,        
//             3 = angles, 4 = dihedrals)                               
//   sysid:    system ID number, -1 for accumulation over all systems   
//-----------------------------------------------------------------------------
static void AccumulateSamplingTable(prmset *mp, imat *S, dmat* Sv, int order,
                                    int sysid)
{
  int i, j, k, cfmin, cfmax, dev, nrec, nvar, currsys;
  int mbinL, mbinH;
  int* isvar;
  int* rowinS;
  double gssfac, gssarg, bondlen, anglwidth, dihetwist;
  double* dtmp;
  bidx *tbond;
  aidx *tangl;
  hidx *tdihe;

  // Compile sampling histograms.  A relevant sysid (zero or greater
  // value) causes sampling for that specific system to get tested. 
  // Sending this function a sysid of -1 causes ALL conformations to
  // contribute to the sampling table.                              
  if (sysid >= 0) {
    cfmin = mp->FirstConf[sysid];
    cfmax = mp->LastConf[sysid];
  }
  else {
    cfmin = 0;
    cfmax = mp->nconf-1;
  }

  // Allocate wisdom tables and pre-compute constants
  isvar = (int*)malloc(sizeof(int));
  rowinS = (int*)malloc(sizeof(int));
  gssfac = 1.0/(sqrt(2.0*PI));
  currsys = -1;

  // Loop over all (relevant) conformations
  nrec = 0;
  for (i = cfmin; i <= cfmax; i++) {
    if (sysid >= 0 && mp->conf[i].GroupNum != sysid) {
      continue;
    }
    if (order == 2) {
      dtmp = mp->conf[i].bmap.val;
      for (j = 0; j < mp->conf[i].bmap.nbond; j++) {
        tbond = &mp->conf[i].bmap.id[j];
        if (mp->bonds[tbond->key].fitcolX < 0) {
          continue;
        }

        // Log this sample in the coarse-grained histogram
        dev = (dtmp[j] / 0.027777) + 9.00001;
        if (dev >= 0 && dev < 18) {
          S->map[mp->bonds[tbond->key].samprow][dev] += 1;
        }

        // Log this sample in the fine-grained
        // histogram for minima checking      
        if (sysid == -1) {
          bondlen = (dtmp[j] + mp->bonds[tbond->key].l0) * 100.0;
          mbinL = bondlen - 6;
          mbinH = mbinL + 12;
          if (mbinL < 0 || mbinL >= 720 || mbinH < 0 || mbinH >= 720) {
            continue;
          }
          for (k = mbinL; k <= mbinH; k++) {          
            Sv->map[mp->bonds[tbond->key].samprow][k] +=
              gssfac*exp(-0.5*(bondlen-k)*(bondlen-k));
          }
        }
      }
    }
    else if (order == 3) {
      dtmp = mp->conf[i].amap.val;
      for (j = 0; j < mp->conf[i].amap.nangl; j++) {
        tangl = &mp->conf[i].amap.id[j];
        if (mp->angls[tangl->key].fitcolX < 0) {
          continue;
        }

        // Log this sample in the coarse-grained histogram
        dev = (dtmp[j] / 0.023271) + 9.00001;
        if (dev >= 0 && dev < 18) {
          S->map[mp->angls[tangl->key].samprow][dev] += 1;
        }

        // Log this sample in the fine-grained
        // histogram for minima checking      
        if (sysid == -1) {
          anglwidth = (dtmp[j] + mp->angls[tangl->key].th0)*180.0/PI + 360.0;
          mbinL = anglwidth - 6;
          mbinH = mbinL + 12;
          for (k = mbinL; k <= mbinH; k++) {          
            Sv->map[mp->angls[tangl->key].samprow][k] +=
              gssfac*exp(-0.5*(anglwidth-k)*(anglwidth-k));
          }
        }
      }
    }
    else {

      // Set pointers and prepare wisdom tables
      dtmp = mp->conf[i].hmap.val;
      if (mp->conf[i].GroupNum != currsys) {
        free(isvar);
        free(rowinS);
        nvar = mp->conf[i].hmap.ntterm;
        isvar = (int*)malloc(nvar*sizeof(int));
        rowinS = (int*)malloc(nvar*sizeof(int));
        currsys = mp->conf[i].GroupNum;
        for (j = 0; j < nvar; j++) {
          tdihe = &mp->conf[i].hmap.id[j];
          if (mp->torsions[tdihe->key].fitcol < 0) {
            isvar[j] = 0;
            rowinS[j] = -1;
          }
          else {
            isvar[j] = 1;
            rowinS[j] = mp->torsions[tdihe->key].samprow;
          }
        }
      }
      for (j = 0; j < nvar; j++) {
        if (isvar[j] == 0) {
          continue;
        }

        // Log this sample in the coarse-grained histogram
        dev = dtmp[j] * 9.0 / PI + 9.00001;
        if (dev >= 0 && dev < 18) {
          S->map[rowinS[j]][dev] += 1;
        }

        // Log this sample in the fine-grained
        // histogram for minima checking      
        if (sysid == -1) {
          dihetwist = dtmp[j] * 180.0 / PI + 180.0;
          mbinL = dihetwist - 6;
          mbinH = mbinL + 12;
          for (k = mbinL; k <= mbinH; k++) {

            // The second index may need to be imaged,
            // particularly in the case of dihedrals. 
            if (k >= 720) {
              k -= 720;
            }
            else if (k < 0) {
              k += 720;
            }
            Sv->map[rowinS[j]][k] +=
              gssfac*exp(-0.5*(dihetwist-k)*(dihetwist-k));
          }
        }
      }
    }
  }

  // Free allocated memory
  free(isvar);
  free(rowinS);
}

//-----------------------------------------------------------------------------
// CharacterHistogram: print a character-keyed histogram, values ranging
//                     from 0 to 10.                                    
//                                                                      
// Arguments:                                                           
//   S:      integer matrix containing the histogram data               
//   ncol:   the number of bins in the histogram                        
//   outp:   the output file                                            
//-----------------------------------------------------------------------------
static void CharacterHistogram(int* S, int ncol, FILE *outp)
{
  int i;
  char occmap[16];

  sprintf(occmap, " .:=eoUO0@");
  for (i = 0; i < ncol; i++) {
    if (S[i] < 10) {
      fprintf(outp, "%c", occmap[S[i]]);
    }
    else {
      fprintf(outp, "X");
    }
  }
}

//-----------------------------------------------------------------------------
// MatchSamprow2Column: search through all adjustable variables of a specified
//                      order and find the variable's index into the fitting   
//                      matrix columns based on its index into the sampling
//                      matrix.                 
//                                                                      
// Arguments:                                                           
//   mp:      fitting control data (contains notes about the fitting matrix)
//   sidx:    the index into the sampling table                         
//   order:   2 = bonds, 3 = angles, 4 = torsions                       
//   XorY:    Flag to indicate whether to return the X or Y column      
//            when fitting bonds and angles                             
//-----------------------------------------------------------------------------
static int MatchSamprow2Column(prmset *mp, int sidx, int order, int XorY)
{
  int i;

  if (order == 2) {
    for (i = 0; i < mp->nbond; i++) {
      if (mp->bonds[i].samprow == sidx && mp->bonds[i].fitcolX != -1) {
        if (XorY == 0) {
          return mp->bonds[i].fitcolX;
        }
        else {
          return mp->bonds[i].fitcolY;
        }
      }
    }
    printf("MatchSamprow2Column >> Unable to find match for bond "
           "%d, %.4s-%.4s.\n", i, mp->bonds[i].atype, mp->bonds[i].btype);
    exit(1);
  }
  else if (order == 3) {
    for (i = 0; i < mp->nangl; i++) {
      if (mp->angls[i].samprow == sidx && mp->angls[i].fitcolX != -1) {
        if (XorY == 0) {
          return mp->angls[i].fitcolX;
        }
        else {
          return mp->angls[i].fitcolY;
        }
      }
    }
    printf("MatchSamprow2Column >> Unable to find match for angle "
           "%d, %.4s-%.4s-%.4s.\n", i, mp->angls[i].atype, mp->angls[i].btype,
           mp->angls[i].ctype);
    exit(1);
  }
  else if (order == 4) {
    for (i = 0; i < mp->ntor; i++) {
      if (mp->torsions[i].samprow == sidx && mp->torsions[i].fitcol != -1) {
        return mp->torsions[i].fitcol;
      }
    }
    printf("MatchSamprow2Column >> Unable to find match for angle "
           "%d, %.4s-%.4s-%.4s-%.4s.\n", i, mp->torsions[i].atype,
           mp->torsions[i].btype, mp->torsions[i].ctype,
           mp->torsions[i].dtype);
    exit(1);
  }

  return -1;
}

//-----------------------------------------------------------------------------
// FindColumnTerm: find the term, be it a bond, angle, or dihedral, 
//                 corresponding to the column of interest.  It is terribly
//                 inefficient to search through all bonds, angles, and
//                 dihedrals just to find the one that feeds into a particular
//                 matrix column when all of that data could be recorded, but
//                 this routine is only called at output and therefore saves
//                 some complexity in the data structures.
//                                                                      
// Arguments:                                                           
//   mp:      fitting control data (contains notes about the fitting    
//            matrix)                                                   
//   cc:      the column to match a term against                        
//   outp:    the output file                                           
//   order:   returned value, identifying the column as associated with 
//            a bond (return 2), angle (return 3), or torsion (return 4)
//   parmidx: the parameter index (the relevant list is implied by      
//            the order)                                                
//   rval:    flag to set the actions taken upon finding the term       
//            corresponding to the column cc                            
//-----------------------------------------------------------------------------
static double FindColumnTerm(prmset *mp, int cc, FILE *outp, int *order,
                             int *parmidx, int rval)
{
  int i;

  *order = -1;

  // Leading spaces for output in correlated columns output
  for (i = 0; i < mp->nbond; i++) {
    if (mp->bonds[i].fitcolX == cc) {
      *order = 2;
      *parmidx = i;
      if (rval == 0) {
        fprintf(outp, "  ");
      }
      if (rval <= 0) {
        fprintf(outp, " BOND %.2s %.2s      ", mp->bonds[i].atype,
                mp->bonds[i].btype);
      }
      else if (rval == 1) {
        return  mp->bonds[i].K;
      }
      else if (rval == 2) {
        return  mp->bonds[i].l0;
      }
    }
  }
  for (i = 0; i < mp->nangl; i++) {
    if (mp->angls[i].fitcolX == cc) {
      *order = 3;
      *parmidx = i;
      if (rval == 0) {
        fprintf(outp, "  ");
      }
      if (rval <= 0) {
        fprintf(outp, " ANGL %.2s %.2s %.2s   ", mp->angls[i].atype,
                mp->angls[i].btype, mp->angls[i].ctype);
      }
      else if (rval == 1) {
        return mp->angls[i].K;
      }
      else if (rval == 2) {
        return mp->angls[i].th0;
      }
    }
  }
  for (i = 0; i < mp->ntor; i++) {
    if (mp->torsions[i].fitcol == cc) {
      *order = 4;
      *parmidx = i;
      if (rval == 0) {
        fprintf(outp, "  ");
      }
      if (rval <= 0) {
        if (mp->torsions[i].impr == 0) {
          fprintf(outp, " DIHE ");
        }
        else {
          fprintf(outp, " IMPR ");
        }
        fprintf(outp, "%.2s %.2s %.2s %.2s", mp->torsions[i].atype,
                mp->torsions[i].btype, mp->torsions[i].ctype,
                mp->torsions[i].dtype);
      }
      else if (rval == 1) {
        return mp->torsions[i].K;
      }
      else if (rval == 2) {
        return mp->torsions[i].phase;
      }
      else if (rval == 3) {
        return mp->torsions[i].pn;
      }
    }
  }

  // Carriage return for output relating to term glossary
  if (rval == -1) {
    fprintf(outp, "\n");
  }

  // By default return zero
  return 0.0;
}

//-----------------------------------------------------------------------------
// SamplingMatrixOutput: print the output of the sampling matrix, in    
//                       human-readable format.  Because of referencing between
//                       parameters in the tables made from the input force
//                       field and the columns of the fitting matrix, this
//                       function also logs the fitted values of each parameter
//                       next to their original values.
//                                                                      
// Arguments:                                                           
//   mp:       fitting control data (contains atom types and examples of
//             each atom in various topologies)                         
//   order:    the order of bonded terms to seek out (2 = bonds,        
//             3 = angles, 4 = dihedrals)                               
//   outp:     the output file                                          
//   S:        the sampling matrix                                      
//   x:        the solution vector                                      
//-----------------------------------------------------------------------------
static void SamplingMatrixOutput(prmset *mp, int order, FILE *outp, imat *S,
                                 double* x)
{
  int i, j, parmidx, nvar, ninst, ifitcolX, ifitcolY;
  double Korig, Eqorig, Knew, Eqnew;

  if (order == 2) {
    if (mp->FitBondEq == 1) {
      fprintf(outp, "                                               "
              "   Fitted  Original\n");
      fprintf(outp, "  Term  Amber Types   Count  -0.25   +0    0.25"
              "   K / Eq   K / Eq \n");
    }
    else {
      fprintf(outp, "  Term  Amber Types   Count  -0.25   +0    0.25"
              "   Fitted  Original\n");
    }
  }
  else if (order == 3) {
    if (mp->FitAnglEq == 1) {
      fprintf(outp, "                                               "
              "   Fitted  Original\n");
      fprintf(outp, "  Term  Amber Types   Count  -15     +0      15"
              "   K / Eq   K / Eq \n");
    }
    else {
      fprintf(outp, "  Term  Amber Types   Count  -15     +0      15"
              "   Fitted  Original\n");
    }
  }
  else if (order == 4) {
    fprintf(outp, "  Term  Amber Types   Count  -PI     +0      PI"
            "   Fitted  Original\n");
  }
  fprintf(outp, "  ----  -----------  ------- "
          "------------------  -------- --------\n");
  nvar = 0;
  for (i = 0; i < S->row; i++) {
    if (order == 2) {
      while (mp->bonds[nvar].fitcolX < 0) {
        nvar++;
      }
      ninst = mp->bonds[nvar].ninst;
      fprintf(outp, "  BOND  %.2s %.2s       ", mp->bonds[nvar].atype,
              mp->bonds[nvar].btype);
    }
    if (order == 3) {
      while (mp->angls[nvar].fitcolX < 0) {
        nvar++;
      }
      ninst = mp->angls[nvar].ninst;
      fprintf(outp, "  ANGL  %.2s %.2s %.2s    ", mp->angls[nvar].atype,
              mp->angls[nvar].btype, mp->angls[nvar].ctype);
    }
    if (order == 4) {
      while (mp->torsions[nvar].fitcol < 0) {
        nvar++;
      }
      ninst = mp->torsions[nvar].ninst;
      if (mp->torsions[nvar].impr == 0) {
        fprintf(outp, "  DIHE  ");
      }
      else {
        fprintf(outp, "  IMPR  ");
      }
      fprintf(outp, "%.2s %.2s %.2s %.2s ", mp->torsions[nvar].atype,
              mp->torsions[nvar].btype, mp->torsions[nvar].ctype,
              mp->torsions[nvar].dtype);
    }
    fprintf(outp, " %7d ", ninst);
    CharacterHistogram(S->map[i], S->col, outp);
    if ((order == 2 && mp->FitBondEq == 0) ||
        (order == 3 && mp->FitAnglEq == 0) || order == 4) {
      ifitcolX = MatchSamprow2Column(mp, i, order, 0);
      Korig = FindColumnTerm(mp, ifitcolX, outp, &j, &parmidx, 1);
      if (order == 2) {
        mp->bonds[nvar].Kfit = x[ifitcolX];
      }
      else if (order == 3) {
        mp->angls[nvar].Kfit = x[ifitcolX];
      }
      else {
        mp->torsions[nvar].Kfit = x[ifitcolX];
      }
      fprintf(outp, "  %8.2lf %8.2lf\n", x[ifitcolX], Korig);
    }
    else if ((order == 2 && mp->FitBondEq == 1) ||
             (order == 3 && mp->FitAnglEq == 1)) {
      ifitcolX = MatchSamprow2Column(mp, i, order, 0);
      Korig = FindColumnTerm(mp, ifitcolX, outp, &j, &parmidx, 1);
      ifitcolY = MatchSamprow2Column(mp, i, order, 1);
      Eqorig = FindColumnTerm(mp, ifitcolX, outp, &j, &parmidx, 2);
      Knew = x[ifitcolX] + x[ifitcolY];
      fprintf(outp, "  %8.2lf %8.2lf\n", Knew, Korig);
      if (order == 2) {
        Eqnew = (mp->bonds[nvar].lbasisX * x[ifitcolX] +
                 mp->bonds[nvar].lbasisY * x[ifitcolY]) / Knew;
        mp->bonds[nvar].Kfit = Knew;
        mp->bonds[nvar].l0fit = Eqnew;
      }
      else {
        Eqnew = (mp->angls[nvar].tbasisX * x[ifitcolX] +
                 mp->angls[nvar].tbasisY * x[ifitcolY]) / Knew;
        mp->angls[nvar].Kfit = Knew;
        mp->angls[nvar].th0fit = Eqnew;
      }
      if (order == 2) {
        fprintf(outp, "                                               "
                "  %8.2lf %8.2lf\n", Eqnew, Eqorig);
      }
      else {
        fprintf(outp, "                                               "
                "  %8.2lf %8.2lf\n", Eqnew * 180.0 / PI, Eqorig * 180.0 / PI);
      }
    }
    nvar++;
  }
  fprintf(outp, "\n");
}

//-----------------------------------------------------------------------------
// DetailParameter: details the instances in which a parameter occurs in the
//                  fitting set.                                    
//                                                                      
// Arguments:                                                           
//   mp:       fitting control data (contains atom types and examples of
//             each atom in various topologies)                         
//   n:        the column of the parameter to detail                    
//   outp:     the output file                                          
//-----------------------------------------------------------------------------
static void DetailParameter(prmset *mp, int n, FILE *outp)
{
  int h, i, j, bkey, akey, hkey, ninst, atma, atmb, atmc, atmd;
  int resa, resb, resc, resd, nspc, rescon, bcol, acol, hcol;
  char *tpatoms, *tpres;
  imat* BondSamp;
  imat* AnglSamp;
  imat* DiheSamp;
  dmat *SampVtmp;
  itrack* inst;

  // Determine the maximum possible number of instances
  j = 0;
  for (i = 0; i < mp->nconf; i++) {
    j += mp->conf[i].bmap.nbond;
    j += mp->conf[i].amap.nangl;
    j += mp->conf[i].hmap.ntterm;
  }
  inst = (itrack*)malloc(j*sizeof(itrack));

  // Check for this term in exactly one example of each system
  ninst = 0;
  for (h = 0; h < mp->nunisys; h++) {
    i = mp->FirstConf[h];
    tpatoms = mp->conf[i].tp->AtomNames;
    tpres = mp->conf[i].tp->ResNames;
    for (j = 0; j < mp->conf[i].bmap.nbond; j++) {
      bkey = mp->conf[i].bmap.id[j].key;
      if (mp->bonds[bkey].fitcolX == n) {

        // Column n pertains to a bond in this molecule
        inst[ninst].sysid = mp->conf[i].GroupNum;
        inst[ninst].sysno = i;
        inst[ninst].order = 2;
        inst[ninst].tnum = j;
        atma = mp->conf[i].bmap.id[j].a;
        atmb = mp->conf[i].bmap.id[j].b;
        resa = LocateResID(mp->conf[i].tp, atma, 0, mp->conf[i].tp->nres);
        resb = LocateResID(mp->conf[i].tp, atmb, 0, mp->conf[i].tp->nres);
        sprintf(inst[ninst].atom, "%.4s%.4s", &tpatoms[4*atma],
                &tpatoms[4*atmb]);
        sprintf(inst[ninst].res, "%.4s%.4s", &tpres[4*resa], &tpres[4*resb]);
        ninst++;
      }
    }
    for (j = 0; j < mp->conf[i].amap.nangl; j++) {
      akey = mp->conf[i].amap.id[j].key;
      if (mp->angls[akey].fitcolX == n) {

        // Column n pertains to an angle in this molecule
        inst[ninst].sysid = mp->conf[i].GroupNum;
        inst[ninst].sysno = i;
        inst[ninst].order = 3;
        inst[ninst].tnum = j;
        atma = mp->conf[i].amap.id[j].a;
        atmb = mp->conf[i].amap.id[j].b;
        atmc = mp->conf[i].amap.id[j].c;
        resa = LocateResID(mp->conf[i].tp, atma, 0, mp->conf[i].tp->nres);
        resb = LocateResID(mp->conf[i].tp, atmb, 0, mp->conf[i].tp->nres);
        resc = LocateResID(mp->conf[i].tp, atmc, 0, mp->conf[i].tp->nres);
        sprintf(inst[ninst].atom, "%.4s%.4s%.4s", &tpatoms[4*atma],
                &tpatoms[4*atmb], &tpatoms[4*atmc]);
        sprintf(inst[ninst].res, "%.4s%.4s%.4s", &tpres[4*resa],
                &tpres[4*resb], &tpres[4*resc]);
        ninst++;
      }
    }
    for (j = 0; j < mp->conf[i].hmap.ntterm; j++) {
      hkey = mp->conf[i].hmap.id[j].key;
      if (mp->torsions[hkey].fitcol == n) {

        // Column n pertains to a dihedral in this molecule
        inst[ninst].sysid = mp->conf[i].GroupNum;
        inst[ninst].sysno = i;
        inst[ninst].order = 4;
        inst[ninst].tnum = j;
        atma = mp->conf[i].hmap.id[j].a;
        atmb = mp->conf[i].hmap.id[j].b;
        atmc = mp->conf[i].hmap.id[j].c;
        atmd = mp->conf[i].hmap.id[j].d;
        resa = LocateResID(mp->conf[i].tp, atma, 0, mp->conf[i].tp->nres);
        resb = LocateResID(mp->conf[i].tp, atmb, 0, mp->conf[i].tp->nres);
        resc = LocateResID(mp->conf[i].tp, atmc, 0, mp->conf[i].tp->nres);
        resd = LocateResID(mp->conf[i].tp, atmd, 0, mp->conf[i].tp->nres);
        sprintf(inst[ninst].atom, "%.4s%.4s%.4s%.4s", &tpatoms[4*atma],
                &tpatoms[4*atmb], &tpatoms[4*atmc], &tpatoms[4*atmd]);
        sprintf(inst[ninst].res, "%.4s%.4s%.4s%.4s", &tpres[4*resa],
                &tpres[4*resb], &tpres[4*resc], &tpres[4*resd]);
        ninst++;
      }
    }
  }

  // Create sampling tables for conformations 
  // sampled by each individual system.       
  if (mp->nbvar > 0) BondSamp = (imat*)malloc(mp->nunisys*sizeof(imat));
  if (mp->navar > 0) AnglSamp = (imat*)malloc(mp->nunisys*sizeof(imat));
  if (mp->nhvar > 0) DiheSamp = (imat*)malloc(mp->nunisys*sizeof(imat));
  MapSamplingIndices(mp, 2);
  MapSamplingIndices(mp, 3);
  MapSamplingIndices(mp, 4);
  for (i = 0; i < mp->nunisys; i++) {
    if (mp->nbvar > 0) {
      BondSamp[i] = CreateImat(mp->nbvar, 18);
      AccumulateSamplingTable(mp, &BondSamp[i], SampVtmp, 2, i);
    }
    if (mp->navar > 0) {
      AnglSamp[i] = CreateImat(mp->navar, 18);
      AccumulateSamplingTable(mp, &AnglSamp[i], SampVtmp, 3, i);
    }
    if (mp->nhvar > 0) {
      DiheSamp[i] = CreateImat(mp->nhvar, 18);
      AccumulateSamplingTable(mp, &DiheSamp[i], SampVtmp, 4, i);
    }
  }

  // Recount the instances
  for (i = 0; i < ninst; i++) {
    fprintf(outp, "  ");
    fprintf(outp, "%.4s(", inst[i].res);
    rescon = 0;
    nspc = 7;
    for (j = 0; j < inst[i].order; j++) {
      if (str4cmp(&inst[i].res[4*j], &inst[i].res[4*rescon]) == 0) {
        fprintf(outp, "%.4s", &inst[i].atom[4*j]);
        nspc += 4;
      }
      else {
        fprintf(outp, ") %.4s(%.4s", &inst[i].res[4*j], &inst[i].atom[4*j]);
        rescon = j;
        nspc += 11;
      }
      if (j < inst[i].order-1 &&
          str4cmp(&inst[i].res[4*(j+1)], &inst[i].res[4*rescon]) == 0) {
        fprintf(outp, " ");
        nspc += 1;
      }
    }
    fprintf(outp, ")");
    nspc += 1;
    for (j = nspc; j < 41; j++) {
      fprintf(outp, " ");
    }
    fprintf(outp, "%5d  ", mp->GroupCount[inst[i].sysid]);
    if (inst[i].order == 2) {
      bkey = mp->conf[inst[i].sysno].bmap.id[inst[i].tnum].key;
      bcol = mp->bonds[bkey].samprow;
      CharacterHistogram(BondSamp[inst[i].sysid].map[bcol], 18, outp);
    }
    else if (inst[i].order == 3) {
      akey = mp->conf[inst[i].sysno].amap.id[inst[i].tnum].key;
      acol = mp->angls[akey].samprow;
      CharacterHistogram(AnglSamp[inst[i].sysid].map[acol], 18, outp);
    }
    else if (inst[i].order == 4) {
      hkey = mp->conf[inst[i].sysno].hmap.id[inst[i].tnum].key;
      hcol = mp->torsions[hkey].samprow;
      CharacterHistogram(DiheSamp[inst[i].sysid].map[hcol], 18, outp);
    }
    fprintf(outp, "\n");
  }

  // Free allocated memory
  for (i = 0; i < mp->nunisys; i++) {
    if (mp->nbvar > 0) DestroyImat(&BondSamp[i]);
    if (mp->navar > 0) DestroyImat(&AnglSamp[i]);
    if (mp->nhvar > 0) DestroyImat(&DiheSamp[i]);
  }
  if (mp->nbvar > 0) free(BondSamp);
  if (mp->navar > 0) free(AnglSamp);
  if (mp->nhvar > 0) free(DiheSamp);
  free(inst);
}

//-----------------------------------------------------------------------------
// WarningsOutput: runs a series of diagnostics on the parameter and its   
//                 possible effects on the resulting model.         
//                                                                      
// Arguments:                                                           
//   mp:       fitting control data (contains atom types and examples of
//             each atom in various topologies)                         
//   sampV:    the sampling histogram, binned with Gaussian quadrature  
//             over a preset range (720 bins, with padding for angle and
//             dihedral angle measurements)                             
//   order:    the order of the parameters we've analyzed (2 = bonds,   
//             3 = angles, 4 = dihedrals, no support yet for impropers) 
//   outp:     the output file                                          
//-----------------------------------------------------------------------------
static cmat WarningsOutput(prmset *mp, dmat *sampV, double* x, int order,
                           FILE *outp)
{
  int i, j, idxbin, jmin, jmax, nwarn, maxln, nln;
  double eqval, hsum, sumbelow, sumabove, sumnear;
  double *dtmp;
  cmat wtext;

  // Pre-allocate wtext to hold whatever output we generate
  maxln = 128;
  wtext = CreateCmat(maxln, 128);

  // Normalize the fine-grained histograms
  for (i = 0; i < sampV->row; i++) {
    dtmp = sampV->map[i];
    hsum = 1.0/DSum(dtmp, sampV->col);
    for (j = 0; j < sampV->col; j++) {
      dtmp[j] *= hsum;
    }
  }

  // Loop over all parameters and check for problems
  nln = 0;
  nwarn = 0;
  if (order == 2) {
    for (i = 0; i < mp->nbond; i++) {

      // Skip if this bond was not optimizable
      if (mp->bonds[i].fitcolX < 0) {
        continue;
      }

      // Pre-emptively allocate new memory
      if (nln > maxln-3) {
        maxln += 32;
        wtext = ReallocCmat(&wtext, maxln, 128);
      }

      // Check for optimization of the equilibrium constant.   
      // If the equilibrium was not optimized, report only     
      // if the new stiffness constant is drastically different
      // from the original.                                    
      if (mp->bonds[i].fitcolY < 0) {
        if (fabs(mp->bonds[i].Kfit -
                 mp->bonds[i].K) / fabs(mp->bonds[i].K) > 0.5) {
          sprintf(wtext.map[nln], "  - Bond %2.2s %2.2s stiffness changes "
                  "from %9.4lf to %9.4lf", mp->bonds[i].atype,
                  mp->bonds[i].btype, mp->bonds[i].K, mp->bonds[i].Kfit);
          nln++;
        }
        nwarn++;
        continue;
      }

      // If we're still here, check the population of     
      // samples at and around the new equilibrium value. 
      idxbin = mp->bonds[i].l0fit * 100.0;
      if (idxbin < 5 || idxbin >= 714) {
        sprintf(wtext.map[nln], "  - Bond %2.2s %2.2s has an extreme "
                "equilibrium value (%9.4lf)", mp->bonds[i].atype,
                mp->bonds[i].btype, mp->bonds[i].l0fit);
        nln++;
        nwarn++;
        continue;
      }
      sumbelow = DSum(sampV->map[mp->bonds[i].samprow], idxbin);
      sumabove = DSum(&sampV->map[mp->bonds[i].samprow][idxbin+1], 719-idxbin);
      sumnear = DSum(&sampV->map[mp->bonds[i].samprow][idxbin-5], 11);
      if (sumbelow < 0.2 || sumabove < 0.2) {
        sprintf(wtext.map[nln], "  - Bond %2.2s %2.2s has lopsided sampling "
                "(%6.2lf%% below, %6.2lf%% above", mp->bonds[i].atype,
                mp->bonds[i].btype, 100.0*sumbelow, 100.0*sumabove);
        nln++;
        sprintf(wtext.map[nln], "    the optimized equilibrium value of "
                "%9.4lf).", mp->bonds[i].l0fit);
        nln++;
        nwarn++;
      }
      if (sumnear < 0.05) {
        sprintf(wtext.map[nln], "  - Bond %2.2s %2.2s has poor sampling "
                "within 0.05A of its optimized", mp->bonds[i].atype,
                mp->bonds[i].btype);
        nln++;
        sprintf(wtext.map[nln], "    equilbrium value (%9.4lf).",
                mp->bonds[i].l0fit);
        nln++;
        nwarn++;
      }
    }
  }
  if (order == 3) {
    for (i = 0; i < mp->nangl; i++) {

      // Skip if this angle was not optimizable
      if (mp->angls[i].fitcolX < 0) {
        continue;
      }

      // Check for optimization of the equilibrium constant.   
      // If the equilibrium was not optimized, report only     
      // if the new stiffness constant is drastically different
      // from the original.                                    
      if (mp->angls[i].fitcolY < 0) {
        if (fabs(mp->angls[i].Kfit -
                 mp->angls[i].K) / fabs(mp->angls[i].K) > 0.25) {
          if (nln > maxln-1) {
            maxln += 32;
            wtext = ReallocCmat(&wtext, maxln, 128);
          }
          sprintf(wtext.map[nln], "  - Angle %2.2s %2.2s %2.2s stiffness "
                  "changes from %9.4f to %9.4f", mp->angls[i].atype,
                  mp->angls[i].btype, mp->angls[i].ctype, mp->angls[i].K,
                  mp->angls[i].Kfit);
          nwarn++;
          nln += 1;
        }
        continue;
      }

      // If we're still here, check the population of     
      // samples at and around the new equilibrium value. 
      idxbin = mp->angls[i].th0fit * 180.0 / PI + 360.0;
      if (idxbin < 5 || idxbin >= 714) {
        sprintf(wtext.map[nln], "  - Angle %2.2s %2.2s %2.2s has an extreme "
                "equilibrium value (%9.4lf)", mp->angls[i].atype,
                mp->angls[i].btype, mp->angls[i].ctype, mp->angls[i].th0fit);
        nln++;
        nwarn++;
        continue;
      }
      sumbelow = DSum(sampV->map[mp->angls[i].samprow], idxbin);
      sumabove = DSum(&sampV->map[mp->angls[i].samprow][idxbin+1], 719-idxbin);
      sumnear = DSum(&sampV->map[mp->angls[i].samprow][idxbin-3], 7);
      if (sumbelow < 0.2 || sumabove < 0.2) {
        sprintf(wtext.map[nln], "  - Angle %2.2s %2.2s %2.2s has lopsided "
                "sampling (%6.2lf%% below, %6.2lf%% above", mp->angls[i].atype,
                mp->angls[i].btype, mp->angls[i].ctype, 100.0*sumbelow,
                100.0*sumabove);
        nln++;
        sprintf(wtext.map[nln], "    the optimized equilibrium value of "
                "%9.4lf).", mp->angls[i].th0fit * 180.0/PI);
        nln++;
        nwarn++;
      }
      if (sumnear < 0.05) {
        sprintf(wtext.map[nln], "  - Angle %2.2s %2.2s %2.2s has poor "
                "sampling within 0.05 radians of its optimized",
                mp->angls[i].atype, mp->angls[i].btype, mp->angls[i].ctype);
        nln++;
        sprintf(wtext.map[nln], "    equilbrium value (%9.4lf).",
                mp->angls[i].th0fit * 180.0/PI);
        nln++;
        nwarn++;
      }
    }
  }
  if (order == 4) {
    for (i = 0; i < mp->ntor; i++) {

      // Skip if this torsion was not optimizable
      if (mp->torsions[i].fitcol < 0) {
        continue;
      }

      // Check for radical changes to the stiffness constants.  
      // This is a bit trickier than with bonds and angles--the 
      // torsions could be sums of multiple Fourier series terms
      // and the amplitudes of individual terms could change in 
      // radical ways which have little effect on the overall   
      // shape of the potential.  So, test the potential over a 
      // wide range, see what happens, and report cases in which
      // the overall potential has changed.  Also report minima 
      // that appear in the parameter set but are not sampled in
      // the data set.  But, we're not done yet: if the torsion 
      // is stretched across an improper dihedral in all cases, 
      // then SO LONG AS its summed potential cannot overwhelm  
      // the impropers minima in the excluded regions can be    
      // ignored, as can large changes in the torsion potential 
      // in these regions.                                      
    }
  }

  // Check whether there were any warnings at all,
  // and print a title for this block of warnings.
  if (nln > 0) {
    wtext = ReallocCmat(&wtext, nln+2, 128);
    for (i = nln-1; i >= 0; i--) {
      for (j = 0; j < 128; j++) {
        wtext.map[i+1][j] = wtext.map[i][j];
      }
    }
    if (order == 2) {
      sprintf(wtext.map[0], " Bond parameter warnings (%d):", nwarn);
    }
    else if (order == 3) {
      sprintf(wtext.map[0], " Angle parameter warnings (%d):", nwarn);
    }
    else if (order == 4) {
      sprintf(wtext.map[0], " Dihedral parameter warnings (%d):", nwarn);
    }
    for (i = 0; i < 128; i++) {
      wtext.map[nln+1][i] = '\0';
    }
  }
  else {
    wtext = ReallocCmat(&wtext, 1, 1);
    wtext.data[0] = '\0';
  }

  return wtext;
}

//-----------------------------------------------------------------------------
// PrintSamplingTable: print a table to describe the sampling of each fitted  
//                     parameter across all conformations.       
//                                                                      
// Arguments:                                                           
//   mp:       fitting control data (contains atom types and examples of
//             each atom in various topologies)                         
//   order:    the order of bonded terms to seek out (2 = bonds,        
//             3 = angles, 4 = dihedrals)                               
//   x:        the solution vector (to print alongside sampling data)   
//   outp:     the output file                                          
//-----------------------------------------------------------------------------
static cmat PrintSamplingTable(prmset *mp, int order, double* x, FILE *outp)
{
  int nvar;
  char desc[128];
  imat samp;
  dmat sampV;
  cmat warnings;

  // Set flags and strings
  if (order == 2) {
    nvar = (mp->FitBondEq == 1) ? mp->nbvar / 2 : mp->nbvar;
    sprintf(desc, " Bond sampling:\n"
            "                             Bins in Angstroms\n");
  }
  else if (order == 3) {
    nvar = (mp->FitAnglEq == 1) ? mp->navar / 2 : mp->navar;
    sprintf(desc, " Angle sampling:\n"
            "                              Bins in Degrees \n");
  }
  else if (order == 4) {
    nvar = mp->nhvar;
    sprintf(desc, " Torsion term sampling:\n"
            "                              Bins in Radians \n");
  }

  // Allocate the histogram matrix, accumulate the histogram,
  // print the results, then destroy the histogram matrix.   
  if (nvar > 0) {
    samp = CreateImat(nvar, 18);
    sampV = CreateDmat(nvar, 720, 0);
    MapSamplingIndices(mp, order);
    AccumulateSamplingTable(mp, &samp, &sampV, order, -1);
    fprintf(outp, "%s", desc);
    SamplingMatrixOutput(mp, order, outp, &samp, x);
    DestroyImat(&samp);
    warnings = WarningsOutput(mp, &sampV, x, order, outp);
    DestroyDmat(&sampV);
  }
  else {
    warnings = CreateCmat(1, 1);
  }

  return warnings;
}

//-----------------------------------------------------------------------------
// PrintCMAPSampling: print the sampling for individual CMAP terms.
//
// Arguments:
//   mp:       fitting control data (contains atom types and examples of
//             each atom in various topologies)
//   outp:     the output file
//-----------------------------------------------------------------------------
static void PrintCMAPSampling(prmset *mp, FILE *outp)
{
  int i, j, k, cmidx, col1, col2;
  double xval, yval;
  char xlabel[64], ylabel[64];
  imat* cstable;
  
  // Accumulate sampling histograms
  cstable = (imat*)malloc(mp->ncmaps * sizeof(imat));
  for (i = 0; i < mp->ncmaps; i++) {
    cstable[i] = CreateImat(30, 54);
  }
  for (i = 0; i < mp->nconf; i++) {
    for (j = 0; j < mp->conf[i].mmap.ncmap; j++) {
      cmidx = mp->conf[i].mmap.id[j].key;
      xval = mp->conf[i].mmap.val[2*j];
      yval = mp->conf[i].mmap.val[2*j + 1];
      xval /= 2.0 * PI;
      yval /= 2.0 * PI;
      if (xval < -0.5) {
	xval -= floor(xval);
      }
      else if (xval >= 0.5) {
	xval -= ceil(xval);
      }
      if (yval < -0.5) {
	yval -= floor(yval);
      }
      else if (yval >= 0.5) {
	yval -= ceil(yval);
      }
      col1 = (int)((yval + 0.5) * 30.0);
      col2 = (int)(-(xval + 0.5) * 54.0);
      col1 += ((col1 < 0) - (col1 >= 30)) * 30;
      col2 += ((col2 < 0) - (col2 >= 54)) * 54;
      cstable[cmidx].map[col1][col2] += 1;
    }
  }

  // Write out the tables
  sprintf(ylabel, "    SECOND DIHEDRAL ANGLE     ");
  sprintf(xlabel, "                              FIRST DIHEDRAL ANGLE");
  for (i = 0; i < mp->ncmaps; i++) {

    // Do not print CMAPs that were not fitted
    if (mp->cmaps[i].fitcoli < 0) {
      continue;
    }

    // Print the description, then the sampling histogram
    fprintf(outp, "  CMAP Surface %2d (number according to input frcmod "
	    "files):\n", i+1);
    fprintf(outp, "    Title: %s\n", mp->cmaps[i].title);
    fprintf(outp, "    Residues:");
    for (j = 0; j < mp->cmaps[i].resnames.row; j++) {
      fprintf(outp, " %4.4s", mp->cmaps[i].resnames.map[j]);
    }
    fprintf(outp, "\n    Atoms (residue #): ");
    for (j = 0; j < 5; j++) {
      for (k = 4*j; k < 4*(j+1); k++) {
	if (mp->cmaps[i].atoms[k] != ' ') {
	  fprintf(outp, "%c", mp->cmaps[i].atoms[k]);
	}
      }
      if (mp->cmaps[i].relResID[j] == 0) {
        fprintf(outp, "(0)");
      }
      else {
        fprintf(outp, "(%+d)",  mp->cmaps[i].relResID[j]);
      }
      if (j < 4) {
	fprintf(outp, " -- ");
      }
    }
    fprintf(outp, "\n\n        180 +");
    for (j = 0; j < 54; j++) {
      fprintf(outp, "-");
    }
    fprintf(outp, "+\n");    
    for (j = 0; j < 30; j++) {
      if (j == 14) {
        fprintf(outp, "       %c  0 |", ylabel[j]);
      }
      else {
        fprintf(outp, "       %c    |", ylabel[j]);
      }
      CharacterHistogram(cstable[i].map[j], 54, outp);
      fprintf(outp, "|\n");
    }
    fprintf(outp, "       -180 +");
    for (j = 0; j < 54; j++) {
      fprintf(outp, "-");
    }
    fprintf(outp, "+\n");
    fprintf(outp, "          -180");
    fprintf(outp, "                          0                         180\n");
    fprintf(outp, "%s\n\n", xlabel);    
  }

  // Free allocated memory
  for (i = 0; i < mp->ncmaps; i++) {
    DestroyImat(&cstable[i]);
  }
}

//-----------------------------------------------------------------------------
// PrintEAConstants: print the energy adjustment constants emerging from the
//                   fit.  Large values of these constants can indicate that
//                   certain fitted parameters are being allowed to contribute
//                   large baseline values to the energy estimates.
//
// Arguments:                                                           
//   mp:       fitting control data (contains atom types and examples of
//             each atom in various topologies)                         
//   x:        the solution vector (to print alongside sampling data)   
//   outp:     the output file                                          
//-----------------------------------------------------------------------------
static void PrintEAConstants(prmset *mp, double* x, FILE *outp)
{
  int h, i, j, k, vstart;
  double maxwt;

  if (mp->nunisys == 1) {
    if (mp->wtfloor > 100.0) {
      fprintf(outp,
              "      System      Adjustment\n"
              " ---------------- ----------");
    }
    else {
      fprintf(outp,
              "      System      Adjustment  Weight\n"
              " ---------------- ---------- -------");
    }
  }
  else {
    if (mp->wtfloor > 100.0) {
      fprintf(outp,
              "      System      Adjustment         System      "
              "Adjustment\n"
              " ---------------- ----------    ---------------- "
              "----------");
    }
    else {
      fprintf(outp,
              "      System      Adjustment  Weight         System      "
              "Adjustment  Weight\n"
              " ---------------- ---------- -------    ---------------- "
              "---------- -------");
    }
  }
  h = 2;
  vstart = mp->nparm;
  for (i = 0; i < mp->nunisys; i++) {
    if (h == 2) {
      h = 0;
      fprintf(outp, "\n ");
    }
    for (j = 0; j < mp->nconf; j++) {
      if (mp->conf[j].GroupNum == i) {
        if (mp->wtfloor > 100.0) {
          fprintf(outp, "%-16.16s %10.4lf    ", mp->conf[j].tp->source,
                  x[vstart]);
        }
        else {
          maxwt = 0.0;
          for (k = 0; k < mp->nconf; k++) {
            if (mp->conf[k].GroupNum == i && mp->conf[k].wt > maxwt) {
              maxwt = mp->conf[k].wt;
            }
          }
          fprintf(outp, "%-16.16s %10.4lf %7.4lf    ", mp->conf[j].tp->source,
                  x[vstart], maxwt);
        }
        h++;
        break;
      }
    }
    vstart++;
  }
  if (h != 0) {
    fprintf(outp, "\n");
  }
  fprintf(outp, "\n");
  if (mp->fitscee == 1 || mp->fitscnb == 1) {
    PrintVADesc(0, " ", 1, " ", 1, "In addition, 1:4 scaling factors were "
                "also fitted, and should be applied in any simulations with "
                "the resulting parameters.\n", 77, 0, outp);
    if (mp->fitscee == 1) {
      fprintf(outp, " - SCEE = %10.6lf\n",
              1.0/x[mp->nbvar+mp->navar+mp->nhvar]);
    }
    if (mp->fitscnb == 1) {
      fprintf(outp, " - SCNB = %10.6lf\n\n",
              1.0/x[mp->nbvar+mp->navar+mp->nhvar+mp->fitscee]);
    }
  }
}

//-----------------------------------------------------------------------------
// PrintMatrixAnalysis: print an analysis of the fitting matrix, indentifying
//                      highly correlated columns as well as columns with no
//                      fitting data in them.           
//
// Arguments:                                                           
//   mp:      fitting control data (contains notes about the fitting    
//            matrix)                                                   
//   x:       the solution vector                                       
//   outp:    the output file                                           
//-----------------------------------------------------------------------------
static void PrintMatrixAnalysis(prmset *mp, double* x, FILE *outp)
{
  int i, order, parmidx;

  // Zero-data columns?
  if (mp->nzerocol == 0) {
    fprintf(outp, " - No columns with zero data detected\n");
  }
  else {
    fprintf(outp, " - %d columns with zero data detected, corresponding to "
            "parameters:\n", mp->nzerocol);
    for (i = 0; i < mp->nzerocol; i++) {
      FindColumnTerm(mp, mp->zerocol[i], outp, &order, &parmidx, 0);
      fprintf(outp, "\n");
    }
  }

  // Correlated columns?
  if (mp->ncorrcol == 0) {
    fprintf(outp, " - No correlated columns detected\n");
  }
  else {
    fprintf(outp, " - %d highly correlated column pairs detected, "
            "corresponding to parameters:\n\n", mp->ncorrcol);
    fprintf(outp, "     First term,        Second term,     Fitted   Fitted   "
            "   Pearson\n");
    fprintf(outp, "   Type    Atoms      Type    Atoms      value 1  value 2  "
            " Correlation\n");
    fprintf(outp, "   ---- -----------   ---- -----------   -------- -------- "
            " -----------\n");
    for (i = 0; i < mp->ncorrcol; i++) {
      FindColumnTerm(mp, mp->corrcol[2*i], outp, &order, &parmidx, 0);
      FindColumnTerm(mp, mp->corrcol[2*i+1], outp, &order, &parmidx, 0);
      fprintf(outp, "   %8.2lf %8.2lf  %11.8lf\n",
              x[mp->corrcol[2*i]], x[mp->corrcol[2*i+1]], mp->corrval[i]);
    }
  }
}

//-----------------------------------------------------------------------------
// ParameterDescriptions: this routine will print information to put all fitted
//                        parameters in context.  It details the residues and 
//                        atoms that each parameter affects and gives counts on
//                        the number of instances of each occurrence.  This
//                        routine will also check for criteria that may
//                        indicate a parameter has been poorly fitted.  Output
//                        is tabulated in a special section of the mdout file.
//                                                                      
// Arguments:                                                           
//   mp:       fitting control data (contains atom types and examples of
//             each atom in various topologies)                         
//   outp:     pointer to the mdout file                                
//-----------------------------------------------------------------------------
static void ParameterDescription(prmset *mp, FILE *outp)
{
  int i, order, parmidx;

  fprintf(outp, " [ Parameter ] [ Atom types ]                   Sampling "
          "Histogram\n"
          "  Residue (Atom Names)                   Count  MIN     +0     MAX "
          " Warnings\n"
          " --------------------------------------- -----  ------------------ "
          " --------\n");
  for (i = 0; i < mp->nparm; i++) {
    if (mp->verbose == 1) {
      fprintf(stderr, "\rmdgx >> Detailing parameter %4d of %4d.", i,
              mp->nparm);
      fflush(stderr);
    }
    FindColumnTerm(mp, i, outp, &order, &parmidx, -1);
    DetailParameter(mp, i, outp);
  }
  if (mp->verbose == 1) {
    printf("\rmdgx >> Descriptions complete for all instances of %4d "
           "parameters.\n", mp->nparm);
  }
}

//-----------------------------------------------------------------------------
// CountFittedTerms: count the number of fitted terms in a system.      
//
// Arguments:                                                           
//   sysno:    the system number                                        
//   mp:       the master parameter and fitting set                     
//-----------------------------------------------------------------------------
static int CountFittedTerms(int sysno, prmset *mp)
{
  int i, exid, psum;
  int* pterm;
  mmsys *exconf;

  // First, find one example of this system
  exid = -1;
  for (i = 0; i < mp->nconf; i++) {
    if (mp->conf[i].GroupNum == sysno) {
      exid = i;
      break;
    }
  }
  if (exid == -1) {
    printf("CountFittedTerms >> Error.  System %d does not exist.\n", sysno);
    exit(1);
  }
  exconf = &mp->conf[exid];

  // Bonds, angles, and dihedrals
  pterm = (int*)calloc(mp->nparm, sizeof(int));
  for (i = 0; i < exconf->bmap.nbond; i++) {
    if (mp->bonds[exconf->bmap.id[i].key].fitcolX >= 0) {
      pterm[mp->bonds[exconf->bmap.id[i].key].fitcolX] = 1;
    }
    if (mp->bonds[exconf->bmap.id[i].key].fitcolY >= 0) {
      pterm[mp->bonds[exconf->bmap.id[i].key].fitcolY] = 1;
    }
  }
  for (i = 0; i < exconf->amap.nangl; i++) {
    if (mp->angls[exconf->amap.id[i].key].fitcolX >= 0) {
      pterm[mp->angls[exconf->amap.id[i].key].fitcolX] = 1;
    }
    if (mp->angls[exconf->amap.id[i].key].fitcolY >= 0) {
      pterm[mp->angls[exconf->amap.id[i].key].fitcolY] = 1;
    }
  }
  for (i = 0; i < exconf->hmap.ntterm; i++) {
    if (mp->torsions[exconf->hmap.id[i].key].fitcol >= 0) {
      pterm[mp->torsions[exconf->hmap.id[i].key].fitcol] = 1;
    }
  }
  psum = ISum(pterm, mp->nparm);

  // Free allocated memory
  free(pterm);

  return psum;
}

//-----------------------------------------------------------------------------
// PrintAccuracy: report the accuracy of this model in relation to the  
//                training set.
//
// Arguments:                                                           
//   mp:       fitting control data (contains atom types and examples of
//             each atom in various topologies)                         
//   A:        the fitting matrix, in its original form before the QR   
//             decomposition                                            
//   x:        the solution vector, containing all fitted parameters    
//   tj:       trajectory control information                           
//   outp:     the output file                                          
//-----------------------------------------------------------------------------
static void PrintAccuracy(prmset *mp, dmat *A, double* x, trajcon *tj,
                          FILE *outp)
{
  int i, j, k, m, ncfg, fsys, natom, noutlier;
  int* nsysvar;
  int* confid;
  double fitnrg, eQQ, eLJ, fitrmse, fintol;
  double *atmp;
  double* enorm;
  double* eorig;
  double* efin;
  double* eelec;
  double* elj;
  double* ebond;
  double* eangl;
  double* edihe;
  double* ecmap;
  double* enmr;
  double* efit;
  double* badnrg;
  char** badconf;
  FILE *mscript;
  time_t ct;

  // If a MatLab readable display is requested, print it
  if (mp->ao[0] != '\0') {
    mscript = FOpenSafe(mp->ao, tj->OverwriteOutput);
    ct = time(&ct);
    fprintf(mscript, "%% Model accuracy report for fitting ordered by\n"
            "%% %s.\n", tj->inpname);
    fprintf(mscript, "%% Written on %s\n", asctime(localtime(&ct)));
  }
  fprintf(outp, 
          "                      RMS Error         Correlation    Fitted   "
          "Fitted   Model\n"
          "      System        Orig    Fitted    Orig    Fitted   Terms    "
          "Energy   Diff.\n"
          " ----------------  -------  -------  -------  -------  ------  "
          "-------  -------\n");

  // Count the number of fitted terms in each system
  nsysvar = (int*)malloc(mp->nunisys*sizeof(int));
  for (i = 0; i < mp->nunisys; i++) {
    nsysvar[i] = CountFittedTerms(i, mp);
  }

  // Compute the accuracy of the model
  confid = (int*)malloc(mp->nconf*sizeof(int));
  enorm = (double*)malloc(mp->nconf*sizeof(double));
  eorig = (double*)malloc(mp->nconf*sizeof(double));
  efin = (double*)malloc(mp->nconf*sizeof(double));
  eelec = (double*)malloc(mp->nconf*sizeof(double));
  elj = (double*)malloc(mp->nconf*sizeof(double));
  ebond = (double*)malloc(mp->nconf*sizeof(double));
  eangl = (double*)malloc(mp->nconf*sizeof(double));
  edihe = (double*)malloc(mp->nconf*sizeof(double));
  ecmap = (double*)malloc(mp->nconf*sizeof(double));
  enmr = (double*)malloc(mp->nconf*sizeof(double));
  efit = (double*)calloc(mp->nconf, sizeof(double));
  badconf = (char**)malloc(mp->nconf*sizeof(char*));
  badnrg = (double*)malloc(mp->nconf*sizeof(double));
  noutlier = 0;
  for (i = 0; i < mp->nunisys; i++) {
    ncfg = 0;
    for (j = 0; j < mp->nconf; j++) {
      if (mp->conf[j].GroupNum != i) {
        continue;
      }
      if (mp->ao[0] != '\0' && ncfg == 0) {
        fprintf(mscript, "%% System %s\ndata%d = [\n"
                "%% Target  Original    Model     Error\n",
                mp->conf[j].tp->source, i);
      }
      fsys = j;
      fitnrg = 0.0;
      atmp = A->map[j];
      for (k = 0; k < A->col; k++) {
        fitnrg += atmp[k] * x[k];
      }
      confid[ncfg] = j;
      ebond[ncfg] = DSum(mp->conf[j].bmap.Ucontrib, mp->conf[j].bmap.nbond);
      eangl[ncfg] = DSum(mp->conf[j].amap.Ucontrib, mp->conf[j].amap.nangl);
      edihe[ncfg] = DSum(mp->conf[j].hmap.Ucontrib, mp->conf[j].hmap.ntterm);
      ecmap[ncfg] = DSum(mp->conf[j].mmap.Ucontrib, mp->conf[j].mmap.ncmap);
      enmr[ncfg] =  DSum(mp->conf[j].rmap.Ucontrib, mp->conf[j].rmap.nops);
      natom = mp->conf[j].tp->natom;
      eLJ = 0.0;
      eQQ = 0.0;
      for (k = 0; k < natom-1; k++) {
        for (m = k+1; m < natom; m++) {
          eQQ += mp->conf[j].nbnrg.map[k][m] * mp->conf[j].excl.map[k][m];
          eLJ += mp->conf[j].nbnrg.map[m][k] * mp->conf[j].excl.map[m][k];
        }
      }
      eelec[ncfg] = eQQ;
      elj[ncfg] = eLJ;      
      efin[ncfg] = mp->conf[j].nonfitmm + fitnrg;
      enorm[ncfg] = mp->conf[j].enorm;
      eorig[ncfg] = mp->conf[j].eorig;
      efit[i] += fitnrg;
      ncfg++;
    }
    efit[i] /= ncfg;

    // Print the Matlab / Octave readable detailed report
    if (mp->ao[0] != '\0') {
      for (j = 0; j < ncfg; j++) {
        fprintf(mscript, "%9.4lf %9.4lf %9.4lf %9.4lf   %% %s\n", enorm[j],
                eorig[j], efin[j], enorm[j]-efin[j],
                mp->conf[confid[j]].crdsrc);
      }
      fprintf(mscript, "];\n\n");
      fprintf(mscript, "%% System %s\n", mp->conf[confid[0]].tp->source);
      fprintf(mscript, "components%d = [\n", i);
      fprintf(mscript, "%%    elec      LJ        bond      angl      dihe");
      if (mp->ncmaps > 0) {
        fprintf(mscript, "     cmap ");
      }
      if (mp->conf[confid[0]].rmap.nops > 0) {
        fprintf(mscript, "      nmr ");
      }
      fprintf(mscript, "\n");
      for (j = 0; j < ncfg; j++) {
        fprintf(mscript, "%9.4lf %9.4lf %9.4lf %9.4lf %9.4lf ",
                eelec[j], elj[j], ebond[j], eangl[j], edihe[j]);
        if (mp->ncmaps > 0) {
          fprintf(mscript, "%9.4lf ", ecmap[j]);
        }
        if (mp->conf[confid[0]].rmap.nops > 0) {
          fprintf(mscript, "%9.4lf ", enmr[j]);
        }
        fprintf(mscript, "   %% %s\n", mp->conf[confid[j]].crdsrc);
      }
      fprintf(mscript, "];\nfigure; hold on;\n"
              "plot(data%d(:,1), data%d(:,2), 'k.');\n"
              "plot(data%d(:,1), data%d(:,3), 'r.');\n", i, i, i, i);
      fprintf(mscript, "legend('Original', 'Fitted');\n"
              "xlabel('Target energy');\nylabel('Model energy');\n");
      fprintf(mscript, "title('Model results for %s');\n",
              mp->conf[fsys].tp->source);
      fprintf(mscript, "n = input('Ready?');\n\n");
    }

    // Print the standard report
    fitrmse = VecRMSD(efin, enorm, ncfg);
    fprintf(outp, " %-16.16s  %7.4lf  %7.4lf  %7.4lf  %7.4lf  %6d  %7.2lf  "
            "%7.3lf\n", mp->conf[fsys].tp->source, VecRMSD(eorig, enorm, ncfg),
            fitrmse, Pearson(eorig, enorm, ncfg), Pearson(efin, enorm, ncfg),
            nsysvar[i], efit[i], VecRMSD(efin, eorig, ncfg));

    // Record information on outliers that were hard to fit
    fintol = mp->fsigtol*fitrmse;
    if (fintol < mp->fdevfloor) {
      fintol = mp->fdevfloor;
    }
    for (j = 0; j < ncfg; j++) {
      if (fabs(efin[j] - enorm[j]) > fintol) {
        badconf[noutlier] = mp->conf[confid[j]].crdsrc;
        badnrg[noutlier] = efin[j] - enorm[j];
        noutlier++;
      }
    }
  }

  // Print information about the outliers to the standard report
  if (noutlier > 0) {
    fprintf(outp, "\n%% Outliers among all of the fitting data were computed "
            "at %6.2lf sigma.\n%% The data points should be scrutinized, "
            "because they were not fitted\n%% well by the final model.\n"
            "%% Eliminating these data points may yield a better result.\n\n",
            mp->fsigtol);
    fprintf(outp, " Data point coordinates                                  "
            "      E(mm)-E(qm)\n"
                  " --------------------------------------------------------"
            "----  -----------\n");
    for (i = 0; i < noutlier; i++) {
      fprintf(outp, " %-60.60s  %11.4f\n", badconf[i], badnrg[i]);
    }
  }

  // Complete the MatLab output if requested
  if (mp->ao[0] != '\0') {  
    fclose(mscript);
  }

  // Free allocated memory
  free(confid);
  free(enorm);
  free(eorig);
  free(efin);
  free(ebond);
  free(eangl);
  free(edihe);
  free(ecmap);
  free(enmr);
  free(efit);
  free(elj);
  free(eelec);
  free(nsysvar);
  free(badnrg);
  free(badconf);
}

//-----------------------------------------------------------------------------
// SortEpacket: function for quicksort of energy packet contributions.  
//                                                                      
// Arguments:                                                           
//-----------------------------------------------------------------------------
static int SortEpacket(const void *pcA, const void *pcB)
{
  double eA = fabs(((epacket*)pcA)[0].eave);
  double eB = fabs(((epacket*)pcB)[0].eave);

  if (eA < eB) {
    return 1;
  }
  else if (eA > eB) {
    return -1;
  }
  else {
    return 0;
  }
}

//-----------------------------------------------------------------------------
// PrintEnergyContrib: print the energetic contributions of each of the fitted
//                     terms, system-by-system.
//
// Arguments:                                                           
//   mp:       fitting control data (contains atom types and examples of
//             each atom in various topologies)                         
//   A:        the fitting matrix, in its original form before the QR   
//             decomposition                                            
//   x:        the solution vector, containing all fitted parameters    
//   outp:     the output file                                          
//-----------------------------------------------------------------------------
static void PrintEnergyContrib(prmset *mp, dmat *A, double* x, FILE *outp)
{
  int i, j, k, order, parmidx;
  double edr, invN;
  double *dtmp;
  epacket* epc;

  // Array to store the values of each parameter's contribution
  epc = (epacket*)malloc(mp->nparm*sizeof(epacket));

  // Loop over all systems
  for (i = 0; i < mp->nunisys; i++) {

    // Print subsection heading
    fprintf(outp, " System: %s\n", mp->tpencyc[i].source);
    fprintf(outp, "                      Energy    Std. Dev.\n");

    // Sum up the energies
    for (j = 0; j < mp->nparm; j++) {
      epc[j].fitcol = j;
      epc[j].eave = 0.0;
      epc[j].estd = 0.0;
    }
    for (j = 0; j < mp->nconf; j++) {
      if (mp->conf[j].GroupNum == i) {
        dtmp = A->map[j];
        for (k = 0; k < mp->nparm; k++) {
          edr = dtmp[k]*x[k];
          epc[k].eave += edr;
          epc[k].estd += edr*edr;
        }
      }
    }

    // Compute average and standard deviation, then print
    invN = 1.0/mp->GroupCount[i];
    for (j = 0; j < mp->nparm; j++) {
      epc[j].eave *= invN;
      epc[j].estd = sqrt(invN*epc[j].estd - epc[j].eave*epc[j].eave);
    }
    qsort(epc, mp->nparm, sizeof(epacket), SortEpacket);
    for (j = 0; j < mp->nparm; j++) {
      if (fabs(epc[j].eave) > 1.0e-8) {
        FindColumnTerm(mp, epc[j].fitcol, outp, &order, &parmidx, 0);
        if (order >= 0) {
          fprintf(outp, "  %9.4lf  %9.4lf\n", epc[j].eave, epc[j].estd);
        }
      }
    }
  }
  fprintf(outp, "\n\n");

  // Free allocated memory
  free(epc);
}

//-----------------------------------------------------------------------------
// PrintChangeLog: print the atom type change log.                      
//                                                                      
// Arguments:                                                           
//   mp:       fitting control data (contains atom types and examples of
//             each atom in various topologies)                         
//   outp:     the output file                                          
//-----------------------------------------------------------------------------
static void PrintChangeLog(prmset *mp, FILE *outp)
{
  int i;

  fprintf(outp,
          "      System       Atom No Res    Change   Oper'n  "
          "Instances (ambmask format)\n"
          " ----------------  ---- -- ----  --------  ------  "
          "---------------------------\n");
  for (i = 0; i < mp->nchng; i++) {
    fprintf(outp, "%s", mp->ChangeLog.map[i]);
  }
}

//-----------------------------------------------------------------------------
// PrintFrcmodFile: print a file of the parm##.dat or frcmod formats.
//
// Arguments:                                                           
//   mp:       fitting control data (contains atom types and examples of
//             each atom in various topologies)                         
//   x:        solution vector containing new parameters, referenced by 
//             indices in parameter arrays of mp                        
//   ovrw:     flag to indicate that output overwriting is acceptable   
//   fname:    the file name to write                                   
//-----------------------------------------------------------------------------
static void PrintFrcmodFile(prmset *mp, double* x, int ovrw, char* fname)
{
  int i, j, nhbonds;
  FILE *outp;

  // Print the frcmod file for this parameter set
  SetParmOutputFlags(mp);
  outp = FOpenSafe(fname, ovrw);
  fprintf(outp, "%s\n", mp->ititl);
  if (mp->reportall == 0) {
    fprintf(outp, "MASS\n");
  }
  PrintParmAtoms(mp, outp);
  if (mp->Hydrophilics.row > 0 && mp->reportall > 0) {
    for (i = 0; i < mp->Hydrophilics.row; i++) {
      fprintf(outp, "%.4s", mp->Hydrophilics.map[i]);
    }
    fprintf(outp, "\n");
  }
  if (mp->reportall == 0) {
    fprintf(outp, "BOND\n");
  }
  PrintParmBond(mp, 2, x, outp);
  if (mp->reportall == 0) {
    fprintf(outp, "ANGL\n");
  }
  PrintParmBond(mp, 3, x, outp);
  if (mp->reportall == 0) {
    fprintf(outp, "DIHE\n");
  }
  PrintParmBond(mp, 4, x, outp);
  if (mp->reportall == 0) {
    fprintf(outp, "IMPROPER\n");
  }
  PrintParmBond(mp, 5, x, outp);
  nhbonds = 0;
  for (i = 0; i < mp->nhb1012; i++) {
    if (mp->hb1012[i].inreport == 1) {
      nhbonds++;
    }
  }
  if (nhbonds > 0) {
    if (mp->reportall == 0) {
      fprintf(outp, "HBON\n");
    }
    PrintParmBond(mp, 9, x, outp);
  }
  if (mp->reportall == 0) {
    fprintf(outp, "NONB\n");
  }
  else {
    if (mp->neqgroups > 0) {
      for (i = 0; i < mp->neqgroups; i++) {
        for (j = 0; j < mp->eqgroups[i].natom; j++) {
          fprintf(outp, "%.4s", &mp->eqgroups[i].types[4*j]);
        }
        fprintf(outp, "\n");
      }
      fprintf(outp, "\n");
    }
    fprintf(outp, "MOD4      RE\n");
  }
  for (i = 0; i < mp->natom; i++) {
    if (mp->atoms[i].inreport == 1) {
      fprintf(outp, "  %.2s             %10.6lf  %10.6lf  %s\n",
              mp->atoms[i].atype, mp->atoms[i].ljsig, mp->atoms[i].ljeps,
              mp->atoms[i].comment);
    }
  }

  // CMAP parameters can only be reported in a frcmod-style file
  if (mp->reportall == 0) {
    j = 0;
    for (i = 0; i < mp->ncmaps; i++) {
      j += (mp->cmaps[i].fitcoli >= 0);
    }
    if (j > 0) {
      fprintf(outp, "\nCMAP\n");
      PrintParmCMAP(mp, x, outp);
    }
  }
  
  if (mp->reportall > 0) {
    fprintf(outp, "\nEND\n");
  }
  fclose(outp);
}

//-----------------------------------------------------------------------------
// PrintNMROperationsFile: print the NMR operations file.  This will reprint
//                         the input file, with the exception that stiffnesses
//                         of each operation may have different values than
//                         they did initially, even if they were not present in
//                         the input at all.
//
// Arguments:
//   mp:       fitting control data (contains the name of the NMR operations
//             input and output files, as well as the cleaned list of NMR
//             operations)
//   tj:       trajectory control data (for the name of the main input file
//             and file overwriting permissions)
//   x:        the result of the least-squares fitting (needed for stiffness
//             coefficients on each restraint)
//-----------------------------------------------------------------------------
static void PrintNMROperationsFile(prmset *mp, trajcon *tj, double* x)
{
  int i, j, k;
  char* msg;
  FILE *foutp;
  time_t ct;

  // Preamble
  msg = (char*)malloc(2048*sizeof(char));
  ct = time(&ct);
  sprintf(msg, "NMR restraints fitted to reproduce energy surfaces for "
          "systems given in %s on %s.  Restraints were initially laid out in "
          "%s.", tj->inpname, asctime(localtime(&ct)), mp->NMROpsFile);
  foutp = FOpenSafe(mp->NMRParmFile, tj->OverwriteOutput);
  PrintParagraph(msg, 79, "%%", foutp);

  // Loop over all operations
  for (i = 0; i < mp->nops; i++) {
    fprintf(foutp, "&nmropt\n");
    fprintf(foutp, "  Atom1  ");
    for (j = 0; j < mp->nmrops[i].amask.row; j++) {
      fprintf(foutp, "'%s'", mp->nmrops[i].amask.map[j]);
      if (j < mp->nmrops[i].amask.row-1) {
        fprintf(foutp, ", ");
      }
      else {
        fprintf(foutp, "\n");
      }
    }
    fprintf(foutp, "  Atom2  ");
    for (j = 0; j < mp->nmrops[i].bmask.row; j++) {
      fprintf(foutp, "'%s'", mp->nmrops[i].bmask.map[j]);
      if (j < mp->nmrops[i].bmask.row-1) {
        fprintf(foutp, ", ");
      }
      else {
        fprintf(foutp, "\n");
      }
    }
    if (mp->nmrops[i].order >= 3) {
      fprintf(foutp, "  Atom3  ");
      for (j = 0; j < mp->nmrops[i].cmask.row; j++) {
        fprintf(foutp, "'%s'", mp->nmrops[i].cmask.map[j]);
        if (j < mp->nmrops[i].cmask.row-1) {
          fprintf(foutp, ", ");
        }
        else {
          fprintf(foutp, "\n");
        }
      }
    }
    if (mp->nmrops[i].order >= 4) {
      fprintf(foutp, "  Atom4  ");
      for (j = 0; j < mp->nmrops[i].dmask.row; j++) {
        fprintf(foutp, "'%s'", mp->nmrops[i].dmask.map[j]);
        if (j < mp->nmrops[i].dmask.row-1) {
          fprintf(foutp, ", ");
        }
        else {
          fprintf(foutp, "\n");
        }
      }
    }
    fprintf(foutp, "  style %s,  label %s,\n&end\n", mp->nmrops[i].style,
            mp->nmrops[i].label);
    if (i < mp->nops-1) {
      fprintf(foutp, "\n");
    }
  }

  // Close the output file
  fclose(foutp);
}

//-----------------------------------------------------------------------------
// PrintSpectrumReport: print a summary of the spectral resampling.     
//                                                                      
// Arguments:                                                           
//   mp:       fitting control data (contains atom types and examples of
//             each atom in various topologies)                         
//   outp:     the output file                                          
//-----------------------------------------------------------------------------
static void PrintSpectrumReport(prmset *mp, sresamp *ms, trajcon *tj,
                                FILE *outp)
{
  int i, j, k, order, parmidx, ival, nactive, nsamp;
  char tval[64];
  char* sname;
  FILE *tmpfi;
  torterm *tmptor;

  fprintf(outp, " Reduced matrix equation consisted of the following "
          "variables:\n\n");
  fprintf(outp, "   Term  Atom Types  Stiffness   Eq/Phs  Nval  Sample"
          "  Minimum  Maximum   Spc\n"
          "   ---- -- -- -- --  ---------  -------  ----  ------  "
          "-------  -------  -----\n");
  nactive = 0;
  nsamp = 0;
  for (i = 0; i < ms->A->col; i++) {
    if (ms->specvar[i].level == 0 || ms->specvar[i].order == 0) {
      continue;
    }
    nactive++;
    FindColumnTerm(mp, ms->specvar[i].colXmain, outp, &order, &parmidx, 0);
    if (order == 2) {
      fprintf(outp, "  %9.4lf  %7.2lf      ", mp->bonds[parmidx].K, 
              mp->bonds[parmidx].l0);
    }
    else if (order == 3) {
      fprintf(outp, "  %9.4lf  %7.2lf      ", mp->angls[parmidx].K,
              mp->angls[parmidx].th0 * 180.0 / PI);
    }
    else if (order == 4) {
      fprintf(outp, "  %9.4lf  %7.2lf  %4.0lf", mp->torsions[parmidx].K,
              mp->torsions[parmidx].phase, mp->torsions[parmidx].pn);
    }
    if (ms->specvar[i].level == 2) {
      fprintf(outp, "    YES   %7.3lf  %7.3lf  %5.3lf", ms->specvar[i].llimX,
              ms->specvar[i].hlimX, ms->specvar[i].gspcX);
      if (order == 2 || order == 3) {
        fprintf(outp, "\n                                                    "
                "   %7.3lf  %7.3lf  %5.3lf", ms->specvar[i].llimY,
                ms->specvar[i].hlimY, ms->specvar[i].gspcY);
      }
      nsamp++;
    }
    else {
      fprintf(outp, "        ");
    }
    fprintf(outp, "\n");
  }
  fprintf(outp, "   --> Total %d parameters (plus %d system offset columns), "
          "with \n       %d parameters explicitly sampled.\n\n",
          nactive - mp->nunisys, mp->nunisys, nsamp);
  fprintf(outp, " Alternative solutions to the parameter set were found:\n\n");
  fprintf(outp, " ID  Parameter File Name  Additional restraint commands (to "
          "reproduce)\n");
  fprintf(outp, " --  -------------------  ----------------------"
          "-------------------------------\n");
  for (i = 1; i < mp->nspectralx; i++) {

    // The easy stuff: print the index of this solution and the  
    // name of the parameter file, then print the parameter file.
    sname = (char*)malloc(1024*sizeof(char));
    if (strlen(tj->dumpname) > 16) {
      for (j = 0; j < 14; j++) {
        sname[j] = tj->dumpname[j];
      }
      for (j = 14; j < 16; j++) {
        sname[j] = '.';
      }
      sname[16] = '\0';
    }
    else {
      sprintf(sname, "%s", tj->dumpname);
    }
    sprintf(&sname[strlen(sname)], "%-2.2d", i);
    fprintf(outp, " %2d  %19.19s  ", i, sname);
    PrintFrcmodFile(mp, ms->xalt[i].x, tj->OverwriteOutput, sname);

    // The hard stuff: print the restraint commands that
    // reproduce this solution as the optimum if run    
    // without spectrum resampling.                     
    for (j = 0; j < ms->xalt[i].rstpos.row; j++) {
      if (ms->xalt[i].rstpos.map[j][0] < 0) {
        fprintf(outp, "\n");
      }
      FindColumnTerm(mp, ms->xalt[i].rstpos.map[j][1], tmpfi, &order,
                     &parmidx, 1);

      // FIX ME!!!  Need to make this work for bonds and angles.
      if (order == 4) {
        tmptor = &mp->torsions[parmidx];
        fprintf(outp, "shrst ");
        sprintf(tval, "%s", tmptor->atype);
        RemoveWhiteSpace(tval, 1, ' ');
        fprintf(outp, "%s ", tval);
        sprintf(tval, "%s", tmptor->btype);
        RemoveWhiteSpace(tval, 1, ' ');
        fprintf(outp, "%s ", tval);
        sprintf(tval, "%s", tmptor->ctype);
        RemoveWhiteSpace(tval, 1, ' ');
        fprintf(outp, "%s ", tval);
        sprintf(tval, "%s", tmptor->dtype);
        RemoveWhiteSpace(tval, 1, ' ');
        fprintf(outp, "%s per %3.1lf ", tval, tmptor->pn);
        k = ms->xalt[i].rstval.col-1;
        sprintf(tval, "%.7lf", ms->xalt[i].rstval.map[j][k] /
                ms->xalt[i].rstval.map[j][0]);
        RemoveWhiteSpace(tval, 1, '0');
        k = strlen(tval);
        if (tval[k-1] == '.') {
          tval[k] = '0';
          tval[k+1] = '\0';
        }
        fprintf(outp, "rwt %.7lf trg %s\n", ms->xalt[i].rstval.map[j][0],
                tval);
      }
    }
    free(sname);
  }

  // Print the matrices relating each solution.
  fprintf(outp, "\n RMSD in energy estimates (as numbered above) between each "
          "solution.\n From the reduced matrix:\n\n");
  fprintf(outp, "     ");
  for (i = 1; i < mp->nspectralx; i++) {
    fprintf(outp, "  %7d", i);
  }
  fprintf(outp, "\n");
  for (i = 0; i < mp->nspectralx-1; i++) {
    fprintf(outp, "%5d", i);
    for (j = 0; j < i; j++) {
      fprintf(outp, "         ");
    }
    for (j = i+1; j < mp->nspectralx; j++) {
      fprintf(outp, "  %7.4f", ms->ssd1.map[i][j]);
    }
    fprintf(outp, "\n");
  }
  fprintf(outp, "\n After recomputing each solution with the full matrix:"
          "\n\n");
  fprintf(outp, "     ");
  for (i = 1; i < mp->nspectralx; i++) {
    fprintf(outp, "  %7d", i);
  }
  fprintf(outp, "\n");
  for (i = 0; i < mp->nspectralx-1; i++) {
    fprintf(outp, "%5d", i);
    for (j = 0; j < i; j++) {
      fprintf(outp, "         ");
    }
    for (j = i+1; j < mp->nspectralx; j++) {
      fprintf(outp, "  %7.4f", ms->ssd2.map[i][j]);
    }
    fprintf(outp, "\n");
  }
  fprintf(outp, "\n");
}

//-----------------------------------------------------------------------------
// ParamReport: report the best parameters and their fit to the target data.
//              New parameters are reported in a frcmod-format file, while   
//              statistics from the run are reported in a text file.
//
// Arguments:
//   mp:       fitting control data (contains atom types and examples of
//             each atom in various topologies)                         
//   ms:       resampling data (also contains the fitting matrix in its 
//             original form and the global optimum solution)           
//   A:        the fitting matrix, in its original form before the QR   
//             decomposition                                            
//   x:        the solution vector, containing all fitted parameters    
//   tj:       trajectory control information                           
//-----------------------------------------------------------------------------
void ParamReport(prmset *mp, sresamp *ms, trajcon *tj)
{
  int i, j, nsysline;
  char* desc;
  FILE *outp;
  time_t ct;
  cmat warn2, warn3, warn4, warn7, lnwords;

  // Input file header
  ct = time(&ct);
  outp = FOpenSafe(tj->outbase, tj->OverwriteOutput);
  PrintSplash(outp);
  fprintf(outp, "Run on %s", asctime(localtime(&ct)));

  // Reprint the input file
  ReprintInputFile(tj, "sys", "System", mp->PrintFitPoints, outp);

  // Print the energies according to the derived model
  HorizontalRule(outp, 1);
  PrintVADesc(0, "(1.)", 4, " ", 1, "Energies of each system, according to "
              "the fitted parameters.  Units on errors are kcal/mol.  The "
              "original model's energies are compared to target energies "
              "after applying an offset to make the average of each set of "
              "energies for any given system equal (see Section 3).\n", 74, 0,
              outp);
  PrintAccuracy(mp, ms->A, ms->x, tj, outp);
  HorizontalRule(outp, 1);

  // Print the sampling in each fitted parameter, across all systems    
  HorizontalRule(outp, 1);
  PrintVADesc(0, "(2.)", 4, " ", 1, "Parameter sampling.  Bonds and angles "
              "are binned over 18 intervals based on deviations, positive or "
              "negative, from the ideal bond length.  Each interval "
              "represents 0.028 Angstrom or 0.023 radian deviation for bonds "
              "or angles, respectively.  Dihedrals and impropers are binned "
              "at 10-degree intervals, ranging from zero to 180 (the range of "
              "the arccosine function).\n", 74, 0, outp);
  fprintf(outp, " Sampling histogram key:  (zero)   . : = e o U O 0 @ X "
          "( > ten)\n\n");
  warn2 = PrintSamplingTable(mp, 2, ms->x, outp);
  warn3 = PrintSamplingTable(mp, 3, ms->x, outp);
  warn4 = PrintSamplingTable(mp, 4, ms->x, outp);
  if (mp->ncmaps > 0) {
    PrintCMAPSampling(mp, outp);
  }
  HorizontalRule(outp, 1);

  // Print warnings stemming from the sampling of each parameter
  HorizontalRule(outp, 1);
  PrintVADesc(0, "(3.)", 4, " ", 1, "Warnings concerning parameter sampling "
              "and the fitted values of new parameters.", 74, 0, outp);
  fprintf(outp, "\n");
  if (warn2.data[0] != '\0') {
    for (i = 0; i < warn2.row; i++) {
      fprintf(outp, "%s\n", warn2.map[i]);
    }
  }
  if (warn3.data[0] != '\0') {
    for (i = 0; i < warn3.row; i++) {
      fprintf(outp, "%s\n", warn3.map[i]);
    }
  }
  if (warn4.data[0] != '\0') {
    for (i = 0; i < warn4.row; i++) {
      fprintf(outp, "%s\n", warn4.map[i]);
    }
  }
  if (warn2.data[0] == '\0' && warn3.data[0] == '\0' &&
      warn4.data[0] == '\0') {
    fprintf(outp, "  - All parameters pass basic sanity checks.\n\n");
  }
  HorizontalRule(outp, 1);

  // Print energy adjustment constants        
  // (otherwise these will never get noticed).
  HorizontalRule(outp, 1);
  PrintVADesc(0, "(4.)", 4, " ", 1, "Energy adjustment constants.  These "
              "constants are included to bring the overall energy of quantum "
              "mechanical and molecular mechanical calculations into general "
              "agreement so that the fitted parameters can address the energy "
              "differences between multiple states of each system.  Gross "
              "differences between quantum and molecular mechanics are first "
              "eliminated by subtracting off the average energy of all states "
              "for each system; these parameters are then included to \"mop "
              "up.\"\n", 74, 0, outp);
  if (mp->wtfloor <= 100.0) {
    PrintVADesc(0, " ", 4, " ", 1, "Before fitting parameters, the importance "
                "of each system was measured as a function of the target "
                "energy, the lowest energy gaining the most importance and "
                "higher energies receiving less importance down to a minimum "
                "weight specified in the input.  This expanded table includes "
                "the mamimum weights applied to any one conformation of each "
                "system.\n", 74, 0, outp);
  }
  PrintEAConstants(mp, ms->x, outp);
  HorizontalRule(outp, 1);

  // Print an analysis of the fitting matrix,
  // checking for bad conditioning.          
  HorizontalRule(outp, 1);
  PrintVADesc(0, "(5.)", 4, " ", 1, "Fitting matrix analysis.  This is a "
              "check against features of the fitting set that might create a "
              "poorly conditioned matrix.\n", 74, 0, outp);
  PrintMatrixAnalysis(mp, ms->x, outp);
  HorizontalRule(outp, 1);

  // Print the contributions of each parameter
  // to the energy of each system             
  HorizontalRule(outp, 1);
  PrintVADesc(0, "(6.)", 4, " ", 1, "Energy contributions of each fitted "
              "parameter.  All parameters that contribute to each system are "
              "printed, this time with respect to the amount of energy that "
              "they contribute to the molecular mechanical system.\n", 74, 0,
              outp);
  PrintEnergyContrib(mp, ms->A, ms->x, outp);
  HorizontalRule(outp, 1);

  // Print the contexts in which each parameter
  // is sampled; this expounds on section (2.).
  HorizontalRule(outp, 1);
  PrintVADesc(0, "(7.)", 4, " ", 1, "Context of each parameter.  Instances of "
              "each parameter found in the fitting conformations are "
              "enumerated below.\n", 74, 0, outp);
  ParameterDescription(mp, outp);
  HorizontalRule(outp, 1);

  // If changes have been made to any atom types, describe them.
  if (mp->nchng > 0) {
    HorizontalRule(outp, 1);
    PrintVADesc(0, "(8.)", 4, " ", 1, "Changes mades to atom types in each "
                "system's topology file.  These changes were made to expand "
                "or reduce the parameters available for fitting.  The "
                "resulting parameter file must be paired with library files "
                "which reflect all of these changes, and a leap source file "
                "which includes all new atom types.\n", 74, 0, outp);
    PrintChangeLog(mp, outp);
    HorizontalRule(outp, 1);
  }

  // If spectral resampling has been done, report the results.
  if (mp->spectrum == 1) {
    HorizontalRule(outp, 1);
    desc = (char*)malloc(1024*sizeof(char));
    sprintf(desc, "Spectral resampling performed on specific parameters.  "
            "Each parameter in the list below was tightly restrained to "
            "particular values while other parameters were reoptimized in a "
            "reduced matrix equation.  Of those putative solutions that came "
            "within %.2f%% of the accuracy of the global optimum, %d were "
            "selected for evaluation in the context of the full matrix "
            "equation.\n", mp->spvtol, mp->nspectralx);
    PrintVADesc(0, "(9.)", 4, " ", 1, desc, 74, 0, outp);
    free(desc);
    PrintSpectrumReport(mp, ms, tj, outp);
    HorizontalRule(outp, 1);
  }
  fclose(outp);

  // Print the frcmod file for this parameter set
  PrintFrcmodFile(mp, ms->x, tj->OverwriteOutput, tj->dumpname);

  // Print NMR operations with their fitted coefficients
  if (mp->nops > 0) {
    PrintNMROperationsFile(mp, tj, ms->x);
  }
}
