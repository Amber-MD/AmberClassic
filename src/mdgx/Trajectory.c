#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "Parse.h"
#include "mdgxVector.h"
#include "CrdManip.h"
#include "CellManip.h"
#include "Constants.h"
#include "Constraints.h"
#include "pmeRecip.h"
#include "Matrix.h"
#include "Trajectory.h"
#include "Manual.h"
#include "Parse.h"
#include "Thermostats.h"
#include "Integrator.h"
#include "Timings.h"
#include "Topology.h"
#include "Random.h"
#include "AmberNetcdf.h"
#include "Barostats.h"
#include "SpecialMath.h"
#include "MPIMap.h"
#include "BroadcastCommand.h"

#include "pmeDirectDS.h"
#include "CompFrcDS.h"

//-----------------------------------------------------------------------------
// UpdateStepNumber: this routine updates the step number, and manages the
//                   related output file number.                    
//                                                                      
// Arguments:                                                           
//   tj:     trajectory control information                             
//-----------------------------------------------------------------------------
void UpdateStepNumber(trajcon *tj, int newstep)
{
  tj->currstep = newstep;
  tj->currtime = tj->currstep*tj->dt + tj->starttime;

  // If there is no file segment length specified,
  // the step count can be updated and the current
  // file remains set to -1.                      
  if (tj->nfistep == 0) {
    return;
  }

  // Update the file number
  tj->currfi = tj->currstep / tj->nfistep;
}

//-----------------------------------------------------------------------------
// SpliceFileName: splice together an output file name.  If multiple output
//                 files are requested from a single run, the file name will   
//                 reflect the appropriate number in this sequence.
//                                                                      
// Arguments:                                                           
//   tj:     trajectory control information                             
//   base:   the base file name                                         
//   suff:   the file name suffix                                       
//   fname:  on output, the spliced file name                           
//   dprc:   flag to indicate that, if this is the last frame of a file,
//           the step count should be deprecated by 1 in order to write 
//           the result to the proper file                              
//-----------------------------------------------------------------------------
void SpliceFileName(trajcon *tj, char* base, char* suff, char* fname,
                    int dprc)
{
  long long int currstep;

  // Deprecate the step count
  if (dprc == 1) {
    currstep = tj->currstep;
    if (tj->currstep == tj->nstep ||
        (tj->currstep > 0 && tj->nfistep > 0 &&
         tj->currstep % tj->nfistep == 0)) {
      UpdateStepNumber(tj, tj->currstep-1);
    }
  }

  if (tj->nfistep > 0) {
    sprintf(fname, "%s%d%s", base, tj->currfi, suff);
  }
  else {
    sprintf(fname, "%s%s", base, suff);
  }

  // Reset the step count
  if (dprc == 1) {
    UpdateStepNumber(tj, currstep);
  }
}

//-----------------------------------------------------------------------------
// InitializeEnergy: initialize the energy accumulators and set flags for  
//                   calculating energies.                          
//                                                                      
// Arguments:                                                           
//   sysUV:  the energy / virial accumulator                            
//   tj:     the trajectory control data                                
//-----------------------------------------------------------------------------
void InitializeEnergy(Energy *sysUV, trajcon *tj, prmtop *tp,
                      int allocBondUdc)
{
  // Initialize energy accumulators
  sysUV->Esummed = 0;
  sysUV->delec = 0.0;
  sysUV->relec = 0.0;
  sysUV->vdw12 = 0.0;
  sysUV->vdw6 = 0.0;
  sysUV->bond = 0.0;
  sysUV->angl = 0.0;
  sysUV->dihe = 0.0;
  sysUV->cmap = 0.0;
  SetDVec(sysUV->Vir, 9, 0.0);
  if (allocBondUdc == 1) {
    sysUV->nUdc = 3*tp->nBAH.nbond;
    sysUV->BondUdc = (double*)malloc(3*tp->nBAH.nbond*sizeof(double));
  }
  SetDVec(sysUV->BondUdc, 3*tp->nBAH.nbond, 0.0);
  sysUV->dVdL = 0.0;

  // Decide whether energies or virials need to be computed
  if (tj->ntp > 0 && tj->barostat == 1) {
    sysUV->updateU = 2;
  }
  else if ((tj->ntp > 0 && tj->barostat == 2 &&
            tj->currstep % tj->MCBarostatFreq == 0) ||
           (tj->ntpr > 0 && tj->currstep % tj->ntpr == 0) || tj->TI == 1 ||
           tj->mode == 1) {
    sysUV->updateU = 1;
  }
  else {
    sysUV->updateU = 0;
  }
}

//-----------------------------------------------------------------------------
// DestroyEnergyTracker: de-allocate memory for an energy tracking struct.
//                                                                      
// Arguments:                                                           
//   sysUV:  the energy / virial accumulator                            
//-----------------------------------------------------------------------------
void DestroyEnergyTracker(Energy *sysUV)
{
  free(sysUV->BondUdc);
}

//-----------------------------------------------------------------------------
// ExtendCoordinates: this function is called in order to bring a set of
//                    coordinates read from a file up to the size of   
//                    coordinates expected by the existing topology.  It exists
//                    because mdgx is able to extend topologies at run time.
//                                                                      
// Arguments:                                                           
//   tc:     the coordinates                                            
//   tp:     the topology                                               
//-----------------------------------------------------------------------------
void ExtendCoordinates(coord *tc, prmtop *tp)
{
  int i, j, oldnum;

  // Check to make sure that the number of atoms in the
  // coordinates is what the topology started out with 
  if (tc->natom != tp->norigatom) {
    printf("ExtendCoordinates >> Error.  Atom count in the coordinates read "
           "from file\nExtendCoordinates >> is %d, topology originally "
           "contained %d.\n", tc->natom, tp->norigatom);
    exit(1);
  }

  // Reallocate arrays
  free(tc->atmid);
  tc->atmid = CountUp(tp->natom);
  tc->loc = (double*)realloc(tc->loc, 3*tp->natom*sizeof(double));
  tc->prvloc = (double*)realloc(tc->prvloc, 3*tp->natom*sizeof(double));
  tc->scrloc = (double*)realloc(tc->scrloc, 3*tp->natom*sizeof(double));
  tc->vel = (double*)realloc(tc->vel, 3*tp->natom*sizeof(double));
  tc->prvvel = (double*)realloc(tc->prvvel, 3*tp->natom*sizeof(double));
  tc->frc = (double*)realloc(tc->frc, 3*tp->natom*sizeof(double));
  tc->prvfrc = (double*)realloc(tc->prvfrc, 3*tp->natom*sizeof(double));
  tc->scrfrc = (double*)realloc(tc->scrfrc, 3*tp->natom*sizeof(double));

  // Spread the atoms into their new registers
  for (i = tp->natom-1; i >= 0; i--) {
    oldnum = tp->OldAtomNum[i];
    if (oldnum > -1) {
      for (j = 0; j < 3; j++) {
        tc->loc[3*i+j] = tc->loc[3*oldnum+j];
        tc->vel[3*i+j] = tc->vel[3*oldnum+j];
      }
    }
    else {
      for (j = 0; j < 3; j++) {
        tc->loc[3*i+j] = 0.0;
        tc->vel[3*i+j] = 0.0;
      }
    }
  }

  // Update the number of atoms
  tc->natom = tp->natom;
}

//-----------------------------------------------------------------------------
// ValidCoord: look at the first few and (in the case of an ascii file) the
//             last few lines of a file to determine whether it is a valid
//             inpcrd format file, restart file, trajectory in the old .crd
//             format, NetCDF restart file, NetCDF coordinates file, or an 
//             IPolQ energy and structure file.  Return 1 for ascii inpcrd, 2
//             for ascii restart, 3 for .crd format trajectory, 4 for NetCDF
//             restart, 5 for a NetCDF trajectory, and 6 for an IPolQ energy
//             and structure file.  Return 0 if no format can be identified.
//
// Arguments:
//   fname:    the name of the file 
//-----------------------------------------------------------------------------
int ValidCoord(char* fname)
{
  int i, j, k, nlnexpect, isascii, testchar, slen, natom, nword;
  long long int il, fsize;
  char nbuff[16];
  FILE *finp;
  cmat contents, LW, LW2;
  cdftrj Acdf;

  // Ascii files could fall into one of three categories.
  if (DetectBinaryFile(fname) == 0) {
    contents = Ascii2Mem(fname, 96, 8, "Coordinate file not found.\n");
    if (contents.row < 2) {

      // Free allocated memory and return
      DestroyCmat(&contents);
      return 0;
    }
    LW = ParseWords(contents.map[1]);
    if (LW.row < 1) {

      // Free allocated memory and return
      DestroyCmat(&LW);
      DestroyCmat(&contents);
      return 0;
    }
    if (WordIsInteger(LW.map[0]) == 1) {
      natom = atoi(LW.map[0]);
      DestroyCmat(&LW);

      // Could be a restart file or inpcrd.  The third line will be used
      // to check that the format is all right, and the total number of
      // lines will determine whether
      slen = strlen(contents.map[2]);
      nword = CountWords(contents.map[2]);      
      if (nword <= 6) {
        for (i = 0; i < slen; i+=12) {
          k = 0;
          for (j = i; (j < i+12 && j < slen); j++) {
            if (contents.map[2][j] != ' ') {
              nbuff[k] = contents.map[2][j];
              k++;
            }
          }
          nbuff[k] = '\0';
          if (WordIsNumber(nbuff) == 0) {
            for (k = 0; k < j-i; k++) {
              if (nbuff[k] != ' ' && nbuff[k] != '\n') {

                // Free allocated memory and return
                DestroyCmat(&contents);
                return 0;
              }
            }
          }
        }

        // The number of lines expected in an inpcrd file is two plus
        // however many are needed to give coordinates of natom atoms
        // at two sets per line, and possibly a final line for the box
        // dimensions.
        nlnexpect = 2 + (natom/2) + (natom % 2);
        if (contents.row == nlnexpect || contents.row == nlnexpect + 1) {

          // Free allocated memory and return
          DestroyCmat(&contents);
          return 1;
        }

        // For restart files, an set of velocities (in addition to
        // the header and tail) is expected.
        nlnexpect = 2 + 2*((natom/2) + (natom % 2));
        if (contents.row == nlnexpect || contents.row == nlnexpect + 1) {

          // Free allocated memory and return
          DestroyCmat(&contents);
          return 2;
        }

        // Otherwise this file format is not recognized.
        DestroyCmat(&contents);
        return 0;
      }
      else {

        // Free allocated memory and return
        DestroyCmat(&contents);
        return 0;
      }
    }
    else {
      DestroyCmat(&LW);


      // Try to detect the IPolQ energy and coordinates format.  This
      // format contains data extracted from QM/MM calculations
      // driven by the implicitly polarized charge framework.
      nword = CountWords(contents.map[0]);
      if (nword >= 5) {
        LW = ParseWords(contents.map[0]);
        if (strcmp(LW.map[0], "%") == 0 && strcmp(LW.map[2], "by") == 0 &&
            strcmp(LW.map[3], "mdgx") == 0 &&
            strcmp(LW.map[1], "Generated") == 0 &&
            strncmp(LW.map[4], "potential", 9) == 0) {

          // So far, so good.  The 5th line should be coordinates.
          nword = CountWords(contents.map[4]);
          if (nword == 3) {
            LW2 = ParseWords(contents.map[4]);
            if (strcmp(LW2.map[0], "coords") == 0 &&
                strcmp(LW2.map[1], "=") == 0 && strcmp(LW2.map[2], "[") == 0) {

              // The IPolQ energy and coordinates format has been detected
              DestroyCmat(&LW2);
              DestroyCmat(&LW);
              
              return 6;
            }
            DestroyCmat(&LW2);
          }
        }
        DestroyCmat(&LW);
      }

      
      // Could be an old .crd format trajectory.
      // Check the third line to be sure.
      slen = strlen(contents.map[2]);
      nword = CountWords(contents.map[2]);
      if (nword <= 10) {

        for (i = 0; i < slen; i+=8) {
          k = 0;
          for (j = i; (j < i+8 && j < slen); j++) {
            if (contents.map[2][j] != ' ') {
              nbuff[k] = contents.map[2][j];
              k++;
            }
          }
          nbuff[k] = '\0';
          if (WordIsNumber(nbuff) == 0) {
            for (k = 0; k < j-i; k++) {
              if (nbuff[k] != ' ' && nbuff[k] != '\n') {

                // Free allocated memory and return
                DestroyCmat(&contents);
                return 0;
              }
            }
          }
        }

        // Free allocated memory and return
        DestroyCmat(&contents);
        return 3;
      }
      else {

        // Free allocated memory and return
        DestroyCmat(&contents);
        return 0;
      }
    }
  }

  // If we are still here, the file must be binary.
  if (netcdfLoad(&Acdf, fname) != 0) {
    return 0;
  }
  if (Acdf.isNCrestart == 1) { 
    netcdfClose(&Acdf);

    // This is a binary restart file
    return 4;
  }
  netcdfClose(&Acdf);

  // This is a binary trajectory file
  return 5;
}

//-----------------------------------------------------------------------------
// ReadRst: reads an AMBER restart or coordinates input file.
//                                                                      
// Arguments:                                                           
//   tp:     the topology                                               
//   source: the name of the input coordinates source file              
//-----------------------------------------------------------------------------
coord ReadRst(prmtop *tp, char* source, double *starttime)
{
  int i, iorig, rsttype, nrow, ncol, hasbox, wcount;
  char *ctmp;
  cmat Crst, lwords;
  coord tc;
  cdftrj Acdf;

  // Binary file?  No problem.
  if (DetectBinaryFile(source) == 1) {
    if (netcdfLoad(&Acdf, source) != 0 || Acdf.isNCrestart == 0) {
      printf("ReadRst >> Attempted to read %s as a NetCDF\nReadRst >> restart "
             "file, but failed.\n", source);
      exit(1);
    }
    tc = CreateCoord(Acdf.ncatom);
    if (netcdfGetFrame(&Acdf, 0, tc.loc, tc.gdim) != 0 ||
        netcdfGetVelocity(&Acdf, 0, tc.vel) != 0) {
      printf("ReadRst >> Attempted to read %s as a NetCDF\nReadRst >> restart "
             "file, but failed.\n", source);
      exit(1);
    }
    netcdfClose(&Acdf);

    // Add details that would have existed for the ascii-format case
    hasbox = (fabs(tc.gdim[0]) > 1.0e-4 && fabs(tc.gdim[1]) > 1.0e-4 &&
	      fabs(tc.gdim[0]) > 1.0e-4);
    if (hasbox) {
      for (i = 3; i < 6; i++) {
        tc.gdim[i] *= PI/180.0;
        tc.hgdim[i-3] = 0.5*tc.gdim[i-3];
      }
      CompXfrm(tc.gdim, tc.U, tc.invU);
      tc.isortho = TestUnitCellOrtho(tc.gdim);
    }
    else {
      tc.isortho = 0;
    }
  }
  else {

    // This is an ascii character file.  Read it into memory.
    Crst = Ascii2Mem(source, 96, 8, "Missing inpcrd or restrt file.\n");

    // Check whether this is a restart file or an input coordinates file
    i = tp->natom / 2 + tp->natom % 2;
    iorig = (tp->EPInserted == 1) ? tp->norigatom / 2 + tp->norigatom % 2 : -1;
    if (Crst.row == i+3 || Crst.row == iorig+3 ||
        Crst.row == i+2 || Crst.row == iorig+2) {
      hasbox = ((tp->EPInserted == 0 && Crst.row == i+2) ||
                (tp->EPInserted == 1 && Crst.row == iorig+2)) ? 0 : 1;
      rsttype = 1;
    }
    else if (Crst.row == 2*i+3 || Crst.row == 2*iorig+3 ||
             Crst.row == 2*i+2 || Crst.row == 2*iorig+2) {
      hasbox = ((tp->EPInserted == 0 && Crst.row == 2*i+2) ||
                (tp->EPInserted == 1 && Crst.row == 2*iorig+2)) ? 0 : 1;
      rsttype = 2;
    }
    else {
      printf("ReadRst >> Error.  Topology %s has %d atoms\nReadRst >> but "
             "restart file %s has %d lines.\nReadRst >> This cannot be a valid "
             "coordinate/velocity set.\n", tp->source, tp->natom,
             source, Crst.row);
      exit(1);
    }

    // Check the atom count
    sscanf(Crst.map[1], "%d", &i);
    if (i != tp->natom && tp->EPInserted == 0) {
      printf("ReadRst >> Error.  Coordinate file %s contains %d atoms, "
	     "topology\nReadRst >> %s contains %d.\n", source, i, tp->source,
	     tp->natom);
      exit(1);
    }
    else if (i > tp->natom && tp->EPInserted == 0) {
      printf("ReadRst >> Error.  Extra points were inserted into this "
	     "topology, but there are\nReadRst >> as many or more coordinates "
	     "than total points in the topology.\n");
      exit(1);
    }

    // Get the starting time
    if (*starttime < 0.0) {
      wcount = CountWords(Crst.map[1]);
      if (wcount == 2) {
        lwords = ParseWords(Crst.map[1]);
        if (WordIsNumber(lwords.map[1]) == 1) {
          sscanf(lwords.map[1], "%lf", starttime);
        }
        else {
          printf("ReadRst >> Warning.  Invalid format for time stamp (%s) in\n"
                 "ReadRst >> file %s.\n", lwords.map[1], source);
          *starttime = 0.0;
        }
        DestroyCmat(&lwords);
      }
      else {
        *starttime = 0.0;
      }
    }
  
    // Read coordinates
    tc = CreateCoord(i);
    nrow = 2;
    ncol = 0;
    for (i = 0; i < 3*tc.natom; i++) {
      tc.loc[i] = RealXpYf(&Crst.map[nrow][ncol], 12, 7);
      ncol += 12;
      if (ncol == 72 || i == 3*tc.natom-1) {
        ncol = 0;
        nrow++;
      }
    }
    if (rsttype == 2) {
      for (i = 0; i < 3*tc.natom; i++) {
        tc.vel[i] = RealXpYf(&Crst.map[nrow][ncol], 12, 7);
        ncol += 12;
        if (ncol == 72 || i == 3*tc.natom-1) {
          ncol = 0;
          nrow++;
        }
      }
    }
    else {
      SetDVec(tc.vel, 3*tc.natom, 0.0);
    }

    // Read the box dimensions and compute transformation matrices
    if (hasbox) {
      for (i = 0; i < 6; i++) {
        tc.gdim[i] = RealXpYf(&Crst.map[nrow][ncol], 12, 7);
        ncol += 12;
      }
      for (i = 3; i < 6; i++) {
        tc.gdim[i] *= PI/180.0;
        tc.hgdim[i-3] = 0.5*tc.gdim[i-3];
      }
      CompXfrm(tc.gdim, tc.U, tc.invU);
      tc.isortho = TestUnitCellOrtho(tc.gdim);
    }
    else {
      for (i = 0; i < 6; i++) {
        tc.gdim[i] = 0.0;
      }
      tc.isortho = 0;
    }

    // Free allocated memory
    DestroyCmat(&Crst);
  }

  // Initialize atom identification numbers
  for (i = 0; i < tc.natom; i++) {
    tc.atmid[i] = i;
  }

  // If extra points were inserted into the topology, we may have to
  // adjust the coordinates in order to accommodate the insertions  
  if (tp->EPInserted == 1 && tc.natom < tp->natom) {
    ExtendCoordinates(&tc, tp);
  }

  return tc;
}

//-----------------------------------------------------------------------------
// ReadCrdTraj: read a series of frames from a trajectory file in .crd format.
//              Check for the existence of box dimensions with each frame.
//              Return the results as a double-precision matrix, with each
//              frame being a column.  If there is box information, then this
//              will be included in the final six rows.
//
// Arguments:
//   tp:       topology that describes the system
//   source:   name of the trajectory file to read
//   readbox:  flag to have box information stored alongside coordinates
//-----------------------------------------------------------------------------
dmat ReadCrdTraj(prmtop *tp, char* source, int readbox)
{
  int i, j, jln, jpos, nfrm, ndim, nval, maxval;
  int lines_per_frame, last_crd_line, last_crd_pos;
  char line[256], crdbuff[9];
  char *ctmp;
  FILE *finp;
  dmat allcrd;
  cmat traj;

  // Open the file
  if ((finp = fopen(source, "r")) == NULL) {
    printf("ReadCrd >> Error.  Trajectory file %s not found.\n", source);
    exit(1);
  }

  // Read the trajectory verbatim into memory
  ndim = (tp->ifbox > 0 && readbox == 1) ? 3*tp->natom + 6 : 3*tp->natom;
  i = 0;
  while (fgets(line, 96, finp) != NULL) {
    i++;
  }
  fclose(finp);
  traj = Ascii2Mem(source, 96, i,
                   "Attempted to parse this as a .crd format trajectory.");

  // The first frame will serve as a check that
  // this trajectory looks like the topology.
  nval = 0;
  maxval = 3*tp->natom;
  lines_per_frame = maxval/10;
  if (maxval % 10 > 0) {
    lines_per_frame += 1;
    last_crd_pos = (maxval % 10 - 1)*8;
  }
  else {
    last_crd_pos = 72;
  }
  last_crd_line = lines_per_frame - 1;
  if (tp->ifbox > 0) {
    lines_per_frame += 1;
  }

  // Count the number of actual frames by checking the last
  // expected coordinate in each of them
  nfrm = 0;
  crdbuff[8] = '\0';
  for (i = 1; i < traj.row; i += lines_per_frame) {
    ctmp = traj.map[i+last_crd_line];
    for (j = 0; j < 8; j++) {
      crdbuff[j] = ctmp[j+last_crd_pos];
    }
    if (WordIsNumber(crdbuff) == 1) {
      nfrm++;
    }
    else {
      break;
    }
  }

  // Read in the frames
  allcrd = CreateDmat(ndim, nfrm, 0);
  for (i = 0; i < nfrm; i++) {
    jln = 1 + i*lines_per_frame;
    jpos = 0;
    for (j = 0; j < maxval; j++) {
      allcrd.map[j][i] = RealXpYf(&traj.map[jln][jpos], 8, 3);
      jpos += 8;
      if (jpos == 80) {
        jpos = 0;
        jln++;
      }
    }
    if (tp->ifbox > 0 && readbox == 1) {
      allcrd.map[maxval][i] = RealXpYf(&traj.map[jln][0], 8, 3);
      allcrd.map[maxval+1][i] = RealXpYf(&traj.map[jln][8], 8, 3);
      allcrd.map[maxval+2][i] = RealXpYf(&traj.map[jln][16], 8, 3);

      // Box angles are not present in th eold .crd format files.
      allcrd.map[maxval+3][i] = 0.0;
      allcrd.map[maxval+4][i] = 0.0;
      allcrd.map[maxval+5][i] = 0.0;
    }
    nfrm++;
  }

  // Free allocated memory
  DestroyCmat(&traj);

  return allcrd;
}

//-----------------------------------------------------------------------------
// ReadCDFTraj: read a trajectory in NetCDF format.  Like ReadCrd above, this
//              returns the trajectory as a double-precision matrix, with each
//              frame occupying a column and the box dimensions, if applicable,
//              in the last six rows.
//
// Arguments:
//   tp:       topology that describes the system
//   source:   name of the trajectory file to read
//   readbox:  flag to have box information stored alongside coordinates
//-----------------------------------------------------------------------------
dmat ReadCDFTraj(prmtop *tp, char* source, int readbox)
{
  int i, j, ndim, natom3;
  double boxdim[6];
  double* crdbuff;
  cdftrj Acdf;
  dmat allcrd;

  if (netcdfLoad(&Acdf, source) != 0) {
    printf("ReadCDFTraj >> Error.  File %s could not be read.\n", source);
    exit(1);
  }
  natom3 = 3*tp->natom;
  ndim = (tp->ifbox > 0 && readbox == 1) ? natom3 + 6 : natom3;
  allcrd = CreateDmat(ndim, Acdf.ncframe, 0);
  crdbuff = (double*)malloc(natom3*sizeof(double));
  for (i = 0; i < Acdf.ncframe; i++) {
    if (netcdfGetFrame(&Acdf, i, crdbuff, boxdim) != 0) {
      printf("ReadCDFTraj >> Error.  File %s has a corrupted frame.\n",
             source);
      exit(1);
    }
    for (j = 0; j < natom3; j++) {
      allcrd.map[j][i] = crdbuff[j];
    }
    if (tp->ifbox > 0 && readbox == 1) {
      for (j = 0; j < 6; j++) {
        allcrd.map[natom3+j][i] = boxdim[j];
      }
    }
  }
  netcdfClose(&Acdf);

  return allcrd;
}

//-----------------------------------------------------------------------------
// InitCoords: initialize coordinates, first trying to open a specified 
//             starting coordinates file and then proceeding to test for
//             the existence of successive (complete) output files, finally
//             beginning the dynamics at the first non-existent or unfinished
//             segment.  This routine also initializes the step number.
//                                                                      
// Arguments:                                                           
//   tp:     the topology                                               
//   tj:     trajectory control information                             
//   n:      the number of the coordinates file to read (coordinates for
//           the nth system)                                            
//-----------------------------------------------------------------------------
coord InitCoords(prmtop *tp, trajcon *tj, int n)
{
  int i, fcmplt, fndstrt, coherent, tid;
  char fname[MAXNAME], line[128];
  FILE *tst;
  coord tc;

  // If this is the master process of the cell grid's
  // communicator, read the restart file from disk.  
  // Otherwise, wait for the master to deliver the   
  // coordinates over MPI.                           
  UpdateStepNumber(tj, 0);
#ifdef MPI
  MPI_Comm_rank(tj->SysComm[n], &tid);
#else
  tid = 0;
#endif
  if (tid == 0) {

    // First try reading the starting coordinates file
    tc = ReadRst(tp, tj->ipcname.map[n], &tj->starttime);

    // Now, start looking for segments of a trajectory
    if (tj->nfistep > 0) {
      fndstrt = 0;
      while (fndstrt == 0) {
        SpliceFileName(tj, tj->outbase, tj->outsuff, fname, 0);

        // Check for the existence, and then the completeness, 
        // of the output diagnostics file.                     
        fcmplt = 0;
        if ((tst = fopen(fname, "r")) != NULL) {
          while(fcmplt == 0 && fgets(line, 128, tst) != NULL) {
            if (line[0] == '$' &&
                strncmp(line, "$! Closing watermark.  DO NOT alter this "
                        "segment of the file. !$", 64) == 0) {
              fcmplt = 1;
            }
          }
          fclose(tst);
        }
        if (fcmplt == 0) {
          fndstrt = 1;
          continue;
        }

        // We check for existence of the restart file here, because  
        // lack of a restart file just means that the trajectory     
        // segment, for whatever reason, was not finished.  It should
        // not prompt the program to abort.  The existence of restart
        // files for all alternate systems is checked here, to ensure
        // that the simulation starts from a coherent point.         
        coherent = 1;
        for (i = 0; i < tj->nsys; i++) {
          SpliceFileName(tj, tj->rstbase.map[i], tj->rstsuff.map[i], fname, 0);
          if ((tst = fopen(fname, "r")) != NULL) {
            fclose(tst);
          }
          else {
            coherent = 0;
            break;
          }
        }
        if (coherent == 1) {
          UpdateStepNumber(tj, tj->currstep+tj->nfistep);
        }
        else {
          fndstrt = 1;
        }
      }

      // If the trajectory has advanced, 
      // read in the latest restart file.
      if (tj->currstep > 0) {

        // Treat this as any other restart
        tj->irest = 1;
        DestroyCoord(&tc);
        SpliceFileName(tj, tj->rstbase.map[n], tj->rstsuff.map[n], fname, 1);
        tc = ReadRst(tp, fname, &tj->starttime);
      }
    }
  }
  else {

    // Allocate memory for a blank set of coordinates
    tc = CreateCoord(tp->natom);
  }

  // If there is more than one thread,
  // broadcast the initial coordinates
#ifdef MPI
  if (tj->SystemCPUs.col > 1) {
    BroadcastCoordinates(&tc, tj, n);
  }
#endif

  // We're done if the final restart segment has been written
  if (tj->currstep == tj->nstep) {
    exit(1);
  }

  // If the trajectory has advanced, we must adjust the random seeds
  if (tj->currstep > 0) {
    tj->igseed = -1;
    tj->rndcon = SetPRNG(tj);
  }

  // Trap against virial computations in non-orthorhombic cells
  if (tj->ntp > 0 && tj->barostat == 1 && tc.isortho == 0) {
    printf("InitCoords >> Virial computations in non-orthorhombic unit cells "
           "are not\nInitCoords >> yet enabled.  Use Monte-Carlo Barostat "
           "instead.\n");
    exit(1);
  }

  return tc;
}

//-----------------------------------------------------------------------------
// SummarizeOutputFiles: print the restart, and all trajectory files to 
//                       the output diagnostics file header.            
//                                                                      
// Arguments:                                                           
//   tj:     trajectory control data                                    
//   c:      used in cases of thermodynamic integration to signify "I"  
//           initial or "F" final states, but white space otherwise     
//   n:      the system number (0 when there is only one system, but    
//           0 or 1 in the case of thermodynamic integration and up to  
//           MAXSYS in replica exchange)                                
//-----------------------------------------------------------------------------
static void SummarizeOutputFiles(trajcon *tj, char c, int n, FILE *outp)
{
  if (tj->ntwr > 0) {
    fprintf(outp, " Restart file %c      %-45s   %-8s\n", c,
            tj->rstbase.map[n], tj->rstsuff.map[n]);
  }
  if (tj->ntwx > 0) {
    fprintf(outp, " Coordinate Traj. %c  %-45s   %-8s\n", c,
            tj->trjbase.map[n], tj->trjsuff.map[n]);
  }
  if (tj->ntwv > 0) {
    fprintf(outp, " Velocity Traj. %c    %-45s   %-8s\n", c,
            tj->velbase.map[n], tj->velsuff.map[n]);
  }
  if (tj->ntwf > 0) {
    fprintf(outp, " Force Traj. %c       %-45s   %-8s\n", c,
            tj->frcbase.map[n], tj->frcsuff.map[n]);
  }
}

//-----------------------------------------------------------------------------
// PrintSystemSummary: print a summary of the system to the output diagnostics
//                     file header.                         
//                                                                      
// Arguments:                                                           
//   tp:     the topology                                               
//   outp:   the output file pointer                                    
//-----------------------------------------------------------------------------
static void PrintSystemSummary(prmtop *tp, FILE *outp)
{
  fprintf(outp, " Number of atoms:        %10d     Number of residues:      "
          "%10d\n Number of rigid waters: %10d     Number of extra points:  "
          "%10d\n\n", tp->natom, tp->nres, tp->nwat, tp->nxtrapt);
  fprintf(outp, " Bonds with hydrogen:    %10d     Bonds without hydrogen:  "
          "%10d\n Number of bond types:   %10d\n\n", tp->withH.nbond,
          tp->woH.nbond, tp->nBAH.nbond);
  fprintf(outp, " Number of angles:       %10d     Number of angle types:   "
          "%10d\n", tp->withH.nangl + tp->woH.nangl, tp->nBAH.nangl);
  fprintf(outp, " Number of dihedrals:    %10d     Number of dihedral types:"
          "%10d\n\n", tp->withH.ndihe + tp->woH.ndihe, tp->nBAH.ndihe);
  fprintf(outp, " Total system charge:    %10.6lf     Number of vdW types:    "
          " %10d\n\n", DSum(tp->Charges, tp->natom), tp->LJutab.row);
  if (tp->initq - DSum(tp->Charges, tp->natom) > 1.0e-8) {
    fprintf(outp, " Note: total system charge was rounded from its initial "
            "value %10.6lf\n by applying a correction charge equally over all "
            "atoms.\n\n", tp->initq);
  }
}

//-----------------------------------------------------------------------------
// GetMostPopulatedCell: this function serves only OpenDiagnosticsFile() at the
//                       moment, and so is in this library.  It scans all cells
//                       in a cell grid and records the largest or smallest
//                       number of atoms.           
//                                                                      
// Arguments:                                                           
//   CG:       the cell grid                                            
//   extnum:   the extremum to seek, -1 for minimum and +1 for maximum  
//   order:    single atoms or pairs
//-----------------------------------------------------------------------------
static int GetMostPopulatedCell(cellgrid *CG, int extnum, int order)
{
  int i, ccon, mcount, mbuff;
#ifdef MPI
  MPI_Op myop;
#endif

  if (order == 1) {
    mcount = (extnum == -1) ? CG->maxatom : 0;
  }
  else if (order == 2) {
    mcount = (extnum == -1) ? 1000000000 : 0;
  }
  for (i = 0; i < CG->MyCellCount; i++) {
    if (order == 1) {
      ccon = CG->data[CG->MyCellDomain[i]].nr[0];
    }
    else if (order == 2) {
      ccon = CalcCellR2Cost(&CG->data[CG->MyCellDomain[i]], CG);
    }
    if ((extnum == -1 && ccon < mcount) || (extnum == 1 && ccon > mcount)) {
      mcount = ccon;
    }
  }
#ifdef MPI
  myop = (extnum == -1) ? MPI_MIN : MPI_MAX;
  MPI_Reduce(&mcount, &mbuff, 1, MPI_INT, myop, 0, CG->dspcomm);
#else
  mbuff = mcount;
#endif

  return mbuff;
}

//-----------------------------------------------------------------------------
// OpenDiagnosticsFile: open a new diagnostics file and print its preamble.
//                      The result is very similar to typical AMBER output,
//                      with a few differences where they make the output
//                      clearer or more meaningful.     
//                                                                      
// Arguments:                                                           
//   tj:     trajectory control information                             
//   tp:     the topology                                               
//   dcinp:  direct space control data                                  
//   rcinp:  reciprocal space control data                              
//   Etab:   the electrostatic direct space spline interpolation table  
//   crd:    the coordinates                                            
//   CG:     the cell grid                                              
//   n:      the system number                                          
//-----------------------------------------------------------------------------
void OpenDiagnosticsFile(trajcon *tj, prmtop *tp, dircon *dcinp, reccon *rcinp,
			 FrcTab *Etab, coord *crd, cellgrid *CG, int n)
{
  int i, iwarning, minCatom, maxCatom, minCpair, maxCpair;
  char fname[MAXNAME];
  FILE *outp;
  time_t ct;

  // Get collective information about the
  // cell grid for performance reporting
  if (tp->ifbox > 0 && CG != NULL) {
    maxCatom = GetMostPopulatedCell(CG, 1, 1);
    minCatom = GetMostPopulatedCell(CG, -1, 1);
    maxCpair = GetMostPopulatedCell(CG, 1, 2);
    minCpair = GetMostPopulatedCell(CG, -1, 2);
  }  
#ifdef MPI
  // Bail right out if this is not the master process
  if (tp->ifbox > 0 && CG != NULL) {
    if (tj->topchk == 1) {
      CoordinateReduction(CG, crd, tp);
    }
    if (CG->tid != 0) {
      return;
    }
  }
  else {
    MPI_Comm_rank(MPI_COMM_WORLD, &i);
    if (i != 0) {
      return;
    }
  }
#endif
  ct = time(&ct);
  SpliceFileName(tj, tj->outbase, tj->outsuff, fname, 0);
  outp = FOpenSafe(fname, tj->OverwriteOutput);
  PrintSplash(outp);
  fprintf(outp, "Run on %s", asctime(localtime(&ct)));

  // Print file assignments
  fprintf(outp, "\nINPUT FILE ASSIGNMENTS:\n"
          " Designation         Base name                                     "
          "  Suffix\n"
          " -----------------   --------------------------------------------- "
          "  ----------\n");
  fprintf(outp, " Input File          %-45s   NONE\n", tj->inpname);
  if (tj->TI != 1) {
    fprintf(outp, " Input Coordinates   %-45s   NONE\n", tj->ipcname.map[n]);
    fprintf(outp, " Topology            %-45s   NONE\n", tp->source);
  }
  else {
    fprintf(outp, " Input Coordinates I %-45s   NONE\n", tj->ipcname.map[0]);
    fprintf(outp, " Topology I          %-45s   NONE\n", tp[0].source);
    fprintf(outp, " Input Coordinates F %-45s   NONE\n", tj->ipcname.map[1]);
    fprintf(outp, " Topology F          %-45s   NONE\n", tp[1].source);
  }
  fprintf(outp, "\n");
  fprintf(outp, "OUTPUT FILE ASSIGNMENTS:\n"
          " Designation         Base name                                     "
          "  Suffix\n"
          " -----------------   --------------------------------------------- "
          "  ----------\n");
  fprintf(outp, " State Information   %-45s   %-8s\n", tj->outbase,
          tj->outsuff);
  if (tj->TI == 1) {
    SummarizeOutputFiles(tj, 'I', 0, outp);
    SummarizeOutputFiles(tj, 'F', 1, outp);
  }
  else {
    SummarizeOutputFiles(tj, ' ', n, outp);
  }
  if (tj->OverwriteOutput == 1) {
    fprintf(outp, "\n!!! Output overwriting is ENABLED !!!\n");
  }

  // Reprint the input file
  HorizontalRule(outp, 0);
  fprintf(outp, "\nINPUT LINE TEXT:\n\n");
  PrintParagraph(tj->inpline, 79, NULL, outp);
  fprintf(outp, "\nINPUT FILE TEXT:\n\n");
  for (i = 0; i < tj->inptext.row; i++) {
    fprintf(outp, "%s", tj->inptext.map[i]);
  }
  HorizontalRule(outp, 1);

  // Compilation conditionals?

  // Print a report on the topology file
  HorizontalRule(outp, 0);
  fprintf(outp, "(1.) SYSTEM SPECIFICATIONS\n\n");
  if (tj->TI == 1) {
    fprintf(outp, "%% System I -->\n");
  }
  PrintSystemSummary(tp, outp);
  if (tj->TI == 1) {
    fprintf(outp, "%% System F -->\n");
    PrintSystemSummary(&tp[1], outp);
  }
  HorizontalRule(outp, 1);

  // Print control data
  HorizontalRule(outp, 0);
  fprintf(outp, "(2.) CONTROL DATA FOR THE RUN\n\n");
  fprintf(outp, "Mode of operation: ");
  if (tj->mode == 0) {
    fprintf(outp, "MOLECULAR DYNAMICS\n");
  }
  else {
    fprintf(outp, "ENERGY MINIMIZATION\n");
  }
  fprintf(outp, "\nNature and format of input:\n");
  if (tj->irest == 0) {
    fprintf(outp, " - Dynamics are beginning from input coordinates.\n");
  }
  else {
    fprintf(outp, " - Dynamics are being restarted from previous "
            "checkpoint.\n");
  }
  fprintf(outp, "\nNature and format of output:\n");
  if (tj->ntpr > 0) {
    fprintf(outp, " - State information / diagnostics: %7d steps\n", tj->ntpr);
  }
  if (tj->ntwr > 0) {
    fprintf(outp, " - Restart / checkpoint file:       %7d steps\n", tj->ntwr);
  }
  if (tj->ntwx > 0) {
    fprintf(outp, " - Coordinate sets output:          %7d steps\n", tj->ntwx);
  }
  if (tj->ntwv > 0) {
    fprintf(outp, " - Velocity sets output:            %7d steps\n", tj->ntwv);
  }
  if (tj->ntwf > 0) {
    fprintf(outp, " - Force sets output:               %7d steps\n", tj->ntwf);
  }
  HorizontalRule(outp, 1);  

  // Print parameters for nonbonded interactions
  HorizontalRule(outp, 0);  
  fprintf(outp, "(3.) NONBONDED INTERACTIONS\n\n");
  if (tp->ifbox > 0) {
    fprintf(outp, "The direct space sum is computed on a grid of linked "
            "cells, each at least\n%.2lfA thick in all dimensions.  Cell grid "
            "details:\n\n", dcinp->Mcut);
    fprintf(outp, " - Cell count:       %6d %6d %6d\n", CG->ng[0], CG->ng[1],
            CG->ng[2]);
    fprintf(outp, " - Cell atom limit:         %6d\n", CG->maxatom);
    fprintf(outp, " - Cell max atoms, pairs:   %6d %6d\n", maxCatom, maxCpair);
    fprintf(outp, " - Cell min atoms, pairs:   %6d %6d\n\n", minCatom,
            minCpair);
    if (rcinp->nlev == 1) {
      fprintf(outp, "The reciprocal space sum is computed by traditional "
              "Smooth Particle Mesh Ewald.\n"
              " - Mesh grid size:       %6d %6d %6d\n"
              " - Interpolation order:  %6d %6d %6d\n", rcinp->ng[0],
              rcinp->ng[1], rcinp->ng[2], rcinp->ordr[0], rcinp->ordr[1],
              rcinp->ordr[2]);
    }
    else {
      fprintf(outp, "The reciprocal space sum is computed by Multi-Level "
	      "Ewald.\n - Interpolation order:  %6d %6d %6d\n", rcinp->ordr[0],
              rcinp->ordr[1], rcinp->ordr[2]);
     for (i = 0; i < rcinp->nlev; i++) {
        fprintf(outp,
                " - Level %2d mesh size:   %6d %6d %6d [ %6.4lf x coarsened, ",
                i+1, rcinp->QL[i].pag, rcinp->QL[i].row, rcinp->QL[i].col,
                rcinp->cfac[i]);
        if (i < rcinp->nlev-1) {
          fprintf(outp, "%d padded ]\n", rcinp->PadYZ[i]);
        }
        else {
          fprintf(outp, "unpadded ]\n");
        }
      }
    }
    fprintf(outp, " - Direct Sum Tolerance: %20.13e\n"
            " - Ewald coefficient:    %20.13e\n"
            " - Gaussian width:       %20.13e\n\n", dcinp->Dtol,
	    dcinp->ewcoeff, rcinp->S);
    PrintParagraph("The direct space electrostatic force is approximated "
		   "using cubic spline interpolation. Accuracy of the "
		   "interpolation table is not tested for r < 1.0 Angstroms; "
		   "the statistics below reflect the inaccuracy that is "
		   "likely to be encountered in a typical nonbonded force "
		   "calculation.  The lookup table is indexed by the square "
		   "of the distance r_ij between charges i and j.", 79, NULL,
		   outp);
    if (Etab->FitType == 0) {
      fprintf(outp, " - Spline fitting method: BEST POSSIBLE FIT\n");
    }
    else {
      fprintf(outp, " - Spline fitting method: CONTINUOUS DERIVATIVE\n");
    }
    fprintf(outp, " - Spline discretization:  %12.5e\n"
            " - Maximum relative error: %12.5e at r_ij = %7.3lf\n"
            " - Maximum absolute error: %12.5e at r_ij = %7.3lf\n\n",
            Etab->dr, Etab->fmaxrelerr, Etab->fmaxrelerrloc, Etab->fmaxerr,
            Etab->fmaxerrloc);
    if (Etab->FitType == 0) {
      PrintParagraph("The direct space electrostatic energy is also "
                     "approximated using cubic spline interpolation. Note "
		     "that, because there are separate tables for force and "
		     "energy, the electrostatic force is not exactly the "
		     "derivative of the energy.  However, this approach "
		     "maximizes the accuracy of the lookup table for a given "
		     "amount of storage, which leads to acceptable energy "
		     "conservation nonetheless.", 79, NULL, outp);
      fprintf(outp, " - Spline fitting method: BEST POSSIBLE FIT\n");
    }
    else {
      fprintf(outp, " - Spline fitting method: CONTINUOUS DERIVATIVE\n");
    }
    fprintf(outp, " - Spline discretization:  %12.5e\n"
            " - Maximum relative error: %12.5e at r_ij = "
            "%7.3lf\n - Maximum absolute error: %12.5e at r_ij = %7.3lf\n\n",
            Etab->dr, Etab->umaxrelerr, Etab->umaxrelerrloc, Etab->umaxerr,
            Etab->umaxerrloc);
    fprintf(outp, "Other direct space parameters:\n"
            " - Electrostatic cutoff: %7.4lf\n"
            " - van der Waals cutoff: %7.4lf\n", dcinp->Ecut, dcinp->Vcut);
    if (dcinp->LRvdw == 1) {
      fprintf(outp, " - Long-ranged vdW interactions approximated by "
              "homogeneity approximation.\n\n");
    }
    else {
      fprintf(outp, " - Long-ranged vdW interactions are not treated.\n\n");
    }
  }
  else {
    fprintf(outp, "Simulation entails a non-periodic system with no "
	    "cutoff.\n\n");
#ifdef CUDA    
    PrintParagraph("The CUDA implicit solvent simulation module in mdgx "
		   "is designed for maximum throughput with many small "
		   "systems running concurrently.", 79, NULL, outp);
#else
    PrintParagraph("Note:  This is Not a CUDA executable."
                   "The CPU-based implicit solvent simulation module in mdgx "
		   "is designed as a means for checking the high throughput "
                   "CUDA engine.  It will not produce significantly better "
                   "trajectories or speeds than sander or pmemd CPU-based "
                   "code.  If you are not using this to debug or check system "
                   "behavior, the CUDA version is free, too, and will run on "
		   "most NVIDIA GPUs.", 79, NULL, outp);
#endif
    PrintParagraph("Forces in this module are computed in 32-bit floating "
                   "point precision and accumulated in 32-bit integers.  This "
                   "places strict limits on the sizes of forces that can be "
                   "represented: +/- 4096 kcal/mol-A.  Forces that fall "
                   "outside this range will wrap the format and introduce "
                   "wild errors in the dynamics, but no well-behaved system "
                   "suitable for simulations in Amber should display such "
                   "energy gradients.", 79, NULL, outp);
#ifdef CUDA
    PrintParagraph("Energies are likewise accumulated in 32-bit integers, and "
                   "therefore rounded at several levels.  Bond, angle, and "
                   "dihedral terms accumulate energy with a precision of "
                   "1/65536 of one kcal/mol per interaction.  The sum of many "
                   "interactions will display error in accordance will all "
                   "such terms' energies being rounded towards zero by such "
                   "an increment.  Non-bonded interactions accumulate in "
                   "local tiles at a precision of 1/16384 of one kcal/mol, "
                   "then get rounded to 1 part in 4096 kcal/mol for summation "
                   "into the complete value.  Energies produced by the CPU "
                   "and GPU versions of the code are expected to match to "
                   "within a few thousandths of a kcal/mol for bonds, angles, "
                   "and dihedrals, and to within a few hundredths of a "
                   "kcal/mol for non-bonded interactions.  The dynamics is "
                   "identical regardless of the precision of energies.", 79,
                   NULL, outp);
#else
    PrintParagraph("Energies are computed in 32-bit floating point precision "
                   "and accumulated in 64-bit double precision.", 79, NULL,
                   outp);
#endif
  }
  HorizontalRule(outp, 1);

  // Check the topology and conformation for problems
  HorizontalRule(outp, 0);
  fprintf(outp, " (4.) SYSTEM TOPOLOGY AND CONFORMATION\n\n");
  if (tj->topchk == 1) {
    if (tj->TI == 1) {
      fprintf(outp, "%% System I -->\n");
    }
    fprintf(outp, " - Checking chirality of standard amino acids.\n");
    iwarning = ProteinChiralityCheck(tp, crd, outp);
    fprintf(outp, " - Checking for omitted disulfides.\n");
    iwarning += FindDisulfides(tp, crd, outp);
    fprintf(outp, " - Checking for (nonstandard) Lennard-Jones combining "
            "rules.\n");
    iwarning += CheckLJRules(tp, outp);
    if (iwarning == 0) {
      fprintf(outp, " - No problems with the topology detected.\n\n");
    }
    if (tj->TI == 1) {
      fprintf(outp, "%% System F -->\n");
      fprintf(outp, " - Checking chirality of standard amino acids.\n");
      iwarning = ProteinChiralityCheck(tp, crd, outp);
      fprintf(outp, " - Checking for omitted disulfides.\n");
      iwarning += FindDisulfides(tp, crd, outp);
      if (iwarning == 0) {
        fprintf(outp, " - No problems with the topology detected.\n\n");
      }
    }
  }
  else {
    fprintf(outp, " - Check skipped.\n");
  }
  HorizontalRule(outp, 1);
  if (tp->ifbox > 0) {
    HorizontalRule(outp, 0);
    fprintf(outp, " (5.) RESULTS\n\n");
  }
  
  // Close the diagnostics file (it will be reopened
  // when more output is ready)                     
  fclose(outp);
}

//-----------------------------------------------------------------------------
// ComputeAverageEnergies: compute averages and standard deviations for 
//                         quantities stored in the Energy struct.      
//                                                                      
// Arguments:                                                           
//   tj:     trajectory control information                             
//   sysUV:  the energy and virial information                          
//-----------------------------------------------------------------------------
static void ComputeAverageEnergies(trajcon *tj, Energy *sysUV)
{
  long long int stepcount;
  double npt, afac;

  stepcount = (tj->nfistep > 0) ? tj->nfistep : tj->nstep;
  if (tj->ntpr > 0) {
    npt = stepcount/tj->ntpr;
    afac = 1.0/npt;
  }
  else {
    npt = 0.0;
    afac = 1.0;
  }
  sysUV->AVET *= afac;
  sysUV->AVEP *= afac;
  sysUV->AVEV *= afac;
  sysUV->AVEetot *= afac;
  sysUV->AVEkine *= afac;
  sysUV->AVEeptot *= afac;
  sysUV->AVEbond *= afac;
  sysUV->AVEangl *= afac;
  sysUV->AVEdihe *= afac;
  sysUV->AVEcmap *= afac;
  sysUV->AVEelec *= afac;
  sysUV->AVEdelec *= afac;
  sysUV->AVErelec *= afac;
  sysUV->AVEvdw *= afac;
  if (tj->TI == 1) {
    sysUV->AVEdVdL /= stepcount;
  }
  afac = (tj->ntpr > 0) ? 1.0/(npt-1.0) : 1.0;
  sysUV->RMST = SafeSqrt((sysUV->RMST - npt*sysUV->AVET*sysUV->AVET)*afac);
  sysUV->RMSP = SafeSqrt((sysUV->RMSP - npt*sysUV->AVEP*sysUV->AVEP)*afac);
  sysUV->RMSV = SafeSqrt((sysUV->RMSV - npt*sysUV->AVEV*sysUV->AVEV)*afac);
  sysUV->RMSetot = 
    SafeSqrt((sysUV->RMSetot - npt*sysUV->AVEetot*sysUV->AVEetot)*afac);
  sysUV->RMSkine = 
    SafeSqrt((sysUV->RMSkine - npt*sysUV->AVEkine*sysUV->AVEkine)*afac);
  sysUV->RMSeptot =
    SafeSqrt((sysUV->RMSeptot - npt*sysUV->AVEeptot*sysUV->AVEeptot)*afac);
  sysUV->RMSbond =
    SafeSqrt((sysUV->RMSbond - npt*sysUV->AVEbond*sysUV->AVEbond)*afac);
  sysUV->RMSangl =
    SafeSqrt((sysUV->RMSangl - npt*sysUV->AVEangl*sysUV->AVEangl)*afac);
  sysUV->RMSdihe =
    SafeSqrt((sysUV->RMSdihe - npt*sysUV->AVEdihe*sysUV->AVEdihe)*afac);
  sysUV->RMScmap =
    SafeSqrt((sysUV->RMScmap - npt*sysUV->AVEcmap*sysUV->AVEcmap)*afac);
  sysUV->RMSelec =
    SafeSqrt((sysUV->RMSelec - npt*sysUV->AVEelec*sysUV->AVEelec)*afac);
  sysUV->RMSdelec =
    SafeSqrt((sysUV->RMSdelec - npt*sysUV->AVEdelec*sysUV->AVEdelec)*afac);
  sysUV->RMSrelec =
    SafeSqrt((sysUV->RMSrelec - npt*sysUV->AVErelec*sysUV->AVErelec)*afac);
  sysUV->RMSvdw =
    SafeSqrt((sysUV->RMSvdw - npt*sysUV->AVEvdw*sysUV->AVEvdw)*afac);
  if (tj->TI == 1) {
    sysUV->RMSdVdL =
      SafeSqrt((sysUV->RMSdVdL -
                (stepcount-1.0)*sysUV->AVEdVdL*sysUV->AVEdVdL) / stepcount);
  }
}

//-----------------------------------------------------------------------------
// PrintMeanEnergies: print out the mean energies.                      
//                                                                      
// Arguments:                                                           
//   sysUV:    the energy and virial information                        
//   density:  the system density                                       
//   tp:      system topology (contains the CMAP count)
//   outp:     the output diagnostics file (it's being closed up)       
//   format:   code for the format to write, 0 = PME, periodic simulation,
//             1 = isolated system simulation (with continuum solvation model)
//-----------------------------------------------------------------------------
static void PrintMeanEnergies(Energy *sysUV, double density, prmtop *tp,
			      FILE *outp, int format)
{
  char cmapline[96];

  // CMAP terms are only printed if they are present
  if (tp->ncmap > 0) {
    sprintf(cmapline, "%% CMAP: %17.4lf\n", sysUV->AVEcmap);
  }
  else {
    cmapline[0] = '\0';
  }
  if (format == 0) {
    fprintf(outp, "%% Temperature:%11.4lf   Pressure: %13.4lf   "
	    "Volume: %15.4lf\n%% Etot: %17.4lf   EKtot: %16.4lf   "
	    "EPtot: %16.4lf\n%% Bond: %17.4lf   Angle: %16.4lf   "
	    "Dihedral: %13.4lf\n%s%% Elec: %17.4lf   vdW:  %17.4lf   "
	    "Virial: %15.4lf\n%% Density: %14.4lf\n", sysUV->AVET,
	    sysUV->AVEP, sysUV->AVEV, sysUV->AVEetot, sysUV->AVEkine,
	    sysUV->AVEeptot, sysUV->AVEbond, sysUV->AVEangl, sysUV->AVEdihe,
	    cmapline, sysUV->AVEelec, sysUV->AVEvdw, sysUV->AVEVir, density);
  }
  else if (format == 1) {
    fprintf(outp, "%% Temperature:%11.4lf\n%% Etot: %17.4lf   EKtot: "
	    "%16.4lf   EPtot: %16.4lf\n%% Bond: %17.4lf   Angle: %16.4lf   "
	    "Dihedral: %13.4lf\n%s%% Elec: %17.4lf   vdW:  %17.4lf   "
	    "Solvent: %14.4lf\n", sysUV->AVET, sysUV->AVEetot, sysUV->AVEkine,
	    sysUV->AVEeptot, sysUV->AVEbond, sysUV->AVEangl, sysUV->AVEdihe,
	    cmapline, sysUV->AVEelec, sysUV->AVEvdw, sysUV->AVErelec);
  }
}

//-----------------------------------------------------------------------------
// PrintRMSEnergies: print out root mean squared deviations in energies.
//                                                                      
// Arguments:                                                           
//   sysUV:   the energy and virial information                          
//   tp:      system topology (contains the CMAP count)
//   outp:    the output diagnostics file (it's being closed up)         
//   format:  code for the format to write, 0 = PME, periodic simulation,
//            1 = isolated system simulation (with continuum solvation model)
//-----------------------------------------------------------------------------
static void PrintRMSEnergies(Energy *sysUV, prmtop *tp, FILE *outp, int format)
{
  char cmapline[96];
  
  // CMAP terms are only printed if they exist
  if (tp->ncmap > 0) {
    sprintf(cmapline, "%% CMAP: %17.4lf\n", sysUV->RMScmap);
  }
  else {
    cmapline[0] = '\0';
  }
  if (format == 0) {
   fprintf(outp, "%% Temperature:%11.4lf   Pressure: %13.4lf   Volume: %15.4lf"
           "\n%% Etot: %17.4lf   EKtot: %16.4lf   EPtot: %16.4lf\n"
           "%% Bond: %17.4lf   Angle: %16.4lf   Dihedral: %13.4lf\n%s"
           "%% Elec: %17.4lf   vdW:  %17.4lf   Virial: %15.4lf\n",
           sysUV->RMST, sysUV->RMSP, sysUV->RMSV, sysUV->RMSetot,
           sysUV->RMSkine, sysUV->RMSeptot, sysUV->RMSbond, sysUV->RMSangl,
           sysUV->RMSdihe, cmapline, sysUV->RMSelec, sysUV->RMSvdw,
	   sysUV->RMSVir);
  }
  else if (format == 1) {
    fprintf(outp, "%% Temperature:%11.4lf\n"
            "%% Etot: %17.4lf   EKtot: %16.4lf   EPtot: %16.4lf\n"
            "%% Bond: %17.4lf   Angle: %16.4lf   Dihedral: %13.4lf\n%s"
            "%% Elec: %17.4lf   vdW:  %17.4lf   Solvent: %14.4lf\n",
            sysUV->RMST, sysUV->RMSetot, sysUV->RMSkine, sysUV->RMSeptot,
	    sysUV->RMSbond, sysUV->RMSangl, sysUV->RMSdihe, cmapline,
	    sysUV->RMSelec, sysUV->RMSvdw, sysUV->RMSrelec);
  }
}

//-----------------------------------------------------------------------------
// CloseDiagnosticsFile: close up a diagnostics file, after printing averages
//                       and timing information.               
//                                                                      
// Arguments:                                                           
//   tj:      trajectory input control data
//   tp:      the system topology
//   rcinp:   reciprocal space calculation control data
//   sysUV:   system energy and virial accumulators
//   etimers: timings data and accumulators
//   CG:      the cell grid for direct space decomposition
//   format:  code for the format to write, 0 = PME, periodic simulation,
//            1 = isolated system simulation (with continuum solvation model)
//-----------------------------------------------------------------------------
void CloseDiagnosticsFile(trajcon *tj, prmtop *tp, reccon *rcinp,
			  Energy *sysUV, execon *etimers, cellgrid* CG,
			  int format)
{
  long long int stepcount;
  double density;
  char fname[MAXNAME];
  FILE *outp;
  time_t ct;
  dmat alltime;

  // Gather timing data; if this is not the master
  // process, free allocated memory and bail out  
  alltime = GatherTimingData(etimers, CG, rcinp);
  if (CG->tid != 0) {
    DestroyDmat(&alltime);
    return;
  }

  // Deprecate the step count to ensure we close the proper file
  SpliceFileName(tj, tj->outbase, tj->outsuff, fname, 1);
  outp = fopen(fname, "a");
  fprintf(outp, "\n");

  // Print out averages
  stepcount = (tj->nfistep > 0) ? tj->nfistep : tj->nstep;
  if (tj->TI == 1) {
    sysUV[1].AVEdVdL = sysUV[0].AVEdVdL;
    sysUV[1].RMSdVdL = sysUV[0].RMSdVdL;
  }
  ComputeAverageEnergies(tj, sysUV);
  density = DSum(tp->Masses, tp->natom)/(sysUV->AVEV*AVOGADRO*1.0e-24);
  fprintf(outp, "\n%% Averages over %lld steps sampled at %lld points:\n%%\n",
          stepcount, ((tj->ntpr > 0) ? stepcount / tj->ntpr : 1));
  if (tj->TI == 1) {
    fprintf(outp, "%% SYSTEM I -->\n");
  }
  PrintMeanEnergies(sysUV, density, tp, outp, format);
  if (tj->TI == 1) {
    ComputeAverageEnergies(tj, &sysUV[1]);
    fprintf(outp, "%% SYSTEM F -->\n");
    PrintMeanEnergies(&sysUV[1], density, tp, outp, format);
    fprintf(outp, "%%\n%% dV/dL: %16.4lf\n", sysUV[0].AVEdVdL);
  }
  fprintf(outp, "\n%% Root mean squared deviations over %lld steps sampled at "
          "%lld points:\n%%\n", stepcount,
          ((tj->ntpr > 0) ? stepcount / tj->ntpr : 1));
  if (tj->TI == 1) {
    fprintf(outp, "%% SYSTEM I -->\n");
  }
  PrintRMSEnergies(sysUV, tp, outp, format);
  if (tj->TI == 1) {
    fprintf(outp, "%% SYSTEM F -->\n");
    PrintRMSEnergies(sysUV, tp, outp, format);
    fprintf(outp, "%%\n%% dV/dL: %16.4lf (computed and averaged over "
            "%lld steps)\n", sysUV[0].RMSdVdL, stepcount);
  }
  fprintf(outp, "\n");
  HorizontalRule(outp, 1);
  HorizontalRule(outp, 0);
  ct = time(&ct);
  if (tj->mode != 5) {
    fprintf(outp, " (%d.) TIMINGS\n\n", 6 + (tp->ifbox == 0));
    fprintf(outp, " Run completed on %s\n", asctime(localtime(&ct)));
  }
  else {
    fprintf(outp, " (%d.) Timings for dynamics\n\n", 6 + (tp->ifbox == 0));
    fprintf(outp, " Dynamics completed on %s\n", asctime(localtime(&ct)));
  }
  PrintTimingData(rcinp, &alltime, outp);

  // Print the closing seal.  This will serve as a flag to tell mdgx
  // on any subsequent runs that this segment was completed.  Do not
  // add a closing seal to runs performed for the purpose of fitting
  // IPolQ charges; there's more to do, and the runs cannot be      
  // checkpointed in any traditional sense.                         
  HorizontalRule(outp, (tj->mode == 5));
  if (tj->mode != 5) {
    fprintf(outp, "$! Closing watermark.  DO NOT alter this segment of the "
            "file. !$\n");
    HorizontalRule(outp, 0);
  }
  fclose(outp);

  // Free allocated memory on the master process
  DestroyDmat(&alltime);
}

//-----------------------------------------------------------------------------
// InitAverageEnergy: (re)initialize the accumulators for mean and standard
//                    deviation energy statistics.  All of the accumulators are
//                    initialized to zero, as a subsequent call to
//                    AccAverageEnergy (below) will add in the current values
//                    of all these quantities.  There is one exception,
//                    however, relating to dV/dL in thermodynamic integration.
//                    The accumuators for dV/dL are initialized to their
//                    current values here, as AccAverageEnergy, which only
//                    accumulates on selected steps as instructed by the ntpr
//                    input parameter, does not handle accumulation of dV/dL. 
//                    When TI is active, dV/dL is accumulated at every  
//                    simulation step in the integration routines.      
//                                                                      
// Arguments:                                                           
//   sysUV:    system energy and virial data (state information)        
//-----------------------------------------------------------------------------
void InitAverageEnergy(Energy *sysUV)
{
  sysUV->AVEvdw = 0.0;
  sysUV->AVEbond = 0.0;
  sysUV->AVEangl = 0.0;
  sysUV->AVEdihe = 0.0;
  sysUV->AVEcmap = 0.0;
  sysUV->AVEkine = 0.0;
  sysUV->AVEelec = 0.0;
  sysUV->AVErelec = 0.0;
  sysUV->AVEdelec = 0.0;
  sysUV->AVEeptot = 0.0;
  sysUV->AVEetot = 0.0;
  sysUV->AVEP = 0.0;
  sysUV->AVEV = 0.0;
  sysUV->AVET = 0.0;
  sysUV->AVEVir = 0.0;
  sysUV->AVEdVdL = sysUV->dVdL;
  sysUV->RMSvdw = 0.0;
  sysUV->RMSbond = 0.0;
  sysUV->RMSangl = 0.0;
  sysUV->RMSdihe = 0.0;
  sysUV->RMScmap = 0.0;
  sysUV->RMSkine = 0.0;
  sysUV->RMSelec = 0.0;
  sysUV->RMSrelec = 0.0;
  sysUV->RMSdelec = 0.0;
  sysUV->RMSeptot = 0.0;
  sysUV->RMSetot = 0.0;
  sysUV->RMSP = 0.0;
  sysUV->RMSV = 0.0;
  sysUV->RMST = 0.0;
  sysUV->RMSVir = 0.0;
  sysUV->RMSdVdL = sysUV->dVdL * sysUV->dVdL;
}

//-----------------------------------------------------------------------------
// AccAverageEnergy: accumulators mean standard deviation energy stats.  This
//                   routine does not accumulate dV/dL statistics required by
//                   thermodynamic integration; those are accumulated at every
//                   step in the dynamics routines.
//                                                                      
// Arguments:                                                           
//   sysUV:    system energy and virial data (state information)        
//-----------------------------------------------------------------------------
void AccAverageEnergy(Energy *sysUV)
{
  sysUV->AVEvdw   += sysUV->vdw12 + sysUV->vdw6;
  sysUV->AVEbond  += sysUV->bond;
  sysUV->AVEangl  += sysUV->angl;
  sysUV->AVEdihe  += sysUV->dihe;
  sysUV->AVEcmap  += sysUV->cmap;
  sysUV->AVEkine  += sysUV->kine;
  sysUV->AVEelec  += sysUV->elec;
  sysUV->AVErelec += sysUV->relec;
  sysUV->AVEdelec += sysUV->delec;
  sysUV->AVEeptot += sysUV->eptot;
  sysUV->AVEetot  += sysUV->etot;
  sysUV->AVEP     += sysUV->P;
  sysUV->AVEV     += sysUV->V;
  sysUV->AVET     += sysUV->T;
  sysUV->AVEVir   += SumTrace3(sysUV->Vir);
  sysUV->RMSvdw   += pow(sysUV->vdw12 + sysUV->vdw6, 2.0);
  sysUV->RMSbond  += sysUV->bond * sysUV->bond;
  sysUV->RMSangl  += sysUV->angl * sysUV->angl;
  sysUV->RMSdihe  += sysUV->dihe * sysUV->dihe;
  sysUV->RMScmap  += sysUV->cmap * sysUV->cmap;
  sysUV->RMSkine  += sysUV->kine * sysUV->kine;
  sysUV->RMSelec  += sysUV->elec * sysUV->elec;
  sysUV->RMSrelec += sysUV->relec * sysUV->relec;
  sysUV->RMSdelec += sysUV->delec * sysUV->delec;
  sysUV->RMSeptot += sysUV->eptot * sysUV->eptot;
  sysUV->RMSetot  += sysUV->etot * sysUV->etot;
  sysUV->RMSP     += sysUV->P * sysUV->P;
  sysUV->RMSV     += sysUV->V * sysUV->V;
  sysUV->RMST     += sysUV->T * sysUV->T;
  sysUV->RMSVir   += pow(SumTrace3(sysUV->Vir), 2.0);
}

//-----------------------------------------------------------------------------
// PrintStateInfo: print the instantaneous information concerning the state
//                 variables of the system.
//
// Arguments:
//   sysUV:    system energy and virial data (state information)
//   tj:       trajectory control data
//   tp:       system topology (tells us whether there are CMAPs)
//   outp:     output diagnostics file
//   format:   code for the format to write, 0 = PME, periodic simulation,
//             1 = isolated system simulation (with continuum solvation model)
//-----------------------------------------------------------------------------
void PrintStateInfo(Energy *sysUV, trajcon *tj, prmtop *tp, FILE *outp,
		    int format)
{
  char cmapline[96];
  
  // CMAP terms are only printed if they exist
  if (tp->ncmap > 0) {
    sprintf(cmapline, " CMAP: %17.4lf\n", sysUV->cmap);
  }
  else {
    cmapline[0] = '\0';
  } 
  if (format == 0) {
    fprintf(outp, " Step:%18lld   Time: %17.4lf\n Temperature:"
            "%11.4lf   Pressure: %13.4lf   Volume: %15.4lf\n"
            " Etot: %17.4lf   EKtot: %16.4lf   EPtot: %16.4lf\n"
            " Bond: %17.4lf   Angle: %16.4lf   Dihedral: %13.4lf\n%s"
            " Elec: %17.4lf   vdW:  %17.4lf   Virial: %15.4lf\n", tj->currstep,
            tj->currtime, sysUV->T, sysUV->P, sysUV->V, sysUV->etot,
	    sysUV->kine, sysUV->eptot, sysUV->bond, sysUV->angl, sysUV->dihe,
	    cmapline, sysUV->elec, sysUV->vdw12 + sysUV->vdw6,
	    SumTrace3(sysUV->Vir));
  }
  else if (format == 1) {
    fprintf(outp, " Step:%18lld   Time: %17.4lf\n Temperature:"
            "%11.4lf\n"
            " Etot: %17.4lf   EKtot: %16.4lf   EPtot: %16.4lf\n"
            " Bond: %17.4lf   Angle: %16.4lf   Dihedral: %13.4lf\n%s"
            " Elec: %17.4lf   vdW:  %17.4lf   Solvent: %14.4lf\n",
	    tj->currstep, tj->currtime, sysUV->T, sysUV->etot,
	    sysUV->kine, sysUV->eptot, sysUV->bond, sysUV->angl, sysUV->dihe,
	    cmapline, sysUV->elec, sysUV->vdw12 + sysUV->vdw6, sysUV->relec);
  }
}

//-----------------------------------------------------------------------------
// OpenNewOutput: this function determines whether to open a new output file.
//                The conditions for doing so are laid out within the function.
//                                                                      
// Arguments:                                                           
//   tj:       trajectory control data                                  
//   ival:     the interval with which the output is written            
//-----------------------------------------------------------------------------
static int OpenNewOutput(trajcon *tj, int ival)
{
  // If this is the beginning of the simulation, open a new file
  if (tj->currstep == 0 || (tj->irest == 1 && tj->currstep == tj->ntwx) ||
      (tj->ntpr == 0 && (tj->currstep == tj->nstep ||
                         (tj->nfistep > 0 && tj->currstep == tj->nfistep)))) {
    return 1;
  }

  // If there are multiple segments to write, a new file may
  // need to be opened under certain circumstances          
  if (tj->nfistep > 0) {

    // One frame per file is a particular case
    if (tj->nfistep == 1) {
      return 1;
    }
    else if (tj->currfi > 0 && tj->currstep % tj->nfistep == ival) {
      return 1;
    }
  }

  // Otherwise a new diagnostics file need not be written
  return 0;
}

//-----------------------------------------------------------------------------
// WriteDiagnostics: write the diagnostic output file in AMBER format.  This
//                   is the file specified by "-o" (default name mdout) in
//                   SANDER/PMEMD.                            
//                                                                      
// Arguments:                                                           
//   tj:       trajectory contol data                                   
//   tp:       the topology                                             
//   dcinp:    direct space control data                                
//   rcinp:    reciprocal space control data                            
//   Etab:     the electrostatic direct space spline interpolation table
//   CG:       the cell grid                                            
//   crd:      the coordinates                                          
//   sysUV:    system energy and virial data (state information)        
//   etimers:  execution time counters                                  
//-----------------------------------------------------------------------------
void WriteDiagnostics(trajcon *tj, prmtop *tp, dircon *dcinp, reccon *rcinp,
                      FrcTab *Etab, cellgrid *CG, coord *crd, Energy *sysUV,
                      execon *etimers, int n)
{
  char fname[MAXNAME];
  FILE *outp;

  // The volume, temperature, and current system pressure are computed
  // here.  A conversion factor of 6.9477e4 puts the pressure in units
  // of bar; if the barostat is Monte-Carlo or non-existent, pressure 
  // is set to zero.  Otherwise, the pressure should reflect only the 
  // kinetic contributions and read way too high.                     
  sysUV->V = crd->invU.data[0]*crd->invU.data[4]*crd->invU.data[8];
  if (tj->ntt == 0 || tj->ntt == 3) {

    // If there is no thermostat (or, if there is a 
    // Langevin integrator), the system kinetic     
    // energy and temperature must be computed here.
    if (tj->TI == 1) {
      sysUV->T = SystemTemperatureTI(CG, crd, tp, sysUV, tj, 1);
    }
    else {
      sysUV->T = SystemTemperature(CG, crd, tp, sysUV, tj, 1);
    }
  }
  sysUV->P = (tj->ntp > 0 && tj->barostat == 1) ? 
    CurrentSystemPressure(sysUV, crd)*PCONVFAC : 0.0;
#ifdef MPI
  SumTotalEnergy(CG, sysUV);
#else
  SumTotalEnergy(sysUV);
#endif
  if (tj->TI == 1) {
    sysUV[1].V = sysUV[0].V;
    sysUV[1].T = sysUV[0].T;
    sysUV[1].P = (tj->ntp > 0 && tj->barostat == 1) ? 
      CurrentSystemPressure(&sysUV[1], &crd[1])*PCONVFAC : 0.0;
#ifdef MPI
    SumTotalEnergy(&CG[1], &sysUV[1]);
#else
    SumTotalEnergy(&sysUV[1]);
#endif
  }

  // Determine whether this is a new file.  New files must be
  // opened at the beginning of a simulation and when it is  
  // time to write the first frame of any additional segments
  if (OpenNewOutput(tj, tj->ntpr) == 1) {

    // Begin a new file
    OpenDiagnosticsFile(tj, tp, dcinp, rcinp, Etab, crd, CG, n);

    // (Re)Initialize running averages and standard deviations
    InitAverageEnergy(sysUV);
    if (tj->TI == 1) {
      InitAverageEnergy(&sysUV[1]);
    }
  }

  // Compute and print running averages and standard deviations
  if ((tj->ntpr > 0 && tj->currstep % tj->ntpr == 0) ||
      (tj->ntpr == 0 && 
       (tj->currstep == tj->nstep ||
        (tj->nfistep > 0 && tj->currstep % tj->nfistep == 0)))) {

    // Accumulate the average energy; the  
    // initial configuration does not count
    if (tj->currstep > 0) {
      AccAverageEnergy(sysUV);
      if (tj->TI == 1) {
        AccAverageEnergy(&sysUV[1]);
      }
    }
#ifdef MPI
    if (CG->tid == 0) {
#endif
      SpliceFileName(tj, tj->outbase, tj->outsuff, fname, 1);
      outp = fopen(fname, "a");
      fprintf(outp, "\\\\\\\n");
      if (tj->TI == 1) {
        fprintf(outp, "%% System I -->\n");
      }
      PrintStateInfo(sysUV, tj, tp, outp, 0);
      if (tj->TI == 1) {
        fprintf(outp, "%% System F -->\n");
        PrintStateInfo(&sysUV[1], tj, tp, outp, 0);
        fprintf(outp, " dV/dL: %16.4lf\n", sysUV[0].dVdL);
      }
      fclose(outp);
#ifdef MPI
    }
#endif
  }

  // Finalize the diagnostics file if all of the relevant   
  // data has been accumulated in this file.  This happens  
  // if all of the output steps that call for printing      
  // diagnostics have passed.  If thermodynamic integration 
  // is active (tj->TI == 1), then the diagnostic file      
  // cannot be closed until the final step of the simulation
  // as dV/dL is accumulated over all steps.                
  if ((tj->TI != 1 && tj->nfistep > 0 && tj->currstep > 0 &&
       tj->currstep == tj->currfi*tj->nfistep) ||
      (tj->TI != 1 && tj->currstep == tj->nstep) ||
      (tj->TI == 1 && tj->nfistep > 0 && tj->currstep > 0 &&
       tj->currstep == tj->currfi*tj->nfistep) ||
      (tj->TI == 1 && tj->currstep == tj->nstep)) {
    CloseDiagnosticsFile(tj, tp, rcinp, sysUV, etimers, CG, 0);

    // Reset the timers, in case this file is but one
    // in a series that is to be printed.            
    InitExecon(etimers);
    mdgxStartTimer(etimers);
  }
}

//-----------------------------------------------------------------------------
// WriteRst: write a restart file in AMBER format.                      
//                                                                      
// Arguments:                                                           
//   CG:     the cell grid (if this is NULL, the behavior will change slightly
//           and an inpcrd file with no velocities will be written)
//   tc:     the coordinates                                            
//   tj:     the trajectory control parameters (if the mode is zero a restart
//           file including velocities must be written)
//   n:      write the restart file for the nth system
//   prbox:  flag to have box dimensions printed to the end of the restart
//           file.  Values include 0, DO NOT print; 1 ALWAYS print, -1 print
//           ONLY if ifbox == 1 in the topology.
//-----------------------------------------------------------------------------
void WriteRst(cellgrid *CG, coord *tc, prmtop *tp, trajcon *tj, int n,
              int prbox)
{
  int i, j, atmid, natm, nrec;
  double* crdcpy;
  double* halfvel;
  double *veltmp, *mtmp;
  char fname[MAXNAME];
  FILE *outp;
  cell *C;
  
  // Merge coordinates
#ifdef MPI
  if (CG != NULL) {
    CoordinateReduction(CG, tc, tp);
  }
#endif

  // Rewind velocities if a cell grid is provided.  If not, as is the case
  // with the multi-simulator, the velocities will have to be rewound by
  // another method.
  if (CG != NULL) {

    // Only the master process writes the trajectory
    if (CG->tid != 0) {
      return;
    }

    // The velocities must be saved after one half step 
    // in order for the restart to preserve the original
    // trajectory.  Therefore, rewind the velocities by 
    // half of the contribution from the latest force   
    // update, which will be reapplied immediately after
    // the restart occurs.                              
    halfvel = CpyDVec(tc->vel, 3*tc->natom);
    mtmp = tp->InvMasses;
    const double hdt = 0.5*418.4*tj->dt;
    for (i = 0; i < CG->ncell; i++) {
      C = &CG->data[i];
      for (j = 0; j < C->nr[0]; j++) {
        atmid = C->data[j].id;
        const double hmdt = hdt*mtmp[atmid];
        halfvel[3*atmid] -= hmdt*C->data[j].frc[0];
        halfvel[3*atmid+1] -= hmdt*C->data[j].frc[1];
        halfvel[3*atmid+2] -= hmdt*C->data[j].frc[2];
      }
    }
    veltmp = tc->vel;
    tc->vel = halfvel;
  }

  // Splice together the name of the restart file  
  // The step number is deprecated and returned    
  // to its original value so that the restart     
  // file of any particular segment reflects the   
  // state of the system at the END of the segment.
  SpliceFileName(tj, tj->rstbase.map[n], tj->rstsuff.map[n], fname, 0);
  outp = FOpenSafe(fname, tj->OverwriteOutput);
  
  // Reposition atom groups to avoid splitting
  // things across cell boundaries and fulfill
  // user-specified re-imaging of coordinates 
  crdcpy = CpyDVec(tc->loc, 3*tc->natom);
  // FIX ME!  Need implementation for iwrap!  

  // Write out positions
  natm = (tp->EPInserted == 0) ? tc->natom : tp->norigatom;
  if (CG != NULL || tj->mode == 0 || tj->mode == 9) {
    fprintf(outp, "\n%6d%18.9e\n", natm, tj->currtime);
  }
  else {
    fprintf(outp, "\n%6d\n", natm);
  }
  nrec = 0;
  for (i = 0; i < tc->natom; i++) {

    // Skip this atom if it is an inserted extra point
    if (tp->EPInserted == 1 && tp->OldAtomNum[i] == -1) {
      continue;
    }

    // If we're still here, this atom was
    // part of the original topology     
    fprintf(outp, "%12.7lf%12.7lf%12.7lf", crdcpy[3*i], crdcpy[3*i+1],
            crdcpy[3*i+2]);
    nrec++;
    if (nrec % 2 == 0) {
      fprintf(outp, "\n");
    }
  }
  if (natm % 2 != 0) {
    fprintf(outp, "\n");
  }

  // Write out velocities if PME is in effect (that's a pretty weird
  // justification!), in standard dynamics mode, or when doing dynamics
  // on arrays of small GB systems.
  if (CG != NULL || tj->mode == 0 || tj->mode == 9) {
    nrec = 0;
    for (i = 0; i < tc->natom; i++) {

      // Skip this atom if it is an inserted extra point
      if (tp->EPInserted == 1 && tp->OldAtomNum[i] == -1) {
        continue;
      }

      // If we're still here, this atom was
      // part of the original topology     
      fprintf(outp, "%12.7lf%12.7lf%12.7lf", tc->vel[3*i], tc->vel[3*i+1],
              tc->vel[3*i+2]);
      nrec++;
      if (nrec % 2 == 0) {
        fprintf(outp, "\n");
      }
    }
    if (natm % 2 != 0) {
      fprintf(outp, "\n");
    }
  }

  // Write box dimensions unless this is a GB array simulation.  This
  // should work off of some deeper parameter that states whether the
  // simulation has periodic boundaries.
  if (prbox == 1 || (prbox == -1 && tp->ifbox > 0)) {
    fprintf(outp, "%12.7lf%12.7lf%12.7lf%12.7lf%12.7lf%12.7lf\n", tc->gdim[0],
            tc->gdim[1], tc->gdim[2], tc->gdim[3]*180.0/PI,
	    tc->gdim[4]*180.0/PI, tc->gdim[5]*180.0/PI);
  }
  fclose(outp);

  // Replace the velocities with their original values
  if (CG != NULL) {
    tc->vel = veltmp;
    free(halfvel);
  }

  // Free the last of the allocated memory
  free(crdcpy);
}

//-----------------------------------------------------------------------------
// WriteCrd: write a set of coordinates to a trajectory file in .crd format.
//                                                                      
// Arguments:                                                           
//   CG:   the cell grid (if NULL, all of the coordinates will be assumed to
//         be contained in coord and no merger of coordinates or forces will
//         be performed)
//   crd:  the coordinates structure                                    
//   cvf:  flag to direct writing of coordinates (1), velocities (2), or
//         forces (3)                                                   
//   tj:   trajectory control information                               
//   tp:   the topology (in case extra points have been inserted)       
//   n:    the trajectory is written for the nth system                 
//-----------------------------------------------------------------------------
void WriteCrd(cellgrid *CG, coord *crd, int cvf, trajcon *tj, prmtop *tp,
              int n)
{
  int h, i, atmcon, ntwcon;
  int *tmpDomain;
  double *dtmp;
  char fname[MAXNAME];
  char *ctmp, *csuf;
  FILE *outp;

  // Merge coordinates if we're serious about writing the trajectory.
  // This routine may have been called merely to initialize the file.
#ifdef MPI
  if (tj->currstep > 0 && CG != NULL) {
    CoordinateReduction(CG, crd, tp);
  }
#endif

  // Bail out if this is not the master process of this cell grid
  if (CG != NULL && CG->tid != 0) {
    return;
  }

  // Set a pointer to avoid conditionals in the loop
  if (cvf == 1) {
    dtmp = crd->loc;
    ctmp = tj->trjbase.map[n];
    csuf = tj->trjsuff.map[n];
    ntwcon = tj->ntwx;
  }
  else if (cvf == 2) {
    dtmp = crd->vel;
    ctmp = tj->velbase.map[n];
    csuf = tj->velsuff.map[n];
    ntwcon = tj->ntwv;
  }
  else if (cvf == 3) {

    // The master process must be put in
    // charge of all cells temporarily  
    if (CG != NULL) {
      tmpDomain = CG->MyCellDomain;
      CG->MyCellDomain = CountUp(CG->ncell);
      h = CG->MyCellCount;
      CG->MyCellCount = CG->ncell;
      MapCellForcesToAtoms(CG, crd);
      free(CG->MyCellDomain);
      CG->MyCellDomain = tmpDomain;
      CG->MyCellCount = h;
    }
    dtmp = crd->frc;
    ctmp = tj->frcbase.map[n];
    csuf = tj->frcsuff.map[n];
    ntwcon = tj->ntwf;
  }
  
  // Open the output file
  SpliceFileName(tj, ctmp, csuf, fname, 1);
  if (OpenNewOutput(tj, ntwcon) == 1) {

    // Open a new coordinates file
    outp = FOpenSafe(fname, tj->OverwriteOutput);
    time_t ct;
    ct = time(&ct);
    fprintf(outp, "Trajectory written by mdgx on %s", asctime(localtime(&ct)));
  }
  else {
    outp = fopen(fname, "a");
  }

  // Do not write trajectory at step 0
  if (tj->currstep == 0) {
    fclose(outp);
    return;
  }

  // Write trajectory output
  h = 0;
  atmcon = 0;
  for (i = 0; i < 3*crd->natom; i++) {
    atmcon = i/3;
    if (tp->EPInserted == 1 && tp->OldAtomNum[atmcon] == -1) {
      continue;
    }
    fprintf(outp, "%8.3lf", dtmp[i]);
    h++;
    if (h == 10) {
      h = 0;
      fprintf(outp, "\n");
    }
  }
  if (h > 0) {
    fprintf(outp, "\n");
  }
  if (tp->ifbox > 0) {
    fprintf(outp, "%8.3lf%8.3lf%8.3lf\n", crd->gdim[0], crd->gdim[1],
            crd->gdim[2]);
  }
  fclose(outp);
}

//-----------------------------------------------------------------------------
// FillNetCDFArray: allocate and populate an array of doubles for NetCDF
//                  writing.  This routine checks the topology to see if
//                  extra points have been inserted into the topology,  
//                  and condenses the data (forces, velocities, or      
//                  positions) to correspond to atoms in the original   
//                  topology if necessary.                              
//                                                                      
// Arguments:                                                           
//  tp:     the topology                                                
//  adata:  data on the atoms (3 * (number of atoms currently in tp)    
//          elements, applicable to forces, velocities, or positions)   
//-----------------------------------------------------------------------------
static double* FillNetCDFArray(prmtop *tp, double* adata)
{
  int i, j, o3i;
  double* cdfdata;

  if (tp->EPInserted == 0) {
    cdfdata = CpyDVec(adata, 3*tp->natom);
  }
  else {
    cdfdata = (double*)malloc(3*tp->norigatom*sizeof(double));
    for (i = 0; i < tp->natom; i++) {
      if (tp->OldAtomNum[i] > -1) {
        o3i = 3*tp->OldAtomNum[i];
        for (j = 0; j < 3; j++) {
          cdfdata[o3i+j] = adata[3*i+j];
        }
      }
    }
  }

  return cdfdata;
}

//-----------------------------------------------------------------------------
// WriteCDF: write a set of coordinates to a trajectory file in NetCDF format.
//                                                                      
// Arguments:                                                           
//   CG:   the cell grid (if NULL, this will change the behavior slightly, and
//         it is assumed that all coordinates or forces to be printed are
//         already in the coord data structure)
//   crd:  the coordinates structure                                    
//   cvf:  flag to direct writing of coordinates (1), velocities (2), or
//         forces (3)                                                   
//   tj:   trajectory control data                                      
//   Acdf: pointer to NetCDF file management struct                     
//   tp:   the topology                                                 
//   n:    write trajectory data for the nth system                     
//-----------------------------------------------------------------------------
void WriteCDF(cellgrid *CG, coord *crd, int cvf, trajcon *tj, cdftrj *Acdf,
              prmtop *tp, int n)
{
  int i, ntwcon, nsysatom, hasbox;
  double savegdim[3];
  double* cdfdata;
  char fname[MAXNAME];
  char *ctmp, *csuf;

  // Merge coordinates
#ifdef MPI
  if (CG != NULL) {
    CoordinateReduction(CG, crd, tp);
  }
#endif

  // Bail out if this is not the master    
  // process of this cell grid communicator
  if (CG != NULL && CG->tid != 0) {
    return;
  }

  // Test box dimensions
  hasbox = (crd->gdim[0] > 0.001 && crd->gdim[1] > 0.001 &&
	    crd->gdim[2] > 0.001);
  
  // Set a pointer to avoid conditionals in the loop
  if (cvf == 1) {
    cdfdata = FillNetCDFArray(tp, crd->loc);
    ctmp = tj->trjbase.map[n];
    csuf = tj->trjsuff.map[n];
    ntwcon = tj->ntwx;
  }
  else if (cvf == 2) {
    cdfdata = FillNetCDFArray(tp, crd->vel);
    ctmp = tj->velbase.map[n];
    csuf = tj->velsuff.map[n];
    ntwcon = tj->ntwv;
  }
  else if (cvf == 3) {
    if (CG != NULL) {
      MapCellForcesToAtoms(CG, crd);
    }
    cdfdata = FillNetCDFArray(tp, crd->frc);
    ctmp = tj->frcbase.map[n];
    csuf = tj->frcsuff.map[n];
    ntwcon = tj->ntwf;
  }

  // The number of atoms contained in each file
  nsysatom = (tp->EPInserted == 0) ? crd->natom : tp->norigatom;

  // Open the output file if needed
  SpliceFileName(tj, ctmp, csuf, fname, 0);
  if (OpenNewOutput(tj, ntwcon) == 1) {
    if (netcdfCreate(Acdf, fname, nsysatom, 1) == 1) {
      printf("WriteCDF >> Error.  Unable to open NetCDF file %s.\n", fname);
      exit(1);      
    }
  }

  // Write to file and increment frame
  if (tj->currstep > 0) {
    if (hasbox == 1) {
      for (i = 3; i < 6; i++) {
	savegdim[i-3] = crd->gdim[i];
	crd->gdim[i] *= 180.0 / PI;
      }
      if (netcdfWriteNextFrame(Acdf, cdfdata, crd->gdim) == 1) {
        printf("WriteCDF >> Error.  Unable to write to NetCDF file %s.\n",
	       fname);
        exit(1);
      }
      for (i = 3; i < 6; i++) {
	crd->gdim[i] = savegdim[i-3];
      }
    }
    else {
      if (netcdfWriteNextFrame(Acdf, cdfdata, NULL) == 1) {
        printf("WriteCDF >> Error.  Unable to write to NetCDF file %s.\n",
	       fname);
        exit(1);
      }
    }
  }

  // Close the output file if this is the end
  if (tj->currstep == tj->nstep || (tj->nfistep > 0 && tj->currstep > 0 &&
                                    tj->currstep == tj->currfi*tj->nfistep)) {
    SpliceFileName(tj, ctmp, csuf, fname, 0);
    if (netcdfClose(Acdf) == 1) {
      printf("WriteCDF >> Error.  Unable to close NetCDF file %s.\n", fname);
      exit(1);
    }
  }

  // Free allocated memory
  free(cdfdata);
}

//-----------------------------------------------------------------------------
// ReprintInputFile: reprint the input file into some developing output file.
//                   This routine will shorten the input file if requested by
//                   eliminating all but the first two and last two lines that
//                   contain particular, repetitive keywords.
//
// Arguments:
//   tj:        trajectory control data, contains a verbatim copy of the input
//   key1:      repetitive keyword to omit from output
//   key2:      alias of the repetitive keyword to omit
//   printkey:  flag to have all the keys printed regardless of repetition
//   outp:      the output file being written
//-----------------------------------------------------------------------------
void ReprintInputFile(trajcon *tj, char* key1, char* key2, int printkey,
                      FILE *outp)
{
  int i, j, nkeyline;
  cmat lnwords;

  HorizontalRule(outp, 0);
  fprintf(outp, "\nINPUT LINE TEXT:\n\n");
  PrintParagraph(tj->inpline, 79, NULL, outp);
  fprintf(outp, "\nINPUT FILE TEXT:\n\n");
  j = 0;
  nkeyline = 0;
  if (printkey == 0) {
    for (i = 0; i < tj->inptext.row; i++) {
      lnwords = ParseWords(tj->inptext.map[i]);
      if (strcmp(lnwords.map[0], "sys") == 0 ||
          strcmp(lnwords.map[0], "System") == 0) {
        nkeyline++;
      }
      DestroyCmat(&lnwords);
    }
  }
  for (i = 0; i < tj->inptext.row; i++) {

    // Omit most lines containing data point specifications,
    // unless this information is specifically requested to
    // be repeated in the output.
    if (printkey == 0 && nkeyline > 5) {
      lnwords = ParseWords(tj->inptext.map[i]);
      if (lnwords.row > 1 && (strcmp(lnwords.map[0], "sys") == 0 ||
                              strcmp(lnwords.map[0], "System") == 0)) {
        j++;
      }
      DestroyCmat(&lnwords);
      if (j == 3 || j == nkeyline-2) {
        fprintf(outp, "  ...\n");
        continue;
      }
      else if (j > 2 && j < nkeyline-2) {
        continue;
      }
    }
    fprintf(outp, "%s", tj->inptext.map[i]);
  }
  fprintf(outp, "\n");
  HorizontalRule(outp, 1);
}

//-----------------------------------------------------------------------------
// DestroyTrajCon: free all memory associated with a trajectory control data
//                 structure.                                      
//-----------------------------------------------------------------------------
void DestroyTrajCon(trajcon *tj)
{
  int i;

  DestroyCmat(&tj->inptext);
  DestroyCmat(&tj->ipcname);
  DestroyCmat(&tj->rstbase);
  DestroyCmat(&tj->trjbase);
  DestroyCmat(&tj->velbase);
  DestroyCmat(&tj->frcbase);
  DestroyCmat(&tj->rstsuff);
  DestroyCmat(&tj->trjsuff);
  DestroyCmat(&tj->velsuff);
  DestroyCmat(&tj->frcsuff);
  DestroyImat(&tj->SystemCPUs);
  for (i = 0; i < tj->nCPUcluster; i++) {
    free(tj->CPUcluster[i].atoms);
  }
  free(tj->CPUcluster);
  free(tj->inpline);
  free(tj->MySystemDomain);
#ifdef MPI
  for (i = 0; i < tj->nsys; i++) {
    MPI_Comm_free(&tj->SysComm[i]);
  }
  free(tj->SysComm);
#endif
  free(tj->Leash.GridFile);
  free(tj->Leash.GridDefsFile);
  free(tj->Leash.BellyMask);
  free(tj->Leash.FrozenMask);
}

//-----------------------------------------------------------------------------
// PositionCompromise: make a compromise between the positions of two atoms
//                     based on the TI mixing parameter.          
//                                                                      
// Arguments:                                                           
//   atm[A,B]:    atoms A and B                                         
//   mx[A,B]:     mixing parameter weights for systems A and B          
//-----------------------------------------------------------------------------
static double PositionCompromise(atomc *atmA, atomc *atmB, double mxA,
                                 double mxB)
{
  double dx, dy, dz;

  dx = atmB->loc[0] - atmA->loc[0];
  dy = atmB->loc[1] - atmA->loc[1];
  dz = atmB->loc[2] - atmA->loc[2];
  atmA->loc[0] += dx*mxB;
  atmB->loc[0] -= dx*mxA;
  atmA->loc[1] += dy*mxB;
  atmB->loc[1] -= dy*mxA;
  atmA->loc[2] += dz*mxB;
  atmB->loc[2] -= dz*mxA;

  return dx*dx + dy*dy + dz*dz;
}

//-----------------------------------------------------------------------------
// SynchronizeCoordinates: this function performs a periodic realignment of
//                         coordinates for corresponding atoms during TI
//                         calculations.                             
//                                                                      
// Arguments:                                                           
//   CG:      array containing the pair of cell grids                   
//   tj:      trajectory control information                            
//-----------------------------------------------------------------------------
void SynchronizeCoordinates(cellgrid* CG, trajcon *tj)
{
  int i, j, acon, bcon;
  double mxA, mxB, acdr;
  cell *cA, *cB;
  atomc *atmA, *atmB;
  prmtop *tpA;
  prmcorr *prc;
  char fname[MAXNAME];
  FILE *outp;

  // Determine the mixing parameters
  mxA = 1.0;
  for (i = 0; i < tj->mxorder; i++) {
    mxA *= (1.0 - tj->lambda);
  }
  mxB = 1.0 - mxA;

  // Pointers to topologies
  prc = &tj->prc;
  tpA = prc->tpA;

  const int ncell = CG->ng[0]*CG->ng[1]*CG->ng[2];
  acdr = 0.0;
  for (i = 0; i < ncell; i++) {
    cA = &CG[0].data[i];
    cB = &CG[1].data[i];
    atmA = cA->data;
    atmB = cB->data;
    if (prc->relate == 0) {

      // Simplest case: the two systems contain only
      // matching atoms in exactly the same order.  
      for (j = 0; j < cA->nr[0]; j++) {
        if (tpA->Masses[atmA[j].id] > 1.0e-8) {
          acdr += PositionCompromise(&atmA[j], &atmB[j], mxA, mxB);
        }
      }
    }
    else if (prc->relate == 1) {

      // More complex case: the two systems have some
      // unique atoms, but otherwise the same order  
      // of common atoms.                            
      acon = 0;
      bcon = 0;
      while (acon < cA->nr[0] && bcon < cB->nr[0]) {
        if (prc->matchA[atmA[acon].id] < 0) {
          acon++;
          continue;
        }
        if (prc->matchB[atmB[bcon].id] < 0) {
          bcon++;
          continue;
        }
        if (tpA->Masses[atmA[acon].id] > 1.0e-8) {
          acdr += PositionCompromise(&atmA[acon], &atmB[bcon], mxA, mxB);
        }
        acon++;
        bcon++;
      }
    }
    else if (prc->relate == 2) {

      // Horrible case: no orderly correspondence
      for (j = 0; j < cA->nr[0]; j++) {
        acon = atmA[j].id;
        if (prc->matchA[acon] < 0 || tpA->Masses[acon] < 1.0e-8) {
          continue;
        }
        bcon = cB->GPSptr[prc->matchA[acon]];
        acdr += PositionCompromise(&atmA[j], &atmB[bcon], mxA, mxB);
      }
    }
  }

  // Write to output
  if (tj->ntpr > 0 && CG->tid == 0) {
    acdr = sqrt(acdr / prc->comAB);
    SpliceFileName(tj, tj->outbase, tj->outsuff, fname, 1);
    outp = fopen(fname, "a");
    fprintf(outp, "%% Synchronized coordinates at step %9lld: RMS deviation "
            "%12.5e\n", tj->currstep, acdr);
    fclose(outp);
  }
}

//-----------------------------------------------------------------------------
// ReadListOfDoubles: read a list of double-precision real numbers, given in
//                    free format but for comment lines denoted by '%' or '#'.
//                    Returns the result as a matrix so that the number of
//                    elements can be tracked.
//
// Arguments:
//   fname:     name of the file to read
//-----------------------------------------------------------------------------
dmat ReadListOfDoubles(char* fname)
{
  int i, j, nval;
  FILE *finp;
  dmat T;
  cmat fmem, L;

  // Open the file and read it into memory
  fmem = Ascii2Mem(fname, 512, 8, "Routine called from ReadListOfDoubles.");
  nval = 0;
  for (i = 0; i < fmem.row; i++) {
    RemoveComments(fmem.map[i]);
    L = ParseWords(fmem.map[i]);
    for (j = 0; j < L.row; j++) {
      if (WordIsNumber(L.map[j]) == 1) {
        nval++;
      }
    }
    DestroyCmat(&L);
  }
  T = CreateDmat(1, nval, 0);
  nval = 0;
  for (i = 0; i < fmem.row; i++) {
    L = ParseWords(fmem.map[i]);
    for (j = 0; j < L.row; j++) {
      if (WordIsNumber(L.map[j]) == 1) {
        T.map[0][nval] = atof(L.map[j]);
        nval++;
      }
    }
    DestroyCmat(&L);
  }
  DestroyCmat(&fmem);

  return T;
}
