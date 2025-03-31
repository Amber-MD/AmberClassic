#include <time.h>
#include <stdio.h>
#include "Timings.h"
#include "Matrix.h"
#include "mdgxVector.h"

#include "pmeRecipDS.h"

//-----------------------------------------------------------------------------
// InitExecon: initialize a structure to keep timing information during an
//             mdgx run.                                             
//-----------------------------------------------------------------------------
void InitExecon(execon *tm)
{
  tm->bonds = 0.0;
  tm->cellcomm = 0.0;
  tm->nbInt = 0.0;
  tm->nbBsp = 0.0;
  tm->nbPtM = 0.0;
  tm->nbFFT = 0.0;
  tm->nbCnv = 0.0;
  tm->nbMtP = 0.0;
  tm->nbMtM = 0.0;
  tm->nbRecAll = 0.0;
  tm->Setup = 0.0;
  tm->Integ = 0.0;
  tm->Write = 0.0;
  tm->Constraints = 0.0;
  tm->Thermostat = 0.0;
  tm->Barostat = 0.0;
#ifdef MPI
  tm->mpiMeshPack = 0.0;
  tm->mpiMeshPullWait = 0.0;
  tm->mpiMeshPushWait = 0.0;
#endif
  gettimeofday(&tm->t0, NULL);
}

//-----------------------------------------------------------------------------
// mdgxStartTimer: record the time just before beginning some piece of code.
//                                                                      
// Arguments:                                                           
//   tm:     the execon struct in which the timings information is kept 
//-----------------------------------------------------------------------------
void mdgxStartTimer(execon *tm)
{
  gettimeofday(&tm->tti, NULL);
}

//-----------------------------------------------------------------------------
// mdgxStopTimer: record the time just after finishing some piece of code.
//                Before exiting, this function automatically sets the initial
//                time field to the final time field, as if mdgxStartTime were
//                implicitly called, so that multiple calls to mdgxStopTimer()
//                can be used to record a number of successive time intervals
//                without redundant calls to mdgxStartTimer().
//                                                                      
// Arguments:                                                           
//   tm:     the execon struct in which the timings information is kept 
//-----------------------------------------------------------------------------
double mdgxStopTimer(execon *tm)
{
  double dt;

  gettimeofday(&tm->ttf, NULL);
  dt = tm->ttf.tv_sec - tm->tti.tv_sec +
    (1.0e-6)*(tm->ttf.tv_usec - tm->tti.tv_usec);
  tm->tti = tm->ttf;

  return dt;
}

//-----------------------------------------------------------------------------
// GatherTimingData: gather timing data from all processes.  The data is
//                   returned as a double-precision matrix that mirrors many
//                   fields of the execon struct.                  
//                                                                      
// Arguments:                                                           
//   tm:     the execon struct in which the timings information is kept 
//   CG:     the cell grid (for communicator to count number of threads)
//   rcinp:  reciprocal space control data (to indicate MLE is active)  
//-----------------------------------------------------------------------------
dmat GatherTimingData(execon *tm, cellgrid *CG, reccon *rcinp)
{
  int nthr;
  dmat alltime;

#ifdef MPI
  MPI_Comm_size(CG->dspcomm, &nthr);
#else
  nthr = 1;
#endif
  alltime = CreateDmat(20, nthr, 0);
#ifdef MPI
  MPI_Gather(&tm->bonds, 1, MPI_DOUBLE, alltime.map[0], 1, MPI_DOUBLE, 0,
	     CG->dspcomm);
  MPI_Gather(&tm->cellcomm, 1, MPI_DOUBLE, alltime.map[1], 1, MPI_DOUBLE, 0,
	     CG->dspcomm);
  MPI_Gather(&tm->nbInt, 1, MPI_DOUBLE, alltime.map[2], 1, MPI_DOUBLE, 0,
	     CG->dspcomm);
  MPI_Gather(&tm->nbBsp, 1, MPI_DOUBLE, alltime.map[3], 1, MPI_DOUBLE, 0,
	     CG->dspcomm);
  MPI_Gather(&tm->nbPtM, 1, MPI_DOUBLE, alltime.map[4], 1, MPI_DOUBLE, 0,
	     CG->dspcomm);
  MPI_Gather(&tm->nbFFT, 1, MPI_DOUBLE, alltime.map[5], 1, MPI_DOUBLE, 0,
	     CG->dspcomm);
  MPI_Gather(&tm->nbCnv, 1, MPI_DOUBLE, alltime.map[6], 1, MPI_DOUBLE, 0,
	     CG->dspcomm);
  MPI_Gather(&tm->nbMtP, 1, MPI_DOUBLE, alltime.map[7], 1, MPI_DOUBLE, 0,
	     CG->dspcomm);
  if (rcinp->nlev > 1) {
    MPI_Gather(&tm->nbMtM, 1, MPI_DOUBLE, alltime.map[8], 1, MPI_DOUBLE, 0,
	       CG->dspcomm);
  }
  MPI_Gather(&tm->nbRecAll, 1, MPI_DOUBLE, alltime.map[9], 1, MPI_DOUBLE, 0,
	     CG->dspcomm);
  MPI_Gather(&tm->Setup, 1, MPI_DOUBLE, alltime.map[10], 1, MPI_DOUBLE, 0,
	     CG->dspcomm);
  MPI_Gather(&tm->Integ, 1, MPI_DOUBLE, alltime.map[11], 1, MPI_DOUBLE, 0,
	     CG->dspcomm);
  MPI_Gather(&tm->Write, 1, MPI_DOUBLE, alltime.map[12], 1, MPI_DOUBLE, 0,
             CG->dspcomm);
  MPI_Gather(&tm->Constraints, 1, MPI_DOUBLE, alltime.map[13], 1, MPI_DOUBLE,
	     0, CG->dspcomm);
  MPI_Gather(&tm->Barostat, 1, MPI_DOUBLE, alltime.map[14], 1, MPI_DOUBLE, 0,
             CG->dspcomm);
  MPI_Gather(&tm->Thermostat, 1, MPI_DOUBLE, alltime.map[15], 1, MPI_DOUBLE, 0,
             CG->dspcomm);
  MPI_Gather(&tm->mpiMeshPack, 1, MPI_DOUBLE, alltime.map[16], 1, MPI_DOUBLE,
	     0, CG->dspcomm);
  MPI_Gather(&tm->mpiMeshPullWait, 1, MPI_DOUBLE, alltime.map[17], 1,
	     MPI_DOUBLE, 0, CG->dspcomm);
  MPI_Gather(&tm->mpiMeshPushWait, 1, MPI_DOUBLE, alltime.map[18], 1,
	     MPI_DOUBLE, 0, CG->dspcomm);
#else
  alltime.map[0][0]  = tm->bonds;
  alltime.map[1][0]  = tm->cellcomm;
  alltime.map[2][0]  = tm->nbInt;
  alltime.map[3][0]  = tm->nbBsp;
  alltime.map[4][0]  = tm->nbPtM;
  alltime.map[5][0]  = tm->nbFFT;
  alltime.map[6][0]  = tm->nbCnv;
  alltime.map[7][0]  = tm->nbMtP;
  if (rcinp->nlev > 1) {
    alltime.map[8][0]  = tm->nbMtM;
  }
  alltime.map[9][0]  = tm->nbRecAll;
  alltime.map[10][0] = tm->Setup;
  alltime.map[11][0] = tm->Integ;
  alltime.map[12][0] = tm->Write;
  alltime.map[13][0] = tm->Constraints;
  alltime.map[14][0] = tm->Barostat;
  alltime.map[15][0] = tm->Thermostat;
#endif

  return alltime;
}

//-----------------------------------------------------------------------------
// PrintTimeLine: print timing data for one category of mdgx execution. 
//                                                                      
// Arguments:                                                           
//   alltime    : vector of time intervals
//   tt         : total time
//   ncpu       : the number of CPUs at work on this problem
//   desc       : description to be printed along with this output
//   outp       : pointer to the file being written
//-----------------------------------------------------------------------------
static void PrintTimeLine(double *alltime, double tt, int ncpu, char* desc,
			  FILE *outp)
{
  int i, mincpu, maxcpu;
  double avet, mint, maxt;
  
  avet = DAverage(alltime, ncpu);
  fprintf(outp, " %-22s  %9.2lf  %6.2lf", desc, avet, 100.0*avet/tt);
  if (ncpu > 1) {
    mint = 2.0*avet;
    maxt = 0.5*avet;
    for (i = 0; i < ncpu; i++) {
      if (alltime[i] < mint) {
	mint = alltime[i];
	mincpu = i;
      }
      if (alltime[i] > maxt) {
	maxt = alltime[i];
	maxcpu = i;
      }
    }
    fprintf(outp, "  %9.2lf %4d %9.2lf %4d\n", mint, mincpu, maxt, maxcpu);
  }
  else {
    fprintf(outp, "\n");
  }
}

//-----------------------------------------------------------------------------
// PrintTimingData: this function prints out all data related to timings
//                  accumulated over the course of the run.             
//                                                                      
// Arguments:                                                           
//   rcinp:    reciprocal space control data (to indicate MLE usage)    
//   alltime:  matrix of timing data for all threads (used in MPI mode) 
//   outp:     the output file                                          
//-----------------------------------------------------------------------------
void PrintTimingData(reccon *rcinp, dmat *alltime, FILE *outp)
{
  int i, ncpu;
  double tt;
  
  // Compute total time spent in this run
  ncpu = alltime->col;
  tt = 0.0;
  for (i = 0; i < alltime->row; i++) {
    tt += alltime->map[i][0];
  }
  fprintf(outp, " Timings are expressed as a percentage of total tracked "
	  "CPU execution time.\n\n"
	  " Segment                 Time / Percentage");
  if (ncpu > 1) {
    fprintf(outp, "   Min Time  CPU  Max Time  CPU");
  }
  fprintf(outp, "\n ----------------------  ---------  ------");
  if (ncpu > 1) {
    fprintf(outp, "  --------- ---- --------- ----");
  }
  fprintf(outp, "\n");
  PrintTimeLine(alltime->map[10], tt, ncpu, "Setup", outp);
  PrintTimeLine(alltime->map[0], tt, ncpu, "Bonded Interactions", outp);
  PrintTimeLine(alltime->map[1], tt, ncpu, "Cell Communication", outp);
  PrintTimeLine(alltime->map[2], tt, ncpu, "Nonbonded Interactions", outp);
  PrintTimeLine(alltime->map[3], tt, ncpu, "B-Spline Computation", outp);
  PrintTimeLine(alltime->map[4], tt, ncpu, "Particle -> Mesh", outp);
#ifdef MPI
  if (ncpu > 1) {
    PrintTimeLine(alltime->map[17], tt, ncpu, "Mesh Transmission I", outp);
    PrintTimeLine(alltime->map[18], tt, ncpu, "Mesh Transmission II", outp);
    PrintTimeLine(alltime->map[16], tt, ncpu, "Mesh Packaging", outp);
  }
#endif
  PrintTimeLine(alltime->map[6], tt, ncpu, "Convolution Kernel", outp);
  PrintTimeLine(alltime->map[5], tt, ncpu, "Fast Fourier Transform", outp);
  if (rcinp->nlev > 1) {
    PrintTimeLine(alltime->map[8], tt, ncpu, "Mesh -> Mesh", outp);
  }
  PrintTimeLine(alltime->map[7], tt, ncpu, "Mesh -> Particle", outp);
  PrintTimeLine(alltime->map[11], tt, ncpu, "Integration", outp);
  PrintTimeLine(alltime->map[13], tt, ncpu, "Constraints", outp);
  PrintTimeLine(alltime->map[15], tt, ncpu, "Thermostat", outp);
  PrintTimeLine(alltime->map[14], tt, ncpu, "Barostat", outp);
  PrintTimeLine(alltime->map[12], tt, ncpu, "Output Printing", outp);
  fprintf(outp, " Total CPU Time          %9.2lf  %6.2lf\n", tt,
	  100.0);
}
