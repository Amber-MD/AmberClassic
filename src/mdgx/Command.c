#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "Constants.h"
#include "Parse.h"
#include "Macros.h"
#include "mdgxVector.h"
#include "Matrix.h"
#include "pmeDirect.h"
#include "Random.h"
#include "Command.h"
#include "ConfigSamp.h"
#include "LoopBuilder.h"
#include "Peptide.h"

//-----------------------------------------------------------------------------
// PreExistingFileNames: for certain output trajectories, it may be possible
//                       to specify multiple file base names even while the
//                       first may be given by a default.  This routine
//                       automates checks for the existence of all such files.
//                                                                       
// Arguments:                                                            
//   C:       matrix of pre-specified file names                         
//   def0:    the default for the first file name                        
//   fspec:   integer array filled with ones if file names have been     
//            specified or changed from the defaults, or zeros otherwise 
//-----------------------------------------------------------------------------
static void PreExistingFileNames(cmat *C, char* def0, int* fspec)
{
  int i;

  // Test the existence of the first file 
  fspec[0] = (strcmp(C->map[0], def0) == 0) ? 0 : 1;
  for (i = 1; i < C->row; i++) {
    fspec[i] = (C->map[i][0] == '\0') ? 0 : 1;
  }
  for (i = C->row; i < MAXSYS; i++) {
    fspec[i] = 0;
  }
}

//-----------------------------------------------------------------------------
// GetDataFileNames: this function reads in the names of all necessary data
//                   (input and output) for a molecular dynamics run.  This
//                   routine also seeks "base" names and the corresponding
//                   suffixes, so mdgx can write numerous output files with
//                   the same base and suffix but different numbers (i.e.
//                   trj101.crd, where the base is "trj", the suffix is 
//                   ".crd", and the number is 101).  Note that the
//                   initialization of file names is done in the
//                   CommandLineControl function.
//
// Arguments:                                                            
//   tp:      the topology structure(s) (each topology stores its own source   
//            file name); tp is an array with two elements for up to two
//            systems
//   tj:      the trajectory structure (stores the restart file name, as well
//            as many other output file names)                      
//   inp:     the input file                                             
//-----------------------------------------------------------------------------
static void GetDataFileNames(prmtop* tp, trajcon *tj, FILE *inp)
{
  int i, collect, EXoutbase, EXforcedump, EXparmfile, EXfmodfile, nfmodfile; 
  int* EXinpcrd;
  int* EXfrcbase;
  int* EXrstbase;
  int* EXtrjbase;
  int* EXvelbase;
  int* EXfrcsuff;
  int* EXrstsuff;
  int* EXtrjsuff;
  int* EXvelsuff;
  int EXtp[2], EXeprules[2];
  char line[MAXLINE], searchtag[32], searchalias[32];
  cmat lwords;

  // Allocate memory for pre-existing file name/suffix flags 
  EXinpcrd = (int*)calloc(MAXSYS, sizeof(int));
  EXfrcbase = (int*)calloc(MAXSYS, sizeof(int));
  EXrstbase = (int*)calloc(MAXSYS, sizeof(int));
  EXtrjbase = (int*)calloc(MAXSYS, sizeof(int));
  EXvelbase = (int*)calloc(MAXSYS, sizeof(int));
  EXfrcsuff = (int*)calloc(MAXSYS, sizeof(int));
  EXrstsuff = (int*)calloc(MAXSYS, sizeof(int));
  EXtrjsuff = (int*)calloc(MAXSYS, sizeof(int));
  EXvelsuff = (int*)calloc(MAXSYS, sizeof(int));

  // First, verify that certain variables which could have    
  // been given on the command line are not already specified 
  EXoutbase = (strcmp(tj->outbase, "mdout") != 0) ? 1 : 0;
  EXforcedump = (strcmp(tj->dumpname, "forcedump.dat") != 0) ? 1 : 0;
  EXparmfile = (tj->parmfile[0] == '\0') ? 0 : 1;
  EXfmodfile = 0;
  while (EXfmodfile < tj->fmodfile.row &&
	 tj->fmodfile.map[EXfmodfile][0] != '\0') {
    EXfmodfile++;
  }
  nfmodfile = 0;
  EXtp[0] = (strcmp(tp[0].source, "prmtop") != 0) ? 1 : 0;
  EXtp[1] = (tp[1].source[0] == '\0') ? 0 : 1;
  for (i = 0; i < 2; i++) {
    EXeprules[i] = (tp[i].eprulesource[0] != '\0') ? 1 : 0;
  }
  PreExistingFileNames(&tj->ipcname, "inpcrd", EXinpcrd);
  PreExistingFileNames(&tj->trjbase, "mdcrd", EXtrjbase);
  PreExistingFileNames(&tj->velbase, "mdvel", EXvelbase);
  PreExistingFileNames(&tj->frcbase, "mdfrc", EXfrcbase);
  PreExistingFileNames(&tj->rstbase, "restrt", EXrstbase);

  // Search through the &files namelist 
  collect = AdvanceToSegment(inp, "files", 1);
  while (collect == 1) {
    collect = ReadNamelistLine(line, &lwords, "GetDataFileNames", inp);
    if (collect == 0) {
      continue;
    }

    // List of single string arguments to search for 
    if (EXtp[0] == 0) {
      SeekString(lwords, tp[0].source, "Topology", "-p");
    }
    if (EXeprules[0] == 0) {
      SeekString(lwords, tp[0].eprulesource, "EPRules", "-xpt");
    }
    for (i = 0; i < 2; i++) {
      sprintf(searchtag, "Topology%d", i+1);
      sprintf(searchalias, "-p%d", i+1);
      if (EXtp[i] == 0) {
	SeekString(lwords, tp[i].source, searchtag, searchalias);
      }
      sprintf(searchtag, "EPRules%d", i+1);
      sprintf(searchalias, "-xpt%d", i+1);
      if (EXeprules[i] == 0) {
	SeekString(lwords, tp[i].eprulesource, searchtag, searchalias);
      }
    }
    if (EXoutbase == 0) {
      SeekString(lwords, tj->outbase, "OutputBase", "-o");
    }
    if (EXforcedump == 0) {
      SeekString(lwords, tj->dumpname, "ForceDump", "-d");
    }
    if (EXparmfile == 0) {
      SeekString(lwords, tj->parmfile, "ParmFile", "-parm");
    }
    if (EXfmodfile == 0) {
      SeekStringInc(lwords, &tj->fmodfile, "FrcmodFile", "-fmod", &nfmodfile);
    }
    if (EXinpcrd[0] == 0) {
      SeekString(lwords, tj->ipcname.map[0], "StartCoord", "-c");
    }
    if (EXtrjbase[0] == 0) {
      SeekString(lwords, tj->trjbase.map[0], "CrdTrajBase", "-x");
    }
    if (EXvelbase[0] == 0) {
      SeekString(lwords, tj->velbase.map[0], "VelTrajBase", "-v");
    }
    if (EXfrcbase[0] == 0) {
      SeekString(lwords, tj->frcbase.map[0], "FrcTrajBase", "-f");
    }
    if (EXrstbase[0] == 0) {
      SeekString(lwords, tj->rstbase.map[0], "Restart", "-r");
    }
    SeekString(lwords, tj->outsuff, "OutputSuff", "-osf");
    SeekString(lwords, tj->rsrptname, "ResReport", "-rrp");
    SeekString(lwords, tj->trjsuff.map[0], "CrdTrajSuff", "-xsf");
    SeekString(lwords, tj->velsuff.map[0], "VelTrajSuff", "-vsf");
    SeekString(lwords, tj->frcsuff.map[0], "FrcTrajSuff", "-fsf");
    SeekString(lwords, tj->rstsuff.map[0], "RestartSuff", "-rsf");

    // List of group string arguments to search for 
    tj->ipcname = SeekNString(lwords, &tj->ipcname, EXinpcrd, "StartCoord",
			      "-c");
    tj->trjbase = SeekNString(lwords, &tj->trjbase, EXtrjbase, "CrdTrajBase",
			      "-x");
    tj->velbase = SeekNString(lwords, &tj->velbase, EXvelbase, "VelTrajBase",
			      "-v");
    tj->frcbase = SeekNString(lwords, &tj->frcbase, EXfrcbase, "FrcTrajBase",
			      "-f");
    tj->rstbase = SeekNString(lwords, &tj->rstbase, EXrstbase, "Restart",
			      "-r");
    tj->trjsuff = SeekNString(lwords, &tj->trjsuff, EXtrjsuff, "CrdTrajSuff",
			      "-xsf");
    tj->velsuff = SeekNString(lwords, &tj->velsuff, EXvelsuff, "VelTrajSuff",
			      "-vsf");
    tj->frcsuff = SeekNString(lwords, &tj->frcsuff, EXfrcsuff, "FrcTrajSuff",
			      "-fsf");
    tj->rstsuff = SeekNString(lwords, &tj->rstsuff, EXrstsuff, "Restart",
			      "-rsf");

    // Free allocated memory 
    DestroyCmat(&lwords);
  }

  // Free allocated memory 
  free(EXinpcrd);
  free(EXrstbase);
  free(EXtrjbase);
  free(EXvelbase);
  free(EXfrcbase);
  free(EXrstsuff);
  free(EXtrjsuff);
  free(EXvelsuff);
  free(EXfrcsuff);
  if (EXfmodfile == 0 && nfmodfile > 0) {
    tj->fmodfile = ReallocCmat(&tj->fmodfile, nfmodfile, MAXNAME);
  }
}

//-----------------------------------------------------------------------------
// CountOutputStreams: this function counts the number of output streams 
//                     specified in a particular character matrix.  Also 
//                     checks to make sure that there are no gaps in the 
//                     list of file names.                               
//                                                                       
// Arguments:                                                            
//   T:       the character matrix containing file (base) names          
//-----------------------------------------------------------------------------
static int CountOutputStreams(cmat *T)
{
  int i, nout, nhigh;

  if (T->row == 0) {
    return 0;
  }

  nout = 0;
  for (i = 0; i < T->row; i++) {
    if (T->map[i][0] != '\0') {
      nout++;
      nhigh = i;
    }
  }
  if (nout != nhigh+1) {
    printf("CountOutputStreams >> Error.  List of output names invalid:\n");
    for (i = 0; i < nhigh+1; i++) {
      printf("CountOutputStreams >>  %s\n", T->map[i]);
    }
    exit(1);
  }

  return nout;
}

//-----------------------------------------------------------------------------
// ExpandSuffixes: this function fills in blank suffix names when more   
//                 base file names have been specified.  If a particular 
//                 suffix name is specified, all blank suffixes after it 
//                 will be assigned the same string until a non-blank    
//                 suffix is encountered, at which point that becomes    
//                 the new replacement for blank suffixes.               
//                                                                       
// Arguments:                                                            
//   S:     the character matrix of suffixes                             
//   B:     the character matrix of base names                           
//-----------------------------------------------------------------------------
static void ExpandSuffixes(cmat *S, cmat *B)
{
  int i;
  char *ctmp;

  // Expand the suffixes to be as numerous as the bases 
  if (S->row < B->row) {
    *S = ReallocCmat(S, B->row, S->col);
  }

  // The first suffix must be specified or else 
  // we just wipe all the suffixes and bail out 
  if (S->map[0][0] == '\0') {
    for (i = 0; i < S->row*S->col; i++) {
      S->data[i] = '\0';
    }
    return;
  }
  for (i = 0; i < S->row; i++) {
    if (S->map[i][0] != '\0') {
      ctmp = S->map[i];
    }
    else {
      strcpy(S->map[i], ctmp);
    }
  }
}

//-----------------------------------------------------------------------------
// SetOutputFrequency: this function ensures that the frequency of each  
//                     output stream is a factor of the total number of  
//                     steps in either the segment or full trajectory.   
//                                                                       
// Arguments:                                                            
//   tj:     trajectory control information                              
//   ival:   the output frequency (also returned, possibly modified)     
//   vname:  the name of the variable being examined (printed in any     
//           warnings)                                                   
//-----------------------------------------------------------------------------
static int SetOutputFrequency(trajcon *tj, int ival, char* vname)
{
  int i;
  long long int stepcount, ficount;

  // Check to make sure that the interval  
  // not greater than the number of steps. 
  if (ival > tj->nstep) {
    printf("GetCntrlNamelist >> Warning.  When the output frequency (%d) is\n"
	   "GetCntrlNamelist >> greater than the anticipated number of steps "
	   "(%lld)\nGetCntrlNamelist >> the output cannot be properly "
	   "averaged.\nGetCntrlNamelist >> Setting %s to %lld.\n", ival,
	   tj->nstep, vname, tj->nstep);
    return tj->nstep;
  }

  // Bail right out if the interval is just zero 
  if (ival == 0 || tj->nfistep == 0) {
    return ival;
  }

  // Make sure this interval is a clean one 
  stepcount = (tj->nfistep > 0) ? tj->nfistep : tj->nstep;
  if (stepcount % ival > 0) {
    printf("GetCntrlNamelist >> Warning.  When the run output is split into "
           "segments,\nGetCntrlNamelist >> all output frequencies must be "
	   "factors of nfistep.");
    if ((double)stepcount / (double)ival < 10.0) {
      for (i = ival+1; i <= stepcount; i++) {
	if (tj->nfistep % i == 0) {
	  ival = i;
	  break;
	}
      }
      printf("%s set to %d.\n", vname, ival);
    }
    else {
      if (tj->nfistep > 0) {
	ficount = tj->nstep / tj->nfistep;
      }
      stepcount = (ival+1) * (stepcount / ival);
      if (tj->nfistep > 0) {
	tj->nfistep = stepcount;
	tj->nstep = stepcount*ficount;
	printf("nfistep set to %lld,\nGetCntrlNamelist >> ", tj->nfistep);
	printf("nstep   set to %lld.\n", tj->nstep);
      }
      else {
	tj->nstep = stepcount;
	printf("nstep set to %lld.\n", stepcount);
      }
    }
  }

  return ival;
}

//-----------------------------------------------------------------------------
// GetCntrlNamelist: this function reads in run parameters, as would be found 
//                   in the &cntrl namelist in PMEMD or SANDER command input
//                   files.                                
//
// Arguments:                                                            
//   dcinp:   direct space command information                           
//   rcinp:   reciprocal space command information                       
//   tp:      the topology structure array                               
//   tj:      the trajectory structure (stores the restart file name, as 
//            well as many other output file names)                      
//   inp:     the input file                                             
//-----------------------------------------------------------------------------
static void GetCntrlNamelist(dircon *dcinp, prmtop* tp, trajcon *tj, FILE *inp)
{
  int i, collect, EXigseed, EXlambda, ntc;
  double EVcut, lampow, lampowm1, mfac, mxA, dmxA;
  double* pval;
  char line[MAXLINE];
  cmat lwords;

  // Default settings for dynamics / Hamiltonian 
  tp->lj14fac = 2.0;
  tp->elec14fac = 1.2;
  dcinp->Ecut = 8.0;
  dcinp->Vcut = 8.0;
  dcinp->MaxDens = 2.5;
  EVcut = -1.0;
  sprintf(tp->WaterName, "WAT ");
  tj->mode = 0;
  tj->starttime = -1.0;
  tj->dt = 0.001;
  tp->rattle = 0;
  tp->settle = 0;

  // Default settings for thermostat and barostat 
  ntc = 0;
  tj->Ttarget = -100.0;
  tj->Tinit = -100.0;
  tj->Ptarget = 1.0;
  tj->rattletol = 1.0e-6;
  tj->MaxRattleIter = 100;
  tj->RemoveMomentum = 0;
  tj->ioutfm = 0;
  tj->ntt = 0;
  tj->ntp = 0;
  tj->barostat = 1;
  tj->vrand = 1000;
  tj->BerendsenTCoupl = 0.4;
  tj->BerendsenPTime = 1.0;
  tj->BerendsenPCoupl = 44.6;
  tj->MCBarostatFac[0] = 2.0e-3;
  tj->MCBarostatFac[1] = -1.0;
  tj->MCBarostatFac[2] = -1.0;
  tj->MCBarostatFreq = 100;
  tj->npth.TauT = 1.0;
  tj->npth.TauP = 1.0;
  tj->lnth.gamma_ln = 0.0;

  // Default settings for trajectory control 
  tj->nstep = 1;
  tj->nfistep = 0;
  tj->currfi = 0;
  tj->ntwr = 0;
  tj->ntwx = 0;
  tj->ntwv = 0;
  tj->ntwf = 0;
  tj->ntpr = 0;
  tj->irest = 0;
  tj->topchk = 1;

  // Default settings for energy minimization 
  tj->EMinStep = 1.0;
  tj->EMinStep0 = 1.0;
  tj->EMinTol = 1.0e-4;

  // Default settings for thermodynamic integration 
  tj->TI = 0;
  tj->mxorder = 1;
  tj->nsynch = 1000;
  EXigseed = (tj->igseed == 72177) ? 0 : 1;
  EXlambda = (fabs(tj->lambda) < 1.0e-8) ? 0 : 1;

  // Allocate space for topology atom masks 
  for (i = 0; i < 2; i++) {
    tp[i].norattlemask = (char*)calloc(MAXNAME, sizeof(char));
    tp[i].rattlemask = (char*)calloc(MAXNAME, sizeof(char));
  }
  collect = AdvanceToSegment(inp, "cntrl", 1);
  while (collect == 1) {
    collect = ReadNamelistLine(line, &lwords, "GetCntrlNamelist", inp);
    if (collect == 0) {
      continue;
    }

    // Dynamics / Hamiltonian directives 
    SeekReal(lwords, &tp->lj14fac, "Vdw14Fac", "scnb");
    SeekReal(lwords, &tp->elec14fac, "Elec14Fac", "scee");
    SeekReal(lwords, &tj->rattletol, "RattleTol", "tol");
    SeekReal(lwords, &dcinp->Ecut, "ElecCut", "es_cutoff");
    SeekReal(lwords, &dcinp->Vcut, "VdwCut", "vdw_cutoff");
    SeekReal(lwords, &EVcut, "DirectCutoff", "cut");
    SeekReal(lwords, &tj->starttime, "StartTime", "t");
    SeekReal(lwords, &tj->dt, "TimeStep", "dt");
    SeekReal(lwords, &dcinp->MaxDens, "MaxDensity", "rho");
    SeekInt(lwords, &tj->irest, "RestartMD", "irest");
    SeekInt(lwords, &tp->rattle, "DoRATTLE", "rigidbond");
    SeekInt(lwords, &tp->settle, "DoSETTLE", "rigidwat");
    SeekInt(lwords, &tj->MaxRattleIter, "MaxRattleIter", "vlimit");
    SeekInt(lwords, &tj->RemoveMomentum, "ZeroMomentum", "nscm");
    SeekInt(lwords, &tj->mode, "RunMode", "imin");
    SeekInt(lwords, &ntc, "ShakeLevel", "ntc");
    if (EXigseed == 0) {
      SeekInt(lwords, &tj->igseed, "RandomSeed", "ig");
    }

    // Thermodynamic integration directives 
    SeekInt(lwords, &tj->TI, "RunTI", "icfe");
    SeekInt(lwords, &tj->mxorder, "MixOrder", "klambda");
    SeekInt(lwords, &tj->nsynch, "SynchTI", "nsynch");
    if (EXlambda == 0) {
      SeekReal(lwords, &tj->lambda, "MixFactor", "clambda");
    }

    // Thermostat and barostat directives 
    SeekReal(lwords, &tj->Ttarget, "Temperature", "temp0");
    SeekReal(lwords, &tj->Ptarget, "Pressure", "pres0");
    SeekReal(lwords, &tj->Tinit, "StartTemp", "tempi");
    SeekReal(lwords, &tj->BerendsenTCoupl, "BerendsenTC", "taup");
    SeekReal(lwords, &tj->BerendsenPCoupl, "BerendsenPC", "comp");
    SeekReal(lwords, &tj->npth.TauT, "HooverTC", "tauthv");
    SeekReal(lwords, &tj->npth.TauP, "HooverPC", "tauphv");
    SeekReal(lwords, &tj->lnth.gamma_ln, "LangevinFreq", "gamma_ln");
    SeekReal(lwords, &tj->MCBarostatFac[0], "MCBarostatPC", "mccomp");
    SeekReal(lwords, &tj->MCBarostatFac[0], "MCBarostatPCX", "mccompx");
    SeekReal(lwords, &tj->MCBarostatFac[1], "MCBarostatPCY", "mccompy");
    SeekReal(lwords, &tj->MCBarostatFac[2], "MCBarostatPCZ", "mccompz");
    SeekInt(lwords, &tj->ntt, "Thermostat", "ntt");
    SeekInt(lwords, &tj->ntp, "CoordRscl", "ntp");
    SeekInt(lwords, &tj->barostat, "Barostat", "barostat");
    SeekInt(lwords, &tj->vrand, "RandomReset", "vrand");
    SeekInt(lwords, &tj->MCBarostatFreq, "MCBarostatFrq", "mcbfrq");

    // Output management directives 
    SeekLLInt(lwords, &tj->nstep, "StepCount", "nstlim");
    SeekLLInt(lwords, &tj->nfistep, "FileStepCount", "nfistep");
    SeekInt(lwords, &tj->ntpr, "WriteDiagnostics", "ntpr");
    SeekInt(lwords, &tj->ntwx, "WriteCoordinates", "ntwx");
    SeekInt(lwords, &tj->ntwv, "WriteVelocities", "ntwv");
    SeekInt(lwords, &tj->ntwf, "WriteForces", "ntwf");
    SeekInt(lwords, &tj->ntwr, "WriteRestart", "ntwr");
    SeekInt(lwords, &tj->topchk, "TopologyCheck", "tchk");
    SeekInt(lwords, &tj->ioutfm, "CoordFormat", "ioutfm");
    SeekString(lwords, tp[0].WaterName, "WaterName1", "watnam");
    SeekString(lwords, tp[1].WaterName, "WaterName2", "watnam");
    SeekString(lwords, tp[0].rattlemask, "RattleMask1", "rattlemask");
    SeekString(lwords, tp[1].rattlemask, "RattleMask2", "rattlemask");
    SeekString(lwords, tp[0].norattlemask, "NoRattleMask1", "norattlemask");
    SeekString(lwords, tp[1].norattlemask, "NoRattleMask2", "norattlemask");

    // Free allocated memory 
    DestroyCmat(&lwords);
  }

  // Restate control data in internal units/conventions 
  if (EVcut > 0.0) {
    dcinp->Ecut = EVcut;
    dcinp->Vcut = EVcut;
  }
  tp->lj14fac = 1.0 - 1.0/tp->lj14fac;
  tp->elec14fac = 1.0 - 1.0/tp->elec14fac;
  dcinp->Mcut = MAX(dcinp->Ecut, dcinp->Vcut);
  dcinp->invMcut = 1.0/dcinp->Mcut;
  dcinp->invEcut = 1.0/dcinp->Ecut;
  tj->BerendsenPCoupl *= 1.0e-6;
  if (tj->Ttarget < 0.0) {
    tj->Ttarget = 298.0;
    if (tj->Tinit < 0.0) {
      tj->Tinit = 0.0;
    }
  }
  else {
    if (tj->Tinit < 0.0) {
      tj->Tinit = tj->Ttarget;
    }
  }
  if (tj->TI == 1) {
    pval = PascalTriangle(tj->mxorder+1);
    lampow = 1.0;
    lampowm1 = 1.0/tj->lambda;
    mfac = 1.0;
    mxA = 0.0;
    dmxA = 0.0;
    for (i = 0; i <= tj->mxorder; i++) {
      mxA += mfac*lampow*pval[i];
      dmxA += i*mfac*lampowm1*pval[i];
      mfac *= -1.0;
      lampow *= tj->lambda;
      lampowm1 *= tj->lambda;
    }
    tj->mxA = mxA;
    tj->mxB = 1.0 - mxA;
    tj->dmxA = dmxA;
    tj->dmxB = -dmxA;
  }
  tj->rattle = tp->rattle;
  if (ntc > 0) {

    // ntc settings can determine the SHAKE and RATTLE 
    // preferences, and ntc will take precedence.      
    if (ntc == 2) {
      if (tp->rattle == 0) {
	tp->rattle = 1;
      }
      if (tp->settle == 0) {
	tp->settle = 1;
      }
    }
    else if (ntc > 2) {
      printf("GetCntrlNamelist >> Error.  ntc > 2 is not supported.\n");
    }
  }

  // Check input: gamma_ln must for Langevin thermostat 
  if (tj->ntt == 3) {
    if (tj->lnth.gamma_ln < 1.0e-8) {
      printf("GetCntrlNamelist >> Error.  Langevin thermostat cannot be "
	     "called with a\nGetCntrlNamelist >> collision frequency of "
	     "%8.4lf / ps.\n", tj->lnth.gamma_ln);
      exit(1);
    }
  }

  // Check input: van-der Waals cutoff 
  // must be .GE. electrostatic cutoff 
  if (dcinp->Vcut < dcinp->Ecut) {
    printf("GetCntrlNamelist >> Error.  Vdw cutoff (%9.4lf) must be >= "
	   "electrostatic\nGetCntrlNamelist >> cutoff (%9.4lf)\n.",
	   dcinp->Vcut, dcinp->Ecut);
    exit(1);
  }

  // Check input: Monte-Carlo barostat is   
  // the only permissible option with ntp=2 
  if (tj->ntp == 2 && tj->barostat != 2) {
    printf("GetCntrlNamelist >> Error.  The Monte-Carlo barostat must be "
	   "used for\nGetCntrlNamelist >> anisotropic rescaling.\n");
    exit(1);
  }

  // Check input: certain RATTLE procedures are not allowed 
  if (tp->rattle == 2) {
    printf("GetCntrlNamelist >> Error.  Constraining all bonds in not "
	   "supported.\nGetCntrlNamelist >> Set DoRATTLE/ntc = 0 or 1.\n");
    exit(1);
  }

  // Check input: problems with file segmentation 
  if (tj->nfistep > 0) {
    if (tj->nstep % tj->nfistep != 0) {
      printf("GetCntrlNamelist >> Warning.  nstlim (mod) nfistep is nonzero.\n"
	     "GetCntrlNamelist >> nfistep / FileStepCount set to zero.\n");
      tj->nfistep = 0;
    }
    if (tj->ntwr > tj->nfistep || tj->ntwx > tj->nfistep ||
	tj->ntwv > tj->nfistep || tj->ntwf > tj->nfistep ||
	tj->ntpr > tj->nfistep) {
      printf("GetCntrlNamelist >> Warning.  When the run output is split into "
	     "segments,\nGetCntrlNamelist >> ntwr, ntpr, ntwx, ntwv, and "
	     "ntwf cannot exceed nfistep.\n");
      if (tj->ntwr > tj->nfistep) {
	tj->ntwr = tj->nfistep;
      }
      if (tj->ntwx > tj->nfistep) {
	tj->ntwx = tj->nfistep;
      }
      if (tj->ntwv > tj->nfistep) {
	tj->ntwv = tj->nfistep;
      }
      if (tj->ntwf > tj->nfistep) {
	tj->ntwf = tj->nfistep;
      }
      if (tj->ntpr > tj->nfistep) {
	tj->ntpr = tj->nfistep;
      }
    }
  }
  tj->ntpr = SetOutputFrequency(tj, tj->ntpr, "ntpr");
  tj->ntwx = SetOutputFrequency(tj, tj->ntwx, "ntwx");
  tj->ntwv = SetOutputFrequency(tj, tj->ntwv, "ntwv");
  tj->ntwf = SetOutputFrequency(tj, tj->ntwf, "ntwf");
  tj->ntwr = SetOutputFrequency(tj, tj->ntwr, "ntwr");

  // After reading the control namelist, we should know 
  // whether TI is active or not, and therefore whether 
  // we are dealing with two systems or just one.  We   
  // also know exactly how many systems have been       
  // specified, and how many topologies are available.  
  tj->nsys = tj->ipcname.row;
  if (tj->TI == 1 && tj->nsys != 2) {
    printf("GetCntrlNamelist >> Error.  Thermodynamic integration has been "
	   "activated.\nGetCntrlNamelist >> Two sets of initial coordinates "
	   "are required.\n");
    exit(1);
  }
  if (tj->TI == 1 && tp[1].source[0] == '\0') {
    printf("GetCntrlNamelist >> Error.  Thermodynamic integration has been "
	   "activated.\nGetCntrlNamelist >> Two topologies are required.\n");
    exit(1);
  }
  tj->ntop = (tj->TI == 1) ? 2 : 1;
  if (tj->ntop == 1) {
    free(tp[1].rattlemask);
    free(tp[1].norattlemask);
  }

  // Check the output stream counts 
  if (tj->ntwr > 0 && CountOutputStreams(&tj->rstbase) != tj->nsys) {
    printf("GetCntrlNamelist >> Error. %d restart file base names "
	   "were specified\nGetCntrlNamelist >> for %d systems.\n",
	   CountOutputStreams(&tj->rstbase), tj->nsys);
    exit(1);
  }
  if (tj->ntwx > 0 && CountOutputStreams(&tj->trjbase) != tj->nsys) {
    printf("GetCntrlNamelist >> Error. %d trajectory file base names "
	   "were specified\nGetCntrlNamelist >> for %d systems.\n",
	   CountOutputStreams(&tj->trjbase), tj->nsys);
    exit(1);
  }
  if (tj->ntwv > 0 && CountOutputStreams(&tj->velbase) != tj->nsys) {
    printf("GetCntrlNamelist >> Error. %d velocity file base names "
	   "were specified\nGetCntrlNamelist >> for %d systems.\n",
	   CountOutputStreams(&tj->velbase), tj->nsys);
    exit(1);
  }
  if (tj->ntwf > 0 && CountOutputStreams(&tj->frcbase) != tj->nsys) {
    printf("GetCntrlNamelist >> Error. %d force file base names "
	   "were specified\nGetCntrlNamelist >> for %d systems.\n",
	   CountOutputStreams(&tj->frcbase), tj->nsys);
    exit(1);
  }

  // Check for output file suffixes 
  ExpandSuffixes(&tj->rstsuff, &tj->rstbase);
  ExpandSuffixes(&tj->trjsuff, &tj->trjbase);
  ExpandSuffixes(&tj->velsuff, &tj->velbase);
  ExpandSuffixes(&tj->frcsuff, &tj->frcbase);
}

//-----------------------------------------------------------------------------
// GetTopolNamelist: this function reads in additional information about a
//                   topology that cannot be specified in the typical AMBER  
//                   topology format.  Essentially, information provided in
//                   the &topol namelist alerts MDGX that it's dealing with a
//                   non-standard topology, and to look for certain things
//                   when it actually does read in the topology file.
//
// Arguments:                                                            
//   tp:      the topology structure                                     
//   inp:     the input file (mdin)                                      
//-----------------------------------------------------------------------------
static void GetTopolNamelist(prmtop *tp, FILE *inp)
{
  int collect;
  char line[MAXLINE];
  cmat lwords;

  // Default settings 
  tp->ljbuck = 0;
  tp->qshape = 0;
  collect = AdvanceToSegment(inp, "cntrl", 1);
  while (collect == 1) {
    collect = ReadNamelistLine(line, &lwords, "GetTopolNamelist", inp);
    if (collect == 0) {
      continue;
    }

    // Information about the topology that cannot go in a 
    // standard AMBER prmtop file                         
    SeekInt(lwords, &tp->ljbuck, "vdWstyle", "ljstyle");
    SeekInt(lwords, &tp->qshape, "ChargeStyle", "qshape");

    // Free allocated memory 
    DestroyCmat(&lwords);
  }
}

//-----------------------------------------------------------------------------
// GetEwaldNamelist: this function reads in Ewald parameters, as would be
//                   found in the &ewald namelist in PMEMD or SANDER command
//                   input files.                                
//
// Arguments:                                                            
//   dcinp:   direct space command information                           
//   rcinp:   reciprocal space command information                       
//   inp:     the input file (mdin)                                      
//-----------------------------------------------------------------------------
static void GetEwaldNamelist(dircon *dcinp, reccon *rcinp, FILE *inp)
{
  int i, collect, genordr;
  char line[MAXLINE], qformstr[MAXNAME];
  cmat lwords;

  // Default settings for SPME 
  rcinp->ng = (int*)malloc(3*sizeof(int));
  rcinp->ng[0] = -1;
  rcinp->ng[1] = -1;
  rcinp->ng[2] = -1;
  rcinp->ordr[0] = 4;
  rcinp->ordr[1] = 4;
  rcinp->ordr[2] = 4;
  sprintf(qformstr, "POINT");
  dcinp->Dtol = 6.0e-6;
  dcinp->lkpspc = 0.0625;
  dcinp->ewcoeff = -1.0;
  dcinp->LRvdw = 1;
  rcinp->S = -1.0;
  genordr = -1;

  // Default settings for MLE 
  rcinp->ggordr = 8;
  rcinp->nlev = 1;
  rcinp->nslab = 1;
  rcinp->nstrip = 1;
  for (i = 0; i < 4; i++) {
    rcinp->cfac[i] = -1.0;
    rcinp->PadYZ[i] = -1;
  }
  rcinp->cfac[0] = 1.0;
  rcinp->PadYZ[0] = 1;

  // Scan the input file 
  collect = AdvanceToSegment(inp, "ewald", 1);
  while (collect == 1) {
    collect = ReadNamelistLine(line, &lwords, "GetEwaldNamelist", inp);
    if (collect == 0) {
      continue;
    }

    // PME directives 
    SeekReal(lwords, &dcinp->Dtol, "DSumTol", "dsum_tol");
    SeekReal(lwords, &dcinp->ewcoeff, "EwaldCof", "ew_coeff");
    SeekReal(lwords, &rcinp->S, "Sigma", "sigma");
    SeekReal(lwords, &dcinp->lkpspc, "SplnSpc", "eedtbdns");
    SeekInt(lwords, &rcinp->ng[0], "MeshDimX", "nfft1");
    SeekInt(lwords, &rcinp->ng[1], "MeshDimY", "nfft2");
    SeekInt(lwords, &rcinp->ng[2], "MeshDimZ", "nfft3");
    SeekInt(lwords, &rcinp->ordr[0], "OrderX", "ordr1");
    SeekInt(lwords, &rcinp->ordr[1], "OrderY", "ordr2");
    SeekInt(lwords, &rcinp->ordr[2], "OrderZ", "ordr3");
    SeekInt(lwords, &genordr, "Order", "order");
    SeekString(lwords, qformstr, "ChargeFormat", "qform");

    // Long-ranged vdW directives 
    SeekInt(lwords, &dcinp->LRvdw, "vdwmeth", "vdwMethod");

    // MLE directives 
    SeekInt(lwords, &rcinp->nlev, "EwaldLevels", "nlev");
    SeekInt(lwords, &rcinp->PadYZ[0], "Padding1", "lpad1");
    SeekInt(lwords, &rcinp->PadYZ[1], "Padding2", "lpad2");
    SeekInt(lwords, &rcinp->PadYZ[2], "Padding3", "lpad3");
    SeekReal(lwords, &rcinp->cfac[1], "Spread2", "cfac2");
    SeekReal(lwords, &rcinp->cfac[2], "Spread3", "cfac3");
    SeekReal(lwords, &rcinp->cfac[3], "Spread4", "cfac4");
    SeekInt(lwords, &rcinp->ggordr, "GridOrder", "ggordr");

    // Free allocated memory 
    DestroyCmat(&lwords);
  }

  // Translate the charge format string 
  for (i = 0; i < strlen(qformstr); i++) {
    qformstr[i] = ToUpper(qformstr[i]);
  }
  if (strcmp(qformstr, "POINT") == 0) {
    rcinp->qform = 0;
  }
  else if (strcmp(qformstr, "GAUSS") == 0) {
    rcinp->qform = 1;

    // Set the electrostatic cutoff to zero so that direct space interactions
    // are not counted.  This will generate other issues for exclusions, but
    // those will be dealt with in the routines for making lookup tables and
    // adjusting bonded interactions.
    dcinp->Ecut = 0.0;
  }

  // If a general order is specified, apply it 
  if (genordr > 0) {
    rcinp->ordr[0] = genordr;
    rcinp->ordr[1] = genordr;
    rcinp->ordr[2] = genordr;
  }

  // Compute the Ewald coefficient and then the Gaussian spread 
  if (rcinp->S <= 0.0 && dcinp->ewcoeff <= 0.0) {
    dcinp->ewcoeff = EwaldCoefficient(dcinp->Ecut, dcinp->Dtol);
    rcinp->S = 0.5/dcinp->ewcoeff;
  }
  else {
    if (rcinp->S > 0.0) {
      dcinp->ewcoeff = 0.5/rcinp->S;
    }
    else {
      rcinp->S = 0.5/dcinp->ewcoeff;
    }
    dcinp->Dtol = (1.0 - erf(dcinp->ewcoeff*dcinp->Ecut))/dcinp->Ecut;
  }

  // Check input 
  if (dcinp->Dtol < 0.0) {
    printf("GetEwaldNamelist >> Error.  Direct sum tolerance must be a "
	   "positive real value.GetEwaldNamelist >> Value of %9.6lf is "
	   "unacceptable.\n", dcinp->Dtol);
    exit(1);
  }
  for (i = 0; i < 3; i++) {
    if (rcinp->ordr[i] < 3) {
      printf("GetEwaldNamelist >> Error.  Order must be at least 3.\n"
	     "GetEwaldNamelist >> Order was specified as %d.\n",
	     rcinp->ordr[i]);
      exit(1);
    }
  }
}

//-----------------------------------------------------------------------------
// GetForceNamelist: this function reads the &force namelist, which has no
//                   analog other than perhaps &debug in SANDER.      
//
// Arguments:                                                            
//   tj:      the trajectory control information (the force report file  
//            is considered a trajectory type of output)                 
//   inp:     the input file (mdin)                                      
//-----------------------------------------------------------------------------
static void GetForceNamelist(trajcon *tj, FILE *inp)
{
  int collect;
  char line[MAXLINE];
  cmat lwords;

  // By default, all output is written to the variable "S" 
  // in the force report                                   
  sprintf(tj->DMPvar, "S");
  tj->DMPcrd = 1;
  tj->DMPbond = 1;
  tj->DMPangl = 1;
  tj->DMPdihe = 1;
  tj->DMPrelec = 1;
  tj->DMPdelec = 1;
  tj->DMPvdw = 1;
  tj->DMPall = 1;

  collect = AdvanceToSegment(inp, "force", 1);
  while (collect == 1) {
    collect = ReadNamelistLine(line, &lwords, "GetForceNamelist", inp);
    if (collect == 0) {
      continue;
    }

    // Directives for output organization 
    SeekString(lwords, tj->DMPvar, "VarName", "var");

    // Directives for output suppression 
    SeekInt(lwords, &tj->DMPcrd, "DumpCoord", "dumpcrd");
    SeekInt(lwords, &tj->DMPbond, "DumpBond", "dumpbond");
    SeekInt(lwords, &tj->DMPangl, "DumpAngl", "dumpangl");
    SeekInt(lwords, &tj->DMPdihe, "DumpDihe", "dumpdihe");
    SeekInt(lwords, &tj->DMPrelec, "DumpRElec", "dumprelec");
    SeekInt(lwords, &tj->DMPdelec, "DumpDElec", "dumpdelec");
    SeekInt(lwords, &tj->DMPvdw, "DumpVdw", "dumpvdw");
    SeekInt(lwords, &tj->DMPall, "DumpAll", "dumpall");

    // Free allocated memory 
    DestroyCmat(&lwords);
  }
}

//-----------------------------------------------------------------------------
// GetParamNamelist: this function reads a &param namelist, which has no 
//                   analog in sander.  It is distinct from the &fit     
//                   namelist reader, because although they both do      
//                   fitting the methods are distinct.                   
//                                                                       
// Arguments:                                                            
//   tp:      the topology struct (stores its own source file)           
//   myfit:   the fitting struct                                         
//   inp:     the input file (mdin)                                      
//-----------------------------------------------------------------------------
static void GetParamNamelist(prmset *myparms, trajcon *tj, FILE *inp)
{
  int collect, maxconf, nlabels;
  int maxbadj, maxaadj, maxhadj, maxmadj, maxrecast, maxcleave;
  int maxbrst, maxarst, maxhrst, maxrrst, maxgeom, maxspectral;
  char line[MAXLINE], ParmOutput[MAXNAME], rscalemeth[MAXNAME];
  cmat lwords, glossary, nmrlabels;

  // If no parameter fitting input is present, bail right out. 
  collect = AdvanceToSegment(inp, "param", 1);
  if (collect == 0) {
    return;
  }

  // Default values 
  myparms->ep[0] = '\0';
  myparms->sr[0] = '\0';
  myparms->ao[0] = '\0';
  sprintf(myparms->WaterName, "WAT");
  sprintf(myparms->NrgUnits, "kcal/mol");
  myparms->ljbuck = 0;
  myparms->lj14fac = 2.0;
  myparms->elec14fac = 1.2;
  myparms->cmapFitResln = 0;
  myparms->nbadj = 0;
  myparms->naadj = 0;
  myparms->nhadj = 0;
  myparms->nmadj = 0;
  myparms->nradj = 0;
  myparms->nuserrstB = 0;
  myparms->nuserrstA = 0;
  myparms->nuserrstH = 0;
  myparms->nuserrstR = 0;
  myparms->nanglsum = 0;
  myparms->grstB = -1.0;          // All restraints on fitted parameters are
  myparms->grstA = -1.0;          // initialized with negative weights.  This
                                  // is so that values <= 0.0 can be detected
                                  // to indicate that restraints DO NOT apply
                                  // to a particular parameter.
  myparms->grstBcpl = 5000.0;     // Penalize a 50 kcal/mol change to a bond
                                  // spring constant by the same amount as a
                                  // 0.01A change to its equilibrium value
  myparms->grstAcpl = 114.5916;   // Penalize a 2 kcal/mol change to an angle
				  // spring constant by the same amount as a
                                  // 1 degree change to its equilibrium value
  myparms->grstH = -1.0;
  myparms->grstM = -1.0;
  myparms->grst14 = -1.0;
  myparms->grstR = -1.0;
  myparms->FitAllBonds = 0;
  myparms->FitAllAngles = 0;
  myparms->FitAllTorsions = 0;
  myparms->FitAllCMAPs = 0;
  myparms->FitAllNMROps = 0;
  myparms->FitBondEq = 0;
  myparms->FitAnglEq = 0;
  myparms->spectrum = 0;
  myparms->lpost = 0.2;
  myparms->thpost = 0.2;
  myparms->mmtol = 30.0;
  myparms->esigtol = 5.0;
  myparms->fsigtol = 4.0;
  myparms->fdevfloor = 5.0;
  myparms->wtfloor = 1000.0;
  myparms->spvtol = 1.05;
  myparms->fitscnb = 0;
  myparms->fitscee = 0;
  myparms->reportall = 1;
  myparms->zeroNonfit = 0;
  myparms->verbose = 1;
  myparms->nrecast = 0;
  myparms->ncleave = 0;
  myparms->nspectralH = 0;
  myparms->nspectralx = 8;
  myparms->RemoveOutliers = 0;
  myparms->PrintFitPoints = 0;
  myparms->RestScaleMeth = 0;
  myparms->ititl = (char*)malloc(MAXLINE*sizeof(char));
  sprintf(myparms->ititl, "Generated by mdgx executing %s.", tj->inpname);
  sprintf(ParmOutput, "standard");
  myparms->icomm = (char*)calloc(MAXLINE, sizeof(char));
  myparms->NMROpsFile = (char*)calloc(MAXLINE, sizeof(char));
  myparms->NMRParmFile = (char*)calloc(MAXLINE, sizeof(char));
  sprintf(myparms->icomm, "Fit by mdgx");

  // Pre-allocate arrays for adjustable terms specifications 
  myparms->badj = (xbonddef*)malloc(32*sizeof(xbonddef));
  myparms->aadj = (xangldef*)malloc(32*sizeof(xangldef));
  myparms->hadj = (torterm*)malloc(32*sizeof(torterm));
  myparms->madj = (xcmapdef*)malloc(32*sizeof(xcmapdef));
  myparms->userrstB = (bondrst*)malloc(32*sizeof(bondrst));
  myparms->userrstA = (anglrst*)malloc(32*sizeof(anglrst));
  myparms->userrstH = (torrst*)malloc(32*sizeof(torrst));
  myparms->anglsum = (geomrst*)malloc(32*sizeof(geomrst));
  myparms->recast = (typeswitch*)malloc(32*sizeof(typeswitch));
  myparms->cleave = (typebranch*)malloc(32*sizeof(typebranch));
  myparms->spectralH = (specreq*)malloc(32*sizeof(specreq));
  nmrlabels = CreateCmat(32, MAXNAME);
  nlabels = 0;
  maxbadj = 32;
  maxaadj = 32;
  maxhadj = 32;
  maxmadj = 32;
  maxbrst = 32;
  maxarst = 32;
  maxhrst = 32;
  maxrrst = 32;
  maxgeom = 32;
  maxrecast = 32;
  maxcleave = 32;
  maxspectral = 32;

  // Create a glossary of terms to help contain open-ended input declarations
  glossary = CreateCmat(32, 32);
  AddToGlossary(&glossary, 10, "System", "sys", "vdWstyle", "ljstyle",
		"FitBonds", "bonds", "FitAngles", "angles", "RestScale",
		"rscaling");
  AddToGlossary(&glossary, 10, "FitTorsions", "torsions", "FitBondEq",
		"bondeq", "FitAnglEq", "angleq", "FitCmaps", "cmaps",
		"FitLJ14", "fitscnb");
  AddToGlossary(&glossary, 10, "FitEE14", "fitscee", "ReportAll", "repall",
		"ParmOutput", "parmout", "ShowProgress", "verbose",
		"ElimOutliers", "elimsig");
  AddToGlossary(&glossary, 8, "FittedMMOnly", "zeromm", "PrintFitPoints",
		"printpts", "WaterName", "watnam", "EnergyUnits", "eunits");
  AddToGlossary(&glossary, 8, "SeriesReport", "srep", "AccReport", "accrep",
		"ParmTitle", "title", "Vdw14Fac", "scnb");
  AddToGlossary(&glossary, 10, "Elec14Fac", "scee", "BondRest", "brst",
		"AngleRest", "arst", "TorsionRest", "hrst", "CmapRest",
		"mrst");
  AddToGlossary(&glossary, 8, "BondCoupling", "brstcpl", "AngleCoupling",
		"arstcpl", "BondBasisSep", "lpost", "AnglBasisSep", "thpost");
  AddToGlossary(&glossary, 8, "Sc14Rest", "rst14", "ConfTol", "ctol",
		"EOutlier", "esigtol", "FOutlier", "fsigtol");
  AddToGlossary(&glossary, 10, "OutlierMinDev", "fdevfloor", "MinWeight",
		"wtmin", "SpectrumTol", "spvtol", "FitB", "fitb", "FitM",
		"fitm");
  AddToGlossary(&glossary, 12, "FitA", "fita", "FitH", "fith", "RestrainB",
		"sbrst", "RestrainA", "sarst", "CMAPDensity", "cmapdens",
		"FitM", "fitm");
  AddToGlossary(&glossary, 12, "NMRFile", "nmrfile", "NMROutput", "nmrout",
                "FitR", "fitr", "RestrainR", "srrst", "NMROpsRest", "rrst",
		"FitNMROps", "nmrops");
  AddToGlossary(&glossary, 10, "RestrainH", "shrst", "Geometry", "geom",
		"ReplaceType", "recast", "BranchType", "branch", "Spectrum",
		"spectrum");

  // The presence of a &param namelist will 
  // override directives for dynamics.      
  tj->mode = 4;
  maxconf = 1024;
  myparms->nconf = 0;
  myparms->conf = (mmsys*)malloc(maxconf*sizeof(mmsys));
  while (collect == 1) {
    collect = ReadNamelistLine(line, &lwords, "GetParamNamelist", inp);
    if (collect == 0) {
      continue;
    }

    // Input directives 
    SeekInt(lwords, &myparms->ljbuck, "vdWstyle", "ljstyle");
    SeekInt(lwords, &myparms->FitAllBonds, "FitBonds", "bonds");
    SeekInt(lwords, &myparms->FitAllAngles, "FitAngles", "angles");
    SeekInt(lwords, &myparms->FitAllTorsions, "FitTorsions", "torsions");
    SeekInt(lwords, &myparms->FitAllCMAPs, "FitCmaps", "cmaps");
    SeekInt(lwords, &myparms->FitAllNMROps, "FitNMROps", "nmrops");
    SeekInt(lwords, &myparms->FitBondEq, "FitBondEq", "bondeq");
    SeekInt(lwords, &myparms->FitAnglEq, "FitAnglEq", "angleq");
    SeekInt(lwords, &myparms->cmapFitResln, "CMAPDensity", "cmapdens");
    SeekInt(lwords, &myparms->fitscnb, "FitLJ14", "fitscnb");
    SeekInt(lwords, &myparms->fitscee, "FitEE14", "fitscee");
    SeekInt(lwords, &myparms->reportall, "ReportAll", "repall");
    SeekInt(lwords, &myparms->verbose, "ShowProgress", "verbose");
    SeekInt(lwords, &myparms->RemoveOutliers, "ElimOutliers", "elimsig");
    SeekInt(lwords, &myparms->zeroNonfit, "FittedMMOnly", "zeromm");
    SeekInt(lwords, &myparms->PrintFitPoints, "PrintFitPoints", "printpts");
    SeekString(lwords, myparms->WaterName, "WaterName", "watnam");
    SeekString(lwords, myparms->NrgUnits, "EnergyUnits", "eunits");
    SeekString(lwords, myparms->sr, "SeriesReport", "srep");
    SeekString(lwords, myparms->ao, "AccReport", "accrep");
    SeekString(lwords, myparms->ititl, "ParmTitle", "title");
    SeekString(lwords, myparms->NMROpsFile, "NMRFile", "nmrfile");
    SeekString(lwords, myparms->NMRParmFile, "NMROutput", "nmrout");
    SeekString(lwords, ParmOutput, "ParmOutput", "parmout");
    SeekString(lwords, rscalemeth, "RestScale", "rscaling");
    SeekReal(lwords, &myparms->lj14fac, "Vdw14Fac", "scnb");
    SeekReal(lwords, &myparms->elec14fac, "Elec14Fac", "scee");
    SeekReal(lwords, &myparms->grstB, "BondRest", "brst");
    SeekReal(lwords, &myparms->grstA, "AngleRest", "arst");
    SeekReal(lwords, &myparms->grstH, "TorsionRest", "hrst");
    SeekReal(lwords, &myparms->grstM, "CmapRest", "mrst");
    SeekReal(lwords, &myparms->grstR, "NMROpsRest", "rrst");
    SeekReal(lwords, &myparms->grstBcpl, "BondCoupling", "brstcpl");
    SeekReal(lwords, &myparms->grstAcpl, "AngleCoupling", "arstcpl");
    SeekReal(lwords, &myparms->lpost, "BondBasisSep", "lpost");
    SeekReal(lwords, &myparms->thpost, "AnglBasisSep", "thpost");
    SeekReal(lwords, &myparms->grst14, "Sc14Rest", "rst14");
    SeekReal(lwords, &myparms->mmtol, "ConfTol", "ctol");
    SeekReal(lwords, &myparms->esigtol, "EOutlier", "esigtol");
    SeekReal(lwords, &myparms->fsigtol, "FOutlier", "fsigtol");
    SeekReal(lwords, &myparms->fdevfloor, "OutlierMinDev", "fdevfloor");
    SeekReal(lwords, &myparms->wtfloor, "MinWeight", "wtmin");
    SeekReal(lwords, &myparms->spvtol, "SpectrumTol", "spvtol");
    SeekSinglePoint(lwords, myparms, "System", "sys", &maxconf, &glossary);
    SeekBondTermID(lwords, myparms, "FitB", "fitb", &maxbadj, 2);
    SeekBondTermID(lwords, myparms, "FitA", "fita", &maxaadj, 3);
    SeekBondTermID(lwords, myparms, "FitH", "fith", &maxhadj, 4);
    SeekBondTermID(lwords, myparms, "FitM", "fitm", &maxmadj, 7);
    SeekStringInc(lwords, &nmrlabels, "FitR", "fitr", &nlabels);
    SeekSpecParmRest(lwords, myparms, "RestrainB", "sbrst", 2, &maxbrst,
		     &glossary);
    SeekSpecParmRest(lwords, myparms, "RestrainA", "sarst", 3, &maxarst,
		     &glossary);
    SeekSpecParmRest(lwords, myparms, "RestrainH", "shrst", 4, &maxhrst,
		     &glossary);
    SeekSpecParmRest(lwords, myparms, "RestrainR", "srrst", 4, &maxrrst,
		     &glossary);
    SeekGeomRest(lwords, myparms, "Geometry", "geom", &maxgeom);
    SeekRecast(lwords, myparms, "ReplaceType", "recast", &maxrecast, 0);
    SeekRecast(lwords, myparms, "BranchType", "branch", &maxcleave, 1);
    SeekSpecReq(lwords, myparms, "Spectrum", "spectrum", &maxspectral,
		&glossary);

    // Free allocated memory 
    DestroyCmat(&lwords);
  }

  // Conversions of input
  if (rscalemeth[0] != '\0') {
    if (strcmp(rscalemeth, "estimate") == 0) {
      myparms->RestScaleMeth = 0;
    }
    else if (strcmp(rscalemeth, "spy") == 0) {
      myparms->RestScaleMeth = 1;
    }
    else {
      printf("GetParamNamelist >> Warning.  Unrecognized restraint scaling "
             "method\nGetParamNamelist >> %s.  Valid methods are 'estimate' "
             "and 'spy.'\nGetParamNamelist >> See the manual for "
             "descriptions.  The default\nGetParamNamelist >> 'estimate' "
	     "method remains in effect.\n", rscalemeth);
    }
  }
  if (CaselessStrcmp(ParmOutput, "standard", -1) == 0) {
    myparms->reportall = 1;
  }
  else if (CaselessStrcmp(ParmOutput, "stdplus", -1) == 0) {
    myparms->reportall = 2;
  }
  else if (CaselessStrcmp(ParmOutput, "frcmod", -1) == 0) {
    myparms->reportall = 0;
  }
  else {
    printf("GetParamNamelist >> Error.  Unrecognized parameter file type\n"
	   "GetParamNamelist >> %s.\n", ParmOutput);
    printf("GetParamNamelist >> Please specify 'standard' for parm.dat "
	   "format,\nGetParamNamelist >> 'stdplus' for a parm.dat including "
	   "new atom types,\nGetParamNamelist >> or 'frcmod' for a force "
	   "field modification file.\n");
    exit(1);
  }
  
  // Check the existence of an NMR output file if fitting is requested
  if (myparms->NMROpsFile[0] != '\0' && myparms->NMRParmFile[0] == '\0') {
    strcpy(myparms->NMRParmFile, myparms->NMROpsFile);
  }

  // Free allocated memory
  DestroyCmat(&glossary);
}

//-----------------------------------------------------------------------------
// ReallocNails: (re)allocate an array of nails.                         
//
// Arguments:                                                            
//   L:      the array of nails                                          
//   n:      the number of elements currently in the array               
//   s:      the size of array to allocate                               
//-----------------------------------------------------------------------------
static nail* ReallocNails(nail* L, int n, int s)
{
  int i, minelem;
  nail* NL;
  
  NL = (nail*)calloc(s, sizeof(nail));
  minelem = MIN(n, s);
  for (i = 0; i < minelem; i++) {
    NL[i] = L[i];
  }
  if (n < s) {
    for (i = n; i < s; i++) {
      NL[i].maskstr = (char*)malloc(MAXNAME*sizeof(char));
    }
  }
  if (s < n) {
    for (i = s; i < n; i++) {
      free(L[i].maskstr);
    }
  }
  if (n > 0) {
    free(L);
  }

  return NL;
}

//-----------------------------------------------------------------------------
// GetFitqNamelist: this function reads a &fitq namelist, which has no analog
//                  whatsoever in sander.                         
//
// Arguments:                                                            
//   tp:      the topology struct (stores its own source file)           
//   tj:      trajectory control information                             
//   myfit:   the fitting struct (stores lots of input data files)       
//   inp:     the input file (mdin)                                      
//-----------------------------------------------------------------------------
static void GetFitqNamelist(fset *myfit, trajcon *tj, FILE *inp)
{
  int i, collect, maxqeq, maxqmin, maxqfix, maxqsum;
  int mxsys, nsys, nsysI, nsysR;
  char line[MAXLINE], maxmemstr[MAXNAME];
  cmat lwords, qeqmask, qminmask, glossary;
  FILE *outfi;
  
  // If no charge fitting input was received, bail out. 
  collect = AdvanceToSegment(inp, "fitq", 1);
  if (collect == 0) {
    return;
  }

  // A glossary of all terms in this namelist, including aliases
  glossary = CreateCmat(32, 32);
  AddToGlossary(&glossary, 10, "FitProb", "fprob", "TestProb", "tprob",
		"Proximity", "flim", "ProbeSig", "psig", "ProbeEps",
		"peps");
  AddToGlossary(&glossary, 10, "ProbeArm", "parm", "StericLimit", "pnrg",
		"MinQWeight", "minqwt", "TetherWeight", "tetherwt",
		"HistogramBin", "hbin");
  AddToGlossary(&glossary, 10, "AcceptAll", "racc", "AcceptMax", "rmax",
		"Exclusive", "excl", "FitPoints", "nfpt", "ShowDipoles",
		"dpall");
  AddToGlossary(&glossary, 10, "Verbose", "verbose", "Tether", "tether",
		"MaxSnap", "snap", "MaxMemory", "maxmem", "EPRulesFile",
		"eprules");
  AddToGlossary(&glossary, 12, "EPExtension", "epext", "CFExtension", "cfext",
		"HistFile", "hist", "PtRecFile", "ptrecord", "AmberLibrary",
		"amblib", "LibMatching", "libmatch");
  AddToGlossary(&glossary, 12, "RespPhi", "resp", "IPolQPhi", "ipolq",
		"EqualizeQ", "equalq", "MinimizeQ", "minq", "FixQ", "fixq",
		"SumQ", "sumq");
  
  // Default input values 
  nsysI = 0;
  nsysR = 0;
  myfit->Rc = 3.0;
  myfit->Rmax = 6.0;
  myfit->psig = 3.16435;
  myfit->peps = 0.16275;
  myfit->prbarm = 0.9572;
  myfit->stericlim = 3.0;
  myfit->nfitpt = 1000;
  myfit->flimit = 0.4;
  myfit->qminwt = 1.0e-3;
  myfit->qtthwt = 1.0e-3;
  myfit->fhistbin = 0.1;
  myfit->DispAllDP = 0;
  myfit->maxsnap = 20;
  myfit->epext[0] = '\0';
  myfit->confext[0] = '\0';
  myfit->histfile[0] = '\0';
  myfit->ptrecfile[0] = '\0';
  myfit->MaxMem = 1073741824;
  myfit->verbose = 1;
  myfit->tether = 0;
  myfit->namblib = 0;
  myfit->maxamblib = 0;
  sprintf(myfit->qlibmatch, "NONE");
  maxmemstr[0] = '\0';
  
  // Allocate restraint and grid name matrices 
  // and weigthing factors for all systems     
  mxsys = 32;
  myfit->gname = CreateCmat(32, MAXNAME);
  myfit->auxgname = CreateCmat(32, MAXNAME);
  myfit->tpname = CreateCmat(32, MAXNAME);
  myfit->eprule = (char*)malloc(MAXNAME*sizeof(char));
  myfit->eprule[0] = '\0';
  myfit->wt = (double*)calloc(32, sizeof(double));
  myfit->nqeq = 0;
  myfit->nqmin = 0;
  myfit->nqfix = 0;
  myfit->nqsum = 0;
  maxqeq = 32;
  maxqmin = 32;
  maxqfix = 32;
  maxqsum = 32;
  myfit->qfixsolv = ReallocNails(myfit->qfixsolv, 0, 32);
  myfit->qfixvacu = ReallocNails(myfit->qfixvacu, 0, 32);
  myfit->qsum = ReallocNails(myfit->qsum, 0, 32);
  qeqmask = CreateCmat(32, MAXNAME);
  qminmask = CreateCmat(32, MAXNAME);

  // The presence of a &fit namelist will 
  // override directives for dynamics.    
  tj->mode = 3;
  while (collect == 1) {
    collect = ReadNamelistLine(line, &lwords, "GetFitqNamelist", inp);
    if (collect == 0) {
      continue;
    }

    // General fitting parameters 
    SeekReal(lwords, &myfit->fitprob, "FitProb", "fprob");
    SeekReal(lwords, &myfit->testprob, "TestProb", "tprob");
    SeekReal(lwords, &myfit->flimit, "Proximity", "flim");
    SeekReal(lwords, &myfit->psig, "ProbeSig", "psig");
    SeekReal(lwords, &myfit->peps, "ProbeEps", "peps");
    SeekReal(lwords, &myfit->prbarm, "ProbeArm", "parm");
    SeekReal(lwords, &myfit->stericlim, "StericLimit", "pnrg");
    SeekReal(lwords, &myfit->qminwt, "MinQWeight", "minqwt");
    SeekReal(lwords, &myfit->qtthwt, "TetherWeight", "tetherwt");
    SeekReal(lwords, &myfit->fhistbin, "HistogramBin", "hbin");
    SeekReal(lwords, &myfit->Rc, "AcceptAll", "racc");
    SeekReal(lwords, &myfit->Rmax, "AcceptMax", "rmacc");
    SeekInt(lwords, &myfit->exclusive, "Exclusive", "excl");
    SeekInt(lwords, &myfit->nfitpt, "FitPoints", "nfpt");
    SeekInt(lwords, &myfit->DispAllDP, "ShowDipoles", "dpall");
    SeekInt(lwords, &myfit->verbose, "Verbose", "verbose");
    SeekInt(lwords, &myfit->tether, "Tether", "tether");
    SeekInt(lwords, &myfit->maxsnap, "MaxSnap", "snap");
    SeekString(lwords, maxmemstr, "MaxMemory", "maxmem");
    SeekString(lwords, myfit->eprule, "EPRuleFile", "eprules");
    
    // Output file parameters 
    SeekString(lwords, myfit->epext, "EPExtension", "epext");
    SeekString(lwords, myfit->confext, "CFExtension", "cfext");
    SeekString(lwords, myfit->histfile, "HistFile", "hist");
    SeekString(lwords, myfit->ptrecfile, "PtRecFile", "ptrecord");
    SeekString(lwords, myfit->qlibmatch, "LibMatching", "libmatch");
    SeekAmberLibrary(lwords, myfit, &glossary);

    // Lists of input data files and their respective weights 
    SeekSSR(lwords, myfit->gname.map[nsysR], myfit->tpname.map[nsysR],
	    &myfit->wt[nsysR], "RespPhi", "resp", &nsysR);
    SeekS3R(lwords, myfit->gname.map[nsysI], myfit->auxgname.map[nsysI],
	    myfit->tpname.map[nsysI], &myfit->wt[nsysI], "IPolQPhi", "ipolq",
	    &nsysI);
    if (nsysI == mxsys || nsysR == mxsys) {
      mxsys += 32;
      myfit->gname = ReallocCmat(&myfit->gname, mxsys, MAXNAME);
      myfit->auxgname = ReallocCmat(&myfit->auxgname, mxsys, MAXNAME);
      myfit->tpname = ReallocCmat(&myfit->tpname, mxsys, MAXNAME);
      myfit->wt = (double*)realloc(myfit->wt, mxsys*sizeof(double));
    }

    // Lists of restraints 
    SeekStringInc(lwords, &qeqmask, "EqualizeQ", "equalq", &myfit->nqeq);
    SeekStringInc(lwords, &qminmask, "MinimizeQ", "minq", &myfit->nqmin);
    SeekFixQ(lwords, myfit->qfixsolv[myfit->nqfix].maskstr,
             &myfit->qfixsolv[myfit->nqfix].target,
             myfit->qfixvacu[myfit->nqfix].maskstr,
             &myfit->qfixvacu[myfit->nqfix].target, "FixQ", "fixq",
             &myfit->nqfix);
    SeekStringPlusVal(lwords, myfit->qsum[myfit->nqsum].maskstr,
		      &myfit->qsum[myfit->nqsum].target, "SumQ", "sumq",
		      &myfit->nqsum);
    if (myfit->nqeq == maxqeq) {
      maxqeq += 32;
      qeqmask = ReallocCmat(&qeqmask, maxqeq, MAXNAME);
    }
    if (myfit->nqmin == maxqmin) {
      maxqmin += 32;
      qminmask = ReallocCmat(&qminmask, maxqmin, MAXNAME);
    }
    if (myfit->nqfix == maxqfix) {
      myfit->qfixsolv = ReallocNails(myfit->qfixsolv, maxqfix, maxqfix+32);
      myfit->qfixvacu = ReallocNails(myfit->qfixvacu, maxqfix, maxqfix+32);
      maxqfix += 32;
    }
    if (myfit->nqsum == maxqsum) {
      myfit->qsum = ReallocNails(myfit->qsum, maxqsum, maxqsum+32);
      maxqsum += 32;
    }

    // Free allocated memory 
    DestroyCmat(&lwords);
  }

  // Allocate memory for charge equalization and minimization,
  // then free the associated buffers
  myfit->qeq = ReallocNails(myfit->qeq, 0, myfit->nqeq);
  for (i = 0; i < myfit->nqeq; i++) {
    strcpy(myfit->qeq[i].maskstr, qeqmask.map[i]);
  }
  myfit->qmin = ReallocNails(myfit->qmin, 0, myfit->nqmin);
  for (i = 0; i < myfit->nqmin; i++) {
    strcpy(myfit->qmin[i].maskstr, qminmask.map[i]);
  }
  DestroyCmat(&qeqmask);
  DestroyCmat(&qminmask);

  // Contract fitting data arrays 
  nsys = MAX(nsysI, nsysR);
  myfit->model = (nsysI > 0) ? 1 : 0;
  myfit->gname = ReallocCmat(&myfit->gname, nsys, MAXNAME);
  myfit->auxgname = ReallocCmat(&myfit->auxgname, nsys, MAXNAME);
  myfit->tpname = ReallocCmat(&myfit->tpname, nsys, MAXNAME);
  myfit->wt = (double*)realloc(myfit->wt, nsys*sizeof(double));
  myfit->flimit *= myfit->flimit;

  // Parse maximum memory string 
  if (maxmemstr[0] != '\0') {
    myfit->MaxMem = ReadNumericalShorthand(maxmemstr);
  }

  // Uppercase the permittivity for unit name matching
  StringToUpper(myfit->qlibmatch);
  
  // Checks on restraints 
  myfit->qfixsolv = ReallocNails(myfit->qfixsolv, maxqfix, myfit->nqfix);
  myfit->qfixvacu = ReallocNails(myfit->qfixvacu, maxqfix, myfit->nqfix);
  myfit->qsum = ReallocNails(myfit->qsum, maxqsum, myfit->nqsum);

  // Checks on other fitting data 
  myfit->ngrd = nsys;
  if (myfit->ngrd == 0) {
    printf("GetFitqNamelist >> Error.  No fitting data was specified.\n");
    exit(1);
  }
  if (nsysI*nsysR != 0) {
    printf("GetFitqNamelist >> Error.  Both IPolQ and standard REsP data were "
	   "specified.\n");
    exit(1);
  }
  if (myfit->Rc < 0.0) {
    printf("GetFitqNamelist >> Error.  Invalid acceptance cutoff %9.6lf "
	   "specified.\n", myfit->Rc);
    exit(1);
  }
  if (myfit->Rmax < 0.0) {
    printf("GetFitqNamelist >> Error.  Invalid acceptance maximum range "
           "%9.6lf.\n", myfit->Rmax);
    exit(1);
  }
  
  // If IPolQ is in effect, check for auxiliary grids in all systems 
  if (myfit->model == 1) {
    for (i = 0; i < myfit->ngrd; i++) {
      if (myfit->auxgname.map[i][0] == '\0') {
	printf("GetFitNamelist >> Error.  Auxiliary grid not found for system"
	       "\nGetFitNamelist >> %s.\n", myfit->gname.map[i]);
	exit(1);
      }
    }
  }

  // Free allocated memory
  DestroyCmat(&glossary);
}

//-----------------------------------------------------------------------------
// AllocRestraintStrings: allocate the strings needed by a particular
//                        restraint structure.
//
// Arguments:                                                            
//   Leash:  the restraint control data struct (probably contained       
//           within a particular trajectory control data struct)         
//-----------------------------------------------------------------------------
static void SetRestraintDefaults(rstrcon *Leash)
{
  // Specifying a &restraint namelist, or using an    
  // &ipolq namelist with a solute mask, will set off  
  // a call to this routing and activate restraints.  
  Leash->active = 1;
  Leash->usegrid = 0;
  Leash->usebelly = 0;
  Leash->XpandGrid = 1;
  Leash->GridFile = (char*)malloc(MAXNAME*sizeof(char));
  Leash->GridDefsFile = (char*)malloc(MAXNAME*sizeof(char));
  Leash->BellyMask = (char*)malloc(MAXLINE*sizeof(char));
  Leash->FrozenMask = (char*)malloc(MAXLINE*sizeof(char));
  Leash->GridFile[0] = '\0';
  Leash->GridDefsFile[0] = '\0';
  Leash->BellyMask[0] = '\0';
  Leash->FrozenMask[0] = '\0';
}

//-----------------------------------------------------------------------------
// GetRestraintNamelist: this function reads restraint information for a
//                       molecular dynamics simulation.  While certain 
//                       restraints have counterparts in Amber, some of  
//                       the restraints available in mdgx do not.
//
// Arguments:                                                            
//   tj:      the trajectory control information (a restraints struct is 
//            included within tj, activated and allocated by this input  
//            routine)                                                   
//   inp:     the input file (mdin)                                      
//-----------------------------------------------------------------------------
static void GetRestraintNamelist(trajcon *tj, FILE *inp)
{
  int collect;
  char line[MAXLINE];
  cmat lwords;

  // Allocate restraint input data 
  SetRestraintDefaults(&tj->Leash);

  // Bail out if there are no restraints 
  collect = AdvanceToSegment(inp, "restraint", 1);
  if (collect == 0) {
    tj->Leash.active = 0;
    return;
  }

  // Set defaults
  while (collect == 1) {
    collect = ReadNamelistLine(line, &lwords, "GetRestraintNamelist", inp);
    if (collect == 0) {
      continue;
    }
    SeekInt(lwords, &tj->Leash.XpandGrid, "ExpandGrid", "focus");
    SeekReal(lwords, &tj->Leash.GridScale, "GridFactor", "gridx");
    SeekString(lwords, tj->Leash.GridFile, "GridFile", "grid");
    SeekString(lwords, tj->Leash.GridDefsFile, "GridDefs", "griddefs");
    SeekString(lwords, tj->Leash.BellyMask, "MobileAtoms", "bellymask");
    SeekString(lwords, tj->Leash.FrozenMask, "FrozenAtoms", "freezemask");

    // Free allocated memory 
    DestroyCmat(&lwords);
  }

  // Set flags 
  tj->Leash.usegrid = (tj->Leash.GridFile[0] != '\0') ? 1 : 0;
  tj->Leash.usebelly = (tj->Leash.BellyMask[0] != '\0' ||
			tj->Leash.FrozenMask[0] != '\0') ? 1 : 0;
}

//-----------------------------------------------------------------------------
// CheckProgramExistence: check the existence of a program, after entering
//                        preparatory commands as part of the system call.
//
// Arguments:                                                            
//   prepcalls:    preparatory commands                                  
//   prog:         the name of the program executable                    
//-----------------------------------------------------------------------------
static void CheckProgramExistence(cmat *prepcalls, char* prog)
{
  int i, slen;
  char firstword[MAXNAME];
  char* syscall;
  struct stat fstt;

  syscall = (char*)malloc(8192*sizeof(char));
  slen = 0;
  for (i = 0; i < prepcalls->row; i++) {
    sscanf(prepcalls->map[i], "%s", firstword);
    if (strcmp(firstword, "source") != 0 && strcmp(firstword, "export") != 0) {
      continue;
    }
    sprintf(&syscall[slen], "%s ; ", prepcalls->map[i]);
    slen = strlen(syscall);
  }
  sprintf(&syscall[slen], "which %s > .mdgx.filetest ; ", prog);
  system(syscall);
  free(syscall);
  if (stat(".mdgx.filetest", &fstt) == -1) {
    printf("CheckProgramExistence >> Error.  File .mdgx.filetest was not "
	   "written.\n");
    exit(1);
  }
  if (fstt.st_size == 0) {
    printf("CheckProgramExistence >> Program %s not found.\n", prog);
    printf("CheckProgramExistence >> Abbreviated preparatory sequence "
	   "follows:\n");
    for (i = 0; i < prepcalls->row; i++) {
      sscanf(prepcalls->map[i], "%s", firstword);
      if (strcmp(firstword, "source") != 0 &&
	  strcmp(firstword, "export") != 0) {
	continue;
      }
      printf("CheckProgramExistence >> %s\n", prepcalls->map[i]);
    }
    printf("CheckProgramExistence >> If you expect that the program will be "
	   "found by\nCheckProgramExistence >> executing the sequence of "
	   "preparatory calls, this\nCheckProgramExistence >> safeguard can "
	   "be disabled by setting the checkex\nCheckProgramExistence >> "
	   "variable of the &ipolq namelist to 0.\n");
    exit(1);
  }
}

//-----------------------------------------------------------------------------
// GetIPolQNamelist: get a namelist of variables for IPolQ processing of a
//                   molecular conformation.  Input in this namelist will   
//                   modify the run parameters in a &cntrl namelist to perform
//                   molecular dynamics, but also intercept the electrostatics
//                   at stated intervals to assemble input files to the ORCA
//                   quantum chemistry package.  mdgx then makes system()
//                   calls to a quantum package and reads the resulting
//                   density grids in order to compute its own electrostatic
//                   potential grids, which then serve as inputs to another
//                   mdgx run with the &fitq namelist.
//
// Arguments:                                                            
//   tj:        trajectory control information                           
//   ipqinp:    IPolQ control data                                       
//   inp:       the input file (mdin)                                    
//-----------------------------------------------------------------------------
static void GetIPolQNamelist(trajcon *tj, ipqcon *ipqinp, FILE *inp)
{
  int i, collect, nprep, npost, nqmod, maxqmod;
  char line[MAXLINE];
  char* SoluteMask;
  FILE *finfi;
  cmat lwords;
  struct stat fstt;

  // If no charge fitting input was received, bail out. 
  collect = AdvanceToSegment(inp, "ipolq", 1);
  if (collect == 0) {
    return;
  }
  tj->mode = 5;

  // Default input variables 
  ipqinp->ntqs = 1000;
  ipqinp->nQshell = 4;
  ipqinp->nVshell = 4;
  ipqinp->nQphpt = 100;
  ipqinp->nVphpt = 20;
  ipqinp->nqframe = 10;
  ipqinp->neqstep = 10000;
  ipqinp->nblock = 4;
  ipqinp->excitation = 0;
  ipqinp->verbose = 0;
  ipqinp->retqminp = 0;
  ipqinp->retqmchk = 0;
  ipqinp->retqmout = 0;
  ipqinp->retptfi = 0;
  ipqinp->checkex = 1;
  ipqinp->CenterGrid = -1;
  ipqinp->MaxCore = 512;
  ipqinp->nQMThreads = 0;
  ipqinp->DoSinglePoint = 0;
  ipqinp->Qshell[0] = 5.0;
  ipqinp->Qshell[1] = 6.0;
  ipqinp->Qshell[2] = 7.0;
  ipqinp->Qshell[3] = -1.0;
  ipqinp->Vshell[0] = 0.0;
  ipqinp->Vshell[1] = 0.3;
  ipqinp->Vshell[2] = 0.5;
  ipqinp->Vshell[3] = 0.7;
  ipqinp->minqfac = 0.01;
  ipqinp->gdim[0] = 101;
  ipqinp->gdim[1] = 101;
  ipqinp->gdim[2] = 101;
  ipqinp->gspc[0] = 0.2;
  ipqinp->gspc[1] = 0.2;
  ipqinp->gspc[2] = 0.2;
  ipqinp->prepcalls = CreateCmat(1, MAXLINE);
  nprep = 0;
  ipqinp->postcalls = CreateCmat(1, MAXLINE);
  npost = 0;
  ipqinp->qmprog = (char*)malloc(MAXNAME*sizeof(char));
  sprintf(ipqinp->qmprog, "orca");
  ipqinp->qmpath = (char*)malloc(MAXNAME*sizeof(char));
  ipqinp->qmpath[0] = '\0';
  ipqinp->uvpath = (char*)malloc(MAXNAME*sizeof(char));
  ipqinp->uvpath[0] = '\0';
  ipqinp->fmpath = (char*)malloc(MAXNAME*sizeof(char));
  ipqinp->fmpath[0] = '\0';
  ipqinp->inpfile = (char*)malloc(MAXNAME*sizeof(char));
  sprintf(ipqinp->inpfile, "IPolQinp");
  ipqinp->outfile = (char*)malloc(MAXNAME*sizeof(char));
  sprintf(ipqinp->outfile, "IPolQout");
  ipqinp->ptqfile = (char*)malloc(MAXNAME*sizeof(char));
  sprintf(ipqinp->ptqfile, "ptqarray.dat");
  ipqinp->finfile = (char*)malloc(MAXNAME*sizeof(char));
  sprintf(ipqinp->finfile, ".mdgx.finqm");
  ipqinp->grdfile = (char*)malloc(MAXNAME*sizeof(char));
  sprintf(ipqinp->grdfile, "IPolQgrd");
  ipqinp->qmmeth = (char*)malloc(MAXNAME*sizeof(char));
  sprintf(ipqinp->qmmeth, "MP2");
  ipqinp->basis = (char*)malloc(MAXNAME*sizeof(char));
  sprintf(ipqinp->basis, "cc-pvTZ");
  ipqinp->scrdir = (char*)malloc(MAXNAME*sizeof(char));
  ipqinp->scrdir[0] = '\0';
  SoluteMask = (char*)malloc(MAXLINE*sizeof(char));
  SoluteMask[0] = '\0';
  nqmod = 0;
  maxqmod = 32;
  ipqinp->QModMask = CreateCmat(32, MAXNAME);
  ipqinp->QModVal = (double*)malloc(32*sizeof(double));
  while (collect == 1) {
    collect = ReadNamelistLine(line, &lwords, "GetIPolQNamelist", inp);
    if (collect == 0) {
      continue;
    }

    // Parameters modifying the MD 
    SeekString(lwords, SoluteMask, "SoluteMol", "solute");

    // Charge density collection 
    SeekInt(lwords, &ipqinp->ntqs, "FrameRate", "ntqs");
    SeekInt(lwords, &ipqinp->nqframe, "FrameCount", "nqframe");
    SeekInt(lwords, &ipqinp->neqstep, "EqStepCount", "nsteqlim");
    SeekInt(lwords, &ipqinp->nblock, "Blocks", "nblock");
    SeekInt(lwords, &ipqinp->verbose, "Verbose", "verbose");
    SeekReal(lwords, &ipqinp->econv, "EConverge", "econv");

    // Long-ranged charge density replication 
    SeekInt(lwords, &ipqinp->nQshell, "QShellCount", "nqshell");
    SeekInt(lwords, &ipqinp->nVshell, "VShellCount", "nvshell");
    SeekInt(lwords, &ipqinp->nQphpt, "QSpherePts", "nqphpt");
    SeekInt(lwords, &ipqinp->nVphpt, "VSpherePts", "nvphpt");
    SeekInt(lwords, &ipqinp->MaxCore, "MaxMemory", "maxcore");
    SeekInt(lwords, &ipqinp->nQMThreads, "QMThreads", "qmthreads");
    SeekInt(lwords, &ipqinp->excitation, "ExcitedState", "excitation");
    SeekReal(lwords, &ipqinp->Qshell[0], "ExpQBoundary", "qshell1");
    SeekReal(lwords, &ipqinp->Qshell[1], "QShell2", "qshell2");
    SeekReal(lwords, &ipqinp->Qshell[2], "QShell3", "qshell3");
    SeekReal(lwords, &ipqinp->Qshell[3], "QShellX", "qshellx");
    SeekReal(lwords, &ipqinp->Vshell[1], "VShell1", "vshell1");
    SeekReal(lwords, &ipqinp->Vshell[2], "VShell2", "vshell2");
    SeekReal(lwords, &ipqinp->Vshell[3], "VShell3", "vshell3");
    SeekReal(lwords, &tj->dt, "TimeStep", "dt");
    SeekReal(lwords, &ipqinp->minqfac, "MinQWeight", "minqwt");
    SeekStringPlusVal(lwords, ipqinp->QModMask.map[nqmod],
		      &ipqinp->QModVal[nqmod], "ModifyQ", "modq", &nqmod);
    if (nqmod == maxqmod) {
      maxqmod += 32;
      ipqinp->QModMask = ReallocCmat(&ipqinp->QModMask, maxqmod, MAXNAME);
      ipqinp->QModVal = (double*)realloc(ipqinp->QModVal,
					 maxqmod*sizeof(double));
    }

    // Quantum calculation 
    SeekRecord(lwords, &ipqinp->prepcalls, "QuantumPrep", "prepqm", &nprep);
    SeekRecord(lwords, &ipqinp->postcalls, "QuantumClean", "cleanqm", &npost);
    SeekString(lwords, ipqinp->qmprog, "QMPackage", "qmprog");
    SeekString(lwords, ipqinp->qmpath, "QMPath", "qmpath");
    SeekString(lwords, ipqinp->inpfile, "QMInputFile", "qmcomm");
    SeekString(lwords, ipqinp->outfile, "QMOutputFile", "qmresult");
    SeekString(lwords, ipqinp->ptqfile, "PointQFile", "ptqfi");
    SeekString(lwords, ipqinp->finfile, "QMSignal", "qmflag");
    SeekString(lwords, ipqinp->qmmeth, "QMTheory", "qmlev");
    SeekString(lwords, ipqinp->basis, "QMBasis", "basis");
    SeekString(lwords, ipqinp->scrdir, "WorkDirectory", "scrdir");
    SeekInt(lwords, &ipqinp->retqminp, "KeepQMInput", "rqminp");
    SeekInt(lwords, &ipqinp->retqmchk, "KeepQMCheckPt", "rqmchk");
    SeekInt(lwords, &ipqinp->retqmout, "KeepQMOutput", "rqmout");
    SeekInt(lwords, &ipqinp->retptfi, "KeepQCloud", "rcloud");
    SeekInt(lwords, &ipqinp->checkex, "CheckExist", "checkex");
    SeekInt(lwords, &ipqinp->DoSinglePoint, "SinglePoint", "singlept");

    // Electrostatic potential evaluation 
    SeekInt(lwords, &ipqinp->gdim[0], "UElecXBin", "unx");
    SeekInt(lwords, &ipqinp->gdim[1], "UElecYBin", "uny");
    SeekInt(lwords, &ipqinp->gdim[2], "UElecZBin", "unz");
    SeekInt(lwords, &ipqinp->CenterGrid, "CenterGrid", "cengrid");
    SeekReal(lwords, &ipqinp->gspc[0], "UElecXSpc", "uhx");
    SeekReal(lwords, &ipqinp->gspc[1], "UElecYSpc", "uhy");
    SeekReal(lwords, &ipqinp->gspc[2], "UElecZSpc", "uhz");
    SeekString(lwords, ipqinp->fmpath, "FormChkPath", "fmpath");
    SeekString(lwords, ipqinp->uvpath, "UEvalPath", "uvpath");
    SeekString(lwords, ipqinp->grdfile, "GridFile", "grid");

    // Free allocated memory
    DestroyCmat(&lwords);
  }

  // Check input data 
  if (ipqinp->ntqs <= 0) {
    printf("GetIPolQNamelist >> Error.  Frame sampling rate %d is invalid.\n",
	   ipqinp->ntqs);
    exit(1);
  }
  if (ipqinp->nqframe <= 0) {
    printf("GetIPolQNamelist >> Error.  At least one frame must be "
	   "sampled.\n");
    exit(1);
  }
  if (ipqinp->nQphpt <= 0) {
    printf("GetIPolQNamelist >> Error.  Sphere surface point count %d is "
	   "invalid.\n", ipqinp->nQphpt);
    exit(1);
  }
  if (ipqinp->nVphpt <= 0) {
    printf("GetIPolQNamelist >> Error.  Sphere surface point count %d is "
	   "invalid.\n", ipqinp->nVphpt);
    exit(1);
  }
  if (ipqinp->nQshell < 1 || ipqinp->Qshell[0] <= 0.0 ||
      (ipqinp->nQshell > 1 && ipqinp->Qshell[1] <= ipqinp->Qshell[0]) ||
      (ipqinp->nQshell > 2 && ipqinp->Qshell[2] <= ipqinp->Qshell[1])) {
    printf("GetIPolQNamelist >> Explicit charge boundary distance must be "
	   "positive.\nGetIPolQNamelist >> and increase with each successive"
	   "shell.\n");
    exit(1);
  }
  if (ipqinp->nVshell < 1 ||
      (ipqinp->nVshell > 1 && ipqinp->Vshell[1] <= 0.0) ||
      (ipqinp->nVshell > 2 && ipqinp->Vshell[2] <= 0.0) ||
      (ipqinp->nVshell > 3 && ipqinp->Vshell[3] <= 0.0)) {
    printf("GetIPolQNamelist >> Electrostatic potential sampling radius must "
	   "be positive.\nGetIPolQNamelist >> Current values are ");
    for (i = 0; i < ipqinp->nVshell; i++) {
      printf("%9.4lf ", ipqinp->Vshell[i]);
    }
    printf(".\n");
    exit(1);
  }
  if (strcmp(ipqinp->qmprog, "orca") != 0 &&
      strcmp(ipqinp->qmprog, "gaussian") != 0) {
    printf("GetIPolQNamelist >> Quantum chemical packages supported are "
	   "'orca' and 'gaussian'.\n");
    exit(1);
  }
  if (ipqinp->qmpath[0] == '\0') {
    printf("GetIPolQNamelist >> Specify path to quantum chemistry "
	   "executable.\n");
    exit(1);
  }
  if (ipqinp->uvpath[0] == '\0') {
    printf("GetIPolQNamelist >> Specify path to electrostatic potential "
	   "evaluator.\n");
    exit(1);
  }
  if (ipqinp->nqframe < 4) {
    printf("GetIPolQNamelist >> Charges from at least 4 frames must be "
	   "collected\nGetIPolQNamelist >> for statistical purposes.\n");
    exit(1);
  }
  if (strcmp(ipqinp->qmprog, "gaussian") == 0 && ipqinp->fmpath[0] == '\0') {
    printf("GetIPolQNamelist >> A path to the formchk program, which writes "
	   "formatted\nGetIPolQNamelist >> checkpoint files from binaries, "
	   "must be input.\n");
    exit(1);
  }
  if (strcmp(ipqinp->qmprog, "gaussian") == 0 && ipqinp->CenterGrid == -1) {
    ipqinp->CenterGrid = 0;
  }
  else if (strcmp(ipqinp->qmprog, "orca") == 0 && ipqinp->CenterGrid == -1) {
    ipqinp->CenterGrid = 1;
  }
  if (ipqinp->CenterGrid != 0 && ipqinp->CenterGrid != 1) {
    printf("GetIPolQNamelist >> Invalid grid centering directive %d.\n",
	   ipqinp->CenterGrid);
    exit(1);
  }
  if (ipqinp->excitation > 0 && strcmp(ipqinp->qmprog, "orca") != 0) {
    printf("GetIPolQNamelist >> Excited state calculations are only valid "
	   "in the ORCA\nGetIPolQNamelist >> quantum chemistry package with "
	   "density functional theory.  Time-dependent DFT will be used.\n");
    exit(1);
  }
  if (ipqinp->DoSinglePoint == 1) {
    if (strcmp(ipqinp->qmprog, "orca") != 0) {
      printf("GetIPolQNamelist >> Single point evaluations are only supported "
             "for the ORCA\nGetIPolQNamelist >> quantum chemistry package for "
             "the time being.\n");
      exit(1);
    }
    if (strcmp(ipqinp->qmmeth, "MP2") == 0) {
      printf("GetIPolQNamelist >> Single point evaluations in IPolQ are not "
	     "supported\nGetIPolQNamelist >> for MP2.\n");
      exit(1);
    }
    if (ipqinp->excitation == 1) {
      printf("GetIPolQNamelist >> Single point evaluations with background "
	     "charge\nGetIPolQNamelist >>  distributions are not supported "
	     "with excited state\nGetIPolQNamelist >> calculations.\n");
      exit(1);
    }
  }
  ipqinp->nqmod = nqmod;

  // Sort buffered input data 
  if (SoluteMask[0] != '\0') {
    if (tj->Leash.active == 0) {
      SetRestraintDefaults(&tj->Leash);
    }
    tj->Leash.usebelly = 1;
    sprintf(tj->Leash.FrozenMask, "%s", SoluteMask);
    tj->Leash.BellyMask[0] = '\0';
  }
  ipqinp->prepcalls = ReallocCmat(&ipqinp->prepcalls, nprep, MAXLINE);
  ipqinp->postcalls = ReallocCmat(&ipqinp->postcalls, npost, MAXLINE);

  // Check that the QM finishing file does not yet exist 
  if ((finfi = fopen(ipqinp->finfile, "r")) != NULL) {
    printf("GetIPolQNamelist >> Error.  The file %s expresses completion\n"
	   "GetIPolQNamelist >> of quantum calculations launched by the "
	   "master process\nGetIPolQNamelist >> but this file already exists."
	   "\n", ipqinp->finfile);
    fclose(finfi);
    exit(1);
  }

  // Check that executables do exist 
  if (ipqinp->checkex == 1) {
    CheckProgramExistence(&ipqinp->prepcalls, ipqinp->qmpath);
    if (strcmp(ipqinp->qmprog, "gaussian") == 0) {
      CheckProgramExistence(&ipqinp->prepcalls, ipqinp->fmpath);
    }
    CheckProgramExistence(&ipqinp->prepcalls, ipqinp->uvpath);
  }

  // Free allocated memory 
  free(SoluteMask);
}

//-----------------------------------------------------------------------------
// GetConfigsNamelist: the entryway to the configuration sampling protocol.
//                     With force field development increasingly becoming the
//                     forte of mdgx, this special run mode will provide the
//                     sampling needed to understand the range of motion that
//                     a molecule can make.  The principal operation (and the
//                     reason that this is a separate module, with its own run
//                     force and energy calculators) is to clone the positions
//                     and rearrange them into many copies of atom 1x, atom 1y,
//                     atom 1z, atom2x, ... all on different rows of a matrix.
//                     Computations of forces and energies then proceed one
//                     term at a time and get done for all copies of the
//                     molecule.  This affords superior vectorization for fast
//                     execution.  As with force field development, the systems
//                     are expected to be small, and all nonbonded interactions
//                     are computed.  The results can then be processed in
//                     batches.
//
// Arguments:
//   tj:        trajectory control information
//   cfsinp:    conformational sampling input control data
//   inp:       the input file (mdin)
//-----------------------------------------------------------------------------
static void GetConfigsNamelist(trajcon *tj, configs *cfsinp, FILE *inp)
{
  int i, j, collect, maxops, maxcombos, EXigseed;
  char line[MAXLINE];
  cmat lwords, glossary;

  // If no configurations input was received, bail out.
  collect = AdvanceToSegment(inp, "configs", 1);
  if (collect == 0) {
    return;
  }
  tj->mode = 6;

  // Allocate memory in the configs data structure and set size limits
  InitConfigs(cfsinp);
  maxops = 32;
  maxcombos = 32;

  // Default input variables
  EXigseed = (tj->igseed == 72177) ? 0 : 1;

  // A glossary of all terms in this namelist, including aliases
  glossary = CreateCmat(32, 32);
  AddToGlossary(&glossary, 10, "Verbose", "verbose", "Replicas", "count",
		"MaxCycles", "maxcyc", "SDsteps", "ncyc", "MaxMemory",
		"maxcore");
  AddToGlossary(&glossary, 10, "CPUCount", "ncpu", "Multiplicity", "spin",
		"ForceConverge", "frctol", "InitialStep", "step0",
		"StepConverge", "steptol");
  AddToGlossary(&glossary, 12, "StrainLimit", "strainlim", "BondStrain",
		"bstrain", "AngleStrain", "astrain", "OutputBase", "outbase",
                "simEtol", "ESimilarity", "rmsdtol", "Similarity");
  AddToGlossary(&glossary, 12, "OutputSuffix", "outsuff", "OutputType",
		"write", "QMTheory", "qmlev", "QMBasis", "basis",
		"MovingAtoms", "belly", "DoRATTLE", "rigidbond");
  AddToGlossary(&glossary, 12, "Checkpoint", "chk", "RandomSeed", "ig",
		"RandomSample", "random", "GridSample", "uniform",
		"RattleTol", "tol", "MaxRattleIter", "vlimit");
  AddToGlossary(&glossary, 10, "RandomPerturb", "rpert", "GridPerturb",
		"gpert", "Set", "set", "Designation", "output", "FreezePolarH",
		"freezeh");
  AddToGlossary(&glossary, 10, "PreventHBonds", "nohbond", "SpinOH", "spinoh",
		"KSpinOH", "KrstSpinOH", "FBHWSpinOH", "fbhwSpinOH",
		"HBondMargin", "hblimit");

  // Loop through the input file to seek out directives
  while (collect == 1) {
    collect = ReadNamelistLine(line, &lwords, "GetIPolQNamelist", inp);
    if (collect == 0) {
      continue;
    }

    // Some of the input can be handled by generic routines
    SeekInt(lwords, &cfsinp->verbose, "Verbose", "verbose");
    SeekInt(lwords, &cfsinp->rattle, "DoRATTLE", "rigidbond");
    SeekInt(lwords, &cfsinp->MaxRattleIter, "MaxRattleIter", "vlimit");
    SeekInt(lwords, &cfsinp->count, "Replicas", "count");
    SeekInt(lwords, &cfsinp->maxcyc, "MaxCycles", "maxcyc");
    SeekInt(lwords, &cfsinp->ncyc, "SDSteps", "ncyc");
    SeekInt(lwords, &cfsinp->QMsettings.MaxCore, "MaxMemory", "maxcore");
    SeekInt(lwords, &cfsinp->QMsettings.ncpu, "CPUCount", "ncpu");
    SeekInt(lwords, &cfsinp->QMsettings.spin, "Multiplicity", "spin");
    SeekInt(lwords, &cfsinp->reshuffle, "ShuffleCount", "nshuffle");
    SeekInt(lwords, &cfsinp->showorigins, "ShowOrigins", "showorig");
    SeekInt(lwords, &cfsinp->atomlimit, "ExclTableSize", "exclmax");
    SeekInt(lwords, &cfsinp->QMsettings.opt, "OptimizeGeom", "optgeom");
    SeekInt(lwords, &cfsinp->QMsettings.xferRst, "ConstrainGeom", "cnstgeom");
    SeekInt(lwords, &cfsinp->FreezePolarH, "FreezePolarH", "freezeh");
    SeekInt(lwords, &cfsinp->PreventHBonds, "PreventHBonds", "nohbond");
    SeekInt(lwords, &cfsinp->SpinHydroxyls, "SpinOH", "spinoh");
    SeekReal(lwords, &cfsinp->fconv, "ForceConverge", "frctol");
    SeekReal(lwords, &cfsinp->stepconv, "StepConverge", "steptol");
    SeekReal(lwords, &cfsinp->step0, "InitialStep", "step0");
    SeekReal(lwords, &cfsinp->rattletol, "RattleTol", "tol");
    SeekReal(lwords, &cfsinp->HBMargin, "HBondMargin", "hblimit");
    SeekReal(lwords, &cfsinp->KSpinOH, "KrstSpinOH", "KSpinOH");
    SeekReal(lwords, &cfsinp->fbhwSpinOH, "FBHWSpinOH", "fbhwSpinOH");
    SeekReal(lwords, &cfsinp->strainlim, "StrainLimit", "strainlim");
    SeekReal(lwords, &cfsinp->maxbstrn, "BondStrain", "bstrain");
    SeekReal(lwords, &cfsinp->maxastrn, "AngleStrain", "astrain");
    SeekReal(lwords, &cfsinp->Ereplace, "ReplacementTol", "erep");
    SeekReal(lwords, &cfsinp->proximity, "ProximateNrg", "eprox");
    SeekReal(lwords, &cfsinp->simEtol, "ESimilarity", "simEtol");
    SeekReal(lwords, &cfsinp->rmsdtol, "Similarity", "rmsdtol");
    SeekMultiString(lwords, &cfsinp->outbase, "OutputBase", "outbase",
		    &glossary, 1);
    SeekMultiString(lwords, &cfsinp->outsuff, "OutputSuffix", "outsuff",
		    &glossary, 1);
    SeekMultiString(lwords, &cfsinp->ostyle, "OutputType", "write",
		    &glossary, 1);
    SeekMultiString(lwords, &cfsinp->belly, "MovingAtoms", "belly",
		    &glossary, 1);
    SeekString(lwords, cfsinp->QMsettings.theory, "QMTheory", "qmlev");
    SeekString(lwords, cfsinp->QMsettings.basis, "QMBasis", "basis");
    SeekString(lwords, cfsinp->QMsettings.checkpoint, "Checkpoint", "chk");
    SeekString(lwords, cfsinp->shuffletype, "ShuffleStyle", "shuffle");
    SeekString(lwords, cfsinp->shfdir, "Direction", "shfdir");
    if (EXigseed == 0) {
      SeekInt(lwords, &tj->igseed, "RandomSeed", "ig");
    }

    // Most of the input variables in the configs namelist will take
    // complex inputs, and so the Parse library grows.
    SeekManipulator(lwords, cfsinp, &maxops, "RandomSample", "random",
		    &glossary);
    SeekManipulator(lwords, cfsinp, &maxops, "GridSample", "uniform",
		    &glossary);
    SeekManipulator(lwords, cfsinp, &maxops, "RandomPerturb", "rpert",
		    &glossary);
    SeekManipulator(lwords, cfsinp, &maxops, "GridPerturb", "gpert",
		    &glossary);
    SeekManipulator(lwords, cfsinp, &maxops, "Set", "set", &glossary);
    SeekOpsCombo(lwords, cfsinp, &maxcombos, "LinkOperations", "combine",
		 &glossary);

    // Free allocated memory
    DestroyCmat(&lwords);
  }

  // Checks on input
  if (cfsinp->reshuffle > 0) {
    if (strcmp(cfsinp->shfdir, "up") != 0 &&
	strcmp(cfsinp->shfdir, "down") != 0) {
      printf("GetConfigsNamelist >> Error.  Invalid reoptimization direction "
	     "%s.\n", cfsinp->shfdir);
      exit(1);
    }
    if (strcmp(cfsinp->shuffletype, "jackknife") != 0 &&
	strcmp(cfsinp->shuffletype, "bootstrap") != 0 &&
	strcmp(cfsinp->shuffletype, "proximity") != 0) {
      printf("GetConfigsNamelist >> Error.  Invalid reoptimization strategy "
	     "%s.\n", cfsinp->shuffletype);
      exit(1);
    }
  }
  if (cfsinp->outbase.row > 1) {
    for (i = 0; i < cfsinp->outbase.row-1; i++) {
      strcpy(cfsinp->outbase.map[i], cfsinp->outbase.map[i+1]);
    }
    cfsinp->outbase = ReallocCmat(&cfsinp->outbase, cfsinp->outbase.row-1,
				  cfsinp->outbase.col);
  }
  if (cfsinp->outsuff.row > 1) {
    for (i = 0; i < cfsinp->outsuff.row-1; i++) {
      strcpy(cfsinp->outsuff.map[i], cfsinp->outsuff.map[i+1]);
    }
    cfsinp->outsuff = ReallocCmat(&cfsinp->outsuff, cfsinp->outsuff.row-1,
				  cfsinp->outsuff.col);
  }
  if (cfsinp->ostyle.row > 1) {
    for (i = 0; i < cfsinp->ostyle.row-1; i++) {
      strcpy(cfsinp->ostyle.map[i], cfsinp->ostyle.map[i+1]);
    }
    cfsinp->ostyle = ReallocCmat(&cfsinp->ostyle, cfsinp->ostyle.row-1,
				 cfsinp->ostyle.col);
  }
  if (cfsinp->outbase.row < cfsinp->ostyle.row) {
    i = cfsinp->outbase.row;
    cfsinp->outbase = ReallocCmat(&cfsinp->outbase, cfsinp->ostyle.row,
				  cfsinp->outbase.col);
    for (j = i; j < cfsinp->outbase.row; j++) {
      strcpy(cfsinp->outbase.map[j], cfsinp->outbase.map[i-1]);
    }
  }
  if (cfsinp->outsuff.row < cfsinp->ostyle.row) {
    i = cfsinp->outsuff.row;
    cfsinp->outsuff = ReallocCmat(&cfsinp->outsuff, cfsinp->ostyle.row,
				  cfsinp->outsuff.col);
    for (j = i; j < cfsinp->outsuff.row; j++) {
      strcpy(cfsinp->outsuff.map[j], cfsinp->outsuff.map[i-1]);
    }
  }
  cfsinp->nbelly = cfsinp->belly.row;
  if (cfsinp->nbelly == 1 && cfsinp->belly.map[0][0] == '\0') {
    cfsinp->nbelly = 0;
  }

  // Clean up the input for interpretation
  for (i = 0; i < cfsinp->ostyle.row; i++) {
    StringToUpper(cfsinp->ostyle.map[i]);

    // Final check on the things to write, to make sure there will be
    // valid output at the end of what could be a significant calculation
    if (strcmp(cfsinp->ostyle.map[i], "INPCRD") != 0 &&
	strcmp(cfsinp->ostyle.map[i], "PDB") != 0 &&
	strcmp(cfsinp->ostyle.map[i], "CRD") != 0 &&
	strcmp(cfsinp->ostyle.map[i], "CDF") != 0 &&
	strcmp(cfsinp->ostyle.map[i], "GAUSSIAN") != 0 &&
	strcmp(cfsinp->ostyle.map[i], "ORCA") != 0 &&
	strcmp(cfsinp->ostyle.map[i], "QRST") != 0 &&
	strcmp(cfsinp->ostyle.map[i], "FREEZE") != 0) {
      printf("GetConfigsNamelist >> Error.  Invalid output format %s.\n",
	     cfsinp->ostyle.map[i]);
      exit(1);
    }
  }
  StringToLower(cfsinp->QMsettings.theory);
  StringToLower(cfsinp->QMsettings.basis);

  // Free allocated memory
  DestroyCmat(&glossary);
}

//-----------------------------------------------------------------------------
// GetSPEvalNamelist: the start of the QM single point reading and catalogging
//                    protocol.  This relatively simple module takes as input
//                    a list of coordinate and energy items (within each item
//                    there can be one or more coordinate sets and energies),
//                    extracts the proper information from quantum calculation
//                    output files if necessary, and reprints the coordinates
//                    for successful calculations with their corresponding
//                    energies.
//
// Arguments:
//   tj:        trajectory control information
//   spelist:   single point evaluation input control data
//   inp:       the input file (mdin)
//-----------------------------------------------------------------------------
static void GetSPEvalNamelist(trajcon *tj, spdata *spelist, FILE *inp)
{
  int i, j, collect, maxitems, maxcombos;
  char line[MAXLINE];
  cmat lwords, glossary;

  // If no configurations input was received, bail out.
  collect = AdvanceToSegment(inp, "speval", 1);
  if (collect == 0) {
    return;
  }
  tj->mode = 7;

  // A glossary of all terms in this namelist, including aliases
  glossary = CreateCmat(32, 32);
  AddToGlossary(&glossary, 6, "Verbose", "verbose", "PrintAllItems",
		"printall", "DataPoint", "data");

  // Default settings
  spelist->PrintItems = 0;
  spelist->enfCluster = 0;
  spelist->clstWidth = 0.0;

  // Loop through the input file to seek out directives
  maxitems = 32;
  spelist->nitems = 0;
  spelist->items = (sptask*)malloc(maxitems*sizeof(sptask));
  while (collect == 1) {
    collect = ReadNamelistLine(line, &lwords, "GetSPEvalNamelist", inp);
    if (collect == 0) {
      continue;
    }

    // Some of the input can be handled by generic routines
    SeekInt(lwords, &spelist->verbose, "Verbose", "verbose");
    SeekInt(lwords, &spelist->PrintItems, "PrintAllItems", "printall");
    SeekReal(lwords, &spelist->clstWidth, "Cluster", "cluster");
    SeekSPItem(lwords, spelist, &maxitems, "DataPoint", "data", &glossary);

    // Free allocated memory
    DestroyCmat(&lwords);
  }
  spelist->items = (sptask*)realloc(spelist->items,
				    spelist->nitems * sizeof(sptask));
  
  // Turn on cluster checking if a width has been given
  if (spelist->clstWidth > 0.0) {
    spelist->enfCluster = 1;
  }

  // Free allcoated memory
  DestroyCmat(&glossary);
}

//-----------------------------------------------------------------------------
// GetLoopNamelist: function for reading the &loop namelist and controlling
//                  loop reconstruction / prediction.
//
// Arguments:
//   tj:       trajectory control information (for random seed specification,
//             among other details)
//   lpbinp:   loop building control data
//   inp:      pointer to the open input file
//-----------------------------------------------------------------------------
static void GetLoopNamelist(trajcon *tj, loopcon *lpbinp, FILE *inp)
{
  int i, collect;
  char* line;
  cmat lwords, glossary;
  void InitLoopControl(loopcon *);
  
  // If no configurations input was received, bail out.
  collect = AdvanceToSegment(inp, "loop", 1);
  if (collect == 0) {
    return;
  }
  tj->mode = 8;

  // A glossary of all terms in this namelist, including aliases
  glossary = CreateCmat(32, 32);
  AddToGlossary(&glossary, 6, "Verbose", "verbose", "MaxStretch", "longbond",
		"ReconLength", "lreconst");
  
  // Default settings
  InitLoopControl(lpbinp);
  line = (char*)malloc(MAXLINE*sizeof(char));
  while (collect == 1) {
    collect = ReadNamelistLine(line, &lwords, "GetLoopNamelist", inp);
    if (collect == 0) {
      continue;
    }

    // Some of the input can be handled by generic routines
    SeekInt(lwords, &lpbinp->verbose, "Verbose", "verbose");
    SeekInt(lwords, &lpbinp->nbundle, "BundleCount", "nbundle");
    SeekInt(lwords, &lpbinp->npath, "BundleSize", "npath");
    SeekReal(lwords, &lpbinp->maxstretch, "MaxStretch", "longbond");
    SeekReal(lwords, &lpbinp->lreconst, "ReconLength", "lreconst");
    SeekMultiString(lwords, &lpbinp->belly, "MovingAtoms", "belly",
		    &glossary, 1);

    // Free allocated memory
    DestroyCmat(&lwords);
  }

  // Free allocated memory
  free(line);
  DestroyCmat(&glossary);
}

//-----------------------------------------------------------------------------
// InitPptdControl: initialize values for a peptide control namelist.
//
// Arguments:
//   ppctrl: the peptide control data to initialize
//-----------------------------------------------------------------------------
static void InitPptdControl(pepcon *ppctrl)
{
  ppctrl->nsys = 0;
  ppctrl->remd = 0;
  ppctrl->igb  = 6;
  ppctrl->dielectric = 80.0;
  ppctrl->nBondSteps = -1;
  ppctrl->rattle = 0;
  ppctrl->gboffset = -1.0e8;
  ppctrl->tpinames = CreateCmat(32, MAXNAME);
  ppctrl->tpfnames = CreateCmat(32, MAXNAME);
  ppctrl->outbases = CreateCmat(32, MAXNAME);
  ppctrl->crdnames = CreateCmat(32, MAXNAME);
  ppctrl->trjbases = CreateCmat(32, MAXNAME);
  ppctrl->rstrtbases = CreateCmat(32, MAXNAME);
  ppctrl->blockreq[0] = '\0';
  ppctrl->Tranges = (double2*)malloc(32 * sizeof(double2));
  ppctrl->Treplicas = (int*)malloc(32 * sizeof(int));
  ppctrl->Preplicas = (int*)malloc(32 * sizeof(int));
  ppctrl->Sreplicas = (int*)malloc(32 * sizeof(int));
}

//-----------------------------------------------------------------------------
// GetPptdNamelist: function for getting peptide implicit simulations input.
//
// Arguments:
//   tj:      trajectory control data (contains &cntrl namelist information)
//   ppctrl:  peptide simulations control data (file names, directives for
//            certain enhanced sampling methods)
//   inp:     pointer to the mdgx command file
//-----------------------------------------------------------------------------
static void GetPptdNamelist(trajcon *tj, pepcon *ppctrl, FILE *inp)
{
  int collect, maxsys;
  char* line;
  cmat lwords, glossary;

  // If no configurations input was received, bail out.
  collect = AdvanceToSegment(inp, "pptd", 1);
  if (collect == 0) {
    return;
  }
  tj->mode = 9;
  
  // A glossary of all terms in this namelist, including aliases
  glossary = CreateCmat(32, 32);
  AddToGlossary(&glossary, 10, "Peptide", "oligomer", "DoREMD", "remd",
		"GBStyle", "igb", "MinorSteps", "bondstep", "Dielectric",
		"diel");
  AddToGlossary(&glossary, 4, "GBOffset", "offset", "BlockSize", "block");

  // Default settings
  maxsys = 32;
  InitPptdControl(ppctrl);
  line = (char*)malloc(MAXLINE*sizeof(char));
  while (collect == 1) {
    collect = ReadNamelistLine(line, &lwords, "GetPptdNamelist", inp);
    if (collect == 0) {
      continue;
    }

    // Some of the input can be handled by generic routines
    SeekPeptideSystem(lwords, ppctrl, &maxsys, "Peptide", "oligomer",
		      &glossary);
    SeekInt(lwords, &ppctrl->remd, "DoREMD", "remd");
    SeekInt(lwords, &ppctrl->igb, "GBStyle", "igb");
    SeekInt(lwords, &ppctrl->nBondSteps, "MinorSteps", "bondstep");
    SeekReal(lwords, &ppctrl->dielectric, "Dielectric", "diel");
    SeekReal(lwords, &ppctrl->gboffset, "GBOffset", "offset");
    SeekString(lwords, ppctrl->blockreq, "BlockSize", "block");
    
    // Free allocated memory
    DestroyCmat(&lwords);
  }

  // Clean up inputs
  if (ppctrl->blockreq[0] != '\0') {
    StringToUpper(ppctrl->blockreq);
    if (strcmp(ppctrl->blockreq, "SMALL") != 0 &&
	strcmp(ppctrl->blockreq, "MEDIUM") != 0 &&
	strcmp(ppctrl->blockreq, "LARGE") != 0) {
      printf("GetPptdNamelist >> Error.  Invalid block size input %s.\n"
	     "GetPptdNamelist >> Valid selections are SMALL, MEDIUM, and "
	     "LARGE.\n", ppctrl->blockreq);
      exit(1);
    }
  }
  
  // Free allocated memory
  free(line);
  DestroyCmat(&glossary);
}

//-----------------------------------------------------------------------------
// GetParmEditNamelist: function for getting topology file editing directives.
//
// Arguments:
//   tj:      trajectory control data (contains the run mode)
//   pmedit:  small data structure to store prmtop editing file names
//   inp:     pointer to the mdgx command file
//-----------------------------------------------------------------------------
static void GetParmEditNamelist(trajcon *tj, parmed *pmedit, FILE *inp)
{
  int collect;
  char* line;
  cmat lwords, glossary;

  // If no such namelist is present, bail out.
  collect = AdvanceToSegment(inp, "parmedit", 1);
  if (collect == 0) {
    return;
  }
  tj->mode = 10;

  // A glossary of all terms in this namelist, including aliases
  glossary = CreateCmat(32, 32);
  AddToGlossary(&glossary, 10, "NewPrmtop", "newprmtop", "NewInpcrd",
                "newinpcrd", "PrmtopTitle", "title", "TopologyCheck", "tchk",
		"ChargeSanity", "qtol");

  // Default settings
  pmedit->topchk = 1;
  pmedit->qtol = 0.00005;
  sprintf(pmedit->prmfile, "prmtop.edit");
  sprintf(pmedit->crdfile, "inpcrd.edit");
  sprintf(pmedit->prmtitle, "Topology produced by mdgx");

  // User input
  line = (char*)malloc(MAXLINE*sizeof(char));
  while (collect == 1) {
    collect = ReadNamelistLine(line, &lwords, "GetParmEditNamelist", inp);
    if (collect == 0) {
      continue;
    }
    SeekInt(lwords, &pmedit->topchk, "TopologyCheck", "tchk");
    SeekReal(lwords, &pmedit->qtol, "ChargeSanity", "qtol");
    SeekString(lwords, pmedit->prmfile, "NewPrmtop", "newprmtop");
    SeekString(lwords, pmedit->crdfile, "NewInpcrd", "newinpcrd");
    SeekString(lwords, pmedit->prmtitle, "PrmtopTitle", "title");
    
    // Free allocated memory
    DestroyCmat(&lwords);
  }

  // Check input
  if (fabs(pmedit->qtol - 0.00005) > 1.0e-8 && pmedit->topchk == 0) {
    printf("GetParmEditNamelist >> A tolerance for residue integral net "
	   "charges was\nGetParmEditNamelist >> entered, but topology "
	   "checking was disabled.\nGetParmEditNamelist >> A full topology "
	   "checkup will be performed.\n");
  }

  // Free allocated memory
  free(line);
  DestroyCmat(&glossary);
}

//-----------------------------------------------------------------------------
// ReadCommandFile: function for reading an mdgx command file.
//
// Arguments:                                                            
//   dcinp:   direct space command information                           
//   rcinp:   reciprocal space command information                       
//   tp:      the topology struct (stores its own source file)           
//   tj:      the trajectory struct (stores the restart file name, as    
//              well as many other output file names)                    
//   myfit:   charge fitting command data                                
//   myparms: parameter fitting control data                             
//   ipqinp:  IPolQ input control data                                   
//   cfsinp:  configuration sampling input control data
//   spelist: single point evaluation control data
//   lpbinp:  loop building control data
//   ppctrl:  peptide implciit solvent simulations control data
//   source:  the command input file (mdin)                              
//-----------------------------------------------------------------------------
void ReadCommFile(dircon *dcinp, reccon *rcinp, prmtop* tp, trajcon *tj,
		  fset *myfit, prmset *myparms, ipqcon *ipqinp,
		  configs *cfsinp, spdata *spelist, loopcon *lpbinp,
		  pepcon *ppctrl, parmed *pmedit, char* source)
{
  int i, maxlinelen;
  char line[MAXLINE];
  FILE *inp;

  // Open the input file 
  if ((inp = fopen(source, "r")) == NULL) {
    printf("ReadCommFile >> Error.  Command file %s not specified.\n",
	   source);
    exit(1);
  }

  // First, scan for file information 
  GetDataFileNames(tp, tj, inp);

  // Next, scan for &cntrl namelist information 
  GetCntrlNamelist(dcinp, tp, tj, inp);

  // If there is a second sytem, it may be necessary to update 
  // certain topology information in the secondary topology.   
  // Note that any such information (copied between topologies 
  // here) is not available to distinguish the initial and     
  // final states in a TI calculation.                         
  if (tj->TI == 1) {
    tp[1].lj14fac = tp[0].lj14fac;
    tp[1].elec14fac = tp[0].elec14fac;
    tp[1].rattle = tp[0].rattle;
    tp[1].settle = tp[0].settle;
    strcpy(tp[1].WaterName, tp[0].WaterName);
  }

  // Scan for any additional information about the topology. 
  // This is done here so that the topology reader can go on 
  // reading standard AMBER topology files, but can also be  
  // queued to look for additional information that is not   
  // found in a standard topology file.                      
  for (i = 0; i < tj->ntop; i++) {
    GetTopolNamelist(&tp[i], inp);
  }

  // Now, scan for other namelist information 
  GetEwaldNamelist(dcinp, rcinp, inp);
  GetRestraintNamelist(tj, inp);
  GetForceNamelist(tj, inp);
  GetParamNamelist(myparms, tj, inp);
  GetFitqNamelist(myfit, tj, inp);
  GetIPolQNamelist(tj, ipqinp, inp);
  GetConfigsNamelist(tj, cfsinp, inp);
  GetSPEvalNamelist(tj, spelist, inp);
  GetLoopNamelist(tj, lpbinp, inp);
  GetPptdNamelist(tj, ppctrl, inp);
  GetParmEditNamelist(tj, pmedit, inp);
  
  // Finally, scan the file verbatim into memory 
  rewind(inp);
  i = 0;
  maxlinelen = 0;
  while (fgets(line, MAXLINE, inp) != NULL) {
    i++;
    if (strlen(line) > maxlinelen) {
      maxlinelen = strlen(line);
    }
  }
  rewind(inp);
  tj->inptext = CreateCmat(i, maxlinelen+1);
  rewind(inp);
  i = 0;
  while (fgets(line, MAXLINE, inp) != NULL) {
    strcpy(tj->inptext.map[i], line);
    i++;
  }
  fclose(inp);
}
