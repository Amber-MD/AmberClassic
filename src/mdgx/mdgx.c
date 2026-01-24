#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "CompFrc.h"
#include "Random.h"
#include "pmeRecip.h"
#include "pmeDirect.h"
#include "BSpline.h"
#include "Topology.h"
#include "ChargeMap.h"
#include "Constants.h"
#include "Macros.h"
#include "mdgxVector.h"
#include "CrdManip.h"
#include "Trajectory.h"
#include "Restraints.h"
#include "Grid.h"
#include "Matrix.h"
#include "CellManip.h"
#include "SpecialMath.h"
#include "Nonbonded.h"
#include "Bonded.h"
#include "Command.h"
#include "Constraints.h"
#include "Integrator.h"
#include "Thermostats.h"
#include "ThermoDyn.h"
#include "Barostats.h"
#include "Manual.h"
#include "mleRecip.h"
#include "Timings.h"
#include "ChargeFit.h"
#include "ParamFit.h"
#include "IPolQ.h"
#include "ConfigSamp.h"
#include "SinglePointEval.h"
#include "LoopBuilder.h"
#include "Peptide.h"
#include "ParmEdit.h"
#include "Debug.h"
#include "BroadcastCommand.h"
#include "MPITypeCast.h"
#include "MPIMap.h"

//-----------------------------------------------------------------------------
// CheckArgumentOverflow: this function is called before accessing the next 
//                        element in the argv[] list.  It checks to make sure
//                        that such an argument exists; if not, the program
//                        exits with an error message.
//                                                                     
// Arguments:                                                           
//   i:     the number of the current element in argv[] minus 1         
//   argc:  the total number of elements in the argv[] list             
//   ev:    the current element in the argv[] list                      
//-----------------------------------------------------------------------------
static void CheckArgumentOverflow(int i, int argc, char* ev)
{
  if (i < argc-2) {
    return;
  }
  else {
    printf("CheckArgumentOverflow >> Error.  No value specified for %s.\n",
	   ev);
    exit(1);
  }
}

//-----------------------------------------------------------------------------
// TakeNumberFromTag: some input stream flags may contain a traditional flag,
//                    such as -c, concatenated with a number.  In these cases,
//                    the number must be safely extracted from the input flag.
//                    A number N typically denotes that the following input
//                    argument pertains to some Nth alternate version of the
//                    system.              
//                                                                      
// Arguments:                                                           
//   tag:   the original input flag                                     
//   npos:  the first character at which the number is expected to start
//   Nmax:  the maximum allowable number for this particular flag (this is a 
//          way to impose limits on the number of systems that may be
//          requested)                                               
//-----------------------------------------------------------------------------
static int TakeNumberFromFlag(char* tag, int npos, int Nmax)
{
  int i, N;

  const int slen = strlen(tag);
  for (i = npos; i < slen; i++) {
    if (tag[i] < '0' || tag[i] > '9') {
      printf("TakeNumberFromFlag >> Error.  Flag %s does not conform to\n"
             "TakeNumberFromFlag >> flag/number format.\n", tag);
      exit(1);
    }
  }
  N = atoi(&tag[npos]);
  if (N > Nmax) {
    printf("TakeNumberFromFlag >> Error.  Flag %s exceeds the maximum "
           "allowed number\nTakeNumberFromFlag >> of alternate inputs.\n",
           tag);
    exit(1);
  }

  return N-1;
}

//-----------------------------------------------------------------------------
// RecordInputString: this function keeps a record of the input string so that
//                    any information specified there can be carried over into
//                    the output diagnostics file.    
//
// Arguments:                                                           
//   argc:     standard C input argc                                    
//   argv:     standard C input argv                                    
//   tj:       trajectory control input                                 
//-----------------------------------------------------------------------------
static void RecordInputString(int argc, char *argv[], trajcon *tj)
{
  int i, next;
  char tag[MAXNAME];

  tj->inpline = (char*)malloc(MAXLINE*sizeof(char));
  next = 0;
  for (i = 0; i < argc; i++) {
    strcpy(tag, *argv++);
    sprintf(&tj->inpline[next], "%s ", tag);
    next += strlen(tag) + 1;
  }
}

//-----------------------------------------------------------------------------
// CommandLineControl: command-line control functions for mdgx.  This function
//                     interprets data on the command line and then prints out
//                     help information or begins dynamics.
//
// Arguments:                                                           
//   argc:     standard C input argc                                    
//   argv:     standard C input argv                                    
//   tj:       trajectory control input                                 
//   tp:       topology information                                     
//-----------------------------------------------------------------------------
static void CommandLineControl(int argc, char *argv[], trajcon *tj, prmtop* tp)
{
  int i, nfmodfile, next;
  char tag[64];

  // Record the input string
  RecordInputString(argc, argv, tj);

  // Allocate memory for trajectory file names
  tj->ipcname = CreateCmat(1, MAXNAME);
  tj->rstbase = CreateCmat(1, MAXNAME);
  tj->trjbase = CreateCmat(1, MAXNAME);
  tj->velbase = CreateCmat(1, MAXNAME);
  tj->frcbase = CreateCmat(1, MAXNAME);
  tj->fmodfile = CreateCmat(1, MAXNAME);
  tj->rstsuff = CreateCmat(1, 32);
  tj->trjsuff = CreateCmat(1, 32);
  tj->velsuff = CreateCmat(1, 32);
  tj->frcsuff = CreateCmat(1, 32);

  // Defaults for critical I/O files
  sprintf(tj->ipcname.map[0], "inpcrd");
  sprintf(tj->dumpname, "forcedump.dat");
  sprintf(tj->rstbase.map[0], "restrt");
  sprintf(tj->trjbase.map[0], "mdcrd");
  sprintf(tj->velbase.map[0], "mdvel");
  sprintf(tj->frcbase.map[0], "mdfrc");
  sprintf(tj->outbase, "mdout");
  tj->outsuff[0] = '\0';
  tj->parmfile[0] = '\0';
  sprintf(tp[0].source, "prmtop");
  tp[1].source[0] = '\0';
  tp[0].eprulesource[0] = '\0';
  tp[1].eprulesource[0] = '\0';
  tj->OverwriteOutput = 0;
  tj->Reckless = 0;
  tj->SyncRNG = 0;

  // Defaults for commonly branched variables
  tj->igseed = 72177;
  tj->lambda = 0.0;
  nfmodfile = 0;

  // If no arguments, print the manual front page
  if (argc == 1 && tj->tid == 0) {
    PrintSplash(stdout);
    PrintUsage();
  }

  // Print additional documentation
  else if (argc == 2 && tj->tid == 0) {
    strcpy(tag, *++argv);
    if (strcmp(tag, "-INPUT") == 0) {
      PrintCommandLineInputOptions();
    }
    else if ((strcmp(tag, "-help") == 0 || strcmp(tag, "-HELP") == 0) ||
             strcmp(tag, "-h") == 0 || strcmp(tag, "--help") == 0)
    {
      PrintSplash(stdout);
      PrintUsage();
      exit(0);
    }
    else if (strcmp(tag, "-IFILE") == 0) {
      PrintInputFormat();
    }
    else if (strcmp(tag, "-FILES") == 0) {
      PrintFilesNamelistVariables();
    }
    else if (strcmp(tag, "-CNTRL") == 0) {
      PrintCntrlNamelistVariables();
    }
    else if (strcmp(tag, "-EWALD") == 0) {
      PrintEwaldNamelistVariables();
    }
    else if (strcmp(tag, "-FORCE") == 0) {
      PrintForceNamelistVariables();
    }
    else if (strcmp(tag, "-FITQ") == 0) {
      PrintFitqNamelistVariables();
    }
    else if (strcmp(tag, "-PARAM") == 0) {
      PrintParamNamelistVariables();
    }
    else if (strcmp(tag, "-IPOLQ") == 0) {
      PrintIPolQNamelistVariables();
    }
    else if (strcmp(tag, "-CONFIGS") == 0) {
      PrintConfigsNamelistVariables();
    }
    else if (strcmp(tag, "-PPTD") == 0) {
      PrintPptdNamelistVariables();
    }
    else if (strcmp(tag, "-RULE") == 0) {
      PrintRuleNamelistVariables();
    }
    else if (strcmp(tag, "-PARMEDIT") == 0) {
      PrintParmEditNamelistVariables();
    }
    else if (strcmp(tag, "-ATTR") == 0) {
      PrintAttributions();
    }
    else {
      printf("CommandLineControl >> Error.  Unrecognized argument %s.\n", tag);
    }
  }
#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  if (argc <= 2) {
    exit(1);
  }

  // Read command line input and execute dynamics
  for (i = 0; i < argc-1; i += 2) {
    strcpy(tag, *++argv);

    // Input command file
    if (strcmp(tag, "-i") == 0) {
      CheckArgumentOverflow(i, argc, "-i");
      strcpy(tj->inpname, *++argv);
    }

    // Variables that get changed to make branches of a run
    else if (strcmp(tag, "-igseed") == 0) {
      CheckArgumentOverflow(i, argc, "-igseed");
      tj->igseed = atoi(*++argv);
    }
    else if (strcmp(tag, "-clambda") == 0) {
      CheckArgumentOverflow(i, argc, "-clambda");
      tj->lambda = atof(*++argv);
    }

    // Topology
    else if (strcmp(tag, "-p") == 0) {
      CheckArgumentOverflow(i, argc, "-p");
      strcpy(tp[0].source, *++argv);
    }
    else if (strncmp(tag, "-p", 2) == 0) {
      CheckArgumentOverflow(i, argc, tag);
      next = TakeNumberFromFlag(tag, 2, 2);
      strcpy(tp[next].source, *++argv);
    }

    // Input coordinates
    else if (strcmp(tag, "-c") == 0) {
      CheckArgumentOverflow(i, argc, "-c");
      strcpy(tj->ipcname.map[0], *++argv);
    }
    else if (strncmp(tag, "-c", 2) == 0) {
      CheckArgumentOverflow(i, argc, tag);
      next = TakeNumberFromFlag(tag, 2, MAXSYS);
      if (next >= tj->ipcname.row) {
	tj->ipcname = ReallocCmat(&tj->ipcname, next+1, MAXNAME);
      }
      strcpy(tj->ipcname.map[next], *++argv);
    }

    // Force / energy evaluation output file
    else if (strcmp(tag, "-d") == 0) {
      CheckArgumentOverflow(i, argc, "-d");
      strcpy(tj->dumpname, *++argv);
    }

    // Output diagnostics file
    else if (strcmp(tag, "-o") == 0) {
      CheckArgumentOverflow(i, argc, "-o");
      strcpy(tj->outbase, *++argv);
    }

    // Coordinates trajectory file
    else if (strcmp(tag, "-x") == 0) {
      CheckArgumentOverflow(i, argc, "-x");
      strcpy(tj->trjbase.map[0], *++argv);
    }
    else if (strncmp(tag, "-x", 2) == 0) {
      CheckArgumentOverflow(i, argc, tag);
      next = TakeNumberFromFlag(tag, 2, MAXSYS);
      if (next >= tj->trjbase.row) {
        tj->trjbase = ReallocCmat(&tj->trjbase, next+1, MAXNAME);
      }
      strcpy(tj->trjbase.map[next], *++argv);
    }

    // Velocity trajectory file
    else if (strcmp(tag, "-v") == 0) {
      CheckArgumentOverflow(i, argc, "-v");
      strcpy(tj->velbase.map[0], *++argv);
    }
    else if (strncmp(tag, "-v", 2) == 0) {
      CheckArgumentOverflow(i, argc, tag);
      next = TakeNumberFromFlag(tag, 2, MAXSYS);
      if (next >= tj->velbase.row) {
        tj->velbase = ReallocCmat(&tj->velbase, next+1, MAXNAME);
      }
      strcpy(tj->velbase.map[next], *++argv);
    }

    // Force trajectory
    else if (strcmp(tag, "-f") == 0) {
      CheckArgumentOverflow(i, argc, "-f");
      strcpy(tj->frcbase.map[0], *++argv);
    }
    else if (strncmp(tag, "-f", 2) == 0) {
      CheckArgumentOverflow(i, argc, tag);
      next = TakeNumberFromFlag(tag, 2, MAXSYS);
      if (next >= tj->frcbase.row) {
        tj->frcbase = ReallocCmat(&tj->frcbase, next+1, MAXNAME);
      }
      strcpy(tj->frcbase.map[next], *++argv);
    }

    // Restart file
    else if (strcmp(tag, "-r") == 0) {
      CheckArgumentOverflow(i, argc, "-r");
      strcpy(tj->rstbase.map[0], *++argv);
    }
    else if (strncmp(tag, "-r", 2) == 0) {
      CheckArgumentOverflow(i, argc, tag);
      next = TakeNumberFromFlag(tag, 2, MAXSYS);
      if (next >= tj->rstbase.row) {
        tj->rstbase = ReallocCmat(&tj->rstbase, next+1, MAXNAME);
      }
      strcpy(tj->rstbase.map[next], *++argv);
    }

    // Extra point rules file
    else if (strcmp(tag, "-xpt") == 0) {
      CheckArgumentOverflow(i, argc, "-xpt");
      strcpy(tp[0].eprulesource, *++argv);
    }
    else if (strncmp(tag, "-xpt", 4) == 0) {
      CheckArgumentOverflow(i, argc, tag);
      next = TakeNumberFromFlag(tag, 3, 2);
      strcpy(tp[next].eprulesource, *++argv);
    }

    // Parameter files
    else if (strcmp(tag, "-parm") == 0) {
      CheckArgumentOverflow(i, argc, "-parm");
      strcpy(tj->parmfile, *++argv);
    }
    else if (strcmp(tag, "-fmod") == 0) {
      CheckArgumentOverflow(i, argc, "-fmod");
      strcpy(tj->fmodfile.map[nfmodfile], *++argv);
      nfmodfile++;
      if (nfmodfile == tj->fmodfile.row) {
	tj->fmodfile = ReallocCmat(&tj->fmodfile, 2*tj->fmodfile.row, MAXNAME);
      }
    }

    // Flag to overwrite existing outputs
    else if (strcmp(tag, "-O") == 0) {
      tj->OverwriteOutput = 1;
      i--;
    }

    // Flag to override case inhibitor
    else if (strcmp(tag, "-Reckless") == 0) {
      tj->Reckless = 1;
      i--;
    }

    // Flag to override desynchronization of random number generators
    else if (strcmp(tag, "-SyncRNG") == 0) {
      tj->SyncRNG = 1;
      i--;
    }

    else {
      printf("CommandLineControl >> Error.  Unrecognized tag %s.\n", tag);
      exit(1);
    }
  }

  // Clean up memory
  if (nfmodfile > 0) {
    tj->fmodfile = ReallocCmat(&tj->fmodfile, nfmodfile, MAXNAME);
  }
}

//-----------------------------------------------------------------------------
// TheGreatCaseInhibitor: this routine serves as a prophylactic against input
//                        combinations that mdgx is not yet ready to handle,
//                        even if certain input directives in isolation are
//                        accetable.  The routine simply checks for a list of
//                        possible combinations, a list which will hopefully
//                        shrink over time, but also which may not be complete.
//
// Arguments:                                                           
//   tj:    the trajectory control information                          
//   tp:    the topology                                                
//   rcinp: the reciprocal space control information                    
//   dcinp: the direct space control information                        
//-----------------------------------------------------------------------------
static void TheGreatCaseInhibitor(trajcon *tj, prmtop *tp, reccon *rcinp)
{
  int problem;

  // Assume no problems to begin
  if (tj->Reckless == 1) {
    return;
  }
  problem = 0;

  // Only the master process will perform case check
  if (tj->tid == 0) {

    // Case: Energy minimization
    if (tj->mode == 1) {
      printf("TheGreatCaseInhibitor >> Energy minimization is not yet "
	     "implemented.\n");
      problem = 1;
    }

    // Case: MLE and NPT, the MLE in mdgx is not yet ready
    // to work with anything but a Monte-Carlo barostat   
    if (rcinp->nlev > 1 && tj->ntp > 0 && tj->barostat == 1) {
      printf("TheGreatCaseInhibitor >> Unable to perform MLE in combination "
	     "with NPT.\nTheGreatCaseInhibitor >> In principle, MLE is "
	     "compatible with constant\nTheGreatCaseInhibitor >> pressure "
	     "simulations, but not as yet in mdgx.\n");
      problem = 1;
    }

    // Case: SETTLE / RATTLE and NPT
    if ((tp->settle == 1 || tp->rattle == 1) &&
	tj->ntp > 0 && tj->barostat == 1) {
      printf("TheGreatCaseInhibitor >> Unable to perform bond constraints in "
	     "combination\nTheGreatCaseInhibitor >> with a Berendsen barostat "
	     "at this time.\n");
      problem = 1;
    }

    // Case: Monte-Carlo barostat and no thermostat
    if (tj->ntp > 0 && tj->barostat == 2 && tj->ntt == 0) {
      printf("TheGreatCaseInhibitor >> Unable to perform isenthalpic ensemble "
	     "with\nTheGreatCaseInhibitor >> Monte-Carlo barostat.\n");
      problem = 1;
    }

    // Case: Thermodynamic Integration and multi-threaded run
    if (tj->TI == 1 && tj->nthreads > 1) {
      printf("TheGreatCaseInhibitor >> Parallel thermmodynamic integration "
	     "is not yet\nTheGreatCaseInhibitor >> stable.\n");
      problem = 1;
    }

    // Case: Grid-based restraints and NPT ensemble
    if (tj->Leash.active == 1 && (tj->Leash.usegrid == 1 && tj->ntp > 0)) {
      printf("TheGreatCaseInhibitor >> Grid-based restraints cannot yet be "
	     "used in conjunction\nTheGreatCaseInhibitor >> with volume "
	     "rescaling.\n");
      problem = 1;
    }
  }

  // Broadcast the problem state and exit if needed
#ifdef MPI
  MPI_Bcast(&problem, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  if (problem == 1) {
    exit(1);
  }
}

//-----------------------------------------------------------------------------
// main
//-----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  int i, isys;
  reccon rcinp;
  dircon dcinp;
  coord* crd;
  cellgrid* CG;
  bckit* PPk;
  prmtop* tp;
  trajcon tj;
  ipqcon ipqinp;
  configs cfsinp;
  spdata spelist;
  loopcon lpbinp;
  pepcon ppctrl;
  FrcTab Etab, EHtab;
  Energy* sysUV;
  execon etimers;
  cdftrj* Acdf;
  fset myfit;
  prmset myparms;
  parmed pmedit;

#ifdef MPI
  // MPI startup
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &tj.tid);
  MPI_Comm_size(MPI_COMM_WORLD, &tj.nthreads);
  DefineMPITypes(&tj);
#else
  tj.tid = 0;
  tj.nthreads = 1;
#endif
  MapProcessors(&tj);

  // Initialize timers
  InitExecon(&etimers);
  mdgxStartTimer(&etimers);

  // Command line input
  tp = (prmtop*)malloc(2*sizeof(prmtop));
  CommandLineControl(argc, argv, &tj, tp);
  if (tj.tid == 0) {

    // Parse input command file
    ReadCommFile(&dcinp, &rcinp, tp, &tj, &myfit, &myparms, &ipqinp, &cfsinp,
		 &spelist, &lpbinp, &ppctrl, &pmedit, tj.inpname);
  }
#ifdef MPI
  if (tj.nthreads > 1) {
    BroadcastInputData(&dcinp, &rcinp, &tj, tp, &ipqinp);
  }
#endif

  // Initialize random number generator
  InitPRNG(&tj);

  // Initialize Langevin thermostat, if requested
  if ((tj.mode == 9 || tj.mode < 3) && tj.ntt == 3) {
    PrepLangevinThermostat(&tj);
  }
  
  // If charge or some other type of fitting is wanted, we  
  // dive right into these routines and never come back out.
  if (tj.mode == 3 && tj.tid == 0) {
    FitCharges(&myfit, &tj, &etimers);
  }
  if (tj.mode == 4 && tj.tid == 0) {
    FitParams(&myparms, &tj);
  }

  // If configuration generation is wanted, we dive right into
  // this routine and never come back out.
  if (tj.mode == 6) {
    SampleConfigurations(&cfsinp, &tj, &tp[0]);
  }
  if (tj.mode == 8) {
    LoopBuilder(&lpbinp, &tj, &tp[0]);
  }
  if (tj.mode == 9) {

    // Although tp[0] is not used at all in these peptide simulations
    // (the topologies for each system are read in apart from the
    // standard topology structures used for PME molecular dynamics),
    // an extra points file could, in theory, be specified in the
    // &files namelist.  It will be some time before this functionality
    // can propagate to real peptide simulations.
    PeptideSimulator(&ppctrl, &tj, tp[0].eprulesource, &etimers);
  }
  
  // Read topology; the master process of MPI_COMM_WORLD will read  
  // one or two topology files, depending on how many are specified,
  // and all processes will then parse these topologies.            
  for (i = 0; i < tj.ntop; i++) {
    GetPrmTop(&tp[i], &tj, 1);
  }

  // If single point data evaluations are wanted, we dive right
  // into this routine and never come back out.
  if (tj.mode == 7) {
    ReadSinglePointEnergies(&spelist, &tj, &tp[0]);
  }

  // If topology and system reprinting is desired, divert into this
  // routine and never come back out.
  if (tj.mode == 10) {
    EditPrmtop(&tp[0], &tj, &pmedit);
  }

  // Compute force / energy lookup tables
  for (i = 0; i < tj.ntop; i++) {  
    LongRangeVDW(&tp[i], &dcinp);
  }
  Etab = DirectSpaceR2(((rcinp.qform == 0) ? dcinp.Ecut : dcinp.Vcut),
		       dcinp.lkpspc, dcinp.ewcoeff, 0, rcinp.qform);
  EHtab = DirectSpaceR2(MINNB, 0.0625*dcinp.lkpspc, dcinp.ewcoeff, 0, 
			rcinp.qform);

  // This is where we will decide which systems of perhaps
  // many this process will actually tend.  This process  
  // will then allocate memory for only those systems, and
  // deallocate only that memory at the end of the run.   
  SelectSystemsToTend(&tj);

  // This is where the an array of topologies might get      
  // created for the purpose of Hamiltonian replica exchange.
  // In that event, two topologies must be specified and the 
  // original two topologies read in will become endpoints   
  // of a much larger array of intermediate topologies.      
  // If (Hamiltonian) replica exchange is the name of the    
  // game, then only the topologies pertaining to systems    
  // tended by this process will be created, and the endpoint
  // topologies will be freed immediately if the endpoint    
  // systems are not tended by this process.  Memory must be 
  // conserved when you have hundreds of systems at once!    
  if (tj.nsys > 2) {
    tp = (prmtop*)realloc(tp, tj.nsys*sizeof(prmtop));
    tp[tj.nsys-1] = tp[1];
    for (i = 0; i < tj.MySystemCount; i++) {
      isys = tj.MySystemDomain[i];
      tp[isys] = InterpolateTopology(&tp[0], &tp[1], (double)isys/tj.nsys);
    }
  }

  // Check for cases which mdgx is not yet prepared to handle
  TheGreatCaseInhibitor(&tj, &tp[0], &rcinp);

  // Allocate system data structures
  crd = (coord*)malloc(tj.nsys*sizeof(coord));
  CG = (cellgrid*)malloc(tj.nsys*sizeof(cellgrid));
  PPk = (bckit*)malloc(tj.nsys*sizeof(bckit));
  sysUV = (Energy*)malloc(tj.nsys*sizeof(Energy));
  Acdf = (cdftrj*)malloc(3*tj.nsys*sizeof(cdftrj));

  // Read starting coordinates
  for (i = 0; i < tj.MySystemCount; i++) {
    isys = tj.MySystemDomain[i];
    crd[isys] = InitCoords(&tp[isys], &tj, isys);
    InitializeEnergy(&sysUV[isys], &tj, &tp[isys], 1);
  }
  
  // Prepare for TI
  if (tj.TI == 1) {
    tj.prc = CompTpCorr(&tp[0], &tp[1], &crd[0], &crd[1]);
  }

  // Preparation for direct space calculation
  for (i = 0; i < tj.MySystemCount; i++) {
    isys = tj.MySystemDomain[i];
    CG[isys] = CreateCellGrid(&crd[isys], &dcinp, &rcinp, &tp[isys], &tj,
			      isys);
  }

  // Preparation for reciprocal space calculation
  if (rcinp.nlev == 1) {
    PrepPME(&CG[0], &rcinp, &crd[0]);
  }
  else {
    PrepMLE(&rcinp, &crd[0], &tp[0]);
  }
  for (i = 0; i < tj.MySystemCount; i++) {
    isys = tj.MySystemDomain[i];
    PPk[isys] = CreateBCKit(&rcinp, &rcinp.QL[0], &crd[isys], &tp[isys],
			    (rcinp.nlev == 1) ? FFTW_MEASURE : FFTW_UNALIGNED);
    if (rcinp.nlev > 1) {
      PPk[isys].forwplan = rcinp.forwplan[0];
      PPk[isys].backplan = rcinp.forwplan[0];
    }
  }

  // Link cell grid, direct space, and reciprocal space calculations
  for (i = 0; i < tj.MySystemCount; i++) {
    isys = tj.MySystemDomain[i];
#ifdef MPI
    LinkCellGrid(&CG[isys], &crd[isys], &rcinp);
#else
    LinkCellGrid(&CG[isys], &rcinp);
#endif
  }

  // Prepare thermo-barostat, if requested
  if (tj.ntt == 4) {
    PrepThermoBarostat(&tp[0], &tj);
  }

  // Prepare Monte-Carlo barostat
  if (tj.ntp > 0 && tj.barostat == 2) {
    PrepMCBarostat(&tj, &crd[0]);
  }

  // Read grid-based restraints if specified
  if (tj.Leash.active == 1 && tj.Leash.usegrid == 1) {
    tj.Leash.Rgrd = ReadRestraintGrid(tj.Leash.GridFile, &crd[0]);

    char fname[64];
    cell *C;
    fbook Tgrd;
    cling Tcr;
    Tcr.ngss = 1;
    Tcr.gss = (edgauss*)malloc(sizeof(edgauss));
    Tcr.gss[0].amp = 1.0;
    Tcr.gss[0].sig = 2.8;
    Tcr.gss[0].xtns = 16.0;
    Tgrd = ScoreOnGrid(&tj.Leash.Rgrd, &Tcr, &tj, &crd[0]);

    CellPrivateGrids(&tj.Leash.Rgrd, 1, &CG[0]);

    for (i = 0; i < CG[0].ncell; i++) {
      C = &CG[0].data[i];
      sprintf(fname, "Cell%d%d%d.xplor", C->gbin[0], C->gbin[1], C->gbin[2]);
      WriteRestraintGrid(fname, &CG[0].data[i].Fscr[0], &tj);
    }

    exit(1);
  }

  // Initialize velocities, or other features of the
  // system in the case of energy minimization      
  if (tj.TI == 1) {
    InitVelocities(crd, CG, tp, &dcinp, &Etab, &EHtab, &rcinp,
		   PPk, &tj, sysUV, &etimers, Acdf, 0);
  }
  else {
    for (i = 0; i < tj.MySystemCount; i++) {
      isys = tj.MySystemDomain[i];
      InitVelocities(&crd[isys], &CG[isys], &tp[isys], &dcinp, &Etab, &EHtab,
		     &rcinp, &PPk[isys], &tj, &sysUV[isys], &etimers,
		     &Acdf[3*isys], isys);
    }
  }

  // Record initialization / startup time
  etimers.Setup = mdgxStopTimer(&etimers);
  
  // Run dynamics or energy minimization
  if (tj.mode == 0 || tj.mode == 1) {

    // Loop for dynamics
    while (tj.currstep < tj.nstep) {
      UpdateStepNumber(&tj, tj.currstep + 1);
      if (tj.mode == 0) {
	if (tj.TI == 1) {
	  Dynamics(crd, CG, tp, &dcinp, &Etab, &EHtab, &rcinp, PPk, &tj,
		   sysUV, &etimers, Acdf, 0);
	}
	else {
	  for (i = 0; i < tj.MySystemCount; i++) {
	    isys = tj.MySystemDomain[i];
	    Dynamics(&crd[isys], &CG[isys], &tp[isys], &dcinp, &Etab, &EHtab,
		     &rcinp, &PPk[isys], &tj, &sysUV[isys], &etimers,
		     &Acdf[3*isys], isys);
	  }
	}
      }
      else if (tj.mode == 1) {
	printf("mdgx >> Error, minimization not yet implemented.\n");
	exit(1);
 
#if 0
	// FIX ME!!!  Minimization does not yet work.
	for (i = 0; i < tj.MySystemCount; i++) {
	  isys = tj.MySystemDomain[i];
	  Minimization(&crd[isys], &CG[isys], &tp[isys], &dcinp, &Etab, &EHtab,
		       &rcinp, &PPk[isys], &tj, &sysUV[isys], &etimers,
		       &Acdf[3*isys]);
	}
#endif
      }
    }

    // Print the final restart file
    for (i = 0; i < tj.MySystemCount; i++) {
      isys = tj.MySystemDomain[i];
      WriteRst(&CG[isys], &crd[isys], &tp[isys], &tj, isys, 1);
    }
  }

  // Run a force calculation
  else if (tj.mode == 2) {
    PrepForceReport(&crd[0], &CG[0], &tp[0], &dcinp, &Etab, &EHtab, &rcinp,
		    &PPk[0], &tj, &etimers);
  }

  // Run an IPolQ data point
  else if (tj.mode == 5) {
    PerformIPolQ(crd, CG, tp, &dcinp, &Etab, &EHtab, &rcinp, PPk, &tj, sysUV,
		 &etimers, Acdf, &ipqinp);
  }

  // Free allocated memory
  for (i = 0; i < tj.ntop; i++) {
    FreeTopology(&tp[i]);
  }
  free(tp);
  free(Acdf);
  FreeFrcTab(&Etab);
  FreeFrcTab(&EHtab);
  DestroyRecCon(&rcinp, &CG[0]);
  for (i = 0; i < tj.MySystemCount; i++) {
    isys = tj.MySystemDomain[i];
    DestroyCoord(&crd[isys]);
    DestroyCellGrid(&CG[isys]);
    DestroyEnergyTracker(&sysUV[isys]);
    DestroyBCKit(&PPk[isys]);
  }
  DestroyTrajCon(&tj);
  free(crd);
  free(CG);
  free(sysUV);
  free(PPk);

#ifdef MPI
  // Cleanup MPI
  FreeMPITypes(&tj);
  MPI_Finalize();
#else
  // Cleanup serial FFTW
  fftw_cleanup();
#endif

  return 0;
}
