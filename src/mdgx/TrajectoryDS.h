#ifndef TrajectoryStructs
#define TrajectoryStructs

#ifndef PREP_API
#ifdef MPI
#include <mpi.h>
#endif

#include "Constants.h"
#include "MatrixDS.h"
#include "TopologyDS.h"
#include "RestraintsDS.h"
#include "MPIMapDS.h"
#endif

//-----------------------------------------------------------------------------
// CoupledHoover: a structure for storing all the necessary parameters to
//                manage a coupled Hoover thermo-barostat.           
//-----------------------------------------------------------------------------
struct CoupledHoover {
  double pmass;     // The mass of the barostat
  double qmass;     // The mass of the thermostat
  double chi;       // The friction coefficient of the thermostat
  double eta;       // The friction coefficient of the barostat
  double sigma;     // Computed from the number of degrees of freedom
  double TauT;      // Time constant for temperature fluctuations
  double TauP;      // Time constant for pressure fluctuations
};
typedef struct CoupledHoover choov;

//-----------------------------------------------------------------------------
// LangevinThermostat: a structure for storing constants related to a Langevin
//                     thermostat.                             
//-----------------------------------------------------------------------------
struct LangevinThermostat {
  double gamma_ln;
  double c_implic;
  double c_explic;
  double sdfac;
};
typedef struct LangevinThermostat lnbath;
  
//-----------------------------------------------------------------------------
// TrajectoryControlData: a structure for storing all parameters for managing 
//                        a trajectory.                         
//-----------------------------------------------------------------------------
struct TrajectoryControlData {

  // Molecular dynamics / minimization run parameters
  int mode;                // The mode for this run (0 = MD, 1 = minimization,
                           //   2 = force and energy evaluation for a single
                           //   snapshot, 3 = charge fitting)
  int nsys;                // The number of systems in play (1 for standard
                           //   MD, minimization, or force evaluation, 2 for
                           //   thermodynamic integration, more for replica
                           //   exchange)
  int MySystemCount;       // The number of systems tended by this process
  int ntop;                // The number of topologies (1 for standard MD,
                           //   minimization, or force evaluation, 2 for
                           //   thermodynamic integration or replica exchange)
  long long int nstep;     // The total number of steps to perform
  long long int nfistep;   // The total number of steps to include in each
                           //   file; additional files may be written if
                           //   nstep > nfistep
  long long int currstep;  // The current step number
  int currfi;              // The current file number
  int RemoveMomentum;      // Remove net system momentum at this interval   
  int irest;               // Flag to tell whether run is a restart or begins
                           //   from initial coordiantes only (no velocities)
  int ntwr;                // The restart file writing frequency
  int ntwx;                // The coordinate trajectory file writing frequency
  int ntwv;                // The velocity trajectory file writing frequency
  int ntwf;                // The force trajectory file writing frequency
  int ntpr;                // The output diagnostics file writing frequency
  int ioutfm;              // Format for writing trajectory coordinates
  int OverwriteOutput;     // Flag to activate overwriting of trajectory and
                           //   other output information
  int Reckless;            // Flag to deactivate case inhibition; default 0, do
                           //   not deactivate inhibition of potentially buggy
                           //   or unphysical cases
  int igseed;              // The random number generator seed
  int SyncRNG;             // Flag to indicate that, in MPI implementations,
                           //   separate processes should (set to 1) or should
                           //   not (default, set to 0) be synchronized
  int rattle;              // Flag to specify whether to constrain bonds to
                           //   hydrogen via RATTLE.
  int MaxRattleIter;       // The maximum number of RATTLE iterations
  int topchk;              // Flag to activate topology checking
  int ntt;                 // Thermostat type (default 0, no thermostating; see
                           //   Manual.c for all options)
  int ntp;                 // Coordinate rescaling type (default 0, no
                           //   coordinate rescaling; set to 1 for isotropic
                           //   and 2 for anisotropic rescaling)
  int barostat;            // Choice of barostat (default 1 for Berendsen, 2
                           //   for Monte Carlo)
  int vrand;               // Interval for random velocity reset
  int MCBarostatFreq;      // The frequency at which the Monte-Carlo barostat
                           //   is applied.  By default, every 100 steps.
  int TI;                  // Flag to activate thermodynamic integration;
                           //   default 0 but set to 1 to activate T.I. along
                           //   a linear mixing path
  int mxorder;             // The order to which the mixing coefficient lambda
                           //   is raised in order to mix forces between the
                           //   initial and final states
  int nsynch;              // In TI calculations, coordinates of the the two
                           //   trajectories will be explicitly synchronized 
                           //   every nsynch steps.  Because forces between
                           //   corresponding particles are synchronized at
                           //   every step and the order of operations hitting
                           //   each particle location should be the same, the
                           //   coordinates should never differ in theory.
                           //   However, this occasional check is just to be on
                           //   the safe side.
  long int rndcon;         // The random number counter

  double mxA;              // The mixing factors for states A and B in
  double mxB;              //   thermodynamic integration calculations
  double dmxA;             // Derivatives of the mixing factors for states A
  double dmxB;             //   and B in thermodynamic integration calculations
  double starttime;        // The starting simulation time
  double currtime;         // The current simulation time
  double dt;               // The time step
  double rattletol;        // The RATTLE tolerance
  double Ttarget;          // The target temperature
  double Tinit;            // The initial temperature
  double BerendsenTCoupl;  // Time constant for Berendsen temperature coupling,
                           //   units of ps^-1, default 0.4 ps^-1
  double BerendsenPTime;   // Time constant for Berendsen pressure coupling,
                           //   units of ps^-1, default 1.0 ps^-1
  double BerendsenPCoupl;  // Pressure constant for Berendsen pressure
                           //   coupling, units of bar^-1, default 44.6 x 10^-6
  double MCBarostatFac[3]; // The rescaling constant multipliers for random
                           //   moves using a Monte-Carlo barostat.  Up to
                           //   three constants may be specified.  Initially,
                           //   only the first has a non-negative value and is
                           //   set to the default of 0.0001.  The default is
                           //   to perform isotropic rescaling, but specifying
                           //   non-negative values for other variables will
                           //   cause the system to rescale anisotropically in
                           //   the other dimensions.
  double mcdVmax;          // Monte-Carlo barostat volume rescaling margin
  double Ptarget;          // The target pressure, units of bar, default 1.0
  double lambda;           // The value of the mixing coefficient in
                           //   thermodynamic integration (TI) calculations
  double EMinStep0;        // The initial energy minimization step size
                           //   (default 0.01 Angstroms)
  double EMinStep;         // The instantaneous energy minimization step size
  double EMinTol;          // The energy minimization convergence criterion
  choov npth;              // Nose-Hoover thermobarostat parameters
  lnbath lnth;             // Langevin thermostat constants
  rstrcon Leash;           // Restraint control structure

  // Topology correlations
  prmcorr prc;             // This is the topology correlation map, created for
                           //   the purposes of running TI

  // Force report parameters
  int DMPcrd;        // Activate coordinate dumping
  int DMPbond;       // Activate bond force dumping
  int DMPangl;       // Activate angle force dumping
  int DMPdihe;       // Activate dihedral force dumping
  int DMPdelec;      // Activate direct sum electrostatic force dumping
  int DMPrelec;      // Activate reciprocal sum electrostatic force dumping
  int DMPvdw;        // Activate van-der Waals force dumping
  int DMPall;        // Activate summed (total) force dumping
  char DMPvar[32];   // Variable name for Matlab output

  // File names
  char inpname[MAXNAME];   // Input command file
  cmat ipcname;            // Input coordinates file(s)
  char dumpname[MAXNAME];  // The force dump file (for printing a comprehensive
                           //   report on all forces and energies)
  char rsrptname[MAXNAME]; // The residue report file (for printing a human-
                           //   readable description of all residues)
  char parmfile[MAXNAME];  // The force field parameter file (for fitting
                           //   torsions and other terms in a new model)
  cmat fmodfile;           // The force field parameter modification file(s)
                           //   (for fitting torsions and other terms)

  // File base names (these serve as the file names in case the
  // corresponding suffixes are NULL strings                   
  cmat rstbase;           // Restart
  cmat trjbase;           // Coordinates trajectory
  cmat velbase;           // Velocity trajectory
  cmat frcbase;           // Force trajectory
  char outbase[MAXNAME];  // Output diagnostics

  // File suffixes
  cmat rstsuff;           // Restart
  cmat trjsuff;           // Coordinates trajectory
  cmat velsuff;           // Velocity trajectory
  cmat frcsuff;           // Force trajectory
  char outsuff[32];       // Output diagnostics

  // Input file text, verbatim
  cmat inptext;
  char* inpline;

  // Parallel control data
  int tid;                   // Thread rank in MPI_COMM_WORLD, if MPI is
                             //   defined, or zero otherwise
  int nthreads;              // The total number of threads in the parallel run
#ifdef MPI
  MPI_Datatype MPI_DIRCON;   // MPI type for dircon structs (see pmeDirectDS.h)
  MPI_Datatype MPI_ATOMB;    // MPI type for atomb buffer (see CellManipDS.h)
  MPI_Datatype MPI_ATOMV;    // MPI type for atombv buffer (see CellManipDS.h)
  MPI_Datatype MPI_ATOMX;    // MPI type for atombx buffer (see CellManipDS.h)
  MPI_Datatype MPI_EPRULE;   // MPI type for eprule struct (see VirtualSites.c)
  MPI_Comm* SysComm;         // Communicators that bundle each system's CPUs
#endif
  int nCPUcluster;           // Number of tightly connected CPU clusters
  lgrp* CPUcluster;          // Clusters of tightly connected CPUs
  imat SystemCPUs;           // List of CPUs devoted to each system, as ranked
                             //   in MPI_COMM_WORLD if MPI is defined or simply
                             //   zeros if the implementation is serial
  int* MySystemDomain;       // List of systems tended by this process
};
typedef struct TrajectoryControlData trajcon;

//-----------------------------------------------------------------------------
// EnergyTracker: a structure for tracking instantaneous energy and related
//                quantities such as pressure and temperature.  
//-----------------------------------------------------------------------------
struct EnergyTracker {

  // Components of energy
  double delec;     // Direct space electrostatic energy
  double relec;     // Reciprocal space electrostatic energy
  double vdw12;     // Total repulsive vdW energy
  double vdw6;      // Total attractive (London dispersion) vdW energy
  double bond;      // Energy due to bonded terms
  double angl;      // Energy due to angle terms
  double dihe;      // Energy due to dihedral terms
  double cmap;      // Energy due to cmap terms
  double kine;      // Total kinetic energy
  double P;         // The current pressure
  double V;         // The current volume
  double T;         // The current temperature
  double dVdL;      // Derivative of the mixed potential energy function with
                    //   respect to the mixing parameter lambda

  // Sums of other energy quantities
  double elec;      // Total electrostatic energy
  double eptot;     // Total potential energy
  double etot;      // Total energy

  // Running averages
  double AVEdelec;
  double AVErelec;
  double AVEvdw;
  double AVEbond;
  double AVEangl;
  double AVEdihe;
  double AVEcmap;
  double AVEkine;
  double AVEelec;
  double AVEeptot;
  double AVEetot;
  double AVEP;
  double AVEV;
  double AVET;
  double AVEVir;
  double AVEdVdL;

  // Running standard deviations
  double RMSdelec;
  double RMSrelec;
  double RMSvdw;
  double RMSbond;
  double RMSangl;
  double RMSdihe;
  double RMScmap;
  double RMSkine;
  double RMSelec;
  double RMSeptot;
  double RMSetot;
  double RMSP;
  double RMSV;
  double RMST;
  double RMSVir;
  double RMSdVdL;

  // Virial tensor
  double Vir[9];

  // Energy decomposition
  double* BondUdc;  // Quantities for recalculating bond energy

  // Directive flags
  int updateU;      // Flag to activate system energy update
  int updateV;      // Flag to activate system virial update
  int nUdc;         // Number of types in the bond energy decomposition
  int Esummed;      // Flag to indicate whether the energy has been
                    //   summed since the last initialization
};
typedef struct EnergyTracker Energy;

#endif
