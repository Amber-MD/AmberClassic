#ifndef MultiSimDataStructures
#define MultiSimDataStructures

#include "GpuDS.h"

#define SM_ATOM_COUNT     256
#define MD_ATOM_COUNT     512
#define LG_ATOM_COUNT     928

//-----------------------------------------------------------------------------
// gpuSpecs: specifications for a GPU device critical to layout of the
//           multi-simulation problem
//-----------------------------------------------------------------------------
struct GpuDeviceConfiguration {
  int major;                 // Major architecture type (i.e. 5 = Maxwell)
  int minor;                 // Minor architecture type (i.e. 6.1 = Pascal GTX)
  int MPcount;               // Number of streaming multiprocessors
  int maxThrPerMP;           // Largest number of threads any one SM can handle
  int maxThrPerBlock;        // Largest size of any one thread block
  long long int cardMemory;  // Memory available on the card, in megabtyes
  char name[128];            // Name of the card
};
typedef struct GpuDeviceConfiguration gpuSpecs;

//-----------------------------------------------------------------------------
// msbath: a structure for storing constants related to a Langevin thermostat
//         in the context of a GPU Multi-Simulator.
//-----------------------------------------------------------------------------
struct MultiSimLangevinThermostat {
  int active;          // Flag to have Langevin dynamics turned on
  int refresh;         // Refresh rate for pseudo-random numbers
  float gamma_ln;      // Langevin collision frequency, from &cntrl namelist
  float c_implic;      // Velocity multiplier in Verlet-V
  float c_explic;      // Velocity multiplier in Verlet-C
  gpuFloat sdfac;      // Temperature scaling factor for Gaussian bumps,
                       //   one value for each system
};
typedef struct MultiSimLangevinThermostat msbath;

//-----------------------------------------------------------------------------
// gpuMultiSim: this is equivalent to base_simulationConst and gpuContext in
//              the pmemd CUDA code, but considerably simpler.  Each of the
//              arrays in this struct (data types gpu<type> are arrays) has
//              information for all systems, back to back, padded by up to 32
//              blank entries for warp coalesced reading.
//-----------------------------------------------------------------------------
struct GpuMultipleSimulator {

  // Trajectory control constants
  int        nsys;
  float      dt;              // Time step for minor (high frequency) updates
  float      dtVF;            // dt scaled by sqrt(418.4) to put it in units of
                              //   Angstroms / ps (velocities are stored in
                              //   the units of the restart files)
  double     invdtVF;         // Inverse dtVF, for constraint velocity updates
  int        igseed;          // The random number seed, passed from &cntrl
  int        ntwx;            // The Trajectory writing frequency
  int        ntpr;            // The output diagnostics writing frequency
  int        igb;             // The Generalized Born style
  int        nsgmdout;        // The number of times diagnostics (mdout) will
                              //   be printed per segment
  long int   rndcon;          // Random number control integer (synchronize
                              //   this with a trajectory control data struct)
  int        slowFrcMult;     // Multiplier for low frequency forces
                              //   (non-bonded and dihedral interactions)
                              //   during dynamics
  int        rattle;          // Flag to have bonds to hydrogen constrained to
                              //   their equilibrium lengths
  int        maxRattleIter;   // Maximum number of iterations for RATTLE
  float      rattleTol;       // Tolerance for converging RATTLE
  float      rtoldt;          // RATTLE tolerance scaled by the time step
  float      hdtm2invm;       // Conversion factor to get the plain inverse
                              //   mass out of the combined half time step and
                              //   inverse mass product
  float      hdtm2mass;       // Conversion factor to get the plain mass out of
                              //   the combined half time step and inverse mass

  // Potential function constants
  float      GBNeckScale;     // Neck GB scaling factor (igb = 7 or 8)
  float      GBNeckCut;       // Neck cutoff for igb ==7 or 8
  float      GBOffset;        // Radius offset used for various GB formulations
  float      kappa;           // Inverse of the Debye-Huckel length
  float      kscale;          // Scaling factor for GB radii in exponentials
  float      dielectric;      // The solvent dielectric constant

  // Langevin thermostat
  msbath     Tstat;           // Copy this data from a trajectory control
                              //   struct once it has been assigned

  // Convenient indexing quantities
  int        ljABoffset;      // The offset for indexing each system's Lennard-
                              //   Jones A and B force or energy tables.  Based
                              //   on the square of the largest atom type count
                              //   in any one system.

  // Atomic indexing
  gpuInt2    atomReadLimits;

  // Bond terms
  gpuInt2    bondReadLimits;
  gpuUInt    bondIDs;
  gpuFloat2  bondBasics;
  gpuFloat4  bondAugs;
  gpuInt2    cnstReadLimits;
  gpuUInt2   cnstInfo;

  // Angle terms
  gpuInt2    anglReadLimits;
  gpuUInt    anglIDs;
  gpuFloat2  anglInfo;

  // Dihedral terms
  gpuInt2    diheReadLimits;
  gpuUInt2   diheIDs;
  gpuFloat2  diheInfo;

  // Non-bonded terms and atom properties
#ifdef AMBER_PLATFORM_AMD
  gpuUInt2    pairBitMasks;   // Exclusion bit masks for GPU tiles
#else
  gpuUInt    pairBitMasks;    // Exclusion bit masks for GPU tiles
#endif
  gpuUInt    attnIDs;         // Atom IDs for all attenuations for partially
                              //   excluded pairs
  gpuFloat2  attnFactors;     // Attenuation factors for all partially excluded
                              //   pairs.  ADDING the non-bonded interactions
                              //   for each pair scaled by these factors will
                              //   produce the correct non-bonded energy.
  gpuInt4    nbReadLimits;    // Start and stop indices for reading non-bonded
                              //   exclusions.  x = start non-bonded tile bit
                              //   masks, y = stop non-bonded tile bit masks,
                              //   z = start non-bonded attenuations, w = stop
                              //   non-bonded attenuations
  gpuFloat   atomQ;           // Atomic partial charges
  gpuFloat   atomMass;        // Atomic masses
  gpuFloat   atomHDTM;        // Atomic inverse masses scaled by the time step
  gpuDouble  invNDF;          // The inverse of "the number of system degrees
                              //   of freedom times the gas constant"
  gpuInt     atomLJID;        // Atomic Lennard-Jones indices into lj[A,B]tab
  gpuInt     typeCounts;      // Numbers of atom types in each system
  gpuFloat2  ljFtab;          // Basic Lennard-Jones A and B coefficients,
  gpuFloat2  ljUtab;          //   both arrays representing very large matrices
                              //   containing tiles that can be read as if they
                              //   were symmetric matrices

  // Generalized Born calculations
  gpuInt     NeckID;          // Index into the static neck parameters array
  gpuFloat   rborn;           // Baseline Generalized Born radii
  gpuFloat   reff;            // Effective Born radii
  gpuFloat   GBalpha;         // Generalized Born alpha parameters
  gpuFloat   GBbeta;          // Generalized Born alpha parameters
  gpuFloat   GBgamma;         // Generalized Born alpha parameters
  gpuFloat   Fscreen;         // Generalized Born alpha parameters
  gpuFloat2  neckFactors;     // Neck GB max value and max position, bundled

  // Energy components
  gpuDouble  Uelec;           // Intramolecular electrostatic energy
  gpuDouble  Uvdw;            // Intramolecular electrostatic energy
  gpuDouble  Usolv;           // Solvent electrostatic contribution
  gpuDouble  Ubond;           // Bond stretching energy
  gpuDouble  Uangl;           // Angle bending energy
  gpuDouble  Udihe;           // Torsion twisting energy
  gpuDouble  Ukine;           // Kinetic energy
  gpuDouble  Temp;            // Temperature

  // Atomic coordinates
  gpuFloat   atomX;           // Atomic X coordinates
  gpuFloat   atomY;           // Atomic Y coordinates
  gpuFloat   atomZ;           // Atomic Z coordinates
  gpuFloat   crdrX;           // X coordinate snapshot for restart printing
  gpuFloat   crdrY;           // Y coordinate snapshot for restart printing
  gpuFloat   crdrZ;           // Z coordinate snapshot for restart printing

  // Atomic velocities
  gpuFloat   velcX;           // Atomic X velocities
  gpuFloat   velcY;           // Atomic Y velocities
  gpuFloat   velcZ;           // Atomic Z velocities
  gpuFloat   velrX;           // Velocity X snapshot for restart printing
  gpuFloat   velrY;           // Velocity Y snapshot for restart printing
  gpuFloat   velrZ;           // Velocity Z snapshot for restart printing

  // Atomic forces
  gpuInt     frcX;            // Atomic X forces
  gpuInt     frcY;            // Atomic Y forces
  gpuInt     frcZ;            // Atomic Z forces
  gpuInt     langX;           //
  gpuInt     langY;           // Langevin bumps in X, Y, and Z
  gpuInt     langZ;           //
#ifdef CUDA
  gpuFloat   prngHeap;        // Pseudo-random number heap
  gpuInt2    prngReadLimits;  // Portion of prngHeap that each system owns

  // Pseudo-random number states
  void       *prngStates;     // Pointer to what is really a curandState_t
                              //   struct array.  This is allocated by code in
                              //   a CUDA unit and sanitized so that no
                              //   C-compiled code has to know what such a
                              //   struct is.
#endif

  // System tracking through the block grid
  gpuInt     systemPos;       // Array of counters to track the block grid's
                              //   progress in processing all simulations.  The
                              //   first half of the array stores counters for
                              //   energy computation steps, the second for
                              //   the runs of ntpr-1 interim dynamics steps.

  // Pointers to the device memory locations of all arrays above
  int2       *DVCatomReadLimits;
  int2       *DVCbondReadLimits;
  uint       *DVCbondIDs;
  float2     *DVCbondBasics;
  float4     *DVCbondAugs;
  int2       *DVCcnstReadLimits;
  uint2      *DVCcnstInfo;
  int2       *DVCanglReadLimits;
  uint       *DVCanglIDs;
  float2     *DVCanglInfo;
  int2       *DVCdiheReadLimits;
  uint2      *DVCdiheIDs;
  float2     *DVCdiheInfo;
#ifdef AMBER_PLATFORM_AMD
  uint2      *DVCpairBitMasks;
#else
  uint       *DVCpairBitMasks;
#endif
  uint       *DVCattnIDs;
  float2     *DVCattnFactors;
  int4       *DVCnbReadLimits;
  float      *DVCatomQ;
  float      *DVCatomMass;
  float      *DVCatomHDTM;
  double     *DVCinvNDF;
  int        *DVCatomLJID;
  int        *DVCtypeCounts;
  float2     *DVCljFtab;
  float2     *DVCljUtab;
  int        *DVCNeckID;
  float      *DVCrborn;
  float      *DVCreff;
  float      *DVCGBalpha;
  float      *DVCGBbeta;
  float      *DVCGBgamma;
  float      *DVCFscreen;
  float2     *DVCneckFactors;
  float      *DVCatomX;
  float      *DVCatomY;
  float      *DVCatomZ;
  float      *DVCcrdrX;
  float      *DVCcrdrY;
  float      *DVCcrdrZ;
  float      *DVCvelcX;
  float      *DVCvelcY;
  float      *DVCvelcZ;
  float      *DVCvelrX;
  float      *DVCvelrY;
  float      *DVCvelrZ;
  int        *DVCfrcX;
  int        *DVCfrcY;
  int        *DVCfrcZ;
  int        *DVClangX;
  int        *DVClangY;
  int        *DVClangZ;
  float      *DVCsdfac;
  float      *DVCprngHeap;
  int2       *DVCprngReadLimits;
  double     *DVCuelec;
  double     *DVCuvdw;
  double     *DVCusolv;
  double     *DVCubond;
  double     *DVCuangl;
  double     *DVCudihe;
  double     *DVCukine;
  double     *DVCtemp;
  int        *DVCsystemPos;
};
typedef struct GpuMultipleSimulator gpuMultiSim;

#endif
