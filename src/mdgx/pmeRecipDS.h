#ifndef pmeRecipStructs
#define pmeRecipStructs

#ifndef PREP_API
#include "fftw3.h"

#include "GridDS.h"
#include "mleRecipDS.h"
#endif

struct pmeRecipControlData {

  // Smooth particle Mesh Ewald
  int qform;           // Charge format: possible values include 0 (standard
                       //   point charges) and 1 (spherical Gaussian charges)
  int ordr[3];         // Particle-to-mesh interpolation order
  int* ng;             // The mesh dimensions
  double S;            // Sigma for Gaussian charges on the mesh (or sigma for
                       //   general charges if they're all spherical Gaussians)

  // Multi-Level (Smooth Particle Mesh) Ewald
  int nlev;            // Number of levels in MLE (maximum 4)
  int nslab;           // The number of slabs to break it down into
  int nstrip;          // The number of strips to break each slab into
  int ggordr;          // Grid-to-grid interpolation order
  int PadYZ[4];        // Padding of grids at each level
  double cfac[4];      // Coarsening factors for each mesh
  g2gmap* SPrv;        
  g2gmap* SPcv;
  dbook* Urec;
  dbook* QL;           
  fftw_plan* forwplan; // Forward FFT plans for each level
  fftw_plan* backplan; // Backward FFT plans for each level
};
typedef struct pmeRecipControlData reccon;

struct BCMeshKit {
  int plans;             // Flag to indicate that FFT plans have been made
                         //   especially for this struct (the plans may be
                         //   borrowed from elsewhere)
  double SelfEcorr;      // Self energy correction for this mesh setup
  double* Bx;            // B mesh prefactors in X, Y, and Z
  double* By;            //
  double* Bz;            //
  double* mx;            // M value (reciprocal space vector index) in X, Y,
  double* my;            //   and Z
  double* mz;            //
  double* mxs;           // Shifted M values in X, Y, or Z
  double* mys;           //
  double* mzs;           //
  double* Ex;            // Exponential tables for X, Y, and Z (used for both
  double* Ey;            //   orthorhombic and non-orthorhombic unit cells)
  double* Ez;            //
  fftw_plan forwplan;    // Forward FFT plan
  fftw_plan backplan;    // Backward FFT plan
};
typedef struct BCMeshKit bckit;

#endif
