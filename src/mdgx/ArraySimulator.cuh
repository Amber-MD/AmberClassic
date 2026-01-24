#ifndef ARRAYSIM_FUNCS
#define ARRAYSIM_FUNCS

#include "MultiSimDS.h"

gpuSpecs GetGpuSpecs(int reckless);

void InitGpuPRNG(gpuMultiSim *gms, int igseed, int nblocks, int blockdim);

void SetGmsImage(gpuMultiSim *gms);

void LaunchDynamics(gpuMultiSim *gms, int blockDim, int nblocks,
		    gpuSpecs *devspc);

#endif
