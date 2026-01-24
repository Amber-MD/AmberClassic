#include <stdio.h>
#include <stdlib.h>
#if defined(AMBER_PLATFORM_AMD)
#  include <hip/hip_runtime.h>
#  include <hip/hip_runtime_api.h>
#  include <hiprand_kernel.h>
#  include "HipDefinitions.h"
#else
#  include <cuda.h>
#  include <nvml.h>
#  include <curand_kernel.h>
#endif
#include "MultiSimDS.h"
#include "KernelMacros.h"

__device__ __constant__ gpuMultiSim cGms;

#if defined(AMBER_PLATFORM_AMD)
#define INC_HEADER "kDynamicsAmd.h"
#else
#define INC_HEADER "kDynamicsNv.h"
#endif
//-----------------------------------------------------------------------------
// GetGpuSpecs: obtain specs for the GPU in the system.  This is here to have
//              all of the CUDA-specific data structures in this one unit,
//              while returning the information that's really needed in a
//              custom struct full of POD types to be incorporated into the
//              main peptide control data structure in Peptide.c.
//
// Arguments:
//   reckless:   taken from the trajectory data structure and fed in from the
//               command line, lets mdgx go ahead and run on GPUs that are
//               already in use
//-----------------------------------------------------------------------------
extern "C" gpuSpecs GetGpuSpecs(int reckless)
{
  int i, ndev, gpucount, seldev, nerrmsg;
#if !defined(AMBER_PLATFORM_AMD)
  unsigned int nvmlItemCount;
  int* validGpus;
  int* gpuList;
  nvmlReturn_t sttNV;
  nvmlProcessInfo_t nvmlInfo[32];
  nvmlDevice_t ntdev;
#endif
  cudaError_t stt;
  cudaDeviceProp devPRP;
  gpuSpecs devspc;
  // Test that there is a GPU in the system
  stt = cudaGetDeviceCount(&gpucount);

  // Unsure why there may be a problem with some machines (with no viable CUDA
  // devices) creating voluminous output, but it may be something wrong with
  // the gpucount variable returned from cudaGetDeviceCount.  Exit with an
  // error if it appears to be too large.
  if (gpucount > 64) {
    printf("mdgx >> Error.  A huge number of CUDA-capable devices were "
           "found.\nmdgx >> this is interpreted as a problem with the "
           "environment.\n");
    cudaDeviceReset();
    exit(1);
  }
  if (gpucount == 0) {
    printf("mdgx >> Error.  No CUDA-capable devices were found.\n");
    cudaDeviceReset();
    exit(1);
  }

  // Activate zero-copy
  cudaSetDeviceFlags(cudaDeviceMapHost);

  // Initialize the NVIDIA Management Library
#if !defined(AMBER_PLATFORM_AMD)
  nvmlInit();

  // Get device properties
  validGpus = (int*)malloc(gpucount * sizeof(int));
  gpuList = (int*)malloc(gpucount * sizeof(int));
  ndev = 0;
  nerrmsg = 0;
  for (i = 0; i < gpucount; i++) {
    cudaGetDeviceProperties(&devPRP, i);
    if (devPRP.major >= 3) {
      nvmlDeviceGetHandleByIndex(i, &ntdev);
      nvmlItemCount = 0;
      sttNV = nvmlDeviceGetComputeRunningProcesses(ntdev, &nvmlItemCount,
                                                   nvmlInfo);
      if (sttNV != NVML_SUCCESS && sttNV != NVML_ERROR_INSUFFICIENT_SIZE) {
	printf("mdgx >> Warning.  Unable to monitor activity on GPU %d "
	       "[error %u]\n", i, sttNV);
        nerrmsg++;
        if (nerrmsg >= 64) {
          printf("mdgx >> Too many warnings have been printed.  There must "
                 "be something wrong\nmdgx >> with the environment.\n");
          cudaDeviceReset();
          exit(1);
        }
      }
      if (nvmlItemCount == 0 || reckless == 1) {
        validGpus[i] = 1;
        gpuList[ndev] = i;
        ndev++;
      }
    }
  }
  if (ndev == 0 && reckless == 0) {
    printf("mdgx >> All GPUs are unavailable, or assisting other customers.  "
           "If you believe\nmdgx >> you have received this message in error "
	   "then you may re-run mdgx with the\nmdgx >> -Reckless flag.  As "
           "the name implies, you should be careful.\n");
    exit(1);
  }

  // Shut down the NVIDIA Management Lbirary
  nvmlShutdown();
  // Select a device from the list
  stt = cudaSetValidDevices(gpuList, ndev);
  if (stt != cudaSuccess) {
    printf("mdgx >> Error searching for CUDA-compatible GPU.\n");
    cudaDeviceReset();
    exit(1);
  }
#endif

  // Establish the CUDA context
  stt = cudaFree(0);
  if (stt != cudaSuccess) {
    printf("mdgx >> Error selecting compatible GPU.\n");
    cudaDeviceReset();
    exit(1);
  }

#if !defined(AMBER_PLATFORM_AMD)
  // Get the device
  stt = cudaGetDevice(&seldev);
  if (stt != cudaSuccess) {
    printf("mdgx >> Error setting cuda device.\n");
    cudaDeviceReset();
    exit(1);
  }
#else
  seldev = 0; // just for reading the properties
#endif
  cudaDeviceSynchronize();
  cudaGetDeviceProperties(&devPRP, seldev);

  // Copy the relevant information for shipment back to the calling function
  devspc.major          = devPRP.major;
  devspc.minor          = devPRP.minor;
  devspc.MPcount        = devPRP.multiProcessorCount;
  devspc.maxThrPerMP    = devPRP.maxThreadsPerMultiProcessor;
  devspc.maxThrPerBlock = devPRP.maxThreadsPerBlock;
  devspc.cardMemory     = devPRP.totalGlobalMem;
  i = strlen(devPRP.name);
  if (i > 127) {
    i = 127;
  }
  strncpy(devspc.name, devPRP.name, i);
  devspc.name[i] = '\0';

#if !defined(AMBER_PLATFORM_AMD)
  // Free allocated memory
  free(gpuList);
  free(validGpus);
#endif

  return devspc;
}

//-----------------------------------------------------------------------------
// kGpuPRNGSetup: kernel for initializing GPU random number generators
//-----------------------------------------------------------------------------
__global__ void kGpuPRNGSetup(curandState_t *states, int igseed)
{
  int tid = threadIdx.x + (blockIdx.x * blockDim.x);
  curand_init(igseed, tid, 0, &states[tid]);
}

//-----------------------------------------------------------------------------
// InitGpuPRNG: initialize pseudo-random number generators on the GPU.
//
// Arguments:
//   gms:       the repository for all parameters amd coordinates
//   igseed:    the random number generator seed
//   nblocks:   the number of blocks that the main dynamics kernels will run
//   blockDim:  dimension of the main dynamics kernel blocks
//-----------------------------------------------------------------------------
extern "C" void InitGpuPRNG(gpuMultiSim *gms, int igseed, int nblocks,
                            int blockdim)
{
  cudaMalloc((void **)&gms->prngStates,
             nblocks * blockdim * sizeof(curandState));
  kGpuPRNGSetup<<<nblocks, blockdim>>>((curandState_t*)gms->prngStates,
                                       igseed);
}

//-----------------------------------------------------------------------------
// SetGmsImage: function to establish a GPU Multi-Simulator on the device,
//              with pointers to all of the device-allocated memory as well as
//              constants describing the simulation conditions.
//
// Arguments:
//   gms:    the repository for all parameters amd coordinates
//-----------------------------------------------------------------------------
extern "C" void SetGmsImage(gpuMultiSim *gms)
{
  cudaError_t status;

  status = cudaMemcpyToSymbol(cGms, gms, sizeof(gpuMultiSim));
  if (status != cudaSuccess) {
    printf("SetGmsImage >> Unable to copy gpuMultiSim struct to the "
           "device (error %d).\n", (int)status);
    exit(1);
  }
}

//----------------------------------------------------------------------------
// kSetSystemCounters: set counters to guide the blocks as they step through
//                     systems during each segment of dynamics, in between
//                     coordinate writes.
//-----------------------------------------------------------------------------
__global__ void kSetSystemCounters(int blocks)
{
  int tid = threadIdx.x + (blockIdx.x * blockDim.x);
  while (tid < 2 * cGms.nsgmdout) {
    cGms.DVCsystemPos[tid] = blocks;
    tid += gridDim.x * blockDim.x;
  }
}

//-----------------------------------------------------------------------------
// Dynamics kernels with RATTLE to compute forces (kDynLoop) or force and
// energies (kEStep).  The concept is that kEStep will launch to get forces and
// energies for the first step, then kDynLoop will fire off for (ntpr - 1)
// steps, then kEStep will launch again, and so on until the maximum number of
// steps has been reached.
//-----------------------------------------------------------------------------
#ifdef AMBER_PLATFORM_AMD
#  define __LAUNCH_BOUNDS__(X, Y) __launch_bounds__(X)
#else
#  define __LAUNCH_BOUNDS__(X, Y) __launch_bounds__(X, Y)
#endif


#define GO_RATTLE
#define ATOM_LIMIT SM_ATOM_COUNT
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 750)
#  define THREAD_COUNT 256
#  define BLOCK_COUNT  4
#elif defined(AMBER_PLATFORM_AMD)
#  define THREAD_COUNT 128
#  define BLOCK_COUNT  8
#else
#  define THREAD_COUNT 288
#  define BLOCK_COUNT  4
#endif
#define SMALL_KERNEL
#define COMPUTE_ENERGY
__global__ void __LAUNCH_BOUNDS__(THREAD_COUNT, BLOCK_COUNT) kEStepSmallRtt(int sgc)
#include INC_HEADER
#undef COMPUTE_ENERGY
__global__ void __LAUNCH_BOUNDS__(THREAD_COUNT, BLOCK_COUNT) kDynLoopSmallRtt(int sgc)
#include INC_HEADER
#define GBSOLVENT
#define COMPUTE_ENERGY
__global__ void __LAUNCH_BOUNDS__(THREAD_COUNT, BLOCK_COUNT) kEStepSmallGBRtt(int sgc)
#include INC_HEADER
#undef COMPUTE_ENERGY
__global__ void __LAUNCH_BOUNDS__(THREAD_COUNT, BLOCK_COUNT) kDynLoopSmallGBRtt(int sgc)
#include INC_HEADER
#undef GBSOLVENT
#undef THREAD_COUNT
#undef BLOCK_COUNT
#undef ATOM_LIMIT
#undef SMALL_KERNEL

#define ATOM_LIMIT MD_ATOM_COUNT
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 750)
#  define THREAD_COUNT 512
#  define BLOCK_COUNT  2
#elif defined(AMBER_PLATFORM_AMD)
#  define THREAD_COUNT 256
#  define BLOCK_COUNT  4
#else
#  define THREAD_COUNT 576
#  define BLOCK_COUNT  2
#endif
#define COMPUTE_ENERGY
__global__ void __LAUNCH_BOUNDS__(THREAD_COUNT, BLOCK_COUNT) kEStepMedRtt(int sgc)
#include INC_HEADER
#undef COMPUTE_ENERGY
__global__ void __LAUNCH_BOUNDS__(THREAD_COUNT, BLOCK_COUNT) kDynLoopMedRtt(int sgc)
#include INC_HEADER
#define GBSOLVENT
#define COMPUTE_ENERGY
__global__ void __LAUNCH_BOUNDS__(THREAD_COUNT, BLOCK_COUNT) kEStepMedGBRtt(int sgc)
#include INC_HEADER
#undef COMPUTE_ENERGY
__global__ void __LAUNCH_BOUNDS__(THREAD_COUNT, BLOCK_COUNT) kDynLoopMedGBRtt(int sgc)
#include INC_HEADER
#undef GBSOLVENT
#undef THREAD_COUNT
#undef BLOCK_COUNT
#undef ATOM_LIMIT

#define ATOM_LIMIT LG_ATOM_COUNT
#if defined(AMBER_PLATFORM_AMD)
#  define THREAD_COUNT 512
#  define BLOCK_COUNT  2
#else
#  define THREAD_COUNT 1024
#  define BLOCK_COUNT  1
#endif
#define COMPUTE_ENERGY
__global__ void __LAUNCH_BOUNDS__(THREAD_COUNT, BLOCK_COUNT) kEStepLargeRtt(int sgc)
#include INC_HEADER
#undef COMPUTE_ENERGY
__global__ void __LAUNCH_BOUNDS__(THREAD_COUNT, BLOCK_COUNT) kDynLoopLargeRtt(int sgc)
#include INC_HEADER
#define GBSOLVENT
#define COMPUTE_ENERGY
__global__ void __LAUNCH_BOUNDS__(THREAD_COUNT, BLOCK_COUNT) kEStepLargeGBRtt(int sgc)
#include INC_HEADER
#undef COMPUTE_ENERGY
__global__ void __LAUNCH_BOUNDS__(THREAD_COUNT, BLOCK_COUNT) kDynLoopLargeGBRtt(int sgc)
#include INC_HEADER
#undef GBSOLVENT
#undef THREAD_COUNT
#undef BLOCK_COUNT
#undef ATOM_LIMIT
#undef GO_RATTLE

//-----------------------------------------------------------------------------
// Kernels without RATTLE.  The concept behind DynLoop and EStep kernels is the
// same as above, but without the register burden of RATTLE the kernels can
// engage additional threads.
//-----------------------------------------------------------------------------
#define ATOM_LIMIT SM_ATOM_COUNT
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 750)
#  define THREAD_COUNT 256
#  define BLOCK_COUNT  4
#elif defined(AMBER_PLATFORM_AMD)
#  define THREAD_COUNT 128
#  define BLOCK_COUNT  8
#else
#  define THREAD_COUNT 320
#  define BLOCK_COUNT  4
#endif
#define COMPUTE_ENERGY
__global__ void __LAUNCH_BOUNDS__(THREAD_COUNT, BLOCK_COUNT) kEStepSmall(int sgc)
#include INC_HEADER
#undef COMPUTE_ENERGY
__global__ void __LAUNCH_BOUNDS__(THREAD_COUNT, BLOCK_COUNT) kDynLoopSmall(int sgc)
#include INC_HEADER
#define GBSOLVENT
#define COMPUTE_ENERGY
__global__ void __LAUNCH_BOUNDS__(THREAD_COUNT, BLOCK_COUNT) kEStepSmallGB(int sgc)
#include INC_HEADER
#undef COMPUTE_ENERGY
__global__ void __LAUNCH_BOUNDS__(THREAD_COUNT, BLOCK_COUNT) kDynLoopSmallGB(int sgc)
#include INC_HEADER
#undef GBSOLVENT
#undef THREAD_COUNT
#undef BLOCK_COUNT
#undef ATOM_LIMIT

#define ATOM_LIMIT MD_ATOM_COUNT
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 750)
#  define THREAD_COUNT 512
#  define BLOCK_COUNT  2
#elif defined(AMBER_PLATFORM_AMD)
#  define THREAD_COUNT 256
#  define BLOCK_COUNT  4
#else
#  define THREAD_COUNT 640
#  define BLOCK_COUNT  2
#endif
#define COMPUTE_ENERGY
__global__ void __LAUNCH_BOUNDS__(THREAD_COUNT, BLOCK_COUNT) kEStepMed(int sgc)
#include INC_HEADER
#undef COMPUTE_ENERGY
__global__ void __LAUNCH_BOUNDS__(THREAD_COUNT, BLOCK_COUNT) kDynLoopMed(int sgc)
#include INC_HEADER
#define GBSOLVENT
#define COMPUTE_ENERGY
__global__ void __LAUNCH_BOUNDS__(THREAD_COUNT, BLOCK_COUNT) kEStepMedGB(int sgc)
#include INC_HEADER
#undef COMPUTE_ENERGY
__global__ void __LAUNCH_BOUNDS__(THREAD_COUNT, BLOCK_COUNT) kDynLoopMedGB(int sgc)
#include INC_HEADER
#undef GBSOLVENT
#undef THREAD_COUNT
#undef BLOCK_COUNT
#undef ATOM_LIMIT

#define ATOM_LIMIT LG_ATOM_COUNT
#if defined(AMBER_PLATFORM_AMD)
#  define THREAD_COUNT 512
#  define BLOCK_COUNT  2
#else
#  define THREAD_COUNT 1024
#  define BLOCK_COUNT  1
#endif
#define COMPUTE_ENERGY
__global__ void __LAUNCH_BOUNDS__(THREAD_COUNT, BLOCK_COUNT) kEStepLarge(int sgc)
#include INC_HEADER
#undef COMPUTE_ENERGY
__global__ void __LAUNCH_BOUNDS__(THREAD_COUNT, BLOCK_COUNT) kDynLoopLarge(int sgc)
#include INC_HEADER
#define GBSOLVENT
#define COMPUTE_ENERGY
__global__ void __LAUNCH_BOUNDS__(THREAD_COUNT, BLOCK_COUNT) kEStepLargeGB(int sgc)
#include INC_HEADER
#undef COMPUTE_ENERGY
__global__ void __LAUNCH_BOUNDS__(THREAD_COUNT, BLOCK_COUNT) kDynLoopLargeGB(int sgc)
#include INC_HEADER
#undef GBSOLVENT
#undef THREAD_COUNT
#undef BLOCK_COUNT
#undef ATOM_LIMIT

//-----------------------------------------------------------------------------
// LaunchDynamics: launch the appropriate kernels for energy and forces.  As
//                 is done in the pmemd code, this function in the CUDA unit
//                 encapsualtes the launch so that the .c libraries can be
//                 built with a standard C compiler.
//
// Arguments:
//   gms:       the repository for all parameters amd coordinates
//   blockDim:  the block size to use, determined by PlanGpuUtilization in
//              Peptide.c
//   devspc:    device specifications
//-----------------------------------------------------------------------------
extern "C" void LaunchDynamics(gpuMultiSim *gms, int blockDim, int nblocks,
			       gpuSpecs *devspc)
{
  int i;

  // Initialize system counters for this portion of dynamics
  kSetSystemCounters<<<nblocks, blockDim>>>(nblocks);

#if defined(AMBER_PLATFORM_AMD)
  int medBlockDim = 256;
  int largeBlockDim = 512;
#else
  int medBlockDim = 512;
  int largeBlockDim = 1024;
#endif

#ifdef AMBER_PLATFORM_AMD
  nblocks = gms->nsys;
#endif

  // Vacuum-phase dynamics
  if (gms->igb == 6) {
    if (gms->rattle == 0) {
      for (i = 0; i < gms->nsgmdout; i++) {
        if (blockDim < medBlockDim) {
          kEStepSmall<<<nblocks, blockDim>>>(i);
          kDynLoopSmall<<<nblocks, blockDim>>>(i);
        }
        else if (blockDim < largeBlockDim) {
          kEStepMed<<<nblocks, blockDim>>>(i);
          kDynLoopMed<<<nblocks, blockDim>>>(i);
        }
        else {
          kEStepLarge<<<nblocks, blockDim>>>(i);
          kDynLoopLarge<<<nblocks, blockDim>>>(i);
        }
      }
    }
    else {
      for (i = 0; i < gms->nsgmdout; i++) {
        if (blockDim < medBlockDim) {
          kEStepSmallRtt<<<nblocks, blockDim>>>(i);
          kDynLoopSmallRtt<<<nblocks, blockDim>>>(i);
        }
        else if (blockDim < largeBlockDim) {
          kEStepMedRtt<<<nblocks, blockDim>>>(i);
          kDynLoopMedRtt<<<nblocks, blockDim>>>(i);
        }
        else {
          kEStepLargeRtt<<<nblocks, blockDim>>>(i);
          kDynLoopLargeRtt<<<nblocks, blockDim>>>(i);
        }
      }
    }
  }

  // Dynamics in Generalized Born solvent
  else {
    if (gms->rattle == 0) {
      for (i = 0; i < gms->nsgmdout; i++) {
        if (blockDim < medBlockDim) {
          kEStepSmallGB<<<nblocks, blockDim>>>(i);
          kDynLoopSmallGB<<<nblocks, blockDim>>>(i);
        }
        else if (blockDim < largeBlockDim) {
          kEStepMedGB<<<nblocks, blockDim>>>(i);
          kDynLoopMedGB<<<nblocks, blockDim>>>(i);
        }
        else {
          kEStepLargeGB<<<nblocks, blockDim>>>(i);
          kDynLoopLargeGB<<<nblocks, blockDim>>>(i);
	}
      }
    }
    else {
      for (i = 0; i < gms->nsgmdout; i++) {
        if (blockDim < medBlockDim) {
          kEStepSmallGBRtt<<<nblocks, blockDim>>>(i);
          kDynLoopSmallGBRtt<<<nblocks, blockDim>>>(i);
        }
        else if (blockDim < largeBlockDim) {
          kEStepMedGBRtt<<<nblocks, blockDim>>>(i);
          kDynLoopMedGBRtt<<<nblocks, blockDim>>>(i);
        }
        else {
          kEStepLargeGBRtt<<<nblocks, blockDim>>>(i);
          kDynLoopLargeGBRtt<<<nblocks, blockDim>>>(i);
	}
      }
    }
  }
}
